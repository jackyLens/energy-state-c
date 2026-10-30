/*
 * energy.c -- energy state model (ODE + EKF), C implementation.
 *
 * Mapping to the Python reference `energy/src/energy_state/` (bit-level parity held up to
 * v1.1.1; the v1.2 additions -- personal observation baselines, gap policy, online load
 * calibration, soft gating -- exist only in C for now):
 *   params.py      -> energy_params_t + energy_phi
 *   ode.py         -> energy_advance / 2x2 Jacobians / analytic discretisation
 *   load.py        -> load_raw_series / energy_load_raw
 *   observation.py -> energy_obs_t + obs_regress (+ personal intercepts, v1.2)
 *   ekf.py         -> ekf_predict_step / ctx_update (robust t, gating, Joseph form)
 *   product.py     -> energy_display / smoothing / daily summary
 */

#include "energy.h"

#include <limits.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define ENERGY_HOURS_PER_STEP (5.0 / 60.0)
#define ENERGY_MAX_SUBSTEP_H (2.5 / 60.0)
#define ENERGY_BINS 288
#define ENERGY_MAX_HOURS (1.0e6)          /* sane time-axis range (~114 years) */
#define ENERGY_MAX_SUBSTEPS 4096          /* sub-step cap for one interval (CPU guard) */
#define ENERGY_MAX_SMOOTH_MIN 1440.0      /* display smoothing cap (1 day) */
#define ENERGY_MAX_SAMPLES 10000000u      /* batch size cap (size_t overflow guard) */

#define ENERGY_DEFAULT_STEPS_REF 500.0    /* ~100 steps/min over 5 min; brisk walking is 500-600 */
#define ENERGY_EKF_MAX_STEP_H 0.5         /* EKF prediction sub-step for long intervals */
#define ENERGY_TRIMP_MAX_DT_H 0.5         /* TRIMP coverage clamp (batch and stream alike) */
#define ENERGY_STRESS_AUTO_MAX 1.5        /* auto stress scale: a value above this means 0-100 */
#define ENERGY_EPISODE_SWITCH 6           /* 30 min of the other regime ends an episode */
#define ENERGY_MIN_NIGHT_SAMPLES 36       /* a night needs >= 3 h of samples to be committed */
#define ENERGY_MIN_DAY_SAMPLES 36         /* a day needs >= 3 h of samples to be committed */
#define ENERGY_EPISODE_MAX 288            /* force-commit an episode after 24 h of one regime */
#define ENERGY_EPISODE_BREAK_H 6.0        /* a data gap longer than this closes open episodes */
#define ENERGY_MAX_GAP_SEGMENTS 120       /* sleep-window segments walked inside one gap (60 days) */
#define ENERGY_LREF_BIN_W 0.125           /* online l_ref histogram bin width (raw units) */
#define ENERGY_LREF_PRIOR 2.0             /* prior l_ref of a fresh stream (typical batch P99.5) */
#define ENERGY_LREF_PRIOR_W 288.0         /* prior weight in samples (one day) */
#define ENERGY_LREF_MEMORY 4032.0         /* histogram is halved beyond ~14 days of samples */
#define ENERGY_PEAK_DECAY_PER_H (1.0 / 24.0) /* online hr_max peak decays 1 bpm per day */

/* ------------------------------------------------------------------ */
/* 小工具                                                              */
/* ------------------------------------------------------------------ */
static double clampd(double v, double lo, double hi)
{
    if (!isfinite(v)) return lo;
    return v < lo ? lo : (v > hi ? hi : v);
}

static double finite_or(double v, double fallback)
{
    return isfinite(v) ? v : fallback;
}

static int cmp_double(const void *pa, const void *pb)
{
    const double a = *(const double *)pa;
    const double b = *(const double *)pb;
    return (a < b) ? -1 : (a > b ? 1 : 0);
}

/* 线性插值分位点，p ∈ [0,100]；会排序 tmp（调用方提供副本）。 */
static double percentile_sorted(double *tmp, size_t n, double p)
{
    if (n == 0) return (double)NAN;
    if (n == 1) return tmp[0];
    qsort(tmp, n, sizeof(double), cmp_double);
    const double pos = (double)(n - 1) * clampd(p, 0.0, 100.0) / 100.0;
    const size_t lo = (size_t)floor(pos);
    const size_t hi = (size_t)ceil(pos);
    if (lo == hi) return tmp[lo];
    const double f = pos - (double)lo;
    return tmp[lo] * (1.0 - f) + tmp[hi] * f;
}

/* 前向/后向线性补缺：两侧有值做线性插值，只有一侧则沿用最近的有限值。 */
static void fill_linear(double *v, size_t n)
{
    if (n == 0) return;
    size_t first = n, last = 0, good = 0;
    for (size_t i = 0; i < n; i++) {
        if (isfinite(v[i])) {
            if (first == n) first = i;
            last = i;
            good++;
        }
    }
    if (good == 0) {
        for (size_t i = 0; i < n; i++) v[i] = 0.0;
        return;
    }
    for (size_t i = 0; i < first; i++) v[i] = v[first];
    for (size_t i = last + 1; i < n; i++) v[i] = v[last];
    size_t i = first;
    while (i <= last) {
        if (isfinite(v[i])) {
            i++;
            continue;
        }
        size_t j = i;
        while (j <= last && !isfinite(v[j])) j++;
        const double a = v[i - 1], b = v[j];
        const double span = (double)(j - (i - 1));
        for (size_t k = i; k < j; k++) {
            const double w = (double)(k - (i - 1)) / span;
            v[k] = a * (1.0 - w) + b * w;
        }
        i = j;
    }
}

static double wrap24(double h)
{
    double r = fmod(h, 24.0);
    if (r < 0.0) r += 24.0;
    return r;
}

/* 钟点 → [0,288) 分箱下标；非有限值按 0 处理，避免 double→int 的未定义行为 */
static int tod_bin(double tod)
{
    if (!isfinite(tod)) return 0;
    int b = (int)(wrap24(tod) / 24.0 * ENERGY_BINS);
    if (b < 0) b = 0;
    if (b >= ENERGY_BINS) b = ENERGY_BINS - 1;
    return b;
}

/* ------------------------------------------------------------------ */
/* Defaults                                                            */
/* ------------------------------------------------------------------ */
double energy_hrv_log_from_ms(double rmssd_ms)
{
    if (!isfinite(rmssd_ms) || rmssd_ms <= 0.0) return (double)NAN;
    return log(rmssd_ms);
}

void energy_sample_init(energy_sample_t *s)
{
    if (!s) return;
    memset(s, 0, sizeof(*s));
    s->hr = NAN;
    s->hrv_log = NAN;
    s->q_hr = 1.0;
    s->q_hrv = 1.0;
    s->load = -1.0;
    s->activity = -1.0;
    s->spo2 = NAN;
    s->spo2_min = NAN;
}

void energy_params_default(energy_params_t *p)
{
    if (!p) return;
    memset(p, 0, sizeof(*p));
    p->kappa_w = 0.100;
    p->kappa_s[ENERGY_WAKE] = 0.050;
    p->kappa_s[ENERGY_LIGHT] = 0.260;
    p->kappa_s[ENERGY_DEEP] = 0.420;
    p->kappa_s[ENERGY_REM] = 0.180;
    p->alpha_sleep = 0.300;
    p->alpha_rest = 0.020;
    p->w_stage[ENERGY_WAKE] = 0.15;
    p->w_stage[ENERGY_LIGHT] = 1.00;
    p->w_stage[ENERGY_DEEP] = 1.25;
    p->w_stage[ENERGY_REM] = 0.95;
    p->eta = 0.60;
    p->s_ref = 0.50;
    p->s_floor = 0.10;
    p->phi_min = 0.30;
    p->phi_max = 2.00;
    p->beta_0 = 0.010;
    p->beta_L = 0.150;
    p->beta_P = 0.030;
    p->beta_S = 0.035;
}

void energy_load_default(energy_load_cfg_t *c)
{
    if (!c) return;
    memset(c, 0, sizeof(*c)); /* hr_base[] = 0 -> automatic */
    c->b_trimp = 1.92;
    c->w_trimp = 1.0;
    c->w_steps = 0.6;
    c->w_dhr = 0.4;
    c->trimp_scale = 10.0;
    c->steps_ref = ENERGY_DEFAULT_STEPS_REF;
    c->l_ref = -1.0;   /* automatic */
    c->hr_rest = -1.0; /* automatic */
    c->hr_max = -1.0;  /* automatic */
}

void energy_display_default(energy_display_cfg_t *c)
{
    if (!c) return;
    memset(c, 0, sizeof(*c));
    c->rho = 0.10;
    c->peak_hour = 10.0;
    c->dip_hour = 14.5;
    c->dip_depth = 0.70;
    c->dip_width_h = 1.6;
    c->smooth_minutes = 30.0;
    c->band_z = 1.96;
}

void energy_ekf_default(energy_ekf_cfg_t *c)
{
    if (!c) return;
    memset(c, 0, sizeof(*c));
    c->q_s = 0.02 * 0.02;
    c->q_e = 0.03 * 0.03;
    c->p0_s = 0.18 * 0.18;
    c->p0_e = 0.22 * 0.22;
    c->gate_sigma = 4.0;
    c->q_floor = 0.05;
    c->robust_nu = 4.0;
    c->gate_max_reject = 6;
}

void energy_obs_default(energy_obs_t *o)
{
    if (!o) return;
    memset(o, 0, sizeof(*o));
    /* Column 0 (intercept) is only a population prior; with baseline.enable=1 (default) it is
     * replaced by the user's own learned level, so only the slopes matter. */
    /* night: [1, E, p_deep, p_rem, p_wake] */
    const double hr_n[5] = {58.0, -12.0, -5.0, 2.0, 9.0};
    const double hrv_n[5] = {3.55, 0.45, 0.30, 0.18, -0.35};
    /* day: [1, E, A, cos1, sin1, cos2, sin2]
     * A = 1 is roughly brisk walking (~550 steps / 5 min) -> about +45 bpm and -0.35 ln(RMSSD)
     * over rest; v1.1 used +105 bpm / -0.55 with A saturating at 120 steps / 5 min. */
    const double hr_d[7] = {80.0, -12.0, 45.0, -2.0, -3.46, 0.75, -1.30};
    const double hrv_d[7] = {3.25, 0.45, -0.35, -0.06, -0.10, 0.025, -0.043};
    memcpy(o->hr_night, hr_n, sizeof(hr_n));
    memcpy(o->hrv_night, hrv_n, sizeof(hrv_n));
    memcpy(o->hr_day, hr_d, sizeof(hr_d));
    memcpy(o->hrv_day, hrv_d, sizeof(hrv_d));
    o->sigma_hr_night = 2.5;
    o->sigma_hr_day = 7.0;
    o->sigma_hrv_night = 0.22;
    o->sigma_hrv_day = 0.40;
}

void energy_baseline_default(energy_baseline_cfg_t *c)
{
    if (!c) return;
    memset(c, 0, sizeof(*c)); /* init: nothing learned yet (all counts 0) */
    c->enable = 1;
    c->cold_nights = 3;
    c->lambda = 0.10;
    c->cold_r_inflate = 4.0;
    c->bootstrap_samples = 6;
    c->min_samples = 12;
    c->rest_activity = 0.10;
    c->anomaly[ENERGY_BL_HR_NIGHT] = 6.0;
    c->anomaly[ENERGY_BL_HRV_NIGHT] = 0.25;
    c->anomaly[ENERGY_BL_HR_DAY] = 8.0;
    c->anomaly[ENERGY_BL_HRV_DAY] = 0.35;
    c->anomaly_weight = 0.25;
}

void energy_gap_default(energy_gap_cfg_t *c)
{
    if (!c) return;
    memset(c, 0, sizeof(*c));
    c->max_hold_h = 1.0;
    c->sleep_start = 23.0;
    c->sleep_end = 7.0;
    c->template_stage[ENERGY_WAKE] = 0.10;
    c->template_stage[ENERGY_LIGHT] = 0.50;
    c->template_stage[ENERGY_DEEP] = 0.20;
    c->template_stage[ENERGY_REM] = 0.20;
}

void energy_config_default(energy_config_t *c)
{
    if (!c) return;
    memset(c, 0, sizeof(*c));
    energy_params_default(&c->params);
    energy_load_default(&c->load);
    energy_display_default(&c->display);
    energy_ekf_default(&c->ekf);
    energy_obs_default(&c->obs);
    energy_baseline_default(&c->baseline);
    energy_gap_default(&c->gap);
    c->stress_scale = 0.0; /* automatic */
    c->x0_s = 0.45;
    c->x0_e = 0.55;
    c->use_ekf = 1;
    c->day_update = 1;
    c->start_year = 0;
    c->start_month = 0;
    c->start_day = 0;
}

const char *energy_version(void)
{
    return ENERGY_VERSION;
}

static int cfg_fail(char *buf, size_t len, const char *msg)
{
    if (buf && len > 0) snprintf(buf, len, "%s", msg);
    return ENERGY_EINVAL;
}

static int baseline_valid(const energy_baseline_t *b)
{
    for (int k = 0; k < 4; k++) {
        if (b->n[k] < 0) return 0;
        if (b->n[k] > 0 && !isfinite(b->level[k])) return 0;
    }
    if (b->n_sleep < 0) return 0;
    if (b->n_sleep > 0) {
        if (!(b->sleep_on >= 0.0 && b->sleep_on < 24.0)) return 0;
        if (!(b->sleep_off >= 0.0 && b->sleep_off < 24.0)) return 0;
        double sum = 0.0;
        for (int i = 0; i < 4; i++) {
            if (!(b->stage_mix[i] >= 0.0 && b->stage_mix[i] <= 1.0)) return 0;
            sum += b->stage_mix[i];
        }
        if (sum > 1.0 + 1e-6) return 0;
    }
    return 1;
}

static int all_finite(const double *v, int n)
{
    for (int i = 0; i < n; i++) {
        if (!isfinite(v[i])) return 0;
    }
    return 1;
}

int energy_config_validate(const energy_config_t *cfg, char *errbuf, size_t errlen)
{
    if (!cfg) return cfg_fail(errbuf, errlen, "cfg 为 NULL");
    const energy_params_t *p = &cfg->params;
    if (!(p->kappa_w > 0.0)) return cfg_fail(errbuf, errlen, "kappa_w 必须 > 0");
    for (int i = 0; i < 4; i++) {
        if (!(p->kappa_s[i] > 0.0)) return cfg_fail(errbuf, errlen, "kappa_s[*] 必须 > 0");
    }
    if (!(p->alpha_sleep > 0.0)) return cfg_fail(errbuf, errlen, "alpha_sleep 必须 > 0");
    if (!(p->alpha_rest >= 0.0)) return cfg_fail(errbuf, errlen, "alpha_rest 必须 >= 0");
    if (!(p->beta_0 >= 0.0) || !(p->beta_L >= 0.0) || !(p->beta_P >= 0.0) || !(p->beta_S >= 0.0)) {
        return cfg_fail(errbuf, errlen, "beta_* 必须 >= 0");
    }
    if (!(p->eta > 0.0)) return cfg_fail(errbuf, errlen, "eta 必须 > 0");
    if (!(p->s_ref > 0.0)) return cfg_fail(errbuf, errlen, "s_ref 必须 > 0");
    if (!(p->s_floor >= 0.0)) return cfg_fail(errbuf, errlen, "s_floor 必须 >= 0");
    if (!(p->phi_min >= 0.0) || !(p->phi_min < p->phi_max) || !isfinite(p->phi_max)) {
        return cfg_fail(errbuf, errlen, "phi_min/phi_max 非法（需 0 <= min < max 且有限）");
    }
    if (!(cfg->x0_s >= 0.0 && cfg->x0_s <= 1.0 && cfg->x0_e >= 0.0 && cfg->x0_e <= 1.0)) {
        return cfg_fail(errbuf, errlen, "x0_s/x0_e 必须在 [0,1]");
    }
    const energy_ekf_cfg_t *e = &cfg->ekf;
    if (!(e->q_s > 0.0) || !(e->q_e > 0.0)) return cfg_fail(errbuf, errlen, "q_s/q_e 必须 > 0");
    if (!(e->p0_s > 0.0) || !(e->p0_e > 0.0)) return cfg_fail(errbuf, errlen, "p0_s/p0_e 必须 > 0");
    if (!(e->gate_sigma > 0.0)) return cfg_fail(errbuf, errlen, "gate_sigma 必须 > 0");
    if (!(e->q_floor > 0.0)) return cfg_fail(errbuf, errlen, "q_floor 必须 > 0");
    if (!(e->robust_nu >= 0.0)) return cfg_fail(errbuf, errlen, "robust_nu 必须 >= 0");
    if (e->gate_max_reject < 0) return cfg_fail(errbuf, errlen, "gate_max_reject 必须 >= 0");
    const energy_display_cfg_t *d = &cfg->display;
    if (!isfinite(d->rho) || d->rho < 0.0 || d->rho > 1.0) {
        return cfg_fail(errbuf, errlen, "rho 必须在 [0,1]");
    }
    if (!isfinite(d->smooth_minutes) || d->smooth_minutes < 0.0) {
        return cfg_fail(errbuf, errlen, "smooth_minutes 必须 >= 0");
    }
    if (!(d->dip_width_h > 0.0) || !(d->band_z > 0.0)) {
        return cfg_fail(errbuf, errlen, "dip_width_h/band_z 必须 > 0");
    }
    const energy_load_cfg_t *l = &cfg->load;
    if (!(l->steps_ref > 0.0)) return cfg_fail(errbuf, errlen, "steps_ref 必须 > 0");
    if (!(l->trimp_scale >= 0.0) || !(l->b_trimp >= 0.0)) {
        return cfg_fail(errbuf, errlen, "trimp_scale/b_trimp 必须 >= 0");
    }
    for (int h = 0; h < 24; h++) {
        if (!(l->hr_base[h] >= 0.0) || !isfinite(l->hr_base[h])) {
            return cfg_fail(errbuf, errlen, "load.hr_base[*] 必须有限且 >= 0（0 = 自动）");
        }
    }
    const energy_obs_t *o = &cfg->obs;
    if (!(o->sigma_hr_night > 0.0) || !(o->sigma_hr_day > 0.0) ||
        !(o->sigma_hrv_night > 0.0) || !(o->sigma_hrv_day > 0.0)) {
        return cfg_fail(errbuf, errlen, "观测噪声 sigma 必须 > 0");
    }
    if (!all_finite(o->hr_night, 5) || !all_finite(o->hrv_night, 5) ||
        !all_finite(o->hr_day, 7) || !all_finite(o->hrv_day, 7)) {
        return cfg_fail(errbuf, errlen, "观测系数必须有限");
    }
    if (!(o->hr_night[1] <= 0.0 && o->hr_day[1] <= 0.0)) {
        return cfg_fail(errbuf, errlen, "心率对 E 的系数必须 <= 0（恢复越好心率越低）");
    }
    if (!(o->hrv_night[1] >= 0.0 && o->hrv_day[1] >= 0.0)) {
        return cfg_fail(errbuf, errlen, "HRV 对 E 的系数必须 >= 0");
    }
    const energy_baseline_cfg_t *b = &cfg->baseline;
    if (b->cold_nights < 0 || b->bootstrap_samples < 1 || b->min_samples < 1) {
        return cfg_fail(errbuf, errlen, "baseline: cold_nights >= 0，bootstrap_samples/min_samples >= 1");
    }
    if (!(b->lambda > 0.0 && b->lambda <= 1.0)) return cfg_fail(errbuf, errlen, "baseline.lambda 必须在 (0,1]");
    if (!(b->cold_r_inflate >= 1.0) || !isfinite(b->cold_r_inflate)) {
        return cfg_fail(errbuf, errlen, "baseline.cold_r_inflate 必须 >= 1");
    }
    if (!(b->rest_activity >= 0.0 && b->rest_activity <= 1.0)) {
        return cfg_fail(errbuf, errlen, "baseline.rest_activity 必须在 [0,1]");
    }
    for (int k = 0; k < 4; k++) {
        if (!(b->anomaly[k] > 0.0) || !isfinite(b->anomaly[k])) {
            return cfg_fail(errbuf, errlen, "baseline.anomaly[*] 必须 > 0");
        }
    }
    if (!(b->anomaly_weight >= 0.0 && b->anomaly_weight <= 1.0)) {
        return cfg_fail(errbuf, errlen, "baseline.anomaly_weight 必须在 [0,1]");
    }
    if (!baseline_valid(&b->init)) return cfg_fail(errbuf, errlen, "baseline.init 非法");
    const energy_gap_cfg_t *g = &cfg->gap;
    if (!(g->max_hold_h > 0.0) || !isfinite(g->max_hold_h)) {
        return cfg_fail(errbuf, errlen, "gap.max_hold_h 必须 > 0");
    }
    if (!(g->sleep_start >= 0.0 && g->sleep_start < 24.0 && g->sleep_end >= 0.0 && g->sleep_end < 24.0)) {
        return cfg_fail(errbuf, errlen, "gap.sleep_start/sleep_end 必须在 [0,24)");
    }
    {
        double sum = 0.0;
        for (int i = 0; i < 4; i++) {
            if (!(g->template_stage[i] >= 0.0 && g->template_stage[i] <= 1.0)) {
                return cfg_fail(errbuf, errlen, "gap.template_stage[*] 必须在 [0,1]");
            }
            sum += g->template_stage[i];
        }
        if (sum > 1.0 + 1e-9) return cfg_fail(errbuf, errlen, "gap.template_stage 之和必须 <= 1");
    }
    if (!(cfg->stress_scale >= 0.0) || !isfinite(cfg->stress_scale)) {
        return cfg_fail(errbuf, errlen, "stress_scale 必须 >= 0（0 = 自动）");
    }
    if (errbuf && errlen > 0) errbuf[0] = '\0';
    return ENERGY_OK;
}

/* ------------------------------------------------------------------ */
/* 状态方程                                                            */
/* ------------------------------------------------------------------ */
double energy_phi(const energy_params_t *p, double s)
{
    const double x = (fmax(0.0, s) + p->s_floor) / (p->s_ref + p->s_floor);
    return clampd(pow(x, p->eta), p->phi_min, p->phi_max);
}

static double phi_dphi(const energy_params_t *p, double s, double *dphi)
{
    const double sm = fmax(0.0, s);
    const double x = (sm + p->s_floor) / (p->s_ref + p->s_floor);
    const double raw = pow(x, p->eta);
    const double phi = clampd(raw, p->phi_min, p->phi_max);
    if (raw <= p->phi_min || raw >= p->phi_max) {
        *dphi = 0.0;
    } else {
        *dphi = p->eta * pow(sm + p->s_floor, p->eta - 1.0) / pow(p->s_ref + p->s_floor, p->eta);
    }
    return phi;
}

double energy_wake_equilibrium(const energy_params_t *p, double s, double load, double stress)
{
    const double alpha = p->alpha_rest;
    const double beta = p->beta_0 + p->beta_L * load + p->beta_P * stress + p->beta_S * s;
    return (alpha + beta) > 0.0 ? alpha / (alpha + beta) : 0.0;
}

double energy_sleep_equilibrium(const energy_params_t *p, double s, double psi)
{
    const double alpha = p->alpha_sleep * energy_phi(p, s) * psi;
    const double beta = p->beta_0;
    return (alpha + beta) > 0.0 ? alpha / (alpha + beta) : 1.0;
}

/*
 * 单步解析推进（与 Python advance_raw 逐行对应）。
 * 输入：m 睡眠概率、lam 分期消散、psi 分期恢复能力、load、stress、h（小时）。
 */
static void advance_raw(const energy_params_t *p,
                        double s, double e, double m, double lam, double psi,
                        double load, double stress, double h,
                        double *s_out, double *e_out)
{
    const double a = p->kappa_w * (1.0 - m);
    const double b = a + m * lam;
    double s1, s_mid;
    if (b > 0.0) {
        const double s_inf = a / b;
        const double d = s - s_inf;
        s1 = s_inf + d * exp(-b * h);
        s_mid = s_inf + d * exp(-0.5 * b * h);
    } else {
        s1 = s;
        s_mid = s;
    }

    const double alpha = (1.0 - m) * p->alpha_rest + m * p->alpha_sleep * energy_phi(p, s_mid) * psi;
    const double beta = p->beta_0 + p->beta_L * load + p->beta_P * stress +
                        (1.0 - m) * p->beta_S * s_mid;
    const double c = alpha + beta;
    double e1;
    if (c > 0.0) {
        const double e_inf = alpha / c;
        e1 = e_inf + (e - e_inf) * exp(-c * h);
    } else {
        e1 = e;
    }
    *s_out = s1;
    *e_out = e1;
}

void energy_advance(const energy_params_t *p,
                    double s, double e, double m, double lam, double psi,
                    double load, double stress, double h,
                    double *s_out, double *e_out)
{
    double s1, e1;
    advance_raw(p, s, e, m, lam, psi, load, stress, h, &s1, &e1);
    *s_out = clampd(s1, 0.0, 1.0);
    *e_out = clampd(e1, 0.0, 1.0);
}

/* 连续时间雅可比 ∂f/∂x（2×2，行优先） */
static void continuous_jacobian(const energy_params_t *p,
                                double s, double e, double m, double lam, double psi,
                                double load, double stress, double F[4])
{
    const double a = p->kappa_w * (1.0 - m);
    const double b = a + m * lam;
    double dphi;
    const double phi = phi_dphi(p, s, &dphi);
    const double alpha = (1.0 - m) * p->alpha_rest + m * p->alpha_sleep * phi * psi;
    const double beta = p->beta_0 + p->beta_L * load + p->beta_P * stress +
                        (1.0 - m) * p->beta_S * s;

    F[0] = -b;
    F[1] = 0.0;
    F[2] = m * p->alpha_sleep * psi * dphi * (1.0 - e) - (1.0 - m) * p->beta_S * e;
    F[3] = -(alpha + beta);
}

/* Φ = expm(F·h)，2×2 闭式（与 Python discrete_jacobian 一致） */
static void discrete_jacobian(const double F[4], double h, double Phi[4])
{
    const double A00 = F[0] * h, A01 = F[1] * h, A10 = F[2] * h, A11 = F[3] * h;
    const double tau = A00 + A11;
    const double det = A00 * A11 - A01 * A10;
    const double disc = 0.25 * tau * tau - det;
    const double N00 = A00 - 0.5 * tau, N01 = A01, N10 = A10, N11 = A11 - 0.5 * tau;
    double c0, c1;
    if (disc > 1e-12) {
        const double delta = sqrt(disc);
        c0 = exp(0.5 * tau) * cosh(delta);
        c1 = exp(0.5 * tau) * sinh(delta) / delta;
    } else if (disc < -1e-12) {
        const double omega = sqrt(-disc);
        c0 = exp(0.5 * tau) * cos(omega);
        c1 = exp(0.5 * tau) * sin(omega) / omega;
    } else {
        c0 = exp(0.5 * tau);
        c1 = c0;
    }
    Phi[0] = c0 + c1 * N00;
    Phi[1] = c1 * N01;
    Phi[2] = c1 * N10;
    Phi[3] = c0 + c1 * N11;
}

/* ------------------------------------------------------------------ */
/* Load / activity                                                      */
/* ------------------------------------------------------------------ */
double energy_activity_step(double steps, double kcal, double steps_ref)
{
    const double s = clampd(steps / (steps_ref > 0 ? steps_ref : ENERGY_DEFAULT_STEPS_REF), 0.0, 1.5);
    double a = 0.8 * s;
    if (isfinite(kcal) && kcal > 0.0) {
        a += 0.2 * clampd(kcal / 6.0, 0.0, 1.5);
    }
    return clampd(a, 0.0, 1.0);
}

double energy_load_raw(const energy_load_cfg_t *c,
                       double hr, double hr_rest, double hr_max, double hr_base,
                       double steps, double dt_h)
{
    const double span = fmax(hr_max - hr_rest, 1.0);
    const double step_frac =
        clampd(steps / (c->steps_ref > 0 ? c->steps_ref : ENERGY_DEFAULT_STEPS_REF), 0.0, 1.5);
    const double r = isfinite(hr) ? clampd((hr - hr_rest) / span, 0.0, 1.0) : 0.0;
    const double r_proxy = clampd(0.6 * step_frac, 0.0, 1.0);
    const double r_eff = isfinite(hr) ? r : r_proxy;
    const double trimp = dt_h * r_eff * exp(c->b_trimp * r_eff) * c->trimp_scale;
    const double dhr = isfinite(hr) ? clampd((hr - hr_base) / span, 0.0, 1.0) : 0.0;
    return c->w_trimp * trimp + c->w_steps * step_frac + c->w_dhr * dhr;
}

/* hr_rest / hr_max / time-of-day baseline from the raw series (as Python load.py).
 * All of them report ENERGY_ENOMEM instead of silently falling back to constants. */
static int estimate_hr_rest(const energy_sample_t *s, size_t n, double *out)
{
    double *buf = (double *)malloc(n * sizeof(double));
    if (!buf) return ENERGY_ENOMEM;
    size_t k = 0;
    for (size_t i = 0; i < n; i++) {
        if (isfinite(s[i].hr) && s[i].stage[ENERGY_DEEP] > 0.5) buf[k++] = s[i].hr;
    }
    if (k >= 12) {
        *out = percentile_sorted(buf, k, 50.0);
        free(buf);
        return ENERGY_OK;
    }
    k = 0;
    for (size_t i = 0; i < n; i++) {
        const double m = s[i].stage[ENERGY_LIGHT] + s[i].stage[ENERGY_DEEP] + s[i].stage[ENERGY_REM];
        if (isfinite(s[i].hr) && m > 0.5) buf[k++] = s[i].hr;
    }
    if (k >= 12) {
        *out = percentile_sorted(buf, k, 20.0);
        free(buf);
        return ENERGY_OK;
    }
    k = 0;
    for (size_t i = 0; i < n; i++) {
        if (isfinite(s[i].hr)) buf[k++] = s[i].hr;
    }
    *out = (k > 0) ? percentile_sorted(buf, k, 5.0) : 60.0;
    free(buf);
    return ENERGY_OK;
}

static int estimate_hr_max(const energy_sample_t *s, size_t n, double *out)
{
    double *buf = (double *)malloc(n * sizeof(double));
    if (!buf) return ENERGY_ENOMEM;
    size_t k = 0;
    for (size_t i = 0; i < n; i++) {
        if (isfinite(s[i].hr)) buf[k++] = s[i].hr;
    }
    *out = (k < 20) ? 185.0 : clampd(percentile_sorted(buf, k, 99.5), 150.0, 210.0);
    free(buf);
    return ENERGY_OK;
}

/* Time-of-day HR baseline: drop > P90 (exercise) -> per-bin median -> circular interpolation
 * -> 7-point circular smoothing. */
static int estimate_hr_base(const energy_sample_t *s, size_t n, double *prof)
{
    double *hr = (double *)malloc(n * sizeof(double));
    if (!hr) return ENERGY_ENOMEM;
    size_t k = 0;
    for (size_t i = 0; i < n; i++) {
        if (isfinite(s[i].hr)) hr[k++] = s[i].hr;
    }
    if (k < (size_t)ENERGY_BINS) {
        const double fill = (k > 0) ? percentile_sorted(hr, k, 50.0) : 70.0;
        for (int i = 0; i < ENERGY_BINS; i++) prof[i] = fill;
        free(hr);
        return ENERGY_OK;
    }
    const double thr = percentile_sorted(hr, k, 90.0);
    free(hr);

    int counts[ENERGY_BINS];
    int pos[ENERGY_BINS];
    int offsets[ENERGY_BINS + 1];
    double *values = (double *)malloc(n * sizeof(double));
    if (!values) return ENERGY_ENOMEM;
    for (int b = 0; b < ENERGY_BINS; b++) {
        counts[b] = 0;
        pos[b] = 0;
    }
    for (size_t i = 0; i < n; i++) {
        if (!isfinite(s[i].hr) || s[i].hr > thr) continue;
        counts[tod_bin(s[i].tod)]++;
    }
    offsets[0] = 0;
    for (int b = 0; b < ENERGY_BINS; b++) offsets[b + 1] = offsets[b] + counts[b];
    for (size_t i = 0; i < n; i++) {
        if (!isfinite(s[i].hr) || s[i].hr > thr) continue;
        const int b = tod_bin(s[i].tod);
        values[offsets[b] + pos[b]++] = s[i].hr;
    }
    int ngood = 0;
    for (int b = 0; b < ENERGY_BINS; b++) {
        prof[b] = (counts[b] >= 3) ? percentile_sorted(values + offsets[b], (size_t)counts[b], 50.0)
                                   : (double)NAN;
        if (isfinite(prof[b])) ngood++;
    }
    free(values);
    if (ngood == 0) {
        for (int b = 0; b < ENERGY_BINS; b++) prof[b] = 70.0;
        return ENERGY_OK;
    }

    /* circular interpolation */
    for (int b = 0; b < ENERGY_BINS; b++) {
        if (isfinite(prof[b])) continue;
        int prev = -1, next = -1;
        for (int j = 1; j <= ENERGY_BINS; j++) {
            const int q = ((b - j) % ENERGY_BINS + ENERGY_BINS) % ENERGY_BINS;
            if (isfinite(prof[q])) {
                prev = q;
                break;
            }
        }
        for (int j = 1; j <= ENERGY_BINS; j++) {
            const int q = (b + j) % ENERGY_BINS;
            if (isfinite(prof[q])) {
                next = q;
                break;
            }
        }
        if (prev < 0 || next < 0) {
            prof[b] = prof[prev >= 0 ? prev : next];
            continue;
        }
        const int d1 = (b - prev + ENERGY_BINS) % ENERGY_BINS;
        const int d2 = (next - b + ENERGY_BINS) % ENERGY_BINS;
        const double w = (d1 + d2) > 0 ? (double)d1 / (double)(d1 + d2) : 0.5;
        prof[b] = prof[prev] * (1.0 - w) + prof[next] * w;
    }
    /* circular 7-point smoothing */
    double smooth[ENERGY_BINS];
    for (int b = 0; b < ENERGY_BINS; b++) {
        double acc = 0.0;
        for (int j = -3; j <= 3; j++) {
            acc += prof[((b + j) % ENERGY_BINS + ENERGY_BINS) % ENERGY_BINS];
        }
        smooth[b] = acc / 7.0;
    }
    memcpy(prof, smooth, sizeof(smooth));
    return ENERGY_OK;
}

/* ------------------------------------------------------------------ */
/* 展示层                                                              */
/* ------------------------------------------------------------------ */
static double circ_raw(const energy_display_cfg_t *c, double tod)
{
    const double w = 2.0 * M_PI / 24.0;
    const double d = wrap24(tod - c->dip_hour + 12.0) - 12.0;
    const double dip = c->dip_depth * exp(-0.5 * (d / c->dip_width_h) * (d / c->dip_width_h));
    return cos(w * (tod - c->peak_hour)) - dip;
}

/* 归一化因子：在 288 点标准网格上取 max|c|。整天覆盖时与 Python 完全一致。 */
static double circ_scale(const energy_display_cfg_t *c)
{
    double mx = 1e-9;
    for (int i = 0; i < ENERGY_BINS; i++) {
        const double v = fabs(circ_raw(c, (double)i * 24.0 / ENERGY_BINS));
        if (v > mx) mx = v;
    }
    return mx;
}

double energy_circadian(const energy_display_cfg_t *c, double tod)
{
    const double scale = circ_scale(c);
    return scale > 1e-9 ? circ_raw(c, wrap24(tod)) / scale : 0.0;
}

double energy_display(const energy_display_cfg_t *c, double e, double tod)
{
    const double mod = clampd(1.0 + c->rho * energy_circadian(c, tod), 0.0, 2.0);
    return clampd(100.0 * clampd(e, 0.0, 1.0) * mod, 0.0, 100.0);
}

/* ------------------------------------------------------------------ */
/* Observation model                                                   */
/* ------------------------------------------------------------------ */
typedef struct {
    double m, lam, psi, sleep_prob;
} stage_mix_t;

static stage_mix_t stage_mix(const energy_params_t *p, const double stage[4])
{
    stage_mix_t out;
    out.sleep_prob = stage[ENERGY_LIGHT] + stage[ENERGY_DEEP] + stage[ENERGY_REM];
    out.m = out.sleep_prob;
    out.lam = stage[ENERGY_WAKE] * p->kappa_s[ENERGY_WAKE] +
              stage[ENERGY_LIGHT] * p->kappa_s[ENERGY_LIGHT] +
              stage[ENERGY_DEEP] * p->kappa_s[ENERGY_DEEP] +
              stage[ENERGY_REM] * p->kappa_s[ENERGY_REM];
    out.psi = stage[ENERGY_WAKE] * p->w_stage[ENERGY_WAKE] +
              stage[ENERGY_LIGHT] * p->w_stage[ENERGY_LIGHT] +
              stage[ENERGY_DEEP] * p->w_stage[ENERGY_DEEP] +
              stage[ENERGY_REM] * p->w_stage[ENERGY_REM];
    return out;
}

/* 返回 1=夜间（有分期记录），0=日间 */
static int is_night(const double stage[4])
{
    return (stage[ENERGY_WAKE] + stage[ENERGY_LIGHT] + stage[ENERGY_DEEP] + stage[ENERGY_REM]) > 0.5;
}

static double sleep_part(const double stage[4])
{
    return stage[ENERGY_LIGHT] + stage[ENERGY_DEEP] + stage[ENERGY_REM];
}

/* Observation model without its intercept: sum of slope x feature (E, stages / activity,
 * circadian harmonics). The intercept is the population prior coef[0] (returned in
 * *intercept) or, by default, the user's personal level. */
static double obs_regress(const energy_config_t *cfg, const energy_sample_t *s,
                          double e_state, int night, int channel,
                          double *slope_e, double *sigma, double *intercept)
{
    const energy_obs_t *o = &cfg->obs;
    double acc;
    const double *coef;
    if (night) {
        coef = (channel == 0) ? o->hr_night : o->hrv_night;
        acc = coef[1] * e_state + coef[2] * s->stage[ENERGY_DEEP] + coef[3] * s->stage[ENERGY_REM] +
              coef[4] * s->stage[ENERGY_WAKE];
        *sigma = (channel == 0) ? o->sigma_hr_night : o->sigma_hrv_night;
    } else {
        coef = (channel == 0) ? o->hr_day : o->hrv_day;
        const double w = 2.0 * M_PI / 24.0;
        const double t = wrap24(s->tod);
        acc = coef[1] * e_state + coef[2] * s->activity + coef[3] * cos(w * t) + coef[4] * sin(w * t) +
              coef[5] * cos(2.0 * w * t) + coef[6] * sin(2.0 * w * t);
        *sigma = (channel == 0) ? o->sigma_hr_day : o->sigma_hrv_day;
    }
    *slope_e = coef[1];
    *intercept = coef[0];
    return acc;
}

/* ------------------------------------------------------------------ */
/* EKF                                                                 */
/* ------------------------------------------------------------------ */
/* One EKF prediction step with explicit inputs (callers keep h <= ENERGY_EKF_MAX_STEP_H). */
static void ekf_predict_step(energy_ctx_t *ctx, const stage_mix_t *mix,
                             double load, double stress, double h)
{
    const energy_config_t *cfg = &ctx->cfg;
    double F[4], Phi[4];
    continuous_jacobian(&cfg->params, ctx->s, ctx->e, mix->m, mix->lam, mix->psi, load, stress, F);
    discrete_jacobian(F, h, Phi);

    double s_new, e_new;
    advance_raw(&cfg->params, ctx->s, ctx->e, mix->m, mix->lam, mix->psi, load, stress, h,
                &s_new, &e_new);

    /* P = Phi P Phi^T + Q h (diagonal Q); symmetric by construction */
    const double p00 = Phi[0] * ctx->p00 + Phi[1] * ctx->p01;
    const double p01 = Phi[0] * ctx->p01 + Phi[1] * ctx->p11;
    const double p10 = Phi[2] * ctx->p00 + Phi[3] * ctx->p01;
    const double p11 = Phi[2] * ctx->p01 + Phi[3] * ctx->p11;
    ctx->p00 = p00 * Phi[0] + p01 * Phi[1] + cfg->ekf.q_s * h;
    ctx->p01 = p00 * Phi[2] + p01 * Phi[3];
    ctx->p11 = p10 * Phi[2] + p11 * Phi[3] + cfg->ekf.q_e * h;
    if (!isfinite(ctx->p00) || !isfinite(ctx->p11) || ctx->p00 <= 0.0 || ctx->p11 <= 0.0) {
        /* covariance lost its meaning (extreme gap / observation): back to cold-start variance */
        ctx->p00 = cfg->ekf.p0_s;
        ctx->p01 = 0.0;
        ctx->p11 = cfg->ekf.p0_e;
    }
    ctx->s = clampd(s_new, 0.0, 1.0);
    ctx->e = clampd(e_new, 0.0, 1.0);
}

/* Cold-start R inflation for a channel with n committed episodes: x cold_r_inflate at n = 0,
 * linearly back to x1 at n = cold_nights. */
static double cold_r_scale(const energy_baseline_cfg_t *bc, int n)
{
    if (bc->cold_nights <= 0 || n >= bc->cold_nights) return 1.0;
    const double frac = (double)(bc->cold_nights - n) / (double)bc->cold_nights;
    return 1.0 + (bc->cold_r_inflate - 1.0) * frac;
}

static void ctx_update(energy_ctx_t *ctx, const energy_sample_t *s, int channel)
{
    const energy_config_t *cfg = &ctx->cfg;
    const energy_baseline_cfg_t *bc = &cfg->baseline;
    const double y = (channel == 0) ? s->hr : s->hrv_log;
    const double q = (channel == 0) ? finite_or(s->q_hr, 0.0) : finite_or(s->q_hrv, 0.0);
    if (!isfinite(y) || q <= 0.0) return;

    const int night = is_night(s->stage);
    const int k = (night ? ENERGY_BL_HR_NIGHT : ENERGY_BL_HR_DAY) + channel;
    double slope, sigma, intercept;
    const double reg = obs_regress(cfg, s, ctx->e, night, channel, &slope, &sigma, &intercept);
    if (!isfinite(reg)) return;

    /* Personal intercept ("your usual level"). A sample's implied level is evaluated on the
     * open-loop reference trajectory, which observations never touch: learning the baseline
     * therefore cannot feed back into E (using the posterior E here would let E and the
     * baseline drift together). Constant offsets of the whole HR series end up in the level,
     * so E only reflects deviations from the user's own norm. */
    double level = intercept;
    double r_scale = 1.0;
    if (bc->enable) {
        double slope_r, sigma_r, icpt_r;
        const double implied =
            y - obs_regress(cfg, s, ctx->ref_e, night, channel, &slope_r, &sigma_r, &icpt_r);
        int have = 0;
        if (ctx->bl.n[k] > 0) {
            level = ctx->bl.level[k];
            have = 1;
        } else if (ctx->acc_n[k] > 0 && ctx->acc_n[k] >= bc->bootstrap_samples) {
            /* first episode: running mean of the episode so far (bootstrap) */
            level = ctx->acc_sum[k] / (double)ctx->acc_n[k];
            have = 1;
        }
        const int learn = (q >= 0.5) && isfinite(implied) &&
                          (night ? sleep_part(s->stage) >= 0.5 : s->activity < bc->rest_activity);
        if (learn && ctx->acc_n[k] < INT_MAX) {
            /* once a level exists, clip artefacts so one spike cannot move the episode mean */
            const double lim = 6.0 * sigma;
            ctx->acc_sum[k] += (ctx->bl.n[k] > 0) ? clampd(implied, level - lim, level + lim) : implied;
            ctx->acc_n[k]++;
        }
        if (!have) return; /* no personal reference yet: "high" or "low" is undefined */
        r_scale = cold_r_scale(bc, ctx->bl.n[k]);
    }

    /* residual against the personal observation model (prior E): consumed by vitality.h */
    const double pred = level + reg;
    ctx->resid[channel] = y - pred;

    if (!cfg->day_update && !night) {
        /* Python: skip_update = (not day_update) and panel.is_day[k] */
        if (s->tod >= 7.0 && s->tod < 22.0) return;
    }

    /* H = [0, slope]: the observation depends on E only */
    const double R = (sigma * sigma) / fmax(q, cfg->ekf.q_floor) * r_scale;
    const double PHt0 = ctx->p01 * slope;
    const double PHt1 = ctx->p11 * slope;
    const double hph = slope * PHt1;
    const double S = hph + R;
    if (!isfinite(S) || S <= 1e-12) return;

    const double innov = y - pred;
    const double z = innov / sqrt(S);
    double R_use = R;
    if (cfg->ekf.robust_nu > 0.0) {
        R_use = R * (cfg->ekf.robust_nu + z * z) / (cfg->ekf.robust_nu + 1.0);
    }
    if (ctx->n_obs < ULONG_MAX) ctx->n_obs++;
    if (fabs(z) > cfg->ekf.gate_sigma) {
        if (cfg->ekf.gate_max_reject <= 0 || ctx->gate_rej[k] < cfg->ekf.gate_max_reject) {
            if (ctx->gate_rej[k] < INT_MAX) ctx->gate_rej[k]++;
            if (ctx->n_gated < ULONG_MAX) ctx->n_gated++;
            /* A rejection means model and data disagree: widen the E band instead of
             * reporting the model's own (false) confidence. */
            if (ctx->p11 < cfg->ekf.p0_e) ctx->p11 = fmin(ctx->p11 * 1.25, cfg->ekf.p0_e);
            return;
        }
        /* Soft acceptance after a run of rejections: inflate R until the normalised innovation
         * equals the gate. Influence stays bounded, but the filter can no longer lock out. */
        const double g = cfg->ekf.gate_sigma;
        const double r_soft = innov * innov / (g * g) - hph;
        if (r_soft > R_use) R_use = r_soft;
    } else {
        ctx->gate_rej[k] = 0;
    }
    const double S_use = hph + R_use;

    const double K0 = PHt0 / S_use;
    const double K1 = PHt1 / S_use;
    ctx->s += K0 * innov;
    ctx->e += K1 * innov;

    /* Joseph form P = (I-KH) P (I-KH)^T + K K^T R with H = [0, slope] */
    const double a00 = 1.0, a01 = -K0 * slope;
    const double a10 = 0.0, a11 = 1.0 - K1 * slope;
    const double t00 = a00 * ctx->p00 + a01 * ctx->p01;
    const double t01 = a00 * ctx->p01 + a01 * ctx->p11;
    const double t10 = a10 * ctx->p00 + a11 * ctx->p01;
    const double t11 = a10 * ctx->p01 + a11 * ctx->p11;
    ctx->p00 = t00 * a00 + t01 * a01 + K0 * K0 * R_use;
    ctx->p01 = t00 * a10 + t01 * a11 + K0 * K1 * R_use;
    ctx->p11 = t10 * a10 + t11 * a11 + K1 * K1 * R_use;
}

/* 开环推进：与 Python simulate 一样切成 <=2.5 分钟的子步 */
static void simulate_advance(const energy_params_t *p, double s, double e,
                             double m, double lam, double psi,
                             double load, double stress, double h,
                             double *s_out, double *e_out)
{
    int n_sub = (int)ceil(h / ENERGY_MAX_SUBSTEP_H - 1e-9);
    if (n_sub < 1) n_sub = 1;
    if (n_sub > ENERGY_MAX_SUBSTEPS) n_sub = ENERGY_MAX_SUBSTEPS;
    const double hs = h / (double)n_sub;
    for (int i = 0; i < n_sub; i++) {
        double s1, e1;
        advance_raw(p, s, e, m, lam, psi, load, stress, hs, &s1, &e1);
        s = clampd(s1, 0.0, 1.0);
        e = clampd(e1, 0.0, 1.0);
    }
    *s_out = s;
    *e_out = e;
}

static void ctx_clamp(energy_ctx_t *ctx)
{
    ctx->s = clampd(ctx->s, 0.0, 1.0);
    ctx->e = clampd(ctx->e, 0.0, 1.0);
}

/* ------------------------------------------------------------------ */
/* Streaming API                                                       */
/* ------------------------------------------------------------------ */
void energy_ctx_init(energy_ctx_t *ctx, const energy_config_t *cfg)
{
    if (!ctx) return;
    memset(ctx, 0, sizeof(*ctx));
    if (cfg) {
        ctx->cfg = *cfg;
    } else {
        energy_config_default(&ctx->cfg);
    }
    ctx->config_error = (energy_config_validate(&ctx->cfg, NULL, 0) != ENERGY_OK);
    if (ctx->config_error) {
        energy_config_default(&ctx->cfg);
    }
    ctx->s = ctx->cfg.x0_s;
    ctx->e = ctx->cfg.x0_e;
    ctx->ref_s = ctx->cfg.x0_s;
    ctx->ref_e = ctx->cfg.x0_e;
    ctx->p00 = ctx->cfg.ekf.p0_s;
    ctx->p01 = 0.0;
    ctx->p11 = ctx->cfg.ekf.p0_e;
    ctx->started = 0;
    ctx->bl = ctx->cfg.baseline.init; /* the profile's baseline, if any */
    ctx->stress_div = 1.0;
    ctx->lc_last_hr = (double)NAN;
    ctx->resid[0] = ctx->resid[1] = (double)NAN;
}

void energy_ctx_get_baseline(const energy_ctx_t *ctx, energy_baseline_t *out)
{
    if (!ctx || !out) return;
    *out = ctx->bl;
}

int energy_ctx_set_baseline(energy_ctx_t *ctx, const energy_baseline_t *bl)
{
    if (!ctx || !bl || !baseline_valid(bl)) return ENERGY_EINVAL;
    ctx->bl = *bl;
    return ENERGY_OK;
}

void energy_ctx_gate_stats(const energy_ctx_t *ctx, unsigned long *n_obs, unsigned long *n_gated)
{
    if (n_obs) *n_obs = ctx ? ctx->n_obs : 0ul;
    if (n_gated) *n_gated = ctx ? ctx->n_gated : 0ul;
}

/* Does the caller's raw sample carry any information? Rows without HR/HRV, stages, steps or a
 * positive load/activity (e.g. a no-wear placeholder) only advance time. */
static int sample_has_data(const energy_sample_t *s)
{
    if (isfinite(s->hr) && isfinite(s->q_hr) && s->q_hr > 0.0) return 1;
    if (isfinite(s->hrv_log) && isfinite(s->q_hrv) && s->q_hrv > 0.0) return 1;
    for (int i = 0; i < 4; i++) {
        if (isfinite(s->stage[i]) && s->stage[i] > 0.0) return 1;
    }
    if (isfinite(s->steps) && s->steps > 0.0) return 1;
    if (isfinite(s->load) && s->load > 0.0) return 1;
    if (isfinite(s->activity) && s->activity > 0.0) return 1;
    return 0;
}

/* ---- online load calibration (stream only, past samples only) ---------------------------- */

/* Stochastic-approximation quantile tracker: running mean for the first 12 samples, then
 * steps of +step*p / -step*(1-p), whose fixed point is the p-quantile of the input. */
static void quantile_track(double *q, int *n, double x, double p, double step)
{
    if (*n < 12) {
        *q = (*n == 0) ? x : *q + (x - *q) / (double)(*n + 1);
    } else {
        *q += (x >= *q) ? step * p : -step * (1.0 - p);
    }
    if (*n < INT_MAX) (*n)++;
}

static int hour_of(double tod)
{
    int h = (int)wrap24(finite_or(tod, 0.0));
    if (h < 0) h = 0;
    if (h > 23) h = 23;
    return h;
}

static int hr_base_given(const energy_load_cfg_t *c)
{
    for (int h = 0; h < 24; h++) {
        if (!(c->hr_base[h] > 0.0)) return 0;
    }
    return 1;
}

/* Same fallback chain as the batch estimate: deep-sleep median -> sleep P20 -> all P5 -> 60. */
static double lc_hr_rest(const energy_ctx_t *ctx)
{
    if (ctx->cfg.load.hr_rest > 0.0) return ctx->cfg.load.hr_rest;
    if (ctx->lc_n_deep >= 12) return ctx->lc_deep;
    if (ctx->lc_n_sleep >= 12) return ctx->lc_sleep;
    if (ctx->lc_n_all >= 1) return ctx->lc_all;
    return 60.0;
}

static double lc_hr_max(const energy_ctx_t *ctx)
{
    if (ctx->cfg.load.hr_max > 0.0) return ctx->cfg.load.hr_max;
    return (ctx->lc_n_hr < 20) ? 185.0 : clampd(ctx->lc_peak, 150.0, 210.0);
}

static double lc_hr_base(const energy_ctx_t *ctx, double tod)
{
    const int h = hour_of(tod);
    if (hr_base_given(&ctx->cfg.load)) return ctx->cfg.load.hr_base[h];
    if (ctx->lc_n_base[h] >= 3) return ctx->lc_base[h];
    return lc_hr_rest(ctx) + 10.0;
}

/* Online P99.5 of the raw load, blended with a one-day prior so a fresh stream starts where a
 * typical batch window ends up (the v1.1 constant 0.25 saturated L at rest). */
static double lc_l_ref(const energy_ctx_t *ctx)
{
    if (ctx->cfg.load.l_ref > 0.0) return ctx->cfg.load.l_ref;
    const double w = ctx->lc_hist_w;
    double data = ENERGY_LREF_PRIOR;
    if (w > 0.0) {
        const double target = 0.995 * w;
        double cum = 0.0;
        data = ENERGY_LREF_BINS * ENERGY_LREF_BIN_W;
        for (int b = 0; b < ENERGY_LREF_BINS; b++) {
            const double hb = ctx->lc_hist[b];
            if (hb > 0.0 && cum + hb >= target) {
                data = ((double)b + (target - cum) / hb) * ENERGY_LREF_BIN_W;
                break;
            }
            cum += hb;
        }
    }
    const double l_ref = (ENERGY_LREF_PRIOR_W * ENERGY_LREF_PRIOR + w * data) / (ENERGY_LREF_PRIOR_W + w);
    return l_ref > 0.25 ? l_ref : 0.25;
}

static double sample_load_raw(const energy_ctx_t *ctx, const energy_sample_t *s, double dt_h)
{
    return energy_load_raw(&ctx->cfg.load, s->hr, lc_hr_rest(ctx), lc_hr_max(ctx),
                           lc_hr_base(ctx, s->tod), s->steps, dt_h);
}

/* Feed one data sample into the online calibration (after its own load was computed). */
static void lc_observe(energy_ctx_t *ctx, const energy_sample_t *s, double raw, double h_since)
{
    ctx->lc_peak = fmax(ctx->lc_peak - ENERGY_PEAK_DECAY_PER_H * fmax(h_since, 0.0), 0.0);
    if (isfinite(s->hr) && s->q_hr > 0.0) {
        const double hr = s->hr;
        if (s->stage[ENERGY_DEEP] > 0.5) quantile_track(&ctx->lc_deep, &ctx->lc_n_deep, hr, 0.50, 1.0);
        if (sleep_part(s->stage) > 0.5) quantile_track(&ctx->lc_sleep, &ctx->lc_n_sleep, hr, 0.20, 1.0);
        quantile_track(&ctx->lc_all, &ctx->lc_n_all, hr, 0.05, 1.0);
        const int hb = hour_of(s->tod);
        quantile_track(&ctx->lc_base[hb], &ctx->lc_n_base[hb], hr, 0.45, 1.0);
        /* peak of min(two adjacent samples): a single artefact spike cannot set hr_max */
        if (isfinite(ctx->lc_last_hr)) {
            const double robust = fmin(hr, ctx->lc_last_hr);
            if (robust > ctx->lc_peak) ctx->lc_peak = robust;
        }
        ctx->lc_last_hr = hr;
        if (ctx->lc_n_hr < INT_MAX) ctx->lc_n_hr++;
    } else {
        ctx->lc_last_hr = (double)NAN;
    }
    if (isfinite(raw) && raw >= 0.0) {
        int b = (int)fmin(raw / ENERGY_LREF_BIN_W, (double)(ENERGY_LREF_BINS - 1));
        if (b < 0) b = 0;
        ctx->lc_hist[b] += 1.0;
        ctx->lc_hist_w += 1.0;
        if (ctx->lc_hist_w > ENERGY_LREF_MEMORY) {
            for (int i = 0; i < ENERGY_LREF_BINS; i++) ctx->lc_hist[i] *= 0.5;
            ctx->lc_hist_w *= 0.5;
        }
    }
}

static double sample_activity(const energy_config_t *cfg, const energy_sample_t *s)
{
    if (isfinite(s->activity) && s->activity >= 0.0) return clampd(s->activity, 0.0, 1.0);
    return energy_activity_step(s->steps, s->kcal, cfg->load.steps_ref);
}

/* ---- personal baseline episodes ----------------------------------------------------------- */

static double baseline_rate(const energy_baseline_cfg_t *bc, int n)
{
    if (n <= 0) return 1.0;
    return (n < bc->cold_nights) ? 1.0 / (double)(n + 1) : bc->lambda;
}

/* Commit the finished episode's mean implied level of channel k into the baseline. */
static void baseline_commit(energy_ctx_t *ctx, int k)
{
    const energy_baseline_cfg_t *bc = &ctx->cfg.baseline;
    if (ctx->acc_n[k] >= bc->min_samples) {
        const double mean = ctx->acc_sum[k] / (double)ctx->acc_n[k];
        const int n = ctx->bl.n[k];
        double lam = baseline_rate(bc, n);
        if (n > 0 && n >= bc->cold_nights) {
            const double dev = mean - ctx->bl.level[k];
            if (fabs(dev) > bc->anomaly[k]) {
                /* Anomalous episode (illness, alcohol, travel): down-weight it. Three anomalies
                 * in a row in the same direction are a real change of the norm (new device,
                 * long-term training) and are learned at the normal rate again. */
                const int dir = dev > 0.0 ? 1 : -1;
                int run = (ctx->anomaly_run[k] * dir > 0) ? ctx->anomaly_run[k] + dir : dir;
                if (run > 3) run = 3;
                if (run < -3) run = -3;
                ctx->anomaly_run[k] = run;
                if (abs(run) < 3) lam *= bc->anomaly_weight;
            } else {
                ctx->anomaly_run[k] = 0;
            }
        }
        ctx->bl.level[k] = (n > 0) ? ctx->bl.level[k] + lam * (mean - ctx->bl.level[k]) : mean;
        if (ctx->bl.n[k] < INT_MAX) ctx->bl.n[k]++;
    }
    ctx->acc_sum[k] = 0.0;
    ctx->acc_n[k] = 0;
}

static double circ_blend(double cur, double obs, double lam)
{
    const double d = wrap24(obs - cur + 12.0) - 12.0; /* shortest signed distance on the clock */
    return wrap24(cur + lam * d);
}

static void episode_close_night(energy_ctx_t *ctx)
{
    if (ctx->night_len >= ENERGY_MIN_NIGHT_SAMPLES) {
        baseline_commit(ctx, ENERGY_BL_HR_NIGHT);
        baseline_commit(ctx, ENERGY_BL_HRV_NIGHT);
        /* sleep habits for the gap policy: usual window and stage composition */
        const double lam = baseline_rate(&ctx->cfg.baseline, ctx->bl.n_sleep);
        const double inv = 1.0 / (double)ctx->night_len;
        if (ctx->bl.n_sleep == 0) {
            ctx->bl.sleep_on = wrap24(ctx->night_on);
            ctx->bl.sleep_off = wrap24(ctx->night_off);
            for (int i = 0; i < 4; i++) ctx->bl.stage_mix[i] = clampd(ctx->night_stage[i] * inv, 0.0, 1.0);
        } else {
            ctx->bl.sleep_on = circ_blend(ctx->bl.sleep_on, ctx->night_on, lam);
            ctx->bl.sleep_off = circ_blend(ctx->bl.sleep_off, ctx->night_off, lam);
            for (int i = 0; i < 4; i++) {
                const double v = clampd(ctx->night_stage[i] * inv, 0.0, 1.0);
                ctx->bl.stage_mix[i] += lam * (v - ctx->bl.stage_mix[i]);
            }
        }
        if (ctx->bl.n_sleep < INT_MAX) ctx->bl.n_sleep++;
    }
    /* too short (nap, stray stage rows): discarded */
    ctx->acc_sum[ENERGY_BL_HR_NIGHT] = ctx->acc_sum[ENERGY_BL_HRV_NIGHT] = 0.0;
    ctx->acc_n[ENERGY_BL_HR_NIGHT] = ctx->acc_n[ENERGY_BL_HRV_NIGHT] = 0;
    ctx->night_len = 0;
    for (int i = 0; i < 4; i++) ctx->night_stage[i] = 0.0;
}

static void episode_close_day(energy_ctx_t *ctx)
{
    if (ctx->day_len >= ENERGY_MIN_DAY_SAMPLES) {
        baseline_commit(ctx, ENERGY_BL_HR_DAY);
        baseline_commit(ctx, ENERGY_BL_HRV_DAY);
    }
    ctx->acc_sum[ENERGY_BL_HR_DAY] = ctx->acc_sum[ENERGY_BL_HRV_DAY] = 0.0;
    ctx->acc_n[ENERGY_BL_HR_DAY] = ctx->acc_n[ENERGY_BL_HRV_DAY] = 0;
    ctx->day_len = 0;
}

/* Commit whatever episodes are still open (end of a batch, or before a long data gap). */
static void baseline_flush(energy_ctx_t *ctx)
{
    episode_close_night(ctx);
    episode_close_day(ctx);
    ctx->run_night = 0;
    ctx->run_day = 0;
}

/* Episode bookkeeping for one data sample. Night = sample with a stage record; an episode
 * ends after ENERGY_EPISODE_SWITCH consecutive samples of the other regime (hysteresis keeps
 * brief awakenings without stages inside the night), or after ENERGY_EPISODE_MAX samples. */
static void episode_track(energy_ctx_t *ctx, const energy_sample_t *s)
{
    if (is_night(s->stage)) {
        if (ctx->run_night < INT_MAX) ctx->run_night++;
        ctx->run_day = 0;
        if (ctx->night_len == 0) ctx->night_on = s->tod;
        ctx->night_off = s->tod;
        ctx->night_len++;
        for (int i = 0; i < 4; i++) ctx->night_stage[i] += s->stage[i];
        if (ctx->run_night >= ENERGY_EPISODE_SWITCH && ctx->day_len > 0) episode_close_day(ctx);
        if (ctx->night_len >= ENERGY_EPISODE_MAX) episode_close_night(ctx);
    } else {
        if (ctx->run_day < INT_MAX) ctx->run_day++;
        ctx->run_night = 0;
        ctx->day_len++;
        if (ctx->run_day >= ENERGY_EPISODE_SWITCH && ctx->night_len > 0) episode_close_night(ctx);
        if (ctx->day_len >= ENERGY_EPISODE_MAX) episode_close_day(ctx);
    }
}

/* ---- time advance and gap policy ---------------------------------------------------------- */

/* Advance state, covariance and the open-loop reference over h hours with constant inputs.
 * The EKF predicts in sub-steps of <= ENERGY_EKF_MAX_STEP_H so the Jacobian is not frozen
 * over a long gap; the reference uses exactly the same discretisation, hence E == E_ref
 * whenever no observation has been applied. */
static void ctx_advance_const(energy_ctx_t *ctx, const double stage[4], double load,
                              double stress, double h)
{
    if (!(h > 0.0)) return;
    const energy_params_t *p = &ctx->cfg.params;
    const stage_mix_t mix = stage_mix(p, stage);
    if (ctx->cfg.use_ekf) {
        int n_sub = (int)fmin(ceil(h / ENERGY_EKF_MAX_STEP_H - 1e-9), (double)ENERGY_MAX_SUBSTEPS);
        if (n_sub < 1) n_sub = 1;
        const double hs = h / (double)n_sub;
        for (int i = 0; i < n_sub; i++) {
            double s1, e1;
            advance_raw(p, ctx->ref_s, ctx->ref_e, mix.m, mix.lam, mix.psi, load, stress, hs, &s1, &e1);
            ctx->ref_s = clampd(s1, 0.0, 1.0);
            ctx->ref_e = clampd(e1, 0.0, 1.0);
            ekf_predict_step(ctx, &mix, load, stress, hs);
        }
    } else {
        double s1, e1;
        simulate_advance(p, ctx->s, ctx->e, mix.m, mix.lam, mix.psi, load, stress, h, &s1, &e1);
        ctx->s = s1;
        ctx->e = e1;
        simulate_advance(p, ctx->ref_s, ctx->ref_e, mix.m, mix.lam, mix.psi, load, stress, h, &s1, &e1);
        ctx->ref_s = s1;
        ctx->ref_e = e1;
    }
}

static int in_window(double t, double on, double off)
{
    t = wrap24(t);
    if (on <= off) return t >= on && t < off;
    return t >= on || t < off;
}

static double hours_to_boundary(double t, double on, double off)
{
    double d_on = wrap24(on - t), d_off = wrap24(off - t);
    if (!(d_on > 1e-9)) d_on = 24.0;
    if (!(d_off > 1e-9)) d_off = 24.0;
    return fmin(d_on, d_off);
}

/* No data for longer than gap.max_hold_h: stop holding the last inputs (v1.1 propagated the
 * last daytime sample through a whole unworn night). Inside the user's usual sleep window
 * (learned, else cfg.gap) advance as soft sleep with the user's typical stage mix (else the
 * template); outside it as quiet wakefulness. L = 0 and stress = 0 in both cases. */
static unsigned ctx_advance_gap(energy_ctx_t *ctx, double tod, double h)
{
    static const double awake[4] = {0.0, 0.0, 0.0, 0.0};
    const energy_gap_cfg_t *g = &ctx->cfg.gap;
    double on = g->sleep_start, off = g->sleep_end;
    const double *tmpl = g->template_stage;
    if (ctx->bl.n_sleep > 0) {
        on = ctx->bl.sleep_on;
        off = ctx->bl.sleep_off;
        tmpl = ctx->bl.stage_mix;
    }
    unsigned flags = ENERGY_FLAG_GAP;
    double t = wrap24(tod);
    for (int seg = 0; h > 1e-12 && seg < ENERGY_MAX_GAP_SEGMENTS; seg++) {
        const int sleeping = in_window(t, on, off);
        double d = hours_to_boundary(t, on, off);
        if (d > h || seg == ENERGY_MAX_GAP_SEGMENTS - 1) d = h;
        ctx_advance_const(ctx, sleeping ? tmpl : awake, 0.0, 0.0, d);
        if (sleeping) flags |= ENERGY_FLAG_GAP_SLEEP;
        t = wrap24(t + d);
        h -= d;
    }
    return flags;
}

/* Advance from the previous sample to t = ctx->t_h + h. Returns ENERGY_FLAG_* for the gap. */
static unsigned ctx_advance(energy_ctx_t *ctx, double h)
{
    const double hold_left = ctx->t_last_data + ctx->cfg.gap.max_hold_h - ctx->t_h;
    const double h_hold = clampd(hold_left, 0.0, h);
    ctx_advance_const(ctx, ctx->prev_stage, ctx->prev_load, ctx->prev_stress, h_hold);
    const double h_gap = h - h_hold;
    if (!(h_gap > 1e-12)) return 0u;
    return ctx_advance_gap(ctx, ctx->tod + h_hold, h_gap);
}

/* ---- per-sample step ---------------------------------------------------------------------- */

static void point_band(const energy_display_cfg_t *d, double e, double sd, double tod,
                       energy_point_t *out)
{
    const double mod = clampd(1.0 + d->rho * energy_circadian(d, tod), 0.0, 2.0);
    const double ec = clampd(e, 0.0, 1.0);
    out->base = 100.0 * ec;
    out->energy = clampd(100.0 * ec * mod, 0.0, 100.0);
    out->lower = clampd(100.0 * clampd(ec - d->band_z * sd, 0.0, 1.0) * mod, 0.0, 100.0);
    out->upper = clampd(100.0 * clampd(ec + d->band_z * sd, 0.0, 1.0) * mod, 0.0, 100.0);
}

static void ctx_emit(const energy_ctx_t *ctx, const energy_sample_t *s, double load,
                     double activity, unsigned flags, energy_point_t *out)
{
    if (!out) return;
    const stage_mix_t mix = stage_mix(&ctx->cfg.params, s->stage);
    memset(out, 0, sizeof(*out));
    out->t_h = s->t_h;
    out->tod = s->tod;
    out->s = ctx->s;
    out->e = ctx->e;
    out->e_sd = ctx->cfg.use_ekf ? sqrt(fmax(ctx->p11, 0.0)) : 0.0;
    point_band(&ctx->cfg.display, ctx->e, out->e_sd, s->tod, out);
    out->load = load;
    out->activity = activity;
    out->sleep_prob = mix.sleep_prob;
    out->flags = flags;
    out->hr_resid = ctx->resid[0];
    out->hrv_resid = ctx->resid[1];
}

/* Field clean-up. Stress must already be scaled to 0-1 by the caller (energy_ctx_step /
 * energy_compute decide the scale per stream / per series, never per point). */
static void sanitize_sample(energy_sample_t *s)
{
    for (int i = 0; i < 4; i++) {
        if (!isfinite(s->stage[i]) || s->stage[i] < 0.0) s->stage[i] = 0.0;
        s->stage[i] = clampd(s->stage[i], 0.0, 1.0);
    }
    const double sum = s->stage[0] + s->stage[1] + s->stage[2] + s->stage[3];
    if (sum > 1.0) {
        for (int i = 0; i < 4; i++) s->stage[i] /= sum;
    }
    if (!isfinite(s->tod) || s->tod < 0.0) s->tod = wrap24(s->t_h);
    if (!isfinite(s->tod)) s->tod = 0.0;
    s->tod = wrap24(s->tod);
    s->stress = clampd(finite_or(s->stress, 0.0), 0.0, 1.0);
    if (!isfinite(s->steps) || s->steps < 0.0) s->steps = 0.0;
    if (!isfinite(s->kcal) || s->kcal < 0.0) s->kcal = 0.0;
    s->q_hr = clampd(finite_or(s->q_hr, 0.0), 0.0, 1.0);
    s->q_hrv = clampd(finite_or(s->q_hrv, 0.0), 0.0, 1.0);
}

/* Shared by energy_ctx_step and energy_compute. has_data is decided on the caller's raw
 * sample (batch interpolation must not turn a no-wear row into data). */
static int ctx_step_internal(energy_ctx_t *ctx, const energy_sample_t *s_in, energy_point_t *out,
                             int has_data)
{
    if (ctx->config_error) return ENERGY_EINVAL;
    energy_sample_t local = *s_in;
    sanitize_sample(&local);
    if (!isfinite(local.t_h) || fabs(local.t_h) > ENERGY_MAX_HOURS) return ENERGY_EINVAL;

    /* interval covered by this sample: 5 min for the first one, the real spacing afterwards */
    double h = ENERGY_HOURS_PER_STEP;
    double h_since_data = 0.0;
    unsigned flags = 0u;
    if (!ctx->started) {
        ctx->started = 1;
        ctx->t_last_data = local.t_h;
    } else {
        h = local.t_h - ctx->t_h;
        if (!(h > 0.0)) return ENERGY_EDOMAIN;
        h_since_data = local.t_h - ctx->t_last_data;
        flags = ctx_advance(ctx, h);
        if (has_data && h_since_data > ENERGY_EPISODE_BREAK_H) baseline_flush(ctx);
    }
    ctx->t_h = local.t_h;
    ctx->tod = local.tod;
    ctx->resid[0] = ctx->resid[1] = (double)NAN;

    double load = 0.0, activity = 0.0;
    if (has_data) {
        /* TRIMP coverage is clamped exactly like the batch path (v1.1 used the raw interval,
         * so the first sample after a 3 h gap saturated L) */
        const double dt = clampd(h, 0.0, ENERGY_TRIMP_MAX_DT_H);
        const double raw = sample_load_raw(ctx, &local, dt);
        load = (isfinite(local.load) && local.load >= 0.0) ? clampd(local.load, 0.0, 1.0)
                                                           : clampd(raw / lc_l_ref(ctx), 0.0, 1.0);
        activity = sample_activity(&ctx->cfg, &local);
        lc_observe(ctx, &local, raw, h_since_data);
        local.load = load;
        local.activity = activity;

        ctx->prev_load = load;
        ctx->prev_stress = local.stress;
        memcpy(ctx->prev_stage, local.stage, sizeof(ctx->prev_stage));
        ctx->t_last_data = local.t_h;

        episode_track(ctx, &local);
        if (ctx->cfg.use_ekf) {
            ctx_update(ctx, &local, 0);
            ctx_update(ctx, &local, 1);
            ctx_clamp(ctx);
        }
    }
    if (ctx->cfg.baseline.enable && ctx->bl.n[ENERGY_BL_HR_NIGHT] < ctx->cfg.baseline.cold_nights) {
        flags |= ENERGY_FLAG_LEARNING;
    }
    ctx_emit(ctx, &local, load, activity, flags, out);
    return ENERGY_OK;
}

int energy_ctx_step(energy_ctx_t *ctx, const energy_sample_t *sample, energy_point_t *out)
{
    if (!ctx || !sample) return ENERGY_EINVAL;
    if (ctx->config_error) return ENERGY_EINVAL;
    energy_sample_t local = *sample;
    double div = ctx->stress_div;
    if (ctx->cfg.stress_scale > 0.0) {
        div = ctx->cfg.stress_scale;
    } else if (isfinite(local.stress) && local.stress > ENERGY_STRESS_AUTO_MAX) {
        div = 100.0; /* the stream turned out to be 0-100: switch for good */
    }
    if (isfinite(local.stress)) local.stress /= div;
    const int rc = ctx_step_internal(ctx, &local, out, sample_has_data(sample));
    if (rc == ENERGY_OK && !(ctx->cfg.stress_scale > 0.0)) ctx->stress_div = div;
    return rc;
}

void energy_ctx_peek(const energy_ctx_t *ctx, double tod, energy_point_t *out)
{
    if (!ctx || !out) return;
    if (!isfinite(tod) || tod < 0.0) {
        tod = ctx->started ? ctx->tod : (isfinite(ctx->t_h) ? wrap24(ctx->t_h) : 0.0);
    }
    memset(out, 0, sizeof(*out));
    out->t_h = ctx->t_h;
    out->tod = tod;
    out->s = ctx->s;
    out->e = ctx->e;
    out->e_sd = ctx->cfg.use_ekf ? sqrt(fmax(ctx->p11, 0.0)) : 0.0;
    out->hr_resid = out->hrv_resid = (double)NAN;
    point_band(&ctx->cfg.display, ctx->e, out->e_sd, tod, out);
    if (ctx->cfg.baseline.enable && ctx->bl.n[ENERGY_BL_HR_NIGHT] < ctx->cfg.baseline.cold_nights) {
        out->flags = ENERGY_FLAG_LEARNING;
    }
}

/* ------------------------------------------------------------------ */
/* Batch                                                               */
/* ------------------------------------------------------------------ */
void energy_output_free(energy_output_t *out)
{
    if (!out) return;
    free(out->points);
    free(out->daily);
    out->points = NULL;
    out->daily = NULL;
    out->n_points = 0;
    out->n_daily = 0;
}

/* Raw load series + calibration (as Python load.py). Every allocation failure is reported as
 * ENERGY_ENOMEM. prof_out (ENERGY_BINS values) may be NULL. */
static int load_raw_series(const energy_sample_t *samples, size_t n,
                           const energy_load_cfg_t *cfg, double *raw,
                           double *hr_rest_out, double *hr_max_out, double *l_ref_out,
                           double *prof_out)
{
    int rc;
    double hr_rest = cfg->hr_rest, hr_max = cfg->hr_max;
    if (!(hr_rest > 0.0) && (rc = estimate_hr_rest(samples, n, &hr_rest)) != ENERGY_OK) return rc;
    if (!(hr_max > 0.0) && (rc = estimate_hr_max(samples, n, &hr_max)) != ENERGY_OK) return rc;
    double prof[ENERGY_BINS];
    if (hr_base_given(cfg)) {
        /* profile value: hourly baseline expanded to the 5-minute grid */
        for (int b = 0; b < ENERGY_BINS; b++) prof[b] = cfg->hr_base[b * 24 / ENERGY_BINS];
    } else if ((rc = estimate_hr_base(samples, n, prof)) != ENERGY_OK) {
        return rc;
    }

    for (size_t i = 0; i < n; i++) {
        const double dt = (i == 0) ? ENERGY_HOURS_PER_STEP
                                   : clampd(samples[i].t_h - samples[i - 1].t_h, 0.0,
                                            ENERGY_TRIMP_MAX_DT_H);
        raw[i] = energy_load_raw(cfg, samples[i].hr, hr_rest, hr_max, prof[tod_bin(samples[i].tod)],
                                 samples[i].steps, dt);
    }
    double l_ref = cfg->l_ref;
    if (!(l_ref > 0.0)) {
        double *tmp = (double *)malloc(n * sizeof(double));
        if (!tmp) return ENERGY_ENOMEM;
        memcpy(tmp, raw, n * sizeof(double));
        l_ref = percentile_sorted(tmp, n, 99.5);
        free(tmp);
        if (!(l_ref > 0.25)) l_ref = 0.25;
    }
    if (hr_rest_out) *hr_rest_out = hr_rest;
    if (hr_max_out) *hr_max_out = hr_max;
    if (l_ref_out) *l_ref_out = l_ref;
    if (prof_out) memcpy(prof_out, prof, sizeof(prof));
    return ENERGY_OK;
}

static int compute_loads(const energy_sample_t *samples, size_t n,
                         const energy_load_cfg_t *cfg,
                         double *load, double *activity)
{
    double *raw = (double *)malloc(n * sizeof(double));
    if (!raw) return ENERGY_ENOMEM;
    double l_ref = 0.0;
    const int rc = load_raw_series(samples, n, cfg, raw, NULL, NULL, &l_ref, NULL);
    if (rc != ENERGY_OK) {
        free(raw);
        return rc;
    }
    for (size_t i = 0; i < n; i++) {
        /* caller-provided values (>= 0) win, otherwise the library's normalised load */
        const double given_load = samples[i].load;
        load[i] = (isfinite(given_load) && given_load >= 0.0)
                      ? clampd(given_load, 0.0, 1.0)
                      : clampd(raw[i] / l_ref, 0.0, 1.0);
        const double given_act = samples[i].activity;
        activity[i] = (isfinite(given_act) && given_act >= 0.0)
                          ? clampd(given_act, 0.0, 1.0)
                          : energy_activity_step(samples[i].steps, samples[i].kcal, cfg->steps_ref);
    }
    free(raw);
    return ENERGY_OK;
}

int energy_load_calibrate(const energy_sample_t *samples, size_t n,
                          const energy_load_cfg_t *cfg, energy_load_stats_t *out)
{
    if (!samples || n == 0 || n > ENERGY_MAX_SAMPLES || !out || !cfg) return ENERGY_EINVAL;
    memset(out, 0, sizeof(*out));
    /* same clean-up as energy_compute (tod inference, stage normalisation, NaN steps) */
    energy_sample_t *work = (energy_sample_t *)malloc(n * sizeof(energy_sample_t));
    double *raw = (double *)malloc(n * sizeof(double));
    int rc = ENERGY_ENOMEM;
    if (work && raw) {
        memcpy(work, samples, n * sizeof(energy_sample_t));
        for (size_t i = 0; i < n; i++) sanitize_sample(&work[i]);
        double prof[ENERGY_BINS];
        rc = load_raw_series(work, n, cfg, raw, &out->hr_rest, &out->hr_max, &out->l_ref, prof);
        if (rc == ENERGY_OK) {
            const int per_hour = ENERGY_BINS / 24;
            for (int h = 0; h < 24; h++) {
                double acc = 0.0;
                for (int j = 0; j < per_hour; j++) acc += prof[h * per_hour + j];
                out->hr_base[h] = acc / (double)per_hour;
            }
        }
    }
    free(work);
    free(raw);
    return rc;
}

static int days_in_month(int y, int m)
{
    static const int d[12] = {31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31};
    if (m == 2) {
        const int leap = (y % 4 == 0 && y % 100 != 0) || (y % 400 == 0);
        return leap ? 29 : 28;
    }
    return (m >= 1 && m <= 12) ? d[m - 1] : 30;
}

static void date_label(const energy_config_t *cfg, int day, char out[11])
{
    out[0] = '\0';
    if (cfg->start_year <= 0) return;
    if (cfg->start_month < 1 || cfg->start_month > 12) return;
    if (cfg->start_day < 1 || cfg->start_day > 31) return;
    int y = cfg->start_year, m = cfg->start_month, d = cfg->start_day + day;
    while (d > days_in_month(y, m)) {
        d -= days_in_month(y, m);
        m++;
        if (m > 12) {
            m = 1;
            y++;
        }
    }
    snprintf(out, 11, "%02d-%02d", m % 100, d % 100); /* bounded: no truncation */
}

static void smooth_display(const double *e_in, size_t n, double minutes, double *e_out)
{
    memcpy(e_out, e_in, n * sizeof(double));
    if (!isfinite(minutes) || minutes <= 0.0) return;
    if (minutes > ENERGY_MAX_SMOOTH_MIN) minutes = ENERGY_MAX_SMOOTH_MIN;
    int k = (int)floor(minutes / 5.0 + 0.5);
    if (k <= 1 || (size_t)k >= n) return;
    const int pad = k / 2;
    for (size_t i = 0; i < n; i++) {
        double acc = 0.0;
        for (int j = 0; j < k; j++) {
            long idx = (long)i - pad + j;
            if (idx < 0) idx = 0;
            if (idx >= (long)n) idx = (long)n - 1;
            acc += e_in[idx];
        }
        e_out[i] = acc / (double)k;
    }
}

/* Interpolated (non-observation) fields of a sample, see energy_compute. */
static double *interp_field(energy_sample_t *s, int f)
{
    switch (f) {
    case 0: return &s->stress;
    case 1: return &s->steps;
    default: return &s->kcal;
    }
}

/* Daily summary over points[i0, i1) (one local calendar day). */
static void daily_summary(const energy_config_t *cfg, const energy_point_t *points,
                          size_t i0, size_t i1, int day, energy_daily_t *d)
{
    const size_t cnt = i1 - i0;
    d->day = day;
    d->n = (cnt > (size_t)INT_MAX) ? INT_MAX : (int)cnt;
    date_label(cfg, day, d->date);

    double dmin = 1e9, dmax = -1e9, dsum = 0.0, debt = 0.0;
    size_t morning = i0, evening = i1 - 1, peak = i0, dip = i0;
    int have_peak = 0, have_dip = 0;
    double peak_v = -1e9, dip_v = 1e9;
    for (size_t k = i0; k < i1; k++) {
        const energy_point_t *p = &points[k];
        if (p->energy < dmin) dmin = p->energy;
        if (p->energy > dmax) dmax = p->energy;
        dsum += p->energy;
        if (p->s > debt) debt = p->s;
        if (p->tod < 9.0) morning = k;  /* last point before 09:00 */
        if (p->tod > 20.0) evening = k; /* last point after 20:00 */
        if (p->tod >= 6.0 && p->tod <= 16.0 && p->energy > peak_v) {
            peak_v = p->energy;
            peak = k;
            have_peak = 1;
        }
        if (p->tod >= 13.0 && p->tod <= 17.0 && p->energy < dip_v) {
            dip_v = p->energy;
            dip = k;
            have_dip = 1;
        }
    }
    if (!have_peak || !have_dip) {
        /* window not covered: fall back to the whole day */
        double mx = -1e9, mn = 1e9;
        for (size_t k = i0; k < i1; k++) {
            if (!have_peak && points[k].energy > mx) {
                mx = points[k].energy;
                peak = k;
            }
            if (!have_dip && points[k].energy < mn) {
                mn = points[k].energy;
                dip = k;
            }
        }
    }
    d->morning = points[morning].energy;
    d->evening = points[evening].energy;
    d->mean = dsum / (double)cnt;
    d->min = dmin;
    d->max = dmax;
    d->peak_tod = points[peak].tod;
    d->dip_tod = points[dip].tod;
    d->sleep_debt = debt;
}

int energy_compute(const energy_sample_t *samples, size_t n,
                   const energy_config_t *cfg_in, energy_output_t *out)
{
    if (!samples || n < 2 || n > ENERGY_MAX_SAMPLES || !out) return ENERGY_EINVAL;
    memset(out, 0, sizeof(*out));

    energy_config_t cfg;
    if (cfg_in) {
        cfg = *cfg_in;
    } else {
        energy_config_default(&cfg);
    }
    if (energy_config_validate(&cfg, NULL, 0) != ENERGY_OK) return ENERGY_EINVAL;

    /* time axis: finite, strictly increasing, sane magnitude */
    for (size_t i = 0; i < n; i++) {
        const double t = samples[i].t_h;
        if (!isfinite(t) || fabs(t) > ENERGY_MAX_HOURS) return ENERGY_EINVAL;
        if (i > 0 && t <= samples[i - 1].t_h) return ENERGY_EDOMAIN;
    }

    /* Single exit: every buffer is released at `cleanup`; on success points/daily are handed
     * over to `out` and the local pointers cleared (v1.1 freed e_raw/e_smooth twice on OOM). */
    int rc = ENERGY_OK;
    energy_sample_t *work = (energy_sample_t *)malloc(n * sizeof(energy_sample_t));
    double *scratch = (double *)malloc(n * sizeof(double));
    double *load = (double *)malloc(n * sizeof(double));
    double *activity = (double *)malloc(n * sizeof(double));
    double *e_smooth = (double *)malloc(n * sizeof(double));
    energy_point_t *points = (energy_point_t *)calloc(n, sizeof(energy_point_t));
    energy_daily_t *daily = (energy_daily_t *)calloc(n / 12 + 2, sizeof(energy_daily_t));
    if (!work || !scratch || !load || !activity || !e_smooth || !points || !daily) {
        rc = ENERGY_ENOMEM;
        goto cleanup;
    }
    memcpy(work, samples, n * sizeof(energy_sample_t));

    /* Stress scale is decided once for the whole series (never per point: a 1.2 on a 0-100
     * device is low stress, not full scale). */
    {
        double div = cfg.stress_scale;
        if (!(div > 0.0)) {
            double mx = -1.0;
            for (size_t i = 0; i < n; i++) {
                if (isfinite(work[i].stress) && work[i].stress > mx) mx = work[i].stress;
            }
            div = (mx > ENERGY_STRESS_AUTO_MAX) ? 100.0 : 1.0;
        }
        for (size_t i = 0; i < n; i++) work[i].stress /= div; /* NaN stays NaN */
    }
    /* linear gap filling of stress / steps / kcal (as the Python API) */
    for (int f = 0; f < 3; f++) {
        for (size_t i = 0; i < n; i++) scratch[i] = *interp_field(&work[i], f);
        fill_linear(scratch, n);
        for (size_t i = 0; i < n; i++) *interp_field(&work[i], f) = scratch[i];
    }
    /* tod back-fill, stage normalisation, q clamping (before the load calibration) */
    for (size_t i = 0; i < n; i++) sanitize_sample(&work[i]);

    rc = compute_loads(work, n, &cfg.load, load, activity);
    if (rc != ENERGY_OK) goto cleanup;

    energy_ctx_t ctx;
    energy_ctx_init(&ctx, &cfg);
    for (size_t i = 0; i < n; i++) {
        energy_sample_t s = work[i];
        s.load = load[i];
        s.activity = activity[i];
        /* data / no-data is judged on the caller's row, not on the interpolated copy */
        rc = ctx_step_internal(&ctx, &s, &points[i], sample_has_data(&samples[i]));
        if (rc != ENERGY_OK) goto cleanup;
    }
    baseline_flush(&ctx); /* commit the last (open) night / day */

    /* display layer: 30-minute centred smoothing of E, then energy and band */
    for (size_t i = 0; i < n; i++) scratch[i] = points[i].e;
    smooth_display(scratch, n, cfg.display.smooth_minutes, e_smooth);
    for (size_t i = 0; i < n; i++) {
        point_band(&cfg.display, e_smooth[i], points[i].e_sd, points[i].tod, &points[i]);
    }

    double mn = 1e9, mx = -1e9, sum = 0.0, band = 0.0;
    for (size_t i = 0; i < n; i++) {
        if (points[i].energy < mn) mn = points[i].energy;
        if (points[i].energy > mx) mx = points[i].energy;
        sum += points[i].energy;
        band += points[i].upper - points[i].lower;
    }

    /* Daily summary per local calendar day: the first sample's tod anchors midnight
     * (v1.1 cut at floor(t_h/24), which is only right when the series starts at 00:00). */
    size_t nd = 0;
    {
        const double t0 = points[0].t_h, tod0 = points[0].tod;
        size_t i = 0;
        while (i < n) {
            const int day = (int)floor((points[i].t_h - t0 + tod0) / 24.0);
            size_t j = i;
            while (j < n && (int)floor((points[j].t_h - t0 + tod0) / 24.0) == day) j++;
            if (j - i >= 12) daily_summary(&cfg, points, i, j, day, &daily[nd++]);
            i = j;
        }
    }

    out->points = points;
    out->n_points = n;
    out->daily = daily;
    out->n_daily = nd;
    out->stats_min = mn;
    out->stats_max = mx;
    out->stats_mean = sum / (double)n;
    out->stats_band_mean = band / (double)n;
    out->baseline = ctx.bl;
    points = NULL; /* ownership moved to out */
    daily = NULL;

cleanup:
    free(work);
    free(scratch);
    free(load);
    free(activity);
    free(e_smooth);
    free(points);
    free(daily);
    return rc;
}
