/*
 * energy.c —— 能量状态模型（ODE + EKF）的 C 实现
 *
 * 与 Python 版 `energy/src/energy_state/` 的对应关系：
 *   params.py      -> energy_params_t + energy_phi
 *   ode.py         -> energy_advance / 2×2 雅可比 / 解析离散化
 *   load.py        -> load_series / energy_load_raw
 *   observation.py -> energy_obs_t + 设计矩阵
 *   ekf.py         -> ctx_predict / ctx_update（稳健 t、门控、Joseph）
 *   product.py     -> energy_display / 平滑 / 每日汇总
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
#define ENERGY_MAX_HOURS (1.0e6)          /* 时间轴合理范围（约 114 年） */
#define ENERGY_MAX_SUBSTEPS 4096          /* 单个缺口的子步上限，防 CPU 爆炸 */
#define ENERGY_MAX_SMOOTH_MIN 1440.0      /* 展示平滑上限（1 天） */
#define ENERGY_MAX_SAMPLES 10000000u      /* 单次批处理上限，防 size_t 溢出 */

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
/* 默认值（与 Python 先验一致）                                          */
/* ------------------------------------------------------------------ */
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
    memset(c, 0, sizeof(*c));
    c->b_trimp = 1.92;
    c->w_trimp = 1.0;
    c->w_steps = 0.6;
    c->w_dhr = 0.4;
    c->trimp_scale = 10.0;
    c->steps_ref = 120.0;
    c->l_ref = -1.0;   /* 自动 */
    c->hr_rest = -1.0; /* 自动 */
    c->hr_max = -1.0;  /* 自动 */
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
}

void energy_obs_default(energy_obs_t *o)
{
    if (!o) return;
    memset(o, 0, sizeof(*o));
    /* 夜间：[1, E, p_deep, p_rem, p_wake] */
    const double hr_n[5] = {58.0, -12.0, -5.0, 2.0, 9.0};
    const double hrv_n[5] = {3.55, 0.45, 0.30, 0.18, -0.35};
    /* 日间：[1, E, A, cos1, sin1, cos2, sin2] */
    const double hr_d[7] = {80.0, -12.0, 105.0, -2.0, -3.46, 0.75, -1.30};
    const double hrv_d[7] = {3.25, 0.45, -0.55, -0.06, -0.10, 0.025, -0.043};
    memcpy(o->hr_night, hr_n, sizeof(hr_n));
    memcpy(o->hrv_night, hrv_n, sizeof(hrv_n));
    memcpy(o->hr_day, hr_d, sizeof(hr_d));
    memcpy(o->hrv_day, hrv_d, sizeof(hrv_d));
    o->sigma_hr_night = 2.5;
    o->sigma_hr_day = 7.0;
    o->sigma_hrv_night = 0.22;
    o->sigma_hrv_day = 0.40;
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
    if (!(p->phi_min >= 0.0) || !(p->phi_min < p->phi_max)) {
        return cfg_fail(errbuf, errlen, "phi_min/phi_max 非法（需 0 <= min < max）");
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
    const energy_obs_t *o = &cfg->obs;
    if (!(o->sigma_hr_night > 0.0) || !(o->sigma_hr_day > 0.0) ||
        !(o->sigma_hrv_night > 0.0) || !(o->sigma_hrv_day > 0.0)) {
        return cfg_fail(errbuf, errlen, "观测噪声 sigma 必须 > 0");
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
/* 负荷 / 活动量                                                        */
/* ------------------------------------------------------------------ */
double energy_activity_step(double steps, double kcal, double steps_ref)
{
    const double s = clampd(steps / (steps_ref > 0 ? steps_ref : 120.0), 0.0, 1.5);
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
    const double step_frac = clampd(steps / (c->steps_ref > 0 ? c->steps_ref : 120.0), 0.0, 1.5);
    const double r = isfinite(hr) ? clampd((hr - hr_rest) / span, 0.0, 1.0) : 0.0;
    const double r_proxy = clampd(0.6 * step_frac, 0.0, 1.0);
    const double r_eff = isfinite(hr) ? r : r_proxy;
    const double trimp = dt_h * r_eff * exp(c->b_trimp * r_eff) * c->trimp_scale;
    const double dhr = isfinite(hr) ? clampd((hr - hr_base) / span, 0.0, 1.0) : 0.0;
    return c->w_trimp * trimp + c->w_steps * step_frac + c->w_dhr * dhr;
}

/* 由原始序列估算 hr_rest / hr_max / 按钟点基线（与 Python load.py 一致）。 */
static double estimate_hr_rest(const energy_sample_t *s, size_t n)
{
    double *buf = (double *)malloc(n * sizeof(double));
    if (!buf) return 60.0;
    size_t k = 0;
    for (size_t i = 0; i < n; i++) {
        if (isfinite(s[i].hr) && s[i].stage[ENERGY_DEEP] > 0.5) buf[k++] = s[i].hr;
    }
    if (k >= 12) {
        const double v = percentile_sorted(buf, k, 50.0);
        free(buf);
        return v;
    }
    k = 0;
    for (size_t i = 0; i < n; i++) {
        const double m = s[i].stage[ENERGY_LIGHT] + s[i].stage[ENERGY_DEEP] + s[i].stage[ENERGY_REM];
        if (isfinite(s[i].hr) && m > 0.5) buf[k++] = s[i].hr;
    }
    if (k >= 12) {
        const double v = percentile_sorted(buf, k, 20.0);
        free(buf);
        return v;
    }
    k = 0;
    for (size_t i = 0; i < n; i++) {
        if (isfinite(s[i].hr)) buf[k++] = s[i].hr;
    }
    const double v = (k > 0) ? percentile_sorted(buf, k, 5.0) : 60.0;
    free(buf);
    return v;
}

static double estimate_hr_max(const energy_sample_t *s, size_t n)
{
    double *buf = (double *)malloc(n * sizeof(double));
    if (!buf) return 185.0;
    size_t k = 0;
    for (size_t i = 0; i < n; i++) {
        if (isfinite(s[i].hr)) buf[k++] = s[i].hr;
    }
    if (k < 20) {
        free(buf);
        return 185.0;
    }
    const double p = percentile_sorted(buf, k, 99.5);
    free(buf);
    return clampd(p, 150.0, 210.0);
}

/* 按钟点分箱的个人心率基线：剔除 >P90 的运动段 → 分箱中位数 → 环形插值 → 7 点环形平滑 */
static void estimate_hr_base(const energy_sample_t *s, size_t n, double *prof)
{
    double *hr = (double *)malloc(n * sizeof(double));
    if (!hr) {
        for (int i = 0; i < ENERGY_BINS; i++) prof[i] = 70.0;
        return;
    }
    size_t k = 0;
    for (size_t i = 0; i < n; i++) {
        if (isfinite(s[i].hr)) hr[k++] = s[i].hr;
    }
    if (k < (size_t)ENERGY_BINS) {
        const double fill = (k > 0) ? percentile_sorted(hr, k, 50.0) : 70.0;
        for (int i = 0; i < ENERGY_BINS; i++) prof[i] = fill;
        free(hr);
        return;
    }
    const double thr = percentile_sorted(hr, k, 90.0);

    int counts[ENERGY_BINS];
    int pos[ENERGY_BINS];
    int offsets[ENERGY_BINS + 1];
    double *values = (double *)malloc(n * sizeof(double));
    if (!values) {
        for (int i = 0; i < ENERGY_BINS; i++) prof[i] = 70.0;
        free(hr);
        return;
    }
    for (int b = 0; b < ENERGY_BINS; b++) {
        counts[b] = 0;
        pos[b] = 0;
    }
    for (size_t i = 0; i < n; i++) {
        if (!isfinite(s[i].hr) || s[i].hr > thr) continue;
        const int b = tod_bin(s[i].tod);
        counts[b]++;
    }
    offsets[0] = 0;
    for (int b = 0; b < ENERGY_BINS; b++) offsets[b + 1] = offsets[b] + counts[b];
    for (size_t i = 0; i < n; i++) {
        if (!isfinite(s[i].hr) || s[i].hr > thr) continue;
        const int b = tod_bin(s[i].tod);
        values[offsets[b] + pos[b]++] = s[i].hr;
    }
    for (int b = 0; b < ENERGY_BINS; b++) {
        prof[b] = (counts[b] >= 3) ? percentile_sorted(values + offsets[b], (size_t)counts[b], 50.0)
                                   : (double)NAN;
    }
    free(values);
    free(hr);

    /* 环形插值 */
    int good_idx[ENERGY_BINS];
    int ngood = 0;
    for (int b = 0; b < ENERGY_BINS; b++) {
        if (isfinite(prof[b])) good_idx[ngood++] = b;
    }
    if (ngood == 0) {
        for (int b = 0; b < ENERGY_BINS; b++) prof[b] = 70.0;
        return;
    }
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
    /* 环形 7 点平滑 */
    double smooth[ENERGY_BINS];
    for (int b = 0; b < ENERGY_BINS; b++) {
        double acc = 0.0;
        for (int j = -3; j <= 3; j++) {
            acc += prof[((b + j) % ENERGY_BINS + ENERGY_BINS) % ENERGY_BINS];
        }
        smooth[b] = acc / 7.0;
    }
    memcpy(prof, smooth, sizeof(smooth));
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
/* 观测模型                                                            */
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

/* 观测预测 + 对状态的雅可比行 + 噪声 σ */
static int obs_predict(const energy_config_t *cfg, const energy_sample_t *s,
                       double e_state, int night, int channel,
                       double *pred, double hvec[2], double *sigma)
{
    const energy_obs_t *o = &cfg->obs;
    if (night) {
        const double *coef = (channel == 0) ? o->hr_night : o->hrv_night;
        const double x[5] = {1.0, e_state, s->stage[ENERGY_DEEP], s->stage[ENERGY_REM],
                             s->stage[ENERGY_WAKE]};
        double acc = 0.0;
        for (int i = 0; i < 5; i++) acc += coef[i] * x[i];
        *pred = acc;
        hvec[0] = 0.0;
        hvec[1] = coef[1];
        *sigma = (channel == 0) ? o->sigma_hr_night : o->sigma_hrv_night;
    } else {
        const double *coef = (channel == 0) ? o->hr_day : o->hrv_day;
        const double w = 2.0 * M_PI / 24.0;
        const double t = wrap24(s->tod);
        const double x[7] = {1.0, e_state, s->activity,
                             cos(w * t), sin(w * t),
                             cos(2.0 * w * t), sin(2.0 * w * t)};
        double acc = 0.0;
        for (int i = 0; i < 7; i++) acc += coef[i] * x[i];
        *pred = acc;
        hvec[0] = 0.0;
        hvec[1] = coef[1];
        *sigma = (channel == 0) ? o->sigma_hr_day : o->sigma_hrv_day;
    }
    return isfinite(*pred) ? 1 : 0;
}

/* ------------------------------------------------------------------ */
/* EKF                                                                 */
/* ------------------------------------------------------------------ */
static void ctx_predict(energy_ctx_t *ctx, double h)
{
    const energy_config_t *cfg = &ctx->cfg;
    const stage_mix_t mix = stage_mix(&cfg->params, ctx->prev_stage);
    double F[4], Phi[4];
    continuous_jacobian(&cfg->params, ctx->s, ctx->e, mix.m, mix.lam, mix.psi,
                        ctx->prev_load, ctx->prev_stress, F);
    discrete_jacobian(F, h, Phi);

    double s_new, e_new;
    advance_raw(&cfg->params, ctx->s, ctx->e, mix.m, mix.lam, mix.psi,
                ctx->prev_load, ctx->prev_stress, h, &s_new, &e_new);

    /* P = Φ P Φᵀ + Q·h（Q 对角） */
    const double p00 = Phi[0] * ctx->p00 + Phi[1] * ctx->p01;
    const double p01 = Phi[0] * ctx->p01 + Phi[1] * ctx->p11;
    const double p10 = Phi[2] * ctx->p00 + Phi[3] * ctx->p01;
    const double p11 = Phi[2] * ctx->p01 + Phi[3] * ctx->p11;
    ctx->p00 = p00 * Phi[0] + p01 * Phi[1] + cfg->ekf.q_s * h;
    ctx->p01 = p00 * Phi[2] + p01 * Phi[3];
    ctx->p11 = p10 * Phi[2] + p11 * Phi[3] + cfg->ekf.q_e * h;
    ctx->p01 = 0.5 * (ctx->p01 + ctx->p01);
    if (!isfinite(ctx->p00) || !isfinite(ctx->p11) || ctx->p00 <= 0.0 || ctx->p11 <= 0.0) {
        /* 极长缺口/极端观测让协方差失去意义时，回到冷启动方差 */
        ctx->p00 = cfg->ekf.p0_s;
        ctx->p01 = 0.0;
        ctx->p11 = cfg->ekf.p0_e;
    }
    ctx->s = clampd(s_new, 0.0, 1.0);
    ctx->e = clampd(e_new, 0.0, 1.0);
}

static void ctx_update(energy_ctx_t *ctx, const energy_sample_t *s, int channel)
{
    const energy_config_t *cfg = &ctx->cfg;
    const double y = (channel == 0) ? s->hr : s->hrv_log;
    const double q = (channel == 0) ? finite_or(s->q_hr, 0.0) : finite_or(s->q_hrv, 0.0);
    if (!isfinite(y) || q <= 0.0) return;

    const int night = is_night(s->stage);
    if (!cfg->day_update && !night) {
        /* Python: skip_update = (not day_update) and panel.is_day[k] */
        const int day = (!night) && s->tod >= 7.0 && s->tod < 22.0;
        if (day) return;
    }

    double pred, hvec[2], sigma;
    if (!obs_predict(cfg, s, ctx->e, night, channel, &pred, hvec, &sigma)) return;

    const double R = (sigma * sigma) / fmax(q, cfg->ekf.q_floor);
    const double PHt0 = ctx->p00 * hvec[0] + ctx->p01 * hvec[1];
    const double PHt1 = ctx->p01 * hvec[0] + ctx->p11 * hvec[1];
    const double S = hvec[0] * PHt0 + hvec[1] * PHt1 + R;
    if (!isfinite(S) || S <= 1e-12) return;

    const double innov = y - pred;
    const double z = innov / sqrt(S);
    double R_use = R, S_use = S;
    if (cfg->ekf.robust_nu > 0.0) {
        const double inflate = (cfg->ekf.robust_nu + z * z) / (cfg->ekf.robust_nu + 1.0);
        R_use = R * inflate;
        S_use = hvec[0] * PHt0 + hvec[1] * PHt1 + R_use;
    }
    if (fabs(z) > cfg->ekf.gate_sigma) return;

    const double K0 = PHt0 / S_use;
    const double K1 = PHt1 / S_use;
    ctx->s += K0 * innov;
    ctx->e += K1 * innov;

    /* Joseph 形式 P = (I-KH)P(I-KH)ᵀ + K Kᵀ R */
    const double a00 = 1.0 - K0 * hvec[0], a01 = -K0 * hvec[1];
    const double a10 = -K1 * hvec[0], a11 = 1.0 - K1 * hvec[1];
    const double t00 = a00 * ctx->p00 + a01 * ctx->p01;
    const double t01 = a00 * ctx->p01 + a01 * ctx->p11;
    const double t10 = a10 * ctx->p00 + a11 * ctx->p01;
    const double t11 = a10 * ctx->p01 + a11 * ctx->p11;
    ctx->p00 = t00 * a00 + t01 * a01 + K0 * K0 * R_use;
    ctx->p01 = t00 * a10 + t01 * a11 + K0 * K1 * R_use;
    ctx->p11 = t10 * a10 + t11 * a11 + K1 * K1 * R_use;
    ctx->p01 = 0.5 * (ctx->p01 + ctx->p01);
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
/* 流式 API                                                            */
/* ------------------------------------------------------------------ */
void energy_ctx_init(energy_ctx_t *ctx, const energy_config_t *cfg)
{
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
    ctx->p00 = ctx->cfg.ekf.p0_s;
    ctx->p01 = 0.0;
    ctx->p11 = ctx->cfg.ekf.p0_e;
    ctx->started = 0;
}

static double sample_load(const energy_config_t *cfg, const energy_sample_t *s, double dt_h)
{
    if (isfinite(s->load) && s->load >= 0.0) return clampd(s->load, 0.0, 1.0);
    const double hr_rest = cfg->load.hr_rest > 0 ? cfg->load.hr_rest : 60.0;
    const double hr_max = cfg->load.hr_max > 0 ? cfg->load.hr_max : 185.0;
    const double raw = energy_load_raw(&cfg->load, s->hr, hr_rest, hr_max, hr_rest, s->steps, dt_h);
    const double l_ref = cfg->load.l_ref > 0 ? cfg->load.l_ref : 0.25;
    return clampd(raw / l_ref, 0.0, 1.0);
}

static double sample_activity(const energy_config_t *cfg, const energy_sample_t *s)
{
    if (isfinite(s->activity) && s->activity >= 0.0) return clampd(s->activity, 0.0, 1.0);
    return energy_activity_step(s->steps, s->kcal, cfg->load.steps_ref);
}

static void ctx_emit(const energy_ctx_t *ctx, const energy_sample_t *s, energy_point_t *out)
{
    if (!out) return;
    const stage_mix_t mix = stage_mix(&ctx->cfg.params, s->stage);
    memset(out, 0, sizeof(*out));
    out->t_h = s->t_h;
    out->tod = s->tod;
    out->s = ctx->s;
    out->e = ctx->e;
    out->e_sd = ctx->cfg.use_ekf ? sqrt(fmax(ctx->p11, 0.0)) : 0.0;
    out->energy = energy_display(&ctx->cfg.display, ctx->e, s->tod);
    out->base = 100.0 * ctx->e;
    const double mod = clampd(1.0 + ctx->cfg.display.rho * energy_circadian(&ctx->cfg.display, s->tod),
                              0.0, 2.0);
    const double z = ctx->cfg.display.band_z;
    out->lower = clampd(100.0 * clampd(ctx->e - z * out->e_sd, 0.0, 1.0) * mod, 0.0, 100.0);
    out->upper = clampd(100.0 * clampd(ctx->e + z * out->e_sd, 0.0, 1.0) * mod, 0.0, 100.0);
    out->load = sample_load(&ctx->cfg, s, ENERGY_HOURS_PER_STEP);
    out->activity = sample_activity(&ctx->cfg, s);
    out->sleep_prob = mix.sleep_prob;
}

static void sanitize_sample(energy_sample_t *s)
{
    for (int i = 0; i < 4; i++) {
        if (!isfinite(s->stage[i]) || s->stage[i] < 0.0) s->stage[i] = 0.0;
        s->stage[i] = clampd(s->stage[i], 0.0, 1.0);
    }
    double sum = s->stage[0] + s->stage[1] + s->stage[2] + s->stage[3];
    if (sum > 1.0) {
        for (int i = 0; i < 4; i++) s->stage[i] /= sum;
    }
    if (!isfinite(s->tod) || s->tod < 0.0) s->tod = wrap24(s->t_h);
    if (!isfinite(s->tod)) s->tod = 0.0;
    if (!isfinite(s->stress)) s->stress = 0.0;
    if (s->stress > 1.5) s->stress /= 100.0;
    s->stress = clampd(s->stress, 0.0, 1.0);
    if (!isfinite(s->steps)) s->steps = 0.0;
    if (!isfinite(s->kcal)) s->kcal = 0.0;
    if (!isfinite(s->q_hr)) s->q_hr = 0.0;
    if (!isfinite(s->q_hrv)) s->q_hrv = 0.0;
    s->q_hr = clampd(s->q_hr, 0.0, 1.0);
    s->q_hrv = clampd(s->q_hrv, 0.0, 1.0);
}

static int ctx_step_internal(energy_ctx_t *ctx, const energy_sample_t *s_in, energy_point_t *out)
{
    if (ctx->config_error) return ENERGY_EINVAL;
    energy_sample_t local = *s_in;
    sanitize_sample(&local);
    if (!isfinite(local.t_h) || fabs(local.t_h) > ENERGY_MAX_HOURS) return ENERGY_EINVAL;

    /* 本步覆盖的时长：首样本用 5 分钟，之后用实际间隔（缺数据时可能更长） */
    double h = ENERGY_HOURS_PER_STEP;
    if (!ctx->started) {
        ctx->t_h = local.t_h;
        ctx->started = 1;
    } else {
        h = local.t_h - ctx->t_h;
        if (h <= 0.0) return ENERGY_EDOMAIN;
        if (ctx->cfg.use_ekf) {
            ctx_predict(ctx, h);
        } else {
            const stage_mix_t mix = stage_mix(&ctx->cfg.params, ctx->prev_stage);
            double s1, e1;
            simulate_advance(&ctx->cfg.params, ctx->s, ctx->e, mix.m, mix.lam, mix.psi,
                             ctx->prev_load, ctx->prev_stress, h, &s1, &e1);
            ctx->s = s1;
            ctx->e = e1;
        }
        ctx->t_h = local.t_h;
    }

    /* 先把“本样本生效的负荷/活动量”落回样本：
     * load/activity 传 <0 表示由库计算，若只把它们算进输出而不写回样本，
     * 观测方程会拿到 -1（活动量）或未归一化值，导致日间 HR/HRV 被错误解释。 */
    local.load = sample_load(&ctx->cfg, &local, h);
    local.activity = sample_activity(&ctx->cfg, &local);

    ctx->prev_load = local.load;
    ctx->prev_stress = clampd(local.stress, 0.0, 1.0);
    memcpy(ctx->prev_stage, local.stage, sizeof(ctx->prev_stage));

    if (ctx->cfg.use_ekf) {
        ctx_update(ctx, &local, 0);
        ctx_update(ctx, &local, 1);
        ctx_clamp(ctx);
    }
    ctx_emit(ctx, &local, out);
    return ENERGY_OK;
}

int energy_ctx_step(energy_ctx_t *ctx, const energy_sample_t *sample, energy_point_t *out)
{
    if (!ctx || !sample) return ENERGY_EINVAL;
    return ctx_step_internal(ctx, sample, out);
}

void energy_ctx_peek(const energy_ctx_t *ctx, double tod, energy_point_t *out)
{
    if (!ctx || !out) return;
    if (!isfinite(tod) || tod < 0.0) {
        tod = isfinite(ctx->t_h) ? wrap24(ctx->t_h) : 0.0;
    }
    memset(out, 0, sizeof(*out));
    out->t_h = ctx->t_h;
    out->tod = tod;
    out->s = ctx->s;
    out->e = ctx->e;
    out->e_sd = ctx->cfg.use_ekf ? sqrt(fmax(ctx->p11, 0.0)) : 0.0;
    out->energy = energy_display(&ctx->cfg.display, ctx->e, tod);
    out->base = 100.0 * ctx->e;
    const double mod = clampd(1.0 + ctx->cfg.display.rho * energy_circadian(&ctx->cfg.display, tod),
                              0.0, 2.0);
    const double z = ctx->cfg.display.band_z;
    out->lower = clampd(100.0 * clampd(ctx->e - z * out->e_sd, 0.0, 1.0) * mod, 0.0, 100.0);
    out->upper = clampd(100.0 * clampd(ctx->e + z * out->e_sd, 0.0, 1.0) * mod, 0.0, 100.0);
}

/* ------------------------------------------------------------------ */
/* 批量                                                                */
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

/* 原始负荷序列 + 标定量（与 Python load.py 对齐） */
static int load_raw_series(const energy_sample_t *samples, size_t n,
                           const energy_load_cfg_t *cfg, double *raw,
                           double *hr_rest_out, double *hr_max_out, double *l_ref_out)
{
    const double hr_rest = cfg->hr_rest > 0 ? cfg->hr_rest : estimate_hr_rest(samples, n);
    const double hr_max = cfg->hr_max > 0 ? cfg->hr_max : estimate_hr_max(samples, n);
    double prof[ENERGY_BINS];
    estimate_hr_base(samples, n, prof);

    for (size_t i = 0; i < n; i++) {
        const double dt = (i == 0) ? ENERGY_HOURS_PER_STEP
                                   : clampd(samples[i].t_h - samples[i - 1].t_h, 0.0, 0.5);
        const int b = tod_bin(samples[i].tod);
        raw[i] = energy_load_raw(cfg, samples[i].hr, hr_rest, hr_max, prof[b],
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
    return ENERGY_OK;
}

static int compute_loads(const energy_sample_t *samples, size_t n,
                         const energy_load_cfg_t *cfg,
                         double *load, double *activity)
{
    double *raw = (double *)malloc(n * sizeof(double));
    if (!raw) return ENERGY_ENOMEM;
    double l_ref = 0.0;
    const int rc = load_raw_series(samples, n, cfg, raw, NULL, NULL, &l_ref);
    if (rc != ENERGY_OK) {
        free(raw);
        return rc;
    }
    for (size_t i = 0; i < n; i++) {
        load[i] = clampd(raw[i] / l_ref, 0.0, 1.0);
        activity[i] = energy_activity_step(samples[i].steps, samples[i].kcal, cfg->steps_ref);
    }
    free(raw);
    return ENERGY_OK;
}

int energy_load_calibrate(const energy_sample_t *samples, size_t n,
                          const energy_load_cfg_t *cfg, energy_load_stats_t *out)
{
    if (!samples || n == 0 || !out || !cfg) return ENERGY_EINVAL;
    double *raw = (double *)malloc(n * sizeof(double));
    if (!raw) return ENERGY_ENOMEM;
    const int rc = load_raw_series(samples, n, cfg, raw, &out->hr_rest, &out->hr_max, &out->l_ref);
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
    snprintf(out, 11, "%02d-%02d", m, d);
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

    /* 时间轴检查：有限、单调、量级合理 */
    for (size_t i = 0; i < n; i++) {
        const double t = samples[i].t_h;
        if (!isfinite(t) || fabs(t) > ENERGY_MAX_HOURS) return ENERGY_EINVAL;
        if (i > 0 && t <= samples[i - 1].t_h) return ENERGY_EDOMAIN;
    }

    /* 拷贝一份可变输入，做缺失值补全（与 Python API 的线性插值一致） */
    energy_sample_t *work = (energy_sample_t *)malloc(n * sizeof(energy_sample_t));
    double *scratch = (double *)malloc(n * sizeof(double));
    energy_point_t *points = (energy_point_t *)calloc(n, sizeof(energy_point_t));
    double *load = (double *)malloc(n * sizeof(double));
    double *activity = (double *)malloc(n * sizeof(double));
    if (!work || !scratch || !points || !load || !activity) {
        free(work);
        free(scratch);
        free(points);
        free(load);
        free(activity);
        return ENERGY_ENOMEM;
    }
    memcpy(work, samples, n * sizeof(energy_sample_t));
    {
        double *fields[3];
        for (size_t i = 0; i < n; i++) {
            scratch[i] = work[i].stress;
        }
        fill_linear(scratch, n);
        for (size_t i = 0; i < n; i++) work[i].stress = scratch[i];
        fields[0] = scratch; /* 仅为消除未使用告警的占位，实际逐列处理 */
        (void)fields;
        for (size_t i = 0; i < n; i++) scratch[i] = work[i].steps;
        fill_linear(scratch, n);
        for (size_t i = 0; i < n; i++) work[i].steps = scratch[i];
        for (size_t i = 0; i < n; i++) scratch[i] = work[i].kcal;
        fill_linear(scratch, n);
        for (size_t i = 0; i < n; i++) work[i].kcal = scratch[i];
    }
    /* 归一化/清洗：tod 回填、分期归一、q 裁剪（必须在负荷计算前完成） */
    for (size_t i = 0; i < n; i++) {
        sanitize_sample(&work[i]);
    }

    int rc = compute_loads(work, n, &cfg.load, load, activity);
    if (rc != ENERGY_OK) {
        free(work);
        free(scratch);
        free(points);
        free(load);
        free(activity);
        return rc;
    }

    energy_ctx_t ctx;
    energy_ctx_init(&ctx, &cfg);
    for (size_t i = 0; i < n; i++) {
        energy_sample_t s = work[i];
        s.load = load[i];
        s.activity = activity[i];
        rc = ctx_step_internal(&ctx, &s, &points[i]);
        if (rc != ENERGY_OK) {
            free(work);
            free(scratch);
            free(points);
            free(load);
            free(activity);
            return rc;
        }
        points[i].load = load[i];
        points[i].activity = activity[i];
    }

    /* 展示层：对 E 做 30 分钟居中平滑后重算能量与置信带 */
    double *e_smooth = (double *)malloc(n * sizeof(double));
    if (!e_smooth) {
        free(work);
        free(scratch);
        free(points);
        free(load);
        free(activity);
        return ENERGY_ENOMEM;
    }
    double *e_raw = (double *)malloc(n * sizeof(double));
    if (!e_raw) {
        free(e_smooth);
        free(work);
        free(scratch);
        free(points);
        free(load);
        free(activity);
        return ENERGY_ENOMEM;
    }
    for (size_t i = 0; i < n; i++) e_raw[i] = points[i].e;
    smooth_display(e_raw, n, cfg.display.smooth_minutes, e_smooth);

    const double zband = cfg.display.band_z;
    for (size_t i = 0; i < n; i++) {
        const double e = clampd(e_smooth[i], 0.0, 1.0);
        const double mod = clampd(1.0 + cfg.display.rho * energy_circadian(&cfg.display, points[i].tod),
                                  0.0, 2.0);
        const double sd = points[i].e_sd;
        points[i].base = 100.0 * e;
        points[i].energy = clampd(100.0 * e * mod, 0.0, 100.0);
        points[i].lower = clampd(100.0 * clampd(e - zband * sd, 0.0, 1.0) * mod, 0.0, 100.0);
        points[i].upper = clampd(100.0 * clampd(e + zband * sd, 0.0, 1.0) * mod, 0.0, 100.0);
    }
    free(e_raw);
    free(e_smooth);

    /* 统计量 */
    double mn = 1e9, mx = -1e9, sum = 0.0, band = 0.0;
    for (size_t i = 0; i < n; i++) {
        if (points[i].energy < mn) mn = points[i].energy;
        if (points[i].energy > mx) mx = points[i].energy;
        sum += points[i].energy;
        band += points[i].upper - points[i].lower;
    }

    /* 每日汇总（自然日 = floor(t_h/24)） */
    energy_daily_t *daily = (energy_daily_t *)calloc(n / 12 + 2, sizeof(energy_daily_t));
    if (!daily) {
        free(e_raw);
        free(e_smooth);
        free(work);
        free(scratch);
        free(points);
        free(load);
        free(activity);
        return ENERGY_ENOMEM;
    }
    size_t nd = 0;
    size_t i = 0;
    while (i < n) {
        const int day = (int)floor(points[i].t_h / 24.0);
        size_t j = i;
        while (j < n && (int)floor(points[j].t_h / 24.0) == day) j++;
        const size_t cnt = j - i;
        if (cnt >= 12) {
            energy_daily_t *d = &daily[nd++];
            d->day = day;
            d->n = (cnt > (size_t)INT_MAX) ? INT_MAX : (int)cnt;
            date_label(&cfg, day, d->date);

            double dmin = 1e9, dmax = -1e9, dsum = 0.0, debt = 0.0;
            long morning = -1, evening = -1;
            long peak = -1, dip = -1;
            double peak_v = -1e9, dip_v = 1e9;
            for (size_t k = i; k < j; k++) {
                const energy_point_t *p = &points[k];
                if (p->energy < dmin) dmin = p->energy;
                if (p->energy > dmax) dmax = p->energy;
                dsum += p->energy;
                if (p->s > debt) debt = p->s;
                if (p->tod < 9.0) morning = (long)k;
                if (p->tod > 20.0) evening = (long)k;
                if (p->tod >= 6.0 && p->tod <= 16.0 && p->energy > peak_v) {
                    peak_v = p->energy;
                    peak = (long)k;
                }
                if (p->tod >= 13.0 && p->tod <= 17.0 && p->energy < dip_v) {
                    dip_v = p->energy;
                    dip = (long)k;
                }
            }
            if (peak < 0) {
                peak_v = -1e9;
                for (size_t k = i; k < j; k++) {
                    if (points[k].energy > peak_v) {
                        peak_v = points[k].energy;
                        peak = (long)k;
                    }
                }
            }
            if (dip < 0) {
                dip_v = 1e9;
                for (size_t k = i; k < j; k++) {
                    if (points[k].energy < dip_v) {
                        dip_v = points[k].energy;
                        dip = (long)k;
                    }
                }
            }
            if (morning < 0) morning = (long)i;
            if (evening < 0) evening = (long)(j - 1);
            d->morning = points[morning].energy;
            d->evening = points[evening].energy;
            d->mean = dsum / (double)cnt;
            d->min = dmin;
            d->max = dmax;
            d->peak_tod = points[peak].tod;
            d->dip_tod = points[dip].tod;
            d->sleep_debt = debt;
        }
        i = j;
    }

    out->points = points;
    out->n_points = n;
    out->daily = daily;
    out->n_daily = nd;
    out->stats_min = mn;
    out->stats_max = mx;
    out->stats_mean = sum / (double)n;
    out->stats_band_mean = band / (double)n;

    free(work);
    free(scratch);
    free(load);
    free(activity);
    return ENERGY_OK;
}
