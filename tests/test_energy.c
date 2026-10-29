/*
 * test_energy.c —— C 版能量模型的单元测试（无第三方依赖）
 *
 * 运行：make test
 */

#include "energy.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static int g_fail = 0;
static int g_checks = 0;

#define CHECK(cond)                                                          \
    do {                                                                     \
        g_checks++;                                                          \
        if (!(cond)) {                                                       \
            printf("FAIL %s:%d  %s\n", __FILE__, __LINE__, #cond);           \
            g_fail++;                                                        \
        }                                                                    \
    } while (0)

#define CHECK_NEAR(a, b, tol)                                                \
    do {                                                                     \
        g_checks++;                                                          \
        const double _a = (a), _b = (b);                                     \
        if (!(fabs(_a - _b) <= (tol))) {                                     \
            printf("FAIL %s:%d  %s=%.12f vs %.12f (tol %g)\n", __FILE__,     \
                   __LINE__, #a, _a, _b, (double)(tol));                     \
            g_fail++;                                                        \
        }                                                                    \
    } while (0)

/* ------------------------------------------------------------------ */
/* 构造 3 天、5 分钟粒度的测试面板（睡眠 23:00–07:00）                    */
/* ------------------------------------------------------------------ */
static energy_sample_t *make_panel(size_t days, size_t *n_out)
{
    const size_t n = 288 * days;
    energy_sample_t *s = (energy_sample_t *)calloc(n, sizeof(energy_sample_t));
    if (!s) return NULL;
    for (size_t i = 0; i < n; i++) {
        const double t = (double)i / 12.0;
        const double tod = fmod(t, 24.0);
        const int asleep = (tod >= 23.0 || tod < 7.0);
        s[i].t_h = t;
        s[i].tod = tod;
        s[i].stage[ENERGY_LIGHT] = asleep ? 1.0 : 0.0;
        s[i].hr = asleep ? 56.0 : 72.0;
        s[i].hrv_log = asleep ? log(45.0) : log(32.0);
        s[i].steps = asleep ? 0.0 : 40.0;
        s[i].kcal = 0.0;
        s[i].stress = 0.25;
        s[i].q_hr = 1.0;
        s[i].q_hrv = 1.0;
        s[i].load = -1.0;
        s[i].activity = -1.0;
    }
    *n_out = n;
    return s;
}

/* ------------------------------------------------------------------ */
static void test_defaults(void)
{
    energy_params_t p;
    energy_params_default(&p);
    CHECK_NEAR(p.kappa_w, 0.10, 1e-12);
    CHECK_NEAR(p.kappa_s[ENERGY_DEEP], 0.42, 1e-12);
    CHECK_NEAR(p.alpha_sleep, 0.30, 1e-12);
    CHECK_NEAR(p.beta_S, 0.035, 1e-12);
    CHECK_NEAR(p.w_stage[ENERGY_DEEP], 1.25, 1e-12);
    CHECK_NEAR(p.s_ref, 0.50, 1e-12);
    CHECK_NEAR(p.s_floor, 0.10, 1e-12);
}

static void test_phi(void)
{
    energy_params_t p;
    energy_params_default(&p);
    /* 与 Python 端一致：φ(0)=0.341, φ(0.5)=1.0, φ(0.8)=1.275, φ(1)=1.439 */
    CHECK_NEAR(energy_phi(&p, 0.0), 0.3413, 5e-4);
    CHECK_NEAR(energy_phi(&p, 0.5), 1.0, 1e-9);
    CHECK_NEAR(energy_phi(&p, 0.8), 1.2754, 5e-4);
    CHECK_NEAR(energy_phi(&p, 1.0), 1.4387, 5e-4);
    /* 单调不减、有上下限 */
    double prev = 0.0;
    for (int i = 0; i <= 100; i++) {
        const double v = energy_phi(&p, (double)i / 100.0);
        CHECK(v >= prev - 1e-12);
        CHECK(v >= p.phi_min - 1e-12 && v <= p.phi_max + 1e-12);
        prev = v;
    }
}

static void test_advance(void)
{
    energy_params_t p;
    energy_params_default(&p);
    double s, e;

    /* 睡眠 8 小时：S 下降、E 上升，且都在 [0,1] */
    energy_advance(&p, 0.6, 0.4, 1.0, 0.26, 1.0, 0.0, 0.0, 8.0, &s, &e);
    CHECK(s < 0.6);
    CHECK(e > 0.4);
    CHECK(s >= 0.0 && s <= 1.0 && e >= 0.0 && e <= 1.0);

    /* 清醒 16 小时：S 上升、E 下降 */
    energy_advance(&p, 0.2, 0.9, 0.0, 0.05, 0.15, 0.3, 0.2, 16.0, &s, &e);
    CHECK(s > 0.2);
    CHECK(e < 0.9);
    CHECK(e >= 0.0);

    /* 极小步长与极大步长都不越界 */
    for (int i = 0; i < 1000; i++) {
        const double h = pow(10.0, -3.0 + 6.0 * (double)i / 999.0);
        energy_advance(&p, 0.99, 0.01, 0.5, 0.26, 1.0, 1.0, 1.0, h, &s, &e);
        CHECK(s >= 0.0 && s <= 1.0);
        CHECK(e >= 0.0 && e <= 1.0);
    }
}

static void test_equilibrium(void)
{
    energy_params_t p;
    energy_params_default(&p);
    CHECK_NEAR(energy_wake_equilibrium(&p, 0.5, 0.0, 0.0), 0.4211, 5e-4);
    CHECK_NEAR(energy_sleep_equilibrium(&p, 0.5, 1.0), 0.9677, 2e-3);
}

static void test_display(void)
{
    energy_display_cfg_t c;
    energy_display_default(&c);
    /* 节律：10:00 接近峰，22:00 接近 -1，归一化后 |c| <= 1 */
    CHECK_NEAR(energy_circadian(&c, 10.0), 0.987, 5e-3);
    CHECK(energy_circadian(&c, 22.0) < -0.99);
    CHECK_NEAR(energy_circadian(&c, 14.5), -0.317, 2e-2);
    for (int i = 0; i < 288; i++) {
        const double v = energy_circadian(&c, (double)i / 12.0);
        CHECK(fabs(v) <= 1.0 + 1e-9);
    }
    /* E=1 时能量不超过 100 */
    CHECK(energy_display(&c, 1.0, 10.0) <= 100.0);
    CHECK(energy_display(&c, 1.0, 22.0) <= 100.0);
    CHECK(energy_display(&c, 0.0, 10.0) == 0.0);
}

static void test_load(void)
{
    energy_load_cfg_t c;
    energy_load_default(&c);
    /* 心率越高、步数越多，原始负荷越大 */
    const double r0 = energy_load_raw(&c, 60.0, 55.0, 185.0, 65.0, 0.0, 5.0 / 60.0);
    const double r1 = energy_load_raw(&c, 120.0, 55.0, 185.0, 65.0, 0.0, 5.0 / 60.0);
    const double r2 = energy_load_raw(&c, 120.0, 55.0, 185.0, 65.0, 200.0, 5.0 / 60.0);
    CHECK(r1 > r0);
    CHECK(r2 > r1);
    /* HR 缺失时退化为步数代理，不应为 NaN */
    const double r3 = energy_load_raw(&c, NAN, 55.0, 185.0, 65.0, 60.0, 5.0 / 60.0);
    CHECK(isfinite(r3) && r3 > 0.0);
    CHECK(energy_load_raw(&c, 300.0, 55.0, 185.0, 65.0, 500.0, 5.0 / 60.0) > 0.0);

    /* 活动量在 [0,1] */
    CHECK_NEAR(energy_activity_step(0.0, 0.0, 120.0), 0.0, 1e-12);
    CHECK(energy_activity_step(120.0, 6.0, 120.0) <= 1.0);
    CHECK(energy_activity_step(1000.0, 1000.0, 120.0) == 1.0);

    /* v1.2：步数参考 500 步/5 min（约 100 步/分钟）。120 步/5 min 只是慢慢走动，
     * 不能再被当成快走（v1.1 在 120 步就饱和到 A=0.8）；真实快走 500–600 步才接近满量程。 */
    CHECK_NEAR(c.steps_ref, 500.0, 1e-12);
    CHECK_NEAR(energy_activity_step(120.0, 0.0, c.steps_ref), 0.8 * 120.0 / 500.0, 1e-12);
    CHECK(energy_activity_step(120.0, 0.0, c.steps_ref) < 0.25);
    CHECK(energy_activity_step(550.0, 0.0, c.steps_ref) > 0.8);
    CHECK_NEAR(energy_activity_step(120.0, 0.0, 0.0), energy_activity_step(120.0, 0.0, 500.0), 1e-12);
    /* 负荷的步数项同样按 500 归一：120 步远小于 550 步 */
    const double w0 = energy_load_raw(&c, NAN, 55.0, 185.0, 65.0, 120.0, 5.0 / 60.0);
    const double w1 = energy_load_raw(&c, NAN, 55.0, 185.0, 65.0, 550.0, 5.0 / 60.0);
    CHECK(w1 > 3.0 * w0);
}

static void test_batch_basic(void)
{
    size_t n = 0;
    energy_sample_t *s = make_panel(3, &n);
    CHECK(s != NULL);

    energy_config_t cfg;
    energy_config_default(&cfg);
    cfg.start_year = 2026;
    cfg.start_month = 9;
    cfg.start_day = 1;
    energy_output_t out;
    const int rc = energy_compute(s, n, &cfg, &out);
    CHECK(rc == ENERGY_OK);
    CHECK(out.n_points == n);
    CHECK(out.n_daily == 3);

    for (size_t i = 0; i < out.n_points; i++) {
        const energy_point_t *p = &out.points[i];
        CHECK(isfinite(p->energy));
        CHECK(p->energy >= 0.0 && p->energy <= 100.0);
        CHECK(p->lower <= p->energy + 1e-9);
        CHECK(p->upper >= p->energy - 1e-9);
        CHECK(p->e >= 0.0 && p->e <= 1.0);
        CHECK(p->s >= 0.0 && p->s <= 1.0);
        CHECK(p->load >= 0.0 && p->load <= 1.0);
        CHECK(p->activity >= 0.0 && p->activity <= 1.0);
    }
    CHECK(out.stats_min >= 0.0 && out.stats_max <= 100.0);
    CHECK(out.stats_mean > 0.0 && out.stats_mean < 100.0);

    for (size_t d = 0; d < out.n_daily; d++) {
        CHECK(out.daily[d].n >= 12);
        CHECK(out.daily[d].min <= out.daily[d].mean + 1e-9);
        CHECK(out.daily[d].mean <= out.daily[d].max + 1e-9);
        CHECK(out.daily[d].morning >= 0.0 && out.daily[d].morning <= 100.0);
        CHECK(out.daily[d].sleep_debt >= 0.0 && out.daily[d].sleep_debt <= 1.0);
        CHECK(strcmp(out.daily[d].date, "") != 0);
    }
    energy_output_free(&out);
    free(s);
}

static void test_stream_equals_batch(void)
{
    size_t n = 0;
    energy_sample_t *s = make_panel(3, &n);
    CHECK(s != NULL);

    energy_config_t cfg;
    energy_config_default(&cfg);
    energy_output_t out;
    CHECK(energy_compute(s, n, &cfg, &out) == ENERGY_OK);

    /* 把批量算出的 load/activity 喂给流式，状态应逐点一致 */
    energy_ctx_t ctx;
    energy_ctx_init(&ctx, &cfg);
    for (size_t i = 0; i < n; i++) {
        energy_sample_t one = s[i];
        one.load = out.points[i].load;
        one.activity = out.points[i].activity;
        energy_point_t p;
        CHECK(energy_ctx_step(&ctx, &one, &p) == ENERGY_OK);
        CHECK_NEAR(p.s, out.points[i].s, 1e-12);
        CHECK_NEAR(p.e, out.points[i].e, 1e-12);
        CHECK_NEAR(p.e_sd, out.points[i].e_sd, 1e-12);
        CHECK(p.flags == out.points[i].flags);
    }
    /* 个体基线学习也是同一套代码：结束时学到的基线一致 */
    energy_baseline_t bl;
    energy_ctx_get_baseline(&ctx, &bl);
    for (int k = 0; k < 4; k++) {
        CHECK(bl.n[k] == out.baseline.n[k]);
        if (bl.n[k] > 0) CHECK_NEAR(bl.level[k], out.baseline.level[k], 1e-12);
    }
    CHECK(out.baseline.n[ENERGY_BL_HR_NIGHT] >= 2);
    energy_output_free(&out);
    free(s);
}

static void test_gap_and_openloop(void)
{
    size_t n = 0;
    energy_sample_t *s = make_panel(1, &n);
    CHECK(s != NULL);

    /* 挖掉中间 6 小时（t=6h..12h），时间轴仍然连续 */
    for (size_t i = 0; i < n; i++) {
        if (s[i].t_h >= 6.0 && s[i].t_h < 12.0) {
            s[i].hr = NAN;
            s[i].hrv_log = NAN;
            s[i].q_hr = 0.0;
            s[i].q_hrv = 0.0;
        }
    }
    energy_config_t cfg;
    energy_config_default(&cfg);
    energy_output_t out;
    CHECK(energy_compute(s, n, &cfg, &out) == ENERGY_OK);
    for (size_t i = 0; i < out.n_points; i++) {
        CHECK(isfinite(out.points[i].energy));
        CHECK(out.points[i].energy >= 0.0 && out.points[i].energy <= 100.0);
    }
    /* 缺口内置信带应变宽（相对无缺口） */
    double band_gap = 0.0;
    for (size_t i = 0; i < out.n_points; i++) {
        if (s[i].t_h >= 6.0 && s[i].t_h < 12.0) {
            band_gap += out.points[i].upper - out.points[i].lower;
        }
    }
    CHECK(band_gap > 0.0);
    energy_output_free(&out);

    /* 开环：没有观测也能算，且带宽为 0 */
    energy_config_default(&cfg);
    cfg.use_ekf = 0;
    for (size_t i = 0; i < n; i++) {
        s[i].hr = NAN;
        s[i].hrv_log = NAN;
    }
    CHECK(energy_compute(s, n, &cfg, &out) == ENERGY_OK);
    for (size_t i = 0; i < out.n_points; i++) {
        CHECK(isfinite(out.points[i].energy));
        CHECK(out.points[i].e_sd == 0.0);
        CHECK_NEAR(out.points[i].lower, out.points[i].energy, 1e-9);
    }
    energy_output_free(&out);
    free(s);
}

static void test_invalid_args(void)
{
    energy_sample_t s[3];
    memset(s, 0, sizeof(s));
    for (int i = 0; i < 3; i++) {
        s[i].t_h = (double)i / 12.0;
        s[i].tod = s[i].t_h;
        s[i].q_hr = 1.0;
        s[i].q_hrv = 1.0;
        s[i].load = -1.0;
        s[i].activity = -1.0;
    }
    energy_output_t out;
    energy_config_t cfg;
    energy_config_default(&cfg);

    CHECK(energy_compute(NULL, 3, &cfg, &out) == ENERGY_EINVAL);
    CHECK(energy_compute(s, 1, &cfg, &out) == ENERGY_EINVAL);
    s[2].t_h = s[1].t_h; /* 非单调 */
    CHECK(energy_compute(s, 3, &cfg, &out) == ENERGY_EDOMAIN);
}

static void test_phi_clamps(void)
{
    energy_params_t p;
    energy_params_default(&p);
    /* 下限：eta 很小 + s_floor=0 时 raw→0，必须被 phi_min 托住 */
    p.eta = 0.1;
    p.s_floor = 0.0;
    CHECK_NEAR(energy_phi(&p, 0.0), p.phi_min, 1e-12);
    /* 上限：eta 很大时 raw 超过 phi_max */
    p.eta = 5.0;
    CHECK_NEAR(energy_phi(&p, 1.0), p.phi_max, 1e-12);
}

static void test_stream_edges(void)
{
    energy_config_t cfg;
    energy_config_default(&cfg);
    energy_ctx_t ctx;
    energy_ctx_init(&ctx, &cfg);

    energy_sample_t a;
    memset(&a, 0, sizeof(a));
    a.t_h = 0.0;
    a.tod = 0.0;
    a.stage[ENERGY_LIGHT] = 1.0;
    a.hr = 55.0;
    a.hrv_log = log(45.0);
    a.q_hr = 1.0;
    a.q_hrv = 1.0;
    a.load = 0.0;
    a.activity = 0.0;

    energy_point_t p;
    CHECK(energy_ctx_step(&ctx, &a, &p) == ENERGY_OK);

    /* 时间不前进：必须报 EDOMAIN，而不是静默算一步 */
    CHECK(energy_ctx_step(&ctx, &a, &p) == ENERGY_EDOMAIN);

    /* 极端观测 + 大初始协方差：状态仍必须被裁剪到 [0,1] */
    energy_config_default(&cfg);
    cfg.ekf.p0_s = 1.0;
    cfg.ekf.p0_e = 1.0;
    energy_ctx_init(&ctx, &cfg);
    energy_sample_t o = a;
    o.hr = 120.0; /* 夜间 120 bpm 的离群值 */
    CHECK(energy_ctx_step(&ctx, &o, &p) == ENERGY_OK);
    CHECK(ctx.e >= 0.0 && ctx.e <= 1.0);
    CHECK(ctx.s >= 0.0 && ctx.s <= 1.0);
    CHECK(p.energy >= 0.0 && p.energy <= 100.0);
}


/* 用观测模型反推一个“与状态一致”的日间 HR，避免测试让滤波器一路撞边界 */
static double synth_day_hr(const energy_config_t *cfg, const energy_sample_t *s, double e)
{
    const double w = 2.0 * M_PI / 24.0;
    const double t = s->tod;
    const double a = (s->activity >= 0.0) ? s->activity
                                          : energy_activity_step(s->steps, s->kcal, cfg->load.steps_ref);
    return cfg->obs.hr_day[0] + cfg->obs.hr_day[1] * e + cfg->obs.hr_day[2] * a +
           cfg->obs.hr_day[3] * cos(w * t) + cfg->obs.hr_day[4] * sin(w * t) +
           cfg->obs.hr_day[5] * cos(2 * w * t) + cfg->obs.hr_day[6] * sin(2 * w * t);
}

static void test_s_clamp(void)
{
    /* HR 与 S 的交叉协方差会让观测更新把 S 推出 [0,1]，必须裁剪 */
    energy_config_t cfg;
    energy_config_default(&cfg);
    energy_ctx_t ctx;
    energy_ctx_init(&ctx, &cfg);
    ctx.p00 = 100.0;
    ctx.p01 = -100.0;
    ctx.p11 = 100.0;

    energy_sample_t s;
    memset(&s, 0, sizeof(s));
    s.t_h = 0.0;
    s.tod = 12.0;
    s.stage[ENERGY_WAKE] = 1.0;
    s.hr = 30.0; /* 极端低心率 */
    s.hrv_log = log(20.0);
    s.q_hr = 1.0;
    s.q_hrv = 1.0;
    s.load = 0.2;
    s.activity = 0.1;

    energy_point_t p;
    CHECK(energy_ctx_step(&ctx, &s, &p) == ENERGY_OK);
    CHECK(ctx.s >= 0.0 && ctx.s <= 1.0);
    CHECK(ctx.e >= 0.0 && ctx.e <= 1.0);
}

/* 生成一条“自洽”的日间序列：开环真值 → 观测模型反推 HR；可选在一个点上加尖峰。
 * 返回：a/b 两个 ctx 的末态差（尖峰的影响）、a 与真值的跟踪误差。 */
static void run_gate_case(double robust_nu, double *diff_out, double *track_out)
{
    enum { K = 40 };
    energy_config_t cfg;
    energy_config_default(&cfg);
    cfg.ekf.robust_nu = robust_nu;
    /* 日间个人基线已建立（画像），且恰好等于合成数据用的截距：门控/稳健项被真正触发，
     * 而不是因为“还没有常态可比”而跳过更新。 */
    cfg.baseline.init.level[ENERGY_BL_HR_DAY] = cfg.obs.hr_day[0];
    cfg.baseline.init.level[ENERGY_BL_HRV_DAY] = cfg.obs.hrv_day[0];
    cfg.baseline.init.n[ENERGY_BL_HR_DAY] = 5;
    cfg.baseline.init.n[ENERGY_BL_HRV_DAY] = 5;

    energy_sample_t base[K];
    for (int k = 0; k < K; k++) {
        memset(&base[k], 0, sizeof(base[k]));
        base[k].t_h = (double)k / 12.0;
        base[k].tod = 10.0 + base[k].t_h;
        base[k].steps = 45.0;
        base[k].q_hr = 1.0;
        base[k].q_hrv = 1.0;
        base[k].hr = NAN;
        base[k].hrv_log = NAN;
        base[k].load = 0.3;  /* 固定负荷，保证开环真值与 EKF 用同一套动力学 */
        base[k].activity = 0.3;
        base[k].stress = 0.2;
    }

    energy_config_t cfg_ol = cfg;
    cfg_ol.use_ekf = 0;
    energy_ctx_t ol;
    energy_ctx_init(&ol, &cfg_ol);
    double truth[K];
    energy_point_t p;
    for (int k = 0; k < K; k++) {
        energy_ctx_step(&ol, &base[k], &p);
        truth[k] = p.e;
    }
    for (int k = 0; k < K; k++) {
        base[k].hr = synth_day_hr(&cfg, &base[k], truth[k]);
        base[k].hrv_log = 3.40;
    }

    energy_ctx_t a, b;
    energy_ctx_init(&a, &cfg);
    energy_ctx_init(&b, &cfg);
    energy_point_t pa, pb;
    for (int k = 0; k < K; k++) {
        energy_ctx_step(&a, &base[k], &pa);
        energy_sample_t sp = base[k];
        if (k == 30) sp.hr = 220.0;
        energy_ctx_step(&b, &sp, &pb);
    }
    unsigned long na = 0, ga = 0, nb = 0, gb = 0;
    energy_ctx_gate_stats(&a, &na, &ga);
    energy_ctx_gate_stats(&b, &nb, &gb);
    CHECK(na >= (unsigned long)K); /* 观测确实参与了更新（基线已建立） */
    CHECK(ga == 0);                /* 自洽数据不应被门控拒绝 */
    CHECK(gb >= 1);                /* 尖峰被拒绝 */
    *diff_out = fabs(pa.e - pb.e);
    *track_out = fabs(pa.e - truth[K - 1]);
}

static void test_gate_with_gaussian_noise(void)
{
    /* 关掉稳健 t，只留门控：孤立的离群点不应改变状态 */
    double diff = 0.0, track = 0.0;
    run_gate_case(0.0, &diff, &track);
    CHECK(track < 0.05);
    CHECK(diff < 2e-3);
}

static void test_gate_rejects_outlier(void)
{
    /* 稳健 t 打开时，同样不应该被单点离群值拽走 */
    double diff = 0.0, track = 0.0;
    run_gate_case(4.0, &diff, &track);
    CHECK(track < 0.05);
    CHECK(diff < 2e-3);
}

static void day_night_case(const energy_config_t *cfg, double *e_night, double *e_day)
{
    energy_sample_t night, day;
    memset(&night, 0, sizeof(night));
    night.t_h = 0.0;
    night.tod = 2.0;
    night.stage[ENERGY_LIGHT] = 1.0;
    night.hr = 70.0;
    night.hrv_log = 3.60;
    night.q_hr = 1.0;
    night.q_hrv = 1.0;
    night.load = 0.2;
    night.activity = 0.0;

    day = night;
    day.tod = 14.0;
    day.stage[ENERGY_LIGHT] = 0.0;
    day.activity = 0.3;

    energy_ctx_t a, b;
    energy_ctx_init(&a, cfg);
    energy_ctx_init(&b, cfg);
    energy_point_t pn, pd;
    CHECK(energy_ctx_step(&a, &night, &pn) == ENERGY_OK);
    CHECK(energy_ctx_step(&b, &day, &pd) == ENERGY_OK);
    *e_night = pn.e;
    *e_day = pd.e;
}

static void test_day_night_regimes_differ(void)
{
    /* 同样的 HR=70：夜间偏高（能量被下调），日间偏低（能量被上调）。
     * 这条断言专门防止“夜间误用日间观测系数”。
     * v1.2：截距来自个人基线。这里给一个画像基线（夜间 52 bpm、日间 78 bpm），
     * 验证 regime 选择；再用 baseline.enable=0（群体截距）复核同一结论。 */
    energy_config_t cfg;
    energy_config_default(&cfg);
    cfg.baseline.init.level[ENERGY_BL_HR_NIGHT] = 52.0;
    cfg.baseline.init.level[ENERGY_BL_HRV_NIGHT] = 3.60;
    cfg.baseline.init.level[ENERGY_BL_HR_DAY] = 78.0;
    cfg.baseline.init.level[ENERGY_BL_HRV_DAY] = 3.60;
    for (int k = 0; k < 4; k++) cfg.baseline.init.n[k] = 5;
    double en = 0.0, ed = 0.0;
    day_night_case(&cfg, &en, &ed);
    CHECK(en < 0.55); /* 夜间 70 bpm 远高于你的夜间常态 → 能量下调 */
    CHECK(ed > 0.55); /* 日间 70 bpm 低于你的日间常态 → 能量上调 */

    energy_config_default(&cfg);
    cfg.baseline.enable = 0;
    day_night_case(&cfg, &en, &ed);
    CHECK(en < 0.55);
    CHECK(ed > 0.55);

    /* 默认配置、没有任何基线：第一点没有“常态”可比，不做校正（E 保持先验） */
    energy_config_default(&cfg);
    day_night_case(&cfg, &en, &ed);
    CHECK_NEAR(en, cfg.x0_e, 1e-12);
    CHECK_NEAR(ed, cfg.x0_e, 1e-12);
}

/* 回归锁：同一面板在 C v1.2 下的输出（v1.1.1 及以前这里锁的是 Python 端的逐位结果；
 * v1.2 引入个体基线、新的日间系数与缺口策略后 C 端领先于 Python，改为锁 C 自身）。
 * 容差 1e-9；用来锁住 ODE/EKF/基线学习/展示层的实现细节。改算法时需要有意识地更新。 */
static void test_golden_regression(void)
{
    size_t n = 0;
    energy_sample_t *s = make_panel(3, &n);
    CHECK(s != NULL);
    energy_config_t cfg;
    energy_config_default(&cfg);
    energy_output_t out;
    CHECK(energy_compute(s, n, &cfg, &out) == ENERGY_OK);
    CHECK(out.n_points == n);

    const size_t idx[5] = {0, 100, 288, 500, 863};
    /* idx 0：还没有任何个人基线 → 不做观测校正，E/方差保持先验 (0.55, 0.22) */
    const double gs[5] = {0.450000000000000, 0.178535253535194, 0.628505190475809,
                          0.689610699937653, 0.649378555741247};
    const double ge[5] = {0.550000000000000, 0.642573358902543, 0.477788886843377,
                          0.207021529770369, 0.505564321131611};
    const double gsd[5] = {0.220000000000000, 0.046412602267840, 0.037905353515403,
                           0.045809858518258, 0.034360090663948};
    const double gen[5] = {50.674498320876395, 70.509072821767603, 42.684251171374882,
                           19.679965735637843, 44.201721705139072};
    for (int k = 0; k < 5; k++) {
        const energy_point_t *p = &out.points[idx[k]];
        CHECK_NEAR(p->s, gs[k], 1e-9);
        CHECK_NEAR(p->e, ge[k], 1e-9);
        CHECK_NEAR(p->e_sd, gsd[k], 1e-9);
        CHECK_NEAR(p->energy, gen[k], 1e-9);
    }
    CHECK(out.stats_min >= 14.5 && out.stats_min <= 14.8);
    CHECK(out.stats_max >= 86.9 && out.stats_max <= 87.2);
    CHECK((out.points[0].flags & ENERGY_FLAG_LEARNING) != 0);
    CHECK(out.points[863].flags == 0);
    energy_output_free(&out);
    free(s);
}

static void test_misc_api(void)
{
    size_t n = 0;
    energy_sample_t *s = make_panel(1, &n);
    CHECK(s != NULL);

    /* cfg 传 NULL：使用默认配置 */
    energy_output_t out;
    CHECK(energy_compute(s, n, NULL, &out) == ENERGY_OK);
    CHECK(out.n_points == n);
    energy_output_free(&out);

    /* 负荷标定：非法参数 */
    energy_load_stats_t st;
    energy_load_cfg_t lc;
    energy_load_default(&lc);
    CHECK(energy_load_calibrate(NULL, n, &lc, &st) == ENERGY_EINVAL);
    CHECK(energy_load_calibrate(s, 0, &lc, &st) == ENERGY_EINVAL);
    CHECK(energy_load_calibrate(s, n, &lc, NULL) == ENERGY_EINVAL);
    CHECK(energy_load_calibrate(s, n, &lc, &st) == ENERGY_OK);
    CHECK(st.hr_rest > 30.0 && st.hr_rest < 120.0);
    CHECK(st.hr_max >= 150.0 && st.hr_max <= 210.0);
    CHECK(st.l_ref > 0.0);

    /* 流式 peek：不推进时间，状态不变 */
    energy_config_t cfg;
    energy_config_default(&cfg);
    energy_ctx_t ctx;
    energy_ctx_init(&ctx, &cfg);
    energy_point_t p;
    CHECK(energy_ctx_step(&ctx, &s[0], &p) == ENERGY_OK);
    const double e_before = ctx.e;
    energy_point_t q;
    energy_ctx_peek(&ctx, 12.0, &q);
    CHECK_NEAR(q.e, e_before, 1e-15);
    CHECK_NEAR(q.tod, 12.0, 1e-15);
    CHECK(q.energy >= 0.0 && q.energy <= 100.0);

    /* 空指针保护 */
    CHECK(energy_ctx_step(NULL, &s[0], &p) == ENERGY_EINVAL);
    CHECK(energy_ctx_step(&ctx, NULL, &p) == ENERGY_EINVAL);

    free(s);
}

static void test_sample_init_and_auto_activity(void)
{
    /* 初始化：q 默认 1、hr/hrv 缺失、load/activity 交给库 */
    energy_sample_t s0;
    memset(&s0, 0xAA, sizeof(s0)); /* 脏内存 */
    energy_sample_init(&s0);
    CHECK(s0.q_hr == 1.0 && s0.q_hrv == 1.0);
    CHECK(isnan(s0.hr) && isnan(s0.hrv_log));
    CHECK(s0.load < 0.0 && s0.activity < 0.0);

    /* 流式：activity=-1（自动）与显式给同一个值，必须逐点一致。
     * 这曾经是个 bug：流式只把自动活动量写进输出，没写回观测方程。 */
    size_t n = 0;
    energy_sample_t *panel = make_panel(2, &n);
    CHECK(panel != NULL);
    energy_config_t cfg;
    energy_config_default(&cfg);
    energy_ctx_t a, b;
    energy_ctx_init(&a, &cfg);
    energy_ctx_init(&b, &cfg);
    energy_point_t pa, pb;
    for (size_t i = 0; i < n; i++) {
        energy_sample_t sa = panel[i];
        energy_sample_t sb = panel[i];
        sa.activity = -1.0;
        sb.activity = energy_activity_step(panel[i].steps, panel[i].kcal, cfg.load.steps_ref);
        CHECK(energy_ctx_step(&a, &sa, &pa) == ENERGY_OK);
        CHECK(energy_ctx_step(&b, &sb, &pb) == ENERGY_OK);
        CHECK_NEAR(pa.e, pb.e, 1e-12);
        CHECK_NEAR(pa.energy, pb.energy, 1e-12);
    }
    free(panel);
}

static void test_batch_honors_provided_load(void)
{
    /* 批量接口必须尊重调用方传入的 load/activity（>=0 直接用） */
    size_t n = 0;
    energy_sample_t *low = make_panel(1, &n);
    energy_sample_t *high = make_panel(1, &n);
    energy_sample_t *auto_panel = make_panel(1, &n);
    CHECK(low && high && auto_panel);
    for (size_t i = 0; i < n; i++) {
        low[i].load = 0.05;
        high[i].load = 0.95;
        low[i].activity = 0.0;
        high[i].activity = 0.0;
        auto_panel[i].load = -1.0;
        auto_panel[i].activity = -1.0;
    }
    energy_config_t cfg;
    energy_config_default(&cfg);
    energy_output_t ol, oh, oa;
    CHECK(energy_compute(low, n, &cfg, &ol) == ENERGY_OK);
    CHECK(energy_compute(high, n, &cfg, &oh) == ENERGY_OK);
    CHECK(energy_compute(auto_panel, n, &cfg, &oa) == ENERGY_OK);

    /* 传 0.05 与 0.95 必须得到不同的曲线，且低负荷能量更高 */
    CHECK(ol.stats_mean > oh.stats_mean + 5.0);
    CHECK_NEAR(ol.points[100].load, 0.05, 1e-12);
    CHECK_NEAR(oh.points[100].load, 0.95, 1e-12);
    CHECK(fabs(oa.stats_mean - oh.stats_mean) > 1e-9); /* 自动负荷与显式 0.95 不同 */

    /* 显式 activity 也必须生效：只改 activity，日间观测不同 → 曲线不同 */
    energy_sample_t *a0 = make_panel(1, &n);
    energy_sample_t *a9 = make_panel(1, &n);
    for (size_t i = 0; i < n; i++) {
        a0[i].load = 0.3; a0[i].activity = 0.0;
        a9[i].load = 0.3; a9[i].activity = 0.9;
    }
    energy_output_t o0, o9;
    CHECK(energy_compute(a0, n, &cfg, &o0) == ENERGY_OK);
    CHECK(energy_compute(a9, n, &cfg, &o9) == ENERGY_OK);
    CHECK_NEAR(o0.points[100].activity, 0.0, 1e-12);
    CHECK_NEAR(o9.points[100].activity, 0.9, 1e-12);
    double maxdiff = 0.0;
    for (size_t i = 0; i < n; i++) {
        const double d = fabs(o0.points[i].e - o9.points[i].e);
        if (d > maxdiff) maxdiff = d;
    }
    CHECK(maxdiff > 1e-3);

    energy_output_free(&ol);
    energy_output_free(&oh);
    energy_output_free(&oa);
    energy_output_free(&o0);
    energy_output_free(&o9);
    free(low); free(high); free(auto_panel); free(a0); free(a9);
}

static void test_hrv_conversion(void)
{
    /* 毫秒 → ln(ms)：必须是自然对数，且拒绝非正/非有限 */
    CHECK_NEAR(energy_hrv_log_from_ms(45.0), log(45.0), 1e-12);
    CHECK_NEAR(energy_hrv_log_from_ms(1.0), 0.0, 1e-12);
    CHECK(isnan(energy_hrv_log_from_ms(0.0)));
    CHECK(isnan(energy_hrv_log_from_ms(-3.0)));
    CHECK(isnan(energy_hrv_log_from_ms(NAN)));
    /* 与直接取对数等价，且落在观测模型的合理区间内 */
    const double v = energy_hrv_log_from_ms(60.0);
    CHECK(v > 3.0 && v < 4.5);

    /* 端到端：用转换函数填 hrv_log，应当与显式 log() 完全一致 */
    size_t n = 0;
    energy_sample_t *a = make_panel(1, &n);
    energy_sample_t *b = make_panel(1, &n);
    for (size_t i = 0; i < n; i++) {
        const double ms = exp(a[i].hrv_log);
        a[i].hrv_log = energy_hrv_log_from_ms(ms);
        b[i].hrv_log = log(ms);
    }
    energy_config_t cfg;
    energy_config_default(&cfg);
    energy_output_t oa, ob;
    CHECK(energy_compute(a, n, &cfg, &oa) == ENERGY_OK);
    CHECK(energy_compute(b, n, &cfg, &ob) == ENERGY_OK);
    for (size_t i = 0; i < n; i++) CHECK_NEAR(oa.points[i].e, ob.points[i].e, 1e-15);
    energy_output_free(&oa);
    energy_output_free(&ob);
    free(a);
    free(b);
}

static void test_config_validation(void)
{
    energy_config_t cfg;
    energy_config_default(&cfg);
    char err[128];
    err[0] = '\0';
    CHECK(energy_config_validate(&cfg, err, sizeof(err)) == ENERGY_OK);
    CHECK(err[0] == '\0');

    energy_config_t bad = cfg;
    bad.params.kappa_w = 0.0;
    CHECK(energy_config_validate(&bad, err, sizeof(err)) == ENERGY_EINVAL);
    CHECK(strstr(err, "kappa_w") != NULL);

    bad = cfg;
    bad.x0_e = 2.0;
    CHECK(energy_config_validate(&bad, err, sizeof(err)) == ENERGY_EINVAL);

    bad = cfg;
    bad.ekf.q_s = 0.0;
    CHECK(energy_config_validate(&bad, err, sizeof(err)) == ENERGY_EINVAL);

    bad = cfg;
    bad.display.rho = 5.0;
    CHECK(energy_config_validate(&bad, err, sizeof(err)) == ENERGY_EINVAL);

    CHECK(energy_config_validate(NULL, err, sizeof(err)) == ENERGY_EINVAL);

    /* 非法配置：ctx 进入错误态，step 直接返回 EINVAL（不会算出垃圾） */
    bad = cfg;
    bad.params.alpha_sleep = -1.0;
    energy_ctx_t ctx;
    energy_ctx_init(&ctx, &bad);
    CHECK(ctx.config_error == 1);
    energy_sample_t s0;
    memset(&s0, 0, sizeof(s0));
    s0.t_h = 0.0;
    s0.tod = 12.0;
    s0.q_hr = 1.0;
    s0.q_hrv = 1.0;
    s0.load = 0.0;
    s0.activity = 0.0;
    energy_point_t pt;
    CHECK(energy_ctx_step(&ctx, &s0, &pt) == ENERGY_EINVAL);

    /* 批量也要拒绝非法配置 */
    size_t n = 0;
    energy_sample_t *panel = make_panel(1, &n);
    energy_output_t out;
    CHECK(energy_compute(panel, n, &bad, &out) == ENERGY_EINVAL);
    free(panel);
}

static void test_degenerate_inputs(void)
{
    size_t n = 0;
    energy_sample_t *s = make_panel(1, &n);
    CHECK(s != NULL);
    energy_config_t cfg;
    energy_config_default(&cfg);
    energy_output_t out;

    /* tod 缺失/非法：必须能跑，且不能有 UB（ASan/UBSan 下验证） */
    for (size_t i = 0; i < n; i += 7) s[i].tod = NAN;
    for (size_t i = 3; i < n; i += 11) s[i].tod = -1.0;
    CHECK(energy_compute(s, n, &cfg, &out) == ENERGY_OK);
    for (size_t i = 0; i < out.n_points; i++) {
        CHECK(isfinite(out.points[i].energy));
        CHECK(out.points[i].energy >= 0.0 && out.points[i].energy <= 100.0);
        CHECK(isfinite(out.points[i].tod));
    }
    energy_output_free(&out);

    /* t_h 非有限 / 越界 */
    energy_sample_t *bad = make_panel(1, &n);
    bad[10].t_h = NAN;
    CHECK(energy_compute(bad, n, &cfg, &out) == ENERGY_EINVAL);
    bad[10].t_h = 2.0e6;
    CHECK(energy_compute(bad, n, &cfg, &out) == ENERGY_EINVAL);
    free(bad);

    /* 标定接口直接吃用户输入：NaN tod 也必须安全（不经过批次清洗） */
    energy_sample_t *rawp = make_panel(1, &n);
    for (size_t i = 0; i < n; i += 5) rawp[i].tod = NAN;
    energy_load_stats_t st;
    CHECK(energy_load_calibrate(rawp, n, &cfg.load, &st) == ENERGY_OK);
    CHECK(isfinite(st.hr_rest) && isfinite(st.l_ref) && st.l_ref > 0.0);
    free(rawp);

    /* q 超过 1：应被裁剪，结果与 q=1 一致 */
    energy_sample_t *a = make_panel(1, &n);
    energy_sample_t *b = make_panel(1, &n);
    for (size_t i = 0; i < n; i++) b[i].q_hr = 5.0;
    energy_output_t oa, ob;
    CHECK(energy_compute(a, n, &cfg, &oa) == ENERGY_OK);
    CHECK(energy_compute(b, n, &cfg, &ob) == ENERGY_OK);
    for (size_t i = 0; i < n; i++) {
        CHECK_NEAR(oa.points[i].e, ob.points[i].e, 1e-12);
    }
    energy_output_free(&oa);
    energy_output_free(&ob);
    free(a);
    free(b);

    /* 极端平滑分钟数：不应有 int 溢出 */
    cfg.display.smooth_minutes = 1.0e9;
    CHECK(energy_compute(s, n, &cfg, &out) == ENERGY_OK);
    energy_output_free(&out);

    /* 非法起始日期：date 留空，不崩溃 */
    cfg.display.smooth_minutes = 30.0;
    cfg.start_year = 2026;
    cfg.start_month = 13;
    cfg.start_day = 1;
    CHECK(energy_compute(s, n, &cfg, &out) == ENERGY_OK);
    CHECK(out.daily[0].date[0] == '\0');
    energy_output_free(&out);
    free(s);

    /* 最小面板 n=2 */
    energy_sample_t two[2];
    memset(two, 0, sizeof(two));
    for (int i = 0; i < 2; i++) {
        two[i].t_h = (double)i / 12.0;
        two[i].tod = two[i].t_h;
        two[i].hr = 60.0;
        two[i].hrv_log = 3.5;
        two[i].q_hr = 1.0;
        two[i].q_hrv = 1.0;
        two[i].load = -1.0;
        two[i].activity = -1.0;
    }
    CHECK(energy_compute(two, 2, NULL, &out) == ENERGY_OK);
    CHECK(out.n_points == 2);
    CHECK(out.n_daily == 0); /* 一天不到 12 点，不产出日汇总 */
    energy_output_free(&out);
}

static void test_nan_fuzz(void)
{
    /* 随机把字段替换成 NaN/极端值：只要求“不崩溃、输出有界” */
    unsigned seed = 987654321u;
    for (int round = 0; round < 20; round++) {
        size_t n = 0;
        energy_sample_t *s = make_panel(2, &n);
        CHECK(s != NULL);
        for (size_t i = 0; i < n; i++) {
            seed = seed * 1664525u + 1013904223u;
            const unsigned r = seed >> 8;
            if (r % 17 == 0) s[i].hr = NAN;
            if (r % 19 == 0) s[i].hrv_log = NAN;
            if (r % 23 == 0) s[i].stress = NAN;
            if (r % 29 == 0) s[i].steps = NAN;
            if (r % 31 == 0) s[i].q_hr = NAN;
            if (r % 37 == 0) s[i].tod = NAN;
            if (r % 41 == 0) s[i].stage[ENERGY_LIGHT] = NAN;
            if (r % 43 == 0) s[i].load = -1.0;
            if (r % 47 == 0) s[i].hr = 400.0;
        }
        energy_config_t cfg;
        energy_config_default(&cfg);
        energy_output_t out;
        CHECK(energy_compute(s, n, &cfg, &out) == ENERGY_OK);
        for (size_t i = 0; i < out.n_points; i++) {
            CHECK(isfinite(out.points[i].energy));
            CHECK(out.points[i].energy >= 0.0 && out.points[i].energy <= 100.0);
            CHECK(out.points[i].s >= 0.0 && out.points[i].s <= 1.0);
            CHECK(out.points[i].e >= 0.0 && out.points[i].e <= 1.0);
        }
        energy_output_free(&out);
        free(s);
    }
}

static void test_fuzz_bounds(void)
{
    /* 80 组随机参数 × 72 小时：状态解析有界 */
    unsigned seed = 12345u;
    for (int c = 0; c < 80; c++) {
        energy_params_t p;
        energy_params_default(&p);
        seed = seed * 1664525u + 1013904223u;
        p.kappa_w = 0.03 + 0.47 * ((seed >> 8) & 0xFFFF) / 65535.0;
        seed = seed * 1664525u + 1013904223u;
        p.alpha_sleep = 0.05 + 1.45 * ((seed >> 8) & 0xFFFF) / 65535.0;
        seed = seed * 1664525u + 1013904223u;
        p.beta_L = 0.9 * ((seed >> 8) & 0xFFFF) / 65535.0;
        seed = seed * 1664525u + 1013904223u;
        p.eta = 0.1 + 1.4 * ((seed >> 8) & 0xFFFF) / 65535.0;

        double s = 0.5, e = 0.5;
        for (int k = 0; k < 72 * 12; k++) {
            const double m = (k % 24 < 8) ? 1.0 : 0.0;
            const double lam = m * 0.26 + (1.0 - m) * 0.05;
            const double psi = m * 1.0 + (1.0 - m) * 0.15;
            energy_advance(&p, s, e, m, lam, psi, 0.4, 0.3, 5.0 / 60.0, &s, &e);
            if (!(s >= 0.0 && s <= 1.0 && e >= 0.0 && e <= 1.0)) {
                printf("FAIL fuzz case %d step %d s=%.6f e=%.6f\n", c, k, s, e);
                g_fail++;
                g_checks++;
                return;
            }
        }
        g_checks++;
    }
}

/* ------------------------------------------------------------------ */
/* v1.2：个体化观测模型 / 缺口 / 压力量纲 / 日历日                        */
/* ------------------------------------------------------------------ */
static unsigned long long g_rng = 1;

static double t_urand(void)
{
    g_rng = g_rng * 6364136223846793005ULL + 1442695040888963407ULL;
    return (double)(g_rng >> 11) / 9007199254740992.0;
}

static double t_nrand(void)
{
    double u1 = t_urand();
    const double u2 = t_urand();
    if (u1 < 1e-12) u1 = 1e-12;
    return sqrt(-2.0 * log(u1)) * cos(2.0 * M_PI * u2);
}

/* 更贴近真实的一段数据：23:00–07:00 睡眠（90 分钟周期的分期），白天有零星短走
 * （450–550 步/5 min）、隔天傍晚一次 45 分钟运动；hr_shift 把整条心率曲线平移
 * （模拟“体质不同、但恢复模式相同”的两个人）。 */
static energy_sample_t *gen_days(int days, double start_tod, double hr_shift,
                                 unsigned long long seed, size_t *n_out)
{
    g_rng = seed;
    const size_t n = (size_t)days * 288;
    energy_sample_t *s = (energy_sample_t *)calloc(n, sizeof(energy_sample_t));
    if (!s) return NULL;
    for (size_t i = 0; i < n; i++) {
        energy_sample_init(&s[i]);
        const double t = (double)i / 12.0;
        const double tod = fmod(start_tod + t, 24.0);
        const int asleep = (tod >= 23.0 || tod < 7.0);
        const int day = (int)floor((t + start_tod) / 24.0);
        s[i].t_h = t;
        s[i].tod = tod;
        s[i].stress = 0.25 + 0.1 * sin(tod / 24.0 * 2.0 * M_PI);
        if (asleep) {
            const double since = tod >= 23.0 ? tod - 23.0 : tod + 1.0;
            const int ph = (int)fmod(since * 60.0, 90.0) / 23;
            const int st = ph == 0 ? ENERGY_LIGHT : ph == 1 ? ENERGY_DEEP : ph == 2 ? ENERGY_LIGHT : ENERGY_REM;
            s[i].stage[st] = 1.0;
            s[i].hr = 54.0 + hr_shift - (st == ENERGY_DEEP ? 3.0 : 0.0) + (st == ENERGY_REM ? 2.0 : 0.0) -
                      2.0 * since / 8.0 + 1.5 * t_nrand();
            s[i].hrv_log = log(50.0) + (st == ENERGY_DEEP ? 0.2 : 0.0) + 0.15 * t_nrand();
        } else {
            double steps = fmax(0.0, 30.0 + 30.0 * t_nrand());
            if (fmod(tod, 3.0) < 0.25) steps = 450.0 + 100.0 * t_urand();
            const int exercise = (day % 2 == 0) && tod >= 18.0 && tod < 18.75;
            if (exercise) steps = 650.0;
            const double a = fmin(steps / 500.0, 1.5);
            s[i].steps = steps;
            s[i].hr = 70.0 + hr_shift + 30.0 * a + (exercise ? 55.0 : 0.0) + 4.0 * t_nrand();
        }
    }
    *n_out = n;
    return s;
}

static const energy_point_t *point_at(const energy_point_t *p, size_t n, double t_h)
{
    for (size_t i = 0; i < n; i++) {
        if (fabs(p[i].t_h - t_h) < 1e-6) return &p[i];
    }
    return NULL;
}

static double e_at(const energy_output_t *o, double t_h)
{
    const energy_point_t *p = point_at(o->points, o->n_points, t_h);
    return p ? p->e : NAN;
}

static void stream_all(const energy_sample_t *s, size_t n, const energy_config_t *cfg,
                       energy_point_t *pts, energy_ctx_t *ctx)
{
    energy_ctx_init(ctx, cfg);
    for (size_t i = 0; i < n; i++) {
        CHECK(energy_ctx_step(ctx, &s[i], pts ? &pts[i] : NULL) == ENERGY_OK);
    }
}

static void test_constitution_invariance(void)
{
    /* 核心原则：能量只跟自己比。两个人恢复模式相同、只是静息心率不同（整条曲线平移），
     * 能量必须相同；v1.1 的群体截距（58 bpm）会把体质差异直接变成能量差异。 */
    const double shifts[4] = {-6.0, 0.0, 6.0, 12.0};
    double mean_bl[4] = {0}, mean_pop[4] = {0}, lvl[4] = {0};
    for (int pop = 0; pop < 2; pop++) {
        energy_config_t cfg;
        energy_config_default(&cfg);
        cfg.baseline.enable = pop ? 0 : 1;
        for (int k = 0; k < 4; k++) {
            size_t n = 0;
            energy_sample_t *s = gen_days(10, 0.0, shifts[k], 777u, &n);
            CHECK(s != NULL);
            if (!s) return;
            energy_output_t out;
            CHECK(energy_compute(s, n, &cfg, &out) == ENERGY_OK);
            double m = 0.0;
            for (size_t i = 288 * 3; i < n; i++) m += out.points[i].e;
            m /= (double)(n - 288 * 3);
            if (pop) {
                mean_pop[k] = m;
            } else {
                mean_bl[k] = m;
                lvl[k] = out.baseline.level[ENERGY_BL_HR_NIGHT];
                CHECK(out.baseline.n[ENERGY_BL_HR_NIGHT] >= 9);
                CHECK(out.baseline.n[ENERGY_BL_HR_DAY] >= 9);
            }
            energy_output_free(&out);
            free(s);
        }
    }
    for (int k = 0; k < 4; k++) {
        CHECK_NEAR(mean_bl[k], mean_bl[1], 0.01);
        /* 平移量全部进了个人基线，而不是进了能量 */
        CHECK_NEAR(lvl[k] - lvl[1], shifts[k], 0.05);
    }
    /* 对照：关掉个体基线时，同样的“平移”造成巨大的能量差（这正是 v1.2 修掉的问题） */
    CHECK(mean_pop[0] - mean_pop[3] > 0.2);
}

static void test_gate_rate_individual(void)
{
    /* 静息心率比群体先验高 10 bpm 的人：观测不应被门控大面积拒绝（v1.1 日间拒绝率约 90%） */
    size_t n = 0;
    energy_sample_t *s = gen_days(10, 0.0, 10.0, 4242u, &n);
    CHECK(s != NULL);
    if (!s) return;
    energy_config_t cfg;
    energy_config_default(&cfg);
    energy_ctx_t ctx;
    stream_all(s, n, &cfg, NULL, &ctx);
    unsigned long n_obs = 0, n_gated = 0;
    energy_ctx_gate_stats(&ctx, &n_obs, &n_gated);
    CHECK(n_obs > 2000);
    CHECK((double)n_gated < 0.05 * (double)n_obs);
    free(s);
}

static void test_daytime_not_saturated(void)
{
    /* 日间观测方程：普通作息 + 零星快走 + 隔天运动，日间 E 不应被推到 1.0 */
    size_t n = 0;
    energy_sample_t *s = gen_days(14, 0.0, 0.0, 99u, &n);
    CHECK(s != NULL);
    if (!s) return;
    energy_config_t cfg;
    energy_config_default(&cfg);
    energy_output_t out;
    CHECK(energy_compute(s, n, &cfg, &out) == ENERGY_OK);
    int day_n = 0, hi = 0;
    double day_mean = 0.0;
    for (size_t i = 288; i < n; i++) {
        if (s[i].tod >= 7.0 && s[i].tod < 23.0) {
            day_n++;
            day_mean += out.points[i].e;
            if (out.points[i].e > 0.98) hi++;
        }
    }
    day_mean /= (double)day_n;
    CHECK(hi < day_n / 100);
    CHECK(day_mean > 0.3 && day_mean < 0.9);
    /* 白天应该有明显的下降：傍晚低于早上 */
    for (int d = 3; d < 14; d++) {
        CHECK(e_at(&out, d * 24 + 22.0) < e_at(&out, d * 24 + 8.0));
    }
    energy_output_free(&out);
    free(s);
}

static void test_individual_deviation(void)
{
    /* 只跟自己比：同一个人某晚心率 +6、HRV −0.3（生病/饮酒），次日早晨能量明显下降，
     * 再下一天恢复；这个坏夜被降权，不会把基线拉走。 */
    size_t n = 0;
    energy_sample_t *s0 = gen_days(14, 0.0, 0.0, 31u, &n);
    energy_sample_t *s1 = gen_days(14, 0.0, 0.0, 31u, &n);
    CHECK(s0 && s1);
    if (!s0 || !s1) return;
    for (size_t i = 0; i < n; i++) {
        if (s1[i].t_h >= 9 * 24 - 1 && s1[i].t_h < 9 * 24 + 7) {
            s1[i].hr += 6.0;
            s1[i].hrv_log -= 0.3;
        }
    }
    energy_config_t cfg;
    energy_config_default(&cfg);
    energy_output_t o0, o1;
    CHECK(energy_compute(s0, n, &cfg, &o0) == ENERGY_OK);
    CHECK(energy_compute(s1, n, &cfg, &o1) == ENERGY_OK);
    const double t_bad = 9 * 24 + 7.5, t_next = 10 * 24 + 7.5;
    CHECK(e_at(&o0, t_bad) - e_at(&o1, t_bad) > 0.15);
    CHECK(fabs(e_at(&o0, t_next) - e_at(&o1, t_next)) < 0.1);
    CHECK(fabs(o0.baseline.level[ENERGY_BL_HR_NIGHT] - o1.baseline.level[ENERGY_BL_HR_NIGHT]) < 1.0);
    energy_output_free(&o0);
    energy_output_free(&o1);
    free(s0);
    free(s1);
}

static void test_baseline_learning_and_api(void)
{
    size_t n = 0;
    energy_sample_t *s = gen_days(5, 0.0, 0.0, 21u, &n);
    CHECK(s != NULL);
    if (!s) return;
    energy_config_t cfg;
    energy_config_default(&cfg);
    energy_output_t out;
    CHECK(energy_compute(s, n, &cfg, &out) == ENERGY_OK);

    /* 学到的基线与作息 */
    const energy_baseline_t *bl = &out.baseline;
    CHECK(bl->n[ENERGY_BL_HR_NIGHT] >= 4 && bl->n[ENERGY_BL_HRV_NIGHT] >= 4);
    CHECK(bl->n[ENERGY_BL_HR_DAY] >= 4);
    CHECK(bl->n[ENERGY_BL_HRV_DAY] == 0); /* 日间没有 HRV：该通道保持“未建立” */
    CHECK(isfinite(bl->level[ENERGY_BL_HR_NIGHT]) && isfinite(bl->level[ENERGY_BL_HRV_NIGHT]));
    CHECK(bl->n_sleep >= 4);
    CHECK(fabs(bl->sleep_on - 23.0) < 0.5);
    CHECK(fabs(bl->sleep_off - 7.0) < 0.5);
    CHECK(bl->stage_mix[ENERGY_DEEP] > 0.1 && bl->stage_mix[ENERGY_LIGHT] > 0.3);

    /* 冷启动标志：前两晚提示“正在了解你的常态”，第 3 晚之后消失 */
    CHECK((out.points[90].flags & ENERGY_FLAG_LEARNING) != 0);       /* day0 07:30 */
    CHECK((out.points[4 * 288 + 90].flags & ENERGY_FLAG_LEARNING) == 0); /* day4 07:30 */

    /* 画像回灌：带着已学基线重新计算，从第一点起就不是冷启动 */
    energy_config_t warm = cfg;
    warm.baseline.init = out.baseline;
    char err[128];
    CHECK(energy_config_validate(&warm, err, sizeof(err)) == ENERGY_OK);
    energy_output_t out2;
    CHECK(energy_compute(s, n, &warm, &out2) == ENERGY_OK);
    CHECK((out2.points[0].flags & ENERGY_FLAG_LEARNING) == 0);
    CHECK(out2.baseline.n[ENERGY_BL_HR_NIGHT] > out.baseline.n[ENERGY_BL_HR_NIGHT]);
    energy_output_free(&out2);

    /* 流式 get/set 往返 */
    energy_ctx_t ctx;
    energy_ctx_init(&ctx, &cfg);
    energy_baseline_t got;
    energy_ctx_get_baseline(&ctx, &got);
    CHECK(got.n[0] == 0 && got.n_sleep == 0);
    CHECK(energy_ctx_set_baseline(&ctx, &out.baseline) == ENERGY_OK);
    energy_ctx_get_baseline(&ctx, &got);
    for (int k = 0; k < 4; k++) {
        CHECK(got.n[k] == out.baseline.n[k]);
        if (got.n[k] > 0) CHECK_NEAR(got.level[k], out.baseline.level[k], 0.0);
    }
    CHECK_NEAR(got.sleep_on, out.baseline.sleep_on, 0.0);
    energy_point_t pk;
    energy_ctx_peek(&ctx, 8.0, &pk);
    CHECK((pk.flags & ENERGY_FLAG_LEARNING) == 0);

    /* 非法基线：拒绝且不改动 ctx */
    energy_baseline_t bad = out.baseline;
    bad.n[0] = -1;
    CHECK(energy_ctx_set_baseline(&ctx, &bad) == ENERGY_EINVAL);
    bad = out.baseline;
    bad.level[ENERGY_BL_HR_NIGHT] = NAN;
    CHECK(energy_ctx_set_baseline(&ctx, &bad) == ENERGY_EINVAL);
    bad = out.baseline;
    bad.sleep_on = 25.0;
    CHECK(energy_ctx_set_baseline(&ctx, &bad) == ENERGY_EINVAL);
    CHECK(energy_ctx_set_baseline(&ctx, NULL) == ENERGY_EINVAL);
    CHECK(energy_ctx_set_baseline(NULL, &out.baseline) == ENERGY_EINVAL);
    energy_ctx_get_baseline(&ctx, &got);
    CHECK(got.n[0] == out.baseline.n[0]);
    warm.baseline.init = bad;
    CHECK(energy_config_validate(&warm, err, sizeof(err)) == ENERGY_EINVAL);

    energy_output_free(&out);
    free(s);
}

static void test_stream_auto_load_matches_batch(void)
{
    /* 不填画像（load 全自动）：流式在线标定应与批量窗口标定接近（v1.1：0.944 vs 0.124） */
    size_t n = 0;
    energy_sample_t *s = gen_days(7, 0.0, 0.0, 5u, &n);
    energy_point_t *p = (energy_point_t *)calloc(n, sizeof(energy_point_t));
    CHECK(s && p);
    if (!s || !p) {
        free(s);
        free(p);
        return;
    }
    energy_config_t cfg;
    energy_config_default(&cfg);
    energy_output_t out;
    CHECK(energy_compute(s, n, &cfg, &out) == ENERGY_OK);
    energy_ctx_t ctx;
    stream_all(s, n, &cfg, p, &ctx);
    double max_de = 0.0, sed_b = 0.0, sed_s = 0.0;
    int sed_n = 0;
    for (size_t i = 0; i < n; i++) {
        if (i >= 288) max_de = fmax(max_de, fabs(out.points[i].e - p[i].e));
        if (s[i].tod >= 9.0 && s[i].tod < 17.0 && s[i].steps < 100.0) {
            sed_b += out.points[i].load;
            sed_s += p[i].load;
            sed_n++;
        }
    }
    CHECK(max_de < 0.06);
    CHECK(sed_s / sed_n < 0.15); /* 久坐时流式负荷不再顶到 1 */
    CHECK(fabs(sed_s / sed_n - sed_b / sed_n) < 0.06);

    /* 用标定结果填画像：流式与批量逐点一致 */
    energy_load_stats_t st;
    CHECK(energy_load_calibrate(s, n, &cfg.load, &st) == ENERGY_OK);
    CHECK(st.hr_base[3] > 0.0 && st.hr_base[14] > st.hr_base[3]);
    energy_config_t prof = cfg;
    prof.load.hr_rest = st.hr_rest;
    prof.load.hr_max = st.hr_max;
    prof.load.l_ref = st.l_ref;
    memcpy(prof.load.hr_base, st.hr_base, sizeof(st.hr_base));
    energy_output_free(&out);
    CHECK(energy_compute(s, n, &prof, &out) == ENERGY_OK);
    stream_all(s, n, &prof, p, &ctx);
    for (size_t i = 0; i < n; i++) {
        CHECK_NEAR(p[i].load, out.points[i].load, 1e-12);
        CHECK_NEAR(p[i].e, out.points[i].e, 1e-12);
    }
    energy_output_free(&out);
    free(p);
    free(s);
}

static void test_stream_load_after_gap(void)
{
    /* 数据中断 3 小时后的第一点：TRIMP 不能按 3 小时累计（v1.1 直接顶到 L=1） */
    energy_config_t cfg;
    energy_config_default(&cfg);
    energy_ctx_t ctx;
    energy_ctx_init(&ctx, &cfg);
    energy_sample_t a;
    energy_sample_init(&a);
    a.t_h = 0.0;
    a.tod = 9.0;
    a.hr = 70.0;
    a.steps = 10.0;
    energy_point_t p, q;
    CHECK(energy_ctx_step(&ctx, &a, &p) == ENERGY_OK);
    a.t_h = 1.0 / 12.0;
    a.tod = 9.0 + 1.0 / 12.0;
    CHECK(energy_ctx_step(&ctx, &a, &p) == ENERGY_OK);
    CHECK((p.flags & ENERGY_FLAG_GAP) == 0);
    a.t_h = 3.0;
    a.tod = 12.0;
    CHECK(energy_ctx_step(&ctx, &a, &q) == ENERGY_OK);
    CHECK(q.load < 0.2);
    CHECK((q.flags & ENERGY_FLAG_GAP) != 0);
    CHECK((q.flags & ENERGY_FLAG_GAP_SLEEP) == 0); /* 09:05–12:00 不在睡眠窗内 */
}

static void test_nowear_night(void)
{
    /* 第 3 晚没戴戒指（8 小时）：次日早晨不能按“整夜清醒”推算（v1.1：E=0.094）。
     * 缺行（没有记录）与空行（有时间戳、全缺失）结果应一致。 */
    size_t n = 0;
    energy_sample_t *s = gen_days(5, 0.0, 0.0, 11u, &n);
    energy_sample_t *absent = (energy_sample_t *)calloc(n, sizeof(energy_sample_t));
    energy_sample_t *empty = (energy_sample_t *)calloc(n, sizeof(energy_sample_t));
    energy_point_t *p = (energy_point_t *)calloc(n, sizeof(energy_point_t));
    CHECK(s && absent && empty && p);
    if (!s || !absent || !empty || !p) {
        free(s); free(absent); free(empty); free(p);
        return;
    }
    size_t m = 0;
    for (size_t i = 0; i < n; i++) {
        empty[i] = s[i];
        if (s[i].t_h >= 71.0 && s[i].t_h < 79.0) {
            energy_sample_init(&empty[i]);
            empty[i].t_h = s[i].t_h;
            empty[i].tod = s[i].tod;
        } else {
            absent[m++] = s[i];
        }
    }
    energy_config_t cfg;
    energy_config_default(&cfg);
    energy_output_t worn, oa, oe;
    CHECK(energy_compute(s, n, &cfg, &worn) == ENERGY_OK);
    CHECK(energy_compute(absent, m, &cfg, &oa) == ENERGY_OK);
    CHECK(energy_compute(empty, n, &cfg, &oe) == ENERGY_OK);
    const double t_m = 79.5;
    CHECK(e_at(&oa, t_m) > 0.6);
    CHECK(fabs(e_at(&oa, t_m) - e_at(&worn, t_m)) < 0.1);
    CHECK(fabs(e_at(&oe, t_m) - e_at(&oa, t_m)) < 0.02);
    const energy_point_t *first = point_at(oa.points, oa.n_points, 79.0);
    CHECK(first != NULL);
    if (first) {
        CHECK((first->flags & ENERGY_FLAG_GAP) != 0);
        CHECK((first->flags & ENERGY_FLAG_GAP_SLEEP) != 0);
    }
    /* 流式（缺行）与批量一致的结论 */
    energy_ctx_t ctx;
    stream_all(absent, m, &cfg, p, &ctx);
    const energy_point_t *ps = point_at(p, m, t_m);
    CHECK(ps != NULL);
    if (ps) CHECK(fabs(ps->e - e_at(&oa, t_m)) < 0.05);

    /* 有效数据中断超过 6 小时也不能把这段时间当成白天学基线：基线仍然稳定 */
    CHECK(fabs(oa.baseline.level[ENERGY_BL_HR_NIGHT] - worn.baseline.level[ENERGY_BL_HR_NIGHT]) < 1.0);

    energy_output_free(&worn);
    energy_output_free(&oa);
    energy_output_free(&oe);
    free(s); free(absent); free(empty); free(p);
}

static void test_stress_scale(void)
{
    /* 压力量纲按整段判断：0–100 设备上的 1.2 分是“几乎没有压力”，不是满分 */
    size_t n = 0;
    energy_sample_t *a = gen_days(2, 0.0, 0.0, 3u, &n);
    energy_sample_t *b = gen_days(2, 0.0, 0.0, 3u, &n);
    energy_sample_t *c = gen_days(2, 0.0, 0.0, 3u, &n);
    CHECK(a && b && c);
    if (!a || !b || !c) {
        free(a); free(b); free(c);
        return;
    }
    for (size_t i = 0; i < n; i++) {
        a[i].stress = (i % 50 == 0) ? 60.0 : 1.2;   /* 0–100 量纲 */
        b[i].stress = (i % 50 == 0) ? 0.60 : 0.012; /* 同一序列的 0–1 量纲 */
        c[i].stress = (i % 50 == 0) ? 1.0 : 1.0;    /* 0–1 量纲的满分压力 */
    }
    energy_config_t cfg;
    energy_config_default(&cfg);
    energy_output_t oa, ob, oc, ox;
    CHECK(energy_compute(a, n, &cfg, &oa) == ENERGY_OK);
    CHECK(energy_compute(b, n, &cfg, &ob) == ENERGY_OK);
    CHECK(energy_compute(c, n, &cfg, &oc) == ENERGY_OK);
    energy_config_t fixed = cfg;
    fixed.stress_scale = 100.0;
    CHECK(energy_compute(a, n, &fixed, &ox) == ENERGY_OK);
    for (size_t i = 0; i < n; i++) {
        CHECK_NEAR(oa.points[i].e, ob.points[i].e, 1e-12);
        CHECK_NEAR(ox.points[i].e, ob.points[i].e, 1e-12);
    }
    CHECK(oa.stats_mean > oc.stats_mean + 1.0);

    /* 流式：显式 stress_scale=100 从第一点起一致；自动模式在首个 >1.5 的值之后一致 */
    energy_ctx_t ca, cb, cx;
    energy_ctx_init(&ca, &cfg);
    energy_ctx_init(&cb, &cfg);
    energy_ctx_init(&cx, &fixed);
    energy_point_t pa, pb, px;
    for (size_t i = 0; i < n; i++) {
        CHECK(energy_ctx_step(&ca, &a[i], &pa) == ENERGY_OK); /* a[0] = 60 → 立即切到 /100 */
        CHECK(energy_ctx_step(&cb, &b[i], &pb) == ENERGY_OK);
        CHECK(energy_ctx_step(&cx, &a[i], &px) == ENERGY_OK);
        CHECK_NEAR(pa.e, pb.e, 1e-12);
        CHECK_NEAR(px.e, pb.e, 1e-12);
    }
    energy_output_free(&oa);
    energy_output_free(&ob);
    energy_output_free(&oc);
    energy_output_free(&ox);
    free(a); free(b); free(c);
}

static void test_daily_calendar(void)
{
    /* 每日汇总按当地日历日切分（起点 22:00），日期跨月正确 */
    size_t n = 0;
    energy_sample_t *s = gen_days(3, 22.0, 0.0, 3u, &n);
    CHECK(s != NULL);
    if (!s) return;
    energy_config_t cfg;
    energy_config_default(&cfg);
    cfg.start_year = 2026;
    cfg.start_month = 10;
    cfg.start_day = 31;
    energy_output_t out;
    CHECK(energy_compute(s, n, &cfg, &out) == ENERGY_OK);
    CHECK(out.n_daily == 4);
    if (out.n_daily == 4) {
        CHECK(out.daily[0].n == 24);
        CHECK(out.daily[1].n == 288);
        CHECK(out.daily[2].n == 288);
        CHECK(out.daily[3].n == 264);
        CHECK(strcmp(out.daily[0].date, "10-31") == 0);
        CHECK(strcmp(out.daily[1].date, "11-01") == 0);
        CHECK(strcmp(out.daily[3].date, "11-03") == 0);
        for (size_t d = 0; d < out.n_daily; d++) CHECK(out.daily[d].day == (int)d);
        /* 早晨值取当日 09:00 前的最后一点：日历日 1 的早晨在第 1 晚之后 */
        CHECK(out.daily[1].morning > 0.0);
    }
    energy_output_free(&out);
    free(s);
}

static void test_config_validation_v12(void)
{
    energy_config_t cfg, bad;
    energy_config_default(&cfg);
    char err[160];
    CHECK(energy_config_validate(&cfg, err, sizeof(err)) == ENERGY_OK);
    CHECK(cfg.baseline.enable == 1);
    CHECK(cfg.stress_scale == 0.0);
    CHECK(cfg.ekf.gate_max_reject == 6);
    CHECK_NEAR(cfg.obs.hr_day[2], 45.0, 1e-12);

    bad = cfg; bad.baseline.lambda = 0.0;
    CHECK(energy_config_validate(&bad, err, sizeof(err)) == ENERGY_EINVAL);
    bad = cfg; bad.baseline.lambda = 1.5;
    CHECK(energy_config_validate(&bad, err, sizeof(err)) == ENERGY_EINVAL);
    bad = cfg; bad.baseline.cold_r_inflate = 0.5;
    CHECK(energy_config_validate(&bad, err, sizeof(err)) == ENERGY_EINVAL);
    bad = cfg; bad.baseline.bootstrap_samples = 0;
    CHECK(energy_config_validate(&bad, err, sizeof(err)) == ENERGY_EINVAL);
    bad = cfg; bad.baseline.anomaly[1] = 0.0;
    CHECK(energy_config_validate(&bad, err, sizeof(err)) == ENERGY_EINVAL);
    bad = cfg; bad.baseline.anomaly_weight = 2.0;
    CHECK(energy_config_validate(&bad, err, sizeof(err)) == ENERGY_EINVAL);
    bad = cfg; bad.baseline.init.n[2] = 1; bad.baseline.init.level[2] = INFINITY;
    CHECK(energy_config_validate(&bad, err, sizeof(err)) == ENERGY_EINVAL);
    bad = cfg; bad.gap.max_hold_h = 0.0;
    CHECK(energy_config_validate(&bad, err, sizeof(err)) == ENERGY_EINVAL);
    bad = cfg; bad.gap.sleep_start = 24.0;
    CHECK(energy_config_validate(&bad, err, sizeof(err)) == ENERGY_EINVAL);
    bad = cfg; bad.gap.template_stage[1] = 0.9;
    CHECK(energy_config_validate(&bad, err, sizeof(err)) == ENERGY_EINVAL);
    bad = cfg; bad.stress_scale = -1.0;
    CHECK(energy_config_validate(&bad, err, sizeof(err)) == ENERGY_EINVAL);
    bad = cfg; bad.ekf.gate_max_reject = -1;
    CHECK(energy_config_validate(&bad, err, sizeof(err)) == ENERGY_EINVAL);
    bad = cfg; bad.obs.hr_day[1] = 5.0; /* 恢复越好心率越高：符号错 */
    CHECK(energy_config_validate(&bad, err, sizeof(err)) == ENERGY_EINVAL);
    bad = cfg; bad.obs.hrv_night[3] = NAN;
    CHECK(energy_config_validate(&bad, err, sizeof(err)) == ENERGY_EINVAL);
    bad = cfg; bad.load.hr_base[5] = NAN;
    CHECK(energy_config_validate(&bad, err, sizeof(err)) == ENERGY_EINVAL);
    bad = cfg; bad.load.steps_ref = 0.0;
    CHECK(energy_config_validate(&bad, err, sizeof(err)) == ENERGY_EINVAL);
    CHECK(strcmp(energy_version(), ENERGY_VERSION) == 0);
}

int main(void)
{
    test_defaults();
    test_phi();
    test_advance();
    test_equilibrium();
    test_display();
    test_load();
    test_batch_basic();
    test_stream_equals_batch();
    test_gap_and_openloop();
    test_invalid_args();
    test_misc_api();
    test_phi_clamps();
    test_stream_edges();
    test_s_clamp();
    test_day_night_regimes_differ();
    test_golden_regression();
    test_gate_with_gaussian_noise();
    test_gate_rejects_outlier();
    test_sample_init_and_auto_activity();
    test_batch_honors_provided_load();
    test_hrv_conversion();
    test_config_validation();
    test_degenerate_inputs();
    test_nan_fuzz();
    test_fuzz_bounds();
    test_constitution_invariance();
    test_gate_rate_individual();
    test_daytime_not_saturated();
    test_individual_deviation();
    test_baseline_learning_and_api();
    test_stream_auto_load_matches_batch();
    test_stream_load_after_gap();
    test_nowear_night();
    test_stress_scale();
    test_daily_calendar();
    test_config_validation_v12();

    printf("%s: %d checks, %d failures\n", g_fail ? "FAILED" : "OK", g_checks, g_fail);
    return g_fail ? 1 : 0;
}
