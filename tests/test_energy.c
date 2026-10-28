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
    }
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
    const double a = 0.8 * (s->steps / 120.0);
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

static void test_day_night_regimes_differ(void)
{
    /* 同样的 HR=70：夜间偏高（能量被下调），日间偏低（能量被上调）。
     * 这条断言专门防止“夜间误用日间观测系数”。 */
    energy_config_t cfg;
    energy_config_default(&cfg);

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
    energy_ctx_init(&a, &cfg);
    energy_ctx_init(&b, &cfg);
    energy_point_t pn, pd;
    CHECK(energy_ctx_step(&a, &night, &pn) == ENERGY_OK);
    CHECK(energy_ctx_step(&b, &day, &pd) == ENERGY_OK);

    CHECK(pn.e < 0.55); /* 夜间 70 bpm 高于模型期望 → 能量下调 */
    CHECK(pd.e > 0.55); /* 日间 70 bpm 低于模型期望 → 能量上调 */
}

/* 金标准：同一面板在 Python 端（energy_api.compute_energy，已与仓库逐位对齐）的输出。
 * 容差 1e-9；用来锁住 ODE/EKF/展示层的实现细节（例如 Joseph 项、观测系数）。 */
static void test_golden_against_python(void)
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
    const double gs[5] = {0.450000000000000, 0.129750033966399, 0.612387642326741,
                          0.676002793326633, 0.634173672123537};
    const double ge[5] = {0.381714615132108, 0.308223666180302, 0.395757730884946,
                          0.216905300515907, 0.386758788459087};
    const double gsd[5] = {0.147048097750264, 0.043930002026221, 0.033332047371451,
                           0.045683659765298, 0.033582207388786};
    const double gen[5] = {33.209909616835517, 33.700015596945711, 35.536163209027023,
                           20.591435931981682, 34.117314068115711};
    for (int k = 0; k < 5; k++) {
        const energy_point_t *p = &out.points[idx[k]];
        CHECK_NEAR(p->s, gs[k], 1e-9);
        CHECK_NEAR(p->e, ge[k], 1e-9);
        CHECK_NEAR(p->e_sd, gsd[k], 1e-9);
        CHECK_NEAR(p->energy, gen[k], 1e-9);
    }
    CHECK(out.stats_min >= 17.0 && out.stats_min <= 17.2);
    CHECK(out.stats_max >= 54.8 && out.stats_max <= 55.0);
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
    test_golden_against_python();
    test_gate_with_gaussian_noise();
    test_gate_rejects_outlier();
    test_config_validation();
    test_degenerate_inputs();
    test_nan_fuzz();
    test_fuzz_bounds();

    printf("%s: %d checks, %d failures\n", g_fail ? "FAILED" : "OK", g_checks, g_fail);
    return g_fail ? 1 : 0;
}
