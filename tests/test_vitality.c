/*
 * test_vitality.c —— 精力层单元 / 回归测试（无第三方依赖）
 *
 * 运行：make test
 */
#include "vitality.h"
#include "../src/energy_internal.h"

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
/* 合成面板：23:00–07:00 睡眠，白天 HR 带节律；可选场景                    */
/* ------------------------------------------------------------------ */
typedef struct {
    int nap_day;        /* 该日 13:00 午睡 20 min（分期 light）；<0 关闭 */
    int exercise_day;   /* 该日 18:00 跑步 30 min；<0 关闭 */
    int arousal_day;    /* 该日 09:30–10:15 心率 +10 bpm 无运动；<0 关闭 */
    int stress_day;     /* 该日 11:00–12:00 心率 +22 bpm 无运动；<0 关闭 */
    int short_day;      /* 该日凌晨 0 点才睡、5 点醒（只睡 5 h）；<0 关闭 */
    int broken_day;     /* 该日凌晨 03:00–03:20 夜醒 20 min（硬 wake 标签）；<0 关闭 */
    double hr_offset;   /* 整体心率平移（体质差异） */
    double spo2;        /* >0：夜间每点血氧均值（恒定）；0 = 不给血氧 */
    int spo2_bad_day;   /* 该日醒来前那一夜血氧不稳：均值 −3、15% 的样本再低 4；<0 关闭 */
    int spo2_sparse_day; /* 该日醒来前那一夜只有 30% 样本有血氧；<0 关闭 */
    int spo2_min_day;   /* 该日醒来前那一夜 spo2 均值不变、但 20% 样本 spo2_min 低 4；<0 关闭 */
} scenario_t;

static void scenario_default(scenario_t *sc)
{
    sc->nap_day = sc->exercise_day = sc->arousal_day = sc->stress_day = sc->short_day = -1;
    sc->broken_day = -1;
    sc->hr_offset = 0.0;
    sc->spo2 = 0.0;
    sc->spo2_bad_day = sc->spo2_sparse_day = sc->spo2_min_day = -1;
}

static energy_sample_t *make_panel(size_t days, const scenario_t *sc, size_t *n_out)
{
    const size_t n = 288 * days;
    energy_sample_t *s = (energy_sample_t *)calloc(n, sizeof(energy_sample_t));
    if (!s) return NULL;
    for (size_t i = 0; i < n; i++) {
        energy_sample_t *p = &s[i];
        energy_sample_init(p);
        const double t = (double)i / 12.0;
        const double tod = fmod(t, 24.0);
        const int day = (int)floor(t / 24.0);
        int asleep = (tod >= 23.0 || tod < 7.0);
        if (sc->short_day == day) asleep = (tod >= 0.0 && tod < 5.0);
        if (sc->short_day >= 0 && sc->short_day == day + 1 && tod >= 23.0) asleep = 0; /* 前一晚熬夜到 0 点 */
        if (sc->broken_day == day && tod >= 3.0 && tod < 3.25) asleep = 0;
        p->t_h = t;
        p->tod = tod;
        double hr, hrv, steps = 0.0;
        if (asleep) {
            const double phase = fmod(tod + 1.0, 24.0) / 8.0;
            const double deep = (phase < 0.5) ? 0.45 : 0.15;
            p->stage[ENERGY_DEEP] = deep;
            p->stage[ENERGY_REM] = (phase < 0.5) ? 0.15 : 0.35;
            p->stage[ENERGY_LIGHT] = 1.0 - deep - p->stage[ENERGY_REM] - 0.05;
            p->stage[ENERGY_WAKE] = 0.05;
            hr = 54.0 + 3.0 * phase;
            hrv = log(48.0 + 4.0 * sin(phase * 6.28));
        } else {
            const double w = 2.0 * M_PI / 24.0;
            hr = 70.0 - 2.0 * cos(w * tod) - 3.46 * sin(w * tod) + 0.75 * cos(2.0 * w * tod) -
                 1.30 * sin(2.0 * w * tod);
            hrv = log(34.0);
            steps = 40.0;
            if (sc->arousal_day == day && tod >= 9.5 && tod < 10.25) hr += 10.0;
            if (sc->stress_day == day && tod >= 11.0 && tod < 12.0) hr += 22.0;
            if (sc->nap_day == day && tod >= 13.0 && tod < 13.3333) {
                p->stage[ENERGY_LIGHT] = 0.9;
                p->stage[ENERGY_WAKE] = 0.1;
                hr = 58.0;
                hrv = log(44.0);
                steps = 0.0;
            }
            if (sc->exercise_day == day && tod >= 18.0 && tod < 18.5) {
                hr = 150.0;
                steps = 650.0;
            }
        }
        p->hr = hr + sc->hr_offset;
        p->hrv_log = hrv;
        p->steps = steps;
        p->stress = 0.25;
        if (asleep && sc->spo2 > 0.0) {
            /* “醒来那天”的夜 = 前一日 23:00 起到当日 07:00 */
            const int wake_day = (tod >= 23.0) ? day + 1 : day;
            const int k = (int)i % 20; /* 确定性的样本抽样，不用随机数 */
            p->spo2 = sc->spo2;
            if (sc->spo2_bad_day == wake_day) {
                p->spo2 = sc->spo2 - 3.0;
                if (k < 3) p->spo2 -= 4.0; /* 15% 的样本瞬时下跌 */
            }
            if (sc->spo2_min_day == wake_day && k < 4) p->spo2_min = sc->spo2 - 4.0; /* 20% 的样本窗口内最低值下跌 */
            if (sc->spo2_sparse_day == wake_day && k >= 6) p->spo2 = NAN; /* 只剩 30% 覆盖 */
        }
    }
    *n_out = n;
    return s;
}

static int run_batch(size_t days, const scenario_t *sc, const vitality_config_t *vc, vitality_output_t *out)
{
    size_t n;
    energy_sample_t *s = make_panel(days, sc, &n);
    CHECK(s != NULL);
    if (!s) return ENERGY_ENOMEM;
    const int rc = vitality_compute(s, n, NULL, vc, out);
    free(s);
    return rc;
}

static const vitality_point_t *point_at(const vitality_output_t *out, int day, double tod)
{
    for (size_t i = 0; i < out->n_points; i++) {
        const vitality_point_t *p = &out->points[i];
        if ((int)floor(p->t_h / 24.0) == day && fabs(p->tod - tod) < 1e-6) return p;
    }
    return NULL;
}

/* ------------------------------------------------------------------ */
static void test_defaults_and_validate(void)
{
    vitality_config_t c;
    vitality_config_default(&c);
    CHECK(vitality_config_validate(&c, NULL, 0) == ENERGY_OK);
    CHECK(strcmp(vitality_version(), VITALITY_VERSION) == 0);
    CHECK_NEAR(c.dyn.tau0_h, 0.75, 1e-12);
    CHECK_NEAR(c.display.rho, 0.15, 1e-12);
    CHECK(c.use_obs == 1 && c.detect_events == 1);

    char buf[128];
    vitality_config_t bad = c;
    bad.dyn.tau0_h = 0.0;
    CHECK(vitality_config_validate(&bad, buf, sizeof(buf)) == ENERGY_EINVAL);
    CHECK(buf[0] != '\0');
    bad = c;
    bad.dyn.d_cap = 1.5;
    CHECK(vitality_config_validate(&bad, NULL, 0) == ENERGY_EINVAL);
    bad = c;
    bad.obs.act_zero = bad.obs.rest_activity;
    CHECK(vitality_config_validate(&bad, NULL, 0) == ENERGY_EINVAL);
    bad = c;
    bad.event.load_off = bad.event.load_on;
    CHECK(vitality_config_validate(&bad, NULL, 0) == ENERGY_EINVAL);
    bad = c;
    bad.event.a_nap = NAN;
    CHECK(vitality_config_validate(&bad, NULL, 0) == ENERGY_EINVAL);
    bad = c;
    bad.morning.lambda = 0.0;
    CHECK(vitality_config_validate(&bad, NULL, 0) == ENERGY_EINVAL);
    bad = c;
    bad.display.score_mode = 3;
    CHECK(vitality_config_validate(&bad, NULL, 0) == ENERGY_EINVAL);
    bad = c;
    bad.display.rho = 0.9;
    CHECK(vitality_config_validate(&bad, NULL, 0) == ENERGY_EINVAL);
    CHECK(vitality_config_validate(NULL, NULL, 0) == ENERGY_EINVAL);

    /* 非法配置：ctx 标记错误，step 返回 EINVAL */
    vitality_ctx_t ctx;
    vitality_ctx_init(&ctx, NULL, &bad);
    energy_sample_t s;
    energy_sample_init(&s);
    s.t_h = 0.0;
    CHECK(vitality_ctx_step(&ctx, &s, NULL, NULL) == ENERGY_EINVAL);
    CHECK(vitality_ctx_event(&ctx, VITALITY_EVENT_CAFFEINE, 0.0) == ENERGY_EINVAL);
}

static void test_subfunctions(void)
{
    vitality_config_t c;
    vitality_config_default(&c);
    /* τ_D：能量越低越慢 */
    CHECK_NEAR(vitality_tau(&c.dyn, 1.0), 0.75, 1e-12);
    CHECK_NEAR(vitality_tau(&c.dyn, 0.3), 0.75 * (1.0 + 1.5 * 0.7), 1e-12);
    CHECK(vitality_tau(&c.dyn, 0.0) > vitality_tau(&c.dyn, 0.5));
    CHECK_NEAR(vitality_tau(&c.dyn, NAN), vitality_tau(&c.dyn, 0.0), 1e-12);

    /* 均衡线：峰值钟点高于低谷，E=0 → 0，E=1 → ≤ 1 */
    const double v_peak = vitality_equilibrium(&c.display, 0.8, 10.0);
    const double v_dip = vitality_equilibrium(&c.display, 0.8, 14.5);
    CHECK(v_peak > v_dip);
    CHECK(v_peak > 0.8 && v_dip < 0.8);
    CHECK_NEAR(vitality_equilibrium(&c.display, 0.0, 10.0), 0.0, 1e-12);
    CHECK(vitality_equilibrium(&c.display, 1.0, 10.0) <= 1.0);
    for (double tod = 0.0; tod < 24.0; tod += 0.25) {
        const double v = vitality_equilibrium(&c.display, 0.6, tod);
        CHECK(v >= 0.6 * (1.0 - c.display.rho) - 1e-9 && v <= 0.6 * (1.0 + c.display.rho) + 1e-9);
    }

    /* 倒 U：峰值在 z_c/2，z_c 处为 0，更高为负；负 z 为犯困 */
    const double zc = c.obs.z_c;
    CHECK(vitality_arousal(&c.obs, zc / 2.0) > vitality_arousal(&c.obs, zc / 4.0));
    CHECK(vitality_arousal(&c.obs, zc / 2.0) > vitality_arousal(&c.obs, 0.9 * zc));
    CHECK_NEAR(vitality_arousal(&c.obs, zc), 0.0, 1e-12);
    CHECK(vitality_arousal(&c.obs, 1.5 * zc) < 0.0);
    CHECK(vitality_arousal(&c.obs, 6.0) >= c.obs.y_min);
    CHECK_NEAR(vitality_arousal(&c.obs, -1.0), -c.obs.drowsy_gain, 1e-12);
    CHECK_NEAR(vitality_arousal(&c.obs, 0.0), 0.0, 1e-12);
    CHECK(isnan(vitality_arousal(&c.obs, NAN)));

    /* 裁剪：上限 d_cap·E，下限 −V* */
    CHECK_NEAR(vitality_d_bounds(&c.dyn, 0.5, 0.5, 1.0), 0.15, 1e-12);
    CHECK_NEAR(vitality_d_bounds(&c.dyn, 0.5, 0.5, -1.0), -0.5, 1e-12);
    CHECK_NEAR(vitality_d_bounds(&c.dyn, 0.5, 0.5, 0.1), 0.1, 1e-12);
    CHECK_NEAR(vitality_d_bounds(&c.dyn, 0.5, 0.5, NAN), 0.0, 1e-12);

    /* 象限 */
    CHECK(vitality_quadrant(&c.display, 0.8, 0.8) == VITALITY_Q_HIGH_HIGH);
    CHECK(vitality_quadrant(&c.display, 0.8, 0.2) == VITALITY_Q_HIGH_LOW);
    CHECK(vitality_quadrant(&c.display, 0.2, 0.8) == VITALITY_Q_LOW_HIGH);
    CHECK(vitality_quadrant(&c.display, 0.2, 0.2) == VITALITY_Q_LOW_LOW);
}

static void test_batch_basic(void)
{
    scenario_t sc;
    scenario_default(&sc);
    vitality_output_t out;
    CHECK(run_batch(4, &sc, NULL, &out) == ENERGY_OK);
    CHECK(out.n_points == 288 * 4);
    CHECK(out.n_daily == 4);
    CHECK(out.energy.n_points == out.n_points);
    int n_wake = 0, n_sleep = 0, n_learning = 0;
    for (size_t i = 0; i < out.n_points; i++) {
        const vitality_point_t *p = &out.points[i];
        CHECK(p->v >= 0.0 && p->v <= 1.0);
        CHECK(p->v_star >= 0.0 && p->v_star <= 1.0);
        CHECK(p->vitality >= 0.0 && p->vitality <= 100.0);
        CHECK(p->score >= 0.0 && p->score <= 100.0);
        CHECK(p->lower <= p->vitality + 1e-9 && p->upper >= p->vitality - 1e-9);
        CHECK(isfinite(p->d) && isfinite(p->d_sd) && p->d_sd >= 0.0);
        CHECK(fabs(p->v - (p->v_star + p->d)) < 1e-9 || p->v == 0.0 || p->v == 1.0);
        CHECK_NEAR(p->e, out.energy.points[i].e, 1e-12);
        CHECK_NEAR(p->energy, out.energy.points[i].energy, 1e-12);
        if (p->flags & VITALITY_FLAG_WAKE) n_wake++;
        if (p->flags & VITALITY_FLAG_SLEEP) n_sleep++;
        if (p->flags & VITALITY_FLAG_LEARNING) n_learning++;
        /* 睡眠中不做观测 */
        if (p->flags & VITALITY_FLAG_SLEEP) CHECK(p->flags & VITALITY_FLAG_NO_OBS);
    }
    CHECK(n_wake == 4);          /* 每天一次晨起初始化 */
    CHECK(n_sleep == 4 * 8 * 12); /* 23:00–07:00 */
    CHECK(n_learning > 0 && n_learning < (int)out.n_points); /* 前 3 晚学习期，之后退出 */
    CHECK(out.habit.n_nights == 4);
    CHECK(out.habit.sleep_h > 7.5 && out.habit.sleep_h < 8.5);
    CHECK(out.habit.deep_frac > 0.2 && out.habit.deep_frac < 0.4);
    CHECK(out.stats_min <= out.stats_mean && out.stats_mean <= out.stats_max);
    /* 每日汇总 */
    for (size_t d = 0; d < out.n_daily; d++) {
        const vitality_daily_t *dd = &out.daily[d];
        CHECK(dd->n > 0);
        CHECK(isfinite(dd->mean) && dd->min <= dd->mean && dd->mean <= dd->max);
        CHECK(isfinite(dd->wake_score));
        CHECK(isfinite(dd->d_init));
    }
    /* 白天的节律：10 点附近高于 14:30 低谷（无事件日） */
    const vitality_point_t *p10 = point_at(&out, 3, 10.0);
    const vitality_point_t *p14 = point_at(&out, 3, 14.5);
    CHECK(p10 && p14 && p10->v_star > p14->v_star);
    vitality_output_free(&out);
    CHECK(out.points == NULL && out.n_points == 0);
}

static void test_stream_equals_batch(void)
{
    scenario_t sc;
    scenario_default(&sc);
    sc.nap_day = 2;
    sc.exercise_day = 2;
    size_t n;
    energy_sample_t *s = make_panel(3, &sc, &n);
    CHECK(s != NULL);
    if (!s) return;
    vitality_output_t out;
    CHECK(vitality_compute(s, n, NULL, NULL, &out) == ENERGY_OK);

    vitality_ctx_t ctx;
    vitality_ctx_init(&ctx, NULL, NULL);
    energy_ctx_t ectx;
    energy_ctx_init(&ectx, NULL);
    double max_dv = 0.0, max_dd = 0.0;
    for (size_t i = 0; i < n; i++) {
        vitality_point_t vp;
        energy_point_t ep;
        /* 与能量层测试同口径：能量层用批量算出的 load/activity（内部测试入口）推进，再用
         * vitality_ctx_feed 喂给精力层（批量负荷用全窗口标定，流式用在线标定，二者本就允许
         * 不同），这样比较的才是精力层本身 */
        CHECK(energy_ctx_step_fixed(&ectx, &s[i], out.energy.points[i].load,
                                    out.energy.points[i].activity, &ep) == ENERGY_OK);
        CHECK(vitality_ctx_feed(&ctx, &s[i], &ep, &vp) == ENERGY_OK);
        CHECK_NEAR(ep.e, out.energy.points[i].e, 1e-12);
        /* 状态量（v/d/v_star）逐点一致；vitality 只差在批量的 15 min 平滑 */
        if (fabs(vp.v - out.points[i].v) > max_dv) max_dv = fabs(vp.v - out.points[i].v);
        if (fabs(vp.d - out.points[i].d) > max_dd) max_dd = fabs(vp.d - out.points[i].d);
        CHECK((vp.flags & ~VITALITY_FLAG_PCT_READY) == (out.points[i].flags & ~VITALITY_FLAG_PCT_READY));
    }
    CHECK(max_dv < 1e-9);
    CHECK(max_dd < 1e-9);
    vitality_habit_t h;
    vitality_ctx_get_habit(&ctx, &h);
    CHECK(h.n_nights == out.habit.n_nights);
    CHECK_NEAR(h.sleep_h, out.habit.sleep_h, 1e-12);
    vitality_output_free(&out);
    free(s);
}

static void test_events(void)
{
    scenario_t sc;
    scenario_default(&sc);
    sc.nap_day = 3;
    sc.exercise_day = 3;
    sc.arousal_day = 3;
    sc.stress_day = 3;
    vitality_output_t out;
    CHECK(run_batch(5, &sc, NULL, &out) == ENERGY_OK);
    int n_nap = 0, n_ex = 0, n_arousal = 0, n_stress = 0;
    double d_before_nap = 0.0, d_after_nap = 0.0, d_before_ex = 0.0, d_after_ex = 0.0;
    for (size_t i = 0; i < out.n_points; i++) {
        const vitality_point_t *p = &out.points[i];
        const int day = (int)floor(p->t_h / 24.0);
        if (p->flags & VITALITY_FLAG_NAP) {
            n_nap++;
            CHECK(day == 3 && p->tod > 13.2 && p->tod < 13.6);
            d_before_nap = out.points[i - 1].d;
            d_after_nap = p->d;
        }
        if (p->flags & VITALITY_FLAG_EXERCISE) {
            n_ex++;
            CHECK(day == 3 && p->tod > 18.4 && p->tod < 18.8);
            d_before_ex = out.points[i - 1].d;
            d_after_ex = p->d;
        }
        if (p->flags & VITALITY_FLAG_AROUSAL) {
            n_arousal++;
            CHECK(day == 3 && p->tod >= 9.5 && p->tod < 10.25);
        }
        if (p->flags & VITALITY_FLAG_STRESS) {
            n_stress++;
            CHECK(day == 3 && p->tod >= 11.0 && p->tod < 12.0);
        }
    }
    CHECK(n_nap == 1);
    CHECK(n_ex == 1);
    CHECK(d_after_nap - d_before_nap > 0.08);   /* 20 min 午睡 → 约 +0.15（观测可能略拉回） */
    CHECK(d_after_ex - d_before_ex > 0.06);     /* 30 min 跑步 → 约 +0.12 */
    CHECK(n_arousal >= 6);                      /* 45 min 的 +10 bpm 至少 6 个点打上唤醒 */
    CHECK(n_stress >= 6);                       /* 1 h 的 +22 bpm，持续 10 min 后打过载 */
    /* 轻度唤醒抬高 D，过载压低 D */
    const vitality_point_t *a0 = point_at(&out, 3, 9.5), *a1 = point_at(&out, 3, 10.25);
    const vitality_point_t *s0 = point_at(&out, 3, 11.0), *s1 = point_at(&out, 3, 12.0);
    CHECK(a0 && a1 && a1->d > a0->d + 0.01);
    CHECK(s0 && s1 && s1->d < s0->d - 0.01);
    /* 过载时段的精力（相对均衡线）应低于唤醒时段 */
    CHECK(a1->v - a1->v_star > s1->v - s1->v_star);
    /* 运动中不用观测 */
    const vitality_point_t *ex = point_at(&out, 3, 18.25);
    CHECK(ex && (ex->flags & VITALITY_FLAG_NO_OBS));
    /* 每日汇总计数 */
    CHECK(out.daily[3].n_nap == 1 && out.daily[3].n_exercise == 1 && out.daily[3].n_stress_samples >= 6);
    CHECK(out.daily[2].n_nap == 0 && out.daily[2].n_exercise == 0);
    vitality_output_free(&out);

    /* 关闭事件检测：没有脉冲标志 */
    vitality_config_t vc;
    vitality_config_default(&vc);
    vc.detect_events = 0;
    CHECK(run_batch(5, &sc, &vc, &out) == ENERGY_OK);
    for (size_t i = 0; i < out.n_points; i++) {
        CHECK((out.points[i].flags & (VITALITY_FLAG_NAP | VITALITY_FLAG_EXERCISE | VITALITY_FLAG_AROUSAL |
                                      VITALITY_FLAG_STRESS | VITALITY_FLAG_FATIGUE)) == 0);
    }
    vitality_output_free(&out);
}

static void test_user_events_and_ceiling(void)
{
    scenario_t sc;
    scenario_default(&sc);
    size_t n;
    energy_sample_t *s = make_panel(2, &sc, &n);
    CHECK(s != NULL);
    if (!s) return;
    vitality_ctx_t ctx;
    vitality_ctx_init(&ctx, NULL, NULL);
    CHECK(vitality_ctx_event(&ctx, VITALITY_EVENT_CUSTOM, 0.0) == ENERGY_EINVAL);
    CHECK(vitality_ctx_event(&ctx, (vitality_event_type_t)99, 0.1) == ENERGY_EINVAL);

    double d_prev = 0.0, d_jump = 0.0, d_caf_first = 0.0, d_caf_full = 0.0;
    int seen_user = 0;
    for (size_t i = 0; i < n; i++) {
        vitality_point_t vp;
        const int day = (int)floor(s[i].t_h / 24.0);
        if (day == 1 && fabs(s[i].tod - 10.0) < 1e-9) {
            CHECK(vitality_ctx_event(&ctx, VITALITY_EVENT_MEDITATION, 0.0) == ENERGY_OK); /* 即时 +0.08 */
        }
        if (day == 1 && fabs(s[i].tod - 15.0) < 1e-9) {
            CHECK(vitality_ctx_event(&ctx, VITALITY_EVENT_CAFFEINE, 0.0) == ENERGY_OK); /* 一阶起效 */
        }
        CHECK(vitality_ctx_step(&ctx, &s[i], NULL, &vp) == ENERGY_OK);
        if (day == 1 && fabs(s[i].tod - 10.0) < 1e-9) {
            CHECK(vp.flags & VITALITY_FLAG_USER);
            d_jump = vp.d;
            seen_user++;
        }
        if (day == 1 && fabs(s[i].tod - 15.0) < 1e-9) d_caf_first = vp.d;
        if (day == 1 && fabs(s[i].tod - 16.0) < 1e-9) d_caf_full = vp.d;
        if (day == 1 && fabs(s[i].tod - (10.0 - 1.0 / 12.0)) < 1e-9) d_prev = vp.d;
    }
    CHECK(seen_user == 1);
    CHECK(d_jump - d_prev > 0.05); /* 冥想即时注入 0.08（观测会略拉回） */
    CHECK(d_caf_first > 0.0 && d_caf_first < 0.08); /* 咖啡第一步只到 1−e^(−5/20) ≈ 22% */
    (void)d_caf_full;

    /* 天花板：能量低时 D 不能超过 d_cap·E */
    vitality_ctx_t c2;
    energy_config_t ec;
    energy_config_default(&ec);
    ec.x0_e = 0.2;
    ec.use_ekf = 0;
    vitality_config_t vc;
    vitality_config_default(&vc);
    vc.use_obs = 0;
    vitality_ctx_init(&c2, &ec, &vc);
    for (int k = 0; k < 5; k++) CHECK(vitality_ctx_event(&c2, VITALITY_EVENT_CUSTOM, 0.5) == ENERGY_OK);
    energy_sample_t a;
    energy_sample_init(&a);
    a.t_h = 0.0;
    a.tod = 10.0;
    a.hr = 70.0;
    vitality_point_t vp;
    CHECK(vitality_ctx_step(&c2, &a, NULL, &vp) == ENERGY_OK);
    CHECK(vp.d <= vc.dyn.d_cap * vp.e + 1e-12);
    CHECK(vp.d > 0.0);
    CHECK(vp.v <= 1.0);
    free(s);
}

static void test_night_resume_merge(void)
{
    /* 夜里醒 20 min 再睡：两段应合并为同一夜（常态只学一次，时长不被截断） */
    scenario_t sc;
    scenario_default(&sc);
    vitality_output_t ref, out;
    CHECK(run_batch(4, &sc, NULL, &ref) == ENERGY_OK);
    sc.broken_day = 2;
    CHECK(run_batch(4, &sc, NULL, &out) == ENERGY_OK);
    int n_wake_day2 = 0;
    double last_d_init = NAN;
    for (size_t i = 0; i < out.n_points; i++) {
        const vitality_point_t *p = &out.points[i];
        if ((p->flags & VITALITY_FLAG_WAKE) && (int)floor(p->t_h / 24.0) == 2) {
            n_wake_day2++;
            last_d_init = p->d_init;
            CHECK(isfinite(p->d_init));
        } else if (!(p->flags & VITALITY_FLAG_WAKE)) {
            CHECK(isnan(p->d_init));
        }
    }
    CHECK(n_wake_day2 == 2);                        /* 03:20 一次、07:10 一次 */
    CHECK_NEAR(out.daily[2].d_init, last_d_init, 0.0); /* 日汇总取早上那次 */
    CHECK(out.habit.n_nights == ref.habit.n_nights); /* 常态仍只按 4 晚学习 */
    CHECK(fabs(out.habit.sleep_h - ref.habit.sleep_h) < 0.3);
    CHECK(out.habit.wake_frac > ref.habit.wake_frac); /* 夜醒计入碎片化 */
    vitality_output_free(&ref);
    vitality_output_free(&out);
}

static void test_morning_init(void)
{
    /* 冷启动后：短夜（5 h）的 D_init 必须低于正常夜 */
    scenario_t sc;
    scenario_default(&sc);
    sc.short_day = 5;
    vitality_output_t out;
    CHECK(run_batch(7, &sc, NULL, &out) == ENERGY_OK);
    double d_init[7];
    for (int d = 0; d < 7; d++) d_init[d] = out.daily[d].d_init;
    CHECK(isfinite(d_init[4]) && isfinite(d_init[5]) && isfinite(d_init[6]));
    CHECK(d_init[5] < d_init[4] - 0.02); /* 短夜 */
    CHECK(d_init[5] < 0.0);
    CHECK(d_init[6] > d_init[5]);         /* 恢复正常后回升 */
    CHECK_NEAR(d_init[0], 0.0, 1e-12);    /* 第一晚没有常态可比：不猜 */
    /* 晨起点带 WAKE，且 LEARNING 只在前 3 晚 */
    int learning_after = 0;
    for (size_t i = 0; i < out.n_points; i++) {
        const vitality_point_t *p = &out.points[i];
        if ((p->flags & VITALITY_FLAG_WAKE) && p->t_h > 24.0 * 4 && (p->flags & VITALITY_FLAG_LEARNING)) learning_after++;
    }
    CHECK(learning_after == 0);
    vitality_output_free(&out);

    /* 画像恢复：set/get habit */
    vitality_ctx_t ctx;
    vitality_ctx_init(&ctx, NULL, NULL);
    vitality_habit_t h;
    memset(&h, 0, sizeof(h));
    h.n_nights = 5;
    h.delta_e = 0.3;
    h.hrv_wake = 3.8;
    h.sleep_h = 7.5;
    h.deep_frac = 0.25;
    h.wake_frac = 0.05;
    CHECK(vitality_ctx_set_habit(&ctx, &h) == ENERGY_OK);
    vitality_habit_t g;
    vitality_ctx_get_habit(&ctx, &g);
    CHECK(g.n_nights == 5 && fabs(g.sleep_h - 7.5) < 1e-12);
    h.deep_frac = 1.5;
    CHECK(vitality_ctx_set_habit(&ctx, &h) == ENERGY_EINVAL);
    h.deep_frac = 0.25;
    h.n_nights = -1;
    CHECK(vitality_ctx_set_habit(&ctx, &h) == ENERGY_EINVAL);
    CHECK(vitality_ctx_set_habit(&ctx, NULL) == ENERGY_EINVAL);
    /* 带画像启动：第一晚就不在学习期 */
    h.n_nights = 5;
    CHECK(vitality_ctx_set_habit(&ctx, &h) == ENERGY_OK);
    energy_sample_t a;
    energy_sample_init(&a);
    a.t_h = 0.0;
    a.tod = 12.0;
    a.hr = 70.0;
    vitality_point_t vp;
    CHECK(vitality_ctx_step(&ctx, &a, NULL, &vp) == ENERGY_OK);
    CHECK((vp.flags & VITALITY_FLAG_LEARNING) == 0);
}

static void test_percentile_and_score_modes(void)
{
    scenario_t sc;
    scenario_default(&sc);
    sc.nap_day = 4;
    vitality_config_t vc;
    vitality_config_default(&vc);
    vc.display.score_mode = 1;
    vitality_output_t out;
    CHECK(run_batch(6, &sc, &vc, &out) == ENERGY_OK);
    int ready_day0 = 0, ready_last = 0;
    for (size_t i = 0; i < out.n_points; i++) {
        const vitality_point_t *p = &out.points[i];
        const int day = (int)floor(p->t_h / 24.0);
        if (p->flags & VITALITY_FLAG_PCT_READY) {
            CHECK(isfinite(p->percentile) && p->percentile >= 0.0 && p->percentile <= 100.0);
            CHECK_NEAR(p->score, p->percentile, 1e-12);
            if (day == 0) ready_day0++;
            if (day == 5) ready_last++;
        } else {
            CHECK(isnan(p->percentile));
            CHECK_NEAR(p->score, p->vitality, 1e-12); /* 学习期退回绝对刻度 */
        }
    }
    CHECK(ready_day0 == 0);
    CHECK(ready_last > 0);
    /* 午睡回弹点的分位应高于其前一点 */
    for (size_t i = 1; i < out.n_points; i++) {
        if (out.points[i].flags & VITALITY_FLAG_NAP) {
            CHECK(out.points[i].percentile > out.points[i - 1].percentile);
        }
    }
    vitality_output_free(&out);

    vc.display.score_mode = 2;
    CHECK(run_batch(4, &sc, &vc, &out) == ENERGY_OK);
    for (size_t i = 0; i < out.n_points; i++) {
        const vitality_point_t *p = &out.points[i];
        if (isfinite(p->percentile)) CHECK_NEAR(p->score, 0.5 * (p->vitality + p->percentile), 1e-9);
    }
    vitality_output_free(&out);
}

static void test_constitution_invariance(void)
{
    /* 同一条曲线整体 +8 bpm：精力应基本不变（能量层个体基线吸收平移） */
    scenario_t a, b;
    scenario_default(&a);
    scenario_default(&b);
    a.nap_day = b.nap_day = 3;
    a.arousal_day = b.arousal_day = 3;
    b.hr_offset = 8.0;
    vitality_output_t oa, ob;
    CHECK(run_batch(5, &a, NULL, &oa) == ENERGY_OK);
    CHECK(run_batch(5, &b, NULL, &ob) == ENERGY_OK);
    double max_diff = 0.0;
    for (size_t i = 0; i < oa.n_points; i++) {
        if (oa.points[i].t_h < 24.0 * 3) continue; /* 学习期后 */
        const double d = fabs(oa.points[i].score - ob.points[i].score);
        if (d > max_diff) max_diff = d;
    }
    CHECK(max_diff < 3.0);
    vitality_output_free(&oa);
    vitality_output_free(&ob);
}

static void test_open_loop(void)
{
    scenario_t sc;
    scenario_default(&sc);
    vitality_config_t vc;
    vitality_config_default(&vc);
    vc.use_obs = 0;
    vc.detect_events = 0;
    vitality_output_t out;
    CHECK(run_batch(3, &sc, &vc, &out) == ENERGY_OK);
    for (size_t i = 0; i < out.n_points; i++) {
        const vitality_point_t *p = &out.points[i];
        CHECK(p->flags & VITALITY_FLAG_NO_OBS);
        /* 没有观测与事件时，D 只能来自晨起初始化并衰减：|D| 单调不增（除晨起点） */
        if (i > 0 && !(p->flags & VITALITY_FLAG_WAKE)) CHECK(fabs(p->d) <= fabs(out.points[i - 1].d) + 1e-12);
    }
    vitality_output_free(&out);
}

static void test_invalid_and_edges(void)
{
    vitality_output_t out;
    CHECK(vitality_compute(NULL, 10, NULL, NULL, &out) == ENERGY_EINVAL);
    energy_sample_t s[3];
    for (int i = 0; i < 3; i++) {
        energy_sample_init(&s[i]);
        s[i].t_h = (double)i / 12.0;
        s[i].hr = 70.0;
    }
    CHECK(vitality_compute(s, 0, NULL, NULL, &out) == ENERGY_EINVAL);
    CHECK(vitality_compute(s, 3, NULL, NULL, NULL) == ENERGY_EINVAL);
    vitality_config_t bad;
    vitality_config_default(&bad);
    bad.dyn.p0 = -1.0;
    CHECK(vitality_compute(s, 3, NULL, &bad, &out) == ENERGY_EINVAL);
    vitality_output_free(&out); /* 失败后也可安全释放 */
    /* 时间不单调 */
    s[2].t_h = s[1].t_h;
    CHECK(vitality_compute(s, 3, NULL, NULL, &out) == ENERGY_EDOMAIN);
    vitality_output_free(&out);

    vitality_ctx_t ctx;
    vitality_ctx_init(&ctx, NULL, NULL);
    CHECK(vitality_ctx_step(&ctx, NULL, NULL, NULL) == ENERGY_EINVAL);
    CHECK(vitality_ctx_step(NULL, &s[0], NULL, NULL) == ENERGY_EINVAL);
    CHECK(vitality_ctx_step(&ctx, &s[0], NULL, NULL) == ENERGY_OK);
    CHECK(vitality_ctx_step(&ctx, &s[0], NULL, NULL) == ENERGY_EDOMAIN);
    energy_point_t ep;
    CHECK(vitality_ctx_feed(&ctx, &s[1], NULL, NULL) == ENERGY_EINVAL);
    memset(&ep, 0, sizeof(ep));
    ep.t_h = NAN;
    CHECK(vitality_ctx_feed(&ctx, &s[1], &ep, NULL) == ENERGY_EINVAL);
    /* peek 不推进时间 */
    vitality_point_t vp;
    vitality_ctx_peek(&ctx, -1.0, &vp);
    CHECK(vp.v >= 0.0 && vp.v <= 1.0 && isfinite(vp.score));
    vitality_ctx_peek(&ctx, 15.0, &vp);
    CHECK_NEAR(vp.tod, 15.0, 1e-12);
    vitality_ctx_peek(NULL, 15.0, &vp);
    vitality_ctx_peek(&ctx, 15.0, NULL);
    /* 未启动的 ctx peek */
    vitality_ctx_t fresh;
    vitality_ctx_init(&fresh, NULL, NULL);
    vitality_ctx_peek(&fresh, 9.0, &vp);
    CHECK(isfinite(vp.v) && vp.v >= 0.0 && vp.v <= 1.0);
    vitality_output_free(NULL);
    vitality_ctx_init(NULL, NULL, NULL);
    vitality_ctx_get_habit(&ctx, NULL);
}

static unsigned rng_state = 12345u;
static double rnd(void)
{
    rng_state = rng_state * 1103515245u + 12345u;
    return (double)((rng_state >> 8) & 0xFFFFFF) / 16777216.0;
}

static void test_nan_fuzz_and_gaps(void)
{
    scenario_t sc;
    scenario_default(&sc);
    sc.nap_day = 1;
    sc.exercise_day = 1;
    size_t n;
    energy_sample_t *s = make_panel(3, &sc, &n);
    CHECK(s != NULL);
    if (!s) return;
    for (size_t i = 0; i < n; i++) {
        const double r = rnd();
        if (r < 0.10) s[i].hr = NAN;
        else if (r < 0.15) s[i].hr = 300.0 * rnd();
        if (rnd() < 0.30) s[i].hrv_log = NAN;
        if (rnd() < 0.05) s[i].steps = NAN;
        if (rnd() < 0.05) s[i].stage[ENERGY_LIGHT] = NAN;
        if (rnd() < 0.03) s[i].q_hr = 0.0;
    }
    /* 挖一个 3 小时缺口 */
    size_t m = 0;
    for (size_t i = 0; i < n; i++) {
        if (s[i].t_h > 24.0 + 12.0 && s[i].t_h < 24.0 + 15.0) continue;
        s[m++] = s[i];
    }
    vitality_output_t out;
    CHECK(vitality_compute(s, m, NULL, NULL, &out) == ENERGY_OK);
    int n_gap = 0;
    for (size_t i = 0; i < out.n_points; i++) {
        const vitality_point_t *p = &out.points[i];
        CHECK(isfinite(p->v) && p->v >= 0.0 && p->v <= 1.0);
        CHECK(isfinite(p->score) && p->score >= 0.0 && p->score <= 100.0);
        CHECK(isfinite(p->d) && isfinite(p->d_sd));
        CHECK(p->lower >= 0.0 && p->upper <= 100.0 && p->lower <= p->upper);
        if (p->flags & VITALITY_FLAG_GAP) n_gap++;
    }
    CHECK(n_gap >= 1);
    vitality_output_free(&out);

    /* 流式：同样的脏数据不崩、有界 */
    vitality_ctx_t ctx;
    vitality_ctx_init(&ctx, NULL, NULL);
    for (size_t i = 0; i < m; i++) {
        vitality_point_t vp;
        CHECK(vitality_ctx_step(&ctx, &s[i], NULL, &vp) == ENERGY_OK);
        CHECK(isfinite(vp.v) && vp.v >= 0.0 && vp.v <= 1.0);
    }
    free(s);
}

static void test_spo2(void)
{
    const size_t days = 6;
    scenario_t sc0;
    scenario_default(&sc0);
    vitality_output_t base, steady, bad, sparse, vmin;
    CHECK(run_batch(days, &sc0, NULL, &base) == ENERGY_OK);

    /* ① 恒定血氧（每晚一样）：常态学到 97 / 0 下跌，x_SpO2 恒为 0 → 与无血氧逐位一致 */
    scenario_t sc1 = sc0;
    sc1.spo2 = 97.0;
    CHECK(run_batch(days, &sc1, NULL, &steady) == ENERGY_OK);
    CHECK(steady.n_points == base.n_points);
    CHECK(steady.habit.n_spo2_nights == (int)days && steady.habit.n_nights == (int)days);
    CHECK_NEAR(steady.habit.spo2_mean, 97.0, 1e-9);
    CHECK_NEAR(steady.habit.spo2_dip, 0.0, 1e-12);
    int n_spo2_pts = 0, n_breath = 0;
    for (size_t i = 0; i < base.n_points; i++) {
        CHECK(steady.points[i].v == base.points[i].v);
        CHECK(steady.points[i].d == base.points[i].d);
        if (steady.points[i].flags & VITALITY_FLAG_SPO2) {
            n_spo2_pts++;
            CHECK(steady.points[i].flags & VITALITY_FLAG_WAKE);
            CHECK_NEAR(steady.points[i].spo2_night, 97.0, 1e-9);
            CHECK(steady.points[i].spo2_dip == 0.0);
        } else {
            CHECK(isnan(steady.points[i].spo2_night) && isnan(steady.points[i].spo2_dip));
        }
        if (steady.points[i].flags & VITALITY_FLAG_BREATH) n_breath++;
    }
    CHECK(n_spo2_pts == (int)days);
    CHECK(n_breath == 0);
    for (size_t d = 0; d < steady.n_daily; d++) {
        CHECK_NEAR(steady.daily[d].spo2_night, 97.0, 1e-9);
        CHECK(steady.daily[d].spo2_dip == 0.0);
    }
    /* 无血氧时 daily 的血氧字段为 NAN */
    CHECK(isnan(base.daily[1].spo2_night) && isnan(base.daily[1].spo2_dip));

    /* ② 第 4 日醒来前那一夜血氧不稳 → 当日 D_init 低于恒定血氧版本，且只在该夜打 BREATH；
     *    次日（夜间恢复正常，但常态已被拉低一点）D_init 回到 ≥ 恒定版本 */
    scenario_t sc2 = sc1;
    sc2.spo2_bad_day = 4;
    CHECK(run_batch(days, &sc2, NULL, &bad) == ENERGY_OK);
    CHECK(isfinite(bad.daily[4].d_init) && isfinite(steady.daily[4].d_init));
    /* x_SpO2 = −0.5（均值 −3 ≥ 1.5 满分、下跌 15% ≥ 10% 满分）→ Δ = gain·w_spo2·0.5 = 0.04 */
    CHECK_NEAR(bad.daily[4].d_init, steady.daily[4].d_init - 0.04, 1e-9);
    CHECK(bad.daily[4].wake_score < steady.daily[4].wake_score - 2.0);
    CHECK(bad.daily[4].spo2_night < 95.0 && bad.daily[4].spo2_night > 93.0);
    CHECK_NEAR(bad.daily[4].spo2_dip, 0.15, 0.02);
    CHECK(bad.daily[5].d_init >= steady.daily[5].d_init);
    CHECK(bad.daily[5].d_init > bad.daily[4].d_init);
    for (size_t d = 0; d < 4; d++) CHECK(bad.daily[d].d_init == steady.daily[d].d_init);
    n_breath = 0;
    for (size_t i = 0; i < bad.n_points; i++) {
        if (bad.points[i].flags & VITALITY_FLAG_BREATH) {
            n_breath++;
            CHECK((bad.points[i].flags & (VITALITY_FLAG_WAKE | VITALITY_FLAG_SPO2)) ==
                  (VITALITY_FLAG_WAKE | VITALITY_FLAG_SPO2));
            CHECK((int)floor(bad.points[i].t_h / 24.0) == 4);
        }
    }
    CHECK(n_breath == 1);
    /* 坏夜之后日间 V 整体低于恒定版本（D 均值回复，差距随时间收窄但不反号） */
    const vitality_point_t *b8 = point_at(&bad, 4, 8.0), *s8 = point_at(&steady, 4, 8.0);
    const vitality_point_t *b12 = point_at(&bad, 4, 12.0), *s12 = point_at(&steady, 4, 12.0);
    CHECK(b8 && s8 && b12 && s12);
    CHECK(b8->v < s8->v && b12->v <= s12->v && (s8->v - b8->v) > (s12->v - b12->v));

    /* ③ 覆盖不足（30%）的夜：不用血氧，D_init 与恒定版本一致，不打 SPO2 标志，常态不计入 */
    scenario_t sc3 = sc1;
    sc3.spo2_sparse_day = 4;
    CHECK(run_batch(days, &sc3, NULL, &sparse) == ENERGY_OK);
    CHECK(sparse.daily[4].d_init == steady.daily[4].d_init);
    CHECK(isnan(sparse.daily[4].spo2_night));
    CHECK(sparse.habit.n_spo2_nights == (int)days - 1);
    const vitality_point_t *w4 = NULL;
    for (size_t i = 0; i < sparse.n_points; i++) {
        if ((sparse.points[i].flags & VITALITY_FLAG_WAKE) && (int)floor(sparse.points[i].t_h / 24.0) == 4)
            w4 = &sparse.points[i];
    }
    CHECK(w4 && (w4->flags & VITALITY_FLAG_SPO2) == 0);

    /* ④ 只有 spo2_min 下跌（均值不变）：下跌占比 20% → x = −0.5·1 = −0.5 → Δ = −0.04（裁剪后同满分） */
    scenario_t sc4 = sc1;
    sc4.spo2_min_day = 4;
    CHECK(run_batch(days, &sc4, NULL, &vmin) == ENERGY_OK);
    CHECK_NEAR(vmin.daily[4].spo2_night, 97.0, 1e-9);
    CHECK_NEAR(vmin.daily[4].spo2_dip, 0.20, 0.02);
    CHECK_NEAR(vmin.daily[4].d_init, steady.daily[4].d_init - 0.04, 1e-9);

    /* ⑤ 流式 = 批量（带血氧） */
    {
        size_t n;
        energy_sample_t *smp = make_panel(days, &sc2, &n);
        CHECK(smp != NULL);
        vitality_ctx_t ctx;
        vitality_ctx_init(&ctx, NULL, NULL);
        energy_ctx_t ectx;
        energy_ctx_init(&ectx, NULL);
        double max_dd = 0.0;
        for (size_t i = 0; i < n; i++) {
            vitality_point_t vp;
            energy_point_t ep; /* 同 test_stream_equals_batch：负荷/活动量用批量标定值 */
            CHECK(energy_ctx_step_fixed(&ectx, &smp[i], bad.energy.points[i].load,
                                        bad.energy.points[i].activity, &ep) == ENERGY_OK);
            CHECK(vitality_ctx_feed(&ctx, &smp[i], &ep, &vp) == ENERGY_OK);
            const double dd = fabs(vp.d - bad.points[i].d);
            if (dd > max_dd) max_dd = dd;
            CHECK((vp.flags & (VITALITY_FLAG_SPO2 | VITALITY_FLAG_BREATH)) ==
                  (bad.points[i].flags & (VITALITY_FLAG_SPO2 | VITALITY_FLAG_BREATH)));
        }
        CHECK(max_dd < 1e-9);
        vitality_habit_t h;
        vitality_ctx_get_habit(&ctx, &h);
        CHECK(h.n_spo2_nights == bad.habit.n_spo2_nights);
        CHECK_NEAR(h.spo2_dip, bad.habit.spo2_dip, 1e-12);
        /* 画像往返：含血氧常态 */
        vitality_ctx_t c2;
        vitality_ctx_init(&c2, NULL, NULL);
        CHECK(vitality_ctx_set_habit(&c2, &h) == ENERGY_OK);
        vitality_habit_t g;
        vitality_ctx_get_habit(&c2, &g);
        CHECK(g.n_spo2_nights == h.n_spo2_nights && g.spo2_mean == h.spo2_mean);
        h.spo2_mean = 120.0;
        CHECK(vitality_ctx_set_habit(&c2, &h) == ENERGY_EINVAL);
        h.spo2_mean = 97.0;
        h.n_spo2_nights = h.n_nights + 1;
        CHECK(vitality_ctx_set_habit(&c2, &h) == ENERGY_EINVAL);
        h.n_spo2_nights = 0;
        CHECK(vitality_ctx_set_habit(&c2, &h) == ENERGY_OK);
        vitality_ctx_get_habit(&c2, &g);
        CHECK(isnan(g.spo2_mean) && isnan(g.spo2_dip));
        free(smp);
    }

    /* ⑥ 非法血氧值（0 / 120 / 负数 / Inf）视为缺失：与无血氧逐位一致 */
    {
        size_t n;
        energy_sample_t *smp = make_panel(days, &sc0, &n);
        CHECK(smp != NULL);
        for (size_t i = 0; i < n; i++) {
            const double junk[4] = {0.0, 120.0, -5.0, INFINITY};
            smp[i].spo2 = junk[i % 4];
            smp[i].spo2_min = junk[(i + 1) % 4];
        }
        vitality_output_t o;
        CHECK(vitality_compute(smp, n, NULL, NULL, &o) == ENERGY_OK);
        CHECK(o.habit.n_spo2_nights == 0 && isnan(o.habit.spo2_mean));
        for (size_t i = 0; i < n; i++) CHECK(o.points[i].v == base.points[i].v);
        vitality_output_free(&o);
        free(smp);
    }

    /* ⑦ 配置校验 */
    {
        vitality_config_t c;
        vitality_config_default(&c);
        c.morning.w_spo2 = NAN;
        CHECK(vitality_config_validate(&c, NULL, 0) == ENERGY_EINVAL);
        vitality_config_default(&c);
        c.morning.spo2_min_cover = 1.5;
        CHECK(vitality_config_validate(&c, NULL, 0) == ENERGY_EINVAL);
        vitality_config_default(&c);
        c.morning.spo2_dip_delta = 0.0;
        CHECK(vitality_config_validate(&c, NULL, 0) == ENERGY_EINVAL);
        /* w_spo2 = 0：有血氧也不影响 D_init */
        vitality_config_default(&c);
        c.morning.w_spo2 = 0.0;
        vitality_output_t o;
        CHECK(run_batch(days, &sc2, &c, &o) == ENERGY_OK);
        CHECK(o.daily[4].d_init == steady.daily[4].d_init);
        CHECK(isfinite(o.daily[4].spo2_night)); /* 统计照样给出 */
        vitality_output_free(&o);
    }

    vitality_output_free(&base);
    vitality_output_free(&steady);
    vitality_output_free(&bad);
    vitality_output_free(&sparse);
    vitality_output_free(&vmin);
}

static void test_energy_residual_exposed(void)
{
    /* energy_point_t.hr_resid：基线建立前 NAN，之后有限；夜间 / 日间都有 */
    scenario_t sc;
    scenario_default(&sc);
    sc.arousal_day = 2;
    size_t n;
    energy_sample_t *s = make_panel(3, &sc, &n);
    CHECK(s != NULL);
    if (!s) return;
    energy_output_t eo;
    CHECK(energy_compute(s, n, NULL, &eo) == ENERGY_OK);
    CHECK(isnan(eo.points[0].hr_resid)); /* 首样本：还没有个人基线 */
    int finite_day = 0, finite_night = 0;
    double resid_arousal = 0.0;
    int n_ar = 0;
    for (size_t i = 0; i < n; i++) {
        const energy_point_t *p = &eo.points[i];
        if (isfinite(p->hr_resid)) {
            if (p->sleep_prob >= 0.5) finite_night++;
            else finite_day++;
        }
        if (p->t_h >= 48.0 + 9.5 && p->t_h < 48.0 + 10.25 && isfinite(p->hr_resid)) {
            resid_arousal += p->hr_resid;
            n_ar++;
        }
    }
    CHECK(finite_day > 100 && finite_night > 100);
    CHECK(n_ar > 0 && resid_arousal / (double)n_ar > 5.0); /* +10 bpm 的升高大部分留在残差里 */
    energy_output_free(&eo);
    free(s);
}

int main(void)
{
    test_defaults_and_validate();
    test_subfunctions();
    test_batch_basic();
    test_stream_equals_batch();
    test_events();
    test_user_events_and_ceiling();
    test_night_resume_merge();
    test_morning_init();
    test_percentile_and_score_modes();
    test_constitution_invariance();
    test_open_loop();
    test_invalid_and_edges();
    test_nan_fuzz_and_gaps();
    test_energy_residual_exposed();
    test_spo2();
    if (g_fail) {
        printf("FAILED: %d of %d checks\n", g_fail, g_checks);
        return 1;
    }
    printf("OK: %d checks, 0 failures\n", g_checks);
    return 0;
}
