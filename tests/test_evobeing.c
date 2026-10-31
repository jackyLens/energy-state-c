/*
 * test_evobeing.c —— 主接口测试：
 *   1. 主接口结果与直接调用 energy_ctx_step / vitality_ctx_step 逐位一致；
 *   2. 调用顺序与错误码；
 *   3. 画像导出 → 重新 init → 导入后状态一致；
 *   4. 输出有界。
 */
#include "evobeing_zhenghuimin.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

static int g_checks = 0, g_fail = 0;
#define CHECK(c) do { g_checks++; if (!(c)) { g_fail++; fprintf(stderr, "FAIL %s:%d %s\n", __FILE__, __LINE__, #c); } } while (0)

static void make_sample(evobeing_sample_t *s, double t)
{
    evobeing_sample_init_zhenghuimin(s);
    const double tod = fmod(t, 24.0);
    const int asleep = (tod >= 23.0 || tod < 7.0);
    s->t_h = t;
    s->tod = tod;
    if (asleep) {
        s->stage[1] = 0.6; s->stage[2] = 0.25; s->stage[3] = 0.1; s->stage[0] = 0.05;
        s->hr = 55.0;
        s->hrv_log = log(46.0);
        s->spo2 = 97.0;
    } else {
        s->hr = 72.0 + 3.0 * sin(t);
        s->steps = (fmod(t, 5.0) < 0.5) ? 300.0 : 30.0;
        s->kcal = 2.0;
    }
    s->stress = 0.2;
}

static void to_energy(const evobeing_sample_t *in, energy_sample_t *out)
{
    energy_sample_init(out);
    out->t_h = in->t_h; out->tod = in->tod;
    memcpy(out->stage, in->stage, sizeof(out->stage));
    out->hr = in->hr; out->hrv_log = in->hrv_log;
    out->steps = in->steps; out->kcal = in->kcal; out->stress = in->stress;
    out->q_hr = in->q_hr; out->q_hrv = in->q_hrv;
    out->spo2 = in->spo2; out->spo2_min = in->spo2_min;
}

static void test_matches_engines(void)
{
    evobeing_ctx_t ctx;
    evobeing_init_zhenghuimin(&ctx, NULL);
    vitality_ctx_t vc;
    vitality_ctx_init(&vc, NULL, NULL);

    const int n = 4 * 288;
    for (int i = 0; i < n; i++) {
        evobeing_sample_t s;
        make_sample(&s, (double)i / 12.0);
        double e = -1.0, v = -1.0;
        CHECK(evobeing_energy_zhenghuimin(&ctx, &s, &e) == EVOBEING_OK);
        CHECK(evobeing_vitality_zhenghuimin(&ctx, &v) == EVOBEING_OK);

        energy_sample_t es;
        to_energy(&s, &es);
        energy_point_t ep; vitality_point_t vp;
        CHECK(vitality_ctx_step(&vc, &es, &ep, &vp) == ENERGY_OK);
        CHECK(e == ep.energy);
        CHECK(v == vp.vitality);
        CHECK(ctx.vitality_point.flags == vp.flags);
        CHECK(e >= 0.0 && e <= 100.0);
        CHECK(v >= 0.0 && v <= 100.0);
    }
}

static void test_call_order(void)
{
    evobeing_ctx_t ctx;
    evobeing_init_zhenghuimin(&ctx, NULL);
    double v = 5.0, e = 5.0;
    CHECK(evobeing_vitality_zhenghuimin(&ctx, &v) == ENERGY_EINVAL); /* 未算能量 */
    CHECK(v == 5.0);
    evobeing_sample_t s;
    make_sample(&s, 0.0);
    CHECK(evobeing_energy_zhenghuimin(&ctx, &s, &e) == EVOBEING_OK);
    CHECK(evobeing_vitality_zhenghuimin(&ctx, &v) == EVOBEING_OK);
    CHECK(evobeing_vitality_zhenghuimin(&ctx, &v) == ENERGY_EINVAL); /* 同一样本只喂一次 */
    CHECK(evobeing_energy_zhenghuimin(&ctx, &s, &e) == ENERGY_EDOMAIN); /* 时间未递增 */
    CHECK(evobeing_energy_zhenghuimin(NULL, &s, &e) == ENERGY_EINVAL);
    CHECK(evobeing_energy_zhenghuimin(&ctx, NULL, &e) == ENERGY_EINVAL);
    CHECK(evobeing_event_zhenghuimin(&ctx, VITALITY_EVENT_CAFFEINE, 0.0) == EVOBEING_OK);
    make_sample(&s, 1.0 / 12.0);
    CHECK(evobeing_energy_zhenghuimin(&ctx, &s, NULL) == EVOBEING_OK);
    CHECK(evobeing_vitality_zhenghuimin(&ctx, NULL) == EVOBEING_OK);
}

static void test_profile_roundtrip(void)
{
    evobeing_ctx_t a;
    evobeing_init_zhenghuimin(&a, NULL);
    const int n = 5 * 288;
    for (int i = 0; i < n; i++) {
        evobeing_sample_t s;
        make_sample(&s, (double)i / 12.0);
        CHECK(evobeing_energy_zhenghuimin(&a, &s, NULL) == EVOBEING_OK);
        CHECK(evobeing_vitality_zhenghuimin(&a, NULL) == EVOBEING_OK);
    }
    evobeing_profile_t p;
    evobeing_get_profile_zhenghuimin(&a, &p);
    CHECK(p.energy.n[ENERGY_BL_HR_NIGHT] >= 3);
    CHECK(p.vitality.n_nights >= 3);

    /* 用画像重新初始化：基线与常态应被沿用，不再处于学习期 */
    evobeing_ctx_t b;
    evobeing_init_zhenghuimin(&b, &p);
    CHECK(b.config_error == 0);
    energy_baseline_t bl;
    energy_ctx_get_baseline(&b.energy, &bl);
    CHECK(memcmp(&bl, &p.energy, sizeof(bl)) == 0);
    vitality_habit_t h;
    vitality_ctx_get_habit(&b.vitality, &h);
    CHECK(h.n_nights == p.vitality.n_nights);
    CHECK(h.delta_e == p.vitality.delta_e);

    evobeing_sample_t s;
    make_sample(&s, 12.0);
    double e = 0.0;
    CHECK(evobeing_energy_zhenghuimin(&b, &s, &e) == EVOBEING_OK);
    CHECK((b.energy_point.flags & ENERGY_FLAG_LEARNING) == 0);

    /* 空画像等价于 NULL */
    evobeing_profile_t zero;
    memset(&zero, 0, sizeof(zero));
    evobeing_ctx_t c, d;
    evobeing_init_zhenghuimin(&c, &zero);
    evobeing_init_zhenghuimin(&d, NULL);
    double ec = 0.0, ed = 0.0;
    CHECK(evobeing_energy_zhenghuimin(&c, &s, &ec) == EVOBEING_OK);
    CHECK(evobeing_energy_zhenghuimin(&d, &s, &ed) == EVOBEING_OK);
    CHECK(ec == ed);
}

static void test_sample_init(void)
{
    evobeing_sample_t s;
    memset(&s, 0xAA, sizeof(s));
    evobeing_sample_init_zhenghuimin(&s);
    CHECK(isnan(s.hr) && isnan(s.hrv_log) && isnan(s.spo2) && isnan(s.spo2_min));
    CHECK(s.q_hr == 1.0 && s.q_hrv == 1.0);
    CHECK(s.t_h == 0.0 && s.steps == 0.0 && s.stage[0] == 0.0);
    evobeing_sample_init_zhenghuimin(NULL);
    evobeing_init_zhenghuimin(NULL, NULL);
}

int main(void)
{
    test_sample_init();
    test_matches_engines();
    test_call_order();
    test_profile_roundtrip();
    printf("%s: %d checks, %d failures\n", g_fail ? "FAILED" : "OK", g_checks, g_fail);
    return g_fail ? 1 : 0;
}
