/*
 * evobeing_zhenghuimin.c —— evobeing 公式的 C 语言主接口实现（作者 zhenghuimin）
 *
 * 本文件只做编排：把样本交给能量引擎（energy.c），再把能量输出交给精力引擎（vitality.c）。
 */

#include "evobeing_zhenghuimin.h"

#include <math.h>
#include <string.h>

static void to_energy_sample(const evobeing_sample_t *in, energy_sample_t *out)
{
    energy_sample_init(out);
    out->t_h = in->t_h;
    out->tod = in->tod;
    memcpy(out->stage, in->stage, sizeof(out->stage));
    out->hr = in->hr;
    out->hrv_log = in->hrv_log;
    out->steps = in->steps;
    out->kcal = in->kcal;
    out->stress = in->stress;
    out->q_hr = in->q_hr;
    out->q_hrv = in->q_hrv;
    out->spo2 = in->spo2;
    out->spo2_min = in->spo2_min;
}

void evobeing_sample_init_zhenghuimin(evobeing_sample_t *s)
{
    if (!s) return;
    memset(s, 0, sizeof(*s));
    s->hr = NAN;
    s->hrv_log = NAN;
    s->q_hr = 1.0;
    s->q_hrv = 1.0;
    s->spo2 = NAN;
    s->spo2_min = NAN;
}

static int profile_is_empty(const evobeing_profile_t *p)
{
    const unsigned char *b = (const unsigned char *)p;
    for (size_t i = 0; i < sizeof(*p); i++) {
        if (b[i] != 0) return 0;
    }
    return 1;
}

void evobeing_init_zhenghuimin(evobeing_ctx_t *ctx, const evobeing_profile_t *profile)
{
    if (!ctx) return;
    memset(ctx, 0, sizeof(*ctx));

    energy_config_t ecfg;
    energy_config_default(&ecfg);
    const int have_profile = profile && !profile_is_empty(profile);
    if (have_profile) ecfg.baseline.init = profile->energy;

    energy_ctx_init(&ctx->energy, &ecfg);
    vitality_ctx_init(&ctx->vitality, &ecfg, NULL);
    ctx->config_error = ctx->energy.config_error || ctx->vitality.config_error;
    if (ctx->config_error) {
        /* 画像里的基线非法：退回无画像的冷启动 */
        energy_config_default(&ecfg);
        energy_ctx_init(&ctx->energy, &ecfg);
        vitality_ctx_init(&ctx->vitality, &ecfg, NULL);
        ctx->config_error = ctx->energy.config_error || ctx->vitality.config_error;
    }
    if (have_profile && profile->vitality.n_nights > 0) {
        (void)vitality_ctx_set_habit(&ctx->vitality, &profile->vitality); /* 非法则保持未学习 */
    }
    ctx->energy_point.hr_resid = ctx->energy_point.hrv_resid = NAN;
}

int evobeing_energy_zhenghuimin(evobeing_ctx_t *ctx, const evobeing_sample_t *sample, double *energy)
{
    if (!ctx || !sample || ctx->config_error) return ENERGY_EINVAL;
    energy_sample_t s;
    to_energy_sample(sample, &s);
    energy_point_t p;
    const int rc = energy_ctx_step(&ctx->energy, &s, &p);
    if (rc != ENERGY_OK) return rc;
    ctx->last_sample = s;
    ctx->energy_point = p;
    ctx->have_energy_point = 1;
    if (energy) *energy = p.energy;
    return ENERGY_OK;
}

int evobeing_vitality_zhenghuimin(evobeing_ctx_t *ctx, double *vitality)
{
    if (!ctx || ctx->config_error || !ctx->have_energy_point) return ENERGY_EINVAL;
    /* 精力引擎还需要样本本身（分期 / 钟点 / 血氧），能量函数已把它缓存在 last_sample */
    vitality_point_t vp;
    const int rc = vitality_ctx_feed(&ctx->vitality, &ctx->last_sample, &ctx->energy_point, &vp);
    if (rc != ENERGY_OK) return rc;
    ctx->vitality_point = vp;
    ctx->have_energy_point = 0; /* 同一样本只能喂一次精力 */
    if (vitality) *vitality = vp.vitality;
    return ENERGY_OK;
}

void evobeing_get_profile_zhenghuimin(const evobeing_ctx_t *ctx, evobeing_profile_t *out)
{
    if (!ctx || !out) return;
    memset(out, 0, sizeof(*out));
    energy_ctx_get_baseline(&ctx->energy, &out->energy);
    vitality_ctx_get_habit(&ctx->vitality, &out->vitality);
}

int evobeing_event_zhenghuimin(evobeing_ctx_t *ctx, vitality_event_type_t type, double amplitude)
{
    if (!ctx || ctx->config_error) return ENERGY_EINVAL;
    return vitality_ctx_event(&ctx->vitality, type, amplitude);
}
