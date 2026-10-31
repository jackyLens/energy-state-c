/*
 * evobeing_demo.c —— 主接口调用示例：每 5 分钟喂一个样本，取能量与精力
 *
 * 合成 3 天数据（23:00–07:00 睡眠，日间静坐，第 2 天 18:00 跑步 30 分钟），
 * 打印每天 08:00 / 12:00 / 18:30 / 22:00 的能量与精力。
 *
 * 运行：make example  或  build/evobeing_demo
 */
#include "evobeing_zhenghuimin.h"

#include <math.h>
#include <stdio.h>

static void make_sample(evobeing_sample_t *s, double t)
{
    evobeing_sample_init_zhenghuimin(s);
    const double tod = fmod(t, 24.0);
    const int day = (int)floor(t / 24.0);
    const int asleep = (tod >= 23.0 || tod < 7.0);
    s->t_h = t;
    s->tod = tod;
    if (asleep) {
        const double phase = fmod(tod + 1.0, 24.0) / 8.0;
        const double deep = (phase < 0.5) ? 0.45 : 0.15;
        s->stage[2] = deep;                                /* 深睡 */
        s->stage[3] = (phase < 0.5) ? 0.15 : 0.35;         /* REM */
        s->stage[1] = 1.0 - deep - s->stage[3] - 0.05;     /* 浅睡 */
        s->stage[0] = 0.05;                                /* 清醒 */
        s->hr = 54.0 + 3.0 * phase;
        s->hrv_log = log(48.0 + 4.0 * sin(phase * 6.28));  /* 夜间逐拍 RR → RMSSD → ln */
        s->spo2 = 97.0;
        s->spo2_min = 96.0;
    } else {
        const double w = 2.0 * M_PI / 24.0;
        s->hr = 70.0 - 2.0 * cos(w * tod) - 3.46 * sin(w * tod);
        s->hrv_log = NAN;                                  /* 日间只有心率 */
        s->steps = 40.0;
        s->kcal = 2.0;
        if (day == 1 && tod >= 18.0 && tod < 18.5) {       /* 跑步 30 分钟 */
            s->hr = 150.0;
            s->steps = 650.0;
            s->kcal = 40.0;
        }
    }
    s->stress = 25.0;                                      /* 设备 0–100 量纲 */
}

int main(void)
{
    evobeing_ctx_t ctx;
    evobeing_init_zhenghuimin(&ctx, NULL);

    printf("evobeing %s (%s)\n", EVOBEING_VERSION, EVOBEING_AUTHOR);
    printf("%4s %6s %8s %8s\n", "day", "tod", "energy", "vitality");
    const int n = 3 * 24 * 12;
    for (int i = 0; i < n; i++) {
        evobeing_sample_t s;
        make_sample(&s, (double)i / 12.0);
        double energy = 0.0, vitality = 0.0;
        if (evobeing_energy_zhenghuimin(&ctx, &s, &energy) != EVOBEING_OK) continue;
        if (evobeing_vitality_zhenghuimin(&ctx, &vitality) != EVOBEING_OK) continue;
        const double tod = s.tod;
        if (fabs(tod - 8.0) < 1e-9 || fabs(tod - 12.0) < 1e-9 || fabs(tod - 18.5) < 1e-9 ||
            fabs(tod - 22.0) < 1e-9) {
            printf("%4d %6.2f %8.1f %8.1f\n", i / 288, tod, energy, vitality);
        }
    }

    evobeing_profile_t profile;
    evobeing_get_profile_zhenghuimin(&ctx, &profile);
    printf("已学习夜数：能量基线 %d 晚，精力常态 %d 晚\n",
           profile.energy.n[ENERGY_BL_HR_NIGHT], profile.vitality.n_nights);
    return 0;
}
