/*
 * vitality_demo.c —— 精力层流式示例：能量 + 精力双轨
 *
 * 合成 6 天数据（23:00–07:00 睡眠），第 4 天安排四个场景：
 *   09:30 一次 10 bpm 的无运动心率升高（开会 / 轻度唤醒）
 *   13:00 午睡 20 分钟（分期给 light）
 *   15:00 用户标记一杯咖啡（手动事件）
 *   18:00 跑步 30 分钟，之后看回弹
 * 第 5 天只睡 5 小时，看晨起初始化与天花板效应。
 * 第 6 天睡够 8 小时，但前一夜血氧不稳（均值 −3 个百分点、15% 样本瞬时下跌），
 *   看“夜间血氧平稳度”对晨起初始化的影响（v1.1）。其余各夜血氧恒定 97%，与不传血氧逐位一致。
 *
 * 运行：make example  或  build/vitality_demo [--csv out.csv]
 */
#include "vitality.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

static void make_sample(energy_sample_t *s, double t, int short_night, double *hr_out)
{
    energy_sample_init(s);
    const double tod = fmod(t, 24.0);
    const int day = (int)floor(t / 24.0);
    const double wake_at = (short_night && day == 4) ? 5.0 : 7.0;
    const double sleep_at = (short_night && day == 4) ? 0.0 : 23.0;
    int asleep = (tod >= 23.0 || tod < wake_at);
    if (short_night && day == 4) asleep = (tod >= sleep_at && tod < wake_at);
    s->t_h = t;
    s->tod = tod;
    double hr, hrv, steps = 0.0;
    if (asleep) {
        const double phase = fmod(tod + 1.0, 24.0) / 8.0; /* 0..1 within the night */
        const double deep = (phase < 0.5) ? 0.45 : 0.15;
        s->stage[ENERGY_DEEP] = deep;
        s->stage[ENERGY_REM] = (phase < 0.5) ? 0.15 : 0.35;
        s->stage[ENERGY_LIGHT] = 1.0 - deep - s->stage[ENERGY_REM] - 0.05;
        s->stage[ENERGY_WAKE] = 0.05;
        hr = 54.0 + 3.0 * phase;
        hrv = log(48.0 + 4.0 * sin(phase * 6.28));
    } else {
        const double w = 2.0 * M_PI / 24.0;
        hr = 70.0 - 2.0 * cos(w * tod) - 3.46 * sin(w * tod) + 0.75 * cos(2.0 * w * tod) - 1.30 * sin(2.0 * w * tod);
        hrv = log(34.0);
        steps = 40.0;
        if (day == 3) {
            if (tod >= 9.5 && tod < 10.25) hr += 10.0;                 /* 轻度唤醒：+10 bpm 无运动 */
            if (tod >= 13.0 && tod < 13.3333) {                        /* 午睡 20 min */
                s->stage[ENERGY_LIGHT] = 0.9;
                s->stage[ENERGY_WAKE] = 0.1;
                hr = 58.0;
                hrv = log(44.0);
                steps = 0.0;
            }
            if (tod >= 18.0 && tod < 18.5) {                           /* 跑步 30 min */
                hr = 150.0;
                steps = 650.0;
            }
        }
    }
    s->hr = hr;
    s->hrv_log = hrv;
    s->steps = steps;
    s->stress = 0.25;
    if (asleep) {
        /* 戒指夜间血氧：平时 97%；第 6 天醒来前那一夜（day 4 的 23:00 起）不稳 */
        const int wake_day = (tod >= 23.0) ? day + 1 : day;
        s->spo2 = 97.0;
        if (wake_day == 5) {
            s->spo2 = 94.0;
            if (((int)floor(t * 12.0 + 0.5)) % 20 < 3) s->spo2_min = 90.0; /* 15% 的窗口内最低值下跌 */
        }
    }
    if (hr_out) *hr_out = hr;
}

int main(int argc, char **argv)
{
    const char *csv = NULL;
    for (int i = 1; i + 1 < argc; i++) {
        if (strcmp(argv[i], "--csv") == 0) csv = argv[i + 1];
    }
    energy_config_t ecfg;
    energy_config_default(&ecfg);
    vitality_config_t vcfg;
    vitality_config_default(&vcfg);
    vitality_ctx_t ctx;
    vitality_ctx_init(&ctx, &ecfg, &vcfg);

    FILE *fp = csv ? fopen(csv, "w") : NULL;
    if (fp) fprintf(fp, "t_h,tod,hr,e,energy,v_star,d,d_init,v,vitality,percentile,score,lower,upper,flags,spo2_night,spo2_dip\n");

    printf("vitality-state-c %s（能量 %s）\n", vitality_version(), energy_version());
    printf("第 4 天（day=3）场景：09:30 轻度唤醒 / 13:00 午睡 20min / 15:00 咖啡标记 / 18:00 跑步 30min\n\n");
    printf(" 钟点   心率   能量  均衡线   D     精力  标志\n");
    const int n = 288 * 6;
    for (int i = 0; i < n; i++) {
        const double t = (double)i / 12.0;
        energy_sample_t s;
        double hr;
        make_sample(&s, t, 1, &hr);
        const int day = (int)floor(t / 24.0);
        if (day == 3 && fabs(fmod(t, 24.0) - 15.0) < 1e-9) {
            vitality_ctx_event(&ctx, VITALITY_EVENT_CAFFEINE, 0.0);
        }
        energy_point_t ep;
        vitality_point_t vp;
        if (vitality_ctx_step(&ctx, &s, &ep, &vp) != ENERGY_OK) {
            fprintf(stderr, "step 失败 @%d\n", i);
            return 1;
        }
        if (fp) {
            fprintf(fp, "%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%.6g,%u,%.6g,%.6g\n", vp.t_h,
                    vp.tod, hr, vp.e, vp.energy, vp.v_star, vp.d, vp.d_init, vp.v, vp.vitality, vp.percentile, vp.score,
                    vp.lower, vp.upper, vp.flags, vp.spo2_night, vp.spo2_dip);
        }
        const double tod = fmod(t, 24.0);
        const int show = (day == 3 || day == 4 || day == 5) && (fabs(tod * 4.0 - floor(tod * 4.0 + 1e-9)) < 1e-9) && tod >= 6.0 && tod < 22.0;
        const unsigned ev = vp.flags & (VITALITY_FLAG_WAKE | VITALITY_FLAG_EXERCISE | VITALITY_FLAG_NAP |
                                        VITALITY_FLAG_USER | VITALITY_FLAG_STRESS);
        if (show || ev) {
            char tag[128] = "";
            if (vp.flags & VITALITY_FLAG_WAKE) strcat(tag, "起床 ");
            if (vp.flags & VITALITY_FLAG_BREATH) strcat(tag, "昨夜血氧不稳 ");
            if (vp.flags & VITALITY_FLAG_AROUSAL) strcat(tag, "唤醒 ");
            if (vp.flags & VITALITY_FLAG_STRESS) strcat(tag, "过载 ");
            if (vp.flags & VITALITY_FLAG_NAP) strcat(tag, "午睡回弹 ");
            if (vp.flags & VITALITY_FLAG_EXERCISE) strcat(tag, "运动回弹 ");
            if (vp.flags & VITALITY_FLAG_USER) strcat(tag, "用户标记 ");
            if (vp.flags & VITALITY_FLAG_SLEEP) strcat(tag, "睡眠 ");
            if (vp.flags & VITALITY_FLAG_LEARNING) strcat(tag, "学习期 ");
            if (day == 4 && tod == 6.0) printf(" ---- 第 5 天：只睡 5 小时 ----\n");
            if (day == 5 && tod == 6.0) printf(" ---- 第 6 天：睡够 8 小时，但夜间血氧不稳 ----\n");
            if ((vp.flags & VITALITY_FLAG_WAKE) && (vp.flags & VITALITY_FLAG_SPO2))
                printf("        [夜间血氧 均值 %.1f%%  瞬时下跌占比 %.0f%%  D_init %+.3f]\n", vp.spo2_night,
                       100.0 * vp.spo2_dip, vp.d_init);
            printf(" %02d:%02d  %5.0f  %5.1f  %5.1f  %+.3f  %5.1f  %s\n", (int)floor(tod),
                   (int)floor((tod - floor(tod)) * 60.0 + 0.5), hr, vp.energy, 100.0 * vp.v_star, vp.d, vp.score, tag);
        }
    }
    if (fp) fclose(fp);
    vitality_habit_t h;
    vitality_ctx_get_habit(&ctx, &h);
    printf("\n学到的个人常态：夜数 %d，夜间回充 ΔE=%.3f，醒前 HRV=%.2f（ln），睡眠 %.1f h，深睡 %.0f%%，碎片 %.0f%%\n",
           h.n_nights, h.delta_e, h.hrv_wake, h.sleep_h, 100.0 * h.deep_frac, 100.0 * h.wake_frac);
    printf("                夜间血氧（%d 夜）：均值 %.1f%%，瞬时下跌占比 %.1f%%\n", h.n_spo2_nights, h.spo2_mean,
           100.0 * h.spo2_dip);
    return 0;
}
