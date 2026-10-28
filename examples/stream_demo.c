/*
 * stream_demo.c —— 流式用法示例：每 5 分钟推一次，打印每小时的能量
 *
 * 这就是“能量值放在一个变量里、每 5 分钟更新”的正规写法：
 *   - 变量不是文件级 static，而是调用方持有的 energy_ctx_t；
 *   - 更新不是 E += delta，而是解析解推进 + EKF 观测校正。
 *
 * 为了演示滤波器确实能跟踪状态，这里先用开环 ODE 生成一条“真值”，
 * 再按观测模型反推出 HR/HRV（相当于合成数据），最后用 EKF 流式估计，
 * 并打印估计误差。
 */

#include "energy.h"

#include <math.h>
#include <stdio.h>

#define N 288 /* 一天，5 分钟一格 */

static void fill_sample(energy_sample_t *s, int k)
{
    const double t_h = (double)k / 12.0;
    const double tod = fmod(t_h, 24.0);
    const int asleep = (tod >= 23.0 || tod < 7.0);
    s->t_h = t_h;
    s->tod = tod;
    s->stage[ENERGY_WAKE] = asleep ? 0.0 : 1.0;
    s->stage[ENERGY_LIGHT] = asleep ? 1.0 : 0.0;
    s->stage[ENERGY_DEEP] = (asleep && tod < 3.0) ? 0.4 : 0.0;
    s->stage[ENERGY_REM] = 0.0;
    s->steps = asleep ? 0.0 : 40.0;
    s->kcal = 0.0;
    s->stress = 0.25;
    s->q_hr = 1.0;
    s->q_hrv = 1.0;
    s->load = -1.0;
    s->activity = -1.0;
    s->hr = NAN;
    s->hrv_log = NAN;
}

/* 观测模型反推：给定状态 E，算出模型期望的 HR / ln(RMSSD) */
static void synth_obs(const energy_config_t *cfg, const energy_sample_t *s, double e,
                      double *hr, double *hrv)
{
    const int night = (s->stage[0] + s->stage[1] + s->stage[2] + s->stage[3]) > 0.5;
    if (night) {
        *hr = cfg->obs.hr_night[0] + cfg->obs.hr_night[1] * e +
              cfg->obs.hr_night[2] * s->stage[ENERGY_DEEP] +
              cfg->obs.hr_night[3] * s->stage[ENERGY_REM] +
              cfg->obs.hr_night[4] * s->stage[ENERGY_WAKE];
        *hrv = cfg->obs.hrv_night[0] + cfg->obs.hrv_night[1] * e +
               cfg->obs.hrv_night[2] * s->stage[ENERGY_DEEP] +
               cfg->obs.hrv_night[3] * s->stage[ENERGY_REM] +
               cfg->obs.hrv_night[4] * s->stage[ENERGY_WAKE];
    } else {
        const double w = 2.0 * M_PI / 24.0;
        const double t = s->tod;
        const double a = 0.8 * (s->steps / 120.0);
        *hr = cfg->obs.hr_day[0] + cfg->obs.hr_day[1] * e + cfg->obs.hr_day[2] * a +
              cfg->obs.hr_day[3] * cos(w * t) + cfg->obs.hr_day[4] * sin(w * t) +
              cfg->obs.hr_day[5] * cos(2 * w * t) + cfg->obs.hr_day[6] * sin(2 * w * t);
        *hrv = cfg->obs.hrv_day[0] + cfg->obs.hrv_day[1] * e + cfg->obs.hrv_day[2] * a +
               cfg->obs.hrv_day[3] * cos(w * t) + cfg->obs.hrv_day[4] * sin(w * t) +
               cfg->obs.hrv_day[5] * cos(2 * w * t) + cfg->obs.hrv_day[6] * sin(2 * w * t);
    }
}

int main(void)
{
    energy_sample_t samples[N];
    double truth_e[N];
    for (int k = 0; k < N; k++) fill_sample(&samples[k], k);

    /* 1) 开环跑一遍，作为“真值” */
    energy_config_t cfg;
    energy_config_default(&cfg);
    /* 生产环境：l_ref 应来自用户画像（长期固定），否则每次窗口的归一化基准不同。
     * 这里固定 1.5，模拟一个“中等活动量”的用户。 */
    cfg.load.l_ref = 1.5;
    cfg.use_ekf = 0;
    energy_ctx_t ctx;
    energy_ctx_init(&ctx, &cfg);
    for (int k = 0; k < N; k++) {
        energy_point_t p;
        energy_ctx_step(&ctx, &samples[k], &p);
        truth_e[k] = p.e;
    }

    /* 2) 按观测模型生成 HR/HRV（相当于戒指测到的值） */
    energy_config_default(&cfg); /* use_ekf = 1 */
    cfg.load.l_ref = 1.5;
    for (int k = 0; k < N; k++) {
        double hr, hrv;
        synth_obs(&cfg, &samples[k], truth_e[k], &hr, &hrv);
        samples[k].hr = hr + 0.5 * sin(0.7 * k); /* 加一点观测噪声 */
        samples[k].hrv_log = hrv - 0.02 * cos(0.3 * k);
    }

    /* 3) 流式 EKF：每 5 分钟一步 */
    energy_ctx_init(&ctx, &cfg);
    double max_err = 0.0;
    printf("时刻   能量    S      E_est   E_true   带宽          \n");
    for (int k = 0; k < N; k++) {
        energy_point_t p;
        if (energy_ctx_step(&ctx, &samples[k], &p) != ENERGY_OK) {
            fprintf(stderr, "step %d failed\n", k);
            return 1;
        }
        const double err = fabs(p.e - truth_e[k]);
        if (err > max_err) max_err = err;
        if (k % 24 == 0) {
            printf("%02d:00  %5.1f  %.3f  %.3f   %.3f    [%5.1f,%5.1f]\n",
                   (int)p.tod, p.energy, p.s, p.e, truth_e[k], p.lower, p.upper);
        }
    }
    printf("\nEKF 对 E 的最大误差: %.4f（%.2f 分）\n", max_err, max_err * 100.0);
    return 0;
}
