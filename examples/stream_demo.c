/*
 * stream_demo.c —— 流式用法示例：每 5 分钟推一次，打印能量
 *
 * 这就是“能量值放在一个变量里、每 5 分钟更新”的正规写法：
 *   - 变量不是文件级 static，而是调用方持有的 energy_ctx_t；
 *   - 更新不是 E += delta，而是解析解推进 + EKF 观测校正（相对个人基线）。
 *
 * 为了演示滤波器确实能跟踪状态，这里先用开环 ODE 生成一条“真值”，
 * 再按观测模型反推出 HR/HRV（相当于合成数据），最后用 EKF 流式估计，
 * 并打印估计误差。最后把整条心率曲线平移 +8 bpm（另一个静息心率更高的人），
 * 演示 v1.2 的核心原则：能量只跟自己比，体质差异进个人基线、不进能量。
 */

#include "energy.h"

#include <math.h>
#include <stdio.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define DAYS 3
#define N (288 * DAYS) /* 5 分钟一格 */

static void fill_sample(energy_sample_t *s, int k)
{
    energy_sample_init(s); /* 不要 memset(0)：q=0 / hr=0 都是错的 */
    const double t_h = (double)k / 12.0;
    const double tod = fmod(t_h, 24.0);
    const int asleep = (tod >= 23.0 || tod < 7.0);
    s->t_h = t_h;
    s->tod = tod;
    /* 白天分期全 0（不要传 wake=1，那会被当成“夜里清醒”） */
    s->stage[ENERGY_LIGHT] = asleep ? 0.6 : 0.0;
    s->stage[ENERGY_DEEP] = (asleep && (tod >= 23.0 || tod < 3.0)) ? 0.4 : 0.0;
    s->stage[ENERGY_REM] = (asleep && !(tod >= 23.0 || tod < 3.0)) ? 0.4 : 0.0;
    /* 白天：普通走动 40 步/5 min，每 3 小时一次 15 分钟快走（550 步/5 min） */
    s->steps = asleep ? 0.0 : (fmod(tod, 3.0) < 0.25 ? 550.0 : 40.0);
    s->stress = 0.25;
}

/* 观测模型反推：给定状态 E，算出模型期望的 HR / ln(RMSSD)（截距用 cfg.obs 的群体先验） */
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
        const double a = energy_activity_step(s->steps, s->kcal, cfg->load.steps_ref);
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
    static energy_sample_t samples[N];
    static double truth_e[N], est_e[N];
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
    energy_config_default(&cfg); /* use_ekf = 1，个体基线默认开启 */
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
    printf("时刻       能量    S      E_est   E_true   带宽           状态\n");
    for (int k = 0; k < N; k++) {
        energy_point_t p;
        if (energy_ctx_step(&ctx, &samples[k], &p) != ENERGY_OK) {
            fprintf(stderr, "step %d failed\n", k);
            return 1;
        }
        est_e[k] = p.e;
        const double err = fabs(p.e - truth_e[k]);
        if (k >= 576 && err > max_err) max_err = err;
        if (k % 36 == 0) {
            printf("D%d %02d:00  %5.1f  %.3f  %.3f   %.3f    [%5.1f,%5.1f]  %s\n", k / 288 + 1,
                   (int)p.tod, p.energy, p.s, p.e, truth_e[k], p.lower, p.upper,
                   (p.flags & ENERGY_FLAG_LEARNING) ? "正在了解你的常态" : "");
        }
    }
    /* Known limitation (v1.2): the daytime baseline level is learned from data and partly
     * absorbs the E signal, so daytime E_est can sit ~0.1 below the open-loop truth.
     * The demo is illustrative only; regression bounds live in tests/test_energy.c. */
    printf("\n第 3 天起（冷启动后）EKF 对 E 的最大误差: %.4f（%.2f 分）\n", max_err, max_err * 100.0);

    energy_baseline_t bl;
    energy_ctx_get_baseline(&ctx, &bl);
    printf("学到的个人基线：夜间 HR %.1f bpm / ln(RMSSD) %.2f（%d 晚），日间 HR %.1f bpm（%d 天），"
           "睡眠窗 %.1f→%.1f\n",
           bl.level[ENERGY_BL_HR_NIGHT], bl.level[ENERGY_BL_HRV_NIGHT], bl.n[ENERGY_BL_HR_NIGHT],
           bl.level[ENERGY_BL_HR_DAY], bl.n[ENERGY_BL_HR_DAY], bl.sleep_on, bl.sleep_off);

    /* 4) 另一个人：恢复模式完全相同，只是心率整体高 8 bpm */
    for (int k = 0; k < N; k++) samples[k].hr += 8.0;
    energy_ctx_init(&ctx, &cfg);
    double max_diff = 0.0;
    for (int k = 0; k < N; k++) {
        energy_point_t p;
        energy_ctx_step(&ctx, &samples[k], &p);
        if (fabs(p.e - est_e[k]) > max_diff) max_diff = fabs(p.e - est_e[k]);
    }
    energy_ctx_get_baseline(&ctx, &bl);
    printf("心率整体 +8 bpm：E 的最大差异 %.2g，夜间 HR 基线变为 %.1f bpm（差异进了基线，没进能量）\n",
           max_diff, bl.level[ENERGY_BL_HR_NIGHT]);
    return 0;
}
