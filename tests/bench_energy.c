/* clock_gettime / struct timespec need POSIX.1b under strict -std=c99. */
#define _POSIX_C_SOURCE 199309L
/*
 * bench_energy.c —— 纯计算性能（不含 CSV I/O）
 *
 * 用法：make bench
 */

#include "energy.h"
#include "vitality.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

static double now_ms(void)
{
    struct timespec ts;
    clock_gettime(CLOCK_MONOTONIC, &ts);
    return (double)ts.tv_sec * 1000.0 + (double)ts.tv_nsec / 1e6;
}

int main(void)
{
    const size_t days = 21;
    const size_t n = 288 * days;
    energy_sample_t *s = (energy_sample_t *)calloc(n, sizeof(energy_sample_t));
    if (!s) {
        fprintf(stderr, "oom\n");
        return 1;
    }
    for (size_t i = 0; i < n; i++) {
        const double t = (double)i / 12.0;
        const double tod = fmod(t, 24.0);
        const int asleep = (tod >= 23.0 || tod < 7.0);
        s[i].t_h = t;
        s[i].tod = tod;
        s[i].stage[ENERGY_LIGHT] = asleep ? 1.0 : 0.0;
        s[i].stage[ENERGY_DEEP] = (asleep && tod < 3.0) ? 0.4 : 0.0;
        s[i].hr = asleep ? 56.0 : 72.0;
        s[i].hrv_log = asleep ? log(45.0) : log(32.0);
        s[i].steps = asleep ? 0.0 : 40.0;
        s[i].stress = 0.25;
        s[i].q_hr = 1.0;
        s[i].q_hrv = 1.0;
        s[i].load = -1.0;
        s[i].activity = -1.0;
    }

    energy_config_t cfg;
    energy_config_default(&cfg);
    energy_output_t out;

    /* 预热 */
    if (energy_compute(s, n, &cfg, &out) != ENERGY_OK) return 1;
    energy_output_free(&out);

    const int reps = 20;
    const double t0 = now_ms();
    for (int i = 0; i < reps; i++) {
        if (energy_compute(s, n, &cfg, &out) != ENERGY_OK) return 1;
        energy_output_free(&out);
    }
    const double dt = (now_ms() - t0) / (double)reps;

    printf("面板: %zu 点（%.0f 天）\n", n, (double)days);
    printf("耗时: %.2f ms/次，%.1f µs/点\n", dt, dt * 1000.0 / (double)n);
    printf("吞吐: %.0f 点/秒（单核）\n", (double)n / (dt / 1000.0));

    /* 流式：每步只做一次解析推进 + 两次观测更新 */
    energy_ctx_t ctx;
    energy_ctx_init(&ctx, &cfg);
    energy_point_t p;
    const double t1 = now_ms();
    for (int r = 0; r < reps; r++) {
        energy_ctx_init(&ctx, &cfg);
        for (size_t i = 0; i < n; i++) {
            energy_ctx_step(&ctx, &s[i], &p);
        }
    }
    const double dts = (now_ms() - t1) / (double)reps;
    printf("流式: %.2f ms/次，%.2f µs/步\n", dts, dts * 1000.0 / (double)n);

    /* 精力层（流式，内含能量层一步）：增量开销 = 本行 − 上一行 */
    vitality_ctx_t *vctx = (vitality_ctx_t *)malloc(sizeof(*vctx));
    if (!vctx) {
        free(s);
        return 1;
    }
    vitality_point_t vp;
    const double t2 = now_ms();
    for (int r = 0; r < reps; r++) {
        vitality_ctx_init(vctx, &cfg, NULL);
        for (size_t i = 0; i < n; i++) {
            vitality_ctx_step(vctx, &s[i], &p, &vp);
        }
    }
    const double dtv = (now_ms() - t2) / (double)reps;
    printf("能量+精力流式: %.2f ms/次，%.2f µs/步（精力层增量约 %.2f µs/步）\n", dtv,
           dtv * 1000.0 / (double)n, (dtv - dts) * 1000.0 / (double)n);
    printf("上下文大小: energy_ctx_t %zu B，vitality_ctx_t %zu B（含 energy_ctx_t）\n", sizeof(energy_ctx_t),
           sizeof(*vctx));
    free(vctx);

    free(s);
    return 0;
}
