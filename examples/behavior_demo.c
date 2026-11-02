/*
 * behavior_demo.c — 行为归因模块的最小用法
 *
 * 合成 28 天的晨起精力和行为记录（饮酒、咖啡、小睡、运动、冥想、晚睡），
 * 打印每种行为对晨起精力的影响百分比。
 */
#include "evobeing_behavior_zhenghuimin.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

static unsigned s = 2024u;
static double urand(void) { s ^= s << 13; s ^= s >> 17; s ^= s << 5; return (double)(s % 1000000u) / 1e6; }
static double nrand(void) { double u1 = urand() + 1e-9, u2 = urand(); return sqrt(-2.0 * log(u1)) * cos(6.283185307 * u2); }

int main(void)
{
    enum { N = 28 };
    double y[N], wake[N];
    int weekday[N];
    evobeing_behavior_event_t ev[256];
    int n_ev = 0, d, i;
    evobeing_behavior_input_t in;
    evobeing_behavior_output_t out;

    /* ---- 合成数据：真实效应 饮酒 −8 分、晚间咖啡 −4、小睡 +3、运动 +3、晚睡 −5 ---- */
    for (d = 0; d < N; ++d) {
        double eff = 0.0, base = (d - 1) * 24.0;
        weekday[d] = (d + 2) % 7;
        wake[d] = d * 24.0 + 7.0 + 0.4 * nrand();
        if (d > 0) {
            if (urand() < 0.25) { ev[n_ev].id = EVOBEING_BEHAVIOR_ALCOHOL;    ev[n_ev++].t_h = base + 20.5; eff -= 8; }
            if (urand() < 0.7)  { double h = urand() < 0.6 ? 9.0 : 16.0; ev[n_ev].id = EVOBEING_BEHAVIOR_CAFFEINE; ev[n_ev++].t_h = base + h; if (h >= 14) eff -= 4; }
            if (urand() < 0.3)  { ev[n_ev].id = EVOBEING_BEHAVIOR_NAP;        ev[n_ev++].t_h = base + 13.5; eff += 3; }
            if (urand() < 0.4)  { ev[n_ev].id = EVOBEING_BEHAVIOR_EXERCISE;   ev[n_ev++].t_h = base + 18.0; eff += 3; }
            if (urand() < 0.2)  { ev[n_ev].id = EVOBEING_BEHAVIOR_MEDITATION; ev[n_ev++].t_h = base + 22.0; }
            if (urand() < 0.25) { ev[n_ev].id = EVOBEING_BEHAVIOR_LATE_BED;   ev[n_ev++].t_h = base + 24.5; eff -= 5; }
        }
        y[d] = 60.0 + eff + 4.0 * nrand();
    }

    memset(&in, 0, sizeof in);
    in.n_days = N; in.wake_score = y; in.wake_t_h = wake; in.weekday = weekday;
    in.events = ev; in.n_events = n_ev;

    if (evobeing_behavior_impact_zhenghuimin(&in, NULL, &out) != 0) {
        printf("有效天数不足\n");
        return 1;
    }
    printf("过去 %d 天，晨起精力均值 %.0f（标准差 %.1f）\n", out.n_days_used, out.mean_score, out.sd_score);
    printf("%-16s %6s %8s %16s %8s  %s\n", "行为", "次数", "影响", "90% 区间", "份额", "可信度");
    for (i = 0; i < out.n_behaviors; ++i) {
        const evobeing_behavior_item_t *it = &out.items[i];
        if (it->status == EVOBEING_BEHAVIOR_INSUFFICIENT) {
            printf("%-16s %3d/%-3d  再记录 %d 次就能看到\n", it->name, it->n_with, it->n_without,
                   it->n_with < 4 ? 4 - it->n_with : 0);
            continue;
        }
        printf("%-16s %3d/%-3d %+7.0f%% [%+5.0f%%, %+5.0f%%] %6.0f%%  %s\n", it->name, it->n_with, it->n_without,
               it->effect_pct, it->ci_lo_pct, it->ci_hi_pct, it->share_pct,
               it->status == EVOBEING_BEHAVIOR_HIGH ? "高" : it->status == EVOBEING_BEHAVIOR_MID ? "中" : "低");
    }
    printf("记录的行为合计解释 %.0f%%，其余 %.0f%% 来自睡眠本身、身体状况等未记录因素\n",
           out.explained_pct, out.unexplained_pct);
    return 0;
}
