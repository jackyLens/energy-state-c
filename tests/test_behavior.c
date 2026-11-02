/*
 * test_behavior.c — 行为归因模块的合成数据测试
 *
 * 合成一个"真实世界"：每天晨起精力 = 本人均值 + 真实行为效应 + 周末效应 + 自相关 + 噪声，
 * 再检查模块能否找回效应的方向与大小、样本不足时是否拒绝给数、混杂是否被控制住等。
 */
#include "evobeing_behavior_zhenghuimin.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

static int n_checks = 0, n_fail = 0;
#define CHECK(cond, msg) do { ++n_checks; if (!(cond)) { ++n_fail; printf("  FAIL %s:%d %s\n", __FILE__, __LINE__, msg); } } while (0)

static unsigned rng_s = 7u;
static double urand(void) { rng_s ^= rng_s << 13; rng_s ^= rng_s >> 17; rng_s ^= rng_s << 5; return (double)(rng_s % 1000000u) / 1000000.0; }
static double nrand(void) { double u1 = urand() + 1e-9, u2 = urand(); return sqrt(-2.0 * log(u1)) * cos(6.283185307 * u2); }

typedef struct {
    double y[EVOBEING_BEHAVIOR_MAX_DAYS];
    double wake[EVOBEING_BEHAVIOR_MAX_DAYS];
    int weekday[EVOBEING_BEHAVIOR_MAX_DAYS];
    evobeing_behavior_event_t ev[EVOBEING_BEHAVIOR_MAX_EVENTS];
    int n_ev;
    int n_days;
} world_t;

/* 真实效应（分）：饮酒 −9，晚间咖啡 −5（白天咖啡 0），小睡 +4，运动 +3，冥想 +2 */
static void gen_world(world_t *w, int n_days, double p_alc, double p_caf, double p_nap, double p_ex, double p_med,
                      double noise, double weekend_eff, int late_caffeine_only)
{
    int d;
    double prev = 0.0;
    w->n_days = n_days; w->n_ev = 0;
    for (d = 0; d < n_days; ++d) {
        double eff = 0.0;
        int wd = d % 7;
        w->weekday[d] = wd;
        w->wake[d] = d * 24.0 + 7.0 + 0.5 * nrand();
        /* 行为发生在第 d-1 天起床之后、第 d 天起床之前 → 影响 y[d]。用第 d-1 天的时刻。 */
        if (d > 0) {
            double base = (d - 1) * 24.0;
            if (urand() < p_alc) { w->ev[w->n_ev].id = EVOBEING_BEHAVIOR_ALCOHOL; w->ev[w->n_ev++].t_h = base + 20.0 + 2.0 * urand(); eff -= 9.0; }
            if (urand() < p_caf) {
                double hr = late_caffeine_only ? 15.0 + 3.0 * urand() : (urand() < 0.5 ? 9.0 + 3.0 * urand() : 15.0 + 3.0 * urand());
                w->ev[w->n_ev].id = EVOBEING_BEHAVIOR_CAFFEINE; w->ev[w->n_ev++].t_h = base + hr;
                if (hr >= 14.0) eff -= 5.0;
            }
            if (urand() < p_nap) { w->ev[w->n_ev].id = EVOBEING_BEHAVIOR_NAP; w->ev[w->n_ev++].t_h = base + 13.5; eff += 4.0; }
            if (urand() < p_ex) { w->ev[w->n_ev].id = EVOBEING_BEHAVIOR_EXERCISE; w->ev[w->n_ev++].t_h = base + 18.0; eff += 3.0; }
            if (urand() < p_med) { w->ev[w->n_ev].id = EVOBEING_BEHAVIOR_MEDITATION; w->ev[w->n_ev++].t_h = base + 22.0; eff += 2.0; }
        }
        if (wd == 5 || wd == 6) eff += weekend_eff;
        w->y[d] = 62.0 + eff + 0.3 * prev + noise * nrand();
        if (w->y[d] < 0) w->y[d] = 0;
        if (w->y[d] > 100) w->y[d] = 100;
        prev = w->y[d] - 62.0;
    }
}

static const evobeing_behavior_item_t *find(const evobeing_behavior_output_t *o, int id, int late)
{
    int i;
    for (i = 0; i < o->n_behaviors; ++i)
        if (o->items[i].id == id && o->items[i].is_late_variant == late) return &o->items[i];
    return NULL;
}

static void run(const world_t *w, const evobeing_behavior_config_t *cfg, evobeing_behavior_output_t *o, int *rc)
{
    evobeing_behavior_input_t in;
    memset(&in, 0, sizeof in);
    in.n_days = w->n_days; in.wake_score = w->y; in.wake_t_h = w->wake; in.weekday = w->weekday;
    in.events = w->ev; in.n_events = w->n_ev;
    *rc = evobeing_behavior_impact_zhenghuimin(&in, cfg, o);
}

static void print_out(const evobeing_behavior_output_t *o)
{
    int i;
    printf("  days=%d mean=%.1f sd=%.1f explained=%.0f%% unexplained=%.0f%%\n",
           o->n_days_used, o->mean_score, o->sd_score, o->explained_pct, o->unexplained_pct);
    for (i = 0; i < o->n_behaviors; ++i) {
        const evobeing_behavior_item_t *it = &o->items[i];
        if (it->status == EVOBEING_BEHAVIOR_INSUFFICIENT)
            printf("    %-16s n=%d/%d  再记录几次\n", it->name, it->n_with, it->n_without);
        else
            printf("    %-16s n=%d/%d  %+6.1f%% [%+.1f, %+.1f]  raw %+5.1f%%  share %4.1f%%  status %d\n",
                   it->name, it->n_with, it->n_without, it->effect_pct, it->ci_lo_pct, it->ci_hi_pct,
                   it->raw_pct, it->share_pct, (int)it->status);
    }
}

int main(void)
{
    world_t w;
    evobeing_behavior_output_t o;
    evobeing_behavior_config_t cfg;
    int rc, i;

    printf("[1] 典型 28 天：方向与大小\n");
    {
        int ok_alc = 0, ok_nap = 0, ok_caf_late = 0, ok_ex = 0, trials = 20, tr, g_alc = 0, g_nap = 0, g_ex = 0;
        double sum_alc = 0, sum_nap = 0, sum_ex = 0;
        for (tr = 0; tr < trials; ++tr) {
            const evobeing_behavior_item_t *a, *np, *cl, *ex;
            rng_s = 100u + (unsigned)tr;
            gen_world(&w, 28, 0.3, 0.6, 0.3, 0.4, 0.2, 4.0, 2.0, 0);
            run(&w, NULL, &o, &rc);
            CHECK(rc == 0, "rc");
            CHECK(o.n_days_used == 28, "n_days_used");
            a = find(&o, EVOBEING_BEHAVIOR_ALCOHOL, 0); np = find(&o, EVOBEING_BEHAVIOR_NAP, 0);
            cl = find(&o, EVOBEING_BEHAVIOR_CAFFEINE, 1); ex = find(&o, EVOBEING_BEHAVIOR_EXERCISE, 0);
            if (tr == 0) print_out(&o);
            if (a && a->status != EVOBEING_BEHAVIOR_INSUFFICIENT) { ++g_alc; sum_alc += a->effect_points; if (a->effect_points < -3.0) ++ok_alc; }
            if (np && np->status != EVOBEING_BEHAVIOR_INSUFFICIENT) { ++g_nap; sum_nap += np->effect_points; if (np->effect_points > 0.5) ++ok_nap; }
            if (ex && ex->status != EVOBEING_BEHAVIOR_INSUFFICIENT) { ++g_ex; sum_ex += ex->effect_points; if (ex->effect_points > 0.0) ++ok_ex; }
            if (cl && cl->status != EVOBEING_BEHAVIOR_INSUFFICIENT && cl->effect_points < -1.0) ++ok_caf_late;
            /* 排序：有数字的在前、按 |effect| 降序 */
            for (i = 1; i < o.n_behaviors; ++i) {
                const evobeing_behavior_item_t *p = &o.items[i - 1], *q = &o.items[i];
                if (p->status != EVOBEING_BEHAVIOR_INSUFFICIENT && q->status != EVOBEING_BEHAVIOR_INSUFFICIENT)
                    CHECK(fabs(p->effect_pct) >= fabs(q->effect_pct) - 1e-9, "sorted");
                if (p->status == EVOBEING_BEHAVIOR_INSUFFICIENT) CHECK(q->status == EVOBEING_BEHAVIOR_INSUFFICIENT, "insufficient last");
            }
            CHECK(o.explained_pct >= 0 && o.explained_pct <= 100.0 + 1e-9, "explained range");
            CHECK(fabs(o.explained_pct + o.unexplained_pct - 100.0) < 1e-9, "sum 100");
            {
                double s = 0; for (i = 0; i < o.n_behaviors; ++i) { CHECK(o.items[i].share_pct >= 0, "share>=0"); s += o.items[i].share_pct; }
                CHECK(fabs(s - o.explained_pct) < 1e-6, "share sum == explained");
            }
            /* 区间包含点估计 */
            for (i = 0; i < o.n_behaviors; ++i) {
                const evobeing_behavior_item_t *it = &o.items[i];
                if (it->status == EVOBEING_BEHAVIOR_INSUFFICIENT) continue;
                CHECK(it->ci_lo_pct <= it->effect_pct + 3.0 && it->ci_hi_pct >= it->effect_pct - 3.0, "ci around estimate");
                CHECK(it->n_with >= 4 && it->n_without >= 4, "gating");
            }
        }
        printf("  饮酒方向正确 %d/%d（均值 %.1f 分，真值 −9），小睡 %d/%d（%.1f，真值 +4），运动 %d/%d（%.1f，真值 +3），晚间咖啡 %d\n",
               ok_alc, g_alc, sum_alc / g_alc, ok_nap, g_nap, sum_nap / g_nap, ok_ex, g_ex, sum_ex / g_ex, ok_caf_late);
        CHECK(g_alc >= 15 && ok_alc >= g_alc - 1, "alcohol negative in nearly all trials given");
        CHECK(sum_alc / g_alc < -6.0 && sum_alc / g_alc > -12.0, "alcohol magnitude near −9");
        CHECK(sum_nap / g_nap > 2.0 && sum_nap / g_nap < 6.0, "nap magnitude near +4");
        CHECK(ok_nap >= g_nap * 3 / 4, "nap positive in most trials");
        CHECK(ok_ex >= g_ex * 3 / 4, "exercise positive in most trials");
    }

    printf("[2] 无效应的行为：平均接近 0、多数为 LOW/MID\n");
    {
        int tr, n_small = 0, n_high = 0, n_given = 0;
        double sum = 0;
        for (tr = 0; tr < 30; ++tr) {
            const evobeing_behavior_item_t *c;
            rng_s = 500u + (unsigned)tr;
            gen_world(&w, 28, 0.0, 0.6, 0.0, 0.0, 0.0, 4.0, 0.0, 0);
            /* 把咖啡全改到上午 → 真实效应 0 */
            for (i = 0; i < w.n_ev; ++i) { double base = floor(w.ev[i].t_h / 24.0) * 24.0; w.ev[i].t_h = base + 9.0; }
            for (i = 0; i < 28; ++i) w.y[i] = 62.0 + 4.0 * nrand();
            run(&w, NULL, &o, &rc);
            c = find(&o, EVOBEING_BEHAVIOR_CAFFEINE, 0);
            CHECK(c != NULL, "caffeine present");
            CHECK(find(&o, EVOBEING_BEHAVIOR_CAFFEINE, 1) == NULL, "no late variant when none late");
            if (c && c->status != EVOBEING_BEHAVIOR_INSUFFICIENT) {
                ++n_given; sum += c->effect_points;
                if (fabs(c->effect_points) < 3.0) ++n_small;
                if (c->status == EVOBEING_BEHAVIOR_HIGH) ++n_high;
            }
        }
        printf("  给出数字 %d 次，|效应|<3 分 %d 次，HIGH %d 次，均值 %.2f 分\n", n_given, n_small, n_high, n_given ? sum / n_given : 0.0);
        CHECK(n_given > 0 && fabs(sum / n_given) < 1.5, "null effect mean near 0");
        CHECK(n_high <= 3, "null effect rarely HIGH");
    }

    printf("[3] 样本不足：出现 < 4 次不给数字\n");
    {
        rng_s = 9u;
        gen_world(&w, 28, 0.0, 0.0, 0.0, 0.0, 0.0, 4.0, 0.0, 0);
        w.n_ev = 0;
        w.ev[w.n_ev].id = EVOBEING_BEHAVIOR_ALCOHOL; w.ev[w.n_ev++].t_h = 3 * 24 + 21.0;
        w.ev[w.n_ev].id = EVOBEING_BEHAVIOR_ALCOHOL; w.ev[w.n_ev++].t_h = 10 * 24 + 21.0;
        w.ev[w.n_ev].id = EVOBEING_BEHAVIOR_ALCOHOL; w.ev[w.n_ev++].t_h = 20 * 24 + 21.0;
        run(&w, NULL, &o, &rc);
        CHECK(rc == 0, "rc");
        CHECK(o.n_behaviors == 1, "one behavior");
        CHECK(o.items[0].status == EVOBEING_BEHAVIOR_INSUFFICIENT, "insufficient");
        CHECK(o.items[0].n_with == 3, "n_with 3");
        CHECK(o.items[0].effect_pct != o.items[0].effect_pct, "effect NAN");
        CHECK(o.items[0].share_pct == 0.0, "share 0");
        CHECK(o.n_insufficient == 1, "n_insufficient");
        CHECK(o.unexplained_pct == 100.0, "all unexplained");
        /* 每天都有 → 没有对照 → 同样不给 */
        w.n_ev = 0;
        for (i = 0; i < 27; ++i) { w.ev[w.n_ev].id = EVOBEING_BEHAVIOR_CAFFEINE; w.ev[w.n_ev++].t_h = i * 24 + 9.0; }
        run(&w, NULL, &o, &rc);
        CHECK(o.n_behaviors == 1 && o.items[0].status == EVOBEING_BEHAVIOR_INSUFFICIENT, "no control days");
        CHECK(o.items[0].n_without < 4, "n_without small");
    }

    printf("[4] 混杂控制：饮酒只在周五/六，周末本身也影响精力\n");
    {
        int tr, better = 0;
        for (tr = 0; tr < 20; ++tr) {
            const evobeing_behavior_item_t *a;
            double raw, adj;
            evobeing_behavior_config_t c2;
            rng_s = 900u + (unsigned)tr;
            gen_world(&w, 28, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0, +6.0, 0);  /* 周末早晨 +6 分 */
            w.n_ev = 0;
            for (i = 1; i < 28; ++i) {
                int wd = i % 7;
                if ((wd == 5 || wd == 6) && urand() < 0.9) {   /* 周五/六晚饮酒，真实效应 −4 */
                    w.ev[w.n_ev].id = EVOBEING_BEHAVIOR_ALCOHOL; w.ev[w.n_ev++].t_h = (i - 1) * 24 + 21.0;
                    w.y[i] -= 4.0;
                }
            }
            run(&w, NULL, &o, &rc);
            a = find(&o, EVOBEING_BEHAVIOR_ALCOHOL, 0);
            if (!a || a->status == EVOBEING_BEHAVIOR_INSUFFICIENT) continue;
            raw = a->raw_pct; adj = a->effect_pct;
            /* 不控制周末时估计会偏正（被周末的 +6 盖住） */
            evobeing_behavior_config_default(&c2); c2.control_weekend = 0;
            run(&w, &c2, &o, &rc);
            if (tr == 0) printf("  raw %+.1f%%  控制周末 %+.1f%%  不控制 %+.1f%%（真值约 −6%%）\n", raw, adj, find(&o, EVOBEING_BEHAVIOR_ALCOHOL, 0)->effect_pct);
            if (adj < find(&o, EVOBEING_BEHAVIOR_ALCOHOL, 0)->effect_pct) ++better;
        }
        printf("  控制周末后估计更偏负（更接近真值）：%d/20\n", better);
        CHECK(better >= 16, "weekend control helps");
    }

    printf("[5] 边界：空输入、缺天、超长输入、窗口截取、事件归属\n");
    {
        evobeing_behavior_input_t in;
        memset(&in, 0, sizeof in);
        CHECK(evobeing_behavior_impact_zhenghuimin(NULL, NULL, &o) == -1, "null in");
        in.n_days = 0; in.wake_score = w.y;
        CHECK(evobeing_behavior_impact_zhenghuimin(&in, NULL, &o) == -1, "n_days 0");
        in.n_days = EVOBEING_BEHAVIOR_MAX_DAYS + 1;
        CHECK(evobeing_behavior_impact_zhenghuimin(&in, NULL, &o) == -1, "too many days");
        in.n_days = 28; in.wake_score = NULL;
        CHECK(evobeing_behavior_impact_zhenghuimin(&in, NULL, &o) == -1, "null scores");

        rng_s = 3u;
        gen_world(&w, 28, 0.3, 0.5, 0.3, 0.3, 0.0, 4.0, 0.0, 0);
        for (i = 0; i < 28; ++i) if (i % 4 == 0) w.y[i] = NAN;   /* 7 天缺晨起值 */
        run(&w, NULL, &o, &rc);
        CHECK(rc == 0, "rc with gaps");
        CHECK(o.n_days_used == 21, "21 valid days");
        for (i = 0; i < 28; ++i) if (i >= 5) w.y[i] = NAN;      /* 只剩 5 天 */
        run(&w, NULL, &o, &rc);
        CHECK(rc == -2, "too few days");
        CHECK(o.n_days_used == 3, "n_days_used = 3 (days 1..3; day 0 and 4 also NAN)");
        CHECK(o.unexplained_pct == 100.0, "unexplained 100");

        /* 56 天输入只用最后 28 天：前 28 天塞满饮酒、后 28 天没有 → 饮酒应该不出现或 insufficient */
        rng_s = 11u;
        gen_world(&w, 56, 0.0, 0.0, 0.3, 0.0, 0.0, 4.0, 0.0, 0);
        for (i = 1; i < 28; ++i) { w.ev[w.n_ev].id = EVOBEING_BEHAVIOR_ALCOHOL; w.ev[w.n_ev++].t_h = (i - 1) * 24 + 21.0; }
        run(&w, NULL, &o, &rc);
        CHECK(rc == 0, "rc 56");
        CHECK(o.n_days_used == 28, "window 28");
        {
            const evobeing_behavior_item_t *a = find(&o, EVOBEING_BEHAVIOR_ALCOHOL, 0);
            CHECK(a == NULL || a->status == EVOBEING_BEHAVIOR_INSUFFICIENT, "old events outside window ignored");
        }

        /* 事件归属：晚上 23:00 和凌晨 01:00 的饮酒都归到下一次起床 */
        rng_s = 12u;
        gen_world(&w, 28, 0.0, 0.0, 0.0, 0.0, 0.0, 2.0, 0.0, 0);
        w.n_ev = 0;
        for (i = 1; i < 28; i += 3) {
            w.ev[w.n_ev].id = EVOBEING_BEHAVIOR_ALCOHOL;
            w.ev[w.n_ev++].t_h = (i % 2) ? (i - 1) * 24 + 23.0 : i * 24 + 1.0;
            w.y[i] -= 10.0;
        }
        run(&w, NULL, &o, &rc);
        {
            const evobeing_behavior_item_t *a = find(&o, EVOBEING_BEHAVIOR_ALCOHOL, 0);
            CHECK(a && a->n_with == 9, "9 nights with alcohol");
            CHECK(a && a->effect_points < -6.0, "effect recovered");
            CHECK(a && a->status == EVOBEING_BEHAVIOR_HIGH, "clear effect → HIGH");
            CHECK(a && a->share_pct > 30.0, "large share");
        }
        /* 最后一次起床之后的事件不计入 */
        w.ev[w.n_ev].id = EVOBEING_BEHAVIOR_NAP; w.ev[w.n_ev++].t_h = 27 * 24 + 14.0;
        run(&w, NULL, &o, &rc);
        CHECK(find(&o, EVOBEING_BEHAVIOR_NAP, 0) == NULL, "event after last wake ignored");
    }

    printf("[6] 自定义行为与配置项\n");
    {
        rng_s = 21u;
        gen_world(&w, 28, 0.0, 0.0, 0.0, 0.0, 0.0, 3.0, 0.0, 0);
        w.n_ev = 0;
        for (i = 1; i < 28; i += 2) { w.ev[w.n_ev].id = EVOBEING_BEHAVIOR_CUSTOM_BASE + 1; w.ev[w.n_ev++].t_h = (i - 1) * 24 + 20.0; w.y[i] += 6.0; }
        evobeing_behavior_config_default(&cfg);
        cfg.custom_name[1] = "sauna";
        cfg.n_boot = 0;
        run(&w, &cfg, &o, &rc);
        CHECK(rc == 0 && o.n_behaviors == 1, "custom behavior");
        CHECK(!strcmp(o.items[0].name, "sauna"), "custom name");
        CHECK(o.items[0].effect_points > 3.0, "custom effect");
        CHECK(o.items[0].ci_lo_pct != o.items[0].ci_lo_pct, "no CI when n_boot=0");
        CHECK(o.items[0].status == EVOBEING_BEHAVIOR_MID, "n_boot=0 → MID by count");
        cfg.custom_name[1] = NULL;
        run(&w, &cfg, &o, &rc);
        CHECK(!strcmp(o.items[0].name, "custom_17"), "default custom name");
        /* 同一输入、同一种子 → 结果完全一致（可复现） */
        {
            evobeing_behavior_output_t o2;
            evobeing_behavior_config_default(&cfg);
            run(&w, &cfg, &o, &rc);
            run(&w, &cfg, &o2, &rc);
            CHECK(memcmp(&o, &o2, sizeof o) == 0, "deterministic");
        }
        /* 无 wake_t_h / 无 weekday 也能跑 */
        {
            evobeing_behavior_input_t in;
            memset(&in, 0, sizeof in);
            in.n_days = 28; in.wake_score = w.y; in.events = w.ev; in.n_events = w.n_ev;
            CHECK(evobeing_behavior_impact_zhenghuimin(&in, NULL, &o) == 0, "no wake/weekday");
            CHECK(o.n_behaviors == 1 && o.items[0].effect_points > 3.0, "default 7:00 wake works");
        }
        CHECK(!strcmp(evobeing_behavior_name_zhenghuimin(EVOBEING_BEHAVIOR_NAP), "nap"), "name");
    }

    printf("[7] 晚间咖啡：14 点后与上午分开\n");
    {
        int tr, ok = 0;
        for (tr = 0; tr < 10; ++tr) {
            const evobeing_behavior_item_t *c0, *c1;
            rng_s = 1300u + (unsigned)tr;
            gen_world(&w, 28, 0.0, 0.8, 0.0, 0.0, 0.0, 3.0, 0.0, 0);
            run(&w, NULL, &o, &rc);
            c0 = find(&o, EVOBEING_BEHAVIOR_CAFFEINE, 0); c1 = find(&o, EVOBEING_BEHAVIOR_CAFFEINE, 1);
            if (tr == 0) print_out(&o);
            if (c0 && c1 && c1->status != EVOBEING_BEHAVIOR_INSUFFICIENT && c0->status != EVOBEING_BEHAVIOR_INSUFFICIENT &&
                c1->effect_points < c0->effect_points - 1.0) ++ok;
        }
        printf("  晚间咖啡比总体咖啡更负：%d/10\n", ok);
        CHECK(ok >= 8, "late caffeine separated");
    }

    printf("%d checks, %d failed\n", n_checks, n_fail);
    return n_fail ? 1 : 0;
}
