/*
 * evobeing_behavior_zhenghuimin.c — 行为对晨起精力的影响（个体内归因）
 *
 * 方法概要（详见 docs/behavior_impact_design.md）：
 *   1. 取最近 window_days 天，y_d = 当天晨起精力，按近期加权（半衰期 half_life_days）。
 *   2. 每条行为记录归到它之后的第一次起床；按行为生成 0/1 列，另按时间点生成"晚间 × 行为"列。
 *   3. 样本够的行为列 + 混杂控制列（周末、前一天精力、趋势）做加权岭回归，
 *      系数 b_j 即"有该行为的夜晚次晨精力高 / 低多少分"，除以本人均值得到 effect_pct。
 *   4. 自助法重采样给出 90% 区间，据此定可信度档位。
 *   5. 口径 B（解释份额）用去一列的加权 R² 差。
 */
#include "evobeing_behavior_zhenghuimin.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAXD EVOBEING_BEHAVIOR_MAX_DAYS
#define MAXT EVOBEING_BEHAVIOR_MAX_TYPES
#define MAXCOL (MAXT * 2 + 4)   /* 行为列 + 晚间列 + 截距/周末/前一天/趋势 */
#define NAN_V (0.0 / 0.0)

/* ------------------------------------------------------------------ */
/* 配置与名字                                                           */
/* ------------------------------------------------------------------ */
void evobeing_behavior_config_default(evobeing_behavior_config_t *c)
{
    memset(c, 0, sizeof *c);
    c->window_days = 28;
    c->half_life_days = 14.0;
    c->min_with = 4;
    c->min_without = 4;
    c->ridge = 2.0;
    c->n_boot = 200;
    c->seed = 12345u;
    c->control_weekend = 1;
    c->control_prev = 1;
    c->control_trend = 1;
    c->high_ci_ratio = 2.0;
    c->late[0].id = EVOBEING_BEHAVIOR_CAFFEINE; c->late[0].late_hour = 14.0;
    c->late[1].id = EVOBEING_BEHAVIOR_EXERCISE; c->late[1].late_hour = 20.0;
    c->late[2].id = EVOBEING_BEHAVIOR_ALCOHOL;  c->late[2].late_hour = 21.0;
    c->n_late = 3;
}

const char *evobeing_behavior_name_zhenghuimin(int id)
{
    switch (id) {
    case EVOBEING_BEHAVIOR_ALCOHOL:    return "alcohol";
    case EVOBEING_BEHAVIOR_CAFFEINE:   return "caffeine";
    case EVOBEING_BEHAVIOR_MEDITATION: return "meditation";
    case EVOBEING_BEHAVIOR_STRESSOR:   return "stressor";
    case EVOBEING_BEHAVIOR_NAP:        return "nap";
    case EVOBEING_BEHAVIOR_EXERCISE:   return "exercise";
    case EVOBEING_BEHAVIOR_OVERLOAD:   return "overload";
    case EVOBEING_BEHAVIOR_LATE_BED:   return "late_bed";
    case EVOBEING_BEHAVIOR_LATE_MEAL:  return "late_meal";
    case EVOBEING_BEHAVIOR_SCREEN:     return "screen";
    default:                           return "custom";
    }
}

static void item_name(const evobeing_behavior_config_t *cfg, int id, int late, char *buf, size_t n)
{
    const char *base = NULL;
    char tmp[16];
    if (id >= EVOBEING_BEHAVIOR_CUSTOM_BASE && id < EVOBEING_BEHAVIOR_CUSTOM_BASE + 8)
        base = cfg->custom_name[id - EVOBEING_BEHAVIOR_CUSTOM_BASE];
    if (!base) base = evobeing_behavior_name_zhenghuimin(id);
    if (!strcmp(base, "custom")) {
        snprintf(tmp, sizeof tmp, "custom_%d", id);
        base = tmp;
    }
    snprintf(buf, n, late ? "%s_late" : "%s", base);
}

/* ------------------------------------------------------------------ */
/* 小型线性代数：加权岭回归，列已中心化                                   */
/* ------------------------------------------------------------------ */
/* 解 (XᵀWX + diag(pen)) b = XᵀWy。X: n×p（行主序，步长 MAXCOL），返回 0 成功。 */
static int ridge_solve(const double *X, const double *y, const double *w, const int *use, int n, int p,
                       const double *pen, double *b)
{
    static double A[MAXCOL][MAXCOL + 1];
    int i, j, k;
    for (i = 0; i < p; ++i) {
        for (j = 0; j <= p; ++j) A[i][j] = 0.0;
        A[i][i] = pen[i];
    }
    for (k = 0; k < n; ++k) {
        const double *xr;
        if (!use[k]) continue;
        xr = X + (size_t)k * MAXCOL;
        for (i = 0; i < p; ++i) {
            double wi = w[k] * xr[i];
            if (wi == 0.0) continue;
            for (j = i; j < p; ++j) A[i][j] += wi * xr[j];
            A[i][p] += wi * y[k];
        }
    }
    for (i = 0; i < p; ++i)
        for (j = 0; j < i; ++j) A[i][j] = A[j][i];
    /* 高斯消元（部分主元） */
    for (i = 0; i < p; ++i) {
        int piv = i;
        double best = fabs(A[i][i]);
        for (k = i + 1; k < p; ++k)
            if (fabs(A[k][i]) > best) { best = fabs(A[k][i]); piv = k; }
        if (best < 1e-12) return -1;
        if (piv != i)
            for (j = 0; j <= p; ++j) { double t = A[i][j]; A[i][j] = A[piv][j]; A[piv][j] = t; }
        for (k = i + 1; k < p; ++k) {
            double f = A[k][i] / A[i][i];
            if (f == 0.0) continue;
            for (j = i; j <= p; ++j) A[k][j] -= f * A[i][j];
        }
    }
    for (i = p - 1; i >= 0; --i) {
        double s = A[i][p];
        for (j = i + 1; j < p; ++j) s -= A[i][j] * b[j];
        b[i] = s / A[i][i];
    }
    return 0;
}

/* 加权 R²（y 已中心化）；skip ≥ 0 时该列系数强制为 0 并重新拟合 */
static double weighted_r2(const double *X, const double *y, const double *w, const int *use, int n, int p,
                          const double *pen, int skip, double *b_out)
{
    double b[MAXCOL], pen2[MAXCOL];
    double sse = 0.0, sst = 0.0;
    int i, k;
    for (i = 0; i < p; ++i) pen2[i] = pen[i];
    if (skip >= 0) pen2[skip] = 1e9;   /* 极大惩罚 ≈ 去掉该列 */
    if (ridge_solve(X, y, w, use, n, p, pen2, b)) return NAN_V;
    for (k = 0; k < n; ++k) {
        const double *xr;
        double yhat = 0.0, r;
        if (!use[k]) continue;
        xr = X + (size_t)k * MAXCOL;
        for (i = 0; i < p; ++i) yhat += xr[i] * b[i];
        r = y[k] - yhat;
        sse += w[k] * r * r;
        sst += w[k] * y[k] * y[k];
    }
    if (b_out) memcpy(b_out, b, sizeof(double) * (size_t)p);
    return sst > 0.0 ? 1.0 - sse / sst : 0.0;
}

static unsigned xorshift(unsigned *s)
{
    unsigned x = *s;
    x ^= x << 13; x ^= x >> 17; x ^= x << 5;
    return *s = x;
}

static int cmp_double(const void *a, const void *b)
{
    double x = *(const double *)a, y = *(const double *)b;
    return x < y ? -1 : x > y ? 1 : 0;
}

/* ------------------------------------------------------------------ */
/* 主函数                                                               */
/* ------------------------------------------------------------------ */
int evobeing_behavior_impact_zhenghuimin(const evobeing_behavior_input_t *in,
                                         const evobeing_behavior_config_t *cfg_in,
                                         evobeing_behavior_output_t *out)
{
    evobeing_behavior_config_t cfg;
    static double X[MAXD * MAXCOL];
    static double occ[MAXT][MAXD], late[MAXT][MAXD];
    static double boot_b[MAXCOL][512];
    double y[MAXD], yc[MAXD], w[MAXD], wake[MAXD], pen[MAXCOL], b[MAXCOL];
    int use[MAXD], has_type[MAXT], has_late[MAXT], col_of[MAXT], late_col_of[MAXT];
    int n, start, d, t, k, p = 0, n_valid = 0, i, nb;
    double wsum = 0.0, ybar = 0.0, ysd = 0.0, lh[MAXT];
    double r2_full, r2_ctrl;
    int ctrl_cols = 0;
    int n_boot;

    if (!in || !out || in->n_days <= 0 || in->n_days > MAXD || !in->wake_score ||
        (in->n_events > 0 && !in->events) || in->n_events > EVOBEING_BEHAVIOR_MAX_EVENTS)
        return -1;
    if (cfg_in) cfg = *cfg_in; else evobeing_behavior_config_default(&cfg);
    if (cfg.window_days <= 0 || cfg.window_days > MAXD) cfg.window_days = 28;
    if (cfg.half_life_days <= 0.0) cfg.half_life_days = 14.0;
    if (cfg.n_boot > 500) cfg.n_boot = 500;
    n_boot = cfg.n_boot;
    memset(out, 0, sizeof *out);

    /* 1. 窗口、目标量、权重 */
    start = in->n_days > cfg.window_days ? in->n_days - cfg.window_days : 0;
    n = in->n_days - start;
    for (d = 0; d < n; ++d) {
        int src = start + d;
        y[d] = in->wake_score[src];
        wake[d] = in->wake_t_h ? in->wake_t_h[src] : (double)src * 24.0 + 7.0;
        use[d] = (y[d] == y[d]) && (y[d] >= 0.0) && (y[d] <= 100.0) && (wake[d] == wake[d]);
        w[d] = use[d] ? pow(0.5, (double)(n - 1 - d) / cfg.half_life_days) : 0.0;
        if (use[d]) { ++n_valid; wsum += w[d]; ybar += w[d] * y[d]; }
    }
    out->n_days_used = n_valid;
    if (n_valid < 7 || wsum <= 0.0) {
        out->mean_score = n_valid ? ybar / (wsum > 0 ? wsum : 1.0) : NAN_V;
        out->unexplained_pct = 100.0;
        return -2;
    }
    ybar /= wsum;
    for (d = 0; d < n; ++d) {
        yc[d] = use[d] ? y[d] - ybar : 0.0;
        if (use[d]) ysd += w[d] * yc[d] * yc[d];
    }
    ysd = sqrt(ysd / wsum);
    out->mean_score = ybar;
    out->sd_score = ysd;

    /* 2. 行为记录 → 每夜 0/1 */
    memset(occ, 0, sizeof occ);
    memset(late, 0, sizeof late);
    memset(has_type, 0, sizeof has_type);
    for (t = 0; t < MAXT; ++t) lh[t] = NAN_V;
    for (i = 0; i < cfg.n_late && i < MAXT; ++i)
        if (cfg.late[i].id >= 0 && cfg.late[i].id < MAXT) lh[cfg.late[i].id] = cfg.late[i].late_hour;
    for (k = 0; k < in->n_events; ++k) {
        const evobeing_behavior_event_t *e = &in->events[k];
        double th = e->t_h, tod;
        if (e->id < 0 || e->id >= MAXT || th != th) continue;
        /* 归到之后第一次起床：找最小的 d，wake[d] > th（跳过没有起床记录的天） */
        for (d = 0; d < n; ++d)
            if (wake[d] == wake[d] && wake[d] > th) break;
        if (d >= n) continue;              /* 在最后一次起床之后，还没有对应的早晨 */
        if (d == 0 && th < wake[0] - 24.0) continue; /* 早于窗口起点一天以上 */
        occ[e->id][d] = 1.0;
        has_type[e->id] = 1;
        tod = th - 24.0 * floor(th / 24.0);
        if (lh[e->id] == lh[e->id] && tod >= lh[e->id]) late[e->id][d] = 1.0;
    }

    /* 3. 设计矩阵：先行为列（样本够的），再晚间列，再控制列 */
    memset(X, 0, sizeof X);
    for (t = 0; t < MAXT; ++t) {
        int nw = 0, nwo = 0, nl = 0;
        col_of[t] = -1; late_col_of[t] = -1; has_late[t] = 0;
        if (!has_type[t]) continue;
        for (d = 0; d < n; ++d) if (use[d]) { if (occ[t][d] > 0) ++nw; else ++nwo; if (late[t][d] > 0) ++nl; }
        if (nw >= cfg.min_with && nwo >= cfg.min_without) {
            col_of[t] = p++;
            /* 晚间列：晚间次数够、且不是全部都在晚间（否则与主列共线，没有额外信息） */
            if (lh[t] == lh[t] && nl >= cfg.min_with && nl < nw) { late_col_of[t] = p++; has_late[t] = 1; }
        }
    }
    nb = p;
    if (cfg.control_weekend && in->weekday) {
        for (d = 0; d < n; ++d) {
            int wd = in->weekday[start + d];
            X[(size_t)d * MAXCOL + (size_t)p] = (wd == 5 || wd == 6) ? 1.0 : 0.0;  /* 周六 / 周日早晨 */
        }
        ++p; ++ctrl_cols;
    }
    if (cfg.control_prev) {
        for (d = 0; d < n; ++d)
            X[(size_t)d * MAXCOL + (size_t)p] = (d > 0 && use[d - 1]) ? yc[d - 1] / (ysd > 1e-9 ? ysd : 1.0) : 0.0;
        ++p; ++ctrl_cols;
    }
    if (cfg.control_trend) {
        for (d = 0; d < n; ++d) X[(size_t)d * MAXCOL + (size_t)p] = ((double)d - 0.5 * (double)(n - 1)) / (double)n;
        ++p; ++ctrl_cols;
    }
    for (t = 0; t < MAXT; ++t) {
        if (col_of[t] >= 0) for (d = 0; d < n; ++d) X[(size_t)d * MAXCOL + (size_t)col_of[t]] = occ[t][d];
        if (late_col_of[t] >= 0) for (d = 0; d < n; ++d) X[(size_t)d * MAXCOL + (size_t)late_col_of[t]] = late[t][d];
    }
    /* 列中心化（加权），截距即 ȳ；惩罚 = ridge · 列加权方差 · Σw */
    for (i = 0; i < p; ++i) {
        double m = 0.0, v = 0.0;
        for (d = 0; d < n; ++d) if (use[d]) m += w[d] * X[(size_t)d * MAXCOL + (size_t)i];
        m /= wsum;
        for (d = 0; d < n; ++d) {
            double *x = &X[(size_t)d * MAXCOL + (size_t)i];
            if (use[d]) { *x -= m; v += w[d] * (*x) * (*x); } else *x = 0.0;
        }
        /* 相当于加 ridge 个"效应为 0 的伪日"：收缩轻微，只为防共线与极端值 */
        pen[i] = cfg.ridge * (v > 1e-9 ? v : 1e-9) / (double)n_valid + 1e-9;
    }

    /* 4. 拟合 */
    r2_full = p > 0 ? weighted_r2(X, yc, w, use, n, p, pen, -1, b) : 0.0;
    if (r2_full != r2_full) r2_full = 0.0;
    /* 只有控制列时的 R²：把全部行为列压掉 */
    {
        double pen2[MAXCOL];
        for (i = 0; i < p; ++i) pen2[i] = i < nb ? 1e9 : pen[i];
        r2_ctrl = (p > nb) ? weighted_r2(X, yc, w, use, n, p, pen2, -1, NULL) : 0.0;
        if (r2_ctrl != r2_ctrl) r2_ctrl = 0.0;
    }

    /* 5. 自助法 */
    if (n_boot > 0 && nb > 0) {
        unsigned s = cfg.seed ? cfg.seed : 12345u;
        static double Xb[MAXD * MAXCOL];
        double yb[MAXD], wb[MAXD], bb[MAXCOL];
        int useb[MAXD], idx[MAXD], nv = 0, r, done = 0;
        for (d = 0; d < n; ++d) if (use[d]) idx[nv++] = d;
        for (r = 0; r < n_boot; ++r) {
            for (k = 0; k < nv; ++k) {
                int src = idx[xorshift(&s) % (unsigned)nv];
                memcpy(&Xb[(size_t)k * MAXCOL], &X[(size_t)src * MAXCOL], sizeof(double) * (size_t)p);
                yb[k] = yc[src]; wb[k] = w[src]; useb[k] = 1;
            }
            if (ridge_solve(Xb, yb, wb, useb, nv, p, pen, bb)) continue;
            for (i = 0; i < p; ++i) boot_b[i][done] = bb[i];
            ++done;
        }
        n_boot = done;
    } else n_boot = 0;

    /* 6. 输出各项 */
    {
        double share_sum = 0.0;
        int m = 0;
        for (t = 0; t < MAXT; ++t) {
            int variant;
            if (!has_type[t]) continue;
            for (variant = 0; variant <= has_late[t]; ++variant) {
                evobeing_behavior_item_t *it = &out->items[m++];
                int c = variant ? late_col_of[t] : col_of[t];
                int nw = 0, nwo = 0;
                double with_sum = 0.0, with_w = 0.0, wo_sum = 0.0, wo_w = 0.0;
                const double *ind = variant ? late[t] : occ[t];
                it->id = t; it->is_late_variant = variant;
                item_name(&cfg, t, variant, it->name, sizeof it->name);
                for (d = 0; d < n; ++d) {
                    if (!use[d]) continue;
                    if (ind[d] > 0) { ++nw; with_sum += w[d] * y[d]; with_w += w[d]; }
                    else { ++nwo; wo_sum += w[d] * y[d]; wo_w += w[d]; }
                }
                it->n_with = nw; it->n_without = nwo;
                it->raw_pct = (with_w > 0 && wo_w > 0 && ybar > 0) ? 100.0 * (with_sum / with_w - wo_sum / wo_w) / ybar : NAN_V;
                it->effect_pct = it->effect_points = it->ci_lo_pct = it->ci_hi_pct = NAN_V;
                it->share_pct = 0.0;
                it->status = EVOBEING_BEHAVIOR_INSUFFICIENT;
                if (c < 0) continue;
                /* 晚间子项的效应 = 主列 + 晚间列（晚间发生时的总效应） */
                it->effect_points = variant ? b[col_of[t]] + b[c] : b[c];
                it->effect_pct = ybar > 0 ? 100.0 * it->effect_points / ybar : NAN_V;
                {
                    double r2d = weighted_r2(X, yc, w, use, n, p, pen, c, NULL);
                    double sh = (r2d == r2d) ? (r2_full - r2d) : 0.0;
                    it->share_pct = sh > 0 ? 100.0 * sh : 0.0;
                    share_sum += it->share_pct;
                }
                if (n_boot >= 20) {
                    double tmp[512];
                    int lo_i = (int)floor(0.05 * (n_boot - 1)), hi_i = (int)ceil(0.95 * (n_boot - 1));
                    for (k = 0; k < n_boot; ++k)
                        tmp[k] = variant ? boot_b[col_of[t]][k] + boot_b[c][k] : boot_b[c][k];
                    qsort(tmp, (size_t)n_boot, sizeof(double), cmp_double);
                    it->ci_lo_pct = 100.0 * tmp[lo_i] / ybar;
                    it->ci_hi_pct = 100.0 * tmp[hi_i] / ybar;
                    {
                        double half = 0.5 * (it->ci_hi_pct - it->ci_lo_pct);
                        double ratio = half > 1e-9 ? fabs(it->effect_pct) / half : 0.0;
                        if (ratio >= cfg.high_ci_ratio && nw >= 6) it->status = EVOBEING_BEHAVIOR_HIGH;
                        else if (ratio >= 1.0) it->status = EVOBEING_BEHAVIOR_MID;
                        else it->status = EVOBEING_BEHAVIOR_LOW;
                    }
                } else {
                    it->status = nw >= 8 ? EVOBEING_BEHAVIOR_MID : EVOBEING_BEHAVIOR_LOW;
                }
            }
        }
        out->n_behaviors = m;
        /* 份额之和不超过"行为合计解释量"，也不超过 100 */
        {
            double cap = 100.0 * (r2_full - r2_ctrl);
            if (cap < 0) cap = 0;
            if (cap > 100) cap = 100;
            if (share_sum > cap && share_sum > 0) {
                for (i = 0; i < m; ++i) out->items[i].share_pct *= cap / share_sum;
                share_sum = cap;
            }
        }
        out->explained_pct = share_sum;
        out->unexplained_pct = 100.0 - share_sum;
        for (i = 0; i < m; ++i) if (out->items[i].status == EVOBEING_BEHAVIOR_INSUFFICIENT) ++out->n_insufficient;
    }

    /* 7. 排序：有数字的按 |effect_pct| 降序，INSUFFICIENT 在最后 */
    for (i = 1; i < out->n_behaviors; ++i) {
        evobeing_behavior_item_t key = out->items[i];
        int j = i - 1;
        double kv = key.status == EVOBEING_BEHAVIOR_INSUFFICIENT ? -1.0 : fabs(key.effect_pct);
        while (j >= 0) {
            double jv = out->items[j].status == EVOBEING_BEHAVIOR_INSUFFICIENT ? -1.0 : fabs(out->items[j].effect_pct);
            if (jv >= kv) break;
            out->items[j + 1] = out->items[j];
            --j;
        }
        out->items[j + 1] = key;
    }
    (void)ctrl_cols;
    return 0;
}
