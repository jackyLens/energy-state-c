/*
 * main.c —— 命令行入口：读 5 分钟面板 CSV → 算能量 → 写逐点 CSV + 打印每日汇总
 *
 * 用法：
 *   energy_cli panel.csv [--out out.csv] [--no-ekf] [--no-smooth] [--json]
 *
 * 输入 CSV（第一行是表头，列顺序无关）：
 *   t_h,tod,wake,light,deep,rem,hr,hrv_log,steps,kcal,stress,q_hr,q_hrv,spo2,spo2_min
 *   缺失值写 nan 或留空。
 *
 * 输出 CSV：t_h,tod,s,e,e_sd,energy,lower,upper,base,load,activity
 */

#include "energy.h"
#include "vitality.h"

#include <locale.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#define MAX_COLS 24
#define LINE_MAX_LEN 8192

typedef struct {
    char name[MAX_COLS][32];
    int n;
} header_t;

static int split_csv(char *line, char **fields, int max_fields)
{
    int n = 0;
    char *p = line;
    while (p && n < max_fields) {
        fields[n++] = p;
        char *comma = strchr(p, ',');
        if (!comma) break;
        *comma = '\0';
        p = comma + 1;
    }
    /* 去掉行尾 \r\n 与首尾空白 */
    for (int i = 0; i < n; i++) {
        char *s = fields[i];
        while (*s == ' ' || *s == '\t' || *s == '"') s++;
        char *e = s + strlen(s);
        while (e > s && (e[-1] == ' ' || e[-1] == '\t' || e[-1] == '\r' || e[-1] == '\n' || e[-1] == '"')) e--;
        *e = '\0';
        fields[i] = s;
    }
    return n;
}

static double to_double(const char *s, double fallback)
{
    if (!s || !*s) return fallback;
    if (strcmp(s, "nan") == 0 || strcmp(s, "NaN") == 0 || strcmp(s, "NA") == 0) return NAN;
    char *end = NULL;
    const double v = strtod(s, &end);
    if (end == s) return fallback;
    return v;
}

static int find_col(const header_t *h, const char *name)
{
    for (int i = 0; i < h->n; i++) {
        if (strcmp(h->name[i], name) == 0) return i;
    }
    return -1;
}

static int load_panel(const char *path, energy_sample_t **out, size_t *n_out)
{
    FILE *fp = fopen(path, "r");
    if (!fp) {
        fprintf(stderr, "打开失败: %s\n", path);
        return -1;
    }
    char line[LINE_MAX_LEN];
    if (!fgets(line, sizeof(line), fp)) {
        fclose(fp);
        return -1;
    }
    header_t h;
    h.n = 0;
    char *fields[MAX_COLS];
    int nf = split_csv(line, fields, MAX_COLS);
    for (int i = 0; i < nf; i++) {
        snprintf(h.name[h.n++], sizeof(h.name[0]), "%s", fields[i]);
    }
    /* 必须有表头；否则按固定列序处理 */
    int idx_t = find_col(&h, "t_h");
    if (idx_t < 0) {
        fprintf(stderr, "表头缺少 t_h 列\n");
        fclose(fp);
        return -1;
    }
    const int it = idx_t;
    const int itod = find_col(&h, "tod");
    const int iwake = find_col(&h, "wake");
    const int ilight = find_col(&h, "light");
    const int ideep = find_col(&h, "deep");
    const int irem = find_col(&h, "rem");
    const int ihr = find_col(&h, "hr");
    const int ihrv = find_col(&h, "hrv_log");
    const int isteps = find_col(&h, "steps");
    const int ikcal = find_col(&h, "kcal");
    const int istress = find_col(&h, "stress");
    const int iqhr = find_col(&h, "q_hr");
    const int iqhrv = find_col(&h, "q_hrv");
    const int ispo2 = find_col(&h, "spo2");
    const int ispo2min = find_col(&h, "spo2_min");

    size_t cap = 2048, n = 0;
    energy_sample_t *arr = (energy_sample_t *)malloc(cap * sizeof(energy_sample_t));
    if (!arr) {
        fclose(fp);
        return -1;
    }
    while (fgets(line, sizeof(line), fp)) {
        if (line[0] == '\n' || line[0] == '\r') continue;
        nf = split_csv(line, fields, MAX_COLS);
        if (nf <= it) continue;
        if (n == cap) {
            cap *= 2;
            energy_sample_t *tmp = (energy_sample_t *)realloc(arr, cap * sizeof(energy_sample_t));
            if (!tmp) {
                free(arr);
                fclose(fp);
                return -1;
            }
            arr = tmp;
        }
        energy_sample_t *s = &arr[n++];
        memset(s, 0, sizeof(*s));
        s->t_h = to_double(fields[it], NAN);
        s->tod = (itod >= 0 && itod < nf) ? to_double(fields[itod], NAN) : -1.0;
        s->stage[ENERGY_WAKE] = (iwake >= 0 && iwake < nf) ? to_double(fields[iwake], 0.0) : 0.0;
        s->stage[ENERGY_LIGHT] = (ilight >= 0 && ilight < nf) ? to_double(fields[ilight], 0.0) : 0.0;
        s->stage[ENERGY_DEEP] = (ideep >= 0 && ideep < nf) ? to_double(fields[ideep], 0.0) : 0.0;
        s->stage[ENERGY_REM] = (irem >= 0 && irem < nf) ? to_double(fields[irem], 0.0) : 0.0;
        s->hr = (ihr >= 0 && ihr < nf) ? to_double(fields[ihr], NAN) : NAN;
        s->hrv_log = (ihrv >= 0 && ihrv < nf) ? to_double(fields[ihrv], NAN) : NAN;
        s->steps = (isteps >= 0 && isteps < nf) ? to_double(fields[isteps], 0.0) : 0.0;
        s->kcal = (ikcal >= 0 && ikcal < nf) ? to_double(fields[ikcal], 0.0) : 0.0;
        s->stress = (istress >= 0 && istress < nf) ? to_double(fields[istress], 0.0) : 0.0;
        s->q_hr = (iqhr >= 0 && iqhr < nf) ? to_double(fields[iqhr], 1.0) : 1.0;
        s->q_hrv = (iqhrv >= 0 && iqhrv < nf) ? to_double(fields[iqhrv], 1.0) : 1.0;
        s->spo2 = (ispo2 >= 0 && ispo2 < nf) ? to_double(fields[ispo2], NAN) : NAN;
        s->spo2_min = (ispo2min >= 0 && ispo2min < nf) ? to_double(fields[ispo2min], NAN) : NAN;
    }
    fclose(fp);
    *out = arr;
    *n_out = n;
    return 0;
}

static int write_points(const char *path, const energy_output_t *out)
{
    FILE *fp = stdout;
    if (path) {
        fp = fopen(path, "w");
        if (!fp) {
            fprintf(stderr, "写入失败: %s\n", path);
            return -1;
        }
    }
    fprintf(fp, "t_h,tod,s,e,e_sd,energy,lower,upper,base,load,activity\n");
    for (size_t i = 0; i < out->n_points; i++) {
        const energy_point_t *p = &out->points[i];
        fprintf(fp, "%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g\n",
                p->t_h, p->tod, p->s, p->e, p->e_sd, p->energy, p->lower, p->upper,
                p->base, p->load, p->activity);
    }
    if (path) {
        const int rc = ferror(fp);
        fclose(fp);
        if (rc) {
            fprintf(stderr, "写入失败（磁盘/权限）: %s\n", path);
            return -1;
        }
    } else if (fflush(fp) != 0) {
        return -1;
    }
    return 0;
}

/* --vitality：在能量列之后追加精力列（同一行，同一时间戳） */
static int write_vpoints(const char *path, const vitality_output_t *vo)
{
    FILE *fp = stdout;
    if (path) {
        fp = fopen(path, "w");
        if (!fp) {
            fprintf(stderr, "写入失败: %s\n", path);
            return -1;
        }
    }
    fprintf(fp, "t_h,tod,s,e,e_sd,energy,lower,upper,base,load,activity,"
                "v_star,d,d_sd,vitality,percentile,vscore,vlower,vupper,quadrant,vflags\n");
    for (size_t i = 0; i < vo->n_points; i++) {
        const energy_point_t *p = &vo->energy.points[i];
        const vitality_point_t *v = &vo->points[i];
        fprintf(fp,
                "%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,"
                "%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%.17g,%d,%u\n",
                p->t_h, p->tod, p->s, p->e, p->e_sd, p->energy, p->lower, p->upper, p->base, p->load,
                p->activity, v->v_star, v->d, v->d_sd, v->vitality, v->percentile, v->score, v->lower,
                v->upper, (int)v->quadrant, v->flags);
    }
    if (path) {
        const int rc = ferror(fp);
        fclose(fp);
        if (rc) {
            fprintf(stderr, "写入失败（磁盘/权限）: %s\n", path);
            return -1;
        }
    } else if (fflush(fp) != 0) {
        return -1;
    }
    return 0;
}

static int run_vitality(const energy_sample_t *samples, size_t n, const energy_config_t *cfg,
                        const char *out_path)
{
    vitality_output_t vo;
    const int rc = vitality_compute(samples, n, cfg, NULL, &vo);
    if (rc != ENERGY_OK) {
        fprintf(stderr, "vitality_compute 失败: %d\n", rc);
        return 1;
    }
    if (write_vpoints(out_path, &vo) != 0) {
        vitality_output_free(&vo);
        return 1;
    }
    fprintf(stderr, "样本 %zu，天数 %zu，能量 %.1f–%.1f（均值 %.1f），精力 %.1f–%.1f（均值 %.1f）\n",
            vo.n_points, vo.n_daily, vo.energy.stats_min, vo.energy.stats_max, vo.energy.stats_mean,
            vo.stats_min, vo.stats_max, vo.stats_mean);
    fprintf(stderr, "\n日期   起床精力  D_init   均值   最低   最高   运动  午睡  过载点  用户标记  夜间SpO2  下跌占比\n");
    for (size_t i = 0; i < vo.n_daily; i++) {
        const vitality_daily_t *d = &vo.daily[i];
        fprintf(stderr, "%-6s %7.1f  %+6.3f  %5.1f  %5.1f  %5.1f  %4d  %4d  %6d  %8d  %8.1f  %8.3f\n",
                d->date[0] ? d->date : "-", d->wake_score, d->d_init, d->mean, d->min, d->max,
                d->n_exercise, d->n_nap, d->n_stress_samples, d->n_user_events, d->spo2_night, d->spo2_dip);
    }
    vitality_output_free(&vo);
    return 0;
}

int main(int argc, char **argv)
{
    setlocale(LC_NUMERIC, "C"); /* CSV 的解析/输出必须是 C locale，避免逗号小数 */
    const char *in_path = NULL;
    const char *out_path = NULL;
    int use_ekf = 1;
    int no_smooth = 0;
    int load_stats = 0;
    int vitality = 0;

    for (int i = 1; i < argc; i++) {
        if (strcmp(argv[i], "--out") == 0 && i + 1 < argc) {
            out_path = argv[++i];
        } else if (strcmp(argv[i], "--no-ekf") == 0) {
            use_ekf = 0;
        } else if (strcmp(argv[i], "--no-smooth") == 0) {
            no_smooth = 1;
        } else if (strcmp(argv[i], "--load-stats") == 0) {
            load_stats = 1;
        } else if (strcmp(argv[i], "--vitality") == 0) {
            vitality = 1;
        } else if (strcmp(argv[i], "--help") == 0 || strcmp(argv[i], "-h") == 0) {
            printf("用法: %s <panel.csv> [--out out.csv] [--no-ekf] [--no-smooth] [--load-stats] [--vitality]\n", argv[0]);
            return 0;
        } else if (!in_path) {
            in_path = argv[i];
        }
    }
    if (!in_path) {
        fprintf(stderr, "用法: %s <panel.csv> [--out out.csv] [--no-ekf] [--no-smooth] [--load-stats] [--vitality]\n", argv[0]);
        return 2;
    }

    energy_sample_t *samples = NULL;
    size_t n = 0;
    if (load_panel(in_path, &samples, &n) != 0 || n < 2) {
        fprintf(stderr, "读取面板失败或样本不足（%zu）\n", n);
        free(samples);
        return 1;
    }

    energy_config_t cfg;
    energy_config_default(&cfg);
    cfg.use_ekf = use_ekf;
    if (no_smooth) cfg.display.smooth_minutes = 0.0;
    cfg.start_year = 2026;
    cfg.start_month = 9;
    cfg.start_day = 1;

    if (load_stats) {
        energy_load_stats_t st;
        if (energy_load_calibrate(samples, n, &cfg.load, &st) == ENERGY_OK) {
            printf("hr_rest=%.17g hr_max=%.17g l_ref=%.17g\n", st.hr_rest, st.hr_max, st.l_ref);
        }
    }

    if (vitality) {
        const int vrc = run_vitality(samples, n, &cfg, out_path);
        free(samples);
        return vrc;
    }

    energy_output_t out;
    const int rc = energy_compute(samples, n, &cfg, &out);
    if (rc != ENERGY_OK) {
        fprintf(stderr, "energy_compute 失败: %d\n", rc);
        free(samples);
        return 1;
    }

    if ((!load_stats || out_path) && write_points(out_path, &out) != 0) {
        energy_output_free(&out);
        free(samples);
        return 1;
    }

    fprintf(stderr, "样本 %zu，天数 %zu，能量 %.1f–%.1f（均值 %.1f，带宽 %.1f）\n",
            out.n_points, out.n_daily, out.stats_min, out.stats_max, out.stats_mean,
            out.stats_band_mean);
    fprintf(stderr, "\n日期   起床   晚间   均值   最低   最高   峰值时刻  低谷时刻  睡眠债\n");
    for (size_t i = 0; i < out.n_daily; i++) {
        const energy_daily_t *d = &out.daily[i];
        fprintf(stderr, "%-6s %5.1f  %5.1f  %5.1f  %5.1f  %5.1f  %5.0f:%02.0f    %5.0f:%02.0f    %.2f\n",
                d->date[0] ? d->date : "-", d->morning, d->evening, d->mean, d->min, d->max,
                floor(d->peak_tod), (d->peak_tod - floor(d->peak_tod)) * 60.0,
                floor(d->dip_tod), (d->dip_tod - floor(d->dip_tod)) * 60.0,
                d->sleep_debt);
    }

    energy_output_free(&out);
    free(samples);
    return 0;
}
