/*
 * evobeing_behavior_zhenghuimin.h — 行为对晨起精力的影响（个体内归因）
 *
 * 独立的离线分析模块：输入过去 28 天每天的晨起精力值和行为记录，
 * 输出每种行为对晨起精力的影响百分比。只跟本人比，不跟别人比。
 * 不依赖、也不改动 evobeing_energy_zhenghuimin / evobeing_vitality_zhenghuimin，
 * 不使用 load / activity 之类的输入。建议放在手机端或云端运行（每天起床后跑一次）。
 *
 * 目标量
 * ------
 *   y_d = 第 d 天起床时刻精力引擎的初始输出，即晨起 WAKE 点的 score（0–100，
 *         vitality_daily_t.wake_score）。一天只有一个值；当天没有主睡眠时为 NAN、不参与分析。
 *
 * 行为记录
 * --------
 *   每条记录 = (行为 id, 发生时刻)。时刻用"自窗口第 0 天 00:00 起的小时数"表示。
 *   一条记录归到它之后的第一次起床：起床之后、下一次起床之前发生的行为，影响下一次起床的 y。
 *   行为分两类，接口上不区分：
 *     - 系统自动识别的：小睡、运动、晚间运动、过载日、晚睡（比本人常态晚）等，
 *       由调用方从每日汇总（vitality_daily_t / energy_daily_t）或自己的检测结果生成记录；
 *     - 用户手动记录的：饮酒、咖啡、冥想、压力事件等。
 *   自带 id 见 EVOBEING_BEHAVIOR_*；EVOBEING_BEHAVIOR_CUSTOM_BASE 起可自定义。
 *
 * 百分比的含义（口径 A：影响幅度）
 * -------------------------------
 *   effect_pct = 100 · b_j / ȳ
 *   b_j 是"有这个行为的夜晚，次晨精力比没有时高 / 低多少分"（已控制其他行为、周末、
 *   前一天的精力和近期趋势之后），ȳ 是本人 28 天晨起精力的加权均值。
 *   例如 −12% 表示"饮酒后的早晨，精力比你平时低约 12%"。各项之和不必为 100%，
 *   也不必为 0。另外给出 share_pct（口径 B：解释份额），各项之和 ≤ 100%，
 *   仅用于排序与可信度判断，不建议直接展示。
 *
 * 可信度
 * ------
 *   28 天内行为至少出现 min_with 次、且至少 min_without 天没出现，才给出数字，否则
 *   status = INSUFFICIENT（"再记录几次就能看到"）。给出数字时按自助法区间宽度分三档：
 *   HIGH（区间不跨 0 且 n 够多）、MID、LOW。越近的日子权重越高（半衰期 half_life_days）。
 *
 * 用法
 * ----
 *   evobeing_behavior_input_t in = {0};
 *   in.n_days = 28;  in.wake_score = y;  in.weekday = wd;   // y[28]，NAN 表示当天没有
 *   in.events = ev;  in.n_events = m;                       // {id, t_h}
 *   in.wake_t_h = wake;                                      // 每天起床时刻（小时，NAN 表示没有）
 *   evobeing_behavior_output_t out;
 *   evobeing_behavior_impact_zhenghuimin(&in, NULL, &out);
 *   for (i = 0; i < out.n_behaviors; ++i)
 *       if (out.items[i].status != EVOBEING_BEHAVIOR_INSUFFICIENT)
 *           printf("%s %+.0f%%\n", out.items[i].name, out.items[i].effect_pct);
 */
#ifndef EVOBEING_BEHAVIOR_ZHENGHUIMIN_H
#define EVOBEING_BEHAVIOR_ZHENGHUIMIN_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

#define EVOBEING_BEHAVIOR_MAX_DAYS   56   /* 输入最多接受的天数（分析窗口默认取最近 28 天） */
#define EVOBEING_BEHAVIOR_MAX_TYPES  24   /* 同时分析的行为种类上限（含自定义） */
#define EVOBEING_BEHAVIOR_MAX_EVENTS 2048

/* 自带的行为 id。调用方可从 CUSTOM_BASE 起自定义（需在 config.custom 里给名字）。 */
enum {
    EVOBEING_BEHAVIOR_ALCOHOL = 0,        /* 用户记录：饮酒 */
    EVOBEING_BEHAVIOR_CAFFEINE = 1,       /* 用户记录：咖啡 / 茶 */
    EVOBEING_BEHAVIOR_MEDITATION = 2,     /* 用户记录：冥想 */
    EVOBEING_BEHAVIOR_STRESSOR = 3,       /* 用户记录：压力事件 */
    EVOBEING_BEHAVIOR_NAP = 4,            /* 自动：小睡（VITALITY_FLAG_NAP 或用户说刚睡了一觉） */
    EVOBEING_BEHAVIOR_EXERCISE = 5,       /* 自动：运动（VITALITY_FLAG_EXERCISE） */
    EVOBEING_BEHAVIOR_OVERLOAD = 6,       /* 自动：过载日（当天 STRESS 样本累计 ≥ 调用方阈值） */
    EVOBEING_BEHAVIOR_LATE_BED = 7,       /* 自动：入睡比本人常态晚 ≥ 1 小时 */
    EVOBEING_BEHAVIOR_LATE_MEAL = 8,      /* 用户记录：晚餐很晚 / 睡前进食 */
    EVOBEING_BEHAVIOR_SCREEN = 9,         /* 用户记录：睡前长时间看屏幕 */
    EVOBEING_BEHAVIOR_CUSTOM_BASE = 16    /* 16..23 自定义 */
};

typedef struct {
    int id;        /* 行为 id */
    double t_h;    /* 发生时刻：自窗口第 0 天 00:00 起的小时数（如第 3 天 21:30 → 3·24+21.5） */
} evobeing_behavior_event_t;

typedef struct {
    int n_days;                   /* 天数 ≤ MAX_DAYS，第 0 天最早 */
    const double *wake_score;     /* [n_days] 每天晨起精力 0–100；无则 NAN */
    const double *wake_t_h;       /* [n_days] 每天起床时刻（小时数）；NULL 则按每天 7:00 处理；NAN 表示当天没起床记录 */
    const int *weekday;           /* [n_days] 0=周一 … 6=周日；NULL 则不做周末控制 */
    const evobeing_behavior_event_t *events;
    int n_events;
} evobeing_behavior_input_t;

/* 时间点细分：某行为在 late_hour（当地时刻）之后发生，另算一个"晚间 × 行为"子项。 */
typedef struct {
    int id;
    double late_hour;   /* NAN = 不细分。默认：咖啡 14:00、运动 20:00、饮酒 21:00 */
} evobeing_behavior_late_rule_t;

typedef struct {
    int window_days;      /* 分析窗口，28；取输入最后 window_days 天 */
    double half_life_days;/* 近期权重半衰期，14：14 天前的一天权重 0.5 */
    int min_with;         /* 行为至少出现次数，4 */
    int min_without;      /* 至少几天没出现，4 */
    double ridge;         /* 岭回归收缩强度：相当于加 ridge 个"效应为 0 的伪日"，2.0 */
    int n_boot;           /* 自助法重采样次数，200；0 = 不算区间 */
    unsigned seed;        /* 自助法随机种子，默认 12345 */
    int control_weekend;  /* 1：加入周五 / 周六晚的指示列 */
    int control_prev;     /* 1：加入前一天晨起精力（自相关） */
    int control_trend;    /* 1：加入线性趋势列 */
    double high_ci_ratio; /* |效应| / 半区间 ≥ 该值记 HIGH，2.0；≥ 1.0 记 MID；否则 LOW */
    evobeing_behavior_late_rule_t late[EVOBEING_BEHAVIOR_MAX_TYPES];
    int n_late;
    const char *custom_name[8];  /* 自定义 id 16..23 的名字；NULL 用 "custom_16" 等 */
} evobeing_behavior_config_t;

typedef enum {
    EVOBEING_BEHAVIOR_INSUFFICIENT = 0, /* 样本不够，不给数字 */
    EVOBEING_BEHAVIOR_LOW = 1,
    EVOBEING_BEHAVIOR_MID = 2,
    EVOBEING_BEHAVIOR_HIGH = 3
} evobeing_behavior_status_t;

typedef struct {
    int id;
    int is_late_variant;     /* 1：这是"晚间 × 行为"子项（晚于 late_hour 发生） */
    char name[32];
    int n_with, n_without;   /* 窗口内有效天中，有 / 无该行为的夜晚数 */
    double effect_pct;       /* 口径 A：次晨精力相对本人均值的变化 %（控制混杂后）；INSUFFICIENT 时 NAN */
    double effect_points;    /* 同上，单位是分（0–100 刻度） */
    double raw_pct;          /* 未控制混杂的简单对比：有 − 无 的加权均值差 / ȳ·100 */
    double ci_lo_pct, ci_hi_pct; /* 自助法 90% 区间（n_boot=0 时 NAN） */
    double share_pct;        /* 口径 B：解释份额 %，各项之和 ≤ 100；INSUFFICIENT 时 0 */
    evobeing_behavior_status_t status;
} evobeing_behavior_item_t;

typedef struct {
    int n_days_used;         /* 窗口内有晨起值的天数 */
    double mean_score;       /* ȳ：加权均值 */
    double sd_score;         /* 加权标准差 */
    double explained_pct;    /* 全部行为列合计解释的份额 %（口径 B，= Σ share_pct） */
    double unexplained_pct;  /* 100 − explained_pct：睡眠本身、天气、身体状况等未记录因素 */
    int n_behaviors;
    evobeing_behavior_item_t items[EVOBEING_BEHAVIOR_MAX_TYPES * 2]; /* 按 |effect_pct| 降序；INSUFFICIENT 排最后 */
    int n_insufficient;
} evobeing_behavior_output_t;

/* 默认配置 */
void evobeing_behavior_config_default(evobeing_behavior_config_t *c);

/*
 * 计算行为对晨起精力的影响。
 *   in   输入（见上）
 *   cfg  NULL 用默认配置
 *   out  结果
 * 返回 0 成功；-1 参数错误（n_days 超限 / 指针为空）；-2 有效天数不足 7 天（out 仍填好
 * n_days_used，items 全部 INSUFFICIENT）。
 */
int evobeing_behavior_impact_zhenghuimin(const evobeing_behavior_input_t *in,
                                         const evobeing_behavior_config_t *cfg,
                                         evobeing_behavior_output_t *out);

/* 自带 id 的名字（自定义 id 返回 "custom_N"）。 */
const char *evobeing_behavior_name_zhenghuimin(int id);

#ifdef __cplusplus
}
#endif
#endif
