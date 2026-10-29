/*
 * vitality.h —— 个体“精力”状态（能量之上的快层），v1.0
 *
 * 能量 E(t) 回答“我还剩多少储备”（电池电量，小时–天级）；
 * 精力 V(t) 回答“我现在状态好不好”（信号强度，分钟–小时级）。
 * 两者是单向级联：能量 → 精力，精力不反馈能量。
 *
 *     V*(t) = E(t) · [1 + ρ · c(t)]                      均衡线：能量 × 日内节律
 *     dD/dt = −D / τ_D(E) + Σ I_m(t)                    偏离量：均值回复 + 急性事件
 *     τ_D(E) = τ0 · [1 + λ (1 − E)]                      能量越低回弹越慢
 *     −V* ≤ D ≤ d_cap · E                                天花板：能量低时再刺激也拉不高
 *     V(t)  = clip(V* + D, 0, 1)   → 精力 0–100
 *
 * 观测（都是能量模型**没解释掉**的那部分，也就是 energy_point_t.hr_resid / hrv_resid）：
 *   - 日间心率残差 → 倒 U 形“唤醒度”变换 u(z)：轻度升高（约 +5–15 bpm，无运动）是“打起精神”，
 *     精力上升；持续大幅升高（> +20 bpm，无运动）是“压力过载”，精力下降；低于常态是“犯困”。
 *   - 日间 HRV 残差 → 线性：比常态高 = 副交感占优 = 状态好。
 *   两路都以个人观测噪声 σ 做 z 化，所以“高/低”只相对这个人自己；不同人的 V 不可比。
 *
 * 事件（从生理效应推断，不猜原因）：
 *   运动结束后的回弹（+）、小睡（+）、持续低 HRV 的疲劳拖拽（−，持续项）；
 *   用户手动标记（咖啡 / 冥想 / 压力事件 …）用 vitality_ctx_event 注入，锦上添花、不是必需。
 *
 * 晨起初始化：每次夜间睡眠结束，D 重置为
 *     D_init = g · [w1·x_ΔE + w2·x_HRV + w3·x_sleep]
 *   三项分别是“今晚能量回充量 / 醒前 1 小时 HRV / 睡眠结构（时长、深睡、碎片化）”相对
 *   **本人**常态的偏差（库内在线学习），因此睡得比平时好 → 醒来高于均衡线，反之低于。
 *
 * 0–100 的两种刻度（都只跟自己比）：
 *   - vitality   = 100·V：“绝对自我刻度”，E 与 D 本身已相对个人基线；
 *   - percentile = 过去约 28 天（指数加权）本人日间 V 的经验分位：“比你 80% 的时候好”。
 *   cfg.score_mode 决定 `score` 取哪一种（学习期分位不可用时自动退回 vitality）。
 *
 * 用法
 * ----
 *   流式（内部自带一个 energy_ctx_t，只需喂样本）：
 *     vitality_ctx_t vc;  vitality_ctx_init(&vc, &ecfg, &vcfg);
 *     for (each 5-min sample) {
 *         vitality_point_t vp;  energy_point_t ep;
 *         vitality_ctx_step(&vc, &sample, &ep, &vp);
 *         printf("能量 %.0f  精力 %.0f\n", ep.energy, vp.score);
 *     }
 *     用户标记：vitality_ctx_event(&vc, VITALITY_EVENT_CAFFEINE, 0.0);
 *
 *   批量：
 *     vitality_output_t out = {0};
 *     vitality_compute(samples, n, &ecfg, &vcfg, &out);   // 内含 energy_compute
 *     vitality_output_free(&out);
 *
 *   已有自己的 energy_ctx_t：每步拿到 energy_point_t 后调 vitality_ctx_feed。
 *
 * 保证与限制
 * ----------
 *   - V、V*、D 有界；无随机性；无全局状态；逐位可复现；
 *   - 单步计算 < 100 次浮点运算，vitality_ctx_t 约 3 KB（含 energy_ctx_t，以 sizeof 为准）；
 *   - 事件检测只依据生理效应（负荷、活动量、残差），不识别原因；咖啡只能通过心率间接体现；
 *   - 夜间（主睡眠）不做观测更新，V 跟随均衡线并带 VITALITY_FLAG_SLEEP；
 *   - 个人常态（夜间回充量、睡眠结构、HRV）至少 1 晚后才用于晨起初始化，3 晚内带 LEARNING 标志；
 *   - 分位刻度至少积累 3 个日间才可用。
 */
#ifndef VITALITY_STATE_H
#define VITALITY_STATE_H

#include "energy.h"

#ifdef __cplusplus
extern "C" {
#endif

#define VITALITY_VERSION "1.0.0"
const char *vitality_version(void);

/* ------------------------------------------------------------------ */
/* 1. 参数                                                             */
/* ------------------------------------------------------------------ */
/* 偏离量 D 的动力学 */
typedef struct {
    double tau0_h;      /* 能量充沛时的基础回归时间常数，0.75 h（45 min） */
    double lambda;      /* 能量低时回归变慢的放大系数，1.5 */
    double d_cap;       /* 正向偏离上限系数：D ≤ d_cap·E，0.30 */
    double q_d;         /* 过程噪声（每小时方差），0.02：与 τ_D 匹配，稳态 P≈q_d·τ/2，使 +10 bpm 持续 30 min 约抬 D 0.05–0.08 */
    double p0;          /* 初始 / 晨起方差，0.04 */
} vitality_dyn_t;

/* 观测模型：y_hr = a_hr·D + ε，y_hrv = a_hrv·D + ε；y 是 z 化残差经变换后的量 */
typedef struct {
    double sigma_hr;    /* 日间心率残差的 z 化尺度（bpm），默认取 energy 的 sigma_hr_day = 7 */
    double sigma_hrv;   /* 日间 HRV 残差的 z 化尺度（ln），默认 0.40 */
    double a_hr;        /* 唤醒度 u(z_hr) 对 D 的斜率，3.0 */
    double a_hrv;       /* z_hrv 对 D 的斜率，2.0 */
    double r_hr;        /* 心率通道噪声方差，1.0 */
    double r_hrv;       /* HRV 通道噪声方差，2.0（日间 HRV 稀疏且噪） */
    double z_c;         /* 倒 U 的零点：z_hr = z_c（≈ z_c·σ bpm）时唤醒收益归零，2.8 */
    double drowsy_gain; /* z_hr < 0 时的斜率（心率低于常态 = 犯困），0.6 */
    double y_min;       /* u(z) 的下限，−2.0 */
    double rest_activity; /* 活动量 A 超过此值观测开始降权，0.10；超过 act_zero 不用 */
    double act_zero;    /* 0.35 */
    double load_zero;   /* 负荷 L 超过此值且 A ≥ rest_activity 视为运动中、不用观测；L 高而 A 低只降权一半，0.60 */
    double gate_sigma;  /* 新息门控，4.0 */
} vitality_obs_t;

/* 事件检测与脉冲 */
typedef struct {
    double a_exercise;      /* 运动结束回弹幅度，0.12 */
    double exercise_min_h;  /* 运动段最短时长，10/60 h；幅度按 min(1, 时长/exercise_full_h) 缩放 */
    double exercise_full_h; /* 0.5 h */
    double load_on;         /* 进入运动段的负荷阈值，0.5 */
    double load_off;        /* 离开运动段的负荷阈值，0.10 */
    double a_nap;           /* 小睡回弹幅度，0.15 */
    double nap_min_h;       /* 小睡最短时长，15/60 h；幅度按 min(1, 时长/nap_full_h) 缩放 */
    double nap_full_h;      /* 20/60 h */
    double nap_max_h;       /* 超过此时长视为主睡眠而非小睡，3.0 h */
    double nap_z_hr;        /* 无分期时的小睡判据：z_hr ≤ 该值且 A≈0 持续 nap_min_h，−0.75 */
    double nap_tod_min, nap_tod_max; /* 小睡只在 [9, 20) 钟点内识别 */
    double a_fatigue;       /* 持续低 HRV 的拖拽速率（每小时、每单位 −z），0.03 */
    double fatigue_z;       /* 触发阈值：z_hrv ≤ −1.0 */
    double fatigue_min_h;   /* 持续 0.5 h 后生效 */
    double arousal_z_lo, arousal_z_hi; /* 标志位用：轻度唤醒 z 区间 [0.7, 2.1]；> hi 为过载 */
    double stress_min_h;    /* 过载需持续 10/60 h 才打标志 */
    /* 用户标记的默认幅度（vitality_ctx_event 传 amplitude<=0 时使用） */
    double a_caffeine;      /* +0.10，经 caffeine_tau_h 一阶起效 */
    double caffeine_tau_h;  /* 20/60 h */
    double a_meditation;    /* +0.08 */
    double a_stressor;      /* −0.10 */
    double a_alcohol;       /* −0.08 */
} vitality_event_t;

/* 晨起初始化 */
typedef struct {
    double gain;          /* D_init = gain·(w1·x1 + w2·x2 + w3·x3)，0.40；x 各在 [−0.5, 0.5]，|D_init| ≤ 0.2 */
    double w_delta_e;     /* 0.4：夜间能量回充量相对本人常态 */
    double w_hrv;         /* 0.3：醒前 1 小时 HRV 相对本人夜间常态 */
    double w_sleep;       /* 0.3：睡眠结构（时长 / 深睡占比 / 碎片化）相对本人常态 */
    double sigma_hrv_night; /* HRV z 化尺度，0.22 */
    double lambda;        /* 常态的 EWMA 系数，0.20（约 5 晚记忆）；前 cold_nights 晚累计平均 */
    int cold_nights;      /* 3 */
    double min_sleep_h;   /* 不足此时长的睡眠段不算“一晚”，不触发晨起初始化（按小睡处理），3.0 */
    double min_sleep_prob; /* 睡眠段判定：sleep_prob ≥ 0.5 */
} vitality_morning_t;

/* 刻度与展示 */
typedef struct {
    double rho;          /* 节律调制幅度，0.15 */
    double peak_hour;    /* 10.0 */
    double dip_hour;     /* 14.5 */
    double dip_depth;    /* 0.70 */
    double dip_width_h;  /* 1.6 */
    double smooth_minutes; /* 批量展示平滑，15 min；0 关闭（流式不平滑） */
    double band_z;       /* 1.96 */
    int score_mode;      /* 0：score = vitality（100·V）；1：score = percentile（学习期退回 vitality）；2：两者平均 */
    double pct_half_life_days; /* 分位直方图的半衰期，14 天（≈ 28 天有效窗口） */
    int pct_min_days;    /* 分位刻度至少需要的日间数，3 */
    double quadrant_e;   /* 象限划分阈值（0–1），0.5 */
    double quadrant_v;   /* 0.5 */
} vitality_display_t;

typedef struct {
    vitality_dyn_t dyn;
    vitality_obs_t obs;
    vitality_event_t event;
    vitality_morning_t morning;
    vitality_display_t display;
    int use_obs;         /* 1：日间残差参与校正；0：纯事件 + 节律开环 */
    int detect_events;   /* 1：自动事件检测；0：只接受手动标记 */
} vitality_config_t;

void vitality_dyn_default(vitality_dyn_t *c);
void vitality_obs_default(vitality_obs_t *c);
void vitality_event_default(vitality_event_t *c);
void vitality_morning_default(vitality_morning_t *c);
void vitality_display_default(vitality_display_t *c);
void vitality_config_default(vitality_config_t *c);

/* 校验配置；非法返回 ENERGY_EINVAL（errbuf 可为 NULL）。 */
int vitality_config_validate(const vitality_config_t *cfg, char *errbuf, size_t errlen);

/* ------------------------------------------------------------------ */
/* 2. 输出                                                             */
/* ------------------------------------------------------------------ */
#define VITALITY_FLAG_SLEEP     1u    /* 主睡眠中：不做观测更新，V 跟随均衡线 */
#define VITALITY_FLAG_LEARNING  2u    /* 个人常态仍在学习期（< cold_nights 晚） */
#define VITALITY_FLAG_GAP       4u    /* 数据缺口（来自能量层） */
#define VITALITY_FLAG_WAKE      8u    /* 本点发生了晨起初始化 */
#define VITALITY_FLAG_EXERCISE  16u   /* 本点注入了运动结束回弹 */
#define VITALITY_FLAG_NAP       32u   /* 本点注入了小睡回弹 */
#define VITALITY_FLAG_AROUSAL   64u   /* 轻度交感激活（心率比常态高 5–15 bpm，无运动） */
#define VITALITY_FLAG_STRESS    128u  /* 过载：心率持续高于常态 > 15 bpm 且无运动 */
#define VITALITY_FLAG_FATIGUE   256u  /* 持续低 HRV 的疲劳拖拽生效中 */
#define VITALITY_FLAG_USER      512u  /* 本点注入了用户标记事件 */
#define VITALITY_FLAG_NO_OBS    1024u /* 本点无可用观测（残差缺失 / 运动中 / 能量基线未建立） */
#define VITALITY_FLAG_PCT_READY 2048u /* 分位刻度可用（≥ pct_min_days 个日间） */

/* 四象限（产品文案用） */
typedef enum {
    VITALITY_Q_HIGH_HIGH = 0, /* 能量高 精力高：状态极佳 */
    VITALITY_Q_HIGH_LOW = 1,  /* 能量高 精力低：暂时低迷，稍事休息即可 */
    VITALITY_Q_LOW_HIGH = 2,  /* 能量低 精力高：状态不错但储备不足，别透支 */
    VITALITY_Q_LOW_LOW = 3    /* 能量低 精力低：建议休息 */
} vitality_quadrant_t;

typedef struct {
    double t_h, tod;
    double e;          /* 能量状态 E（透传） */
    double v_star;     /* 均衡线 [0,1] */
    double d;          /* 偏离量（本点观测更新后） */
    double d_init;     /* 晨起初始化值 D_init（仅 WAKE 点有效，其余 NAN） */
    double d_sd;       /* D 的后验标准差 */
    double v;          /* V = clip(V* + D) [0,1]，未平滑 */
    double vitality;   /* 100·V（批量模式已做 smooth_minutes 平滑） */
    double percentile; /* 本人历史分位 0–100；不可用时 NAN */
    double score;      /* 按 cfg.display.score_mode 选定的 0–100 分 */
    double lower, upper; /* 95% 带（0–100，vitality 刻度） */
    double energy;     /* 能量 0–100（透传，便于双轨展示） */
    double z_hr, z_hrv; /* 本点用于观测的 z 化残差（无则 NAN） */
    vitality_quadrant_t quadrant;
    unsigned flags;    /* VITALITY_FLAG_* */
} vitality_point_t;

/* 每日汇总（与 energy_daily_t 同样按当地日历日切分） */
typedef struct {
    int day;
    char date[11];
    double wake_score;   /* 晨起初始化后第一点的 score；无则 NAN */
    double d_init;       /* 当日晨起 D_init；无则 NAN */
    double mean, min, max; /* 日间（sleep_prob<0.5）score 的统计 */
    double peak_tod, dip_tod;
    int n_exercise, n_nap, n_stress_samples, n_user_events;
    int n;
} vitality_daily_t;

/* 个人常态（晨起初始化与分位刻度的依据；建议落盘到画像，下次传 cfg 沿用） */
typedef struct {
    double delta_e;     /* 典型夜间能量回充量 */
    double hrv_wake;    /* 典型醒前 1 小时 ln(RMSSD) */
    double sleep_h;     /* 典型主睡眠时长 */
    double deep_frac;   /* 典型深睡占比 */
    double wake_frac;   /* 典型夜间清醒占比（碎片化） */
    int n_nights;       /* 已学习夜数 */
} vitality_habit_t;

typedef struct {
    vitality_point_t *points;
    size_t n_points;
    vitality_daily_t *daily;
    size_t n_daily;
    energy_output_t energy;   /* 内部 energy_compute 的完整输出（由 vitality_output_free 释放） */
    vitality_habit_t habit;   /* 结束时学到的个人常态 */
    double stats_min, stats_max, stats_mean;
} vitality_output_t;

void vitality_output_free(vitality_output_t *out);

/* ------------------------------------------------------------------ */
/* 3. 批量                                                             */
/* ------------------------------------------------------------------ */
/* ecfg / vcfg 可为 NULL（用默认）。内部先 energy_compute，再逐点推精力。 */
int vitality_compute(const energy_sample_t *samples, size_t n,
                     const energy_config_t *ecfg, const vitality_config_t *vcfg,
                     vitality_output_t *out);

/* ------------------------------------------------------------------ */
/* 4. 流式                                                             */
/* ------------------------------------------------------------------ */
#define VITALITY_PCT_BINS 50

typedef enum {
    VITALITY_EVENT_CAFFEINE = 0,
    VITALITY_EVENT_MEDITATION = 1,
    VITALITY_EVENT_NAP = 2,      /* 用户说“刚睡了一觉”（自动检测没抓到时） */
    VITALITY_EVENT_STRESSOR = 3, /* 用户标记的压力事件 */
    VITALITY_EVENT_ALCOHOL = 4,
    VITALITY_EVENT_CUSTOM = 5    /* amplitude 必须给（可正可负） */
} vitality_event_type_t;

/* 字段视为内部实现细节 */
typedef struct {
    vitality_config_t cfg;
    energy_ctx_t ectx;        /* 内部能量上下文（vitality_ctx_feed 不使用它） */
    int own_energy;
    int config_error;
    int started;
    double t_h, tod;
    double t0, tod0;          /* 首样本时间与钟点（日历日切分） */
    double d, p;              /* 偏离量及其方差 */
    double e_last;            /* 最近一次能量 E */
    vitality_habit_t habit;
    /* 主睡眠段跟踪 */
    int in_sleep, sleep_n;
    double sleep_t0, sleep_e0;
    double sleep_stage_acc[4];
    double hrv_recent[12];    /* 睡眠中最近 1 小时的 HRV 环形缓冲 */
    int hrv_recent_n, hrv_recent_i;
    int awake_run;            /* 睡眠段内连续清醒样本数（用于结束判定） */
    double last_close_t;      /* 上次主睡眠结束时刻（续睡合并用），无则 NAN */
    int last_close_main;      /* 上次结束的是否为主睡眠（≥ min_sleep_h） */
    vitality_habit_t habit_before_close; /* 上次主睡眠更新常态前的快照（续睡时回滚重算） */
    int awake_between_n;      /* 两段之间清醒样本数（合并时计入碎片化） */
    /* 事件检测 */
    int ex_run, ex_cool, ex_idle; /* 运动段样本数 / 结束后已平静连续样本数 / 结束后总样本数 */
    double ex_h;              /* 运动段累计时长 */
    int nap_run; double nap_h; int nap_cool;
    double last_d_init;        /* 最近一次晨起初始化值（供输出） */
    int fat_run; double fat_h;
    int stress_run;
    double pend_amp, pend_tau_h; /* 延迟起效的用户脉冲 */
    unsigned pending_flags;
    /* 分位直方图（日间 V） */
    double pct_hist[VITALITY_PCT_BINS];
    double pct_w;             /* 总权重 */
    double pct_last_day;      /* 上次衰减的日序 */
    int pct_days;             /* 已积累的日间数 */
    int day_seen;             /* 当前日历日是否已记入 pct_days */
    double last_v;
} vitality_ctx_t;

/* ecfg / vcfg 可为 NULL。cfg 非法时后续 step 返回 ENERGY_EINVAL。 */
void vitality_ctx_init(vitality_ctx_t *ctx, const energy_config_t *ecfg, const vitality_config_t *vcfg);

/* 从画像恢复 / 取出个人常态（不影响当前 D）。 */
int vitality_ctx_set_habit(vitality_ctx_t *ctx, const vitality_habit_t *h);
void vitality_ctx_get_habit(const vitality_ctx_t *ctx, vitality_habit_t *out);

/*
 * 推进一步（内部先 energy_ctx_step）。ep / vp 都可为 NULL。
 * 返回 energy_ctx_step 的返回码；能量层出错时精力不推进。
 */
int vitality_ctx_step(vitality_ctx_t *ctx, const energy_sample_t *sample,
                      energy_point_t *ep, vitality_point_t *vp);

/*
 * 用外部算好的 energy_point_t 推进一步（自己管理 energy_ctx_t 时使用）。
 * 要求 ep 来自同一条样本流、时间严格递增。
 */
int vitality_ctx_feed(vitality_ctx_t *ctx, const energy_sample_t *sample,
                      const energy_point_t *ep, vitality_point_t *vp);

/*
 * 注入用户标记事件；在下一次 step 时生效（咖啡按一阶起效，其余即时）。
 * amplitude <= 0 用默认幅度（CUSTOM 必须给非零值，可为负）。非法返回 ENERGY_EINVAL。
 */
int vitality_ctx_event(vitality_ctx_t *ctx, vitality_event_type_t type, double amplitude);

/* 不推进时间，写出当前状态（tod<0 用上次钟点）。 */
void vitality_ctx_peek(const vitality_ctx_t *ctx, double tod, vitality_point_t *out);

/* ------------------------------------------------------------------ */
/* 5. 子函数（便于测试 / 自检）                                          */
/* ------------------------------------------------------------------ */
double vitality_tau(const vitality_dyn_t *c, double e);                 /* τ_D(E)，小时 */
double vitality_equilibrium(const vitality_display_t *c, double e, double tod); /* V*(E, tod) */
double vitality_arousal(const vitality_obs_t *c, double z_hr);          /* 倒 U 变换 u(z) */
double vitality_d_bounds(const vitality_dyn_t *c, double e, double v_star, double d); /* 裁剪 D */
vitality_quadrant_t vitality_quadrant(const vitality_display_t *c, double e, double v);

#ifdef __cplusplus
}
#endif

#endif /* VITALITY_STATE_H */
