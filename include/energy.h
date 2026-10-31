/*
 * energy.h —— 能量状态引擎（ODE + EKF），evobeing 公式，作者 zhenghuimin
 *
 * 能量 E 表示身体的储备水平（0–100），只跟本人的常态比，不做人与人之间的比较。
 * 观测方程的截距是每个人自己的基线（夜间 / 日间心率与 HRV 的常态水平），由库在线学习：
 * 首晚引导 → 前 3 晚累计平均 → 之后 EWMA λ=0.1，异常夜降权。
 *
 * 状态推进使用方程的解析解（指数更新），支持可变步长：
 *     S(t+h) = S∞ + (S−S∞)·e^(−b·h)
 *     E(t+h) = E∞ + (E−E∞)·e^(−(α+β)·h)
 *
 * 两种用法：
 *   批量：energy_compute(samples, n, cfg, &out)          一次算一整段；
 *   流式：energy_ctx_init + energy_ctx_step               每 5 分钟推一次。
 * 一般调用方不需要直接使用本文件，请用 evobeing_zhenghuimin.h 里的两个主函数。
 *
 * 约定：
 *   - 状态 S/E 恒在 [0,1]，能量恒在 [0,100]；
 *   - 无全局状态，不同 energy_ctx_t 互不影响；无随机性，同一输入结果逐位可复现；
 *   - |t_h| ≤ 1e6 小时；样本数 ≤ 1e7；时间必须严格递增，否则返回 ENERGY_EINVAL / ENERGY_EDOMAIN；
 *   - 缺失观测用 NAN；
 *   - 数据缺口超过 cfg.gap.max_hold_h（默认 1 h）的部分按个人常规睡眠窗 / 清醒静息推进，
 *     输出带 ENERGY_FLAG_GAP；
 *   - 负荷 L 与活动量 A 由库根据心率、步数、热量计算，调用方不传。
 */
#ifndef ENERGY_STATE_H
#define ENERGY_STATE_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

#define ENERGY_VERSION "1.2.0"

/* 运行期版本字符串（用于日志/上报）。 */
const char *energy_version(void);

/* 返回码 */
typedef enum {
    ENERGY_OK = 0,
    ENERGY_EINVAL = -1, /* 参数/输入非法 */
    ENERGY_ENOMEM = -2, /* 内存不足 */
    ENERGY_EDOMAIN = -3 /* 时间轴不单调等 */
} energy_status_t;

/* 分期顺序是全局约定：wake / light / deep / rem */
enum { ENERGY_WAKE = 0, ENERGY_LIGHT = 1, ENERGY_DEEP = 2, ENERGY_REM = 3 };

/* ------------------------------------------------------------------ */
/* 1. 参数                                                             */
/* ------------------------------------------------------------------ */
/* 状态方程参数，单位 1/h（无量纲的已注明）。默认值为先验均值。 */
typedef struct {
    double kappa_w;      /* 清醒期 S 积累，0.100/h */
    double kappa_s[4];   /* 分期 S 消散：wake/light/deep/rem */
    double alpha_sleep;  /* 睡眠期基础充电，0.300/h */
    double alpha_rest;   /* 清醒静息恢复，0.020/h */
    double w_stage[4];   /* 分期恢复权重：0.15/1.00/1.25/0.95 */
    double eta;          /* φ(S) 幂指数，0.60 */
    double s_ref;        /* φ=1 的参考 S，0.50 */
    double s_floor;      /* φ 的平滑下限，0.10 */
    double phi_min;      /* φ 硬下限，0.30 */
    double phi_max;      /* φ 硬上限，2.00 */
    double beta_0;       /* 基础消耗，0.010/h */
    double beta_L;       /* 单位负荷消耗，0.150 */
    double beta_P;       /* 单位压力消耗，0.030 */
    double beta_S;       /* 单位睡眠债的清醒期消耗，0.035 */
} energy_params_t;

/* 负荷 L 与活动量 A 的构造参数 */
typedef struct {
    double b_trimp;     /* Banister 系数，1.92（男）/1.67（女） */
    double w_trimp;     /* 1.0 */
    double w_steps;     /* 0.6 */
    double w_dhr;       /* 0.4 */
    double trimp_scale; /* 10.0，把 TRIMP 提到与步数项同量级 */
    double steps_ref;   /* 500 步 / 5 min ≈ 100 步/分钟（普通快走约 500–600 步 / 5 min） */
    double l_ref;       /* >0 固定参考值（用户画像）；<=0 自动：批量用本窗口 P99.5，流式在线估计 */
    double hr_rest;     /* <=0 自动：深睡期心率中位数（流式：深睡/睡眠心率的滑动均值） */
    double hr_max;      /* <=0 自动：历史 P99.5，限幅 [150,210]（流式：衰减峰值） */
    double hr_base[24]; /* 按小时的个人心率基线（用户画像，来自 energy_load_calibrate）；
                         * 24 个值全部 >0 时使用，否则自动（批量：288 点钟点基线；流式：在线分箱） */
} energy_load_cfg_t;

/* 展示层：节律、平滑、置信带 */
typedef struct {
    double rho;          /* 节律调制幅度，0.10（±10%） */
    double peak_hour;    /* 主峰钟点，10.0 */
    double dip_hour;     /* 午后低谷钟点，14.5 */
    double dip_depth;    /* 低谷深度，0.70 */
    double dip_width_h;  /* 低谷宽度，1.6 h */
    double smooth_minutes; /* 展示平滑，30 min；0 关闭 */
    double band_z;       /* 置信带倍数，1.96 */
} energy_display_cfg_t;

/* EKF */
typedef struct {
    double q_s, q_e;       /* 过程噪声（每小时方差） */
    double p0_s, p0_e;     /* 初始方差 */
    double gate_sigma;     /* 新息门控，4.0 */
    double q_floor;        /* 质量下限，0.05 */
    double robust_nu;      /* Student-t 自由度，4.0；0 表示高斯 */
    int gate_max_reject;   /* 同一通道连续被门控拒绝 N 次后改为“降权接收”（防止模型长期失锁），
                            * 6（30 分钟）；0 = 永远硬拒绝 */
} energy_ekf_cfg_t;

/* 观测模型：夜间 [1,E,p_deep,p_rem,p_wake]，日间 [1,E,A,cos1,sin1,cos2,sin2]
 * 第 0 列（截距）仅在 cfg.baseline.enable=0 时使用；默认由个体基线替代。 */
typedef struct {
    double hr_night[5];
    double hr_day[7];
    double hrv_night[5];
    double hrv_day[7];
    double sigma_hr_night, sigma_hr_day;
    double sigma_hrv_night, sigma_hrv_day;
} energy_obs_t;

/*
 * 个体基线（“你的常态”）。
 *
 * level[k] 是通道 k 的个人观测截距：沿着同样输入下的开环参考轨迹 E_ref（只由 ODE 推进、
 * 从不被观测校正），去掉 E/分期/活动/节律项之后，这个人该通道的平均水平。
 * 观测预测 = level[k] + Σ 斜率×特征，于是 EKF 看到的只是“相对你自己常态的偏离”：
 * 心率比你平时高 → 恢复差；同一条曲线整体平移不改变能量。
 *
 * 学习方式（库内自动，批量/流式同一套代码）：
 *   - 一个“段” = 一晚（连续有分期的样本，只用睡眠概率 >=0.5 的点）或一个白天（静息样本，A < rest_activity）；
 *   - 首段：段内满 bootstrap_samples 个样本后用段内运行均值作临时基线（之前不做校正）；
 *   - 段结束时提交：前 cold_nights 段累计平均（λ=1/(n+1)），之后 EWMA（λ=lambda≈0.1，约 10 段记忆）；
 *   - 提交值偏离当前基线超过 anomaly[k]（生病/饮酒/出差）时，λ 乘 anomaly_weight 降权；
 *   - n[k] < cold_nights 时观测噪声 R 按 cold_r_inflate 线性放大（学习期少信观测）。
 * 同时学习个人常规睡眠窗（sleep_on/sleep_off）与典型分期构成（stage_mix），供缺口策略使用。
 */
enum {
    ENERGY_BL_HR_NIGHT = 0,  /* 夜间心率，bpm */
    ENERGY_BL_HRV_NIGHT = 1, /* 夜间 ln(RMSSD) */
    ENERGY_BL_HR_DAY = 2,    /* 日间心率，bpm */
    ENERGY_BL_HRV_DAY = 3    /* 日间 ln(RMSSD) */
};

typedef struct {
    double level[4];     /* 个人观测截距，下标见 ENERGY_BL_*；n[k]==0 时无意义 */
    int n[4];            /* 已提交的段数（晚 / 天）；0 = 尚未建立 */
    double sleep_on;     /* 常规入睡钟点（0–24），n_sleep>0 时有效 */
    double sleep_off;    /* 常规醒来钟点（0–24） */
    double stage_mix[4]; /* 典型一晚的分期构成 [wake, light, deep, rem] */
    int n_sleep;         /* 已学习的夜数（睡眠窗/分期构成） */
} energy_baseline_t;

typedef struct {
    int enable;              /* 1：使用个体基线（默认）；0：用 cfg.obs 群体截距（仅供对照） */
    int cold_nights;         /* 冷启动段数，3（第 1–3 晚累计平均 + R 放大） */
    double lambda;           /* 冷启动后的 EWMA 系数，0.10（约 10 晚记忆） */
    double cold_r_inflate;   /* 冷启动首段的 R 放大倍数，4.0；随 n 线性回到 1 */
    int bootstrap_samples;   /* 首段引导所需样本数，6（30 分钟） */
    int min_samples;         /* 一段至少多少个有效样本才提交，12 */
    double rest_activity;    /* 日间基线只用 A < 该值的静息样本，0.10 */
    double anomaly[4];       /* 异常段阈值：6 bpm / 0.25 / 8 bpm / 0.35 */
    double anomaly_weight;   /* 异常段的 λ 倍率，0.25（降权而不是完全冻结，持续变化最终会被跟上） */
    energy_baseline_t init;  /* 画像里保存的基线（n[k]>0 的通道 / n_sleep>0 的习惯直接沿用） */
} energy_baseline_cfg_t;

/* 缺口（未佩戴 / 没电 / 数据断流）策略 */
typedef struct {
    double max_hold_h;        /* 缺口前输入最多沿用多久，1.0 h；超出部分按下面的策略推进 */
    double sleep_start;       /* 尚未学到个人睡眠窗时的默认入睡钟点，23.0 */
    double sleep_end;         /* 默认醒来钟点，7.0 */
    double template_stage[4]; /* 尚未学到个人分期构成时的软睡眠模板 {0.10,0.50,0.20,0.20} */
} energy_gap_cfg_t;

typedef struct {
    energy_params_t params;
    energy_load_cfg_t load;
    energy_display_cfg_t display;
    energy_ekf_cfg_t ekf;
    energy_obs_t obs;
    energy_baseline_cfg_t baseline; /* 个体化观测模型 */
    energy_gap_cfg_t gap;           /* 缺口策略 */
    double stress_scale; /* 压力量纲除数：>0 固定（0–100 设备填 100）；0 自动：
                          * 批量看整段有限值最大值（>1.5 则 /100），流式一旦见到 >1.5 的值就切到 /100 */
    double x0_s, x0_e;  /* 状态初值 (0.45, 0.55) */
    int use_ekf;        /* 1: HR/HRV 参与校正；0: 纯 ODE 开环 */
    int day_update;     /* 1: 日间观测也更新状态 */
    int start_year, start_month, start_day; /* >0 时用于生成 MM-DD 日期标签 */
} energy_config_t;

void energy_params_default(energy_params_t *p);
void energy_load_default(energy_load_cfg_t *c);
void energy_display_default(energy_display_cfg_t *c);
void energy_ekf_default(energy_ekf_cfg_t *c);
void energy_obs_default(energy_obs_t *o);
void energy_baseline_default(energy_baseline_cfg_t *c);
void energy_gap_default(energy_gap_cfg_t *c);
void energy_config_default(energy_config_t *c);

/*
 * 校验配置。返回 ENERGY_OK 或 ENERGY_EINVAL；errbuf 可为 NULL。
 * 校验项：速率参数为正、概率/限幅参数在合法区间、EKF 噪声与初始方差为正、
 *         展示参数有限、观测噪声为正。energy_compute 会自动调用；
 *         流式用法在 energy_ctx_init 里检查，非法配置会让后续 step 返回 EINVAL。
 */
int energy_config_validate(const energy_config_t *cfg, char *errbuf, size_t errlen);

/* ------------------------------------------------------------------ */
/* 2. 数据                                                             */
/* ------------------------------------------------------------------ */
/*
 * 一个 5 分钟样本（批量接口传数组，流式接口一次传一个）。
 * 用 energy_sample_init() 初始化后再填字段。
 */
typedef struct {
    /* 时间：距起点的小时数，严格递增；5 分钟一格时相邻样本差 1/12。 */
    double t_h;

    /* 当地钟点，0–24（含小数，例如 13.5 表示 13:30）。
     * 用于昼夜节律、日间观测窗口（07–22 点）与心率基线按钟点分箱。 */
    double tod;

    /* 睡眠分期概率，顺序固定为 [wake, light, deep, rem]，每个值 0–1，四项之和 ≤ 1。
     * 夜间按分期算法输出填写；日间（无分期记录）四项全部填 0。 */
    double stage[4];

    /* 心率，单位 bpm；本样本无心率填 NAN。 */
    double hr;

    /* HRV = ln(RMSSD / ms)，由本样本内的逐拍 RR 间期算出；无 RR 数据填 NAN。
     * 设备给的是 RMSSD 毫秒值时用 energy_hrv_log_from_ms() 换算（45 ms → 3.81）。 */
    double hrv_log;

    /* 本样本时长内的步数。 */
    double steps;

    /* 本样本时长内的活动热量，单位 kcal。 */
    double kcal;

    /* 压力评分，0–1 或 0–100（整条数据流使用同一量纲；库按出现过的最大值自动判断）。 */
    double stress;

    /* 心率观测质量，0–1。1 = 正常；0 = 该点心率不参与估计（伪迹 / 未佩戴）。 */
    double q_hr;

    /* HRV 观测质量，0–1，含义同 q_hr。 */
    double q_hrv;

    /* 血氧饱和度 SpO₂，单位 %，本样本时长内的均值；无血氧填 NAN。
     * 只有夜间血氧被精力层使用；能量层不使用该字段。 */
    double spo2;

    /* 本样本时长内的血氧最低值，单位 %；无则填 NAN。 */
    double spo2_min;
} energy_sample_t;

/* RMSSD(ms) → ln(RMSSD/ms)。<=0 或非有限返回 NAN（与“缺失观测”语义一致）。
 * C 端只有 hrv_log 一个字段，没有 hrv_ms；这是唯一需要的换算。 */
double energy_hrv_log_from_ms(double rmssd_ms);

/*
 * 把样本初始化成默认状态：hr / hrv_log / spo2 / spo2_min = NAN，q_hr / q_hrv = 1，其余 0。
 * 不要用 memset(0) 代替：q=0 会让该样本不参与估计，hr=0 会被当成 0 bpm。
 */
void energy_sample_init(energy_sample_t *s);

/* 输出点标志位（energy_point_t.flags） */
#define ENERGY_FLAG_GAP 1u        /* 本点之前有超过 gap.max_hold_h 的数据缺口（前端提示“数据不足”） */
#define ENERGY_FLAG_GAP_SLEEP 2u  /* 缺口的一部分落在个人睡眠窗内，按软睡眠推进 */
#define ENERGY_FLAG_LEARNING 4u   /* 夜间个体基线仍在冷启动期（前端提示“正在了解你的常态”） */

/* 一个 5 分钟输出点 */
typedef struct {
    double t_h, tod;
    double s, e;      /* 状态 [0,1] */
    double e_sd;      /* E 的后验标准差（开环时为 0） */
    double energy;    /* 0–100，已含展示层节律；批量模式下已做 30 分钟平滑 */
    double lower, upper; /* 95% 置信带（0–100） */
    double base;      /* 100×E（不含节律/平滑） */
    double load, activity, sleep_prob;
    unsigned flags;   /* ENERGY_FLAG_* */
    /* 观测残差（观测值 − 个体化观测模型在校正前的预测），单位 bpm / ln(RMSSD)。
     * 该点没有观测、质量为 0、或个体基线尚未建立时为 NAN。能量模型已经解释了 E、分期、
     * 活动量和节律项，剩下的是“未被能量解释的快变化”——这是精力层（vitality.h）的输入。 */
    double hr_resid, hrv_resid;
} energy_point_t;

/* 每日汇总（按当地日历日切分：以首个样本的钟点 tod 对齐午夜） */
typedef struct {
    int day;           /* 0,1,2... */
    char date[11];     /* "MM-DD"；未配置起始日期时为空 */
    double morning;    /* 当日 tod<9 的最后一点 */
    double evening;    /* 当日 tod>20 的最后一点 */
    double mean, min, max;
    double peak_tod;   /* 06:00–16:00 内的最大值钟点 */
    double dip_tod;    /* 13:00–17:00 内的最小值钟点 */
    double sleep_debt; /* 当日 max S（0–1） */
    int n;             /* 当日样本数 */
} energy_daily_t;

typedef struct {
    energy_point_t *points;
    size_t n_points;
    energy_daily_t *daily;
    size_t n_daily;
    double stats_min, stats_max, stats_mean, stats_band_mean;
    energy_baseline_t baseline; /* 本次计算结束时学到的个体基线（写回画像，下次传 cfg.baseline.init） */
} energy_output_t;

void energy_output_free(energy_output_t *out);

/* 负荷标定结果：用户画像里应固定保存这些量，保证跨窗口可比、流式与批量一致。 */
typedef struct {
    double hr_rest;     /* 静息心率估计 */
    double hr_max;      /* 最大心率估计 */
    double l_ref;       /* 归一化参考值（本次窗口 P99.5 或配置值） */
    double hr_base[24]; /* 按小时的心率基线（288 点钟点基线按小时平均），可直接填 cfg.load.hr_base */
} energy_load_stats_t;

/* 只做负荷标定，不返回序列；用于诊断或写入用户画像。 */
int energy_load_calibrate(const energy_sample_t *samples, size_t n,
                          const energy_load_cfg_t *cfg, energy_load_stats_t *out);

/* ------------------------------------------------------------------ */
/* 3. 主函数：批量计算                                                  */
/* ------------------------------------------------------------------ */
/*
 * 输入 n 个 5 分钟样本，输出逐点曲线 + 每日汇总。
 * - 负荷 L 与活动量 A 由库内部计算；
 * - out 必须是未使用或已 energy_output_free 过的结构（函数内部会先清零，
 *   复用一个仍持有指针的 out 会泄漏）；
 * - 任何返回码下都可以安全调用 energy_output_free。
 * 返回 ENERGY_OK 或 ENERGY_EINVAL / ENERGY_EDOMAIN / ENERGY_ENOMEM。
 */
int energy_compute(const energy_sample_t *samples, size_t n,
                   const energy_config_t *cfg, energy_output_t *out);

/* ------------------------------------------------------------------ */
/* 4. 流式 API（在线每 5 分钟推一次）                                    */
/* ------------------------------------------------------------------ */
#define ENERGY_LREF_BINS 64 /* 流式在线 l_ref 直方图的分箱数（raw ∈ [0,8)，步长 0.125） */

/* 字段视为内部实现细节；需要的量请用下面的访问函数，不要直接读写。 */
typedef struct {
    energy_config_t cfg;
    double s, e;          /* 当前状态 */
    double p00, p01, p11; /* 2×2 协方差 */
    double t_h;           /* 上一次样本时间 */
    double tod;           /* 上一次样本钟点 */
    double t_last_data;   /* 上一个“有数据”样本的时间（缺口判定） */
    double prev_stage[4];
    double prev_load, prev_stress;
    double stress_div;    /* 流式自动量纲：1 或 100（一旦见到 >1.5 即切换） */
    int started;
    int config_error; /* 1 = energy_ctx_init 收到非法配置，后续 step 返回 ENERGY_EINVAL */

    /* 个体基线 */
    double ref_s, ref_e;      /* 开环参考轨迹（同样输入、从不被观测校正） */
    energy_baseline_t bl;     /* 已提交的基线 */
    double acc_sum[4];        /* 当前段的累加器（下标 ENERGY_BL_*） */
    int acc_n[4];
    int run_night, run_day;   /* consecutive night / day samples (episode switch hysteresis) */
    int night_len, day_len;   /* 当前夜/日段的样本数 */
    double night_on, night_off; /* 当前夜段的首/末样本钟点 */
    double night_stage[4];    /* 当前夜段的分期累加 */
    int gate_rej[4];          /* 各通道连续被门控拒绝次数 */
    int anomaly_run[4];       /* 各通道连续同向异常段计数（带符号；|n|>=3 视为常态变化） */
    unsigned long n_obs, n_gated; /* 遥测：参与门控判断的观测数 / 被拒绝数 */
    double resid[2];          /* 本步 HR/HRV 残差（供输出点与精力层使用），无则 NAN */
    double cur_activity;      /* 本步活动量 A（观测方程与日间基线学习使用） */

    /* 流式在线负荷标定（画像缺项时使用） */
    double lc_deep, lc_sleep, lc_all; /* 深睡心率中位数 / 睡眠心率 P20 / 全部心率 P5（在线分位跟踪） */
    int lc_n_deep, lc_n_sleep, lc_n_all, lc_n_hr;
    double lc_peak, lc_last_hr;       /* 衰减峰值心率（相邻两点取小，抗尖峰） */
    double lc_base[24];               /* 按小时的心率 P45（在线分位跟踪，近似批量的去运动中位数） */
    int lc_n_base[24];
    double lc_hist[ENERGY_LREF_BINS]; /* raw 负荷直方图（衰减），用于在线 P99.5 */
    double lc_hist_w;
} energy_ctx_t;

/* 初始化（可在任意时刻重置，等价于“重新冷启动”；cfg->baseline.init 里的基线会被沿用）。 */
void energy_ctx_init(energy_ctx_t *ctx, const energy_config_t *cfg);

/* 取出当前已提交的个体基线（建议每天早上睡眠段结束后落盘到用户画像）。 */
void energy_ctx_get_baseline(const energy_ctx_t *ctx, energy_baseline_t *out);

/* 覆盖个体基线（例如换设备后从画像恢复）；不影响状态 S/E。非法值返回 ENERGY_EINVAL。 */
int energy_ctx_set_baseline(energy_ctx_t *ctx, const energy_baseline_t *bl);

/* 遥测：参与门控判断的观测数与被拒绝数（拒绝率持续 >30% 说明观测模型与该用户不匹配）。 */
void energy_ctx_gate_stats(const energy_ctx_t *ctx, unsigned long *n_obs, unsigned long *n_gated);

/*
 * 推进一步。
 * - 第一次调用：只做观测更新，dt=0；
 * - 之后：用**上一次**有数据样本的分期/负荷/压力按 dt 做解析推进（最多沿用 gap.max_hold_h，
 *   超出部分按缺口策略），再用当前样本观测校正（相对个体基线）。
 * - “无数据”样本（HR/HRV 缺失或 q=0、无分期、无步数）只推进时间，
 *   不覆盖缺口前的输入，因此“整夜空行”与“整夜没有行”结果一致。
 * 返回 ENERGY_OK 或错误码；out 可为 NULL（只更新内部状态）。
 *
 * 注意：流式输出不做展示层的 30 分钟居中平滑（那是离线后处理），
 * 因此 `energy` 与批量结果可能相差 ±1 分。
 */
int energy_ctx_step(energy_ctx_t *ctx, const energy_sample_t *sample, energy_point_t *out);

/* 把当前状态写进 out（不推进时间），用于查询/落盘。 */
void energy_ctx_peek(const energy_ctx_t *ctx, double tod, energy_point_t *out);

/* ------------------------------------------------------------------ */
/* 5. 子函数（可单独调用，便于测试/自检）                                */
/* ------------------------------------------------------------------ */
double energy_phi(const energy_params_t *p, double s);
double energy_wake_equilibrium(const energy_params_t *p, double s, double load, double stress);
double energy_sleep_equilibrium(const energy_params_t *p, double s, double psi);

/* 单步解析推进（不做观测校正），h 单位小时。输出裁剪到 [0,1]。 */
void energy_advance(const energy_params_t *p,
                    double s, double e, double m, double lam, double psi,
                    double load, double stress, double h,
                    double *s_out, double *e_out);

/* 负荷原始量（单点，未归一化）。L = raw / l_ref，再把结果裁剪到 [0,1]。
 * dt_h 为该点覆盖时长；hr 缺失传 NAN。 */
double energy_load_raw(const energy_load_cfg_t *c, double hr, double hr_rest,
                       double hr_max, double hr_base, double steps, double dt_h);
/* 活动量 A ∈ [0,1]（只由步数/热量构成，完全外生）。 */
double energy_activity_step(double steps, double kcal, double steps_ref);

/* 展示层：节律形状（已归一化到 [-1,1]）与 0–100 能量。 */
double energy_circadian(const energy_display_cfg_t *c, double tod);
double energy_display(const energy_display_cfg_t *c, double e, double tod);

#ifdef __cplusplus
}
#endif

#endif /* ENERGY_STATE_H */
