/*
 * energy.h —— 个体能量状态（ODE + EKF）的 C 实现，v1.0
 *
 * 这是 Python 版 `energy/`（energy-state v1.1）的 C 移植，用于嵌入式 / 手机端 / 服务端
 * 的在线计算。模型与参数默认值完全一致，数值对照见 README 的“与 Python 的一致性”。
 *
 * 设计要点
 * --------
 * 1. 状态是显式的 `energy_ctx_t`，不是文件级 static：多用户、多线程、回放、测试
 *    都不会互相污染，也不需要为了换个参数而重启进程。
 * 2. 每次推进用方程的**解析解**（指数更新），不是 `E += delta`：
 *        S(t+h) = S∞ + (S-S∞)·e^(−b·h)
 *        E(t+h) = E∞ + (E-E∞)·e^(−(α+β)·h)
 *    同一个 exp 调用，但精确、恒在 [0,1]、不累积误差，且天然支持可变步长（缺数据）。
 * 3. 两种用法：
 *      - 批量：`energy_compute(samples, n, cfg, &out)`   一次算一整段，适合离线/回放；
 *      - 流式：`energy_ctx_init` + `energy_ctx_step`      每 5 分钟推一次，适合在线。
 *    两者共用同一套子函数，结果逐点一致（流式不做展示层平滑）。
 *
 * 保证与限制
 * ----------
 * 保证：
 *   - 状态 S/E 恒在 [0,1]（解析解 + 两处裁剪），能量恒在 [0,100]；
 *   - 无全局状态：不同 energy_ctx_t 互不影响，可多用户/多线程（每线程一个 ctx）；
 *   - 无随机性；同一输入逐位可复现；
 *   - 所有 double→int 转换都有范围保护（NaN/Inf 不会触发 UB）。
 * 限制（超出返回 ENERGY_EINVAL）：
 *   - |t_h| ≤ 1e6 小时；样本数 ≤ 1e7；时间必须严格递增；
 *   - 缺失值用 NAN（观测类字段）；分期/质量/压力等非观测字段的 NaN 会被清洗；
 *   - 单个缺口最多 4096 个子步（约 1.4 天步长以上按上限收敛，结果仍在 [0,1]）。
 * 资源（Apple clang -O2, 64-bit）：
 *   - energy_ctx_t 712 B；批量接口堆内存约 256 B/样本；
 *   - 最大栈帧 8.3 KB（负荷标定）；流式单步栈帧 432 B。
 *
 * 快速上手
 * --------
 *     energy_config_t cfg;  energy_config_default(&cfg);
 *     energy_ctx_t ctx;     energy_ctx_init(&ctx, &cfg);
 *     for (each 5-min sample) {
 *         energy_point_t p;
 *         energy_ctx_step(&ctx, &sample, &p);
 *         printf("%.1f\n", p.energy);        // 0-100
 *     }
 *
 * 批量：
 *     energy_output_t out = {0};
 *     energy_compute(samples, n, &cfg, &out);
 *     printf("morning %.1f\n", out.daily[0].morning);
 *     energy_output_free(&out);
 */
#ifndef ENERGY_STATE_H
#define ENERGY_STATE_H

#include <stddef.h>

#ifdef __cplusplus
extern "C" {
#endif

#define ENERGY_VERSION "1.0.1"

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
/* 状态方程参数，单位 1/h（无量纲的已注明）。默认值 = Python 端先验均值。 */
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
    double steps_ref;   /* 120 步 / 5 min ≈ 快走 */
    double l_ref;       /* >0 固定参考值（用户画像）；<=0 用本次 P99.5 */
    double hr_rest;     /* <=0 自动：深睡期心率中位数 */
    double hr_max;      /* <=0 自动：历史 P99.5，限幅 [150,210] */
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
} energy_ekf_cfg_t;

/* 观测模型：夜间 [1,E,p_deep,p_rem,p_wake]，日间 [1,E,A,cos1,sin1,cos2,sin2] */
typedef struct {
    double hr_night[5];
    double hr_day[7];
    double hrv_night[5];
    double hrv_day[7];
    double sigma_hr_night, sigma_hr_day;
    double sigma_hrv_night, sigma_hrv_day;
} energy_obs_t;

typedef struct {
    energy_params_t params;
    energy_load_cfg_t load;
    energy_display_cfg_t display;
    energy_ekf_cfg_t ekf;
    energy_obs_t obs;
    double x0_s, x0_e;  /* 状态初值 (0.45, 0.55) */
    int use_ekf;        /* 1: HR/HRV 参与校正；0: 纯 ODE 开环 */
    int day_update;     /* 1: 日间观测也更新状态（与 Python 默认一致） */
    int start_year, start_month, start_day; /* >0 时用于生成 MM-DD 日期标签 */
} energy_config_t;

void energy_params_default(energy_params_t *p);
void energy_load_default(energy_load_cfg_t *c);
void energy_display_default(energy_display_cfg_t *c);
void energy_ekf_default(energy_ekf_cfg_t *c);
void energy_obs_default(energy_obs_t *o);
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
/* 一个 5 分钟样本。缺失值用 NAN；`load`/`activity` < 0 表示由库内部计算。 */
typedef struct {
    double t_h;      /* 距起点小时数，必须单调递增 */
    double tod;      /* 当地钟点 0–24；<0/NaN 时按 t_h % 24 推断 */
    double stage[4]; /* wake/light/deep/rem 概率，和 ≤1。
                      * 约定：全 0 = 无分期记录（日间）；和 >0.5 = 有记录（夜间观测方程）。
                      * 白天请传全 0，不要传 wake=1——否则会走夜间观测方程。 */
    double hr;       /* bpm，缺失 NAN */
    double hrv_log;  /* ln(RMSSD/ms)，缺失 NAN */
    double steps;    /* 每 5 分钟步数，可选（算负荷/活动量用） */
    double kcal;     /* 每 5 分钟热量，可选 */
    double stress;   /* 0–1；若传 0–100 会被 /100 */
    double q_hr;     /* 观测质量 0–1；0 = 该点不参与更新；NaN 视为 0。不知道就传 1 */
    double q_hrv;    /* 同上 */
    double load;     /* 归一化负荷 L∈[0,1]，进状态方程；<0 = 由库计算 */
    double activity; /* 活动量 A∈[0,1]，只进日间观测方程；<0 = 由库计算 */
} energy_sample_t;

/*
 * 把样本初始化成“全缺失、质量正常、负荷/活动量交给库算”的默认状态：
 *   hr/hrv_log = NAN, q_hr/q_hrv = 1, load/activity = -1, 其余 0。
 * 直接 memset(0) 是**错的**：q=0 会让该样本被当成“无效观测”，hr=0 会被当成 0 bpm。
 */
void energy_sample_init(energy_sample_t *s);

/* 一个 5 分钟输出点 */
typedef struct {
    double t_h, tod;
    double s, e;      /* 状态 [0,1] */
    double e_sd;      /* E 的后验标准差（开环时为 0） */
    double energy;    /* 0–100，已含展示层节律；批量模式下已做 30 分钟平滑 */
    double lower, upper; /* 95% 置信带（0–100） */
    double base;      /* 100×E（不含节律/平滑） */
    double load, activity, sleep_prob;
} energy_point_t;

/* 每日汇总（自然日 = 从起点开始的 24h 块） */
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
} energy_output_t;

void energy_output_free(energy_output_t *out);

/* 负荷标定结果：用户画像里应固定保存这三个量，保证跨窗口可比。 */
typedef struct {
    double hr_rest; /* 静息心率估计 */
    double hr_max;  /* 最大心率估计 */
    double l_ref;   /* 归一化参考值（本次窗口 P99.5 或配置值） */
} energy_load_stats_t;

/* 只做负荷标定，不返回序列；用于诊断或写入用户画像。 */
int energy_load_calibrate(const energy_sample_t *samples, size_t n,
                          const energy_load_cfg_t *cfg, energy_load_stats_t *out);

/* ------------------------------------------------------------------ */
/* 3. 主函数：批量计算                                                  */
/* ------------------------------------------------------------------ */
/*
 * 输入 n 个 5 分钟样本，输出逐点曲线 + 每日汇总。
 * - 内部会补算 load/activity（除非样本里给了 >=0 的值）；
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
typedef struct {
    energy_config_t cfg;
    double s, e;          /* 当前状态 */
    double p00, p01, p11; /* 2×2 协方差 */
    double t_h;           /* 上一次样本时间 */
    double prev_stage[4];
    double prev_load, prev_stress;
    int started;
    int config_error; /* 1 = energy_ctx_init 收到非法配置，后续 step 返回 ENERGY_EINVAL */
} energy_ctx_t;

/* 初始化（可在任意时刻重置，等价于“重新冷启动”）。 */
void energy_ctx_init(energy_ctx_t *ctx, const energy_config_t *cfg);

/*
 * 推进一步。
 * - 第一次调用：只做观测更新，dt=0；
 * - 之后：用**上一次**样本的分期/负荷/压力按 dt 做解析推进，再用当前样本观测校正。
 * 返回 ENERGY_OK 或错误码；out 可为 NULL（只更新内部状态）。
 *
 * 注意：流式输出不做展示层的 30 分钟居中平滑（那是离线后处理），
 * 因此 `energy` 与批量结果可能相差 ±1 分；状态 s/e 完全一致。
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
