# energy-state-c —— 个体能量状态（ODE + EKF）的 C 实现

Python 参考实现（ring_personal_agent 仓库的 `energy/`，energy-state v1.1）的 C 移植：一个头文件、一个主函数、
无第三方依赖，用于手机端 / 嵌入式 / 服务端的在线能量计算。

```
输入（5 分钟面板） → 负荷 L / 活动量 A → 解析 ODE 推进 → EKF 用 HR/HRV 校正 → 0–100 曲线 + 置信带 + 每日汇总
```

v1.1.1（修复批量接口忽略调用方 `load`/`activity` 的 bug；v1.1.0 新增 `energy_hrv_log_from_ms()`；v1.0.1 修复流式自动活动量未写回观测方程）。与 Python 的数值一致性：**逐点最大偏差 1e-13**（21 天、6048 点，EKF 与开环两条路径都验证过）。

---

## 0. 能量值该不该做成 static，每 5 分钟 `+=` 一下？

**方向对，做法要改。** “能量是一个随每 5 分钟更新的状态”这个判断是对的，但有三处必须纠正：

| 设想 | 问题 | 做法 |
|---|---|---|
| 只存一个 `energy` 变量 | 状态不止能量：`S`（睡眠压力）决定充电效率 `φ(S)` 和睡眠债消耗 `β_S·S`；`P`（协方差）决定置信带。只存 E 算不回去 | `energy_ctx_t { s, e, p00, p01, p11 }` |
| `E += delta` | 这是欧拉法。模型有解析解 `E(t+h)=E∞+(E−E∞)e^{−(α+β)h}`，同一量级开销，但**精确、恒在 [0,1]、不累积误差**，且天然支持可变步长（断点 6 小时也能一步算对） | 闭式指数更新（`advance_raw`） |
| 文件级 `static` | 多用户/多线程/回放/测试互相污染；无法回滚重算；改参数要重启 | 状态在调用方持有的 ctx 里，可 `energy_ctx_init` 随时重置 |
| 只做加法 | 观测校正不是加常数：EKF 的修正量依赖当前协方差、观测质量、新息大小，还会被门控/稳健化截断 | `ctx_update()` 用 2×2 Kalman + Student-t + 4σ 门控 |

代价只有一个：调用方要多传一个 `energy_ctx_t*`。换来的是可重入、可测试、可回放、可多用户。
如果确实需要全局单例，在应用层包一层即可（本库不提供，避免隐式共享）。

流式写法：

```c
energy_ctx_t ctx;
energy_ctx_init(&ctx, &cfg);
for (;;) {
    energy_sample_t s = read_next_5min_sample();
    energy_point_t p;
    energy_ctx_step(&ctx, &s, &p);      /* 解析推进 + EKF 校正 */
    ui_set_energy(p.energy);            /* 0–100 */
}
```

---

## 1. 目录

```
include/energy.h        公共 API（唯一需要 include 的文件）
src/energy.c            实现：状态方程 / 负荷 / 观测 / EKF / 展示 / 每日汇总
src/main.c              命令行：CSV 进 → CSV 出 + 每日汇总
examples/stream_demo.c  流式用法示例（开环生成真值 → 合成观测 → EKF 跟踪）
tests/test_energy.c     单元测试（无第三方依赖）
tests/bench_energy.c    性能基准
scripts/gen_panel.py    用 Python 生成对照面板与期望输出
scripts/compare.py      比较 C 与 Python
scripts/check_parity.py 一键：生成 → 跑 C → 比较
```

## 2. 构建与运行

```bash
make                 # 生成 build/energy_cli、build/test_energy、build/stream_demo
make test            # 单元测试（61800 个断言）
make sanitize        # AddressSanitizer + UBSan
make coverage        # gcov 行覆盖率
make bench           # 纯计算性能
make example         # 流式示例
make check           # 单测 + 与 Python 的对照（需要 Python 环境）

./build/energy_cli data/sample_panel.csv --out out.csv   # 自带 1 天样例面板
./build/energy_cli data/sample_panel.csv --load-stats    # 只打印 hr_rest/hr_max/l_ref
```

要求：C99 编译器 + libm。已在 Apple clang 21 / `-Wall -Wextra -Wpedantic -Wconversion` 下零告警。

## 3. API

### 3.1 批量（主函数）

```c
energy_config_t cfg;
energy_config_default(&cfg);          /* 群体先验；生产环境覆盖 params/obs/x0 */

energy_output_t out = {0};
int rc = energy_compute(samples, n, &cfg, &out);
if (rc != ENERGY_OK) { /* ENERGY_EINVAL / ENERGY_EDOMAIN / ENERGY_ENOMEM */ }

for (size_t i = 0; i < out.n_points; i++)
    printf("%.1f\n", out.points[i].energy);      /* 0–100，已含展示层节律与平滑 */

for (size_t i = 0; i < out.n_daily; i++)
    printf("%s morning=%.1f debt=%.2f\n",
           out.daily[i].date, out.daily[i].morning, out.daily[i].sleep_debt);

energy_output_free(&out);
```

`energy_compute` 内部流程（主函数只做编排，计算都在子函数里）：

```
energy_compute
├─ fill_linear / sanitize_sample      输入清洗（NaN 补缺、分期归一）
├─ estimate_hr_rest / estimate_hr_max / estimate_hr_base   负荷标定
│   └─ energy_load_raw                每点原始负荷 → P99.5 归一化成 L
│   └─ energy_activity_step           活动量 A（只由步数/热量构成）
├─ energy_ctx_init                    状态初值 (0.45, 0.55)
├─ for each sample: energy_ctx_step
│   ├─ stage_mix                      m / λ / ψ
│   ├─ ctx_predict → advance_raw      解析推进 + 连续雅可比 + expm(F·h)
│   ├─ ctx_update  → obs_predict      夜间/日间观测方程 + 稳健 EKF
│   └─ ctx_emit    → energy_display   节律调制 + 置信带
├─ smooth_display                     30 分钟居中平滑（只影响展示）
└─ 每日汇总                           morning/evening/peak/dip/sleep_debt
```

### 3.2 流式

```c
energy_ctx_t ctx;
energy_ctx_init(&ctx, &cfg);                 /* 可在任意时刻重置 = 冷启动 */
energy_point_t p;
energy_ctx_step(&ctx, &sample, &p);          /* 每个 5 分钟样本调用一次 */
energy_ctx_peek(&ctx, tod, &p);              /* 只查当前状态，不推进时间 */
```

语义：第一次调用只做观测更新；之后每次用**上一个**样本的分期/负荷/压力按 `dt` 推进，
再用当前样本的 HR/HRV 校正。`dt` 可以不是 5 分钟（缺数据、断点都行）。
流式不做展示层的 30 分钟居中平滑（那是离线后处理），所以 `energy` 可能与批量差 ±1 分；
状态 `s/e` 与批量逐点一致（单测里断言 1e-12）。

### 3.3 负荷标定（写入用户画像）

```c
energy_load_stats_t st;
energy_load_calibrate(samples, n, &cfg.load, &st);
/* st.hr_rest, st.hr_max, st.l_ref → 存进用户画像，之后固定传入 cfg.load */
```

不固定 `l_ref` 的话，负荷是**窗口内相对量**：整段活动恒定会被归一成 ~0.8，跨周不可比。

### 3.4 配置校验与版本

```c
char err[128];
if (energy_config_validate(&cfg, err, sizeof(err)) != ENERGY_OK)
    fprintf(stderr, "配置非法: %s\n", err);

printf("energy-state-c %s\n", energy_version());   /* "1.0.0" */
```

`energy_compute` 会自动校验（非法配置返回 `ENERGY_EINVAL`）；
流式用法在 `energy_ctx_init` 里校验，非法配置会把 `ctx.config_error` 置 1，之后每次
`energy_ctx_step` 返回 `ENERGY_EINVAL`——不会静默算出垃圾。

### 3.5 可单独调用的子函数

| 函数 | 用途 |
|---|---|
| `energy_phi(p, s)` | 恢复驱动 φ(S)，含下限与上下限 |
| `energy_advance(...)` | 单步解析推进（无观测），测试/仿真用 |
| `energy_load_raw(...)` / `energy_activity_step(...)` | 单点负荷/活动量 |
| `energy_wake_equilibrium` / `energy_sleep_equilibrium` | 平衡点（量纲检查） |
| `energy_circadian` / `energy_display` | 展示层节律与 0–100 映射 |

## 4. 输入

```c
energy_sample_t s;
energy_sample_init(&s);   /* 推荐：hr/hrv=NAN, q=1, load/activity=-1 */
s.t_h = t; s.tod = tod;
s.stage[...] = ...; s.hr = ...; s.hrv_log = ...;
s.steps = ...; s.stress = ...;
```

> **不要 `memset(s, 0, ...)`**：`q=0` 会被当成“无效观测”而跳过滤波，`hr=0` 会被当成 0 bpm。
> 用 `energy_sample_init()`。

| 字段 | 必须？ | 含义 | 不传/缺省时 |
|---|---|---|---|
| `t_h` | **必须** | 距起点小时，严格递增 | — |
| `tod` | 可选 | 当地钟点 0–24 | `<0`/NaN → 由 `t_h%24` 推断 |
| `stage[4]` | 建议 | wake/light/deep/rem 概率 | 全 0 = 无分期记录（**日间**） |
| `hr` | 建议 | 心率 bpm | NaN = 该点无观测 |
| `hrv_log` | 可选 | **ln(RMSSD/ms)**，每样本一个标量 | NaN = 无观测 |
| `steps`/`kcal` | 可选 | 每 5 分钟步数/热量 | 0 |
| `stress` | 可选 | 0–1（传 0–100 自动 /100） | 0 |
| `q_hr`/`q_hrv` | 可选 | 观测质量 0–1 | **1**（`energy_sample_init` 设好） |
| `load` | 可选 | 归一化负荷 L∈[0,1]，进**状态方程** | `<0` → 由库计算 |
| `activity` | 可选 | 活动量 A∈[0,1]，只进**日间观测方程** | `<0` → 由库计算 |

**四个“可选”字段到底管什么：**

- **`q_hr` / `q_hrv`（观测质量）**：只影响 EKF。`R = σ² / max(q, 0.05)`，q 越小越不信任这个点；
  `q = 0` 直接跳过该通道（只预测、不更新）。设备给出质量位就透传（5 秒质量在 5 分钟窗口内取均值）；
  没有质量信息就保持 1。NaN 视为 0。
- **`load`（负荷 L）**：进状态方程的消耗项 `β = β₀ + β_L·L + β_P·stress + (1−m)·β_S·S`，
  是白天掉电快慢的主因（实测 L=0 与 L=0.6 的日均能量差 14 分）。不传时由库用
  `心率 TRIMP + 步数 + 心率超基线` 算，再按 P99.5 归一化；流式接口没有全局标定，
  会用 `cfg.load` 里的 `hr_rest/hr_max/l_ref`（默认 60/185/0.25），
  **生产环境应把 `l_ref` 固定成用户画像里的值**，否则不同窗口不可比。
- **`activity`（活动量 A）**：只进**日间观测方程**（解释掉运动导致的心率升高），不直接改能量。
  不传时按 `0.8·steps/120 + 0.2·kcal/6` 算。它与 `load` 分开是刻意的：
  `load` 由心率算出，再拿它回归心率就是内生性（会把 `a_E` 从 −12 拉到 −38）。
- 二者都传时，库不再计算；只传一个时，另一个仍会自动算。

模型内部派生的三个量：`m=Σp_sleep`、`λ=Σp_s·κ_s`、`ψ=Σp_s·w_s`。

> **stage 的约定**：是否有分期记录决定用夜间还是日间观测方程（`Σstage > 0.5` 即夜间）。
> 白天请传全 0，**不要传 `wake=1`**——那会被当成“夜里清醒”，走夜间系数。

### 4.1 调用者输入清单（最优 / 最小）

**一句话**：只有 `t_h` 是绝对必须；`stage` 和 `hr` 决定“这是不是一个能量模型”；
其余都有默认值，但省略会损失信息。配置层（`cfg.obs`、`l_ref`）的影响通常比样本字段更大。

**样本字段**

| 优先级 | 字段 | 作用 | 省略后会发生什么 |
|---|---|---|---|
| **必须** | `t_h` | 时间轴 | 无法计算 |
| **必须**（否则模型退化） | `stage[4]` | 决定睡眠概率 m、消散 λ、恢复权重 ψ | m=0，全程清醒：同一天 07:00 的 E 从 **0.711** 掉到 **0.365**（夜里不充电） |
| **必须**（否则只剩先验） | `hr` | 主观测通道 | EKF 退化开环，置信带宽=0，曲线只由模型先验决定 |
| 强烈建议 | `hrv_log` | 第二观测通道（夜间最重要） | 只剩 HR 校正（HRV 本身信息弱，但少一个约束） |
| 强烈建议 | `steps_5m`（或 `load`） | 负荷 L 的主要来源 | L 只由 HR 的 TRIMP/ΔHR 算；HR 也缺时 L=0，曲线几乎不下降 |
| 强烈建议 | `stress` | β_P 消耗项（满量程约 3×β0） | 少一项消耗，能量略偏高 |
| 强烈建议 | `tod`（或 `start`） | 昼夜节律、日间观测窗口、基线分箱 | 按 `t_h%24` 推断；起点不是 00:00 时整体错位 |
| 强烈建议 | `q_hr`/`q_hrv` | 伪迹/摘表时关掉该点 | 默认 1（全部相信）：持续 +6 bpm 偏差可把 E 拉偏 **0.36** |
| 可选 | `kcal` | 只影响活动量 A | A 只由步数算 |
| 可选 | `load` / `activity` | 上游已算好时直接传（批量/流式都尊重，>=0 即采用） | 库内部按固定规则计算 |

**配置层（比样本字段更影响质量）**

| 配置 | 最优做法 | 不配时 |
|---|---|---|
| `cfg.obs` | 传离线拟合的 `fit.json` | 群体先验（夜 58 / 日 80 bpm 基线），设备基线不同会拖偏 |
| `cfg.params` | 同上（个体化） | 群体先验 |
| `cfg.load.l_ref/hr_rest/hr_max` | 用户画像里的固定值 | 批量按本次窗口估；流式用 60/185/0.25 |
| `cfg.x0` | 接续计算时传上次末状态 | (0.45, 0.55) 冷启动 |
| `cfg.display` | 产品统一配置 | ρ=0.10、10:00 峰、14:30 谷、30 分钟平滑 |

**三档调用**

```c
/* ① 最小可用：只有睡眠分期 + 心率（能出曲线，但没有 HRV/活动/压力信息） */
energy_sample_t s; energy_sample_init(&s);
s.t_h = t; s.tod = tod;
s.stage[ENERGY_LIGHT] = p_light;  /* 夜里；白天全 0 */
s.hr = hr;

/* ② 推荐（生产默认）：加 HRV、步数、压力、质量位，并传离线拟合的观测模型 */
s.hrv_log = hrv;  s.steps = steps;  s.stress = stress;
s.q_hr = q_hr;    s.q_hrv = q_hrv;
cfg.obs = load_fit_json("fit.json");      /* 或自己填 cfg.obs */
cfg.load.l_ref = user_profile.l_ref;      /* 跨窗口可比的关键 */

/* ③ 最全：上游已有同款管线时，直接给归一化后的 load/activity，保证与 Python 曲线逐点一致 */
s.load = load;  s.activity = activity;
```

**不要这样做**

- `memset(&s, 0, sizeof(s))`：`q=0` 会被当成“无效观测”跳过滤波，`hr=0` 会被当成 0 bpm。
- 白天传 `stage[WAKE]=1`：判据是“有没有分期记录”，白天全 0 才走日间观测方程。
- 一半字段传、一半让库算：`load`/`activity` 的归一化路径不同，混用会让曲线不可比。
- 把 `load` 当成“绝对运动量”：它是相对个人历史 P99.5 的归一化值。

### 4.2 HRV 的格式：不是数组，是每 5 分钟一个标量

C 端 `energy_sample_t` 里只有 **`hrv_log`** 一个字段（Python 接口的 `hrv_ms` 是便捷参数，C 没有）：

```c
/* 设备给 RMSSD(ms) → 取自然对数再填 */
s.hrv_log = energy_hrv_log_from_ms(rmssd_ms);   /* 45 ms → 3.8067 */
/* 或自己算 */
s.hrv_log = log(rmssd);                          /* rmssd <= 0 时填 NAN */
```

| 项 | 说明 |
|---|---|
| 类型 | `double`，**每个 5 分钟样本一个值**，不是数组（批量时才是一组样本的数组） |
| 单位 | `ln(RMSSD/ms)`；典型 3.0–4.5（RMSSD 20–90 ms），全范围约 1.6–5.3 |
| 缺失 | `NAN`（不是 0，也不要用 0 代替） |
| 质量 | 用 `q_hrv ∈ [0,1]` 表达“这一点可不可信”；0 = 跳过 |
| 日夜 | **不限定夜间**：模型有夜间和日间两套 HRV 系数，按“有没有分期记录”自动选。设备只在夜间给 HRV 时，白天填 `NAN` 即可 |

> 最常见的集成错误：把毫秒直接填进 `hrv_log`（45 而不是 3.81）。量级差 10 倍，观测会被 4σ 门控全部拒掉——不报错，但 HRV 通道等于失效。
> `energy_hrv_log_from_ms()` 就是为这个准备的；它会拒绝 ≤0 / 非有限值并返回 `NAN`。

日/夜两套 HRV 系数（默认值）：

| 通道 | 特征 | 系数 |
|---|---|---|
| 夜间 | `[1, E, p_deep, p_rem, p_wake]` | `[3.55, 0.45, 0.30, 0.18, −0.35]`，σ=0.22 |
| 日间 | `[1, E, A, cos1, sin1, cos2, sin2]` | `[3.25, 0.45, −0.55, −0.06, −0.10, 0.025, −0.043]`，σ=0.40 |

## 5. 输出

| 结构 | 字段 |
|---|---|
| `energy_point_t` | `t_h, tod, s, e, e_sd, energy, lower, upper, base, load, activity, sleep_prob` |
| `energy_daily_t` | `day, date, morning, evening, mean, min, max, peak_tod, dip_tod, sleep_debt, n` |
| `energy_output_t` | `points/n_points, daily/n_daily, stats_min/max/mean/band_mean` |

`energy` 恒在 0–100；`s/e` 恒在 0–1；`lower/upper` 是 95% 置信带（开环时退化为等于 `energy`）。

## 6. 参数默认值

与 Python 完全一致（`energy_config_default`）。

| 参数 | 值 | 含义 |
|---|---|---|
| `kappa_w` | 0.100 /h | 清醒期 S 积累（由 16h→S=0.8 反解） |
| `kappa_s[4]` | 0.05 / 0.26 / 0.42 / 0.18 | wake/light/deep/rem 的 S 消散 |
| `alpha_sleep` | 0.300 /h | 睡眠充电 |
| `alpha_rest` | 0.020 /h | 清醒静息恢复 |
| `w_stage[4]` | 0.15 / 1.00 / 1.25 / 0.95 | 分期恢复权重（light=1 归一化） |
| `eta / s_ref / s_floor` | 0.60 / 0.50 / 0.10 | φ(S) 形状 |
| `phi_min / phi_max` | 0.30 / 2.00 | φ 上下限 |
| `beta_0 / beta_L / beta_P / beta_S` | 0.010 / 0.150 / 0.030 / 0.035 | 基础/负荷/压力/睡眠债消耗 |
| `x0_s, x0_e` | 0.45, 0.55 | 状态初值 |
| `q_s, q_e` | 0.02², 0.03² | 过程噪声 |
| `p0_s, p0_e` | 0.18², 0.22² | 初始协方差 |
| `gate_sigma / q_floor / robust_nu` | 4.0 / 0.05 / 4.0 | 门控 / 质量下限 / Student-t |
| `rho / peak / dip` | 0.10 / 10:00 / 14:30（深 0.7、宽 1.6h） | 展示层节律 |
| `smooth_minutes / band_z` | 30 / 1.96 | 展示平滑 / 置信带 |
| `b_trimp / w_* / trimp_scale / steps_ref` | 1.92 / 1.0,0.6,0.4 / 10 / 120 | 负荷 |

**依据**：只有 `κ_w` 是由“16h 后 S=0.80”反解出来的（`−ln(0.2)/16=0.1006`）；
其余是“生理合理性 + 工程试错”的群体先验，用合成数据验证可辨识性。11 个参数里只有 5 个有清晰剖面峰，
所以**不要单独解释某个参数**，要解释它们组合出的可观测量（充电 τ≈3h、消耗 τ≈5–10h、8h 睡眠 ΔE≈0.4–0.55）。
详见 Python 参考仓库的《能量计算分析》§6.4（参数依据与证据强度）。

## 7. 与 Python 的一致性

```bash
# 需要 Python3 + pandas + 参考实现（energy_api.py / energy_state）。
# 参考实现不在此仓库时，用环境变量指定位置：
export ENERGY_PY_ROOT=/path/to/python/reference      # 含 energy_api.py
export ENERGY_STATE_SRC=/path/to/reference/energy/src # 含 energy_state/
python scripts/check_parity.py 21
```

21 天（6048 点）、EKF 与开环两条路径：

| 列 | C vs Python 最大偏差 |
|---|---|
| `s` / `e` | 1.2e-15 / 7.8e-16 |
| `energy` | 9.9e-14 |
| `lower` / `upper` | 8.5e-14 / 9.9e-14 |
| `hr_rest` / `hr_max` / `l_ref` | 0 / 0 / 0（相对误差） |

偏差来自 CSV 往返的十进制舍入，不是算法差异。

## 8. 性能

`make bench`（Apple clang -O2，单核）：

| 路径 | 21 天 6048 点 | 每点 |
|---|---|---|
| 批量 `energy_compute`（含负荷标定 + EKF + 展示 + 日汇总） | 27 ms | 4.5 µs |
| 流式 `energy_ctx_step` | 18 ms | 3.0 µs |

对比 Python 版 `compute_energy`：762 ms / 次（约 28 倍）；C 命令行含 CSV I/O 约 40 ms。

## 9. 移植范围

已移植：负荷/活动量、状态方程、观测模型（夜间/日间、固定系数）、EKF（稳健 t、门控、Joseph、
2.5 分钟子步的开环积分）、展示层（节律、平滑、置信带）、每日汇总。

**未移植**（Python 端的上层功能）：

| 功能 | 影响 | 替代做法 |
|---|---|---|
| 异常夜冻结 + 二次回放 | Python 的头条 RMSE 依赖它（0.131 → 0.045） | 在输入层把异常夜/次日的 `q_hr/q_hrv` 置 0，再跑一遍 |
| 睡眠报告（episode / recovery_ratio） | 少一组产品字段 | C 端可用 S/E 轨迹自行分段计算 |
| 冷启动阶段 / 文案 | 纯产品层 | 上层实现 |
| 离线拟合（MAP/剖面似然） | C 端不拟合 | 用 Python 离线拟合，把 `params/obs/x0` 填进 `energy_config_t` |

## 10. 嵌入式注意事项

- 流式路径**不分配内存**（只有 `energy_compute` 用 malloc），可在中断/任务里调用；
- 全部 `double`；若要用 `float`，需要把 `energy.h` 里的字段与实现一起改，并重新跑对照；
- 无全局状态，可多实例；一个 ctx 约 150 字节；
- `energy_circadian` 每次会算 288 点归一化因子（约 600 次浮点运算），流式每 5 分钟一次可忽略；
  高频调用可自己缓存 `energy_display` 的结果；
- 观测模型是**群体先验**（夜间 58 bpm、日间 80 bpm 基线）。设备基线不同必须覆盖 `cfg.obs`，
  否则滤波器会把固定偏移当成能量偏差（Python 端同样问题，实测 +9 bpm 能把 E 拖低 0.2）。

## 11. 生产使用

**输入契约**

| 项 | 约束 | 违反时 |
|---|---|---|
| `t_h` | 有限、严格递增、`|t_h| ≤ 1e6` 小时 | `ENERGY_EINVAL` / `ENERGY_EDOMAIN` |
| 样本数 | `2 ≤ n ≤ 1e7` | `ENERGY_EINVAL` |
| 缺失值 | `hr`/`hrv_log` 用 `NAN`；其它字段的 NaN 会被清洗 | 自动处理 |
| 时间步长 | 推荐 5 分钟；其它步长可用（解析解对任意 dt 稳定） | — |
| 单缺口 | 最多 4096 个子步（约 1.4 天以上按上限收敛） | 自动处理 |

**资源占用**（Apple clang -O2，64-bit）

| 项 | 数值 |
|---|---|
| `energy_ctx_t` | 712 B（含整份 config） |
| 批量堆内存 | 约 256 B/样本（6048 点 ≈ 1.5 MB） |
| 最大栈帧 | 8.3 KB（负荷标定的 288 点分箱）；流式单步 432 B |

嵌入式建议用流式接口（无 malloc、512 B 级栈）；批量接口留给手机/服务端/离线回放。

**并发与重入**

无全局变量、无静态可写状态；一个 `energy_ctx_t` 只被一个线程使用即可。负荷标定是纯函数。

**数值健壮性**

- 所有 `double → int` 转换都有有限性/范围保护（NaN 的钟点不会再触发 UB；已用 UBSan 验证）；
- 极长缺口导致协方差失去意义时，自动回到冷启动方差，不会永久卡死；
- `q` 裁剪到 [0,1]，`smooth_minutes` 上限 1440，避免非法值放大成异常结果。

**构建**

- C99 + libm；GCC/Clang 已验证（`-Wall -Wextra -Wpedantic -Wshadow -Wconversion` 零告警，
  clang static analyzer 干净）；
- MSVC 需要 C11 模式与 `isfinite/_finite` 适配（未验证，不在支持列表）；
- 仓库带 `.gitlab-ci.yml`：`make lib`（-Werror）+ `make test` + `make sanitize`。

## 12. 测试

```
make test        # 61800 个断言
make sanitize    # ASan + UBSan，0 error
make coverage    # src/energy.c 行覆盖率 88.54%
make check       # 单测 + 3 天 C/Python 对照
python scripts/check_parity.py 21
```

覆盖的内容：φ 形状与上下限、解析推进、平衡点、展示节律、负荷/活动量、批处理、流式与批量逐点一致、
6 小时缺口、开环、非法参数、S/E 裁剪、门控（分别关掉稳健 t 验证）、昼夜观测系数差异、
金标准回归（5 个参考点，容差 1e-9）、80 组随机参数 × 72h 有界性。

覆盖范围还包括：配置校验、非法/缺失 tod、t≥1e6h、q>1、极端平滑、非法起始日期、
n=2 最小面板，以及 20 轮随机 NaN 注入的模糊测试（断言输出有界、无 UB）。

**变异测试（手工，9 个变异体）**：7 个被杀死；2 个存活是**等价变异**——`ctx_predict` 与 `ctx_clamp`
两处都做了 [0,1] 裁剪，去掉任意一处行为不变（同时去掉两处则 58 个断言失败，已验证）。
去掉 Joseph 项、夜间误用日间系数、日间漏掉活动量 A 这三个变异都被金标准回归杀死。

关键断言：状态恒在 [0,1]；能量恒在 [0,100]；置信带包住点估计；流式与批量 `s/e` 逐点一致；
80 组随机参数 × 72h 零越界；非法输入返回 `ENERGY_EINVAL/EDOMAIN`。

## 13. 决策日志

- **状态用 ctx，不用 static**：多用户与测试隔离；代价是调用方多传一个指针。拒绝的理由（“调用麻烦”）不足以抵消全局可变状态的污染。
- **解析解，不用欧拉 `E+=Δ`**：同量级开销，换来恒有界、可跳步、无数值漂移。欧拉在 6 小时断点上误差会到 1e-3 量级。
- **开环用 2.5 分钟子步，EKF 用单步**：与 Python `simulate` / `run_ekf` 的既有实现对齐（前者二阶积分，后者单步 + 中点 φ），保证 1e-13 对照；不是随意选择。
- **流式不做 30 分钟居中平滑**：居中窗口需要“未来”数据，在线不可得；平滑只在批量做，并在文档里标明 ±1 分差异。
- **负荷标定单独暴露**：`l_ref` 必须进用户画像，否则窗口间不可比；这是 Python 端踩过的坑（评审 #4）。
- **未移植异常夜冻结**：它需要两遍回放与跨夜基线，属于上层编排；在 C 端用 `q=0` 表达同样的“不信任这段观测”，更简单。
- **拒绝把 `energy_circadian` 的归一化固定成 1.0**：会与 Python 差 ~1e-5；改用 288 点网格取 max，整天覆盖时逐位一致，短面板差异可忽略。
