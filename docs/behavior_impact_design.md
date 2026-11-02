# 行为对晨起精力的影响：设计方案

模块：`include/evobeing_behavior_zhenghuimin.h` / `src/evobeing_behavior_zhenghuimin.c`
测试：`tests/test_behavior.c`（707 项断言）　示例：`examples/behavior_demo.c`

## 1. 要解决的问题

每天早上起床，精力引擎给出一个 0–100 的初始值。用户想知道：过去 28 天里，我做过的哪些事（饮酒、咖啡、小睡、运动、晚睡……）让第二天早上的精力变好或变差了，各占多少。

约束：

- 只跟自己比，不跟别人比。
- 每种行为输出一个百分比，各项之和不必是 100%，因为晨起精力还受睡眠本身、身体状况、天气等没有记录的因素影响。
- 行为有两类来源：系统自动识别（小睡、运动、过载、晚睡），用户手动记录（饮酒、咖啡、冥想、压力事件、晚餐很晚、睡前看屏幕，以及自定义）。记录形式统一为"行为 + 时间点"。
- 这是消费电子产品里的自我观察工具，不是医学结论。文案用"关联""有这个行为的早晨"，不用"导致""诊断"。
- 不改动 `evobeing_energy_zhenghuimin` / `evobeing_vitality_zhenghuimin` 两个主函数，不引入 load / activity 输入。本模块是独立的离线分析，建议在手机端或云端每天起床后跑一次。

## 2. 业内做法与本方案的取舍

| 产品 | 做法 | 借鉴 |
|---|---|---|
| WHOOP Journal | 用户勾选行为，30 天后显示"有该行为时恢复分平均高 / 低 X%"，附次数 | 展示口径：有 / 无对比的百分比，带次数 |
| Oura Tags | 标签与睡眠 / 准备度的个人关联，需要足够次数 | 最少次数门槛，否则提示继续记录 |
| Garmin / Fitbit | 不做逐行为归因，只给因素分解（睡眠、压力、活动） | "未解释部分"的呈现 |
| Apple | 不做 | — |

头部产品的共同短板是只做**简单对比**（有行为的日子 vs 没有的日子），不控制混杂。典型问题：饮酒集中在周五、周六晚上，而周末早晨本来就睡得久、精力高，简单对比会把饮酒的负面效应抵消甚至算成正的。本方案在简单对比之上加了一层**个体内加权岭回归**，把其他行为、周末、前一天的精力、近期趋势一起放进模型，再读行为系数；同时保留简单对比作为参考值输出。

## 3. 数据定义

### 3.1 目标量

`y_d` = 第 d 天起床时刻精力引擎的初始输出，即晨起 WAKE 点的 `score`（`vitality_daily_t.wake_score`，0–100）。一天只有一个值；当天没有主睡眠（`d_init` 为 NAN）的日子不参与分析。

### 3.2 行为记录

每条记录 `(id, t_h)`，`t_h` 是自窗口第 0 天 00:00 起的小时数。

**归属规则**：一条记录归到它之后的第一次起床。起床后到下一次起床之间发生的所有行为，影响下一次起床的 y。晚上 23:00 的饮酒和凌晨 01:00 的饮酒都归到接下来的那个早晨。最后一次起床之后的记录（今天白天发生的）还没有对应的早晨，暂不计入。

**自动行为**由调用方从引擎每日汇总生成，建议规则：

| id | 来源 | 建议规则 |
|---|---|---|
| NAP | `vitality_daily_t.n_nap > 0` 或 `VITALITY_FLAG_NAP` | 小睡结束时刻 |
| EXERCISE | `n_exercise > 0` 或 `VITALITY_FLAG_EXERCISE` | 运动结束时刻 |
| OVERLOAD | `n_stress_samples ≥ 12`（即 STRESS 累计 ≥ 1 小时） | 当天 20:00 |
| LATE_BED | 入睡时刻比本人近 14 天中位数晚 ≥ 1 h | 入睡时刻 |

**手动行为**：ALCOHOL、CAFFEINE、MEDITATION、STRESSOR、LATE_MEAL、SCREEN，自定义 id 16–23。引擎里已有 `evobeing_event_zhenghuimin()` 的用户标记，可以同时记一条到行为日志。

### 3.3 采集频率对可用性的影响

按《健康指标日间夜间采集频率》表：心率日间主要每 5 秒一条、夜间每 5 秒或 1 分钟一条；HRV 每 2–20 分钟一条；睡眠按阶段保存、每段睡眠一份汇总；步数 / 锻炼时长每 15 分钟一个汇总点。这些频率对本模块足够：本模块只用每天一个 `wake_score` 和行为时刻，不直接碰原始采样。需要注意的是"多固件混合，不能视为统一硬件采样配置"——不同固件的 `wake_score` 来自同一个精力引擎，刻度一致；但若某天夜间 HRV 稀疏到引擎没有做晨起初始化，那天就是 NAN，模块会自动跳过。

## 4. 算法

### 4.1 窗口与权重

取最近 `window_days = 28` 天。有效天 = 有 `wake_score` 且有起床时刻的天。有效天不足 7 天时不分析（返回 −2）。

每天权重 `w_d = 0.5^((N−1−d) / half_life)`，`half_life = 14` 天：昨天权重 1，两周前 0.5，四周前 0.25。近期的行为变化更快反映到结果里。

`ȳ` = 加权均值，`s_y` = 加权标准差。

### 4.2 特征列

对每种行为生成 0/1 列 `x_j[d]`（第 d 个早晨之前有没有这个行为）。

**时间点细分**：对有"晚间阈值"的行为（默认咖啡 14:00、运动 20:00、饮酒 21:00），若该行为在阈值之后发生，再生成一列 `x_j^late[d]`。这样可以区分"上午一杯咖啡"和"下午三点之后的咖啡"。细分列只在晚间次数 ≥ `min_with` 且不是全部都在晚间时才加（否则与主列完全重合，没有额外信息）。

**样本门槛**：行为列进入回归要求 28 天内 `n_with ≥ 4` 且 `n_without ≥ 4`。不满足的行为标 `INSUFFICIENT`，不给数字，前端显示"再记录 N 次就能看到"。

**控制列**（默认全开，可关）：

- 周末：周六 / 周日早晨 = 1（即周五 / 周六晚）。
- 前一天精力：`(y_{d−1} − ȳ) / s_y`，缺则 0。吸收精力的自相关（连着几天累的、或者刚度假回来）。
- 趋势：`(d − (N−1)/2) / N`。吸收四周内的缓慢变化（季节、工作强度）。

所有列按加权均值中心化，截距即 `ȳ`。

### 4.3 加权岭回归

最小化

    Σ_d w_d (y_d − ȳ − Σ_j x_jd b_j)² + Σ_j λ_j b_j²
    λ_j = ridge · Var_w(x_j) · Σw / n_valid

`ridge = 2.0` 相当于给每一列加 2 个"效应为 0 的伪日"：在 28 天里收缩约 7%，足以防共线和极端值，又不明显压低效应。用部分主元高斯消元求解（p ≤ 51）。

`b_j` 的含义：控制其他列后，有行为 j 的早晨比没有时高 / 低多少分。

### 4.4 输出口径

**口径 A（主显示）**：`effect_pct_j = 100 · b_j / ȳ`。

- 例：ȳ = 60，b_alcohol = −9 → −15%，读作"饮酒后的早晨，精力比你平时低约 15%"。
- 晚间细分项的效应 = `b_j + b_j^late`，即"晚间发生时的总效应"，方便直接显示"14 点后的咖啡 −9%"。
- 各项之和不必为 0 也不必为 100。
- 同时给 `raw_pct`：未控制混杂的简单加权对比，`100·(mean_with − mean_without)/ȳ`，供调试和与 WHOOP 式口径对照。

**口径 B（内部）**：`share_pct_j = 100 · (R²_full − R²_{−j})`，去掉列 j 后加权 R² 的下降量。所有项之和被封顶在 `100·(R²_full − R²_controls)`，即行为合计解释的份额，也不超过 100。`unexplained_pct = 100 − Σ share`。用于排序辅助和"记录的行为合计解释了 X%，其余来自未记录因素"这句话。晚间细分项的份额是它在主列之上的**增量**，所以可能很小。

### 4.5 不确定性与可信度

自助法：对有效天按天有放回重采样 `n_boot = 200` 次，权重随行，重新拟合，取系数的 5%–95% 分位作 90% 区间。

档位：`ratio = |effect| / 半区间宽`。

- HIGH：ratio ≥ 2 且 n_with ≥ 6
- MID：ratio ≥ 1
- LOW：其他（区间跨 0，方向还不确定）
- INSUFFICIENT：没进回归

`n_boot = 0` 时不算区间，按 n_with ≥ 8 给 MID，否则 LOW。固定种子，同样输入结果完全可复现。

### 4.6 排序

有数字的按 |effect_pct| 降序，INSUFFICIENT 排最后。

## 5. 函数的输入与输出

### 5.1 `evobeing_behavior_impact_zhenghuimin`

```c
int evobeing_behavior_impact_zhenghuimin(const evobeing_behavior_input_t *in,
                                         const evobeing_behavior_config_t *cfg,   /* NULL = 默认 */
                                         evobeing_behavior_output_t *out);
```

返回 `0` 成功；`-1` 参数错误（指针为空、`n_days` 不在 1–56、事件数超 2048）；`-2` 有效天数不足 7 天（`out->n_days_used` 已填，各项全部 INSUFFICIENT，`unexplained_pct = 100`）。

**输入 `evobeing_behavior_input_t`**

| 字段 | 类型 | 含义 |
|---|---|---|
| `n_days` | int | 天数，1–56，第 0 天最早；模块只用最后 `window_days`（28）天 |
| `wake_score` | const double[n_days] | 每天晨起精力 0–100；没有则 NAN |
| `wake_t_h` | const double[n_days] | 每天起床时刻（自第 0 天 00:00 起的小时数）；NULL 按每天 7:00；NAN = 当天没起床记录 |
| `weekday` | const int[n_days] | 0 = 周一 … 6 = 周日；NULL 则不做周末控制 |
| `events` | const evobeing_behavior_event_t[n_events] | `{id, t_h}`，id 见 `EVOBEING_BEHAVIOR_*`，t_h 同上 |
| `n_events` | int | ≤ 2048 |

**配置 `evobeing_behavior_config_t`**（`evobeing_behavior_config_default()` 取默认）

| 字段 | 默认 | 含义 |
|---|---|---|
| `window_days` | 28 | 分析窗口 |
| `half_life_days` | 14 | 近期权重半衰期 |
| `min_with` / `min_without` | 4 / 4 | 给数字的最少有 / 无天数 |
| `ridge` | 2.0 | 收缩强度（伪日数） |
| `n_boot` | 200 | 自助法次数，0 = 不算区间，上限 500 |
| `seed` | 12345 | 随机种子 |
| `control_weekend` / `control_prev` / `control_trend` | 1 / 1 / 1 | 三个控制列开关 |
| `high_ci_ratio` | 2.0 | HIGH 档门槛 |
| `late[]`, `n_late` | 咖啡 14:00、运动 20:00、饮酒 21:00 | 晚间细分规则，`late_hour = NAN` 表示不细分 |
| `custom_name[8]` | NULL | 自定义 id 16–23 的名字 |

**输出 `evobeing_behavior_output_t`**

| 字段 | 含义 |
|---|---|
| `n_days_used` | 窗口内有晨起值的天数 |
| `mean_score`, `sd_score` | ȳ、s_y |
| `explained_pct` | 记录的行为合计解释的份额 %（= Σ share_pct） |
| `unexplained_pct` | 100 − explained_pct |
| `n_behaviors` | 条目数（含晚间细分项） |
| `items[]` | 每种行为一条，见下表 |
| `n_insufficient` | INSUFFICIENT 条数 |

**每条 `evobeing_behavior_item_t`**

| 字段 | 含义 |
|---|---|
| `id` | 行为 id |
| `is_late_variant` | 1 = 晚间细分项 |
| `name` | 如 `alcohol`、`caffeine_late`、自定义名 |
| `n_with`, `n_without` | 有 / 无该行为的早晨数 |
| `effect_pct` | **口径 A**：次晨精力相对本人均值的变化 %，带正负；INSUFFICIENT 时 NAN |
| `effect_points` | 同上，单位是分 |
| `raw_pct` | 简单对比 %（未控制混杂） |
| `ci_lo_pct`, `ci_hi_pct` | 90% 区间（`n_boot = 0` 时 NAN） |
| `share_pct` | 口径 B：解释份额 %，各项之和 ≤ 100 |
| `status` | 0 INSUFFICIENT / 1 LOW / 2 MID / 3 HIGH |

### 5.2 `evobeing_behavior_config_default(cfg)`、`evobeing_behavior_name_zhenghuimin(id)`

填默认配置；按 id 取英文名（自定义返回 `"custom"`）。

### 5.3 调用示例

```c
evobeing_behavior_input_t in = {0};
in.n_days = 28; in.wake_score = y; in.wake_t_h = wake; in.weekday = wd;
in.events = ev; in.n_events = n_ev;
evobeing_behavior_output_t out;
if (evobeing_behavior_impact_zhenghuimin(&in, NULL, &out) == 0)
    for (int i = 0; i < out.n_behaviors; ++i)
        if (out.items[i].status != EVOBEING_BEHAVIOR_INSUFFICIENT)
            printf("%s %+.0f%%\n", out.items[i].name, out.items[i].effect_pct);
```

`examples/behavior_demo.c` 的输出（合成数据，真实效应：饮酒 −8 分、晚睡 −5、晚间咖啡 −4、小睡 +3、运动 +3）：

```
过去 28 天，晨起精力均值 57（标准差 6.1）
行为           次数   影响       90% 区间   份额  可信度
alcohol            8/20      -16% [  -19%,   -10%]     47%  高
late_bed           7/21       -6% [   -9%,    -2%]      5%  中
nap                7/21       +6% [   -0%,   +10%]      5%  中
meditation         6/22       +6% [   +0%,    +9%]      4%  中
caffeine          20/8        +3% [   -2%,    +6%]      1%  低
caffeine_late      7/21       -2% [   -6%,    +4%]      5%  低
exercise           4/24       +2% [   -4%,    +6%]      0%  低
记录的行为合计解释 68%，其余 32% 来自睡眠本身、身体状况等未记录因素
```

## 6. 展示建议

- 列表按 |影响| 降序；每行：名字、`+6%` / `−16%`、次数 `8 次`、可信度三档用颜色深浅而不是文字。
- INSUFFICIENT 的行显示"再记录 N 次就能看到"，N = 4 − n_with（或 4 − n_without）。
- LOW 档的数字显示但置灰，附一句"方向还不确定"。
- 底部一句：记录的行为合计解释了 X%，其余来自睡眠本身、身体状况等。
- 文案用"饮酒后的早晨，精力比你平时低 16%"，不用"饮酒导致"。

## 7. 已知边界

- 行为与某个控制项几乎重合（比如 90% 的饮酒都在周五 / 周六晚），模型分不开两者，效应会被压向 0、落在 LOW 档。测试 [4] 覆盖了这种情况：控制周末后估计比不控制更接近真值，但仍不准。解决办法只有让用户多几次在工作日的记录。
- 两种行为总是同时出现（喝酒必配晚餐很晚），同样分不开，份额会在两者之间摊。
- 28 天、4 次的门槛是下限，HIGH 档一般要 6 次以上。
- 效应是"相对本人过去 28 天的均值"，用户习惯改变后（比如戒咖啡）数字会在两周左右跟着变。
- 全部数字来自合成数据测试，未在真实用户数据上验证。

## 8. 测试覆盖（`tests/test_behavior.c`）

1. 典型 28 天、20 次随机世界：饮酒方向正确 18/19，均值 −7.5 分（真值 −9）；小睡 +3.5（真值 +4）；运动 +2.4（真值 +3）；排序、份额之和、区间包含点估计、门槛。
2. 零效应行为：30 次里 27 次 |效应| < 3 分，HIGH 仅 1 次，均值 0.15 分。
3. 样本不足：出现 3 次 / 每天都出现 → INSUFFICIENT。
4. 混杂控制：饮酒只在周末、周末本身 +6 分 → 控制周末后估计更接近真值 20/20。
5. 边界：空指针、天数越界、缺天、56 天输入只用后 28 天、事件归属（23:00 与 01:00 都归下一个早晨、最后一次起床之后的不计）。
6. 自定义行为、名字、`n_boot = 0`、可复现性、不给 `wake_t_h` / `weekday`。
7. 晚间咖啡与上午咖啡分开 10/10。
