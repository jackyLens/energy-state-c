# evobeing —— 能量与精力的 C 语言实现

作者：zhenghuimin。C99，无第三方依赖，只需链接 libm。

智能戒指每 5 分钟给出一个样本（睡眠分期、心率、夜间 HRV、步数、热量、压力、血氧），
本库据此持续计算两个 0–100 的值：

- **能量**：身体的储备水平，小时到天级变化。夜间睡眠充电，日间消耗，用心率与 HRV 校正。
- **精力**：此刻的状态好坏，分钟到小时级变化。叠在能量之上，由日间心率 / HRV 相对本人常态的偏离、
  运动后回弹、小睡、以及昨夜睡眠质量（含夜间血氧平稳度）决定。

两个值都只跟本人的历史比。库在运行中学习这个人的常态（心率、HRV 基线，睡眠窗，夜间回充量，夜间血氧），
不做人与人之间的比较，也不含任何群体统计。面向消费电子场景，不是医疗判断。

## 文件

```
include/evobeing_zhenghuimin.h   主接口：两个计算函数 + 初始化 + 画像导入导出（调用方只需看这个文件）
src/evobeing_zhenghuimin.c       主接口实现（编排）
include/energy.h  src/energy.c   能量引擎（ODE + EKF，个人基线学习，缺口处理）
include/vitality.h src/vitality.c 精力引擎（偏离量 D 的标量卡尔曼滤波，事件，晨起初始化）
src/energy_internal.h            仅测试使用
src/main.c                       命令行：CSV 进 → CSV 出
examples/evobeing_demo.c         主接口调用示例（3 天合成数据）
examples/stream_demo.c           能量引擎流式示例
examples/vitality_demo.c         精力引擎示例（事件日、短睡夜、血氧不稳夜）
tests/                           单元测试
docs/                            算法与产品说明
```

## 集成

把 `src/energy.c`、`src/vitality.c`、`src/evobeing_zhenghuimin.c` 加入工程，头文件路径加 `include/`，链接 `-lm`。

```bash
cc -std=c99 -O2 -Iinclude your_app.c src/energy.c src/vitality.c src/evobeing_zhenghuimin.c -lm
```

glibc 在 `-std=c99` 下需要加 `-D_DEFAULT_SOURCE` 才能看到 `M_PI`（Makefile 已带）。

## 调用

```c
#include "evobeing_zhenghuimin.h"

evobeing_ctx_t ctx;
evobeing_init_zhenghuimin(&ctx, NULL);               /* 首次使用；有画像时传 &profile */

for (每个 5 分钟样本) {
    evobeing_sample_t s;
    evobeing_sample_init_zhenghuimin(&s);
    s.t_h      = t_hours;                            /* 距首样本的小时数，严格递增 */
    s.tod      = 13.5;                               /* 当地钟点，13:30 */
    s.stage[0] = 0; s.stage[1] = 0; s.stage[2] = 0; s.stage[3] = 0;   /* 日间无分期；夜间填 [清醒,浅睡,深睡,REM] 概率 */
    s.hr       = 72.0;                               /* bpm */
    s.hrv_log  = NAN;                                /* ln(RMSSD/ms)，日间无 RR 间期填 NAN */
    s.steps    = 120;
    s.kcal     = 6.0;
    s.stress   = 35;                                 /* 0–1 或 0–100，整条流统一 */
    s.q_hr     = 1.0;                                /* 心率质量 0–1 */
    s.q_hrv    = 1.0;
    s.spo2     = NAN;                                /* 血氧 %，日间可不测 */
    s.spo2_min = NAN;

    double energy = 0, vitality = 0;
    if (evobeing_energy_zhenghuimin(&ctx, &s, &energy) == EVOBEING_OK &&
        evobeing_vitality_zhenghuimin(&ctx, &vitality) == EVOBEING_OK) {
        printf("能量 %.0f  精力 %.0f\n", energy, vitality);
    }
}

evobeing_profile_t profile;
evobeing_get_profile_zhenghuimin(&ctx, &profile);    /* 保存；下次 init 传回即可跳过约 3 晚的学习期 */
```

每个样本先调能量函数，再调精力函数；只需要能量时不调精力函数。
用户手动标记（咖啡、冥想、小睡、压力事件、饮酒）用 `evobeing_event_zhenghuimin()`。
更多输出（置信带、四象限、标志位、晨起初始化值、本夜血氧统计）在 `ctx.energy_point` 与 `ctx.vitality_point`，
字段含义见 `energy.h` / `vitality.h`。

### 输入字段

| 字段 | 定义 | 格式 |
|---|---|---|
| `t_h` | 距首样本的小时数 | 严格递增；5 分钟一格时相差 1/12 |
| `tod` | 当地钟点 | 0–24，含小数 |
| `stage[4]` | 睡眠分期概率 [清醒, 浅睡, 深睡, REM] | 每项 0–1，和 ≤ 1；日间全 0 |
| `hr` | 心率均值 | bpm；未测到 NAN |
| `hrv_log` | ln(RMSSD / ms)，由逐拍 RR 间期算出 | 无 RR 填 NAN；毫秒值用 `energy_hrv_log_from_ms()` |
| `steps` | 本样本步数 | 数值 |
| `kcal` | 本样本活动热量 | kcal |
| `stress` | 压力评分 | 0–1 或 0–100 |
| `q_hr` / `q_hrv` | 观测质量 | 0–1；0 = 不参与计算 |
| `spo2` / `spo2_min` | 血氧均值 / 最低值 | %；未测到 NAN；只用夜间 |

负荷与活动量由库根据心率、步数、热量计算，调用方不传。

### 返回码

`EVOBEING_OK`（0）成功；`ENERGY_EDOMAIN` 时间未递增；`ENERGY_EINVAL` 参数非法或调用顺序错误。
出错时输出变量不变。

## 构建与测试

```bash
make            # 库、命令行、测试、示例
make test       # 三套单元测试
make example    # 运行示例
./build/energy_cli data/sample_panel.csv --vitality --out out.csv
```

测试用合成数据核对：状态有界、无 NaN、流式与批量一致、个人基线学习、体质平移不变性（整条心率曲线平移 ±10 bpm 能量不变）、
缺口与未佩戴、夜间血氧对晨起初始化的影响、主接口与引擎逐位一致、画像导入导出。

## 算法概要

**能量** 用两个隐状态：睡眠压力 S 与能量 E，按解析解推进（支持任意步长）：

```
S(t+h) = S∞ + (S − S∞)·e^(−b·h)
E(t+h) = E∞ + (E − E∞)·e^(−(α+β)·h)
```

睡眠期按分期权重充电，清醒期按负荷、压力、睡眠债消耗。心率与 HRV 作为观测，
用 2×2 扩展卡尔曼滤波校正（Student-t 稳健化、4σ 门控）。观测方程的截距是个人基线，由库在线学习
（首晚引导 → 前 3 晚累计平均 → 之后 EWMA λ=0.1，异常夜降权），所以同一条心率曲线整体平移不改变能量。
输出 0–100 时叠加昼夜节律调制。

**精力** `V = clip(V* + D)`，均衡线 `V* = E·(1 + ρ·c(t))`，偏离量 D 按 `dD/dt = −D/τ(E) + 事件脉冲` 推进，
能量越低回弹越慢，且 `D ≤ 0.30·E`。观测是能量模型没解释掉的日间心率 / HRV 残差：
心率略高于常态是"打起精神"（精力升），持续大幅高于常态是"过载"（精力降），低于常态是"犯困"。
每天醒来时 D 重置为晨起初始化值，由本夜能量回充量、醒前 HRV、睡眠结构、夜间血氧平稳度相对本人常态的偏差给出，
`|D_init| ≤ 0.2`。日间血氧不参与计算：健康人静息血氧稳定在 96–99%，变化在传感器噪声之内。

详细推导、参数表与图见 `docs/`。

## 文档

| 文档 | 内容 |
|---|---|
| [docs/technical.html](docs/technical.html) · [PDF](docs/technical.pdf) | 能量：数据流、收支方程、观测与 EKF、流式时序、个人基线 |
| [docs/product.html](docs/product.html) · [PDF](docs/product.pdf) | 能量：曲线怎么读、刻度语义、能说 / 不能说 |
| [docs/vitality_technical.html](docs/vitality_technical.html) · [PDF](docs/vitality_technical.pdf) | 精力：级联结构、倒 U 观测、事件、晨起初始化、夜间血氧 |
| [docs/vitality_product.html](docs/vitality_product.html) · [PDF](docs/vitality_product.pdf) | 精力：能量与精力并排怎么读、四象限、手动标记 |
| [docs/vitality_design.md](docs/vitality_design.md) | 精力：模型推导、参数表、决策记录 |

## 约定

- 状态 S/E 恒在 [0,1]，能量与精力恒在 [0,100]；无随机性，同一输入结果逐位可复现。
- 无全局状态：每个用户一个 `evobeing_ctx_t`（约 6 KB，以 `sizeof` 为准），可多实例、多线程（每线程各自的 ctx）。
- 流式路径不分配堆内存。
- 数据中断（摘下、没电）后直接喂下一个样本即可；超过 1 小时的缺口按个人常规睡眠窗 / 清醒静息推进，
  输出带 `ENERGY_FLAG_GAP`。
- 学习期（前 3 晚）输出带 `ENERGY_FLAG_LEARNING` / `VITALITY_FLAG_LEARNING`，可提示"正在了解你的常态"。
- 所有数字来自合成数据的测试，参数为先验值；真实用户数据接入后应按 `docs/vitality_design.md` 末尾列出的问题复核。
