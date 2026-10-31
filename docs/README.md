# 说明文档

| 文档 | 面向 | HTML（浏览器打开） | PDF（GitLab 可直接预览） |
|---|---|---|---|
| 产品说明 | 产品 / 设计 / 运营 | [product.html](product.html) | [product.pdf](product.pdf) |
| 技术说明 | 算法 / 后端 / 客户端 | [technical.html](technical.html) | [technical.pdf](technical.pdf) |
| 精力 · 产品说明 | 产品 / 设计 / 运营 | [vitality_product.html](vitality_product.html) | [vitality_product.pdf](vitality_product.pdf) |
| 精力 · 技术说明 | 算法 / 后端 / 客户端 | [vitality_technical.html](vitality_technical.html) | [vitality_technical.pdf](vitality_technical.pdf) |
| 精力 · 算法设计 | 算法 | [vitality_design.md](vitality_design.md)（Mermaid 流程图，GitHub 直接渲染） | — |

这些文档都是**自包含 HTML**（内联 SVG、无外部图片/字体请求），直接双击用浏览器打开即可；PDF 由同一份 HTML 打印生成，适合在 GitLab 里直接预览或转发。

## 图

技术说明里的四张图：

1. **输入 → 输出数据流**（七个阶段：采集 / 装配 / 派生 / 状态 / 校正 / 呈现 / 交付）
2. **能量收支**（α 充电、β 消耗、S 同时作用于两端）
3. **流式调用时序**（一次 `energy_ctx_step` 内的预测与观测更新）
4. **运行状态机**（冷启动 / 稳定 / 降级 / 复位）

产品说明里的三张图：

1. **一天的能量曲线**（真实输出：夜间充电、晨间峰值、日间消耗、95% 带）
2. **睡眠分 → 每日能量 → 行为关联**（三层价值链）
3. **0–100 刻度**（低 / 常见 / 接近满格，以及“0 分结构上不可达”）

配色是按语义分配的六种色相（采集=靛蓝、装配=紫、派生=琥珀、状态=赤陶（焦点）、滤波=玫红、呈现=绿），不是装饰性彩虹；正文表格里的色点与图例一致。

精力技术说明里的五张图：级联结构、一步流水线、倒 U 变换 u(z)、两种 0–100 刻度；产品说明里的五张图：事件日、短睡日与血氧不稳夜的能量/精力双轨曲线（`build/vitality_demo.csv` 真实输出）、四象限、刻度。


## 维护约定

- 文档里的数字（测试断言数等）与仓库实测保持一致；改动后同步更新技术说明与 `gen_vitality_docs.py` 的 `NUMBERS`。
- 图里的数值来自 `data/sample_panel.csv` 的真实输出；换样例数据后需要重新生成曲线路径。
