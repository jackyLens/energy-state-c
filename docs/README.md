# 说明文档

| 文档 | 面向 | HTML（浏览器打开） | PDF（GitLab 可直接预览） |
|---|---|---|---|
| 产品说明 | 产品 / 设计 / 运营 | [product.html](product.html) | [product.pdf](product.pdf) |
| 技术说明 | 算法 / 后端 / 客户端 | [technical.html](technical.html) | [technical.pdf](technical.pdf) |

两份文档都是**自包含 HTML**（内联 SVG、无外部图片/字体请求），直接双击用浏览器打开即可；PDF 由同一份 HTML 打印生成，适合在 GitLab 里直接预览或转发。

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

## 重新生成 PDF

需要 Chrome / Chromium：

```bash
make docs                     # 自动找 Chrome/Chromium
make docs CHROME=/path/to/chrome
```

或手动：

```bash
"/Applications/Google Chrome.app/Contents/MacOS/Google Chrome" \
  --headless=new --no-pdf-header-footer \
  --print-to-pdf=docs/technical.pdf docs/technical.html
```

## v1.2 更新

技术说明 §04 新增"个体化观测模型"（个人基线的引导、冷启动、EWMA 更新、异常段降权、持久化）以及日间方程、缺口、流式负荷、压力量纲的修正；产品说明的"边界"一节改为"只跟自己比"的口径。文档中的 v1.2 数字来自合成数据，未在真实用户数据上验证。

## 维护约定

- 文档里的数字（测试断言数、覆盖率、一致性偏差、性能）必须与仓库实测一致；改动后同步更新 **技术说明 §7** 与首页统计。
- 图里的数值来自 `data/sample_panel.csv` 的真实输出；换样例数据后需要重新生成曲线路径。
- 修改文案时不要出现“测量 / 诊断 / 治疗”等措辞（见产品说明 §6）。
