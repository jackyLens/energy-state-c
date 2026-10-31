#!/usr/bin/env python3
"""生成精力层的两份自包含 HTML 文档（技术说明 / 产品说明）。

用法：
    make build/vitality_demo
    ./build/vitality_demo --csv build/vitality_demo.csv > /dev/null
    python3 scripts/gen_vitality_docs.py build/vitality_demo.csv

- 样式块直接复用 docs/technical.html / docs/product.html 的 <style>，保证两代文档视觉一致；
- 曲线图里的路径来自 vitality_demo 的真实 CSV 输出（第 4 天事件日 / 第 5 天短睡日）；
- 数字（测试断言数、结构体大小）由下方 NUMBERS 字典集中维护，改动后同步更新 docs/README.md。
"""
import csv
import html
import math
import os
import re
import sys

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
DOCS = os.path.join(ROOT, "docs")

# 实测数字（make test，合成数据）
NUMBERS = {
    "energy_checks": "81,794",
    "vitality_checks": "40,489",
    "evobeing_checks": "12,121",
    "ctx_bytes": "3984",
    "ectx_bytes": "2352",
    "evobeing_ctx_bytes": "6736",
    "version": "2.0.0",
}

FLAG = dict(SLEEP=1, LEARNING=2, GAP=4, WAKE=8, EXERCISE=16, NAP=32, AROUSAL=64, STRESS=128,
            FATIGUE=256, USER=512, NO_OBS=1024, PCT_READY=2048, BREATH=4096, SPO2=8192)


def read_style(name):
    text = open(os.path.join(DOCS, name), encoding="utf-8").read()
    m = re.search(r"<style>(.*?)</style>", text, re.S)
    return m.group(1)


def load_csv(path):
    rows = []
    with open(path, encoding="utf-8") as f:
        for r in csv.DictReader(f):
            rows.append({k: (float(v) if k != "flags" else int(v)) for k, v in r.items() if v != ""})
    return rows


def fmt(v, nd=1):
    return f"{v:.{nd}f}"


# ---------------------------------------------------------------- SVG helpers
def box(x, y, w, h, fill, stroke, title, lines=(), r=6, tw=600):
    out = [f'<rect x="{x}" y="{y}" width="{w}" height="{h}" rx="{r}" fill="{fill}" stroke="{stroke}"/>']
    out.append(f'<text x="{x + 12}" y="{y + 20}" font-size="13" font-weight="{tw}" fill="#23252f">{html.escape(title)}</text>')
    for i, ln in enumerate(lines):
        out.append(f'<text x="{x + 12}" y="{y + 38 + i * 15}" font-size="11.5" fill="#5b6070">{html.escape(ln)}</text>')
    return "\n".join(out)


def arrow(x1, y1, x2, y2, color="#7a7f8c", dash=None, marker="arr"):
    d = f' stroke-dasharray="{dash}"' if dash else ""
    return f'<path d="M{x1},{y1} L{x2},{y2}" fill="none" stroke="{color}" stroke-width="1.3"{d} marker-end="url(#{marker})"/>'


def path_arrow(d, color="#7a7f8c", dash=None, marker="arr"):
    dd = f' stroke-dasharray="{dash}"' if dash else ""
    return f'<path d="{d}" fill="none" stroke="{color}" stroke-width="1.3"{dd} marker-end="url(#{marker})"/>'


def defs():
    return ('<defs>'
            '<marker id="arr" markerWidth="8" markerHeight="6" refX="7" refY="3" orient="auto"><polygon points="0 0, 8 3, 0 6" fill="#7a7f8c"/></marker>'
            '<marker id="arr-acc" markerWidth="8" markerHeight="6" refX="7" refY="3" orient="auto"><polygon points="0 0, 8 3, 0 6" fill="#bf4520"/></marker>'
            '</defs>')


FONT = "font-family=\"'PingFang SC','Microsoft YaHei','Noto Sans SC',sans-serif\""


def svg_wrap(w, h, body, title, desc):
    return (f'<svg viewBox="0 0 {w} {h}" role="img" aria-label="{html.escape(title)}" {FONT}>'
            f'<title>{html.escape(title)}</title><desc>{html.escape(desc)}</desc>{defs()}'
            f'<rect width="{w}" height="{h}" fill="#ffffff"/>{body}</svg>')


# ---------------------------------------------------------------- FIG: 级联结构
def fig_cascade():
    b = []
    # 输入
    b.append(box(30, 40, 150, 120, "#eef1ff", "#3f5bd9", "每 5 分钟样本",
                 ["心率 / HRV", "步数 / 分期", "质量 / 压力"]))
    # 能量层
    b.append(f'<rect x="220" y="20" width="300" height="290" rx="8" fill="none" stroke="#bf4520" stroke-dasharray="4,3"/>')
    b.append(f'<text x="232" y="40" font-size="11" fill="#bf4520" letter-spacing=".12em">能量层 energy.c（不改）</text>')
    b.append(box(240, 56, 260, 60, "#fdeee8", "#bf4520", "状态方程 S / E", ["解析解推进，睡眠恢复、日间消耗"]))
    b.append(box(240, 132, 260, 60, "#fbe9ef", "#a83f63", "个体化观测模型 + EKF", ["预测 HR / HRV，用残差校正 E"]))
    b.append(box(240, 208, 120, 84, "#fbe9ef", "#a83f63", "输出（残差）", ["hr_resid", "hrv_resid"]))
    b.append(box(380, 208, 120, 84, "#fdeee8", "#bf4520", "输出", ["E, sleep_prob", "load, activity（内部计算）"]))
    # 精力层
    b.append(f'<rect x="560" y="20" width="330" height="290" rx="8" fill="none" stroke="#0e7a83" stroke-dasharray="4,3"/>')
    b.append(f'<text x="572" y="40" font-size="11" fill="#0e7a83" letter-spacing=".12em">精力层 vitality.c</text>')
    b.append(box(580, 56, 290, 52, "#e6f4f5", "#0e7a83", "① 预测  D ← e^(−h/τ_D)·D", ["τ_D = 0.75 h·(1+1.5(1−E))"]))
    b.append(box(580, 118, 290, 52, "#e6f4f5", "#0e7a83", "② 睡眠段 / 晨起初始化", ["D ← D_init（相对本人常态）"]))
    b.append(box(580, 180, 290, 52, "#e6f4f5", "#0e7a83", "③④ 残差 z 化 → 事件 → 标量 KF", ["倒 U 变换；运动 / 小睡 / 拖拽"]))
    b.append(box(580, 242, 290, 52, "#e8f4ec", "#2c7a4b", "⑤⑥ 裁剪 → 分位 → 输出", ["V = clip(V* + D)，percentile"]))
    # 用户标记
    b.append(box(30, 220, 150, 70, "#f2eeff", "#6d4fc4", "用户标记（可选）", ["咖啡 / 冥想 / 压力", "vitality_ctx_event"]))
    # arrows
    b.append(arrow(180, 100, 240, 86))
    b.append(arrow(370, 116, 370, 132))
    b.append(arrow(300, 192, 300, 208))
    b.append(arrow(440, 192, 440, 208))
    b.append(path_arrow("M500,250 H540 V206 H580", "#a83f63"))
    b.append(path_arrow("M500,270 H548 V82 H580", "#bf4520", marker="arr-acc"))
    b.append(path_arrow("M180,255 H200 V330 H720 V294", "#6d4fc4", dash="5,4"))
    b.append(f'<text x="206" y="345" font-size="11" fill="#6d4fc4">注入到下一步的预测阶段（咖啡 20 min 一阶起效）</text>')
    return svg_wrap(920, 360, "\n".join(b), "级联结构：能量层不变，精力层只消费它的输出与残差",
                    "样本进入能量层得到 E 与观测残差；精力层按预测、晨起初始化、观测更新、裁剪输出四步计算 V。")


# ---------------------------------------------------------------- FIG: 一步流程
def fig_step():
    b = []
    x0, w, h = 40, 250, 54
    col2 = 340
    col3 = 640
    rows = [
        (x0, 30, "energy_ctx_step", ["E · sleep_prob · load · activity · 残差"], "#fdeee8", "#bf4520"),
        (x0, 110, "① 预测", ["D ← k·D，P ← k²P + q_d·h", "待注入的用户脉冲在此加入"], "#e6f4f5", "#0e7a83"),
        (x0, 190, "② 睡眠段跟踪", ["sleep_prob ≥ 0.5 → 累计；无观测", "连续 3 个清醒样本 → 段结束"], "#e6f4f5", "#0e7a83"),
        (x0, 270, "段结束分流", ["≥ 3 h：晨起初始化 WAKE", "15 min–3 h 且 9–20 点：小睡 NAP", "< 2 h 内再睡：续睡合并"], "#f2eeff", "#6d4fc4"),
        (col2, 110, "③ z 化", ["z_hr = hr_resid / 7", "z_hrv = hrv_resid / 0.40"], "#fbf1e2", "#a8650f"),
        (col2, 190, "观测权重 w_act", ["A ≤ 0.10 → 1；A ≥ 0.35 → 0", "L > 0.60：A 高则 0，A 低则 ×0.5"], "#fbf1e2", "#a8650f"),
        (col2, 270, "事件检测", ["运动段（L≥0.5 且 A≥0.35）", "无分期小睡 / 低 HRV 拖拽", "唤醒·过载标志"], "#f2eeff", "#6d4fc4"),
        (col3, 110, "④ 标量 KF ×2", ["y = u(z_hr)，H = 3.0，R = 1.0 / w_act", "y = z_hrv，H = 2.0，R = 2.0 / w_act", "门控 |ν|/√S > 4 丢弃"], "#fbe9ef", "#a83f63"),
        (col3, 190, "⑤ 裁剪", ["−V* ≤ D ≤ 0.30·E"], "#e8f4ec", "#2c7a4b"),
        (col3, 270, "⑥ 分位 + 输出", ["日间 V 记入 50 桶直方图", "vitality / percentile / score / 带 / 象限"], "#e8f4ec", "#2c7a4b"),
    ]
    for (x, y, t, ls, f, s) in rows:
        hh = 38 + 15 * len(ls)
        b.append(box(x, y, w, hh, f, s, t, ls))
    b.append(arrow(x0 + w / 2, 84, x0 + w / 2, 110))
    b.append(arrow(x0 + w / 2, 177, x0 + w / 2, 190))
    b.append(arrow(x0 + w / 2, 257, x0 + w / 2, 270))
    b.append(path_arrow(f"M{x0 + w},140 H{col2}"))
    b.append(arrow(col2 + w / 2, 177, col2 + w / 2, 190))
    b.append(arrow(col2 + w / 2, 257, col2 + w / 2, 270))
    b.append(path_arrow(f"M{col2 + w},140 H{col3}"))
    b.append(arrow(col3 + w / 2, 192, col3 + w / 2, 190 + 0))
    b.append(arrow(col3 + w / 2, 177 + 15, col3 + w / 2, 190))
    b.append(arrow(col3 + w / 2, 243, col3 + w / 2, 270))
    b.append(path_arrow(f"M{x0 + w},300 H{col2 - 20} V 380 H{col3 + w / 2} V 350", "#6d4fc4", dash="5,4"))
    b.append(f'<text x="{col2}" y="372" font-size="11" fill="#6d4fc4">晨起初始化 / 小睡脉冲直接改 D，再进入裁剪与输出</text>')
    b.append(f'<text x="{x0}" y="395" font-size="11" fill="#5b6070">睡眠中（SLEEP）与缺口（GAP）样本跳过 ③④，直接到 ⑤⑥：没有观测时 V 退化为均衡线。</text>')
    return svg_wrap(920, 410, "\n".join(b), "一次 vitality_ctx_step 的六个阶段",
                    "预测、睡眠段跟踪、z 化、权重、事件检测、两次标量卡尔曼更新、裁剪、分位与输出。")


# ---------------------------------------------------------------- FIG: 倒 U
def arousal(z, z_c=2.8, drowsy=0.6, y_min=-2.0):
    u = z * (1 - z / z_c) if z >= 0 else drowsy * z
    return max(y_min, min(u, z_c / 4))


def fig_inverted_u():
    W, H = 920, 300
    px0, px1, py0, py1 = 70, 880, 30, 250
    zmin, zmax, umin, umax = -3.5, 4.5, -2.2, 0.9

    def X(z):
        return px0 + (z - zmin) / (zmax - zmin) * (px1 - px0)

    def Y(u):
        return py1 - (u - umin) / (umax - umin) * (py1 - py0)

    b = [f'<line x1="{px0}" y1="{Y(0)}" x2="{px1}" y2="{Y(0)}" stroke="#c9c4ba"/>',
         f'<line x1="{X(0)}" y1="{py0}" x2="{X(0)}" y2="{py1}" stroke="#c9c4ba"/>']
    pts = []
    z = zmin
    while z <= zmax + 1e-9:
        pts.append(f"{X(z):.1f},{Y(arousal(z)):.1f}")
        z += 0.05
    b.append(f'<polyline points="{" ".join(pts)}" fill="none" stroke="#0e7a83" stroke-width="2.2"/>')
    for z, lab in [(-2, "−14 bpm"), (-1, "−7"), (0, "0"), (1, "+7"), (1.4, ""), (2, "+14"), (2.8, "+20 bpm"), (4, "+28")]:
        b.append(f'<line x1="{X(z)}" y1="{Y(0) - 3}" x2="{X(z)}" y2="{Y(0) + 3}" stroke="#5b6070"/>')
        if lab:
            b.append(f'<text x="{X(z)}" y="{Y(0) + 16}" font-size="11" fill="#5b6070" text-anchor="middle">{lab}</text>')
    b.append(f'<text x="{px1}" y="{Y(0) + 32}" font-size="11" fill="#5b6070" text-anchor="end">日间心率残差 z_hr（σ = 7 bpm）</text>')
    # annotations
    b.append(f'<rect x="{X(0.7)}" y="{py0}" width="{X(2.1) - X(0.7)}" height="{py1 - py0}" fill="#e6f4f5" opacity=".6"/>')
    b.append(f'<text x="{(X(0.7) + X(2.1)) / 2}" y="{py1 - 8}" font-size="11" fill="#0e7a83" text-anchor="middle">轻度唤醒 AROUSAL [0.7, 2.1]</text>')
    b.append(f'<rect x="{X(2.1)}" y="{py0}" width="{px1 - X(2.1)}" height="{py1 - py0}" fill="#fdeee8" opacity=".6"/>')
    b.append(f'<text x="{(X(2.1) + px1) / 2}" y="{py1 - 8}" font-size="11" fill="#bf4520" text-anchor="middle">过载 STRESS（z &gt; 2.1，持续 10 min）</text>')
    b.append(f'<text x="{X(-2.6)}" y="{py0 + 14}" font-size="11" fill="#a8650f" text-anchor="middle">犯困（斜率 0.6）</text>')
    b.append(f'<circle cx="{X(1.4)}" cy="{Y(0.7)}" r="4" fill="#0e7a83"/>')
    b.append(f'<text x="{X(1.4)}" y="{Y(0.7) - 10}" font-size="11" fill="#0e7a83" text-anchor="middle">峰值 z = 1.4（≈ +10 bpm），u = 0.7</text>')
    b.append(f'<circle cx="{X(2.8)}" cy="{Y(0)}" r="4" fill="#bf4520"/>')
    b.append(f'<text x="{X(2.8) + 8}" y="{Y(0) - 8}" font-size="11" fill="#bf4520">零点 z_c = 2.8：再高就是负收益</text>')
    b.append(f'<text x="{px0 - 4}" y="{Y(0.8)}" font-size="11" fill="#5b6070" text-anchor="end">u(z)</text>')
    b.append(f'<text x="{px0 - 4}" y="{Y(-2) + 4}" font-size="11" fill="#5b6070" text-anchor="end">−2</text>')
    return svg_wrap(W, H, "\n".join(b), "倒 U 变换：心率略高于本人常态抬精力，高太多转为过载",
                    "横轴是日间心率残差的 z 值，纵轴是进入卡尔曼更新的观测量 u。")


# ---------------------------------------------------------------- FIG: 双轨曲线（真实输出）
def fig_dual_track(rows, day, title, annotate):
    W, H = 920, 360
    px0, px1, py0, py1 = 60, 890, 24, 290
    day_rows = [r for r in rows if day * 24 <= r["t_h"] < (day + 1) * 24]

    def X(tod):
        return px0 + tod / 24 * (px1 - px0)

    def Y(v):
        return py1 - v / 100 * (py1 - py0)

    b = []
    # 睡眠底色
    segs, cur = [], None
    for r in day_rows:
        s = bool(r["flags"] & FLAG["SLEEP"])
        if s and cur is None:
            cur = r["tod"]
        if not s and cur is not None:
            segs.append((cur, r["tod"]))
            cur = None
    if cur is not None:
        segs.append((cur, 24))
    for a, c in segs:
        b.append(f'<rect x="{X(a):.1f}" y="{py0}" width="{X(c) - X(a):.1f}" height="{py1 - py0}" fill="#eef1ff"/>')
    # grid
    for v in (0, 25, 50, 75, 100):
        b.append(f'<line x1="{px0}" y1="{Y(v)}" x2="{px1}" y2="{Y(v)}" stroke="#ece7dd"/>')
        b.append(f'<text x="{px0 - 8}" y="{Y(v) + 4}" font-size="11" fill="#71768a" text-anchor="end">{v}</text>')
    for hh in range(0, 25, 3):
        b.append(f'<text x="{X(hh)}" y="{py1 + 16}" font-size="11" fill="#71768a" text-anchor="middle">{hh:02d}:00</text>')
    # band
    up = " ".join(f"{X(r['tod']):.1f},{Y(r['upper']):.1f}" for r in day_rows)
    lo = " ".join(f"{X(r['tod']):.1f},{Y(r['lower']):.1f}" for r in reversed(day_rows))
    b.append(f'<polygon points="{up} {lo}" fill="#0e7a83" opacity=".08"/>')
    # curves
    for key, color, width, dash in (("energy", "#bf4520", 1.8, ""), ("v_star", "#71768a", 1.2, ' stroke-dasharray="4,3"'), ("vitality", "#0e7a83", 2.2, "")):
        pts = " ".join(f"{X(r['tod']):.1f},{Y(100 * r[key] if key == 'v_star' else r[key]):.1f}" for r in day_rows)
        b.append(f'<polyline points="{pts}" fill="none" stroke="{color}" stroke-width="{width}"{dash}/>')
    # annotations
    for tod, label, color in annotate:
        r = min(day_rows, key=lambda q: abs(q["tod"] - tod))
        b.append(f'<circle cx="{X(r["tod"]):.1f}" cy="{Y(r["vitality"]):.1f}" r="4" fill="{color}"/>')
        b.append(f'<text x="{X(r["tod"]) + 6:.1f}" y="{Y(r["vitality"]) - 8:.1f}" font-size="11" fill="{color}">{html.escape(label)} {r["vitality"]:.0f}</text>')
    # legend
    lx = px0 + 10
    b.append(f'<rect x="{lx}" y="{py1 + 26}" width="18" height="3" fill="#bf4520"/><text x="{lx + 24}" y="{py1 + 31}" font-size="11" fill="#5b6070">能量</text>')
    b.append(f'<rect x="{lx + 70}" y="{py1 + 26}" width="18" height="3" fill="#71768a"/><text x="{lx + 94}" y="{py1 + 31}" font-size="11" fill="#5b6070">均衡线 V*</text>')
    b.append(f'<rect x="{lx + 180}" y="{py1 + 26}" width="18" height="3" fill="#0e7a83"/><text x="{lx + 204}" y="{py1 + 31}" font-size="11" fill="#5b6070">精力（含 95% 带）</text>')
    b.append(f'<rect x="{lx + 330}" y="{py1 + 23}" width="18" height="9" fill="#eef1ff"/><text x="{lx + 354}" y="{py1 + 31}" font-size="11" fill="#5b6070">睡眠</text>')
    b.append(f'<text x="{px1}" y="{py1 + 31}" font-size="11" fill="#71768a" text-anchor="end">{html.escape(title)}</text>')
    return svg_wrap(W, H, "\n".join(b), title, "来自 examples/vitality_demo.c 的真实输出：能量、均衡线与精力三条曲线。")


# ---------------------------------------------------------------- FIG: 四象限
def fig_quadrants():
    W, H = 920, 330
    b = []
    cells = [
        (260, 30, "#e8f4ec", "#2c7a4b", "能量高 · 精力高", "状态极佳", "适合安排需要专注或体力的事。"),
        (580, 30, "#fbf1e2", "#a8650f", "能量高 · 精力低", "暂时低迷", "储备还在，只是此刻提不起劲：走动、喝水、小睡 15 分钟。"),
        (260, 170, "#e6f4f5", "#0e7a83", "能量低 · 精力高", "别透支", "状态不错但储备不足，今晚早点睡。"),
        (580, 170, "#fdeee8", "#bf4520", "能量低 · 精力低", "建议休息", "两条线都低，别硬撑。"),
    ]
    for x, y, f, s, t, h2, d in cells:
        b.append(f'<rect x="{x}" y="{y}" width="300" height="120" rx="8" fill="{f}" stroke="{s}"/>')
        b.append(f'<text x="{x + 16}" y="{y + 26}" font-size="11" fill="{s}" letter-spacing=".1em">{t}</text>')
        b.append(f'<text x="{x + 16}" y="{y + 56}" font-size="20" fill="#23252f" font-family="Songti SC, Noto Serif SC, serif">{h2}</text>')
        b.append(f'<text x="{x + 16}" y="{y + 86}" font-size="11.5" fill="#5b6070">{d}</text>')
    b.append('<text x="60" y="95" font-size="12" fill="#5b6070">精力 ≥ 50</text>')
    b.append('<text x="60" y="235" font-size="12" fill="#5b6070">精力 &lt; 50</text>')
    b.append('<text x="410" y="312" font-size="12" fill="#5b6070" text-anchor="middle">能量 ≥ 50</text>')
    b.append('<text x="730" y="312" font-size="12" fill="#5b6070" text-anchor="middle">能量 &lt; 50</text>')
    return svg_wrap(W, H, "\n".join(b), "四象限：能量 × 精力的四种组合与文案",
                    "横轴能量高低，纵轴精力高低；每格一句建议。阈值 50 可配置。")


# ---------------------------------------------------------------- FIG: 两种刻度
def fig_scales():
    W, H = 920, 230
    b = []
    b.append(box(40, 30, 400, 170, "#e6f4f5", "#0e7a83", "vitality = 100 × V（绝对自我刻度）",
                 ["立即可用", "E 与 D 本身已相对个人基线", "适合曲线、带宽、与能量同图", "示例：午睡后 55.6 → 68.6"]))
    b.append(box(480, 30, 400, 170, "#f2eeff", "#6d4fc4", "percentile（本人历史分位）",
                 ["需要 ≥ 3 个日间（PCT_READY）", "过去约 28 天日间 V 的经验分位（半衰期 14 天）", "文案：“比你 80% 的时候好”", "score_mode 选择 score 列取哪个"]))
    b.append('<text x="460" y="120" font-size="22" fill="#71768a" text-anchor="middle">×</text>')
    b.append('<text x="460" y="215" font-size="11.5" fill="#bf4520" text-anchor="middle">两种刻度都只跟自己比，库内不存在任何跨人统计</text>')
    return svg_wrap(W, H, "\n".join(b), "两种 0–100 刻度", "绝对自我刻度与个人分位刻度的定义与适用场景。")


# ---------------------------------------------------------------- 文档骨架
def page(title_tag, style, eyebrow, h1, lede, chips, toc, sections, footer_l, footer_r, cover=""):
    toc_html = "".join(f'<a href="#{i}">{t}</a>' for i, t in toc)
    chips_html = "".join(f'<span class="chip">{c}</span>' for c in chips)
    return f"""<!DOCTYPE html>
<html lang="zh-CN">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>{html.escape(title_tag)}</title>
<style>{style}</style>
</head>
<body>
<div class="wrap">
  <header>
    <div class="eyebrow">{eyebrow}</div>
    <h1>{h1}</h1>
    <p class="lede">{lede}</p>
    <div class="meta">{chips_html}</div>
  </header>
{cover}
  <nav class="toc"><span>CONTENTS</span>{toc_html}</nav>
{sections}
  <footer>
    <span>{footer_l}</span>
    <span>{footer_r}</span>
  </footer>
</div>
</body>
</html>
"""


def section(sid, no, title, body):
    return f"""
  <section id="{sid}">
    <div class="sec-head"><span class="sec-no">{no}</span><h2>{title}</h2></div>
    <div class="rule"></div>
{body}
  </section>
"""


def figure(no, title, svg, caption):
    return f"""
    <figure>
      <div class="fig-head"><span class="fig-no">FIG {no}</span><span class="fig-title">{title}</span></div>
      <div class="diagram-container">{svg}</div>
      <figcaption>{caption}</figcaption>
    </figure>
"""


def pick(rows, day, tod):
    return min((r for r in rows if day * 24 <= r["t_h"] < (day + 1) * 24), key=lambda r: abs(r["tod"] - tod))


# ---------------------------------------------------------------- 技术说明
def technical(rows):
    N = NUMBERS
    nap_before, nap_after = pick(rows, 3, 13.25), pick(rows, 3, 13.5)
    ex_before, ex_after = pick(rows, 3, 18.5), pick(rows, 3, 18.5833)
    ar0, ar1 = pick(rows, 3, 9.25), pick(rows, 3, 10.25)
    caf0, caf1 = pick(rows, 3, 14.75), pick(rows, 3, 15.25)
    wake4 = next(r for r in rows if r["flags"] & FLAG["WAKE"] and 96 <= r["t_h"] < 120)
    wake3 = next(r for r in rows if r["flags"] & FLAG["WAKE"] and 72 <= r["t_h"] < 96)
    wake5 = next(r for r in rows if r["flags"] & FLAG["WAKE"] and 120 <= r["t_h"] < 144)

    s1 = section("why", "01", "为什么要第二个状态量", f"""
    <p>能量 E 是慢变量：它回答"储备还剩多少"，由夜间恢复和日间消耗决定，一杯咖啡或一次午睡不会改变储备。但用户此刻的主观状态会变。精力层把这类<b>短时的、会自己消退的</b>偏离单独建模为 D，叠加在"能量 × 昼夜节律"的均衡线上：</p>
    <pre>V(t)  = clip( V*(t) + D(t), 0, 1 )
V*(t) = E(t) · [ 1 + ρ · c(t) ]                       <span class="c">ρ = 0.15，c(t) 复用能量层的节律形状</span>
dD/dt = − D / τ_D(E) + Σ 事件脉冲 + 观测修正          <span class="c">τ_D = 0.75 h · (1 + 1.5 · (1 − E))</span></pre>
    <table>
      <thead><tr><th>现象</th><th>能量 E</th><th>精力 V</th></tr></thead>
      <tbody>
        <tr><td>下午咖啡，心率略升</td><td>继续下降</td><td>短暂上浮后约 1 h 内衰回</td></tr>
        <tr><td>20 min 午睡</td><td>略升</td><td>立即回弹（示例 {fmt(nap_before['vitality'])} → {fmt(nap_after['vitality'])}）</td></tr>
        <tr><td>1 h 会议，心率 +22 bpm 未动</td><td>消耗加快，滞后</td><td>过载，下压</td></tr>
        <tr><td>比平时睡得好 / 差</td><td>晨起 E 高 / 低</td><td>醒来高于 / 低于均衡线（D_init）</td></tr>
      </tbody>
    </table>
    <div class="note"><b>没有观测时它退化为均衡线。</b>D 的开环解是指数衰减，没有事件、没有残差时 |D| 单调不增（测试 <code class="inl">test_open_loop</code>），因此精力层不会"无中生有"。</div>
""")
    s2 = section("arch", "02", "级联结构：能量层不变，精力层只吃它的输出", f"""
    <p>D 没有并入能量层的状态向量。能量层已经在学个人基线，如果再加一维 D，"心率高是基线偏了还是 D 高"不可辨识。级联后能量层代码不变，只新增两列输出：观测更新前的 HR / HRV 残差 <code class="inl">hr_resid / hrv_resid</code>——能量层解释不了的部分，正是精力层要解释的部分。</p>
{figure(1, "级联：能量层 → 残差 → 精力层", fig_cascade(), "能量层的残差已扣除个人基线、昼夜节律、活动负荷与 E 自身的贡献，所以体质差异不进入精力：整体心率 +8 bpm，第 4 天起精力分差 < 3（test_constitution_invariance）。")}
    <div class="cols">
      <div class="card"><h4><span class="dot" style="background:#0e7a83"></span>精力层新增的量</h4>
        <ul><li>隐状态 D 与其方差 P（标量）</li><li>个人常态 <code>vitality_habit_t</code>：夜间回充量、醒前 HRV、睡眠时长 / 深睡 / 碎片化、夜间血氧均值 / 下跌占比（EWMA，可持久化）</li><li>本夜血氧 1% 分箱直方图（31 箱）</li><li>50 桶的日间 V 分位直方图（按日历日指数衰减）</li></ul>
      </div>
      <div class="card"><h4><span class="dot" style="background:#bf4520"></span>从能量层拿的量</h4>
        <ul><li>E、e_sd（折进精力带宽）</li><li>sleep_prob（睡眠段判定）、load、activity（观测权重 / 事件）</li><li>hr_resid、hrv_resid（观测）</li><li>昼夜节律形状 <code>energy_circadian()</code>（均衡线）</li></ul>
      </div>
    </div>
""")
    s3 = section("step", "03", "一步之内：预测、晨起、观测、裁剪、分位", f"""
    <p>流式接口 <code class="inl">vitality_ctx_step</code> 内部先走一次 <code class="inl">energy_ctx_step</code>，然后按六个阶段推进 D。批量接口 <code class="inl">vitality_compute</code> 是同一套代码逐点喂入，再加 15 min 的展示平滑（只平滑 vitality / score / 带，不改状态 v / d）；两者的 v / d 逐点一致到 1e-9（<code class="inl">test_stream_equals_batch</code>）。</p>
{figure(2, "一次 vitality_ctx_step 的六个阶段", fig_step(), "睡眠段与缺口样本没有观测，直接进入裁剪与输出。晨起初始化与小睡脉冲直接改写 D。")}
    <h3>预测</h3>
    <pre>k = exp(−h / τ_D(E))        D ← k·D        P ← k²·P + q_d·h        <span class="c">q_d = 0.02 /h，P0 = 0.01（晨起先验标准差 0.1）</span></pre>
    <p class="small">稳态 P ≈ q_d·τ_D/2 ≈ 0.0075，D 的先验标准差约 0.09；能量低时 τ_D 变长（最长 112 min），疲惫状态下任何波动都恢复得慢。</p>
    <h3>有界</h3>
    <pre>−V* ≤ D ≤ d_cap · E          <span class="c">d_cap = 0.30：E = 0.4 时精力最多比均衡线高 0.12</span></pre>
""")
    s4 = section("obs", "04", "观测：残差 z 化、倒 U 变换、标量卡尔曼", f"""
    <pre>z_hr  = clip(hr_resid / 7, −6, 6)                 z_hrv = clip(hrv_resid / 0.40, −3, 3)
u(z)  = z · (1 − z / 2.8)    z ≥ 0                u(z) = 0.6 · z    z &lt; 0         u ∈ [−2, 0.7]
y = H·D + ε      H = 3.0（心率）/ 2.0（HRV）      R = 1.0 / w_act（心率）/ 2.0 / w_act（HRV）
S = H²P + R      K = P·H / S      |y − H·D| / √S &gt; 4 → 丢弃      D ← D + K·(y − H·D)      P ← (1 − K·H)·P</pre>
{figure(3, "倒 U 变换 u(z_hr)", fig_inverted_u(), "静坐时心率比本人常态高 5–10 bpm 通常是兴奋 / 专注 / 咖啡，高 20 bpm 以上而没在动则是紧张 / 过载。标志位只是文案提示，效应本身走观测通道。")}
    <h3>观测权重：为什么用活动量 A 而不是负荷 L</h3>
    <p>负荷 L 含心率项：心率一高 L 就高，用 L 排除观测会恰好把我们想测的"心率高但没动"排除掉。活动量 A 来自步数，是外生的，所以 A 是主判据：</p>
    <pre>w_act = clip( 1 − (A − 0.10) / (0.35 − 0.10), 0, 1 )
L &gt; 0.60 且 A ≥ 0.10 → w_act = 0（在运动，不用）      L &gt; 0.60 且 A &lt; 0.10 → w_act × 0.5（保留但少信）</pre>
    <p class="small">示例：09:30–10:15 心率 +10 bpm 未运动，D 由 {ar0['d']:+.3f} 升至 {ar1['d']:+.3f}；测试里 +22 bpm 持续 1 h 使 D 下降 &gt; 0.01。</p>
""")
    s5 = section("events", "05", "事件：只认生理效应，不猜原因", f"""
    <table>
      <thead><tr><th>事件</th><th>判据</th><th>效应</th><th>示例实测（第 4 天）</th></tr></thead>
      <tbody>
        <tr><td><span class="dot" style="background:#2c7a4b"></span>运动回弹</td><td>L ≥ 0.5 <b>且</b> A ≥ 0.35 进入；结束后连续 2 点（L ≤ 0.10 或 A ≤ 0.10）；30 min 未平复则放弃</td><td>+0.12 × min(1, 时长/30 min)，≥ 10 min 起算</td><td>18:00–18:30 跑步，18:35 D {ex_before['d']:+.3f} → <b>{ex_after['d']:+.3f}</b>，精力 {fmt(ex_before['vitality'])} → {fmt(ex_after['vitality'])}</td></tr>
        <tr><td><span class="dot" style="background:#3f5bd9"></span>小睡（有分期）</td><td>睡眠段 15 min–3 h，且在 9–20 点</td><td>+0.15 × min(1, 时长/20 min)</td><td>13:00–13:20 午睡，13:30 D {nap_before['d']:+.3f} → <b>{nap_after['d']:+.3f}</b>，精力 {fmt(nap_before['vitality'])} → {fmt(nap_after['vitality'])}</td></tr>
        <tr><td><span class="dot" style="background:#3f5bd9"></span>小睡（无分期）</td><td>z_hr ≤ −0.75 且 A &lt; 0.05 持续 ≥ 15 min 后恢复</td><td>同上</td><td>—</td></tr>
        <tr><td><span class="dot" style="background:#a83f63"></span>疲劳拖拽</td><td>z_hrv ≤ −1 且 w_act &gt; 0 持续 ≥ 30 min</td><td>持续项 dD = −0.03·(−z_hrv)·dt</td><td>—</td></tr>
        <tr><td><span class="dot" style="background:#0e7a83"></span>唤醒 / 过载标志</td><td>0.7 ≤ z_hr ≤ 2.1 / z_hr &gt; 2.1 持续 10 min</td><td>只打标志，效应走观测</td><td>09:30–10:15 AROUSAL</td></tr>
        <tr><td><span class="dot" style="background:#6d4fc4"></span>用户标记</td><td><code>vitality_ctx_event</code></td><td>咖啡 +0.10（20 min 一阶起效）、冥想 +0.08、压力事件 −0.10、酒精 −0.08</td><td>15:00 咖啡标记，D {caf0['d']:+.3f} → 15:15 <b>{caf1['d']:+.3f}</b></td></tr>
      </tbody>
    </table>
    <div class="note warn"><b>关于咖啡：</b>示例里咖啡标记只抬了约 0.04 而不是 0.10——合成数据里心率没有相应升高，观测把它拉回来了。真实用户喝咖啡心率通常升 3–8 bpm，此时观测与标记同向。库不从数据里"识别"咖啡。</div>
""")
    s6 = section("morning", "06", "晨起初始化：四项都相对本人常态", f"""
    <pre>x1 = clip( (ΔE − ΔE̅) / 0.30, ±0.5 )                          <span class="c">今晚能量回充量 vs 本人常态</span>
x2 = clip( (HRV_wake − HRV̅_wake) / 0.44, ±0.5 )               <span class="c">醒前 1 h 的 HRV vs 本人夜间常态</span>
x3 = clip( 0.5·q_时长 + 0.25·q_深睡 − 0.25·q_碎片, ±0.5 )       <span class="c">时长差 2 h / 深睡与碎片化差 10 个百分点记满分</span>
x4 = clip( 0.5·q_血氧均值 − 0.5·q_瞬时下跌, ±0.5 )              <span class="c">夜间血氧平稳度（§7）；没有血氧时为 0</span>
D_init = 0.40 · clip( 0.4·x1 + 0.3·x2 + 0.3·x3 + 0.2·x4, ±0.5 )  <span class="c">|D_init| ≤ 0.2；没有常态的项记 0</span></pre>
    <table>
      <thead><tr><th>规则</th><th>说明</th></tr></thead>
      <tbody>
        <tr><td>主睡眠判定</td><td>sleep_prob ≥ 0.5 为睡眠段；连续 3 个清醒样本（15 min）结束；≥ 3 h 才算"一晚"，否则按小睡</td></tr>
        <tr><td>常态学习</td><td>前 3 晚累计平均（LEARNING），之后 EWMA λ = 0.2（约 5 晚记忆）；<code>get/set_habit</code> 持久化</td></tr>
        <tr><td>续睡合并</td><td>主睡眠结束后 2 h 内再次入睡 → 回滚上次常态更新，整夜重算（夜醒 20 min 不会拆成两晚，<code>test_night_resume_merge</code>）</td></tr>
        <tr><td>不知道就不猜</td><td>第一晚没有常态，D_init = 0</td></tr>
      </tbody>
    </table>
    <p class="small">示例实测：第 4 天 07:10 起床 D_init = {wake3['d_init']:+.3f}；第 5 天只睡 5 h，05:10 起床 D_init = <b>{wake4['d_init']:+.3f}</b>。幅度没有到 −0.2 是因为 x1 项部分抵消：熬到 0 点才睡的人入睡时 E 更低，5 h 回充量并不比平时少多少；真实数据若证实，可改用"醒时 E 相对常态"或提高 w_sleep。</p>
""")
    s6b = section("spo2", "07", "夜间血氧：只进晨起初始化，不进日间观测", f"""
    <p>戒指在指端做反射式 PPG，灌注好、夜间手不动，血氧的可用性明显高于腕戴设备。但"测得准"不等于"该用"——先看生理上它和精力有没有关系：</p>
    <table>
      <thead><tr><th></th><th>日间静息 SpO₂</th><th>夜间 SpO₂</th></tr></thead>
      <tbody>
        <tr><td>和精力的生理关系</td><td>几乎没有。疲劳、压力、咖啡改变的是交感 / 副交感平衡与心率，不改变动脉血氧；健康人静息 SpO₂ 稳定在 96–99%，氧合曲线在这一段是平台，PaO₂ 变化几个 mmHg 都看不出来</td><td>有。夜间通气不稳 → 短暂缺氧 → 化学感受器驱动的交感激活、心率突升、微觉醒。微觉醒常短于分期算法的分辨率，分期仍标"睡着"，但恢复质量已经打折。饮酒、鼻塞、仰睡、高海拔、感冒都先体现在这里</td></tr>
        <tr><td>信噪比</td><td>差：变化幅度（&lt; 1%）小于传感器噪声（约 ±2%），手部活动伪迹大</td><td>好：手不动、灌注稳定；只看相对本夜中位数的下跌，与绝对精度无关</td></tr>
        <tr><td>与 HR / HRV 的关系</td><td>—</td><td>HR / HRV 看到的是交感激活这个<b>结果</b>；血氧看到的是通气不稳这个<b>原因</b>。两者相关但不冗余：心率突升可被分期 / 翻身解释掉，血氧下跌不会</td></tr>
        <tr><td>结论</td><td><b>不用</b></td><td><b>用</b>，作为晨起初始化第四项，权重 0.2</td></tr>
      </tbody>
    </table>
    <h3>本夜统计（常数内存）</h3>
    <pre>每个睡眠样本：low = min(spo2, spo2_min)  → 1% 分箱直方图 [70, 100]；mean 累加
整夜：median = 直方图中位数     dip = 占比( low ≤ median − 3 )     mean = 均值
覆盖：有效样本 ≥ 50% 睡眠样本 且 ≥ 1 h，否则本夜不用血氧（x4 = 0，不打 SPO2 标志，常态不更新）</pre>
    <h3>相对本人常态</h3>
    <pre>q_均值 = clip( (mean − mean̅) / 1.5, ±1 )     高于常态的正向收益减半（SpO₂ 有天花板，"比平时更高"信息量小）
q_下跌 = clip( (dip − dip̅) / 0.10, ±1 )       下跌占比比常态多 10 个百分点记满分
x4 = clip( 0.5·q_均值 − 0.5·q_下跌, ±0.5 )     x4 ≤ −0.25 时晨起点打 BREATH 标志（"昨夜呼吸不如平时平稳"）</pre>
    <table>
      <thead><tr><th>规则</th><th>说明</th></tr></thead>
      <tbody>
        <tr><td>只看夜间</td><td>日间样本的 spo2 字段被完全忽略；不做任何日间观测更新</td></tr>
        <tr><td>只跟自己比</td><td>均值与下跌占比都减去本人常态（前 3 个有血氧的夜累计平均，之后 EWMA λ = 0.2）；高原常住者的常态就是 92%，不会被当成"差"</td></tr>
        <tr><td>不猜原因</td><td>库只说"昨夜血氧比你平时不稳"，不区分饮酒 / 鼻塞 / 海拔 / 任何健康状况，也不设任何绝对阈值（如 90%）</td></tr>
        <tr><td>可退化</td><td>样本不传血氧（NAN / 0）、或 w_spo2 = 0 时，所有输出与不带血氧时逐位一致（<code>test_spo2</code>）；旧画像 <code>vitality_habit_t</code>（memset 0）可直接恢复</td></tr>
        <tr><td>续睡合并</td><td>夜醒后 2 h 内续睡，血氧直方图随分期统计一起保留，按整夜重算</td></tr>
      </tbody>
    </table>
    <p class="small">示例实测：第 6 天睡够 8 h，但前一夜血氧均值 {fmt(wake5['spo2_night'])}%（平时 97%）、瞬时下跌占比 {fmt(100*wake5['spo2_dip'],0)}%（平时 0%）→ x4 = −0.5，D_init = <b>{wake5['d_init']:+.3f}</b>（同样的夜不带血氧时 D_init = 0.000），起床精力 {fmt(wake5['vitality'])}，均衡线 {fmt(100*wake5['v_star'])}。幅度上限是 0.04（0.40 × 0.2 × 0.5）：血氧被刻意放在"修正项"的位置，不会单独决定一天的起点。</p>
    <div class="note warn"><b>为什么不把血氧做成日间观测通道？</b>一旦加进去，卡尔曼更新会把传感器噪声当成状态变化：日间 ±2% 的抖动折算成 z 后与心率通道同量级，而它不携带任何精力信息，只会稀释真正的观测。用了不会更准，只会更吵。</div>
""")
    s7 = section("scale", "08", "两种 0–100 刻度，都只跟自己比", f"""
{figure(4, "vitality 与 percentile", fig_scales(), "score_mode = 0 取 vitality，1 取 percentile（学习期自动退回 vitality），2 取两者平均。")}
    <pre>d_sd = sqrt( P + (V*/E)² · e_sd² )             lower / upper = 100 · clip( V ∓ 1.96 · d_sd )</pre>
    <p class="small">能量层的不确定性按 V 对 E 的灵敏度折进精力带宽；冷启动时带宽大，随能量层基线建立而收窄。</p>
""")
    s8 = section("api", "09", "接口", f"""
    <pre><span class="c">/* 流式：内部自带 energy_ctx_t */</span>
vitality_ctx_t vc;
vitality_ctx_init(&amp;vc, &amp;ecfg, &amp;vcfg);                 <span class="c">/* 配置可为 NULL */</span>
vitality_ctx_set_habit(&amp;vc, &amp;saved_habit);            <span class="c">/* 恢复个人常态；首次使用可不调用 */</span>
for (each 5-min sample) {{
    energy_point_t ep; vitality_point_t vp;
    vitality_ctx_step(&amp;vc, &amp;sample, &amp;ep, &amp;vp);       <span class="c">/* ep.energy, vp.score, vp.quadrant, vp.flags */</span>
}}
vitality_ctx_event(&amp;vc, VITALITY_EVENT_CAFFEINE, 0.0);  <span class="c">/* 用户标记，amplitude ≤ 0 用默认幅度 */</span>
vitality_ctx_get_habit(&amp;vc, &amp;saved_habit);

<span class="c">/* 已有自己的 energy_ctx_t：拿到 energy_point_t 后 */</span>
vitality_ctx_feed(&amp;vc, &amp;sample, &amp;ep, &amp;vp);

<span class="c">/* 批量（内含 energy_compute）*/</span>
vitality_output_t out;
vitality_compute(samples, n, &amp;ecfg, &amp;vcfg, &amp;out);     <span class="c">/* out.points / out.daily / out.energy / out.habit */</span>
vitality_output_free(&amp;out);

<span class="c">/* 命令行：能量列之后追加精力列 */</span>
./build/energy_cli panel.csv --vitality --out out.csv</pre>
    <table>
      <thead><tr><th>标志位</th><th>含义</th></tr></thead>
      <tbody>
        <tr><td><code>SLEEP</code> / <code>GAP</code> / <code>NO_OBS</code></td><td>睡眠中 / 数据缺口 / 本点无可用观测（含运动中、能量基线未建立）</td></tr>
        <tr><td><code>WAKE</code> / <code>LEARNING</code></td><td>本点做了晨起初始化（<code>d_init</code> 字段有效）/ 前 3 晚学习期</td></tr>
        <tr><td><code>EXERCISE</code> / <code>NAP</code> / <code>FATIGUE</code> / <code>USER</code></td><td>本点注入了对应脉冲 / 拖拽</td></tr>
        <tr><td><code>AROUSAL</code> / <code>STRESS</code></td><td>轻度唤醒 / 过载（文案用）</td></tr>
        <tr><td><code>PCT_READY</code></td><td>分位刻度可用（≥ 3 个日间）</td></tr>
        <tr><td><code>SPO2</code> / <code>BREATH</code></td><td>晨起点：本夜血氧覆盖充分、已参与初始化（<code>spo2_night / spo2_dip</code> 字段有效）/ 昨夜血氧平稳度明显低于本人常态</td></tr>
      </tbody>
    </table>
    <p class="small">样本中的 <code class="inl">spo2</code>（窗口均值 %）与 <code class="inl">spo2_min</code>（窗口最低值 %）只在夜间参与计算，没有读数时填 NAN；能量层不读它们。CSV 对应列名 <code class="inl">spo2, spo2_min</code>。</p>
    <p class="small">非法配置：<code class="inl">vitality_config_validate</code> 逐项范围检查，<code class="inl">vitality_compute</code> 返回 <code class="inl">ENERGY_EINVAL</code>，流式上下文每步返回错误；时间倒退返回 <code class="inl">ENERGY_EDOMAIN</code>；NaN 输入按缺失处理（<code class="inl">test_nan_fuzz</code>）。</p>
""")
    s9 = section("verify", "10", "验证：数字有出处", f"""
    <div class="stats">
      <div class="stat g"><div class="n">{N['vitality_checks']}</div><div class="l">精力层测试断言（能量层另有 {N['energy_checks']}）</div></div>
      <div class="stat t"><div class="n">{N['evobeing_checks']}</div><div class="l">主接口层测试断言（evobeing_zhenghuimin）</div></div>
      <div class="stat i"><div class="n">{N['evobeing_ctx_bytes']} B</div><div class="l">evobeing_ctx_t（能量 + 精力上下文）</div></div>
      <div class="stat a"><div class="n">{N['ctx_bytes']} B</div><div class="l">vitality_ctx_t（含 energy_ctx_t {N['ectx_bytes']} B）</div></div>
    </div>
    <table>
      <thead><tr><th>命令</th><th>结果</th></tr></thead>
      <tbody>
        <tr><td><code>make test</code></td><td>能量 OK: {N['energy_checks']} checks, 0 failures；精力 OK: {N['vitality_checks']} checks, 0 failures；主接口 OK: {N['evobeing_checks']} checks, 0 failures</td></tr>
      </tbody>
    </table>
    <p class="small">测试覆盖：默认参数与非法参数；子函数（τ、均衡线、倒 U、裁剪、象限）；有界与睡眠中无观测；流式 == 批量；事件位置 / 方向 / 幅度与日汇总；用户标记与上限；晨起初始化（短睡为负、首晚为 0、学习期退出、常态持久化）；续睡合并；三种 score_mode；体质不变性；开环单调；NaN 与 3 h 缺口；能量层残差暴露；夜间血氧（恒定血氧与无血氧逐位一致、坏夜 D_init −0.04 且只打一次 BREATH、覆盖 30% 的夜被忽略、仅 spo2_min 下跌同样生效、非法值视为缺失、流式 == 批量、画像往返、w_spo2 = 0 退化）。</p>
    <div class="note hot"><b>全部数字来自合成数据。</b>幅度参数（a_nap 0.15、a_exercise 0.12、z_c 2.8 等）需要在真实佩戴数据上用即时自评（"此刻精神如何 1–5"）校准；文档中的曲线是模型行为的展示，不是对真实用户的验证结果。</div>
""")
    s10 = section("limits", "11", "边界", """
    <div class="cols">
      <div class="card"><h4><span class="dot" style="background:#a83f63"></span>没做的</h4>
        <ul><li>识别咖啡 / 冥想 / 压力事件的原因（只接受用户标记）</li><li>跨用户的任何统计或排名</li><li>离线参数拟合（不在设备端做）</li><li>医学用途</li></ul>
      </div>
      <div class="card"><h4><span class="dot" style="background:#a8650f"></span>需要真实数据决定的</h4>
        <ul><li>短睡场景下 D_init 的 ΔE 项是否该换成"醒时 E 相对常态"</li><li>清醒后首小时的残差是否要降权（能量层"心率随 E"假设的伴随现象）</li><li>倒 U 零点 z_c 与各事件幅度</li><li>分位窗口 28 天是否合适（轮班 / 出差）</li><li>血氧项权重 0.2 与下跌判据 3 个百分点（取决于设备的 5 min 聚合方式；给得出 spo2_min 时更灵敏）</li></ul>
      </div>
    </div>
""")
    cover = f"""
  <div class="cover">
    <h3>这份文档回答四个问题</h3>
    <p>为什么能量之外还要一个精力？它怎么从能量层的残差算出来？怎么接入？凭什么相信这些数字？</p>
    <div class="cover-grid">
      <div class="cover-cell"><div class="n" style="color:#0e7a83">1</div><div class="l">个新增隐状态 D</div></div>
      <div class="cover-cell"><div class="n" style="color:#bf4520">0</div><div class="l">处能量层改动（只加两列输出）</div></div>
      <div class="cover-cell"><div class="n" style="color:#2c7a4b">{N['vitality_checks']}</div><div class="l">精力层测试断言</div></div>
      <div class="cover-cell"><div class="n" style="color:#a8650f">{N['evobeing_ctx_bytes']} B</div><div class="l">能量 + 精力上下文内存</div></div>
    </div>
  </div>
"""
    toc = [("why", "1 为什么"), ("arch", "2 级联结构"), ("step", "3 一步之内"), ("obs", "4 观测"), ("events", "5 事件"),
           ("morning", "6 晨起初始化"), ("spo2", "7 夜间血氧"), ("scale", "8 刻度"), ("api", "9 接口"), ("verify", "10 验证"), ("limits", "11 边界")]
    return page("evobeing 技术说明 · 精力", read_style("technical.html"),
                f"evobeing · zhenghuimin · 精力 · 技术说明 v{N['version']}",
                "在能量之上，再算一条精力",
                "能量回答“储备还剩多少”，精力回答“此刻用得出来多少”。精力层把能量层解释不了的心率 / HRV 残差、运动与小睡的回弹、睡眠质量与夜间血氧平稳度相对本人常态的偏差，建模成一个会自己消退的偏离量 D，叠加在“能量 × 昼夜节律”的均衡线上。",
                ["<b>C99</b> + libm", "级联于 energy.c", "<b>1</b> 个标量隐状态", "批量 <b>+</b> 流式"],
                toc, s1 + s2 + s3 + s4 + s5 + s6 + s6b + s7 + s8 + s9 + s10,
                f"evobeing v{N['version']} · C99 · 无第三方依赖",
                "技术说明 · 图表基于真实实现与实测（合成）数据", cover)


# ---------------------------------------------------------------- 产品说明
def product(rows):
    N = NUMBERS
    nap_before, nap_after = pick(rows, 3, 13.25), pick(rows, 3, 13.5)
    ex_before, ex_after = pick(rows, 3, 18.5), pick(rows, 3, 18.5833)
    ann3 = [(7.1667, "起床", "#71768a"), (10.0, "上午高点", "#0e7a83"), (13.5, "午睡后", "#3f5bd9"),
            (15.25, "咖啡标记", "#6d4fc4"), (18.5833, "跑步后", "#2c7a4b")]
    ann4 = [(5.1667, "只睡 5 小时，起床", "#71768a"), (14.5, "午后低谷", "#bf4520")]
    w4 = next(r for r in rows if r["flags"] & FLAG["WAKE"] and 96 <= r["t_h"] < 120)
    w5 = next(r for r in rows if r["flags"] & FLAG["WAKE"] and 120 <= r["t_h"] < 144)
    ann5 = [(7.1667, "起床（昨夜血氧不稳）", "#a83f63"), (10.0, "上午", "#0e7a83")]
    p3, p4 = pick(rows, 3, 14.5), pick(rows, 4, 14.5)
    s1 = section("what", "01", "它是什么，不是什么", """
    <div class="cols">
      <div class="card"><h4><span class="dot" style="background:#0e7a83"></span>它是</h4>
        <ul><li>一条和能量同图的第二条曲线：<b>此刻用得出来多少</b></li><li>由你自己的心率、HRV、活动、睡眠、夜间血氧相对<b>你自己的常态</b>推算</li><li>会对午睡、运动后的回弹、下午的提神、连续紧张的下压作出反应，然后自己消退</li><li>一个 0–100 的数，外加一句四象限文案</li></ul>
      </div>
      <div class="card"><h4><span class="dot" style="background:#bf4520"></span>它不是</h4>
        <ul><li>不是和别人比的排名——库里没有任何跨人数据</li><li>不是情绪、注意力或任何健康状况的判断</li><li>不知道你喝了什么、做了什么：咖啡 / 冥想 / 压力事件需要你手动标记，不标记就只通过心率间接体现</li><li>没有数据时不会编：睡觉、缺口、运动中都回到均衡线</li></ul>
      </div>
    </div>
""")
    s2 = section("day", "02", "一天怎么走：能量与精力并排看", f"""
    <p>下面两张图是示例程序 <code class="inl">make example</code> 的真实输出（合成数据）。赤陶色是能量，青色是精力，灰色虚线是精力的均衡线——没有任何事情发生时精力应该在的位置。</p>
{figure(1, "事件日：唤醒、午睡、咖啡标记、跑步", fig_dual_track(rows, 3, "第 4 天（事件日）", ann3), f"午睡 20 分钟后精力由 {fmt(nap_before['vitality'])} 回到 {fmt(nap_after['vitality'])}，随后约 1 小时衰回均衡线；傍晚跑步期间不用心率做观测，结束平复后回弹 {fmt(ex_before['vitality'])} → {fmt(ex_after['vitality'])}。能量这条线对这些都几乎没有反应——这正是分成两条线的理由。")}
{figure(2, "短睡日：只睡 5 小时", fig_dual_track(rows, 4, "第 5 天（只睡 5 h）", ann4), f"睡得少，能量起点低，精力的均衡线也随之整体下移。晨起 {int(w4['tod']):02d}:{int(round((w4['tod']%1)*60)):02d} 精力 {fmt(w4['vitality'])}，均衡线 {fmt(100*w4['v_star'])}，起床偏置 D_init = {w4['d_init']:+.3f}——略低于均衡线但幅度不大；14:30 能量 {fmt(p4['energy'])}（前一日同刻 {fmt(p3['energy'])}）、精力 {fmt(p4['vitality'])}（前一日 {fmt(p3['vitality'])}）。精力这天整体偏低，主要是均衡线跟着能量一起下移；起床偏置本身只占其中一小部分，对应的体感是起床时比平时更没劲一点。")}
{figure(3, "血氧不稳的夜：睡够了，但起床偏低一点", fig_dual_track(rows, 5, "第 6 天（夜间血氧不稳）", ann5), f"这一夜睡满 8 小时，分期和平时一样，但戒指记录的夜间血氧均值 {fmt(w5['spo2_night'])}%（平时 97%），有 {fmt(100*w5['spo2_dip'],0)}% 的时段出现瞬时下跌（平时没有）。晨起偏置 D_init = {w5['d_init']:+.3f}，起床精力 {fmt(w5['vitality'])}（均衡线 {fmt(100*w5['v_star'])}）；同样的夜若不带血氧，偏置为 0。这是一个修正项：最多只影响起床起点约 4 分，之后照常由白天的数据接管。")}
""")
    s3 = section("quad", "03", "两条线怎么一起读：四象限", f"""
{figure(4, "能量 × 精力", fig_quadrants(), "阈值默认 50 / 50，可按产品需要调整。文案只描述状态与建议，不下任何判断。")}
    <table>
      <thead><tr><th>用户的一天</th><th>能量</th><th>精力</th><th>可以说的话</th></tr></thead>
      <tbody>
        <tr><td>睡够了，上午专注工作</td><td>高</td><td>高</td><td>"状态极佳，适合安排要紧的事。"</td></tr>
        <tr><td>午饭后犯困</td><td>还高</td><td>低</td><td>"储备还在，只是这会儿提不起劲——走动一下，或者小睡 15 分钟。"</td></tr>
        <tr><td>连续熬夜后靠咖啡撑着</td><td>低</td><td>暂时高</td><td>"现在精神不错，但储备不足，别透支，今晚早点睡。"</td></tr>
        <tr><td>熬夜 + 下午 3 点</td><td>低</td><td>低</td><td>"两条线都低了，建议休息。"</td></tr>
      </tbody>
    </table>
""")
    s4 = section("scale", "04", "0–100 怎么读", f"""
{figure(5, "两种刻度", fig_scales(), "第一种立即可用；第二种需要先积累 3 个白天，之后文案可以说“比你 80% 的时候好”。")}
    <div class="note"><b>只跟自己比。</b>两种刻度的参照都是这个人自己：能量相对自己的基线，精力相对自己的均衡线，分位相对自己过去约四周的白天。换一个静息心率 55 的人和一个 75 的人，各自的分数口径是一样的——这是模型假设，也是产品口径。</div>
    <table>
      <thead><tr><th>阶段</th><th>精力能显示什么</th></tr></thead>
      <tbody>
        <tr><td>第 1 晚之前</td><td>精力 = 均衡线（能量 × 节律），事件回弹已可用；晨起不做个性化初始化</td></tr>
        <tr><td>第 1–3 晚（学习期）</td><td>晨起初始化开始相对常态；带 LEARNING 标志，文案宜保守</td></tr>
        <tr><td>第 3 个白天之后</td><td>分位刻度可用，可以说"比你 80% 的时候好"</td></tr>
        <tr><td>换手机 / 重装</td><td>把 <code>vitality_habit_t</code> 存回去即可，不用重新学习</td></tr>
      </tbody>
    </table>
""")
    s5 = section("mark", "05", "手动标记：锦上添花，不是必需", """
    <table>
      <thead><tr><th>标记</th><th>默认效应</th><th>怎么解释给用户</th></tr></thead>
      <tbody>
        <tr><td>咖啡 / 茶</td><td>约 20 分钟内逐渐 +10 分，1 小时左右消退</td><td>"记一下，曲线会把这阵提神画出来。"</td></tr>
        <tr><td>冥想 / 放松</td><td>+8 分</td><td>"放松后的回升。"</td></tr>
        <tr><td>压力事件</td><td>−10 分</td><td>"记下来，之后回看更容易对上。"</td></tr>
        <tr><td>酒精</td><td>−8 分</td><td>—</td></tr>
      </tbody>
    </table>
    <p class="small">标记只是把一个脉冲注入模型，之后仍由心率 / HRV 的观测校正：如果标记了咖啡但心率完全没变，曲线会把它拉回来一部分。运动与午睡不需要标记，库自己能从负荷、活动量和分期里看出来。</p>
""")
    s6 = section("say", "06", "能说与不能说", """
    <div class="cols">
      <div class="card"><h4><span class="dot" style="background:#2c7a4b"></span>可以说</h4>
        <ul><li>"此刻的精力大约 68，比你 80% 的时候好"</li><li>"午睡后回升了 13 分"</li><li>"精神不错，但储备不足，别透支"</li><li>"最近三天的起床精力都偏低"</li><li>"昨夜呼吸不如平时平稳，起床起点低了约 4 分"</li></ul>
      </div>
      <div class="card"><h4><span class="dot" style="background:#bf4520"></span>不能说</h4>
        <ul><li>"你比 70% 的同龄人精力好"（没有跨人数据，也不该有）</li><li>"你的注意力 / 情绪 / 健康如何"</li><li>"检测到你喝了咖啡"（不识别原因）</li><li>任何"测量 / 判定 / 预警"口吻</li></ul>
      </div>
    </div>
    <div class="note hot"><b>一句话：</b>它是一条相对你自己的、会自己消退的状态曲线，用来解释"明明储备还在为什么提不起劲"和"明明没睡够为什么现在还挺精神"。文案保持在这个范围内。</div>
""")
    s7 = section("faq", "07", "常见问题", """
    <table>
      <tbody>
        <tr><td><b>精力和能量为什么有时方向相反？</b></td><td>能量是储备，精力是此刻的可用度。午睡、运动后、咖啡都会让精力短暂高于均衡线，而储备没变；连续紧张则相反。</td></tr>
        <tr><td><b>为什么睡觉时精力曲线是平的？</b></td><td>睡觉时不做观测，精力跟随均衡线。醒来那一刻才根据这一夜相对你常态的好坏，给一个起点。</td></tr>
        <tr><td><b>第一天为什么看不到"比你 x% 的时候好"？</b></td><td>分位需要至少 3 个白天的数据；之前显示的是绝对刻度。</td></tr>
        <tr><td><b>我跑步时精力为什么没有涨？</b></td><td>运动中心率不能当作状态信号，这段不观测；平复后才给回弹。</td></tr>
        <tr><td><b>戒指测血氧很准，为什么白天的血氧不参与？</b></td><td>白天安静状态下血氧几乎是一条直线，困不困、累不累都不改变它——它不携带精力信息。夜里不一样：呼吸不平稳的时段会伴随心跳突升和短暂的浅醒，分期可能仍然显示"睡着"，但恢复质量已经打折。所以只用夜间血氧，而且只跟你自己平时的夜比。</td></tr>
        <tr><td><b>"昨夜血氧不如平时平稳"是在说健康问题吗？</b></td><td>不是。它只是一个相对你自己常态的比较，不设任何绝对阈值，也不区分原因——喝了酒、鼻塞、仰睡、换了海拔都会这样。文案保持在"昨晚恢复打了折扣，今天起点低一点"的范围内。</td></tr>
        <tr><td><b>这些数字在真实用户上验证过吗？</b></td><td>还没有。目前所有参数和图都来自合成数据；上线前需要用即时自评校准幅度。</td></tr>
      </tbody>
    </table>
""")
    toc = [("what", "1 是什么"), ("day", "2 一天怎么走"), ("quad", "3 四象限"), ("scale", "4 刻度"), ("mark", "5 手动标记"), ("say", "6 能说与不能说"), ("faq", "7 常见问题")]
    return page("evobeing 产品说明 · 精力", read_style("product.html"),
                "evobeing · zhenghuimin · 精力 · 产品说明",
                "现在，用得出来多少",
                "能量告诉你今天还剩多少余量，精力告诉你这一刻能用出来多少。它只和你自己比，会对午睡、运动、咖啡和连续紧张作出反应，然后自己消退。",
                ["和能量同图", "只跟自己比", "0–100 + 四象限", "手动标记可选"],
                toc, s1 + s2 + s3 + s4 + s5 + s6 + s7,
                f"evobeing v{N['version']} · 精力 · 面向产品 / 设计 / 运营",
                "产品说明 · 曲线来自示例程序的真实输出（合成数据）")


def main():
    csv_path = sys.argv[1] if len(sys.argv) > 1 else os.path.join(ROOT, "build", "vitality_demo.csv")
    rows = load_csv(csv_path)
    with open(os.path.join(DOCS, "vitality_technical.html"), "w", encoding="utf-8") as f:
        f.write(technical(rows))
    with open(os.path.join(DOCS, "vitality_product.html"), "w", encoding="utf-8") as f:
        f.write(product(rows))
    print("docs/vitality_technical.html, docs/vitality_product.html")


if __name__ == "__main__":
    main()
