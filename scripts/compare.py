"""比较 C 输出与 Python 期望输出（EKF / 开环 / 负荷标定 三项）。

用法：
    python scripts/compare.py            # 需要先有 data/c_out.csv / c_open.csv
    python scripts/check_parity.py       # 一键：生成 → 跑 C → 比较
退出码 0 表示在容差内一致。
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
DATA = HERE.parent / "data"
COLS = ["s", "e", "e_sd", "energy", "lower", "upper", "base"]
TOL = 1e-9


def compare_points(c_path: Path, py_path: Path, label: str) -> bool:
    c = pd.read_csv(c_path)
    py = pd.read_csv(py_path)
    if len(c) != len(py):
        print(f"FAIL {label}: 行数不一致 {len(c)} vs {len(py)}")
        return False
    print(f"\n[{label}] {len(c)} 行")
    print(f"{'列':8s} {'max|diff|':>11s} {'mean|diff|':>11s}")
    ok = True
    for col in COLS:
        d = np.abs(c[col].to_numpy(float) - py[col].to_numpy(float))
        mx = float(np.nanmax(d))
        print(f"{col:8s} {mx:11.3e} {float(np.nanmean(d)):11.3e}" + ("" if mx <= TOL else "  <-- 超差"))
        ok = ok and mx <= TOL
    return ok


def parse_stats(text: str) -> dict[str, float]:
    out = {}
    for key in ("hr_rest", "hr_max", "l_ref"):
        m = re.search(rf"{key}=([-0-9.eE+]+)", text)
        if m:
            out[key] = float(m.group(1))
    return out


def compare_stats(c_text: str) -> bool:
    py = parse_stats((DATA / "expected_stats.txt").read_text(encoding="utf-8"))
    c = parse_stats(c_text)
    print("\n[负荷标定]")
    ok = True
    for key in ("hr_rest", "hr_max", "l_ref"):
        a, b = c.get(key, float("nan")), py.get(key, float("nan"))
        rel = abs(a - b) / max(abs(b), 1e-12)
        print(f"{key:8s} C={a:14.8f}  Py={b:14.8f}  rel={rel:.2e}")
        ok = ok and rel <= 1e-9
    return ok


def main() -> int:
    ok = compare_points(DATA / "c_out.csv", DATA / "expected.csv", "EKF")
    ok = compare_points(DATA / "c_open.csv", DATA / "expected_open.csv", "开环 ODE") and ok
    stats_file = DATA / "c_stats.txt"
    if stats_file.is_file():
        ok = compare_stats(stats_file.read_text(encoding="utf-8")) and ok
    print("\n结果:", "PASS" if ok else "FAIL")
    return 0 if ok else 1


if __name__ == "__main__":
    sys.exit(main())
