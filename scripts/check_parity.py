"""一键对照：生成面板 → 跑 C → 比较 C/Python。

用法：
    python scripts/check_parity.py [days]
"""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
DATA = ROOT / "data"
CLI = ROOT / "build" / "energy_cli"


def main(days: int = 3) -> int:
    if not CLI.is_file():
        print(f"缺少 {CLI}，先执行 make")
        return 1
    subprocess.run([sys.executable, str(HERE / "gen_panel.py"), str(days)], check=True)

    subprocess.run([str(CLI), str(DATA / "panel.csv"), "--out", str(DATA / "c_out.csv")],
                   check=True, stderr=subprocess.DEVNULL)
    subprocess.run([str(CLI), str(DATA / "panel.csv"), "--out", str(DATA / "c_open.csv"), "--no-ekf"],
                   check=True, stderr=subprocess.DEVNULL)
    stats = subprocess.run([str(CLI), str(DATA / "panel.csv"), "--load-stats"],
                           check=True, capture_output=True, text=True)
    (DATA / "c_stats.txt").write_text(stats.stdout, encoding="utf-8")

    return subprocess.run([sys.executable, str(HERE / "compare.py")], check=False).returncode


if __name__ == "__main__":
    sys.exit(main(int(sys.argv[1]) if len(sys.argv) > 1 else 3))
