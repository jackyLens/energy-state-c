"""生成 C/Python 对照数据：3 天 5 分钟面板 + Python 端期望输出。

用法：
    python scripts/gen_panel.py [days]

产物：
    data/panel.csv     C 端输入（含 Python 算好的 load/activity）
    data/expected.csv  Python 端输出（s/e/e_sd/energy/lower/upper/base）
"""

from __future__ import annotations

import os
import sys
import warnings
from pathlib import Path

import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
# 需要 Python 参考实现（energy_api.py + energy_state 包）。可用环境变量指定：
#   ENERGY_PY_ROOT   包含 energy_api.py 的目录（默认 ../../）
#   ENERGY_STATE_SRC 包含 energy_state/ 的 src 目录（默认 $ENERGY_PY_ROOT/repo/energy/src）
ROOT = Path(os.environ.get("ENERGY_PY_ROOT", HERE.parents[1])).expanduser().resolve()
STATE_SRC = Path(os.environ.get("ENERGY_STATE_SRC", ROOT / "repo" / "energy" / "src")).expanduser()
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(STATE_SRC))

from energy_api import compute_energy  # noqa: E402
from energy_state.data.synthetic import generate  # noqa: E402
from energy_state.load import compute_load  # noqa: E402
from energy_state.observation import default_model  # noqa: E402

DATA = HERE.parent / "data"


def main(days: int = 3) -> None:
    DATA.mkdir(parents=True, exist_ok=True)
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", RuntimeWarning)
        ds = generate(days=days, seed=7)
    panel = ds.panel

    frame = pd.DataFrame(
        {
            "t_h": panel.t_h,
            "tod": panel.tod,
            "wake": panel.stage[:, 0],
            "light": panel.stage[:, 1],
            "deep": panel.stage[:, 2],
            "rem": panel.stage[:, 3],
            "hr": panel.hr,
            "hrv_log": panel.hrv_log,
            "steps": panel.steps_5m,
            "kcal": np.zeros(panel.n),
            "stress": panel.stress,
            "q_hr": panel.q_hr,
            "q_hrv": panel.q_hrv,
            "load": panel.load,
            "activity": panel.activity,
        }
    )
    frame.to_csv(DATA / "panel.csv", index=False, float_format="%.17g")

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res = compute_energy(
            panel.t_h,
            start=panel.meta["start_iso"],
            stage_probs=panel.stage,
            hr=panel.hr,
            hrv_log=panel.hrv_log,
            steps_5m=panel.steps_5m,
            stress=panel.stress,
            q_hr=panel.q_hr,
            q_hrv=panel.q_hrv,
            load=panel.load,
            activity=panel.activity,
            obs=default_model(),
            freeze_anomaly_nights=False,
        )

    cols = ["t_h", "s", "e", "e_sd", "energy", "lower", "upper", "base"]
    res.samples[cols].to_csv(DATA / "expected.csv", index=False, float_format="%.17g")

    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        res_open = compute_energy(
            panel.t_h,
            start=panel.meta["start_iso"],
            stage_probs=panel.stage,
            hr=panel.hr,
            hrv_log=panel.hrv_log,
            steps_5m=panel.steps_5m,
            stress=panel.stress,
            q_hr=panel.q_hr,
            q_hrv=panel.q_hrv,
            load=panel.load,
            activity=panel.activity,
            obs=default_model(),
            use_ekf=False,
            freeze_anomaly_nights=False,
        )
    res_open.samples[cols].to_csv(DATA / "expected_open.csv", index=False, float_format="%.17g")

    lres = compute_load(panel)
    (DATA / "expected_stats.txt").write_text(
        f"hr_rest={lres.hr_rest:.17g}\nhr_max={lres.hr_max:.17g}\nl_ref={lres.l_ref:.17g}\n",
        encoding="utf-8",
    )

    print(f"panel: {panel.n} 行 → {DATA / 'panel.csv'}")
    print(f"expected: {len(res.samples)} 行（EKF） / {len(res_open.samples)} 行（开环）")
    print(f"load stats: hr_rest={lres.hr_rest:.3f} hr_max={lres.hr_max:.3f} l_ref={lres.l_ref:.6f}")
    print(f"python energy: {res.energy.min():.2f}–{res.energy.max():.2f} (mean {res.energy.mean():.2f})")


if __name__ == "__main__":
    main(int(sys.argv[1]) if len(sys.argv) > 1 else 3)
