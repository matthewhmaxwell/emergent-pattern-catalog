"""Shared helpers for the repro_v2 group B reproductions (P5, P6, P10).

Plan: docs/validation_rebuild/repro_v2/REPRO_V2_PLAN.md; specs: group_B.md.
"""
from __future__ import annotations

import json
import os
import subprocess
import sys
from pathlib import Path
from typing import Any

import numpy as np

REPO = Path(__file__).resolve().parents[3]
if str(REPO) not in sys.path:
    sys.path.insert(0, str(REPO))

OUT = REPO / "analysis" / "outputs" / "repro_v2"
RAW = OUT / "raw"
SUPPORT = REPO / "docs" / "validation_rebuild" / "repro_v2" / "group_B_support"


def git_info() -> dict[str, Any]:
    """Commit of the worktree and whether model/driver code differs from it."""
    def _git(*args: str) -> str:
        return subprocess.check_output(["git", *args], cwd=REPO, text=True).strip()
    try:
        head = _git("rev-parse", "HEAD")
        dirty = _git("status", "--porcelain", "--", "epc", "analysis/reproductions/repro_v2")
        return {"commit": head, "code_dirty": bool(dirty)}
    except Exception as exc:  # pragma: no cover
        return {"commit": None, "error": repr(exc)}


def to_jsonable(o: Any) -> Any:
    if isinstance(o, dict):
        return {str(k): to_jsonable(v) for k, v in o.items()}
    if isinstance(o, (list, tuple)):
        return [to_jsonable(v) for v in o]
    if isinstance(o, np.ndarray):
        return to_jsonable(o.tolist())
    if isinstance(o, (np.bool_,)):
        return bool(o)
    if isinstance(o, np.integer):
        return int(o)
    if isinstance(o, np.floating):
        f = float(o)
        return f if np.isfinite(f) else str(f)
    if isinstance(o, float) and not np.isfinite(o):
        return str(o)
    return o


def write_json(path: Path, obj: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with open(tmp, "w") as fh:
        json.dump(to_jsonable(obj), fh, indent=1)
    os.replace(tmp, path)


def read_json(path: Path) -> Any:
    with open(path) as fh:
        return json.load(fh)


def t_interval(x: np.ndarray, level: float = 0.95) -> dict[str, float]:
    from scipy import stats
    x = np.asarray(x, dtype=float)
    n = len(x)
    mean = float(np.mean(x))
    sd = float(np.std(x, ddof=1)) if n > 1 else float("nan")
    se = sd / np.sqrt(n) if n > 1 else float("nan")
    tq = float(stats.t.ppf(0.5 + level / 2, n - 1)) if n > 1 else float("nan")
    return {"n": n, "mean": mean, "sd": sd, "se": se,
            "ci95_lo": mean - tq * se, "ci95_hi": mean + tq * se}
