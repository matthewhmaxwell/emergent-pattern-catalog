"""Shared helpers for the repro v2 group A runs (P23, P27, P18).

Spec: docs/validation_rebuild/repro_v2/group_A.md (frozen at commit 39688e3).
"""
from __future__ import annotations

import json
import os
import platform
import socket
import subprocess
import time
from pathlib import Path
from typing import Any

import numpy as np

REPO = Path(__file__).resolve().parents[3]
OUT_DIR = REPO / "analysis" / "outputs" / "repro_v2"
RAW_DIR = OUT_DIR / "raw"
SPEC_DIR = REPO / "docs" / "validation_rebuild" / "repro_v2"
SUPPORT_DIR = SPEC_DIR / "group_A_support"
FROZEN_PLAN_COMMIT = "39688e3"


def git_info() -> dict[str, Any]:
    def _git(*args: str) -> str:
        return subprocess.run(["git", *args], cwd=REPO, capture_output=True,
                              text=True, check=True).stdout.strip()
    status = _git("status", "--porcelain", "--untracked-files=no")
    return {
        "commit": _git("rev-parse", "HEAD"),
        "branch": _git("rev-parse", "--abbrev-ref", "HEAD"),
        "tracked_tree_clean": status == "",
        "frozen_plan_commit": FROZEN_PLAN_COMMIT,
    }


def env_info() -> dict[str, Any]:
    return {
        "host": socket.gethostname(),
        "python": platform.python_version(),
        "numpy": np.__version__,
        "pid": os.getpid(),
    }


def append_jsonl(path: Path, rec: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open(path, "a") as fh:
        fh.write(json.dumps(rec, default=_np_default) + "\n")
        fh.flush()
        os.fsync(fh.fileno())


def read_jsonl(path: Path) -> list[dict[str, Any]]:
    if not path.exists():
        return []
    with open(path) as fh:
        return [json.loads(line) for line in fh if line.strip()]


def write_json(path: Path, obj: Any) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    tmp = path.with_suffix(path.suffix + ".tmp")
    with open(tmp, "w") as fh:
        json.dump(obj, fh, indent=2, default=_np_default)
        fh.write("\n")
    os.replace(tmp, path)


def _np_default(o: Any) -> Any:
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    if isinstance(o, (np.bool_,)):
        return bool(o)
    raise TypeError(f"not JSON serializable: {type(o)}")


DEVIATIONS_FILE = Path(__file__).resolve().parent / "deviations_group_A.json"


def load_deviations(pattern: str) -> list[dict[str, Any]]:
    """Protocol deviations, written (with reasons) before the affected rule is evaluated."""
    if not DEVIATIONS_FILE.exists():
        return []
    return json.loads(DEVIATIONS_FILE.read_text()).get(pattern, [])


def now_iso() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%S%z")
