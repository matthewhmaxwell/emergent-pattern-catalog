"""Extreme-value clustering lens — the extremal index (a documented missing axis). The covered
heavy-tail lens measures the BULK tail shape (how heavy); it says nothing about whether extremes
arrive INDEPENDENTLY or in CLUSTERS in time. Two series with identical tails can differ entirely:
iid heavy-tailed exceedances (theta~1) vs a persistent process whose extremes come in bursts
(theta<1, mean cluster size 1/theta).

  extremal_index — theta via the runs declustering estimator: n_clusters / n_exceedances over a
                   high threshold. ~1 = extremes independent (Poisson); <1 = extremes CLUSTER.
                   THE discriminator (temporal clustering of extremes, orthogonal to tail weight).

Reads an explicit 'scalar'/'order_parameter' series, else derives one (field mean / Kuramoto R /
radius of gyration). Needs enough exceedances; returns None otherwise.
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional

import numpy as np


def _series(history: List[Dict[str, Any]]) -> Optional[np.ndarray]:
    s = []
    for f in history:
        if not isinstance(f, dict):
            continue
        if "scalar" in f:
            s.append(float(f["scalar"])); continue
        if "order_parameter" in f:
            s.append(float(f["order_parameter"])); continue
        if "phases" in f:
            th = np.asarray(f["phases"], dtype=float).ravel()
            s.append(float(abs(np.exp(1j * th).mean()))); continue
        done = False
        for k in ("field", "u", "concentration"):
            if k in f:
                s.append(float(np.asarray(f[k], dtype=float).mean())); done = True; break
        if done:
            continue
        if f.get("positions") is not None:
            p = np.asarray(f["positions"], dtype=float)
            if p.ndim == 2:
                s.append(float(np.hypot(*(p - p.mean(0)).T).mean()))
    s = np.asarray(s, dtype=float)
    return s if s.size >= 80 and np.ptp(s) > 0 else None


def extreme_value(history: List[Dict[str, Any]], q: float = 0.9, run_gap: int = 2) -> Optional[Dict[str, float]]:
    x = _series(history)
    if x is None:
        return None
    u = np.quantile(x, q)
    idx = np.where(x > u)[0]
    if idx.size < 10:
        return None
    gaps = np.diff(idx)
    n_clusters = 1 + int((gaps > run_gap).sum())             # new cluster when gap exceeds run_gap
    theta = float(n_clusters / idx.size)
    return {"extremal_index": round(min(1.0, theta), 3)}
