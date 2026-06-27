"""Anomalous-transport / diffusion-class lens — a documented missing axis. The battery's
heavy-tail lens is on event SIZES (avalanches, wealth); none reads particle TRANSPORT. This
lens reports the diffusion class from trajectories:

  msd_exponent     — beta from time-averaged MSD(Delta) ~ Delta^beta on log-log:
                     beta<1 subdiffusive (caging/CTRW), ~1 normal, 1<beta<2 superdiffusive
                     (Levy/persistent), ~2 ballistic. THE discriminator (the exponent, not the
                     diffusion coefficient / magnitude).
  ergodicity_break — |log(TAMSD/EAMSD)| at matched lags. ~0 ergodic (Brownian); >0 = aging /
                     ergodicity breaking (heavy-tailed-waiting CTRW), where time- and ensemble-
                     averages disagree.

Reads index-aligned 'positions' (Nx2) trajectories over >=8 frames. Returns None otherwise.
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional

import numpy as np


def _trajectories(history: List[Dict[str, Any]]) -> Optional[np.ndarray]:
    seq = []
    for f in history:
        if isinstance(f, dict) and f.get("positions") is not None:
            p = np.asarray(f["positions"], dtype=float)
            if p.ndim == 2 and p.shape[1] == 2:
                seq.append(p)
    if len(seq) < 8:
        return None
    n = min(len(p) for p in seq)
    if n < 4:
        return None
    return np.stack([p[:n] for p in seq])                    # (T, N, 2)


def anomalous_transport(history: List[Dict[str, Any]]) -> Optional[Dict[str, float]]:
    X = _trajectories(history)
    if X is None:
        return None
    T, N, _ = X.shape
    maxlag = max(3, T // 4)
    lags = np.arange(1, maxlag + 1)
    tamsd = np.array([float(((X[l:] - X[:-l]) ** 2).sum(-1).mean()) for l in lags])
    if not np.all(tamsd > 0):
        return None
    eamsd = ((X - X[0:1]) ** 2).sum(-1).mean(1)              # length T
    x = np.log(lags); y = np.log(tamsd)
    if x.std() < 1e-9:
        return None
    beta = float(np.polyfit(x, y, 1)[0])
    ea = eamsd[lags]
    valid = ea > 0
    eb = float(np.mean(np.abs(np.log((tamsd[valid]) / (ea[valid] + 1e-12))))) if valid.any() else 0.0
    return {"msd_exponent": round(beta, 3), "ergodicity_break": round(eb, 3)}
