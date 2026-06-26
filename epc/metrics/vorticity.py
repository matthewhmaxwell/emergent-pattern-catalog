"""Rotational-order (vorticity) lens — coherent CIRCULATION / milling / vortex flow, a
documented under-covered morphology axis (rotation/vortices).

Distinct from the orientation channel (polar/nematic ALIGNMENT — a milling swarm has LOW
polar order, the velocities point every direction around the ring) and from clustering: this
measures net angular momentum about the centre of mass,

    L = | Σ_i (x_i v_yi - y_i v_xi) | / Σ_i |r_i| |v_i|    (r_i relative to the COM)

so |L| ~ 1 for coherent rotation (milling, vortex, spiral flow), ~0 for translation
(flocking), random motion, or a static configuration. Reads per-agent 'velocities' (Nx2),
or derives a velocity from consecutive 'positions' (same-index correspondence). Returns None
for substrates with no agent motion.
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional

import numpy as np


def _pos_vel(history: List[Dict[str, Any]]):
    frames = [f for f in history if isinstance(f, dict)]
    late = frames[max(0, len(frames) // 2 - 1):] or frames        # one extra for differencing
    out = []
    prev = None
    for f in late:
        p = np.asarray(f["positions"], dtype=float)[:, :2] if "positions" in f and \
            np.asarray(f["positions"]).ndim == 2 else None
        if p is None or len(p) < 4:
            prev = p
            continue
        if "velocities" in f:
            v = np.asarray(f["velocities"], dtype=float)
            if v.ndim == 2 and v.shape[0] == p.shape[0] and v.shape[1] >= 2:
                out.append((p, v[:, :2]))
                prev = p
                continue
        if prev is not None and prev.shape == p.shape:
            out.append((p, p - prev))
        prev = p
    return out


def _ang_mom(p: np.ndarray, v: np.ndarray) -> float:
    com = p.mean(0)
    r = p - com
    rxv = r[:, 0] * v[:, 1] - r[:, 1] * v[:, 0]
    denom = float((np.linalg.norm(r, axis=1) * np.linalg.norm(v, axis=1)).sum())
    if denom < 1e-12:
        return 0.0
    return float(abs(rxv.sum()) / denom)


def rotational_order(history: List[Dict[str, Any]]) -> Optional[Dict[str, float]]:
    pairs = _pos_vel(history)
    if len(pairs) < 1:
        return None
    Ls = [_ang_mom(p, v) for p, v in pairs]
    Ls = [x for x in Ls if x == x]
    if not Ls:
        return None
    return {"ang_mom": round(float(np.mean(Ls)), 4)}
