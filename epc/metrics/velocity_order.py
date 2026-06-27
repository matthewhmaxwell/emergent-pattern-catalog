"""Collective velocity-field order from BARE POSITIONS (a documented missing axis). On an
Nx2-only stream the battery never forms velocities, so a milling vortex and a polar flock look
identical to S(k). This lens finite-differences positions into velocities and reports the
collective-motion triple — but its NOVEL coordinate (not owned by vorticity = net angular
momentum, nor by polar order = mean direction) is the CONNECTED VELOCITY-FLUCTUATION
correlation length xi (Cavagna et al.): C(r) = <dv_i . dv_j> over pairs at separation r, with
dv = v - <v>. xi = first zero-crossing of C(r) / system size.

Why xi is the discriminator: a uniformly DRIFTING gas has high polarization yet independent
fluctuations -> xi ~ 0. A real flock has spatially correlated fluctuations -> xi large. So xi
separates genuine collective order from mere drift, which polarization cannot.

Scalars (late-window mean):
  velcorr_length — xi / L, connected velocity-fluctuation correlation length (THE discriminator).
  polarization   — |mean unit velocity| (drift/flock alignment; fooled by drift alone).
  mill_strength  — |angular momentum about COM| / (|r||v|) (rotation; overlaps vorticity).
Reads 'positions' (Nx2) across >=3 frames. Returns None otherwise.
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional

import numpy as np


def _positions(history: List[Dict[str, Any]]) -> Optional[np.ndarray]:
    seq = []
    for f in history:
        if isinstance(f, dict) and f.get("positions") is not None:
            p = np.asarray(f["positions"], dtype=float)
            if p.ndim == 2 and p.shape[1] == 2:
                seq.append(p)
    if len(seq) < 3:
        return None
    n = min(len(p) for p in seq)
    return np.stack([p[:n] for p in seq])                    # (T, N, 2)


def _velcorr_length(pos: np.ndarray, vel: np.ndarray) -> float:
    dv = vel - vel.mean(0)                                   # fluctuations about mean flock velocity
    norm = float((dv * dv).sum(1).mean())
    if norm < 1e-12:
        return 0.0
    d = pos[:, None, :] - pos[None, :, :]
    r = np.hypot(d[..., 0], d[..., 1])
    cij = (dv[:, None, :] * dv[None, :, :]).sum(-1) / norm
    L = float(r.max()) + 1e-9
    nb = 24
    iu = np.triu_indices(len(pos), k=1)
    rb = (r[iu] / L * nb).astype(int).clip(0, nb - 1)
    cv = cij[iu]
    C = np.array([cv[rb == b].mean() if np.any(rb == b) else np.nan for b in range(nb)])
    edges = (np.arange(nb) + 0.5) / nb
    # first zero-crossing of C(r)
    xi = 1.0
    for b in range(1, nb):
        if np.isnan(C[b]):
            continue
        if C[b] <= 0:
            xi = edges[b]
            break
    return float(xi)


def velocity_order(history: List[Dict[str, Any]], gap: int = 1) -> Optional[Dict[str, float]]:
    seq = _positions(history)
    if seq is None:
        return None
    T = seq.shape[0]
    pol, mill, xi = [], [], []
    for t in range(T - gap):
        pos = seq[t]
        vel = seq[t + gap] - seq[t]
        sp = np.hypot(vel[:, 0], vel[:, 1])
        if np.median(sp) < 1e-9:
            continue
        vhat = vel / (sp[:, None] + 1e-12)
        pol.append(float(np.hypot(*vhat.mean(0))))
        c = pos - pos.mean(0)
        cross = c[:, 0] * vel[:, 1] - c[:, 1] * vel[:, 0]
        denom = (np.hypot(c[:, 0], c[:, 1]) * sp).sum() + 1e-12
        mill.append(float(abs(cross.sum()) / denom))
        xi.append(_velcorr_length(pos, vel))
    if not pol:
        return None
    return {"velcorr_length": round(float(np.mean(xi)), 4),
            "polarization": round(float(np.mean(pol)), 4),
            "mill_strength": round(float(np.mean(mill)), 4)}
