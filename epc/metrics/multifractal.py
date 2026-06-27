"""Multifractal-spectrum lens — closes the DENSE-fractal gap the fractal_dimension lens deferred
(it was scoped to sparse fractals via lacunarity; a dense heterogeneous measure has box-counting
D ~ 2 and is invisible to it). The multifractal spectrum width Delta-alpha captures heterogeneous
SCALING of the measure: a monofractal / uniform measure has tau(q) linear -> Delta-alpha ~ 0; a
multiplicative cascade has tau(q) nonlinear -> Delta-alpha wide.

  mf_width — Delta-alpha(field) - Delta-alpha(spatially shuffled field). SURROGATE-GATED so an
             iid heavy-tailed measure (same value histogram, NO spatial structure) does not
             masquerade as multifractal: shuffling leaves its Delta-alpha unchanged -> mf_width~0,
             while a real spatial cascade loses its multifractality under shuffle -> mf_width>0.
             THE discriminator (structural multifractality beyond the value distribution).
  alpha_span — raw Delta-alpha (context).

Reads a 2D 'field'/'u'/'concentration' (-> nonneg measure) or a density grid from 'positions'.
Returns None otherwise. (The queue-flagged fat-tail artifact is exactly what the surrogate removes.)
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional

import numpy as np


def _measure(frame: Dict[str, Any], G: int = 96) -> Optional[np.ndarray]:
    if not isinstance(frame, dict):
        return None
    for k in ("field", "u", "concentration"):
        if k in frame:
            a = np.asarray(frame[k], dtype=float)
            if a.ndim == 2 and min(a.shape) >= 16:
                return a
    pos = frame.get("positions")
    if pos is not None:
        pos = np.asarray(pos, dtype=float)
        if pos.ndim == 2 and pos.shape[1] == 2 and len(pos) >= 64:
            lo = pos.min(0); span = np.ptp(pos, axis=0); span[span == 0] = 1.0
            ij = np.clip(((pos - lo) / span * (G - 1)).astype(int), 0, G - 1)
            grid = np.zeros((G, G)); np.add.at(grid, (ij[:, 1], ij[:, 0]), 1.0)
            return grid
    return None


def _alpha_span(mu: np.ndarray) -> float:
    mu = mu - mu.min()
    s = mu.sum()
    if s <= 0:
        return 0.0
    mu = mu / s
    G = min(mu.shape)
    sizes = [s for s in (2, 4, 8, 16) if s <= G // 4]
    if len(sizes) < 3:
        return 0.0
    qs = np.array([-3, -2, -1, 0.0, 1.0001, 2, 3], dtype=float)   # cap moments (heavy tails unstable at |q|>=4)
    tau = []
    for q in qs:
        logZ, logE = [], []
        for s in sizes:
            ny, nx = (mu.shape[0] // s) * s, (mu.shape[1] // s) * s
            blocks = mu[:ny, :nx].reshape(ny // s, s, nx // s, s).sum(axis=(1, 3))
            p = blocks[blocks > 0]
            if p.size == 0:
                continue
            Z = np.sum(p ** q)
            if Z > 0:
                logZ.append(np.log(Z)); logE.append(np.log(s / G))
        if len(logZ) >= 3:
            tau.append(np.polyfit(logE, logZ, 1)[0])
        else:
            tau.append(np.nan)
    tau = np.array(tau)
    ok = ~np.isnan(tau)
    if ok.sum() < 4:
        return 0.0
    alpha = np.gradient(tau[ok], qs[ok])
    return float(np.nanmax(alpha) - np.nanmin(alpha))


def multifractal(history: List[Dict[str, Any]]) -> Optional[Dict[str, float]]:
    frames = [f for f in history if _measure(f) is not None]
    if not frames:
        return None
    late = frames[len(frames) // 2:] or frames
    raw, net = [], []
    rng = np.random.default_rng(0)
    for f in late[:3]:
        mu = _measure(f)
        a = _alpha_span(mu)
        sh = []
        for _ in range(5):                                   # multi-shuffle baseline (stable for heavy tails)
            flat = mu.ravel().copy(); rng.shuffle(flat)
            sh.append(_alpha_span(flat.reshape(mu.shape)))
        raw.append(a); net.append(a - float(np.mean(sh)))
    if not raw:
        return None
    return {"mf_width": round(float(np.mean(net)), 4),
            "alpha_span": round(float(np.mean(raw)), 4)}
