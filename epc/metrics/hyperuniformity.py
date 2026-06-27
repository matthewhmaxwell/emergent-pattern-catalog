"""Hyperuniformity / giant-number-fluctuations lens — long-wavelength density structure (a
documented missing axis). The covered density/clustering lenses see local structure; this reads
the SCALING of density fluctuations at large scales. Sub-box number variance sigma^2(R) ~ R^alpha:

  number_variance_exponent — alpha (2D):
      alpha ~ 2 (= d)      Poisson / generic disorder
      alpha -> 1 (= d-1)   HYPERUNIFORM: suppressed long-wavelength fluctuations (hidden order:
                           lattices, jammed packings) -- S(k->0)->0
      alpha > 2            clustered / GIANT number fluctuations (active matter) -- S(k->0) large
  THE discriminator is the long-wavelength variance SCALING EXPONENT, not the local density.

Reads 'positions' (Nx2) or a thresholded 2D field as a point pattern. Aggregates the late frames
for a stable estimate. Returns None if too few points. (Small-R/large-scale estimates are the
finite-size-sensitive part -- box sizes capped well below the domain.)
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional

import numpy as np


def _point_set(history: List[Dict[str, Any]]) -> Optional[np.ndarray]:
    late = [f for f in history if isinstance(f, dict)]
    late = late[len(late) // 2:] or late
    pts = []
    for f in late:
        if f.get("positions") is not None:
            p = np.asarray(f["positions"], dtype=float)
            if p.ndim == 2 and p.shape[1] == 2:
                pts.append(p)
    if pts:
        P = np.vstack(pts)
        return P if len(P) >= 64 else None
    # thresholded field -> point pattern (cell centres above median)
    for f in reversed(late):
        for k in ("field", "u", "concentration"):
            if k in f:
                a = np.asarray(f[k], dtype=float)
                if a.ndim == 2:
                    ys, xs = np.where(a > np.median(a))
                    if len(xs) >= 64:
                        return np.column_stack([xs.astype(float), ys.astype(float)])
    return None


def hyperuniformity(history: List[Dict[str, Any]], n_boxes: int = 400) -> Optional[Dict[str, float]]:
    P = _point_set(history)
    if P is None:
        return None
    lo = P.min(0); hi = P.max(0); L = hi - lo
    if np.any(L <= 0):
        return None
    rng = np.random.default_rng(0)
    fracs = np.array([0.08, 0.12, 0.18, 0.26, 0.35])
    Rs = fracs * L.min()
    var, mean = [], []
    for R in Rs:
        room = L - R
        if np.any(room <= 0):
            continue
        corners = lo + rng.uniform(0, 1, (n_boxes, 2)) * room
        counts = np.array([
            np.sum((P[:, 0] >= cx) & (P[:, 0] < cx + R) & (P[:, 1] >= cy) & (P[:, 1] < cy + R))
            for cx, cy in corners], dtype=float)
        var.append(counts.var()); mean.append(counts.mean())
    var = np.asarray(var); used = Rs[:len(var)]
    ok = var > 0
    if ok.sum() < 3:
        return None
    alpha = float(np.polyfit(np.log(used[ok]), np.log(var[ok]), 1)[0])
    return {"number_variance_exponent": round(alpha, 3)}
