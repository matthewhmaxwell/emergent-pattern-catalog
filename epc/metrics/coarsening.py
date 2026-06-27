"""Coarsening-dynamics lens — domain-growth exponent (a documented missing TEMPORAL axis).
structure_factor reports the characteristic scale at ONE snapshot; this reports its time
derivative. Track the principal radial peak k_peak(t) -> domain scale L(t) = 1/k_peak(t), fit
L(t) ~ t^n on log-log. The exponent n is the discriminator:
  n ~ 1/2  non-conserved (Allen-Cahn / curvature-driven) coarsening
  n ~ 1/3  conserved (Cahn-Hilliard / LSW) coarsening
  n ~ 0    a STATIC pattern (Turing, fixed wavelength) or an ARRESTED/glassy/pinned state
A static pattern and a coarsening one can share an instantaneous L; only the growth exponent
tells them apart, and arrest (n->0) is the glassiness signature.

Scalars:
  growth_exponent — n from the log L vs log t fit (THE discriminator; ~0 static, >~0.2 coarsening).
  fit_quality     — R^2 of the log-log fit (clean power law vs noise).
Reads a 2D 'field'/'u'/'concentration', a 'theta_field' (-> cos), or a density grid from
'positions', over >=6 frames. Returns None otherwise.
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional

import numpy as np


def _field_2d(frame: Dict[str, Any], G: int = 96) -> Optional[np.ndarray]:
    if not isinstance(frame, dict):
        return None
    for k in ("field", "u", "concentration"):
        if k in frame:
            a = np.asarray(frame[k], dtype=float)
            if a.ndim == 2 and min(a.shape) >= 12:
                return a
    if "theta_field" in frame:
        a = np.asarray(frame["theta_field"], dtype=float)
        if a.ndim == 2 and min(a.shape) >= 12:
            return np.cos(a)
    pos = frame.get("positions")
    if pos is not None:
        pos = np.asarray(pos, dtype=float)
        if pos.ndim == 2 and pos.shape[1] == 2 and len(pos) >= 16:
            lo = pos.min(0); span = np.ptp(pos, axis=0); span[span == 0] = 1.0
            ij = np.clip(((pos - lo) / span * (G - 1)).astype(int), 0, G - 1)
            grid = np.zeros((G, G)); np.add.at(grid, (ij[:, 1], ij[:, 0]), 1.0)
            return grid
    return None


def _kpeak(field: np.ndarray) -> Optional[float]:
    f = field - field.mean()
    if not np.any(f):
        return None
    P = np.abs(np.fft.fftshift(np.fft.fft2(f))) ** 2
    ny, nx = P.shape; cy, cx = ny // 2, nx // 2
    Y, X = np.indices(P.shape); R = np.hypot(X - cx, Y - cy).astype(int)
    rmax = min(cy, cx) - 1
    if rmax < 4:
        return None
    radial = np.array([P[R == rr].mean() for rr in range(1, rmax)])
    if not np.any(radial):
        return None
    pk = int(np.argmax(radial))                              # index -> r = pk+1
    lo = max(0, pk - 2); hi = min(len(radial), pk + 3)
    w = radial[lo:hi]; rs = np.arange(lo, hi) + 1.0
    if w.sum() <= 0:
        return float(pk + 1)
    return float((w * rs).sum() / w.sum())                   # sub-bin weighted centroid of k_peak


def coarsening(history: List[Dict[str, Any]]) -> Optional[Dict[str, float]]:
    frames = [f for f in history if isinstance(f, dict)]
    ks, ts = [], []
    for t, f in enumerate(frames):
        fld = _field_2d(f)
        if fld is None:
            continue
        kp = _kpeak(fld)
        if kp is not None and kp >= 1:
            ks.append(kp); ts.append(t + 1)
    if len(ks) < 6:
        return None
    L = 1.0 / np.asarray(ks, dtype=float)                     # domain scale ~ 1/k_peak
    t = np.asarray(ts, dtype=float)
    x = np.log(t); y = np.log(L)
    if x.std() < 1e-9 or y.std() < 1e-9:
        return {"growth_exponent": 0.0, "fit_quality": 0.0}
    n, b = np.polyfit(x, y, 1)
    yhat = n * x + b
    ss_res = float(((y - yhat) ** 2).sum()); ss_tot = float(((y - y.mean()) ** 2).sum())
    r2 = 1.0 - ss_res / (ss_tot + 1e-12)
    return {"growth_exponent": round(float(n), 3), "fit_quality": round(r2, 3)}
