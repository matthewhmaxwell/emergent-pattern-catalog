"""Spatial correlation-length / criticality lens — long-range, scale-free spatial
correlations (a documented under-covered axis: correlation-length divergence).

STATUS: DEFERRED (not wired into the descriptor). Validation showed the criticality
discriminator (`plaw_gain`, power-law vs exponential correlation decay) does NOT recover
cleanly from a single ~96^2 snapshot's ACF — it inverted (power-law fields scored negative,
white/exponential positive), because power-law-vs-exponential separation is a finite-size-
scaling problem (multi-resolution), not a one-frame fit. The part that does work, `xi_norm`
(correlation length), largely duplicates the existing Moran's-I autocorrelation channel, so
it adds no clean new axis. Revisit with a proper finite-size-scaling estimator (correlation
length across grid sizes) if the criticality axis is needed. Kept as raw material.

Distinct from structure_factor (which keys on a PEAK = characteristic spacing) and from
Moran's I (a single short-range autocorrelation): this measures HOW FAR order persists and
whether the spatial correlation decays as a POWER LAW (scale-free / critical) versus
exponentially (disordered) or not at all (long-range order).

Radial autocorrelation C(r) via Wiener-Khinchin (C = IFFT(|FFT(field-mean)|^2), radially
averaged, normalized to C(0)=1). Scalars (late-window mean):
  xi_norm   — correlation length / grid half-size: the r at which C(r) falls to 1/e. White
              noise -> ~0; long-range order -> ~1.
  plaw_gain — R^2(power-law fit) - R^2(exponential fit) of C(r) over its decaying range.
              > 0 => scale-free (critical) correlations; <= 0 => exponential/disordered.
              THE discriminator for the criticality axis.
Reads a 'field'/'grid', or rasterizes 'positions'. Returns None otherwise.
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional

import numpy as np


def _grid(frame: Dict[str, Any], grid_n: int = 96) -> Optional[np.ndarray]:
    if "field" in frame:
        a = np.asarray(frame["field"], dtype=float)
        return a if a.ndim == 2 else None
    if "grid" in frame:
        a = np.asarray(frame["grid"], dtype=float)
        return a if a.ndim == 2 else None
    if "positions" in frame:
        p = np.asarray(frame["positions"], dtype=float)
        if p.ndim != 2 or p.shape[1] < 2 or len(p) < 5:
            return None
        p = p[:, :2]
        lo, span = p.min(0), (p.max(0) - p.min(0)); span[span == 0] = 1.0
        ij = np.clip(np.floor((p - lo) / span * (grid_n - 1)).astype(int), 0, grid_n - 1)
        g = np.zeros((grid_n, grid_n)); g[ij[:, 0], ij[:, 1]] = 1.0
        return g
    return None


def _radial_acf(a: np.ndarray) -> np.ndarray:
    a = a - a.mean()
    if a.std() < 1e-12:
        return np.array([1.0])
    p = np.abs(np.fft.fft2(a)) ** 2
    c = np.fft.ifft2(p).real
    c = np.fft.fftshift(c)
    c = c / c.max()
    ny, nx = c.shape; cy, cx = ny // 2, nx // 2
    y, x = np.indices(c.shape)
    r = np.sqrt((y - cy) ** 2 + (x - cx) ** 2).astype(int)
    rmax = min(cy, cx)
    return np.array([c[r == rr].mean() for rr in range(rmax)])


def _fit_r2(x: np.ndarray, y: np.ndarray) -> float:
    if len(x) < 3:
        return -1.0
    s, b = np.polyfit(x, y, 1)
    yhat = s * x + b
    ss_res = float(((y - yhat) ** 2).sum()); ss_tot = float(((y - y.mean()) ** 2).sum())
    return 1.0 - ss_res / ss_tot if ss_tot > 0 else -1.0


def correlation_length(history: List[Dict[str, Any]], grid_n: int = 96) -> Optional[Dict[str, float]]:
    frames = [f for f in history if isinstance(f, dict)]
    if not frames:
        return None
    late = frames[len(frames) // 2:] or frames
    xis, gains = [], []
    for f in late:
        g = _grid(f, grid_n)
        if g is None:
            continue
        acf = _radial_acf(g)
        if acf.size < 4:
            continue
        # correlation length: first r where ACF <= 1/e
        below = np.where(acf <= np.exp(-1.0))[0]
        xi = float(below[0]) if below.size else float(acf.size)
        xis.append(xi / acf.size)
        # power-law vs exponential fit over the positive decaying range (r>=1, C>0)
        r = np.arange(1, acf.size)
        c = acf[1:]
        m = c > 1e-3
        if m.sum() >= 4:
            rr, cc = r[m], c[m]
            r2_pl = _fit_r2(np.log(rr), np.log(cc))            # power law: logC vs logr
            r2_ex = _fit_r2(rr.astype(float), np.log(cc))      # exponential: logC vs r
            gains.append(r2_pl - r2_ex)
    if not xis:
        return None
    return {"xi_norm": round(float(np.mean(xis)), 4),
            "plaw_gain": round(float(np.mean(gains)), 4) if gains else 0.0}
