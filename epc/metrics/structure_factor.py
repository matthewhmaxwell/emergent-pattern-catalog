"""Structure-factor peak — a Tier-2 lens imported from condensed matter.

The static structure factor S(k) = |FFT(density - mean)|^2 reveals characteristic
LENGTH SCALES: a sharp peak at nonzero k means a repeating spacing (Turing stripes,
a lattice, a selected wavelength), while a disordered gas / white-noise density has
a flat S(k). The lens scalar is the PRINCIPAL-PEAK PROMINENCE = max S(k>0) / mean
S(k>0), averaged over the late window.

This adds the "characteristic-scale order" axis the catalog lacks (the clustering-CV
morphology lens sees clumping but not a selected wavelength). It is a SPATIAL measure
(per-frame density field), so it does not inherit the trajectory-smoothness confound
that limited the recurrence lens. Works on positions (binned) or a field/grid.
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional

import numpy as np


def _density_grid(frame: Dict[str, Any], bins: int = 24) -> Optional[np.ndarray]:
    if not isinstance(frame, dict):
        return None
    if "positions" in frame:
        p = np.asarray(frame["positions"], dtype=float)
        if p.ndim != 2 or p.shape[1] < 2:
            return None
        p = p[:, :2]
        rng = [[float(p[:, 0].min()), float(p[:, 0].max()) + 1e-9],
               [float(p[:, 1].min()), float(p[:, 1].max()) + 1e-9]]
        return np.histogram2d(p[:, 0], p[:, 1], bins=bins, range=rng)[0]
    for k in ("field", "grid"):
        if k in frame:
            a = np.asarray(frame[k], dtype=float)
            return a if a.ndim == 2 else None
    return None


def _sk_peak(grid: np.ndarray) -> float:
    a = grid - grid.mean()
    if a.std() < 1e-12:
        return 1.0
    S = np.abs(np.fft.fftshift(np.fft.fft2(a))) ** 2
    ny, nx = S.shape
    S[ny // 2, nx // 2] = 0.0                     # remove k=0 (mean)
    pos = S[S > 0]
    if pos.size == 0:
        return 1.0
    return float(S.max() / (pos.mean() + 1e-12))  # principal-peak prominence


def structure_factor_peak(history: List[Dict[str, Any]], bins: int = 24) -> Optional[Dict[str, float]]:
    grids = [_density_grid(f, bins) for f in history if isinstance(f, dict)]
    grids = [g for g in grids if g is not None]
    if len(grids) < 2:
        return None
    late = grids[len(grids) // 2:]
    vals = [_sk_peak(g) for g in late]
    return {"sk_peak": round(float(np.mean(vals)), 3)}
