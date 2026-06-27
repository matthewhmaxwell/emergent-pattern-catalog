"""Synchronization-texture lens — chimera / cluster sync (a documented missing axis). The
covered global Kuramoto R conflates a chimera (coherent and incoherent domains coexisting in
space) with uniform partial sync: both can give the same R. This lens computes the LOCAL
order-parameter field R(x) = |<e^{i theta}>_neighbourhood| and reports its spatial spread.

  local_R_spread — std of the local-R field. THE discriminator: a chimera has R(x) spanning
                   ~0..1 (incoherent + coherent regions) -> HIGH spread; uniform sync (R~1
                   everywhere) and uniform incoherence (R~0 everywhere) -> LOW spread, even
                   though their GLOBAL R differs wildly.
  global_sync    — |<e^{i theta}>| over all sites (context; what the covered channel sees).
  local_R_bimod  — Sarle bimodality of the R(x) distribution (coherent+incoherent = bimodal).

Reads 2D (or 1D ring) 'phases'. Returns None otherwise. (theta_field = nematic director, not a
phase-sync substrate -> not read here; defect_census handles that.)
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional

import numpy as np


def _phases(frame: Dict[str, Any]) -> Optional[np.ndarray]:
    if isinstance(frame, dict) and "phases" in frame:
        a = np.asarray(frame["phases"], dtype=float)
        if a.ndim in (1, 2) and a.size >= 16:
            return a
    return None


def _local_R(theta: np.ndarray, win: int) -> np.ndarray:
    from scipy.ndimage import uniform_filter, uniform_filter1d
    c, s = np.cos(theta), np.sin(theta)
    if theta.ndim == 2:
        cc = uniform_filter(c, win, mode="wrap"); ss = uniform_filter(s, win, mode="wrap")
    else:
        cc = uniform_filter1d(c, win, mode="wrap"); ss = uniform_filter1d(s, win, mode="wrap")
    return np.sqrt(cc ** 2 + ss ** 2)


def sync_texture(history: List[Dict[str, Any]]) -> Optional[Dict[str, float]]:
    frames = [f for f in history if _phases(f) is not None]
    if not frames:
        return None
    late = frames[len(frames) // 2:] or frames
    spreads, glob, bim = [], [], []
    for f in late:
        th = _phases(f)
        win = max(3, min(th.shape) // 8) if th.ndim == 2 else max(3, th.size // 8)
        R = _local_R(th, win)
        spreads.append(float(R.std()))
        glob.append(float(abs(np.exp(1j * th).mean())))
        z = (R.ravel() - R.mean()) / (R.std() + 1e-12); n = z.size
        g = float((z ** 3).mean()); k = float((z ** 4).mean() - 3)
        bc = (g ** 2 + 1) / (k + 3.0 * (n - 1) ** 2 / ((n - 2) * (n - 3))) if n > 3 else 0.0
        bim.append(bc)
    return {"local_R_spread": round(float(np.mean(spreads)), 4),
            "global_sync": round(float(np.mean(glob)), 4),
            "local_R_bimod": round(float(np.mean(bim)), 3)}
