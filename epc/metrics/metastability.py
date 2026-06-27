"""Metastability / regime-switching lens — a documented missing TEMPORAL axis. A system that
sits in macrostate A, then jumps to B, dwells, jumps back, shares its spectral / complexity /
tail signatures with a unimodal-stationary system; the difference lives in the macrostate TIME
SERIES. This lens reduces each frame to a scalar macrostate, then measures:

  n_macrostates  — number of separated modes in the macrostate distribution (KDE peaks split by
                   a deep valley). >=2 = the system visits distinct metastable states.
  dwell_cv       — coefficient of variation of dwell times (run lengths between mode crossings).
                   THE separator from a periodic oscillation: metastable switching has IRREGULAR
                   (Poisson-like, cv~1) dwells; an oscillation is multimodal too but with REGULAR
                   (cv~0) dwells, and is already owned by the spectral lens.
  bimodality     — Sarle's bimodality coefficient (continuous backup, >0.555 ~ multimodal).

Discriminator = (n_macrostates>=2) WITH high dwell_cv, NOT mere variance or a spectral peak.
A slow drift gives 1 broad mode (not switching); an oscillation gives regular dwells (not meta-
stable). Critical slowing down (early warning) is a separate item — it needs a tipping control.

Reads an explicit 'order_parameter'/'scalar' per frame, else derives one (Kuramoto R from
'phases', spatial mean of a 'field', radius of gyration from 'positions'). Needs >=12 frames.
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional

import numpy as np


def _macrostate_series(history: List[Dict[str, Any]]) -> Optional[np.ndarray]:
    s = []
    for f in history:
        if not isinstance(f, dict):
            continue
        if "order_parameter" in f:
            s.append(float(f["order_parameter"])); continue
        if "scalar" in f:
            s.append(float(f["scalar"])); continue
        if "phases" in f:
            th = np.asarray(f["phases"], dtype=float).ravel()
            s.append(float(abs(np.exp(1j * th).mean()))); continue
        got = False
        for k in ("field", "u", "concentration"):
            if k in f:
                s.append(float(np.asarray(f[k], dtype=float).mean())); got = True; break
        if got:
            continue
        if f.get("positions") is not None:
            p = np.asarray(f["positions"], dtype=float)
            if p.ndim == 2:
                s.append(float(np.hypot(*(p - p.mean(0)).T).mean()))
    s = np.asarray(s, dtype=float)
    return s if s.size >= 12 and np.ptp(s) > 0 else None


def _n_modes(z: np.ndarray) -> int:
    from scipy.ndimage import gaussian_filter1d
    hist, _ = np.histogram(z, bins=40, range=(-3, 3), density=True)
    h = gaussian_filter1d(hist, 2.0)
    peaks = [i for i in range(1, len(h) - 1)
             if h[i] > h[i - 1] and h[i] >= h[i + 1] and h[i] > 0.05 * h.max()]
    if len(peaks) <= 1:
        return max(1, len(peaks))
    modes = 1
    for j in range(1, len(peaks)):
        valley = h[peaks[j - 1]:peaks[j] + 1].min()
        if valley < 0.6 * min(h[peaks[j - 1]], h[peaks[j]]):
            modes += 1
    return modes


def _dwell_cv(z: np.ndarray) -> float:
    b = (z > np.median(z)).astype(int)
    runs, cur = [], 1
    for i in range(1, len(b)):
        if b[i] == b[i - 1]:
            cur += 1
        else:
            runs.append(cur); cur = 1
    runs.append(cur)
    runs = np.asarray(runs, dtype=float)
    if runs.size < 3:
        return 0.0
    return float(runs.std() / (runs.mean() + 1e-12))


def metastability(history: List[Dict[str, Any]]) -> Optional[Dict[str, float]]:
    s = _macrostate_series(history)
    if s is None:
        return None
    z = (s - s.mean()) / (s.std() + 1e-12)
    n = len(z)
    g = float(((z ** 3).mean()))                              # skewness (z already standardized)
    k = float(((z ** 4).mean()) - 3.0)                        # excess kurtosis
    bc = (g ** 2 + 1) / (k + 3.0 * (n - 1) ** 2 / ((n - 2) * (n - 3))) if n > 3 else 0.0
    return {"n_macrostates": _n_modes(z),
            "dwell_cv": round(_dwell_cv(z), 3),
            "bimodality": round(float(bc), 3)}
