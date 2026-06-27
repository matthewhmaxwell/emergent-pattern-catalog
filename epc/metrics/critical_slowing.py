"""Critical-slowing-down / early-warning lens (spun out of metastability). Approaching a tipping
point (saddle-node / fold), the recovery rate -> 0, so on the macrostate series the lag-1
autocorrelation rises toward 1 and the variance rises -- BEFORE the mean shifts. This is the
distinctive pre-tipping signal, invisible to static level or spectral measures.

  ews_ar1_trend — Kendall-tau trend of windowed lag-1 autocorrelation of the DETRENDED series.
                  Rising (->1) = critical slowing down. THE discriminator. Detrending (subtract a
                  slow gaussian-smoothed component) defuses the queue-flagged confound: a slow
                  non-critical DRIFT in the mean leaves residual AR(1) flat -> no false trip.
  ews_var_trend — Kendall-tau trend of windowed variance of residuals (rising = CSD, context).

Reads an explicit 'scalar'/'order_parameter' series, else derives one (field mean / Kuramoto R /
radius of gyration). Needs >=60 frames for windowing. Returns None otherwise.
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional

import numpy as np


def _series(history: List[Dict[str, Any]]) -> Optional[np.ndarray]:
    s = []
    for f in history:
        if not isinstance(f, dict):
            continue
        if "scalar" in f:
            s.append(float(f["scalar"])); continue
        if "order_parameter" in f:
            s.append(float(f["order_parameter"])); continue
        if "phases" in f:
            th = np.asarray(f["phases"], dtype=float).ravel()
            s.append(float(abs(np.exp(1j * th).mean()))); continue
        done = False
        for k in ("field", "u", "concentration"):
            if k in f:
                s.append(float(np.asarray(f[k], dtype=float).mean())); done = True; break
        if done:
            continue
        if f.get("positions") is not None:
            p = np.asarray(f["positions"], dtype=float)
            if p.ndim == 2:
                s.append(float(np.hypot(*(p - p.mean(0)).T).mean()))
    s = np.asarray(s, dtype=float)
    return s if s.size >= 60 and np.ptp(s) > 0 else None


def _detrend(seg: np.ndarray) -> np.ndarray:
    t = np.arange(len(seg))
    return seg - np.polyval(np.polyfit(t, seg, 1), t)        # per-segment LINEAR detrend (kills a drift ramp)


def _ar1(seg: np.ndarray) -> float:
    seg = _detrend(seg)
    if seg.std() < 1e-12 or len(seg) < 4:
        return 0.0
    return float(np.corrcoef(seg[:-1], seg[1:])[0, 1])


def critical_slowing(history: List[Dict[str, Any]]) -> Optional[Dict[str, float]]:
    x = _series(history)
    if x is None:
        return None
    n = len(x); k = n // 3
    early, late = x[:k], x[-k:]
    ar1_rise = _ar1(late) - _ar1(early)
    ve = _detrend(early).var(); vl = _detrend(late).var()
    var_rise = float(np.log((vl + 1e-12) / (ve + 1e-12)))
    return {"ews_ar1_rise": round(ar1_rise, 3), "ews_var_rise": round(var_rise, 3)}
