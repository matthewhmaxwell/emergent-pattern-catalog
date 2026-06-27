"""Chaos-dimension / determinism lens — a documented missing axis. Low-dimensional chaos,
spatiotemporal chaos, and noise all give a BROADBAND spectrum, so the covered spectral lens
cannot tell them apart. Two robust, complementary discriminators from delay embedding:

  zero_one_K   — Gottwald-Melbourne 0-1 test. K ~ 1 for chaotic/stochastic (mean-square
                 displacement grows with time), K ~ 0 for periodic/quasiperiodic (bounded).
  determinism  — nonlinear prediction skill (Sugihara-May): each point's next value is predicted
                 from its delay-embedding nearest neighbour's next value (Theiler-excluded).
                 Deterministic series (chaos, periodic) are predictable -> high corr; stochastic
                 noise is not -> ~0. THIS is what separates chaos from noise where the spectrum
                 cannot (GP correlation-dimension is too fragile on finite data -- the reason
                 recurrence was retired).

Joint readout: periodic = (K~0, determinism high); low-dim chaos = (K~1, determinism high);
noise = (K~1, determinism ~0). Reads an explicit 'scalar'/'order_parameter' series, else derives
one (field mean / Kuramoto R / radius of gyration). Needs >=400 samples; returns None when too
short so the estimate never fires on inadequate data.
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
    return s if s.size >= 400 and np.ptp(s) > 0 else None


def _zero_one(x: np.ndarray, n_c: int = 20) -> float:
    x = (x - x.mean())
    N = len(x); t = np.arange(1, N + 1)
    rng = np.random.default_rng(0)
    Ks = []
    for _ in range(n_c):
        c = float(rng.uniform(0.5, 2 * np.pi - 0.5))
        p = np.cumsum(x * np.cos(c * t)); q = np.cumsum(x * np.sin(c * t))
        ncut = N // 10
        M = np.array([np.mean((p[n:] - p[:-n]) ** 2 + (q[n:] - q[:-n]) ** 2) for n in range(1, ncut)])
        nn = np.arange(1, ncut)
        if M.std() > 0:
            Ks.append(float(np.corrcoef(nn, M)[0, 1]))
    return float(np.median(Ks)) if Ks else 0.0


def _embed(x: np.ndarray, m: int, tau: int) -> np.ndarray:
    M = len(x) - (m - 1) * tau
    return np.stack([x[i * tau:i * tau + M] for i in range(m)], axis=1)


def _pred_skill(x: np.ndarray, m: int = 3, tau: int = 1, theiler: int = 12, sub: int = 600) -> float:
    end = (m - 1) * tau
    M = len(x) - end - 1
    if M < 80:
        return 0.0
    Y = _embed(x, m, tau)[:M]                                 # each row -> next value x[i+end+1]
    tgt = x[end + 1:end + 1 + M]
    base = np.arange(M)
    if M > sub:
        sel = np.linspace(0, M - 1, sub).astype(int); Ys = Y[sel]; ts = tgt[sel]; bi = base[sel]
    else:
        Ys, ts, bi = Y, tgt, base
    d = np.sqrt(((Ys[:, None, :] - Ys[None, :, :]) ** 2).sum(-1))
    far = np.abs(bi[:, None] - bi[None, :]) <= theiler        # Theiler window (+ self)
    d[far] = np.inf
    nn = np.argmin(d, axis=1)
    pred = ts[nn]
    if ts.std() < 1e-9 or pred.std() < 1e-9:
        return 0.0
    return float(np.clip(np.corrcoef(pred, ts)[0, 1], 0.0, 1.0))


def chaos_dimension(history: List[Dict[str, Any]], tau: int = 1) -> Optional[Dict[str, float]]:
    x = _series(history)
    if x is None:
        return None
    x = (x - x.mean()) / (x.std() + 1e-12)
    K = _zero_one(x)
    det = _pred_skill(x, m=3, tau=tau)
    return {"zero_one_K": round(K, 3), "determinism": round(det, 3)}
