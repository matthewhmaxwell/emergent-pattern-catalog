"""Substrate-agnostic information-theoretic detection channels for the emergence
instrument (core layer). Estimators validated on synthetic signals 2026-06-23
(noise→Psi_CE≈0, redundant→≪0, XOR→≈+1; complexity-entropy plane correct).

  psi_ce / psi_ce_best  — Rosas/Mediano synergy / causal-emergence criterion
                          Psi_CE = I(V_t;V_t') - sum_i I(X_i,t; V_t'). >0 ⇒
                          synergistic ("greater-than-the-parts") emergence that
                          clustering / autocorrelation / order-parameter channels
                          cannot see. (PLOS Comp Biol 2020; Phil Trans R Soc A 2022.)
  mpr_complexity        — Martin-Plastino-Rosso statistical complexity on
                          Bandt-Pompe ordinal patterns (structure ⟂ entropy).
  micro_macro           — generic per-component time-series + candidate macro
                          features extracted from any observation history, so the
                          channels apply without a per-system adapter.
"""
from __future__ import annotations

import itertools
import math
from typing import Any, Dict, List, Optional, Sequence, Tuple

import numpy as np


def _shannon_bits(counts: np.ndarray) -> float:
    c = np.asarray(counts, dtype=float)
    tot = c.sum()
    if tot <= 0:
        return 0.0
    p = c[c > 0] / tot
    return float(-(p * np.log2(p)).sum())


def _discretize(a: Sequence[float], bins: int) -> np.ndarray:
    """Integer codes. Low-cardinality (≤ bins distinct, e.g. binary) is FACTORIZED
    one-code-per-value — quantile-binning would collapse binary data to one bin and
    silently zero its mutual information. Higher cardinality is quantile-binned."""
    a = np.asarray(a, dtype=float).ravel()
    if a.size == 0:
        return a.astype(int)
    u = np.unique(a)
    if u.size <= bins:
        return np.searchsorted(u, a)
    edges = np.quantile(a, np.linspace(0, 1, bins + 1))
    edges = np.unique(edges)
    if edges.size < 3:
        return np.searchsorted(u, a)
    return np.clip(np.digitize(a, edges[1:-1]), 0, len(edges) - 2)


def _mi_bits(x: np.ndarray, y: np.ndarray, bx: int, by: int) -> float:
    n = len(x)
    if n < 2:
        return 0.0
    cxy = np.zeros((bx, by), dtype=float)
    np.add.at(cxy, (x, y), 1.0)
    return _shannon_bits(cxy.sum(1)) + _shannon_bits(cxy.sum(0)) - _shannon_bits(cxy.ravel())


def psi_ce(micro: np.ndarray, macro: Sequence[float], bins: int = 4,
           tau: int = 1, max_micro: int = 16) -> float:
    """Practical causal-emergence criterion (Rosas Eq. 10a), in bits. >0 ⇒ macro
    feature V is causally-emergent (synergistic). Sub-samples to max_micro parts."""
    micro = np.asarray(micro, dtype=float)
    macro = np.asarray(macro, dtype=float).ravel()
    T = macro.shape[0]
    if T - tau < 8 or micro.ndim != 2:
        return float("nan")
    n = micro.shape[1]
    idx = np.linspace(0, n - 1, min(n, max_micro)).astype(int) if n > max_micro else range(n)
    Vt = _discretize(macro[:-tau], bins)
    Vt2 = _discretize(macro[tau:], bins)
    bV = int(Vt2.max()) + 1
    i_vv = _mi_bits(Vt, Vt2, int(Vt.max()) + 1, bV)
    s = 0.0
    for i in idx:
        Xi = _discretize(micro[:-tau, i], bins)
        s += _mi_bits(Xi, Vt2, int(Xi.max()) + 1, bV)
    return float(i_vv - s)


def psi_ce_best(micro: np.ndarray, macro_candidates: Dict[str, Any], bins: int = 4,
                tau: int = 1) -> Tuple[float, Optional[str]]:
    """Max Psi_CE over a bank of candidate macro features (the instrument does not
    know V a priori). Returns (best_psi, feature_name)."""
    best, name = float("-inf"), None
    for nm, V in macro_candidates.items():
        v = psi_ce(micro, V, bins=bins, tau=tau)
        if not math.isnan(v) and v > best:
            best, name = v, nm
    return (best if best != float("-inf") else float("nan")), name


def _ordinal_distribution(series: np.ndarray, d: int, delay: int) -> np.ndarray:
    s = np.asarray(series, dtype=float).ravel()
    nperm = math.factorial(d)
    perms = {p: k for k, p in enumerate(itertools.permutations(range(d)))}
    counts = np.zeros(nperm, dtype=float)
    m = s.size - (d - 1) * delay
    for i in range(max(0, m)):
        window = s[i:i + d * delay:delay]
        counts[perms[tuple(np.argsort(window, kind="stable"))]] += 1.0
    return counts


def mpr_complexity(series: Sequence[float], d: int = 4, delay: int = 1) -> Dict[str, float]:
    """Martin-Plastino-Rosso statistical complexity via Bandt-Pompe patterns.
    Returns {H: normalized permutation entropy, C: statistical complexity}."""
    s = np.asarray(series, dtype=float).ravel()
    if s.size < (d - 1) * delay + 6:
        d = 3
    counts = _ordinal_distribution(s, d, delay)
    nperm = math.factorial(d)
    if counts.sum() <= 0:
        return {"H": 0.0, "C": 0.0, "d": d}
    P = counts / counts.sum()
    H = _shannon_bits(counts) / math.log2(nperm)
    Pe = np.full(nperm, 1.0 / nperm)
    M = 0.5 * (P + Pe)

    def _S(p):
        p = p[p > 0]
        return float(-(p * np.log2(p)).sum())
    js = _S(M) - 0.5 * _S(P) - 0.5 * _S(Pe)
    N = nperm
    q0 = -0.5 * (((N + 1) / N) * math.log2(N + 1) - 2 * math.log2(2 * N) + math.log2(N))
    return {"H": float(H), "C": float((js / q0 if q0 > 0 else 0.0) * H), "d": d}


def micro_macro(history: List[Dict[str, Any]], max_parts: int = 24
                ) -> Tuple[Optional[np.ndarray], Optional[Dict[str, Any]]]:
    """Generic per-component time-series + candidate macro features from any
    observation history. Returns (micro [T,n], candidates) or (None, None) when
    the frames carry no usable per-component observable / too few frames."""
    series: List[np.ndarray] = []
    for f in history:
        if not isinstance(f, dict):
            return None, None
        a = None
        if "grid" in f or "field" in f:
            a = np.asarray(f.get("grid", f.get("field")), dtype=float).ravel()
        elif "velocities" in f:
            v = np.asarray(f["velocities"], dtype=float)
            if v.ndim == 2 and v.shape[1] >= 2:
                a = np.arctan2(v[:, 1], v[:, 0])
        elif "positions" in f:
            p = np.asarray(f["positions"], dtype=float)
            if p.ndim == 2 and p.shape[1] >= 1:
                a = p[:, 0]
        elif "phases" in f or "theta" in f:
            a = np.asarray(f.get("phases", f.get("theta")), dtype=float).ravel()
        else:
            for k in ("state", "opinions", "wealth", "attendance", "x",
                      "fraction_on", "task_assignments", "cell_types"):
                if k in f:
                    a = np.asarray(f[k], dtype=float).ravel()
                    break
        if a is None or a.size == 0:
            return None, None
        series.append(a)
    if len(series) < 8:
        return None, None
    n = min(a.size for a in series)
    if n < 2:
        return None, None
    M = np.array([a[:n] for a in series], dtype=float)         # [T, n]
    if n > max_parts:
        M = M[:, np.linspace(0, n - 1, max_parts).astype(int)]
    cands: Dict[str, Any] = {"mean": M.mean(1), "std": M.std(1)}
    u = np.unique(np.round(M, 6))
    if u.size <= 2 and set(u.tolist()).issubset({0.0, 1.0}):
        cands["parity"] = M.astype(int).sum(1) % 2
        cands["majority"] = (M.mean(1) > 0.5).astype(int)
    return M, cands
