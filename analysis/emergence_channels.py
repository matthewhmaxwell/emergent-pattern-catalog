"""Candidate general-purpose (substrate-agnostic) emergence-detection channels,
to extend the morphology-specific generic-emergence indicator.

Two channels from the 2026-06-23 emergence-channel research (25/25 verified):

  psi_ce          — Rosas/Mediano practical causal-emergence criterion
                    Psi_CE = I(V_t; V_t') - sum_i I(X_i,t; V_t').  > 0 is a
                    SUFFICIENT (not necessary) signature of SYNERGISTIC emergence
                    ("greater-than-the-parts" / XOR-type collective effects that
                    correlation, clustering, Moran's I, Kuramoto-r and nematic
                    order all miss). Model-free: needs only micro time-series X_i
                    and a candidate macro feature V. (PLOS Comp Biol 2020; ΦID
                    Phil Trans R Soc A 2022.)

  mpr_complexity  — Martin-Plastino-Rosso statistical complexity on Bandt-Pompe
                    ordinal patterns: C_JS = Q_J[P,P_e] * H_S[P]. A member of the
                    computational-mechanics / structural-complexity family that is
                    ORTHOGONAL to entropy (zero at both order and randomness
                    extremes, peaks in between) — the axis a pure-entropy channel
                    cannot see. (Rosso et al. 2007; Lopez-Ruiz-Mancini-Calbet.)

Both operate on dynamics, not spatial shape, so they catch organization without
assuming what it looks like.
"""
from __future__ import annotations

import itertools
import math
from typing import List, Optional, Sequence

import numpy as np


# ---------------------------------------------------------------------------
# discrete information helpers
# ---------------------------------------------------------------------------
def _shannon_bits(counts: np.ndarray) -> float:
    c = np.asarray(counts, dtype=float)
    tot = c.sum()
    if tot <= 0:
        return 0.0
    p = c[c > 0] / tot
    return float(-(p * np.log2(p)).sum())


def _discretize(a: Sequence[float], bins: int) -> np.ndarray:
    """Map a 1-D series to integer codes. Low-cardinality data (≤ bins distinct
    values, e.g. binary/categorical) is FACTORIZED to one code per value — this
    is essential: quantile-binning collapses binary data to a single bin and
    silently zeros its mutual information. Higher-cardinality data is quantile-
    binned to `bins` codes (robust to scale)."""
    a = np.asarray(a, dtype=float).ravel()
    if a.size == 0:
        return a.astype(int)
    u = np.unique(a)
    if u.size <= bins:
        return np.searchsorted(u, a)            # factorize categorical/binary
    edges = np.quantile(a, np.linspace(0, 1, bins + 1))
    edges = np.unique(edges)
    if edges.size < 3:
        return np.searchsorted(u, a)
    return np.clip(np.digitize(a, edges[1:-1]), 0, len(edges) - 2)


def _mi_bits(x: np.ndarray, y: np.ndarray, bx: int, by: int) -> float:
    """Mutual information (bits) of two integer-coded series via the plug-in
    estimator with a Miller-Madow-style bias note (kept simple; thresholds are
    calibrated against nulls downstream)."""
    n = len(x)
    if n < 2:
        return 0.0
    cxy = np.zeros((bx, by), dtype=float)
    np.add.at(cxy, (x, y), 1.0)
    return _shannon_bits(cxy.sum(1)) + _shannon_bits(cxy.sum(0)) - _shannon_bits(cxy.ravel())


# ---------------------------------------------------------------------------
# Channel 1 — Ψ_CE (synergy / causal emergence)
# ---------------------------------------------------------------------------
def psi_ce(micro: np.ndarray, macro: Sequence[float], bins: int = 4,
           tau: int = 1, max_micro: int = 16) -> float:
    """Practical causal-emergence criterion (Rosas Eq. 10a).

    micro : array [T, n]   per-part time series.
    macro : array [T]      candidate macroscopic feature V.
    Returns Psi_CE in bits. > 0  ⇒ V is a causally-emergent (synergistic)
    feature — its self-predictive information exceeds the sum the parts carry
    about its future. Sub-samples to max_micro parts to keep the sum estimable.
    """
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


def psi_ce_best(micro: np.ndarray, macro_candidates: dict, bins: int = 4,
                tau: int = 1) -> tuple:
    """Max Ψ_CE over a bank of candidate macro features (the instrument does not
    know V a priori, so it tries generic candidates). Returns (best_psi, name)."""
    best, name = float("-inf"), None
    for nm, V in macro_candidates.items():
        v = psi_ce(micro, V, bins=bins, tau=tau)
        if not math.isnan(v) and v > best:
            best, name = v, nm
    return (best if best != float("-inf") else float("nan")), name


# ---------------------------------------------------------------------------
# Channel 2 — MPR statistical complexity (structure orthogonal to entropy)
# ---------------------------------------------------------------------------
def _ordinal_distribution(series: np.ndarray, d: int, delay: int) -> np.ndarray:
    s = np.asarray(series, dtype=float).ravel()
    nperm = math.factorial(d)
    perms = {p: k for k, p in enumerate(itertools.permutations(range(d)))}
    counts = np.zeros(nperm, dtype=float)
    m = s.size - (d - 1) * delay
    if m <= 0:
        return counts
    for i in range(m):
        window = s[i:i + d * delay:delay]
        pat = tuple(np.argsort(window, kind="stable"))
        counts[perms[pat]] += 1.0
    return counts


def mpr_complexity(series: Sequence[float], d: int = 4, delay: int = 1) -> dict:
    """Martin-Plastino-Rosso statistical complexity via Bandt-Pompe patterns.

    Returns {H: normalized permutation entropy in [0,1], C: statistical
    complexity C_JS in [0,~1]}. C is high only for series that are structured
    but not trivially ordered or random (the organization axis). Falls back to
    d=3 when the series is too short for d=4.
    """
    s = np.asarray(series, dtype=float).ravel()
    if s.size < (d - 1) * delay + 6:
        d = 3
    counts = _ordinal_distribution(s, d, delay)
    nperm = math.factorial(d)
    tot = counts.sum()
    if tot <= 0:
        return {"H": 0.0, "C": 0.0, "d": d}
    P = counts / tot
    H_s = _shannon_bits(counts)                       # bits
    H = H_s / math.log2(nperm)                         # normalized 0..1
    # Jensen-Shannon disequilibrium to the uniform distribution
    Pe = np.full(nperm, 1.0 / nperm)
    M = 0.5 * (P + Pe)

    def _S(p):
        p = p[p > 0]
        return float(-(p * np.log2(p)).sum())
    js = _S(M) - 0.5 * _S(P) - 0.5 * _S(Pe)
    # normalization constant Q0 (Martin-Plastino-Rosso) so Q_J in [0,1]
    N = nperm
    q0 = -0.5 * (((N + 1) / N) * math.log2(N + 1) - 2 * math.log2(2 * N) + math.log2(N))
    Q_J = js / q0 if q0 > 0 else 0.0
    return {"H": float(H), "C": float(Q_J * H), "d": d}
