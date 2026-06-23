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


def mpr_emergence(micro: np.ndarray, n_surr: int = 24, seed: int = 0, d: int = 4
                  ) -> Tuple[float, float]:
    """Coordination-gated statistical-complexity channel for COLLECTIVE temporal
    structure (synchronized oscillation / limit cycles).

    Tests whether the collective signal (mean of components) has ordinal-pattern
    complexity exceeding surrogates that break cross-component coordination by
    circularly time-shifting each component independently (preserving each part's
    own dynamics). A synchronized oscillation collapses under desync → flagged; a
    smooth-but-uncoordinated signal (e.g. Brownian drift of a mean over independent
    random walkers) survives → NOT flagged (fixes the plain-shuffle false positive).
    Returns (z, observed_C). NOTE: the per-component/chaos sub-case (collective
    complexity vs independent units) is deferred — see blind_spot_audit report.
    """
    micro = np.asarray(micro, dtype=float)
    if micro.ndim != 2 or micro.shape[0] < 16 or micro.shape[1] < 2:
        return float("nan"), float("nan")
    rng = np.random.default_rng(seed)
    T, n = micro.shape
    c_obs = mpr_complexity(micro.mean(1), d=d)["C"]
    cs = np.empty(n_surr)
    for s in range(n_surr):
        sh = np.empty_like(micro)
        for j in range(n):
            sh[:, j] = np.roll(micro[:, j], int(rng.integers(0, T)))
        cs[s] = mpr_complexity(sh.mean(1), d=d)["C"]
    mu, sd = float(cs.mean()), float(cs.std())
    z = (c_obs - mu) / sd if sd > 1e-9 else (5.0 if c_obs - mu > 1e-6 else 0.0)
    return float(z), float(c_obs)


def plaw_llr(sizes: Sequence[float], xmin: float = 1.0) -> Optional[Tuple[float, float, float]]:
    """Power-law vs exponential log-likelihood ratio (Clauset-style, continuous
    approx). Returns (alpha, LLR, decades); LLR>0 ⇒ power-law preferred. None when
    too few events. A heavy-tailed (SOC) distribution gives LLR>0 over ≥~1.3 decades;
    exponential/uniform distributions give LLR<0."""
    x = np.asarray([s for s in sizes if s is not None and s >= xmin], dtype=float)
    if x.size < 30:
        return None
    n = x.size
    denom = np.sum(np.log(x / (xmin * 0.999)))
    if denom <= 0:
        return None
    alpha = 1.0 + n / denom
    lam = 1.0 / (x.mean() - xmin + 1e-9)
    ll_pl = np.sum(np.log((alpha - 1.0) / xmin) - alpha * np.log(x / xmin))
    ll_exp = np.sum(np.log(max(lam, 1e-12)) - lam * (x - xmin))
    decades = float(np.log10(x.max() / max(x.min(), 1e-9)))
    return float(alpha), float(ll_pl - ll_exp), decades


def heavy_tail_score(history: List[Dict[str, Any]]) -> Optional[Tuple[float, float]]:
    """Self-organized-criticality / heavy-tail channel: is the system's event-size
    distribution power-law (heavy-tailed across ≥~1.3 decades) rather than
    exponential? Derives a size distribution from explicit size observables
    (avalanche_sizes / activity) or, failing that, from per-step total change of a
    grid/field. Returns (LLR, decades) or None. Reads RAW frames (these observables
    are invisible to the spatial/phase channels)."""
    sizes: List[float] = []
    for f in history:
        if isinstance(f, dict) and "avalanche_sizes" in f:
            v = f["avalanche_sizes"]
            sizes.extend(list(v) if hasattr(v, "__len__") else [v])
    if len(sizes) < 30:
        act = [float(f["activity"]) for f in history if isinstance(f, dict) and "activity" in f]
        if len(act) >= 30:
            sizes = act
    if len(sizes) < 30:
        arrs = [np.asarray(f.get("grid", f.get("field")), dtype=float)
                for f in history if isinstance(f, dict) and ("grid" in f or "field" in f)]
        if len(arrs) >= 8:
            sizes = [float(np.abs(arrs[i] - arrs[i - 1]).sum()) for i in range(1, len(arrs))]
    sizes = [float(s) for s in sizes if s is not None and s > 0]
    if len(sizes) < 30:
        return None
    r = plaw_llr(sizes)
    if r is None:
        return None
    _, llr, decades = r
    return llr, decades


def oscillation_score(series: Sequence[float], kmin: int = 2) -> float:
    """Sustained-oscillation / limit-cycle detector: peak-to-mean of the power
    spectrum of the linearly-detrended series, EXCLUDING the lowest `kmin` bins.

    A limit cycle has a prominent interior spectral peak (→ large); Brownian drift
    has only a 1/f² low-frequency ramp removed by detrending + bin exclusion (→ ~1);
    white noise is spectrally flat (→ ~few). Distinguishes oscillation from the
    smooth-drift / noise nulls that fooled raw statistical complexity.
    """
    s = np.asarray(series, dtype=float).ravel()
    if s.size < 16 or np.ptp(s) < 1e-12:
        return 0.0
    t = np.arange(s.size)
    s = s - np.polyval(np.polyfit(t, s, 1), t)        # remove linear trend (drift)
    ps = np.abs(np.fft.rfft(s)) ** 2
    hi = max(kmin + 2, len(ps))
    ps = ps[kmin:hi]
    if ps.size < 3 or ps.sum() <= 0:
        return 0.0
    return float(ps.max() / (ps.mean() + 1e-12))


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
