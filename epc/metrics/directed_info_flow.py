"""Directed information-flow lens — transfer entropy with DIRECTION (asymmetry).

The existing transfer_entropy_global lens measures the MAGNITUDE of information
transfer on a grid CA. This lens measures its DIRECTIONALITY among component time
series: who drives whom. Two systems can share the same total TE but opposite
coordination structure — a leader->follower cascade (asymmetric flow) vs an
egalitarian diffusive mesh (symmetric flow) vs independence (no flow). Nothing else
in the battery sees that axis.

Pairwise boundary-free transfer entropy on quantile-symbolized series:
    TE(i->j) = sum p(j',j,i) log[ p(j'|j,i) / p(j'|j) ]    (j' = j_{t+1})

Scalars:
  mean_te       — average TE over all ordered pairs (coupling magnitude). Separates
                  coupled systems (cascade, mesh) from independent ones (~0).
  net_asymmetry — mean |TE(i->j) - TE(j->i)| over unordered pairs. HIGH for directed
                  cascades, ~0 for symmetric coupling and for independence.
  directionality— net_asymmetry / mean_te: the pure direction axis, scale-free in
                  coupling strength. HIGH cascade (~1.7), LOW mesh (~0.1). GATE: only
                  interpretable when mean_te shows real coupling — independent series
                  give ratio-noise ~0.3-0.4 on a ~0 denominator, so read directionality
                  only when mean_te is above the independence floor (~0.01).

Validated (coupled-logistic controls, T=2000): mean_te coupled 0.055-0.19 vs
independent ~0.004 (gap +0.05); directionality cascade 1.72-1.74 vs mesh 0.05-0.13
(gap +1.60) — the two axes are independent (mesh = high coupling, low direction).

Scope: needs enough samples for the symbol histograms — returns None when T < min_T.
Input is a [T, N] array of component time series (e.g., per-agent observables over
time); pass the relevant component matrix, not a flattened grid.
"""
from __future__ import annotations

from typing import Optional, Dict

import numpy as np


def _symbolize(x: np.ndarray, b: int) -> np.ndarray:
    ranks = np.argsort(np.argsort(x))               # rank-transform: robust to shape
    return np.clip(ranks * b // len(x), 0, b - 1).astype(int)


def _te(src: np.ndarray, tgt: np.ndarray, b: int) -> float:
    jp, j, i = tgt[1:], tgt[:-1], src[:-1]
    idx = (jp * b + j) * b + i
    C3 = np.bincount(idx, minlength=b ** 3).reshape(b, b, b).astype(float)  # [jp, j, i]
    N = C3.sum()
    if N == 0:
        return 0.0
    C_ji = C3.sum(0)                                 # [j, i]
    C2 = C3.sum(2)                                   # [jp, j]
    C1 = C3.sum((0, 2))                              # [j]
    nz = C3 > 0
    P3 = C3 / N
    with np.errstate(divide="ignore", invalid="ignore"):
        cond1 = C3 / C_ji[None, :, :]               # p(jp|j,i)
        cond2 = C2 / C1[None, :]                     # p(jp|j)
        ratio = cond1 / cond2[:, :, None]
    te = float(np.sum(P3[nz] * np.log(ratio[nz])))
    return max(te, 0.0)


def directed_transfer_entropy(series: np.ndarray, bins: int = 4,
                              max_components: int = 24, min_T: int = 60) -> Optional[Dict[str, float]]:
    X = np.asarray(series, dtype=float)
    if X.ndim != 2:
        return None
    T = X.shape[0]
    if T < min_T:
        return None
    X = X[:, X.std(0) > 1e-9]                        # drop constant series
    if X.shape[1] < 2:
        return None
    if X.shape[1] > max_components:
        X = X[:, np.linspace(0, X.shape[1] - 1, max_components).astype(int)]
    n = X.shape[1]
    S = np.column_stack([_symbolize(X[:, k], bins) for k in range(n)])
    te_mat = np.zeros((n, n))
    for i in range(n):
        for j in range(n):
            if i != j:
                te_mat[i, j] = _te(S[:, i], S[:, j], bins)
    off = ~np.eye(n, dtype=bool)
    mean_te = float(te_mat[off].mean())
    tri = np.triu_indices(n, 1)
    net_asym = float(np.abs(te_mat - te_mat.T)[tri].mean())
    return {"mean_te": round(mean_te, 5),
            "net_asymmetry": round(net_asym, 5),
            "directionality": round(net_asym / (mean_te + 1e-9), 4),
            "n_components": n}
