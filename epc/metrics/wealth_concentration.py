"""Wealth concentration and inequality metrics.

These metrics quantify the distribution of a conserved scalar resource
(wealth) across a population of agents.  They are substrate-agnostic
and operate on a single (N,) non-negative float array, or on a
trajectory of such arrays.

Core metrics:
  - gini(w)                            Gini coefficient in [0, 1]
  - top_p_share(w, p)                  Fraction held by top p of agents
  - hill_tail_alpha(w, k)              Hill estimator of Pareto tail α
  - pareto_ks_distance(w, alpha, x_m)  KS distance to fitted Pareto
  - monotonic_growth_fraction(ginis)   Fraction of checkpoints where
                                       Gini does not decrease
  - well_mixed_gini_null(N, mean_w, n) Sample Ginis under the
                                       Boltzmann-Gibbs (Exp) null

The canonical P28 signature is simultaneous:
  - Gini -> 1 at long time
  - top-p shares concentrate monotonically
  - Gini exceeds the well-mixed Exp(mean_w) null by >> 1 null-std

Reference:
    Dragulescu, A. & Yakovenko, V. M. (2001). Statistical mechanics of
        money. European Physical Journal B, 17, 723-729.
    Boghosian, B. M. (2014). Kinetics of wealth and the Pareto law.
        Physical Review E, 89, 042804.
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray


def gini(w: NDArray[np.float64]) -> float:
    """Gini coefficient via sorted-order formula.

    For a non-negative vector w of length N,
        G = (2 * sum_{i=1}^{N} i * w_(i)) / (N * sum(w)) - (N+1)/N
    where w_(i) are the sorted values (ascending). O(N log N).

    Returns 0 for perfect equality, 1 for maximal concentration.
    Robust to zero-sum inputs (returns 0.0).
    """
    w = np.asarray(w, dtype=np.float64)
    if w.size == 0:
        return 0.0
    total = float(w.sum())
    if total <= 0:
        return 0.0
    sw = np.sort(w)
    N = len(sw)
    idx = np.arange(1, N + 1, dtype=np.float64)
    return float((2.0 * (idx * sw).sum()) / (N * total) - (N + 1.0) / N)


def top_p_share(w: NDArray[np.float64], p: float) -> float:
    """Fraction of total wealth held by the top fraction p of agents.

    p must be in (0, 1]. Returns 0 for zero-sum inputs. Uses
    k = max(1, round(p*N)) so p * N < 1 returns the single-top agent
    share.
    """
    w = np.asarray(w, dtype=np.float64)
    if w.size == 0 or p <= 0 or p > 1:
        return 0.0
    total = float(w.sum())
    if total <= 0:
        return 0.0
    N = len(w)
    k = max(1, int(round(p * N)))
    sorted_desc = np.sort(w)[::-1]
    return float(sorted_desc[:k].sum() / total)


def hill_tail_alpha(
    w: NDArray[np.float64],
    k: int | None = None,
    min_k: int = 20,
) -> tuple[float, float, int]:
    """Hill estimator of Pareto tail exponent alpha.

    Given sorted-descending values w_(1) >= ... >= w_(k) >= w_(k+1),
        alpha_hat = 1 / [(1/k) * sum_{i=1}^k ln(w_(i) / w_(k+1))]

    Parameters
    ----------
    w : array-like
        Non-negative values (wealth, rank sizes, etc.).
    k : int or None
        Number of upper order statistics to use. If None, defaults
        to max(min_k, N // 10) — the top 10% (bounded below).
    min_k : int
        Floor on k. Below this, the estimator is unreliable.

    Returns
    -------
    (alpha, w_min, k_used) : tuple
        alpha = NaN if fit degenerate (zero w_min, zero-mean log ratio,
        insufficient data).
    """
    w = np.asarray(w, dtype=np.float64)
    w = w[w > 0]
    n = len(w)
    if n < min_k + 1:
        return (float("nan"), float("nan"), 0)
    if k is None:
        k = max(min_k, n // 10)
    k = min(k, n - 1)
    w_sorted_desc = np.sort(w)[::-1]
    top_k = w_sorted_desc[:k]
    w_min = float(w_sorted_desc[k])
    if w_min <= 0:
        return (float("nan"), w_min, k)
    log_ratio = np.log(top_k / w_min)
    mean_log = float(log_ratio.mean())
    if mean_log <= 0:
        return (float("nan"), w_min, k)
    return (1.0 / mean_log, w_min, k)


def pareto_ks_distance(
    w: NDArray[np.float64],
    alpha: float,
    w_min: float,
) -> float:
    """Kolmogorov-Smirnov distance between upper tail of w and Pareto(alpha, w_min).

    Returns NaN if alpha or w_min is not finite or tail sample too small.
    """
    if not (np.isfinite(alpha) and np.isfinite(w_min)) or w_min <= 0 or alpha <= 0:
        return float("nan")
    w = np.asarray(w, dtype=np.float64)
    tail = np.sort(w[w >= w_min])
    if len(tail) < 10:
        return float("nan")
    n = len(tail)
    empirical = np.arange(1, n + 1, dtype=np.float64) / n
    theoretical = 1.0 - (w_min / tail) ** alpha
    return float(np.max(np.abs(empirical - theoretical)))


def monotonic_growth_fraction(gini_series: NDArray[np.float64]) -> float:
    """Fraction of consecutive checkpoints where Gini does NOT decrease.

    Under pure Yard-Sale condensation, Gini is monotonically
    non-decreasing up to stochastic noise. A value near 1.0 is a
    confirmation that the system is in condensation mode; values near
    0.5 indicate equilibrium fluctuations (consistent with a saturated
    or redistributive regime).

    Uses a small tolerance (1e-4) to ignore floating-point / sampling
    jitter.
    """
    g = np.asarray(gini_series, dtype=np.float64)
    if g.size < 2:
        return float("nan")
    diffs = g[1:] - g[:-1]
    non_decreasing = (diffs >= -1e-4).sum()
    return float(non_decreasing) / float(len(diffs))


def relative_gini_growth(
    gini_series: NDArray[np.float64],
    elapsed_steps_series: NDArray[np.int64] | None = None,
) -> float:
    """Relative Gini growth per checkpoint: (G_end - G_start) / (1 - G_start).

    Normalized so 0 = no growth, 1 = full condensation from start.
    Negative = inequality decreased.

    Returns NaN for zero/one-point series.
    """
    g = np.asarray(gini_series, dtype=np.float64)
    if g.size < 2:
        return float("nan")
    g0 = float(g[0])
    g1 = float(g[-1])
    denom = max(1.0 - g0, 1e-6)
    return (g1 - g0) / denom


def well_mixed_gini_null(
    N: int,
    mean_w: float,
    n_samples: int,
    rng: np.random.Generator | None = None,
) -> NDArray[np.float64]:
    """Sample Gini values from the Boltzmann-Gibbs well-mixed null.

    The Dragulescu-Yakovenko (2001) result: a conserved scalar resource
    exchanged by an ergodic, well-mixing dynamics equilibrates to an
    exponential distribution P(w) = (1/<w>) exp(-w/<w>), analogous to
    the Boltzmann distribution of energy in a gas. The Gini of Exp(β)
    is exactly 0.5 in the large-N limit.

    We sample N values from Exp(mean_w) and compute Gini, repeated
    n_samples times. This gives the null distribution of Gini under
    "symmetric exchange produced a well-mixed equilibrium" — the naive
    prediction for a fair-exchange economy. A Gini > null + 3*null_std
    is evidence that the system is CONDENSING beyond Boltzmann-Gibbs.

    Parameters
    ----------
    N : int
        Population size matching the observed data.
    mean_w : float
        Mean wealth in the observed data.
    n_samples : int
        Number of null Gini values to draw.
    rng : np.random.Generator or None
        Random number generator. If None, uses default_rng().

    Returns
    -------
    ginis : np.ndarray shape (n_samples,)
    """
    if rng is None:
        rng = np.random.default_rng()
    if N < 1 or mean_w <= 0 or n_samples < 1:
        return np.empty(0)
    out = np.empty(n_samples, dtype=np.float64)
    for s in range(n_samples):
        samp = rng.exponential(scale=mean_w, size=N)
        out[s] = gini(samp)
    return out


def concentration_trajectory(
    wealth_series: list[NDArray[np.float64]],
) -> dict[str, NDArray[np.float64]]:
    """Compute concentration metrics along a trajectory.

    Parameters
    ----------
    wealth_series : list of (N,) arrays
        Wealth snapshots at each checkpoint.

    Returns
    -------
    dict with keys:
      - gini            (T,)
      - max_share       (T,)
      - top_1pct_share  (T,)
      - top_10pct_share (T,)
      - mean_wealth     (T,)
      - alpha_hill      (T,)
    """
    T = len(wealth_series)
    out = {
        "gini": np.zeros(T),
        "max_share": np.zeros(T),
        "top_1pct_share": np.zeros(T),
        "top_10pct_share": np.zeros(T),
        "mean_wealth": np.zeros(T),
        "alpha_hill": np.full(T, np.nan),
    }
    for t, w in enumerate(wealth_series):
        w = np.asarray(w, dtype=np.float64)
        if w.size == 0:
            continue
        out["gini"][t] = gini(w)
        total = float(w.sum())
        out["max_share"][t] = float(w.max() / total) if total > 0 else 0.0
        out["top_1pct_share"][t] = top_p_share(w, 0.01)
        out["top_10pct_share"][t] = top_p_share(w, 0.10)
        out["mean_wealth"][t] = float(w.mean())
        alpha, _, _ = hill_tail_alpha(w)
        out["alpha_hill"][t] = alpha
    return out
