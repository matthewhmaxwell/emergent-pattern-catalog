"""Metrics for the P8 traffic-jamming detector.

All metrics operate on a sequence of NS state snapshots. The primary
metric is the stopped-fraction <1[v_i(t) = 0]>_t post-burn-in. Secondary
metrics are jam-lifetime distribution statistics (per-car consecutive-v=0
runs), gap-distribution statistics, and a temporal-shuffle null.

References:
  Bette, H. M., Habel, L., Emig, T. & Schreckenberg, M. (2017).
  "Mechanisms of jamming in the Nagel-Schreckenberg model for traffic
  flow." Phys. Rev. E 95, 012311. — Uses P(v=0) as order parameter;
  decomposes jammed fraction into (jamming rate) x (jam lifetime) x
  (jam size).
"""

from __future__ import annotations

from typing import Any

import numpy as np


def collect_velocity_history(
    state_history: list[dict[str, Any]],
    burn_in: int = 0,
) -> np.ndarray | None:
    """Stack per-step velocity arrays into a (T, n_cars) integer matrix.

    Expects every snapshot past ``burn_in`` to expose a 'velocities' key
    with identical length. Returns None if the observable is missing or
    shape-inconsistent (e.g., n_cars varies across snapshots — which
    should never happen in NS but could in a misconfigured caller).
    """
    if not state_history:
        return None
    snapshots = state_history[burn_in:]
    if not snapshots:
        return None
    per_step: list[np.ndarray] = []
    first_N: int | None = None
    for snap in snapshots:
        v = snap.get("velocities")
        if v is None:
            return None
        arr = np.asarray(v)
        if arr.ndim != 1:
            return None
        if first_N is None:
            first_N = arr.shape[0]
        elif arr.shape[0] != first_N:
            return None
        per_step.append(arr)
    if not per_step:
        return None
    return np.stack(per_step, axis=0)


def collect_gap_history(
    state_history: list[dict[str, Any]],
    burn_in: int = 0,
) -> np.ndarray | None:
    """Stack per-step gap arrays into a (T, n_cars) integer matrix."""
    if not state_history:
        return None
    snapshots = state_history[burn_in:]
    if not snapshots:
        return None
    per_step: list[np.ndarray] = []
    first_N: int | None = None
    for snap in snapshots:
        g = snap.get("gaps")
        if g is None:
            return None
        arr = np.asarray(g)
        if arr.ndim != 1:
            return None
        if first_N is None:
            first_N = arr.shape[0]
        elif arr.shape[0] != first_N:
            return None
        per_step.append(arr)
    if not per_step:
        return None
    return np.stack(per_step, axis=0)


def stopped_fraction(velocity_history: np.ndarray) -> float:
    """Time-and-car-averaged fraction of cars with v = 0.

    Primary metric for P8. Bette et al. (2017) use this as the order
    parameter for the NS jamming transition.
    """
    if velocity_history.size == 0:
        return 0.0
    return float((velocity_history == 0).mean())


def mean_velocity(velocity_history: np.ndarray) -> float:
    """Time-and-car-averaged mean velocity."""
    if velocity_history.size == 0:
        return 0.0
    return float(velocity_history.mean())


def flow(velocity_history: np.ndarray, density: float) -> float:
    """Flow = density * mean_velocity."""
    return density * mean_velocity(velocity_history)


def jam_lifetime_stats(velocity_history: np.ndarray) -> dict[str, float]:
    """Per-car consecutive-v=0 runs: distribution statistics.

    For each car, a "jam event" is a maximal consecutive interval where
    that car's velocity equals zero. We collect all such intervals
    across all cars and report distribution statistics.

    Returns
    -------
    dict with keys:
      n_jam_events : int   total number of stopped-runs
      mean         : float mean run length
      median       : float median run length
      p95          : float 95th percentile
      max          : int   longest run
    All zero if no car is ever stopped.
    """
    if velocity_history.size == 0:
        return {"n_jam_events": 0, "mean": 0.0, "median": 0.0, "p95": 0.0, "max": 0}
    T, N = velocity_history.shape
    is_stopped = velocity_history == 0
    lifetimes: list[int] = []
    for car in range(N):
        s = is_stopped[:, car]
        if not s.any():
            continue
        # Run-length encoding of True segments
        # Pad with False at both ends to catch edge runs
        padded = np.concatenate(([False], s, [False]))
        diffs = np.diff(padded.astype(np.int8))
        starts = np.where(diffs == 1)[0]
        ends = np.where(diffs == -1)[0]
        for a, b in zip(starts, ends):
            lifetimes.append(int(b - a))
    if not lifetimes:
        return {"n_jam_events": 0, "mean": 0.0, "median": 0.0, "p95": 0.0, "max": 0}
    arr = np.asarray(lifetimes)
    return {
        "n_jam_events": int(arr.size),
        "mean": float(arr.mean()),
        "median": float(np.median(arr)),
        "p95": float(np.percentile(arr, 95)),
        "max": int(arr.max()),
    }


def gap_cv(gap_history: np.ndarray) -> float:
    """Coefficient of variation of the gap distribution.

    Pooled across all cars and all timesteps. In free flow, gaps are
    near-uniform (low CV ~ 0.15-0.6). In jammed regimes gaps become
    bimodal (high CV ~ 0.8-1.4).
    """
    if gap_history.size == 0:
        return 0.0
    mu = float(gap_history.mean())
    if mu <= 0.0:
        return 0.0
    return float(gap_history.std() / mu)


def zero_gap_fraction(gap_history: np.ndarray) -> float:
    """Fraction of (car, time) samples with gap = 0 (bumper-to-bumper)."""
    if gap_history.size == 0:
        return 0.0
    return float((gap_history == 0).mean())


def temporal_shuffle_null_stopped(
    velocity_history: np.ndarray,
    n_permutations: int,
    seed: int,
) -> np.ndarray:
    """Null distribution for stopped_fraction by independent per-car
    temporal shuffle.

    WARNING: this null preserves the marginal velocity distribution per
    car AND the global stopped-fraction. It is useful as a null for
    secondary metrics that are sensitive to TEMPORAL CORRELATION
    (e.g., jam-lifetime distribution), not for the stopped-fraction
    primary itself. For the primary, the relevant null is a counterfactual
    low-density / deterministic run, which is handled separately in the
    detector's tier logic.

    For symmetry with other detectors' ``*_null_distribution`` metrics,
    we still expose this function; its output is expected to be tightly
    concentrated around the observed stopped_fraction.
    """
    if velocity_history.size == 0:
        return np.zeros(n_permutations, dtype=np.float64)
    rng = np.random.default_rng(seed)
    T, N = velocity_history.shape
    nulls = np.zeros(n_permutations, dtype=np.float64)
    for k in range(n_permutations):
        shuffled = velocity_history.copy()
        # Shuffle each car's v-trace independently
        for car in range(N):
            rng.shuffle(shuffled[:, car])
        nulls[k] = float((shuffled == 0).mean())
    return nulls


def temporal_shuffle_null_jam_p95(
    velocity_history: np.ndarray,
    n_permutations: int,
    seed: int,
) -> np.ndarray:
    """Null distribution for jam-lifetime p95 under independent per-car
    temporal shuffle.

    This is the detector's primary null. The shuffle preserves the per-
    car marginal distribution of v (so stopped_fraction is unchanged),
    but destroys the persistence structure that makes jam-lifetimes
    heavy-tailed. Under the null, lifetimes are geometric with parameter
    equal to the stopped marginal, so p95 is small (typically < 5 cells
    at NS canonical parameters).

    An observed p95 >> null-distribution p95 is the signature of genuine
    NS jamming as opposed to independent per-car v-fluctuations.
    """
    if velocity_history.size == 0:
        return np.zeros(n_permutations, dtype=np.float64)
    rng = np.random.default_rng(seed)
    T, N = velocity_history.shape
    nulls = np.zeros(n_permutations, dtype=np.float64)
    for k in range(n_permutations):
        shuffled = velocity_history.copy()
        for car in range(N):
            rng.shuffle(shuffled[:, car])
        stats = jam_lifetime_stats(shuffled)
        nulls[k] = stats["p95"]
    return nulls
