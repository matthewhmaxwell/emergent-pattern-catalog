"""
Optimized boundary-conditioned Transfer Entropy (global aggregate).

Sprint 5 implementation: vectorized global boundary TE matching the semantics
of p13_p15_discriminator._boundary_te but with np.add.at batch counting.
3.5× faster than per-cell approach at 60×60.

The key insight: boundary-conditioned TE restricts measurement to cells
at state-transition boundaries (where adjacent cells have different states).
This isolates the information-transfer signal at structure interactions.

Raw average TE gives WRONG ordering (GH > GoL) because GH's deterministic
wave fronts create trivially high TE everywhere. Boundary conditioning
correctly gives GoL >> GH (15-16× ratio at 60×60).

Sprint 5 benchmark results (60×60, 300 steps, 99 perms):
  GH spiral (κ=5): TE = 0.000628  →  ratio vs GH control = 1.0×  → P13
  GH random (κ=3): TE = 0.000444  →  ratio vs GH control = 0.7×  → P13
  GoL random:      TE = 0.009511  →  ratio vs GH control = 15.1× → P15_candidate
  GoL R-pentomino: TE = 0.010131  →  ratio vs GH control = 16.1× → P15_candidate

Architecture decision #16 (boundary-conditioned TE).
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class GlobalBoundaryTEResult:
    """Result of global boundary-conditioned TE computation."""
    te_mean: float              # Global boundary TE (bits)
    null_mean: float            # Mean TE under temporal permutation null
    null_std: float             # Std of null TE
    ratio_vs_null: float        # te_mean / null_mean
    p_value: float              # Permutation test p-value
    n_boundary_obs: int         # Total boundary cell observations
    n_permutations: int         # Number of permutations


def compute_global_boundary_te(
    grids: np.ndarray,
    n_states: int = 2,
    n_permutations: int = 99,
    seed: int = 42,
) -> GlobalBoundaryTEResult:
    """Compute boundary-conditioned Transfer Entropy using global aggregate counts.

    Accumulates frequency tables across ALL boundary cells and ALL timesteps,
    then computes a single global TE. This matches the semantics of
    P13P15Discriminator._boundary_te but uses np.add.at for batch counting.

    Parameters
    ----------
    grids : (T, rows, cols) uint8 array
        Time series of grid states.
    n_states : int
        Number of discrete states.
    n_permutations : int
        Permutations for significance test.
    seed : int
        RNG seed.

    Returns
    -------
    GlobalBoundaryTEResult
    """
    T, rows, cols = grids.shape
    grids_int = grids.astype(np.intp)

    te_obs, n_obs = _global_te_pass(grids_int, grids_int, n_states)

    if n_obs == 0:
        return GlobalBoundaryTEResult(
            te_mean=0.0, null_mean=0.0, null_std=0.0,
            ratio_vs_null=0.0, p_value=1.0,
            n_boundary_obs=0, n_permutations=n_permutations,
        )

    # Permutation null: shuffle time axis of source grids
    rng = np.random.default_rng(seed)
    null_tes = np.zeros(n_permutations)

    for perm_idx in range(n_permutations):
        perm = rng.permutation(T)
        grids_shuffled = grids_int[perm]
        null_te, _ = _global_te_pass(grids_int, grids_shuffled, n_states)
        null_tes[perm_idx] = null_te

    null_mean = float(np.mean(null_tes))
    null_std = float(np.std(null_tes))
    ratio = te_obs / null_mean if null_mean > 0 else 0.0
    p_value = float(np.mean(null_tes >= te_obs))
    if p_value == 0:
        p_value = 1.0 / (n_permutations + 1)

    return GlobalBoundaryTEResult(
        te_mean=te_obs,
        null_mean=null_mean,
        null_std=null_std,
        ratio_vs_null=ratio,
        p_value=p_value,
        n_boundary_obs=n_obs,
        n_permutations=n_permutations,
    )


def _global_te_pass(
    grids_target: np.ndarray,
    grids_source: np.ndarray,
    n_states: int,
) -> tuple[float, int]:
    """Compute global boundary TE in one pass over the trajectory.

    Parameters
    ----------
    grids_target : (T, rows, cols) int array
        Target grid (boundary detection + target variable).
    grids_source : (T, rows, cols) int array
        Source grid (neighbor values for TE). May be temporally shuffled.
    n_states : int
        Number of discrete states.

    Returns
    -------
    (te, n_observations) : float, int
    """
    T, rows, cols = grids_target.shape
    dirs = [(0, 1), (1, 0), (0, -1), (-1, 0)]  # VN neighbors

    # Global frequency tables per direction
    joint = [np.zeros((n_states, n_states, n_states), dtype=np.float64) for _ in range(4)]
    yt_xt = [np.zeros((n_states, n_states), dtype=np.float64) for _ in range(4)]
    ynext_yt = np.zeros((n_states, n_states), dtype=np.float64)
    yt_count = np.zeros(n_states, dtype=np.float64)
    total = 0

    for t in range(T - 1):
        g = grids_target[t]
        g_next = grids_target[t + 1]
        g_src = grids_source[t]

        # Vectorized boundary detection (Moore neighborhood, periodic)
        is_boundary = np.zeros((rows, cols), dtype=bool)
        for dr in (-1, 0, 1):
            for dc in (-1, 0, 1):
                if dr == 0 and dc == 0:
                    continue
                shifted = np.roll(np.roll(g, -dr, axis=0), -dc, axis=1)
                is_boundary |= (g != shifted)

        br, bc = np.where(is_boundary)
        if len(br) == 0:
            continue

        y_vals = g[br, bc]
        y_next_vals = g_next[br, bc]

        # Accumulate marginal counts
        np.add.at(ynext_yt, (y_next_vals, y_vals), 1)
        np.add.at(yt_count, (y_vals,), 1)
        total += len(br)

        # VN neighbor values and joint counts per direction
        for d_idx, (dr, dc) in enumerate(dirs):
            x_vals = np.roll(np.roll(g_src, -dr, axis=0), -dc, axis=1)[br, bc]
            np.add.at(joint[d_idx], (y_next_vals, y_vals, x_vals), 1)
            np.add.at(yt_xt[d_idx], (y_vals, x_vals), 1)

    if total == 0:
        return 0.0, 0

    # H(Y_next | Y_curr) for boundary cells
    h_y_yt = 0.0
    for yn in range(n_states):
        for y in range(n_states):
            n = ynext_yt[yn, y]
            if n <= 0:
                continue
            p = n / total
            p_cond = n / yt_count[y] if yt_count[y] > 0 else 0
            if p_cond > 0:
                h_y_yt -= p * np.log2(p_cond)

    # TE per direction, then average
    te_dirs = []
    for d_idx in range(4):
        h_y_yt_xt = 0.0
        for yn in range(n_states):
            for y in range(n_states):
                for x in range(n_states):
                    n = joint[d_idx][yn, y, x]
                    if n <= 0:
                        continue
                    n_yx = yt_xt[d_idx][y, x]
                    if n_yx <= 0:
                        continue
                    p = n / total
                    p_cond = n / n_yx
                    if p_cond > 0:
                        h_y_yt_xt -= p * np.log2(p_cond)
        te = max(0.0, h_y_yt - h_y_yt_xt)
        te_dirs.append(te)

    return float(np.mean(te_dirs)), total
