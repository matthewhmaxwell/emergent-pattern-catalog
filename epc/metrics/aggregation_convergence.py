"""
Aggregation Convergence Metric — Temporal guard for P1 false positives.

Problem: GoL produces genuine spatial autocorrelation (Moran's I is significant
vs label-shuffle null) because B3/S23 dynamics create live-cell clusters. But
there is no attraction mechanism — it's rule-driven structural clustering. The
P1 detector fires as a conceptual false positive.

Solution: Real similarity-driven aggregation (Schelling, Zhang sorting) shows
monotonic convergence of the aggregation index toward a steady state. GoL and
other rule-driven systems show oscillatory or chaotic aggregation trajectories.

This metric computes:
1. Spearman ρ(I(t), t) — trend direction and strength
2. Coefficient of variation of I(t) in the measurement window — stability
3. Plateau detection — whether I(t) reaches and holds a steady state

Usage: Add as a secondary check in the P1 detector. If convergence fails,
demote P1 to "structural autocorrelation" and do not confirm.

Architecture decision: This is a temporal signature test, not a model-metadata
check. It works on state-history only, preserving P1's "state-history only"
observable scope.
"""

from __future__ import annotations

import numpy as np
from scipy import stats
from dataclasses import dataclass
from typing import List, Optional


@dataclass
class ConvergenceResult:
    """Result of aggregation convergence analysis."""
    spearman_rho: float          # Trend direction: +1 = monotonic increase
    spearman_p: float            # Significance of trend
    cv_late: float               # CV of I(t) in last 20% — low = plateau
    is_monotonic: bool           # ρ > threshold AND p < 0.01
    is_plateaued: bool           # CV in late window < threshold
    is_oscillatory: bool         # High variance, no trend
    has_gain: bool               # I_late > I_initial + threshold (convergence from below)
    i_initial: float             # I at first snapshot
    i_early: float               # Mean I in first 20%
    i_late: float                # Mean I in last 20%
    convergence_score: float     # Combined score [0, 1]
    n_windows: int               # Number of time windows used
    i_trajectory: np.ndarray     # The Moran's I time series


def compute_morans_i_2d(grid: np.ndarray) -> float:
    """Compute Moran's I for a binary 2D grid with queen contiguity.
    
    Parameters
    ----------
    grid : np.ndarray
        2D binary array (0/1 states).
    
    Returns
    -------
    float
        Moran's I value. Range approx [-1, 1].
        Positive = spatial clustering of like values.
    """
    rows, cols = grid.shape
    x = grid.astype(float).ravel()
    n = len(x)
    x_bar = x.mean()
    
    if np.std(x) < 1e-12:
        return 0.0  # Uniform grid — no autocorrelation
    
    # Build queen-contiguity weight sum using vectorized shifts
    # For each cell, sum of (x_i - x_bar)(x_j - x_bar) over neighbors
    padded = np.pad(grid.astype(float), 1, mode='wrap')  # periodic BC
    
    deviations = grid.astype(float) - x_bar
    
    # 8 neighbors (queen contiguity)
    neighbor_sum = np.zeros_like(deviations)
    for di in [-1, 0, 1]:
        for dj in [-1, 0, 1]:
            if di == 0 and dj == 0:
                continue
            shifted = padded[1+di:rows+1+di, 1+dj:cols+1+dj]
            neighbor_sum += (shifted - x_bar)
    
    numerator = np.sum(deviations * neighbor_sum)
    denominator = np.sum(deviations ** 2)
    
    if abs(denominator) < 1e-12:
        return 0.0
    
    # W = total weight (8 neighbors per cell for periodic BC interior)
    W = n * 8  # Each cell has 8 neighbors with periodic BC
    
    I = (n / W) * (numerator / denominator)
    return float(I)


def compute_morans_i_1d(array: np.ndarray, types: np.ndarray) -> float:
    """Compute Moran's I for a 1D array with type labels.
    
    Parameters
    ----------
    array : np.ndarray
        1D array of agent positions (indices 0..N-1).
    types : np.ndarray
        1D array of type labels (integer).
    
    Returns
    -------
    float
        Moran's I.
    """
    n = len(types)
    x = types.astype(float)
    x_bar = x.mean()
    
    if np.std(x) < 1e-12:
        return 0.0
    
    deviations = x - x_bar
    denominator = np.sum(deviations ** 2)
    
    if abs(denominator) < 1e-12:
        return 0.0
    
    # Adjacent neighbors in 1D
    numerator = np.sum(deviations[:-1] * deviations[1:]) * 2  # symmetric
    W = 2 * (n - 1)  # Each pair counted twice
    
    I = (n / W) * (numerator / denominator)
    return float(I)


def compute_convergence(
    history: List[dict],
    grid_key: str = 'grid',
    types_key: Optional[str] = None,
    n_windows: int = 20,
    rho_threshold: float = 0.3,
    cv_plateau_threshold: float = 0.10,
    late_fraction: float = 0.2,
) -> ConvergenceResult:
    """Compute temporal convergence of aggregation metric.
    
    Divides the state history into n_windows temporal bins, computes
    Moran's I for the last snapshot in each bin, then analyzes the
    resulting time series for monotonic convergence vs oscillation.
    
    Parameters
    ----------
    history : list of dict
        State history. Each dict must contain either:
        - grid_key: 2D binary array (for lattice models like GoL, GH), OR
        - types_key: 1D type labels (for 1D models like sorting)
    grid_key : str
        Key for 2D grid in state dict.
    types_key : str or None
        Key for 1D type labels. If provided, uses 1D Moran's I.
    n_windows : int
        Number of temporal windows to sample. More = finer resolution.
    rho_threshold : float
        Minimum Spearman ρ for monotonic convergence.
    cv_plateau_threshold : float
        Maximum CV in late window for plateau detection.
    late_fraction : float
        Fraction of trajectory considered "late" for plateau test.
    
    Returns
    -------
    ConvergenceResult
    """
    T = len(history)
    if T < n_windows * 2:
        n_windows = max(5, T // 2)
    
    # Sample one snapshot per window, ALWAYS including the very first timestep
    window_size = T // n_windows
    indices = [0]  # Always start with the first timestep
    indices += [min((i + 1) * window_size - 1, T - 1) for i in range(n_windows)]
    # Remove duplicates and sort
    indices = sorted(set(indices))
    n_samples = len(indices)
    
    i_values = np.zeros(n_samples)
    for w, idx in enumerate(indices):
        state = history[idx]
        if types_key is not None and types_key in state:
            i_values[w] = compute_morans_i_1d(
                np.arange(len(state[types_key])), state[types_key]
            )
        elif grid_key in state:
            i_values[w] = compute_morans_i_2d(state[grid_key])
        else:
            raise KeyError(
                f"State dict must contain '{grid_key}' or '{types_key}'. "
                f"Available keys: {list(state.keys())}"
            )
    
    # 1. Spearman trend test
    time_indices = np.arange(n_samples)
    rho, p_val = stats.spearmanr(time_indices, i_values)
    
    # 2. Late-window CV (plateau test)
    late_start = int(n_samples * (1 - late_fraction))
    late_values = i_values[late_start:]
    late_mean = np.mean(late_values)
    if abs(late_mean) > 1e-12:
        cv_late = float(np.std(late_values) / abs(late_mean))
    else:
        cv_late = float('inf')
    
    # 3. Early-vs-late I comparison (convergence from below)
    # Use the FIRST snapshot's I, not the mean of early windows.
    # This catches fast convergence (e.g., Schelling settles in 10 steps
    # out of 150 — the first 20% of windows already shows the plateau).
    i_initial = float(i_values[0])     # Very first sampled I
    early_end = int(n_samples * late_fraction)
    early_end = max(early_end, 1)
    early_values = i_values[:early_end]
    i_early = float(np.mean(early_values))  # Still compute for reporting
    i_late = float(late_mean)
    gain_threshold = 0.10  # I must increase by at least 0.10
    # Use i_initial (not i_early) for the gain check
    has_gain = (i_late - i_initial) > gain_threshold
    
    # 4. Classification
    is_monotonic = (rho > rho_threshold) and (p_val < 0.01)
    is_plateaued = cv_late < cv_plateau_threshold
    
    # Oscillatory: no clear trend AND high variance
    overall_cv = float(np.std(i_values) / max(abs(np.mean(i_values)), 1e-12))
    is_oscillatory = (abs(rho) < rho_threshold) and (overall_cv > 0.15)
    
    # 5. Convergence score
    trend_component = np.clip(rho, 0, 1) * 0.33
    plateau_component = np.clip(1.0 - cv_late / cv_plateau_threshold, 0, 1) * 0.33
    gain_component = np.clip((i_late - i_initial) / 0.3, 0, 1) * 0.34
    convergence_score = float(trend_component + plateau_component + gain_component)
    
    return ConvergenceResult(
        spearman_rho=float(rho),
        spearman_p=float(p_val),
        cv_late=cv_late,
        is_monotonic=is_monotonic,
        is_plateaued=is_plateaued,
        is_oscillatory=is_oscillatory,
        has_gain=has_gain,
        i_initial=i_initial,
        i_early=i_early,
        i_late=i_late,
        convergence_score=convergence_score,
        n_windows=n_samples,
        i_trajectory=i_values,
    )


def p1_temporal_guard(
    history: List[dict],
    grid_key: str = 'grid',
    types_key: Optional[str] = None,
    n_windows: int = 20,
    types_are_constant: bool = True,
) -> tuple[bool, ConvergenceResult]:
    """Check whether P1 detection should be allowed based on temporal signature.
    
    Returns (passes, result) where passes=True means the aggregation shows
    convergent dynamics consistent with genuine similarity-driven aggregation.
    
    P1 requires TWO conditions:
    1. TYPE CONSTANCY: The "type" labels must be persistent agent properties
       that don't change over time. If grid states change every step (like GoL
       alive/dead), they are NOT types — they are dynamic states. This is
       checked via the types_are_constant flag.
    2. CONVERGENCE FROM BELOW: Moran's I must increase from the initial
       condition, indicating that spatial sorting occurred.
    
    Parameters
    ----------
    history : list of dict
        Full state history.
    grid_key : str
        Key for 2D grid data.
    types_key : str or None
        Key for 1D type labels.
    n_windows : int
        Temporal resolution.
    types_are_constant : bool
        Whether the type labels are constant over time. If False (e.g., GoL
        alive/dead states change every step), P1 is structurally inapplicable.
        Default True for backward compatibility with sorting/Schelling.
    
    Returns
    -------
    passes : bool
        True if temporal signature is consistent with P1.
    result : ConvergenceResult
        Full convergence analysis.
    """
    result = compute_convergence(
        history, grid_key=grid_key, types_key=types_key, n_windows=n_windows
    )
    
    # Condition 1: Type constancy
    # P1 requires constant type labels (agent identity). If types change
    # every step (like GoL cell states), this is NOT similarity-driven
    # aggregation — it's structural autocorrelation from dynamics.
    if not types_are_constant:
        return False, result
    
    # Condition 2: Convergence from below
    # I must increase from the initial condition to the final state.
    # Systems that start with structural autocorrelation (no initial random
    # mixing) don't show this gain.
    passes = result.has_gain and (result.is_monotonic or result.is_plateaued)
    
    return passes, result
