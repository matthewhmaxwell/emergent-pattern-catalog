"""Metrics for chimera-state detection (P10).

Chimera states are coexisting coherent and incoherent subpopulations in a
system of phase oscillators with non-local coupling. The ring topology of
the Abrams–Strogatz (2004) / Kuramoto–Battogtokh (2002) model organizes
phases by position; a chimera has a spatially-contiguous coherent arc
(local order parameter ~1) and a spatially-contiguous incoherent arc
(local order parameter ~0.5).

DISCRIMINATION CHALLENGE

The NAIVE chimera signatures — "gap between max and min local r", "local
r CV", "persistence of window ordering" — all give false positives on
ordinary all-to-all Kuramoto near its critical coupling K_c. Phase 1 of
Sprint 18 quantified this:

  Random uniform phases N=128: gap = 0.54 ± 0.08
  Chimera (canonical):          gap = 0.52
  Kuramoto K=Kc=1.0 (partial):  gap = 0.68, "coexistence" = True
  Chimera persistence_corr:     +0.55 to +0.68
  Kuramoto K=2.0 persistence:   +0.93 (higher than chimera!)

The reason: Kuramoto all-to-all with heterogeneous ω produces persistent
partial sync — oscillators near the center of ω entrain, tails drift — and
after the model's internal ω-sort, this looks window-wise indistinguishable
from a chimera's arc structure. Any window-level static signature fails.

THE CLEAN DISCRIMINATOR: `pos_vel_ac`

Time-average per-oscillator phase velocity, then compute its SPATIAL
autocorrelation at lag k on the ring. A chimera has neighbors on the ring
that drift at similar effective rates (both locked in the coherent arc, or
both driven by the common non-local coupling field in the incoherent arc);
this autocorrelation is high (+0.85 to +0.95). Ordinary all-to-all Kuramoto
has each oscillator drifting at its own ω_i (randomly assigned to ring
positions by the ω-sort), so neighbor velocities are not spatially
correlated beyond the uniform mean-field term; this autocorrelation is
low (+0.1 to +0.45).

Phase 1j confirmed the separation at lag 4 on N=128 rings:
  Chimera (coex-passing, multi-seed β=0.18): 0.929 ± 0.007
  Kuramoto K=1.0, 6 seeds:                   0.312 ± 0.130
  Separation gap: chim min 0.919 - Kur max 0.448 = +0.47

See ADR 50 in REPLICATION_NOTES.md.

PRIMARY METRIC: `pos_vel_ac(vel_matrix, lag=4)`

SECONDARY GATES:
- `coexistence_gate`: windows persistently coherent AND persistently
  incoherent. Keeps full-sync and full-incoherent at screening rejection.
- `persistence_corr`: temporal correlation of local_r vectors (used as
  a diagnostic; NOT the primary).
- `local_r_time_matrix`: (T_rec, n_windows) cached for null.

NULL:
- SURROGATE null: circular-permute oscillator indices. This destroys
  ring adjacency while preserving the marginal distribution of phases
  per frame. Under the null, pos_vel_ac collapses to ~0.
- n_permutations >= 199 for CONFIRMATION (floor p = 0.005).
"""

from __future__ import annotations

from typing import Any

import numpy as np
from numpy.typing import NDArray


# ----------------------------------------------------------------------
# Local order parameter and window helpers
# ----------------------------------------------------------------------


def compute_local_r(theta: NDArray[np.float64], n_windows: int = 16) -> NDArray[np.float64]:
    """Local Kuramoto order parameter in contiguous position-based windows.

    theta : (N,) array of phases on a ring.
    n_windows : number of contiguous arc windows.
    Returns (n_windows,) array of local |<exp(iθ)>|.

    For a ring of N oscillators the windows tile the index axis in equal
    chunks; the last window absorbs the remainder if N is not divisible.
    """
    N = len(theta)
    win = max(1, N // n_windows)
    out = np.empty(n_windows, dtype=np.float64)
    for k in range(n_windows):
        start = k * win
        stop = (k + 1) * win if k < n_windows - 1 else N
        chunk = theta[start:stop]
        if len(chunk) >= 2:
            out[k] = float(np.abs(np.mean(np.exp(1j * chunk))))
        else:
            out[k] = 1.0  # singleton window is trivially coherent
    return out


def local_r_time_matrix(
    theta_hist: NDArray[np.float64], n_windows: int = 16
) -> NDArray[np.float64]:
    """(T_rec, n_windows) matrix of local r over time.

    theta_hist : (T_rec, N) array of phase snapshots.
    """
    T = len(theta_hist)
    out = np.empty((T, n_windows), dtype=np.float64)
    for t in range(T):
        out[t] = compute_local_r(theta_hist[t], n_windows=n_windows)
    return out


# ----------------------------------------------------------------------
# Coexistence gate (screening prerequisite)
# ----------------------------------------------------------------------


def coexistence_stats(
    lr_matrix: NDArray[np.float64],
    coh_thresh: float = 0.85,
    incoh_thresh: float = 0.70,
    persistence_fraction: float = 0.90,
) -> dict[str, Any]:
    """Test for coherent/incoherent coexistence in local-r over time.

    Two complementary measures:

    1. ``per_frame_coexistence_fraction``: fraction of frames in which
       AT LEAST ONE window exceeds ``coh_thresh`` AND AT LEAST ONE
       window is below ``incoh_thresh`` in the SAME frame. A drifting
       chimera satisfies this at ~100% of frames even though the
       coherent arc moves around the ring.

    2. ``n_persistent_coh`` / ``n_persistent_incoh``: number of windows
       that spend at least ``persistence_fraction`` of frames above /
       below the thresholds. These measure spatial stationarity of the
       chimera arc. For a stationary chimera (short runs), both are
       >= 1. For a drifting chimera (long runs), both may be 0 despite
       genuine coexistence.

    The ``coexistence`` flag uses the per-frame measure (drift-invariant).
    The persistence counts are kept as diagnostics.

    Returns dict with
        n_persistent_coh : int (spatial stationarity diagnostic)
        n_persistent_incoh : int
        per_frame_coexistence_fraction : float in [0, 1]
        coexistence : bool (per-frame fraction >= persistence_fraction)
        gap_timeavg : float (time-averaged local r max - min)
        local_r_min : float
        local_r_max : float
        above_frac_max : float
        below_frac_max : float
    """
    T, W = lr_matrix.shape
    above = np.mean(lr_matrix > coh_thresh, axis=0)
    below = np.mean(lr_matrix < incoh_thresh, axis=0)
    n_coh = int(np.sum(above >= persistence_fraction))
    n_incoh = int(np.sum(below >= persistence_fraction))

    # Drift-invariant per-frame coexistence
    any_coh_per_frame = np.any(lr_matrix > coh_thresh, axis=1)
    any_incoh_per_frame = np.any(lr_matrix < incoh_thresh, axis=1)
    per_frame_coexistence = np.logical_and(any_coh_per_frame,
                                            any_incoh_per_frame)
    per_frame_coex_frac = float(np.mean(per_frame_coexistence))

    time_mean = np.mean(lr_matrix, axis=0)
    return {
        "n_persistent_coh": n_coh,
        "n_persistent_incoh": n_incoh,
        "per_frame_coexistence_fraction": per_frame_coex_frac,
        "coexistence": (per_frame_coex_frac >= persistence_fraction),
        "gap_timeavg": float(time_mean.max() - time_mean.min()),
        "local_r_min": float(time_mean.min()),
        "local_r_max": float(time_mean.max()),
        "above_frac_max": float(above.max()),
        "below_frac_max": float(below.max()),
    }


# ----------------------------------------------------------------------
# Primary: position-velocity spatial autocorrelation
# ----------------------------------------------------------------------


def phase_velocities(
    theta_hist: NDArray[np.float64], record_dt: float = 1.0
) -> NDArray[np.float64]:
    """Instantaneous per-oscillator phase velocity.

    Unwraps phases along the time axis and finite-differences.
    Returns (T_rec - 1, N) array of dθ/dt.

    theta_hist must have T_rec >= 2.
    """
    if len(theta_hist) < 2:
        raise ValueError(
            f"phase_velocities: need >= 2 frames, got {len(theta_hist)}"
        )
    theta_unwrap = np.unwrap(theta_hist, axis=0)
    dt = float(record_dt) if record_dt > 0 else 1.0
    return np.diff(theta_unwrap, axis=0) / dt


def pos_vel_ac(
    vel: NDArray[np.float64], lag: int = 4
) -> float:
    """Spatial autocorrelation of time-averaged phase velocity on a ring.

    vel : (T, N) per-oscillator velocities.
    lag : ring-index lag at which to compute spatial autocorrelation.

    Returns a scalar in [-1, 1].

    Mechanism: a chimera has spatially-correlated velocity structure
    (neighbors drift together because the non-local coupling kernel mixes
    them); ordinary Kuramoto has each oscillator drifting at its own ω,
    so spatial correlation is weak.

    DEGENERATE CASE: if std(time-averaged velocity) == 0 (full sync, every
    oscillator at identical mean frequency), the autocorrelation is
    undefined and we return 1.0. This does NOT promote full-sync to
    chimera because the coexistence gate rejects full-sync at screening.
    """
    v_mean = np.mean(vel, axis=0)
    N = len(v_mean)
    if N <= lag:
        return 0.0
    centered = v_mean - np.mean(v_mean)
    std = float(np.std(v_mean))
    if std < 1e-9:
        return 1.0  # degenerate, handled by coexistence gate upstream
    centered = centered / std
    # Circular autocorrelation at lag `lag`
    return float(np.mean(centered * np.roll(centered, -lag)))


# ----------------------------------------------------------------------
# Surrogate null (position-shuffle)
# ----------------------------------------------------------------------


def shuffle_null_pos_vel_ac(
    theta_hist: NDArray[np.float64],
    record_dt: float,
    n_permutations: int,
    lag: int,
    rng: np.random.Generator,
) -> NDArray[np.float64]:
    """Surrogate null distribution for pos_vel_ac.

    For each permutation, randomly relabel oscillator indices (permute the
    second axis of theta_hist), recompute velocities, and recompute
    pos_vel_ac. This destroys ring adjacency while preserving
    per-frame distributions of phases and per-oscillator trajectories.

    Under the null (no spatial structure), pos_vel_ac is ~0 with small
    variance set by finite N.

    Returns (n_permutations,) array of null pos_vel_ac values.
    """
    T, N = theta_hist.shape
    out = np.empty(n_permutations, dtype=np.float64)
    for i in range(n_permutations):
        perm = rng.permutation(N)
        theta_perm = theta_hist[:, perm]
        vel_perm = phase_velocities(theta_perm, record_dt=record_dt)
        out[i] = pos_vel_ac(vel_perm, lag=lag)
    return out


# ----------------------------------------------------------------------
# Secondary diagnostics
# ----------------------------------------------------------------------


def persistence_corr(lr_matrix: NDArray[np.float64]) -> float:
    """Mean Pearson correlation of consecutive local_r vectors.

    Diagnostic only; NOT a primary (Kuramoto near K_c also scores high).
    """
    T = len(lr_matrix)
    if T < 2:
        return 0.0
    corrs: list[float] = []
    for t in range(T - 1):
        a, b = lr_matrix[t], lr_matrix[t + 1]
        if np.std(a) > 1e-6 and np.std(b) > 1e-6:
            corrs.append(float(np.corrcoef(a, b)[0, 1]))
    return float(np.mean(corrs)) if corrs else 0.0


def global_order_parameter(
    theta_hist: NDArray[np.float64]
) -> tuple[float, float]:
    """Mean and std of the global Kuramoto order parameter over time."""
    z = np.mean(np.exp(1j * theta_hist), axis=1)
    r = np.abs(z)
    return float(np.mean(r)), float(np.std(r))


def velocity_spatial_neighbor_corr(
    vel: NDArray[np.float64]
) -> float:
    """Mean Pearson correlation of velocity time series between ring neighbors.

    Secondary diagnostic: a chimera's incoherent arc has neighbors driven
    by the common coupling field, so their velocity TIME SERIES are
    correlated (even if the time-averaged values aren't). Complements
    pos_vel_ac which uses only the time-averaged spatial structure.
    """
    T, N = vel.shape
    if T < 3:
        return 0.0
    corrs: list[float] = []
    for i in range(N - 1):
        a, b = vel[:, i], vel[:, i + 1]
        if np.std(a) > 1e-6 and np.std(b) > 1e-6:
            corrs.append(float(np.corrcoef(a, b)[0, 1]))
    # Wrap around the ring
    a, b = vel[:, -1], vel[:, 0]
    if np.std(a) > 1e-6 and np.std(b) > 1e-6:
        corrs.append(float(np.corrcoef(a, b)[0, 1]))
    return float(np.mean(corrs)) if corrs else 0.0
