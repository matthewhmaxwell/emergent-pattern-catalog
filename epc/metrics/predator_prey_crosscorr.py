"""P11 — Predator-prey oscillation cross-correlation metrics.

Detects the signature of a bilateral predator-prey relationship in two-
species population time series:

  - Strong anti-correlation between species densities at nonzero lag
    (predator density lags prey density by a quarter period of the
    underlying oscillation; at that lag rho is near zero; at ~half
    a period later rho reaches its minimum value).
  - Oscillation in each species' density (FFT peak-to-mean ratio
    significantly above the noise floor of ~2-3).

Reference physics: Mobilia-Georgiev-Täuber (2007), Täuber (2024) — in
finite lattice LV systems, resonant amplification of demographic noise
produces persistent erratic oscillations. Global density amplitudes
shrink as O(1/√L²) with system size, but the FFT peak-to-mean ratio
remains robustly > 10 across the coexistence phase.

Primary metric: rho_anti = min_tau in excluded-zero-lag band of
    Pearson(species_A(t), species_B(t+tau))

Secondary metrics:
  - tau_anti: lag at which rho_anti is achieved
  - fft_peak_to_mean: predator density FFT peak-to-mean (skip DC)
  - rho_at_zero_lag: instantaneous anti-correlation (for contrast)

Empirically measured ranges (Sprint 11 characterization):
  LV (coexistence phase):  rho_anti in [-0.90, -0.73],
                            |tau_anti| in [11, 18],
                            fft_peak_to_mean in [18, 32]
  RPS species pairs:        rho_anti in [-0.96, -0.93]  <- fires if
                            gated only on rho_anti! Must use n_species
                            prerequisite at detector level.
  SIR S-vs-I:               rho_anti in [-0.37, -0.35]  <- borderline
  Nowak-May (c vs d):       rho_anti = -0.98 at lag +3  <- TRAP
                            caused by strict A+B = const conservation.
                            Detector must check A+B variance > 0.
  Schelling:                rho_anti = 0.00 (zero species variance;
                            agent identity conserved per-agent)
  White noise:              rho_anti in [-0.09, -0.07]
"""

from __future__ import annotations

import numpy as np


def cross_correlation_lag_range(
    x: np.ndarray,
    y: np.ndarray,
    max_lag: int,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute Pearson cross-correlation between x and y over lags
    from -max_lag to +max_lag, inclusive.

    Convention: positive lag k means y lags x — i.e. we correlate
    x[0:n-k] with y[k:n].

    Parameters
    ----------
    x, y : np.ndarray
        1D arrays of equal length.
    max_lag : int
        Maximum absolute lag to evaluate. The output has length
        2*max_lag + 1 with lag axis np.arange(-max_lag, max_lag+1).

    Returns
    -------
    lags : np.ndarray
        Integer array of lags from -max_lag to +max_lag.
    rhos : np.ndarray
        Pearson correlations, same length as lags.
    """
    x = np.asarray(x, dtype=np.float64)
    y = np.asarray(y, dtype=np.float64)
    n = len(x)
    if n != len(y):
        raise ValueError(
            f"x and y must have equal length, got {n} vs {len(y)}"
        )
    if max_lag < 0:
        raise ValueError(f"max_lag must be >= 0, got {max_lag}")
    if max_lag >= n:
        max_lag = n - 1

    lags = np.arange(-max_lag, max_lag + 1)
    rhos = np.zeros(2 * max_lag + 1, dtype=np.float64)

    for idx, k in enumerate(lags):
        if k >= 0:
            # Positive lag: y lags x. Correlate x[0:n-k] vs y[k:n].
            xa = x[: n - k]
            yb = y[k:]
        else:
            # Negative lag: x lags y. Correlate x[-k:] vs y[:n+k].
            xa = x[-k:]
            yb = y[: n + k]

        if len(xa) < 3:
            rhos[idx] = 0.0
            continue

        sx = xa.std()
        sy = yb.std()
        if sx < 1e-12 or sy < 1e-12:
            rhos[idx] = 0.0
            continue

        rhos[idx] = float(
            np.mean((xa - xa.mean()) * (yb - yb.mean())) / (sx * sy)
        )

    return lags, rhos


def predator_prey_rho_anti(
    prey: np.ndarray,
    predator: np.ndarray,
    max_lag: int = 100,
    min_abs_lag: int = 5,
) -> dict[str, float]:
    """Compute the anti-correlation signature of a predator-prey pair.

    The LV predator-prey signature is a strong negative Pearson
    correlation between the two species' densities at nonzero lag.
    We search over |lag| >= min_abs_lag to avoid near-instantaneous
    anti-correlation (which would be a conservation artifact in
    strictly-conserved 2-species systems like Nowak-May, where
    prey_frac + predator_frac = 1 exactly, producing rho=-1 at lag 0).

    Parameters
    ----------
    prey : np.ndarray
        1D time series of prey density (or species-A fraction).
    predator : np.ndarray
        1D time series of predator density (or species-B fraction).
    max_lag : int
        Maximum absolute lag to search.
    min_abs_lag : int
        Minimum |lag| to consider. Lags with |tau| < min_abs_lag are
        excluded from the search (avoids conservation-artifact
        near-zero-lag anti-correlation).

    Returns
    -------
    dict with:
      rho_anti         : float, min cross-correlation over the excluded-
                         zero-lag band (most negative value). Expected
                         negative for LV.
      tau_anti         : int, lag at which rho_anti is achieved.
      rho_at_zero_lag  : float, rho at lag = 0 (diagnostic).
      n_samples        : int, length of the time series.
    """
    prey = np.asarray(prey, dtype=np.float64)
    predator = np.asarray(predator, dtype=np.float64)

    if len(prey) < 3 * min_abs_lag + 3:
        return {
            "rho_anti": 0.0,
            "tau_anti": 0,
            "rho_at_zero_lag": 0.0,
            "n_samples": int(len(prey)),
        }

    # Cap max_lag at n//3 to guarantee >= n/3 samples in the correlation
    effective_max_lag = min(max_lag, max(1, len(prey) // 3))

    lags, rhos = cross_correlation_lag_range(prey, predator, effective_max_lag)

    # Lag-zero cross-correlation (diagnostic)
    zero_idx = np.where(lags == 0)[0][0]
    rho_zero = float(rhos[zero_idx])

    # Restrict to |lag| >= min_abs_lag for the anti-correlation search
    mask = np.abs(lags) >= min_abs_lag
    if not np.any(mask):
        return {
            "rho_anti": 0.0,
            "tau_anti": 0,
            "rho_at_zero_lag": rho_zero,
            "n_samples": int(len(prey)),
        }

    rhos_band = rhos[mask]
    lags_band = lags[mask]
    arg_min = int(np.argmin(rhos_band))

    return {
        "rho_anti": float(rhos_band[arg_min]),
        "tau_anti": int(lags_band[arg_min]),
        "rho_at_zero_lag": rho_zero,
        "n_samples": int(len(prey)),
    }


def fft_peak_to_mean(
    x: np.ndarray,
    min_period: int = 3,
) -> dict[str, float]:
    """FFT-based oscillation metric.

    Measures the prominence of the strongest non-DC FFT component
    against the mean of all non-DC components in the band with period
    >= min_period.

    LV (coexistence phase) empirically gives peak-to-mean in [18, 32].
    White noise gives ~2-3. The threshold for "has an oscillation" is
    about 8-10 based on the SIR borderline case.

    Parameters
    ----------
    x : np.ndarray
        1D time series (typically predator density).
    min_period : int
        Exclude FFT frequencies corresponding to period < min_period.
        Prevents aliased high-frequency noise from dominating.

    Returns
    -------
    dict with:
      fft_peak_to_mean : float, ratio of max |FFT| to mean |FFT|
                         in the (DC-excluded) band.
      dominant_period  : float, period (in samples) at the FFT peak.
                         -1.0 if no meaningful peak found.
    """
    x = np.asarray(x, dtype=np.float64)
    n = len(x)
    if n < 8:
        return {"fft_peak_to_mean": 0.0, "dominant_period": -1.0}

    x0 = x - x.mean()
    if x0.std() < 1e-12:
        return {"fft_peak_to_mean": 0.0, "dominant_period": -1.0}

    F = np.abs(np.fft.rfft(x0))
    # Skip DC (k=0). Cap at period >= min_period.
    k_max = n // min_period
    if k_max < 2:
        return {"fft_peak_to_mean": 0.0, "dominant_period": -1.0}

    band = F[1 : k_max + 1]
    band_mean = band.mean()
    if band_mean < 1e-12:
        return {"fft_peak_to_mean": 0.0, "dominant_period": -1.0}

    peak_to_mean = float(band.max() / band_mean)
    k_best = int(np.argmax(band)) + 1  # offset for skipped DC
    dominant_period = float(n / k_best)

    return {
        "fft_peak_to_mean": peak_to_mean,
        "dominant_period": dominant_period,
    }


def species_time_series_variance_check(
    prey: np.ndarray,
    predator: np.ndarray,
) -> dict[str, float]:
    """Prerequisite-style variance checks for P11.

    Three checks are reported:

    1. std(prey) and std(predator) — if either is ~0, there's no
       time variation and the system cannot be oscillating. Catches
       Schelling (identity is conserved per-agent, so species
       fractions never change).

    2. std(prey + predator) — if this is ~0, the two species are
       trivially conserved against each other, i.e. there's no
       third reservoir (like empty sites). In that case any anti-
       correlation is an algebraic artifact of prey = 1 - predator.
       Catches Nowak-May.

    Returns
    -------
    dict with:
      prey_std          : float
      predator_std      : float
      total_std         : float, std(prey + predator)
      prey_mean         : float
      predator_mean     : float
    """
    prey = np.asarray(prey, dtype=np.float64)
    predator = np.asarray(predator, dtype=np.float64)
    total = prey + predator

    return {
        "prey_std": float(prey.std()),
        "predator_std": float(predator.std()),
        "total_std": float(total.std()),
        "prey_mean": float(prey.mean()),
        "predator_mean": float(predator.mean()),
    }


def circular_shift_null(
    prey: np.ndarray,
    predator: np.ndarray,
    n_permutations: int,
    max_lag: int,
    min_abs_lag: int,
    seed: int = 42,
) -> dict[str, float]:
    """Null distribution for rho_anti via circular time-series shift.

    At each permutation, we circularly shift the predator time series
    by a uniformly random offset in [1, n-1] (excluding 0, the
    identity). The prey series is unchanged. The circular shift
    preserves:
      - each series' marginal distribution (exactly),
      - each series' autocorrelation (exactly),
      - the FFT magnitude spectrum of each series (exactly).
    It destroys the phase relationship between the two series.

    For an LV-like system, the observed rho_anti should be much more
    negative than 99th-percentile of the null distribution.

    Parameters
    ----------
    prey, predator : np.ndarray
        1D time series.
    n_permutations : int
        Number of permutations (>= 99 for p < 0.01, >= 199 for p < 0.005).
    max_lag : int
        Max lag for rho_anti computation (same as the observed stat).
    min_abs_lag : int
        Min |lag| band (same as observed stat).
    seed : int
        RNG seed.

    Returns
    -------
    dict with:
      p_value         : (n_less_or_equal + 1) / (n_permutations + 1).
                        Fraction of null rho_anti <= observed (one-sided
                        lower tail since we want rho_anti more negative
                        than null).
      null_mean       : float
      null_std        : float
      observed        : float, observed rho_anti
      n_permutations  : int
    """
    prey = np.asarray(prey, dtype=np.float64)
    predator = np.asarray(predator, dtype=np.float64)
    n = len(prey)

    observed = predator_prey_rho_anti(
        prey, predator, max_lag=max_lag, min_abs_lag=min_abs_lag
    )["rho_anti"]

    if n < 3 * min_abs_lag + 3 or n_permutations < 1:
        return {
            "p_value": 1.0,
            "null_mean": 0.0,
            "null_std": 0.0,
            "observed": observed,
            "n_permutations": 0,
        }

    rng = np.random.default_rng(seed)
    null_rhos = np.zeros(n_permutations, dtype=np.float64)

    # Draw n_permutations distinct-ish shifts in [1, n-1].
    # We allow duplicates (it's a small bias but standard practice).
    shifts = rng.integers(1, n, size=n_permutations)

    for i, shift in enumerate(shifts):
        pred_shifted = np.roll(predator, int(shift))
        result = predator_prey_rho_anti(
            prey, pred_shifted, max_lag=max_lag, min_abs_lag=min_abs_lag
        )
        null_rhos[i] = result["rho_anti"]

    # One-sided lower tail: count nulls <= observed (more negative or equal)
    n_le = int(np.sum(null_rhos <= observed))
    p_value = (n_le + 1) / (n_permutations + 1)

    return {
        "p_value": float(p_value),
        "null_mean": float(null_rhos.mean()),
        "null_std": float(null_rhos.std()),
        "observed": float(observed),
        "n_permutations": int(n_permutations),
    }


# ---------------------------------------------------------------------------
# Time series extraction from state history
# ---------------------------------------------------------------------------

def extract_species_fractions(
    state_history: list[dict],
    prey_state: int = 1,
    predator_state: int = 2,
) -> tuple[np.ndarray, np.ndarray]:
    """Extract two species' fraction time series from a state history.

    Works for any lattice model where state dicts contain a 'grid' key
    with integer-encoded species labels. For Lotka-Volterra lattice
    (the canonical positive), prey=1 and predator=2. For generic
    two-species lattice models, caller should specify the correct
    state integers.

    Parameters
    ----------
    state_history : list[dict]
        List of state snapshots, each with a 'grid' key (2D integer
        np.ndarray).
    prey_state : int
        Integer label of the prey species on the grid.
    predator_state : int
        Integer label of the predator species on the grid.

    Returns
    -------
    prey_fraction, predator_fraction : np.ndarray
        1D arrays of length len(state_history) giving each species'
        fraction of total cells at each timestep.
    """
    if not state_history:
        return np.zeros(0), np.zeros(0)
    if "grid" not in state_history[0]:
        raise ValueError(
            "state_history entries must contain 'grid' key for P11 metric"
        )

    prey = np.zeros(len(state_history), dtype=np.float64)
    pred = np.zeros(len(state_history), dtype=np.float64)
    for i, s in enumerate(state_history):
        g = s["grid"]
        total = g.size
        prey[i] = float((g == prey_state).sum()) / total
        pred[i] = float((g == predator_state).sum()) / total

    return prey, pred
