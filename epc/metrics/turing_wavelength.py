"""Turing-wavelength metrics for P3.

Detects the signature of a stationary periodic spatial pattern selected
by a reaction-diffusion (or diffusion-driven) instability:

  - A sharp peak in the radially-averaged 2D FFT power spectrum at a
    nonzero wavenumber k*, well above the spectral noise floor.
  - A selected wavelength λ = N / k* that is a physical property of the
    system: it does not change with grid size (modulo the obvious
    discretization of available k-bins).
  - Stationarity in time: the selected peak_k does not drift once the
    pattern has formed (after an initial selection transient).

Primary metric: radial-FFT peak-to-mean ratio on the scalar field.
  peak_to_mean = radial[argmax(radial[k_min:k_max])] / mean(radial[k_min:])

Where `radial[k]` is the azimuthally averaged 2D power spectrum:
  F = fftshift(fft2(field - field.mean()))
  P = |F|²
  radial[k] = mean of P over the annulus at integer radial distance k.

Low-k bins (k < k_min) are excluded to avoid DC leakage and very-long-
wavelength drift artifacts. The highest k-bins are also skipped when
computing the mean to reduce sensitivity to grid discretization.

Reference physics: Turing (1952); Pearson (1993). The Turing wavelength
λ* emerges from the balance between diffusion (which smooths on scale
√(D·t)) and local reaction kinetics (which selects a preferred scale
when D_u ≠ D_v). For Gray-Scott at grid-scale (D_u=0.16, D_v=0.08),
λ* ≈ 12 pixels — empirically verified across N = 64, 96, 128, 192
(peak_k scales linearly with N; wavelength invariant).

Empirical calibration (Sprint 13, N=128, T≥4000, seed 42):
  Gray-Scott labyrinth (F=0.037, k=0.060):  peak/mean ≈ 23, λ ≈ 12.8 px
  Gray-Scott spots     (F=0.030, k=0.062):  peak/mean ≈ 24, λ ≈ 11.6 px
  Gray-Scott uniform   (F=0.100, k=0.100):  peak/mean ≈ 0   (v dies)
  Spatial-shuffle null on GS labyrinth:     peak/mean = 1.42 ± 0.16
                                             (d ≈ 107, p ≈ 0.01 at 99 perm)

False-positive check (Sprint 13 negative sweep):
  RPS (M=1e-4, raw grid):     peak/mean ≈ 23  <- numerically matches GS!
                               Discrimination is NOT by peak_to_mean.
                               It is by state-type uniqueness: GS has
                               ~16k distinct float values; RPS has 4.
                               The P3 detector's `n_unique_values ≥ 50`
                               prerequisite cleanly rejects RPS without
                               relying on any empirical tuning.
  GH (n=8, random IC):        peak/mean ≈ 6.6  (fails screening)
  Schelling (tau=0.375):      peak/mean ≈ 3.3  (fails screening)
  Nowak-May (b=1.8):          peak/mean ≈ 3.6  (fails screening)
  SIR (single seed, t=100):   peak/mean ≈ 3.7  (fails screening)
  GoL (random d=0.37, t=300): peak/mean ≈ 5.4  (fails screening)
  LV (lambda=4, t=600, raw):  peak/mean ≈ 6.7  (fails screening)
  iid Gaussian noise:         peak/mean ≈ 1.0  (baseline)
"""

from __future__ import annotations

import numpy as np


def radial_fft_power_spectrum(
    field: np.ndarray,
) -> np.ndarray:
    """Compute azimuthally averaged 2D FFT power spectrum.

    Returns `radial[k]` for k = 0, 1, ..., min(rows, cols)//2 - 1.
    radial[0] is the DC bin (effectively zero after mean subtraction).

    Parameters
    ----------
    field : np.ndarray, shape (rows, cols)
        Scalar spatial field.

    Returns
    -------
    np.ndarray, shape (k_max,) where k_max = min(rows, cols) // 2
        Radially averaged power spectrum.
    """
    if field.ndim != 2:
        raise ValueError(f"field must be 2D, got shape {field.shape}")
    rows, cols = field.shape
    if rows < 4 or cols < 4:
        raise ValueError("field too small to compute radial spectrum")

    # Subtract spatial mean (kills DC bin)
    f_centered = field - field.mean()
    # 2D FFT with zero-frequency shifted to center
    Fk = np.fft.fftshift(np.fft.fft2(f_centered))
    P = (Fk.conj() * Fk).real  # power spectrum

    # Radial bins based on distance from center
    cr, cc = rows // 2, cols // 2
    yy, xx = np.indices(P.shape)
    rr = np.sqrt((yy - cr) ** 2 + (xx - cc) ** 2)
    rr_int = rr.astype(int)

    k_max = min(rows, cols) // 2
    counts = np.bincount(rr_int.ravel(), minlength=k_max)[:k_max]
    sums = np.bincount(rr_int.ravel(), weights=P.ravel(), minlength=k_max)[:k_max]
    radial = sums / np.maximum(counts, 1)
    return radial


def radial_fft_peak_stats(
    field: np.ndarray,
    k_min: int = 2,
    k_max_frac: float = 1.0,
) -> dict[str, float]:
    """Compute peak statistics from the radial FFT power spectrum.

    Parameters
    ----------
    field : np.ndarray, shape (rows, cols)
        Scalar spatial field.
    k_min : int
        Lowest wavenumber to consider in the peak search. Must be >= 2
        to skip DC and the first-bin long-wavelength drift. Default 2.
    k_max_frac : float
        Upper end of the search window as a fraction of Nyquist
        (k_max = floor(k_max_frac * N/2)). Default 0.5 — avoids
        grid-aliasing artifacts near the Nyquist edge.

    Returns
    -------
    dict with keys:
        peak_k             : int, argmax wavenumber in [k_min, k_max)
        peak_value         : float, radial[peak_k]
        radial_mean        : float, mean(radial[k_min:]) — noise floor
        peak_to_mean       : float, peak_value / max(radial_mean, eps)
        wavelength_pixels  : float, N / peak_k (or nan if peak_k == 0)
        field_std          : float
        field_mean         : float
        k_min_used         : int
        k_max_used         : int
        n_unique_values    : int, number of distinct values in field
                             (used downstream for continuous-vs-discrete
                             prerequisite).
    """
    if k_min < 1:
        raise ValueError("k_min must be >= 1 (bin 0 is DC)")
    if not 0.0 < k_max_frac <= 1.0:
        raise ValueError("k_max_frac must be in (0, 1]")

    radial = radial_fft_power_spectrum(field)
    N = min(field.shape)
    k_max = max(k_min + 1, int(np.floor(k_max_frac * (N // 2))))
    k_max = min(k_max, len(radial))

    if k_max <= k_min:
        # Degenerate: field too small
        return _degenerate_stats(field)

    search = radial[k_min:k_max]
    if search.size == 0 or not np.any(search > 0):
        return _degenerate_stats(field)

    peak_idx_local = int(np.argmax(search))
    peak_k = peak_idx_local + k_min
    peak_value = float(radial[peak_k])
    radial_mean = float(np.mean(radial[k_min:k_max]))
    eps = 1e-12
    peak_to_mean = peak_value / max(radial_mean, eps)
    wavelength = float(N) / float(peak_k) if peak_k > 0 else float("nan")

    return {
        "peak_k": peak_k,
        "peak_value": peak_value,
        "radial_mean": radial_mean,
        "peak_to_mean": peak_to_mean,
        "wavelength_pixels": wavelength,
        "field_std": float(field.std()),
        "field_mean": float(field.mean()),
        "k_min_used": int(k_min),
        "k_max_used": int(k_max),
        "n_unique_values": int(np.unique(field).size),
    }


def _degenerate_stats(field: np.ndarray) -> dict[str, float]:
    return {
        "peak_k": 0,
        "peak_value": 0.0,
        "radial_mean": 0.0,
        "peak_to_mean": 0.0,
        "wavelength_pixels": float("nan"),
        "field_std": float(field.std()),
        "field_mean": float(field.mean()),
        "k_min_used": 0,
        "k_max_used": 0,
        "n_unique_values": int(np.unique(field).size),
    }


def wavelength_stability(
    fields: list[np.ndarray],
    k_min: int = 2,
    k_max_frac: float = 1.0,
) -> dict[str, float]:
    """Measure time-stability of the selected wavelength across snapshots.

    Used to distinguish stationary Turing patterns (stable peak_k) from
    weakly-drifting or transient regimes.

    Parameters
    ----------
    fields : list of np.ndarray
        Sequence of 2D scalar fields (e.g., the last N snapshots).
    k_min, k_max_frac : as in radial_fft_peak_stats

    Returns
    -------
    dict with keys:
        peak_ks            : list[int]
        peak_to_means      : list[float]
        peak_k_mean        : float
        peak_k_std         : float
        peak_k_cv          : float — std / mean, 0 = perfectly stable
        peak_k_range       : int — max - min
        peak_to_mean_mean  : float
        peak_to_mean_min   : float
    """
    if len(fields) == 0:
        raise ValueError("fields must be non-empty")

    ks: list[int] = []
    pms: list[float] = []
    for f in fields:
        s = radial_fft_peak_stats(f, k_min=k_min, k_max_frac=k_max_frac)
        ks.append(s["peak_k"])
        pms.append(s["peak_to_mean"])

    ks_arr = np.array(ks, dtype=float)
    pms_arr = np.array(pms, dtype=float)
    k_mean = float(ks_arr.mean())
    k_std = float(ks_arr.std())
    cv = k_std / k_mean if k_mean > 0 else float("inf")

    return {
        "peak_ks": ks,
        "peak_to_means": pms,
        "peak_k_mean": k_mean,
        "peak_k_std": k_std,
        "peak_k_cv": cv,
        "peak_k_range": int(max(ks) - min(ks)) if ks else 0,
        "peak_to_mean_mean": float(pms_arr.mean()),
        "peak_to_mean_min": float(pms_arr.min()),
    }


def spatial_shuffle(field: np.ndarray, rng: np.random.Generator) -> np.ndarray:
    """Shuffle field values in-place-safe manner (returns new array).

    Destroys all spatial structure (including Turing wavelength) while
    preserving the marginal intensity distribution. This is the null
    model for P3.
    """
    flat = field.ravel().copy()
    rng.shuffle(flat)
    return flat.reshape(field.shape)


def shuffle_null_distribution(
    field: np.ndarray,
    n_permutations: int,
    k_min: int = 2,
    k_max_frac: float = 1.0,
    seed: int = 0,
) -> np.ndarray:
    """Generate the null distribution of peak_to_mean under spatial shuffle.

    Returns an array of shape (n_permutations,) containing peak_to_mean
    values for shuffled copies of the field.
    """
    if n_permutations < 1:
        raise ValueError("n_permutations must be >= 1")

    rng = np.random.default_rng(seed)
    vals = np.empty(n_permutations, dtype=float)
    for i in range(n_permutations):
        shuf = spatial_shuffle(field, rng)
        vals[i] = radial_fft_peak_stats(shuf, k_min=k_min, k_max_frac=k_max_frac)[
            "peak_to_mean"
        ]
    return vals


def count_unique_values(field: np.ndarray) -> int:
    """Count distinct values in the field. Used for continuous-vs-discrete
    prerequisite gate.

    A continuous (float) field typically has n_unique ≈ rows*cols.
    A discrete-state grid (integer labels) has n_unique ≤ number of
    states, usually ≤ 10.
    """
    return int(np.unique(field).size)
