"""Reproduction of Reichenbach, Mobilia & Frey (2007) Fig 2c: spiral wavelength λ ~ M^(1/2).

Figure identity: λ(M) log-log scaling in the coexistence regime, from Fig. 2c of
Reichenbach, T., Mobilia, M. & Frey, E. (2007). "Mobility promotes and jeopardizes
biodiversity in rock-paper-scissors games." Nature 448(7157), 1046–1049.
arXiv: q-bio/0702032 (open access).

Sprint history:
  Sprint 54 (2026-05-26): L=100, narrow M sweep [3e-4..5e-4], ACF estimator.
    slope=0.366, outside [0.4,0.6]. Root cause: M range too narrow.
  Sprint 58 (2026-05-29): L=100, wide M sweep [1e-5..5e-4], ACF estimator.
    slope=0.107. Root cause: formula breaks down at M ≪ M_c.
  Sprint 59 (2026-05-29): L=100, dense near-M_c sweep [2e-4..5e-4], ACF estimator.
    Fit restricted to formula-valid sub-range. Result carried forward.
  Sprint 63 (2026-05-30): L=200, FFT structure-factor ring-peak estimator.
    Target: resolve low-M spirals that exceeded L=100 measurement window.
    M values: 5 points in [5e-5, 5e-4] (wider range enabled by L=200).

Wavelength estimator — FFT structure-factor ring peak (Sprint 63):
  Compute 2D power spectrum P(k) = |FFT(field)|² of species-A density field.
  Azimuthally average to get S(k) vs |k|. Find dominant non-zero-k peak.
  λ = L / k_peak (in lattice units, where k is in units of 2π/L).
  Validated on synthetic sinusoidal field of known wavelength before use.

  Motivation for FFT-ring over ACF at L=200: at larger L, the radial ACF
  first-zero method becomes noisy at large lag distances. The FFT ring-peak
  directly measures the dominant spatial frequency and is more robust when
  λ is a large fraction of L.

M_c ≈ 4.5e-4 (Reichenbach 2007, L-independent in continuum limit, σ=μ=1).
  λ = 0.8·L·√(M/M_c): at L=200, wavelengths are 2× larger in lattice units
  than at L=100, giving low-M spirals room to fit within the measurement window.

Published result: λ ∝ M^(1/2) (Reichenbach 2007 Eq. 2 / Fig. 2c).
Tolerance: log-log fit slope within [0.4, 0.6], R² ≥ 0.90.
"""
from __future__ import annotations

import json
import os
import sys
import time
from multiprocessing import Pool

import numpy as np

sys.path.insert(0, ".")
from epc.models.rps_spatial import RPSSpatialModel

# --- Sprint 63 parameters: L=200, FFT-ring estimator ---
L: int = 200
# M_c ≈ 4.5e-4 (L-independent in continuum limit, σ=μ=1; Reichenbach 2007)
M_C: float = 4.5e-4

# 5 M values spanning ~1 decade in coexistence regime
# At L=200, λ = 0.8·200·√(M/M_c) ranges from ~53 (M=5e-5) to ~134 (M=5e-4)
# All fit within L=200 measurement window (λ/L ≤ 0.67)
M_VALUES: list[float] = [
    5e-5,    # M/M_c = 0.111, λ_formula ≈ 53
    1e-4,    # M/M_c = 0.222, λ_formula ≈ 75
    2e-4,    # M/M_c = 0.444, λ_formula ≈ 107
    3.5e-4,  # M/M_c = 0.778, λ_formula ≈ 141
    5e-4,    # M/M_c = 1.111, λ_formula ≈ 169 — at/above M_c, may extinct
]
N_SEEDS: int = 15
T_EQ: int = 2500       # longer equilibration for L=200 spirals
T_MEASURE: int = 300   # measurement window (generations)
MEASURE_EVERY: int = 30  # snapshot stride → 10 snapshots per run
N_SNAPSHOTS: int = T_MEASURE // MEASURE_EVERY  # = 10


def fft_ring_wavelength(field: np.ndarray, pad_factor: int = 4) -> float:
    """Estimate spiral wavelength via zero-padded FFT structure-factor ring peak.

    Zero-pads the field to (pad_factor × L)² before computing the 2D power
    spectrum, giving k-resolution of 1/pad_factor pixel. Azimuthally averages
    S(k) and finds the dominant non-zero-k peak with parabolic interpolation.
    Returns λ = L / k_peak.

    Parameters
    ----------
    field : np.ndarray
        2D real-valued density field (L×L).
    pad_factor : int
        Zero-padding multiplier (default 4 → k resolution of 0.25 pixel).

    Returns
    -------
    float
        Estimated wavelength in lattice units, or np.nan.
    """
    Lf = field.shape[0]
    field_c = field - field.mean()

    # Zero-pad for finer k resolution
    Lp = Lf * pad_factor
    padded = np.zeros((Lp, Lp))
    padded[:Lf, :Lf] = field_c

    # 2D power spectrum of zero-padded field
    fft_f = np.fft.fft2(padded)
    power = np.abs(fft_f) ** 2
    power = np.fft.fftshift(power)
    cy, cx = Lp // 2, Lp // 2

    # Radial distances (in units of original k pixels, not padded pixels)
    # One original k-unit = pad_factor padded pixels
    y, x = np.ogrid[-cy:Lp - cy, -cx:Lp - cx]
    r = np.sqrt(x.astype(float) ** 2 + y.astype(float) ** 2) / pad_factor

    # Azimuthal average in fine radial bins (width = 1/pad_factor original k-units)
    max_k = Lf // 2  # Nyquist in original k-units
    bin_width = 1.0 / pad_factor
    n_bins = int(max_k / bin_width)
    k_centers = (np.arange(n_bins) + 0.5) * bin_width
    s_k = np.zeros(n_bins)
    for bi in range(n_bins):
        k_lo = bi * bin_width
        k_hi = (bi + 1) * bin_width
        mask = (r >= k_lo) & (r < k_hi)
        if np.any(mask):
            s_k[bi] = float(np.mean(power[mask]))

    # Find dominant peak for k >= 0.75 (skip DC neighborhood)
    start_idx = max(1, int(0.75 / bin_width))
    if start_idx >= n_bins or np.all(s_k[start_idx:] == 0):
        return np.nan

    peak_idx = int(np.argmax(s_k[start_idx:]) + start_idx)
    k_peak = k_centers[peak_idx]

    if k_peak >= max_k - 1 or s_k[peak_idx] <= 0:
        return np.nan

    # Parabolic sub-bin interpolation
    if 1 <= peak_idx < n_bins - 1 and s_k[peak_idx - 1] > 0 and s_k[peak_idx + 1] > 0:
        y_m = np.log(s_k[peak_idx - 1])
        y_0 = np.log(s_k[peak_idx])
        y_p = np.log(s_k[peak_idx + 1])
        denom = 2.0 * (2.0 * y_0 - y_m - y_p)
        if abs(denom) > 1e-12:
            delta = (y_m - y_p) / denom
            k_peak = k_centers[peak_idx] + delta * bin_width

    # λ = L / k_peak (in original lattice units)
    wavelength = float(Lf) / float(k_peak)
    return wavelength


def validate_fft_estimator() -> None:
    """Validate FFT-ring estimator on synthetic fields with known wavelength.

    Uses ring-pattern fields (sum of x- and y-axis sinusoids) that produce
    an isotropic ring in k-space, matching the spectral structure of spirals.
    """
    print("=== FFT-ring estimator validation ===")
    test_L = 200
    rng = np.random.default_rng(12345)
    for true_lambda in [20.0, 40.0, 80.0, 120.0]:
        # Ring pattern: sin(kx) + sin(ky) + noise → isotropic ring at |k|=L/λ
        x = np.arange(test_L)
        y = np.arange(test_L)
        X, Y = np.meshgrid(x, y)
        k = 2.0 * np.pi / true_lambda
        field = np.sin(k * X) + np.sin(k * Y) + 0.1 * rng.standard_normal((test_L, test_L))

        measured = fft_ring_wavelength(field)
        if np.isnan(measured):
            print(f"  λ_true={true_lambda:6.1f}: measured=NaN  FAIL")
        else:
            rel_err = abs(measured - true_lambda) / true_lambda
            status = "OK" if rel_err < 0.15 else "WARN"
            print(f"  λ_true={true_lambda:6.1f}: measured={measured:6.1f}  "
                  f"rel_err={rel_err*100:.1f}%  {status}")
    print("=== End validation ===\n")


def run_one_seed(M: float, seed: int) -> float:
    """Run one (M, seed) trial; return mean wavelength over measurement window."""
    model = RPSSpatialModel(
        rows=L,
        cols=L,
        mobility=M,
        selection_rate=1.0,
        reproduction_rate=1.0,
        neighborhood="von_neumann",
        seed=seed,
    )
    model.run(T_EQ, record_every=T_EQ)
    full_history = model.run(T_MEASURE, record_every=MEASURE_EVERY)
    snapshots = full_history[2:]  # skip initial state and early transient

    lambdas: list[float] = []
    for snap in snapshots:
        grid: np.ndarray = snap["grid"]
        field = (grid == RPSSpatialModel.SPECIES_A).astype(float)
        lam = fft_ring_wavelength(field)
        lambdas.append(lam)

    valid = [l for l in lambdas if not np.isnan(l)]
    if not valid:
        return np.nan
    return float(np.mean(valid))


def _worker(args: tuple) -> tuple:
    """Multiprocessing worker: run one (M_index, seed_index) trial.

    Seed convention: seed = 6300 + M_index * 100 + seed_index (Sprint 63 offset).
    """
    M_index, seed_index = args
    M = M_VALUES[M_index]
    seed = 6300 + M_index * 100 + seed_index
    lam = run_one_seed(M, seed)
    return M_index, seed_index, lam


def main() -> None:
    """Run L=200 FFT-ring λ(M) sweep and save JSON artifact."""
    validate_fft_estimator()

    t_start = time.time()

    print(f"Sprint 63 — L={L} FFT-ring wavelength sweep")
    print(f"M_c(L={L}) = {M_C:.4e}")
    print(f"M values: {M_VALUES}")
    print(f"N_SEEDS={N_SEEDS}, T_EQ={T_EQ}, T_MEASURE={T_MEASURE}")
    print()

    # Dispatch 5×15 = 75 independent runs in parallel
    tasks = [
        (M_index, seed_index)
        for M_index in range(len(M_VALUES))
        for seed_index in range(N_SEEDS)
    ]
    n_workers = min(os.cpu_count() or 1, 8)
    print(f"Running {len(tasks)} simulations on {n_workers} workers ...\n")

    wavelengths: list[list[float]] = [[np.nan] * N_SEEDS for _ in range(len(M_VALUES))]
    with Pool(processes=n_workers) as pool:
        for M_index, seed_index, lam in pool.map(_worker, tasks):
            wavelengths[M_index][seed_index] = lam

    # Per-M statistics
    per_point: list[dict] = []
    for M_index, M in enumerate(M_VALUES):
        seeds = wavelengths[M_index]
        valid = [l for l in seeds if not np.isnan(l)]
        n_valid = len(valid)
        lam_mean = float(np.mean(valid)) if valid else float("nan")
        lam_std = float(np.std(valid, ddof=1)) if len(valid) >= 2 else float("nan")
        lam_sem = lam_std / float(np.sqrt(n_valid)) if n_valid >= 2 else float("nan")
        per_point.append({
            "M": M,
            "lambda_mean": round(lam_mean, 4) if not np.isnan(lam_mean) else None,
            "lambda_sem": round(lam_sem, 4) if not np.isnan(lam_sem) else None,
            "n_valid": n_valid,
            "n_seeds": N_SEEDS,
        })

    # Log-log fit on all points with n_valid >= 5
    fit_points = [
        p for p in per_point
        if p["lambda_mean"] is not None and p["n_valid"] >= 5
    ]
    n_fit_points = len(fit_points)

    if n_fit_points < 3:
        slope = None
        intercept = None
        r_squared = None
        tolerance_pass = False
        r_squared_pass = False
        overall_pass = False
    else:
        log_M = np.array([np.log(p["M"]) for p in fit_points])
        log_lam = np.array([np.log(p["lambda_mean"]) for p in fit_points])
        slope_arr, intercept_arr = np.polyfit(log_M, log_lam, 1)
        slope = float(slope_arr)
        intercept = float(intercept_arr)

        y_pred = slope * log_M + intercept
        ss_res = float(np.sum((log_lam - y_pred) ** 2))
        ss_tot = float(np.sum((log_lam - log_lam.mean()) ** 2))
        r_squared = float(1.0 - ss_res / ss_tot) if ss_tot > 0 else 0.0

        tolerance_pass = bool(0.40 <= slope <= 0.60)
        r_squared_pass = bool(r_squared >= 0.90)
        overall_pass = bool(tolerance_pass and r_squared_pass)

    elapsed = time.time() - t_start

    # Console summary
    for p in per_point:
        M = p["M"]
        lam_m = f"{p['lambda_mean']:7.2f}" if p["lambda_mean"] is not None else "    NaN"
        lam_s = f"{p['lambda_sem']:5.2f}" if p["lambda_sem"] is not None else "  NaN"
        nv = f"{p['n_valid']}/{N_SEEDS}"
        print(f"M={M:.2e}: λ={lam_m} ± {lam_s}  n_valid={nv}")

    print()
    print(f"Fit points used: {n_fit_points}")
    if slope is not None:
        pass_str = "PASS" if tolerance_pass else "FAIL"
        r2_str = "PASS" if r_squared_pass else "FAIL"
        print(f"Log-log slope: {slope:.4f}  (target 0.500, band [0.40, 0.60])  {pass_str}")
        print(f"R²:            {r_squared:.4f}  (target ≥ 0.90)                     {r2_str}")
    else:
        print("Log-log slope: None  (insufficient fit points)")
        print("R²:            None")
    ov_str = "PASS" if overall_pass else "FAIL"
    print(f"Overall:       {ov_str}")

    # JSON output
    out: dict = {
        "sprint": 63,
        "L": L,
        "estimator": "fft_ring_peak",
        "m_c": M_C,
        "m_c_note": f"M_c = {M_C:.4e} (L-independent, Reichenbach 2007, sigma=mu=1)",
        "M_values": M_VALUES,
        "n_fit_points": n_fit_points,
        "per_point": per_point,
        "log_log_slope": round(slope, 6) if slope is not None else None,
        "log_log_intercept": round(intercept, 6) if intercept is not None else None,
        "r_squared": round(r_squared, 6) if r_squared is not None else None,
        "tolerance_pass": tolerance_pass,
        "r_squared_pass": r_squared_pass,
        "overall_pass": overall_pass,
        "n_seeds": N_SEEDS,
        "T_eq": T_EQ,
        "T_measure": T_MEASURE,
        "elapsed_seconds": round(elapsed, 1),
        "sprint_history": [
            {"sprint": 54, "L": 100, "estimator": "acf_first_zero", "slope": 0.366, "note": "M range too narrow"},
            {"sprint": 58, "L": 100, "estimator": "acf_first_zero", "slope": 0.107, "note": "formula breaks down at M << M_c"},
            {"sprint": 59, "L": 100, "estimator": "acf_first_zero", "slope": None, "note": "dense near-M_c sweep"},
        ],
    }

    os.makedirs("analysis/outputs", exist_ok=True)
    output_path = "analysis/outputs/p12_reichenbach2007_reproduction.json"
    with open(output_path, "w") as f:
        json.dump(out, f, indent=2)

    print(f"\nOutput: {output_path}  [{elapsed:.0f}s]")


if __name__ == "__main__":
    main()
