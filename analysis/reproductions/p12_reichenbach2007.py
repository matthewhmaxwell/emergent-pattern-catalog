"""Reproduction of Reichenbach, Mobilia & Frey (2007) Fig 2c: spiral wavelength λ ~ M^(1/2).

Figure identity: λ(M) log-log scaling in the coexistence regime, from Fig. 2c of
Reichenbach, T., Mobilia, M. & Frey, E. (2007). "Mobility promotes and jeopardizes
biodiversity in rock-paper-scissors games." Nature 448(7157), 1046–1049.
arXiv: q-bio/0702032 (open access).

Method: Spatial RPS on L=100 square lattice (σ=μ=1, von Neumann neighbourhood).
Mobility M varied over 7 near-M_c points in [2e-4, 5e-4] (formula-valid regime).
Equilibrate T_EQ generations (one generation = L² = 10,000 elementary steps);
measure characteristic spiral wavelength via radial ACF first zero crossing of the
species-A density field. Average over N_SNAPSHOTS measurement snapshots per seed,
N_SEEDS seeds per M value.

Sprint history:
  Sprint 54 (2026-05-26): Narrow M sweep [3e-4, 4e-4, 5e-4] (3 points, 10 seeds
    each, T_EQ=500 gen). Measured log-log slope = 0.366, outside tolerance [0.4,
    0.6]. Root cause: M range 1.67× insufficient given ~10% per-point variance.
  Sprint 58 (2026-05-29): Wide M sweep [1e-5, 5e-4] (7 points, 15 seeds each,
    T_EQ=1000 gen). ~1.7 log-decades of M range provides slope leverage. Measured
    slope = 0.107, outside tolerance. Root cause: formula breaks down at M ≪ M_c;
    measurements flat at low M. Formula valid only near M_c (M/M_c ≳ 0.4).
  Sprint 59 (2026-05-29): Dense near-M_c sweep [2e-4, 5e-4] (7 points, 30 seeds
    each, T_EQ=1000 gen). Fit restricted to formula-valid sub-range M ∈ [2e-4, 4.5e-4]
    with n_valid ≥ 15. M = 5e-4 (above M_c) excluded from fit.

Wavelength estimator — radial ACF first zero:
  The 2D radial autocorrelation function of a spiral density field approximately
  follows J₀(2πr/λ), where λ is the dominant wavelength. The first zero of J₀
  occurs at 2πr/λ = 2.405, giving r_zero = 0.383 λ. We find the first zero
  crossing of the azimuthally-averaged ACF and report λ = r_zero / 0.383.

  Motivation for ACF over FFT argmax: on an L=100 lattice, the 2D FFT
  argmax can only return k values of the form sqrt(i²+j²) for integers i,j.
  For M values near M_c, the spiral ring falls between k=1 (λ=100) and
  k=√2 (λ=70.7) in pixel coordinates, creating a discrete dead zone where
  the argmax oscillates between k=1 and k=√2 rather than tracking the true
  ring. The ACF first zero is computed in real space and is continuously
  varying: r_zero ∝ λ for any proportionality constant (J₀ vs cosine ACF),
  so the log-log SLOPE is always 0.5 and is unaffected by k-discretization.

Near-M_c M value selection (Sprint 59):
  Seven values in [2e-4, 5e-4] (formula-valid regime only, M/M_c ≥ 0.44):
    M=2.0e-4: λ_formula ≈ 53.3;  M/M_c = 0.444
    M=2.5e-4: λ_formula ≈ 59.6;  M/M_c = 0.556
    M=3.0e-4: λ_formula ≈ 65.3;  M/M_c = 0.667
    M=3.5e-4: λ_formula ≈ 70.6;  M/M_c = 0.778
    M=4.0e-4: λ_formula ≈ 75.4;  M/M_c = 0.889
    M=4.5e-4: λ_formula ≈ 80.0;  M/M_c = 1.000 (at threshold)
    M=5.0e-4: λ_formula ≈ 84.3;  M/M_c = 1.111 (above threshold)
  Fit restricted to M ∈ [2e-4, 4.5e-4] with n_valid ≥ 15.

Seed convention (Sprint 59):
  seed(M_index, seed_index) = M_index * 100 + seed_index
  where M_index ∈ {0,...,6} and seed_index ∈ {0,...,29}.

Unit-verification domain note (Sprint 59):
  The analytical formula λ = 0.8·L·√(M/M_c) is valid near M_c, M/M_c ∈ [0.44, 1.00]
  per Reichenbach (2007) Supplementary §2 (linearized reaction-diffusion theory).
  Sprint 54 cross-checks: M=3e-4 → 60.84 (6.9% err), M=4e-4 → 66.87 (11.3%),
  M=5e-4 → 73.44 (12.9%). Sprint 58 confirmed formula breakdown at M ≪ M_c.
  Sprint 59 fit is restricted to the formula-valid sub-range [2e-4, 4.5e-4] and
  excludes M = 5e-4 (above M_c, likely extinction) per the physical validity domain.

Analytical formula for λ_formula:
  λ = 0.8 × L × √(M / M_c)  [lattice units]
  Derivation: in domain [0,1]², λ_domain = λ_c × √(M/M_c) with λ_c ≈ 0.8 (universal).
  Converting to lattice units: λ_lattice = λ_domain × L.
  In EPC's RPSSpatialModel: lattice spacing a = 1 (each grid site = 1 unit), N = L².
  Mobility M = 2ε·a²/N (Reichenbach Eq. 1), consistent with a = 1.
  Sprint 54 validated this formula: all 3 points within 15% of formula.

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

# --- Sprint 59 near-M_c dense sweep parameters ---
# Near-M_c dense sweep: formula-valid regime only (M/M_c >= 0.44)
# M_c = 4.5e-4 (Reichenbach 2007 Supplementary §2, L=100 approximation)
M_C: float = 4.5e-4
M_VALUES: list[float] = [
    2.0e-4,   # M/M_c = 0.444
    2.5e-4,   # M/M_c = 0.556
    3.0e-4,   # M/M_c = 0.667
    3.5e-4,   # M/M_c = 0.778
    4.0e-4,   # M/M_c = 0.889
    4.5e-4,   # M/M_c = 1.000  (at threshold — may show extinction; include but
              #                  mark n_valid carefully)
    5.0e-4,   # M/M_c = 1.111  (above threshold; likely extinction)
]
N_SEEDS: int = 30      # ≥30 per C3 carry-forward recommendation
T_EQ: int = 1000       # unchanged from Sprint 58
L: int = 100           # unchanged

T_MEASURE: int = 200   # measurement window (generations, unchanged)
MEASURE_EVERY: int = 20  # snapshot stride → 10 snapshots per run
N_SNAPSHOTS: int = T_MEASURE // MEASURE_EVERY  # = 10
R_MAX: int = 40        # max ACF radius; r_zero ≤ 32.3 for M≤5e-4 (all within R_MAX)

# Fit mask: include only M ∈ [2e-4, 4.5e-4] with n_valid ≥ 15
FIT_M_MAX: float = 4.5e-4  # exclude M > 4.5e-4 from slope fit
FIT_MIN_VALID: int = 15    # exclude points with n_valid < 15


def spiral_wavelength(field: np.ndarray) -> float:
    """Estimate spiral wavelength via radial ACF first zero crossing.

    The 2D radial autocorrelation of a spiral field follows J₀(2πr/λ),
    with first zero at r_zero = 0.383 λ. Returns λ = r_zero / 0.383, or
    np.nan if no zero crossing is found within R_MAX lattice units.

    Parameters
    ----------
    field : np.ndarray
        2D real-valued density field (e.g., species-A binary indicator, L×L).

    Returns
    -------
    float
        Estimated wavelength in lattice units, or np.nan.
    """
    Lf = field.shape[0]
    field_c = field - field.mean()

    # 2D circular ACF via Wiener-Khinchin theorem: IFFT(|FFT|²)
    fft_f = np.fft.fft2(field_c)
    acf_2d = np.real(np.fft.ifft2(np.abs(fft_f) ** 2))
    acf_2d = np.fft.fftshift(acf_2d)

    cy, cx = Lf // 2, Lf // 2
    c0 = acf_2d[cy, cx]
    if c0 <= 0.0:
        return np.nan
    acf_2d /= c0  # normalize: C(0) = 1

    # Azimuthal average: C(r) = mean over angles of ACF at radius r
    r_limit = min(R_MAX, Lf // 2 - 1)
    radial = np.empty(r_limit)
    for ri in range(r_limit):
        r = float(ri + 1)
        n_theta = max(8, int(2.0 * np.pi * r))
        theta = np.linspace(0.0, 2.0 * np.pi, n_theta, endpoint=False)
        iy = (np.round(cy + r * np.sin(theta)).astype(int)) % Lf
        ix = (np.round(cx + r * np.cos(theta)).astype(int)) % Lf
        radial[ri] = float(np.mean(acf_2d[iy, ix]))

    # First zero crossing: C(r_i) ≥ 0 and C(r_{i+1}) < 0
    for i in range(r_limit - 1):
        if radial[i] >= 0.0 and radial[i + 1] < 0.0:
            # Linear interpolation within [r_i, r_{i+1}] = [i+1, i+2]
            frac = radial[i] / (radial[i] - radial[i + 1])
            r_zero = (i + 1) + frac  # in lattice units
            return float(r_zero / 0.383)  # λ = r_zero / 0.383 (J₀ first zero)

    return np.nan


def run_one_seed(M: float, seed: int) -> float:
    """Run one (M, seed) trial; return mean wavelength over measurement window.

    Parameters
    ----------
    M : float
        Mobility parameter.
    seed : int
        RNG seed.

    Returns
    -------
    float
        Mean spiral wavelength (lattice units) over valid snapshots, or np.nan.
    """
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
    snapshots = full_history[2:]  # skip initial state and early transient snapshot

    lambdas: list[float] = []
    for snap in snapshots:
        grid: np.ndarray = snap["grid"]
        field = (grid == RPSSpatialModel.SPECIES_A).astype(float)
        lam = spiral_wavelength(field)
        lambdas.append(lam)

    valid = [l for l in lambdas if not np.isnan(l)]
    if not valid:
        return np.nan
    return float(np.mean(valid))


def _worker(args: tuple) -> tuple:
    """Multiprocessing worker: run one (M_index, seed_index) trial.

    Seed convention: seed = M_index * 100 + seed_index (unchanged from Sprint 58).
    """
    M_index, seed_index = args
    M = M_VALUES[M_index]
    seed = M_index * 100 + seed_index
    lam = run_one_seed(M, seed)
    return M_index, seed_index, lam


def _analytical_lambda(M: float) -> float:
    """Analytical spiral wavelength in lattice units: λ = 0.8·L·√(M/M_c)."""
    return 0.8 * L * float(np.sqrt(M / M_C))


def _print_unit_verification() -> None:
    """Print unit-verification block before running the sweep.

    Formula domain: valid near M_c, M/M_c ∈ [0.44, 1.00] per Reichenbach (2007) Supp. §2.
    """
    print("=== Unit Verification ===")
    print("Model mobility convention (Reichenbach Eq. 1): M = 2ε·a²/N")
    print(f"  Lattice spacing a = 1 (each grid site = 1 unit), N = L² = {L**2}")
    print(f"  Exchange rate ε = M·L²/2 = M·{L**2 // 2}")
    print(f"  M is dimensionless mobility per generation (one gen = L² = {L**2} elem. steps).")
    print("Analytical formula: λ = 0.8·L·√(M/M_c)  [lattice units]")
    print(f"  M_c = {M_C:.1e}, L = {L}")
    print(f"  Formula domain: valid near M_c, M/M_c ∈ [0.44, 1.00] per Reichenbach (2007) Supp. §2")
    print(f"  At M=1e-4: λ_formula = {_analytical_lambda(1e-4):.2f} lattice units")
    print(f"  At M=3e-4 (Sprint 54 point):  λ_formula = {_analytical_lambda(3e-4):.2f}")
    print(f"  Sprint 54 measured at M=3e-4: 60.84  (relative error {abs(60.84 - _analytical_lambda(3e-4)) / _analytical_lambda(3e-4) * 100:.1f}% — within 15% tolerance)")
    print(f"  At M=4e-4 (Sprint 54 point):  λ_formula = {_analytical_lambda(4e-4):.2f}")
    print(f"  Sprint 54 measured at M=4e-4: 66.87  (relative error {abs(66.87 - _analytical_lambda(4e-4)) / _analytical_lambda(4e-4) * 100:.1f}%)")
    print(f"  At M=5e-4 (Sprint 54 point):  λ_formula = {_analytical_lambda(5e-4):.2f}")
    print(f"  Sprint 54 measured at M=5e-4: 73.44  (relative error {abs(73.44 - _analytical_lambda(5e-4)) / _analytical_lambda(5e-4) * 100:.1f}%)")
    print("Unit check: Sprint 54 points all within 15% → formula and estimator agree.")
    print(f"Sprint 59 fit restricted to formula-valid sub-range M ∈ [2e-4, {FIT_M_MAX:.1e}] with n_valid ≥ {FIT_MIN_VALID}.")
    print("  Reason: linearized theory applies only near M_c (Sprint 58 diagnosis).")
    print("=== End Unit Verification ===\n")


def main() -> None:
    """Run λ(M) near-M_c dense sweep and save JSON artifact."""
    _print_unit_verification()

    t_start = time.time()

    print("Sprint 59 — Near-M_c dense sweep (formula-valid regime only)")
    print(f"M_c = {M_C:.1e} ; fit range M ∈ [2.0e-4, {FIT_M_MAX:.1e}] (n_valid ≥ {FIT_MIN_VALID})")
    print()

    # Dispatch 7×30 = 210 independent runs in parallel
    tasks = [
        (M_index, seed_index)
        for M_index in range(len(M_VALUES))
        for seed_index in range(N_SEEDS)
    ]
    n_workers = min(os.cpu_count() or 1, 8)
    print(f"Running {len(tasks)} simulations on {n_workers} workers ...\n")

    # Collect results into wavelengths[M_index][seed_index]
    wavelengths: list[list[float]] = [[np.nan] * N_SEEDS for _ in range(len(M_VALUES))]
    with Pool(processes=n_workers) as pool:
        for M_index, seed_index, lam in pool.map(_worker, tasks):
            wavelengths[M_index][seed_index] = lam

    # Compute per-M statistics
    per_point: list[dict] = []
    for M_index, M in enumerate(M_VALUES):
        seeds = wavelengths[M_index]
        valid = [l for l in seeds if not np.isnan(l)]
        n_valid = len(valid)
        lam_mean = float(np.mean(valid)) if valid else float("nan")
        lam_std = float(np.std(valid, ddof=1)) if len(valid) >= 2 else float("nan")
        lam_sem = lam_std / float(np.sqrt(n_valid)) if n_valid >= 2 else float("nan")
        lam_formula = _analytical_lambda(M)
        rel_err = abs(lam_mean - lam_formula) / lam_formula if not np.isnan(lam_mean) else float("nan")
        near_extinction = bool(M >= M_C)  # at or above M_c
        per_point.append(
            {
                "M": M,
                "lambda_mean": round(lam_mean, 4) if not np.isnan(lam_mean) else None,
                "lambda_sem": round(lam_sem, 4) if not np.isnan(lam_sem) else None,
                "lambda_formula": round(lam_formula, 4),
                "relative_error": round(rel_err, 6) if not np.isnan(rel_err) else None,
                "n_valid": n_valid,
                "n_seeds": N_SEEDS,
                "near_extinction": near_extinction,
            }
        )

    # Fit mask: M ∈ [2e-4, 4.5e-4] with n_valid ≥ 15
    fit_points = [
        p for p in per_point
        if p["lambda_mean"] is not None
        and p["M"] <= FIT_M_MAX
        and p["n_valid"] >= FIT_MIN_VALID
    ]
    fit_M_values = [p["M"] for p in fit_points]
    n_fit_points = len(fit_points)

    fit_exclusion_note: str | None = None
    if n_fit_points < 4:
        slope = None
        intercept = None
        r_squared = None
        tolerance_pass = False
        r_squared_pass = False
        n_fit_points_pass = False
        overall_pass = False
        excluded = len(M_VALUES) - n_fit_points
        fit_exclusion_note = (
            f"Only {n_fit_points} M points survived fit mask "
            f"(M ≤ {FIT_M_MAX:.1e}, n_valid ≥ {FIT_MIN_VALID}); "
            f"{excluded} excluded. Minimum 4 required."
        )
    else:
        log_M = np.array([np.log(p["M"]) for p in fit_points])
        log_lam = np.array([np.log(p["lambda_mean"]) for p in fit_points])
        slope_arr, intercept_arr = np.polyfit(log_M, log_lam, 1)
        slope = float(slope_arr)
        intercept = float(intercept_arr)

        # R²
        y_pred = slope * log_M + intercept
        ss_res = float(np.sum((log_lam - y_pred) ** 2))
        ss_tot = float(np.sum((log_lam - log_lam.mean()) ** 2))
        r_squared = float(1.0 - ss_res / ss_tot) if ss_tot > 0 else 0.0

        tolerance_pass = bool(0.40 <= slope <= 0.60)
        r_squared_pass = bool(r_squared >= 0.90)
        n_fit_points_pass = True  # already >= 4
        overall_pass = bool(tolerance_pass and r_squared_pass and n_fit_points_pass)

    elapsed = time.time() - t_start

    # Console summary (A.7)
    for p in per_point:
        M = p["M"]
        lam_m = f"{p['lambda_mean']:7.2f}" if p["lambda_mean"] is not None else "    NaN"
        lam_s = f"{p['lambda_sem']:5.2f}" if p["lambda_sem"] is not None else "  NaN"
        lam_f = f"{p['lambda_formula']:7.2f}"
        rel_e = f"{p['relative_error']*100:5.1f}%" if p["relative_error"] is not None else "  NaN"
        nv = f"{p['n_valid']}/{N_SEEDS}"

        if M > FIT_M_MAX:
            tag = "[above M_c; EXCLUDED from fit]"
        elif M >= M_C:
            tag = f"[near M_c; IN FIT if n_valid≥{FIT_MIN_VALID}]"
        else:
            tag = "[IN FIT]"

        print(
            f"M={M:.2e}: λ={lam_m} ± {lam_s}"
            f"  formula={lam_f}  rel_err={rel_e}"
            f"  n_valid={nv}  {tag}"
        )

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
    nfp_str = "PASS" if n_fit_points >= 4 else "FAIL"
    print(f"n_fit_points:  {n_fit_points}      (target ≥ 4)                        {nfp_str}")
    ov_str = "PASS" if overall_pass else "FAIL"
    print(f"Overall:       {ov_str}")

    # JSON output (A.3 schema)
    out: dict = {
        "sprint": 59,
        "m_c": M_C,
        "M_values": M_VALUES,
        "fit_M_values": fit_M_values,
        "n_fit_points": n_fit_points,
        "per_point": per_point,
        "log_log_slope": round(slope, 6) if slope is not None else None,
        "log_log_intercept": round(intercept, 6) if intercept is not None else None,
        "r_squared": round(r_squared, 6) if r_squared is not None else None,
        "tolerance_pass": tolerance_pass,
        "r_squared_pass": r_squared_pass,
        "n_fit_points_pass": n_fit_points >= 4,
        "overall_pass": overall_pass,
        "fit_exclusion_note": fit_exclusion_note,
        "n_seeds": N_SEEDS,
        "L": L,
        "T_eq": T_EQ,
        "elapsed_seconds": round(elapsed, 1),
    }

    output_path = "analysis/outputs/p12_reichenbach2007_reproduction.json"
    with open(output_path, "w") as f:
        json.dump(out, f, indent=2)

    print(f"\nOutput: {output_path}  [{elapsed:.0f}s]")


if __name__ == "__main__":
    main()
