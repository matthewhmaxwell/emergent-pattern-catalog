"""
P12 (spatial RPS) multi-seed extension — Sprint 56 dim2 closure.

Runs N_SEEDS independent seeds of the Reichenbach-Mobilia-Frey (2007)
spatial RPS model at one canonical mobility M=1e-4 (coexistence regime,
M < M_c ≈ 4.5e-4), measures spiral wavelength λ per seed via the radial
ACF first-zero estimator, and reports aggregate (mean, std, CV).

Parameters: L=100, σ=μ=1, M=1e-4 (canonical coexistence regime),
T_EQ=500 generations (equilibration), T_MEASURE=200 generations
(measurement window), snapshot stride=20 (10 snapshots per seed).

Wavelength estimator — radial ACF first zero (same as p12_reichenbach2007.py):
  λ = r_zero / 0.383  where r_zero is the first zero crossing of the
  azimuthally-averaged 2D ACF of the species-A density field.

At M=1e-4: λ_formula = 0.8 × 100 × √(1e-4/4.5e-4) ≈ 37.7 lattice units;
r_zero ≈ 14.4 (well within R_MAX=40).

Dim2 closure criterion: CV of per-seed mean λ values < 15% (typical
spiral characterisation tolerance).

Output: analysis/outputs/p12_multiseed.json
"""

from __future__ import annotations

import json
import time
from typing import Any, Dict, List, Optional

import numpy as np

from epc.models.rps_spatial import RPSSpatialModel

# --- Parameters ---
L: int = 100                   # lattice size
M: float = 1e-4                # canonical mobility (coexistence regime)
M_C: float = 4.5e-4           # critical mobility (Reichenbach 2007)
T_EQ: int = 500               # equilibration (generations)
T_MEASURE: int = 200          # measurement window (generations)
MEASURE_EVERY: int = 20       # snapshot stride → 10 snapshots per seed
R_MAX: int = 40               # max ACF radius

N_SEEDS: int = 20
SEED_BASE: int = 100          # seeds 100–119; avoids reproduction seeds 0–9


def spiral_wavelength(field: np.ndarray) -> float:
    """Estimate spiral wavelength via radial ACF first zero crossing.

    Same estimator as analysis/reproductions/p12_reichenbach2007.py.
    λ = r_zero / 0.383 (J₀ first zero at 2.405 rad).

    Returns np.nan if no zero crossing found within R_MAX.
    """
    Lf = field.shape[0]
    field_c = field - field.mean()
    fft_f = np.fft.fft2(field_c)
    acf_2d = np.real(np.fft.ifft2(np.abs(fft_f) ** 2))
    acf_2d = np.fft.fftshift(acf_2d)
    cy, cx = Lf // 2, Lf // 2
    c0 = acf_2d[cy, cx]
    if c0 <= 0.0:
        return np.nan
    acf_2d /= c0

    r_limit = min(R_MAX, Lf // 2 - 1)
    radial = np.empty(r_limit)
    for ri in range(r_limit):
        r = float(ri + 1)
        n_theta = max(8, int(2.0 * np.pi * r))
        theta = np.linspace(0.0, 2.0 * np.pi, n_theta, endpoint=False)
        iy = (np.round(cy + r * np.sin(theta)).astype(int)) % Lf
        ix = (np.round(cx + r * np.cos(theta)).astype(int)) % Lf
        radial[ri] = float(np.mean(acf_2d[iy, ix]))

    for i in range(r_limit - 1):
        if radial[i] >= 0.0 and radial[i + 1] < 0.0:
            frac = radial[i] / (radial[i] - radial[i + 1])
            r_zero = (i + 1) + frac
            return float(r_zero / 0.383)

    return np.nan


def run_seed(seed: int) -> Dict[str, Any]:
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
    # Equilibration (T_EQ generations); no snapshots needed
    model.run(T_EQ, record_every=T_EQ)
    # Measurement window
    full_history = model.run(T_MEASURE, record_every=MEASURE_EVERY)
    snapshots = full_history[2:]  # skip initial state and equilibration endpoint

    lambdas: List[Optional[float]] = []
    for snap in snapshots:
        grid: np.ndarray = snap["grid"]
        field = (grid == RPSSpatialModel.SPECIES_A).astype(float)
        lam = spiral_wavelength(field)
        lambdas.append(None if np.isnan(lam) else float(lam))

    valid = [l for l in lambdas if l is not None]
    mean_lam = float(np.mean(valid)) if valid else float("nan")
    return {
        "seed": seed,
        "mean_wavelength": mean_lam,
        "n_valid_snapshots": len(valid),
        "n_total_snapshots": len(lambdas),
        "per_snapshot_wavelengths": lambdas,
    }


def main() -> None:
    lam_expected = 0.8 * L * float(np.sqrt(M / M_C))
    t0 = time.time()
    print(
        f"P12 RPS multi-seed: L={L}, M={M:.1e} (M_c={M_C:.1e}), "
        f"T_eq={T_EQ} gen, T_measure={T_MEASURE} gen, "
        f"N_seeds={N_SEEDS}, λ_formula={lam_expected:.1f}"
    )

    per_seed: List[Dict[str, Any]] = []
    for i in range(N_SEEDS):
        seed = SEED_BASE + i
        t_seed = time.time()
        result = run_seed(seed)
        elapsed = time.time() - t_seed
        per_seed.append(result)
        print(
            f"  seed {seed:3d}: λ={result['mean_wavelength']:.2f}  "
            f"n_valid={result['n_valid_snapshots']}/{result['n_total_snapshots']}  "
            f"({elapsed:.1f}s)"
        )

    valid_means = [r["mean_wavelength"] for r in per_seed
                   if not np.isnan(r["mean_wavelength"])]
    lam_mean = float(np.mean(valid_means)) if valid_means else float("nan")
    lam_std = float(np.std(valid_means, ddof=1)) if len(valid_means) > 1 else float("nan")
    lam_cv = float(lam_std / lam_mean) if lam_mean > 0 else float("nan")
    n_valid_seeds = len(valid_means)

    print(f"\nλ aggregate ({n_valid_seeds}/{N_SEEDS} valid seeds):")
    print(f"  mean = {lam_mean:.2f}")
    print(f"  std  = {lam_std:.2f}")
    print(f"  CV   = {lam_cv:.4f}")
    print(f"  λ_formula = {lam_expected:.2f}")
    print(f"  rel_error = {abs(lam_mean - lam_expected) / lam_expected:.3f}")

    output: Dict[str, Any] = {
        "sprint": 56,
        "description": "P12 spatial RPS multi-seed spiral wavelength at M=1e-4",
        "params": {
            "L": L,
            "M": M,
            "M_c": M_C,
            "T_eq_generations": T_EQ,
            "T_measure_generations": T_MEASURE,
            "measure_stride": MEASURE_EVERY,
            "r_max_acf": R_MAX,
            "n_seeds": N_SEEDS,
            "seed_base": SEED_BASE,
        },
        "primary_observable": "spiral wavelength λ (radial ACF first-zero estimator)",
        "expected_wavelength_formula": round(lam_expected, 2),
        "per_seed": per_seed,
        "aggregate": {
            "lambda_mean": lam_mean,
            "lambda_std": lam_std,
            "lambda_cv": lam_cv,
            "lambda_min": float(np.nanmin([r["mean_wavelength"] for r in per_seed])),
            "lambda_max": float(np.nanmax([r["mean_wavelength"] for r in per_seed])),
            "n_valid_seeds": n_valid_seeds,
        },
        "elapsed_seconds": time.time() - t0,
    }

    out_path = "analysis/outputs/p12_multiseed.json"
    with open(out_path, "w") as f:
        json.dump(output, f, indent=2)
    print(f"\nSaved → {out_path}  (elapsed {output['elapsed_seconds']:.1f}s)")


if __name__ == "__main__":
    main()
