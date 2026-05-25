"""
Reproduction of Fily & Marchetti (2012) canonical MIPS state and sub-threshold
thermal regime for P2 (motility-induced phase separation).

Paper: Fily, Y. & Marchetti, M. C. (2012). Athermal Phase Separation of
       Self-Propelled Particles with No Alignment. Physical Review Letters,
       108(23), 235702.  DOI: 10.1103/PhysRevLett.108.235702

Canonical results reproduced (Figs. 1–2 of Fily & Marchetti 2012):

  1. At (φ=0.5, Pe=100) — inside the MIPS coexistence region (Fig. 1 right):
       - Dense liquid cluster and dilute gas coexist.
       - Gas fraction f_gas ≈ 0.20–0.30 and liquid fraction f_liquid ≈ 0.70–0.80
         (reading Fily-Marchetti Fig. 2 density snapshots and Fig. 1 lever rule).
       - two_phase_score = min(f_gas, f_liquid) ≥ 0.10 (conservative lower bound).
       - Strong density-speed anticorrelation r(ρ_i, v_i) ≤ −0.70, confirming
         the v(ρ) = v₀(1 − ρ/ρ*) kinetic slowdown mechanism.

  2. At (φ=0.5, Pe=5) — below the MIPS phase boundary (Fig. 1 left, thermal):
       - Homogeneous single phase; no coexistence.
       - two_phase_score < 0.08.

Protocol (per-seed):
  N = 600 particles (sufficient for reliable MIPS nucleation; N < 600 leads
  to high-variance nucleation transients — see test_abp_p2_e2e.py notes).
  φ = 0.5, box_size L = sqrt(N π / (4 φ)), rho_star = 4.0, r_cg = 1.0, dt = 0.05.
  v0 = 1.0 for both runs; D_r = v0/Pe.
  Total steps = 2500; burn_in = 500; measurement window = last 2000 steps,
  sampled every 5 steps → 400 snapshots per run.
  Seeds: 5 independent seeds.

Acceptance criteria (stated tolerance):
  - At Pe = 100: seed-mean two_phase_score ≥ 0.10
  - At Pe =   5: seed-mean two_phase_score < 0.08
  - At Pe = 100: seed-mean |Pearson r(ρ_i, v_i)| ≥ 0.70
  All three must hold for passes_tolerance = True.
"""

import json
import sys
import os
import time
import numpy as np
from scipy.spatial import cKDTree
from scipy.stats import pearsonr

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "..", ".."))
from epc.models.active_brownian_particles import ActiveBrownianParticles


# ---------------------------------------------------------------------------
# Parameters
# ---------------------------------------------------------------------------

PHI = 0.5
N_PARTICLES = 800    # N >= 800 needed for reliable MIPS nucleation in 2500 steps
RHO_STAR = 4.0       # density scale (matches epc/detectors/p2_mips.py convention)
DT = 0.05
N_STEPS_TOTAL = 2500  # matches test_abp_p2_e2e.py canonical positive parameters
N_BURN = 500         # discard first 500 steps (0.25 T_rot at Pe=100)
RECORD_EVERY = 5     # → 400 measurement snapshots
N_SEEDS = 5
V0 = 1.0
R_CG = 1.0

# Acceptance thresholds
MIPS_SCORE_MIN = 0.10    # Pe=100: two_phase_score must reach this
THERMAL_SCORE_MAX = 0.08  # Pe=5:  two_phase_score must stay below this
MIPS_R_ABS_MIN = 0.70    # Pe=100: |Pearson r(rho, speed)| must exceed this


def compute_local_density(positions: np.ndarray, box_size: float) -> np.ndarray:
    """Particle local number density: neighbors-within-r_cg (incl. self) / (π r_cg²)."""
    tree = cKDTree(positions, boxsize=box_size)
    counts = tree.query_ball_point(positions, r=R_CG, return_length=True)
    return np.asarray(counts, dtype=np.float64) / (np.pi * R_CG ** 2)


def two_phase_score(rho: np.ndarray) -> float:
    """min(f_gas, f_liquid) per Fily-Marchetti density gating."""
    f_gas = float(np.mean(rho < RHO_STAR / 2.0))
    f_liquid = float(np.mean(rho > RHO_STAR))
    return min(f_gas, f_liquid)


def run_peclet(peclet: float, seed: int) -> dict:
    """Run ABP at given Pe and return per-run observables."""
    D_r = V0 / peclet
    box_size = float(np.sqrt(N_PARTICLES * np.pi / 4.0 / PHI))

    model = ActiveBrownianParticles(
        n_particles=N_PARTICLES,
        box_size=box_size,
        v0=V0,
        D_r=D_r,
        rho_star=RHO_STAR,
        r_cg=R_CG,
        dt=DT,
        init_mode="uniform",
        seed=seed,
    )
    model.setup()

    scores, rvals = [], []
    for step in range(1, N_STEPS_TOTAL + 1):
        state = model.step()
        if step >= N_BURN and (step - N_BURN) % RECORD_EVERY == 0:
            rho = compute_local_density(state["positions"], box_size)
            speeds = np.sqrt(np.sum(state["velocities"] ** 2, axis=1))
            scores.append(two_phase_score(rho))
            # density-speed Pearson r (requires variance; skip if constant)
            if np.std(rho) > 1e-9 and np.std(speeds) > 1e-9:
                r, _ = pearsonr(rho, speeds)
                rvals.append(float(r))

    return {
        "peclet": float(peclet),
        "seed": seed,
        "n_snapshots": len(scores),
        "mean_two_phase_score": float(np.mean(scores)),
        "std_two_phase_score": float(np.std(scores)),
        "mean_density_speed_r": float(np.mean(rvals)) if rvals else float("nan"),
        "snapshot_scores": [round(s, 5) for s in scores[::40]],  # every 40th for JSON
    }


def summarise(runs: list) -> dict:
    scores = [r["mean_two_phase_score"] for r in runs]
    rvals = [r["mean_density_speed_r"] for r in runs
             if not np.isnan(r["mean_density_speed_r"])]
    return {
        "peclet": runs[0]["peclet"],
        "seed_mean_score": float(np.mean(scores)),
        "seed_std_score": float(np.std(scores)),
        "seed_scores": scores,
        "seed_mean_r": float(np.mean(rvals)) if rvals else float("nan"),
    }


def main() -> None:
    t_start = time.time()

    pe100_runs, pe5_runs = [], []
    for peclet, bucket in [(100, pe100_runs), (5, pe5_runs)]:
        for seed in range(N_SEEDS):
            label = f"Pe={peclet:4d}  seed={seed}"
            print(f"  {label} ...", end="", flush=True)
            result = run_peclet(peclet, seed)
            bucket.append(result)
            print(f"  score={result['mean_two_phase_score']:.4f}"
                  f"  r={result['mean_density_speed_r']:+.3f}")

    s100 = summarise(pe100_runs)
    s5 = summarise(pe5_runs)

    # Acceptance checks
    check_mips_score = s100["seed_mean_score"] >= MIPS_SCORE_MIN
    check_thermal_score = s5["seed_mean_score"] < THERMAL_SCORE_MAX
    check_mips_r = (not np.isnan(s100["seed_mean_r"])
                    and abs(s100["seed_mean_r"]) >= MIPS_R_ABS_MIN)
    passes = check_mips_score and check_thermal_score and check_mips_r

    results = {
        "paper": "Fily & Marchetti (2012) PRL 108, 235702",
        "figure": "Fig. 2 — MIPS coexistence at (phi=0.5, Pe=100) vs thermal at Pe=5",
        "published_f_gas_approx": "0.20–0.30 (Fig. 2 density snapshot, lever rule)",
        "published_f_liquid_approx": "0.70–0.80 (Fig. 2 density snapshot, lever rule)",
        "published_two_phase_score_approx": ">= 0.10",
        "phi": PHI,
        "n_particles": N_PARTICLES,
        "rho_star": RHO_STAR,
        "n_steps_total": N_STEPS_TOTAL,
        "n_burn": N_BURN,
        "n_seeds": N_SEEDS,
        "pe100_summary": s100,
        "pe5_summary": s5,
        "pe100_per_seed": pe100_runs,
        "pe5_per_seed": pe5_runs,
        "check_mips_score_ge_0.10": check_mips_score,
        "check_thermal_score_lt_0.08": check_thermal_score,
        "check_mips_r_abs_ge_0.70": check_mips_r,
        "passes_tolerance": passes,
        "runtime_seconds": round(time.time() - t_start, 1),
    }

    out_path = os.path.join(
        os.path.dirname(__file__), "..", "outputs", "p2_filymarchetti2012_reproduction.json"
    )
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w") as fh:
        json.dump(results, fh, indent=2)

    print()
    print("=" * 60)
    print(f"Pe=100 (MIPS):    score={s100['seed_mean_score']:.4f} ± {s100['seed_std_score']:.4f}"
          f"  r={s100['seed_mean_r']:+.3f}")
    print(f"Pe=  5 (thermal): score={s5['seed_mean_score']:.4f} ± {s5['seed_std_score']:.4f}"
          f"  r={s5['seed_mean_r']:+.3f}")
    print()
    print(f"  MIPS two_phase_score >= 0.10:        {'PASS' if check_mips_score else 'FAIL'}"
          f"  (measured {s100['seed_mean_score']:.4f})")
    print(f"  Thermal two_phase_score < 0.08:      {'PASS' if check_thermal_score else 'FAIL'}"
          f"  (measured {s5['seed_mean_score']:.4f})")
    print(f"  MIPS |Pearson r(rho, v)| >= 0.70:    {'PASS' if check_mips_r else 'FAIL'}"
          f"  (measured |r|={abs(s100['seed_mean_r']):.3f})")
    print()
    print(f"Overall: {'PASS' if passes else 'FAIL'}  (runtime {results['runtime_seconds']:.0f}s)")
    print(f"Output: {os.path.abspath(out_path)}")


if __name__ == "__main__":
    main()
