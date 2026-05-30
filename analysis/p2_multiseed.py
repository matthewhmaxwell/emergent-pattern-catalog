"""
P2 (MIPS / active Brownian particles) multi-seed extension — Sprint 60 dim2 closure.

Runs N_SEEDS independent seeds at the canonical MIPS state (φ=0.5, Pe=100)
and reports per-seed two_phase_score plus aggregate (mean, std, CV).

Primary observable: two_phase_score = min(f_gas, f_liquid) per Fily-Marchetti
density gating. Secondary: Pearson r(ρ_i, v_i) density-speed anticorrelation.

Parameters match the existing reproduction script (p2_filymarchetti2012.py):
  N=800, φ=0.5, rho_star=4.0, dt=0.05, v0=1.0, r_cg=1.0
  Pe=100, T_total=2500, burn_in=500, record_every=5 → 400 snapshots

Output: analysis/outputs/p2_multiseed.json
"""

from __future__ import annotations

import json
import os
import sys
import time
from typing import Any, Dict, List

import numpy as np
from scipy.spatial import cKDTree
from scipy.stats import pearsonr

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from epc.models.active_brownian_particles import ActiveBrownianParticles

# Parameters (match p2_filymarchetti2012.py canonical positive)
PHI = 0.5
N_PARTICLES = 800
RHO_STAR = 4.0
DT = 0.05
V0 = 1.0
R_CG = 1.0
PECLET = 100
N_STEPS_TOTAL = 2500
N_BURN = 500
RECORD_EVERY = 5

N_SEEDS = 20
SEED_BASE = 100


def compute_local_density(positions: np.ndarray, box_size: float) -> np.ndarray:
    """Particle local number density: neighbors-within-r_cg / (π r_cg²)."""
    tree = cKDTree(positions, boxsize=box_size)
    counts = tree.query_ball_point(positions, r=R_CG, return_length=True)
    return np.asarray(counts, dtype=np.float64) / (np.pi * R_CG ** 2)


def two_phase_score(rho: np.ndarray) -> float:
    """min(f_gas, f_liquid) per Fily-Marchetti density gating."""
    f_gas = float(np.mean(rho < RHO_STAR / 2.0))
    f_liquid = float(np.mean(rho > RHO_STAR))
    return min(f_gas, f_liquid)


def run_seed(seed: int) -> Dict[str, Any]:
    """Run one ABP seed at Pe=100 and return per-seed observables."""
    D_r = V0 / PECLET
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

    scores: List[float] = []
    rvals: List[float] = []
    for step in range(1, N_STEPS_TOTAL + 1):
        state = model.step()
        if step >= N_BURN and (step - N_BURN) % RECORD_EVERY == 0:
            rho = compute_local_density(state["positions"], box_size)
            speeds = np.sqrt(np.sum(state["velocities"] ** 2, axis=1))
            scores.append(two_phase_score(rho))
            if np.std(rho) > 1e-9 and np.std(speeds) > 1e-9:
                r, _ = pearsonr(rho, speeds)
                rvals.append(float(r))

    return {
        "seed": seed,
        "two_phase_score": float(np.mean(scores)),
        "two_phase_score_std": float(np.std(scores)),
        "density_speed_r": float(np.mean(rvals)) if rvals else float("nan"),
        "n_snapshots": len(scores),
    }


def main() -> None:
    t0 = time.time()
    print(f"P2 MIPS multi-seed: N={N_PARTICLES}, φ={PHI}, Pe={PECLET}, "
          f"N_seeds={N_SEEDS}")

    per_seed: List[Dict[str, Any]] = []
    for i in range(N_SEEDS):
        seed = SEED_BASE + i
        t_seed = time.time()
        result = run_seed(seed)
        elapsed = time.time() - t_seed
        per_seed.append(result)
        print(f"  seed {seed:3d}: score={result['two_phase_score']:.4f}  "
              f"r={result['density_speed_r']:+.3f}  ({elapsed:.1f}s)")

    scores = np.array([r["two_phase_score"] for r in per_seed])
    score_mean = float(np.mean(scores))
    score_std = float(np.std(scores, ddof=1))
    score_cv = float(score_std / score_mean) if score_mean > 0 else float("nan")

    rvals = np.array([r["density_speed_r"] for r in per_seed
                      if not np.isnan(r["density_speed_r"])])
    r_mean = float(np.mean(rvals)) if len(rvals) > 0 else float("nan")
    r_std = float(np.std(rvals, ddof=1)) if len(rvals) > 1 else float("nan")
    r_cv = float(r_std / abs(r_mean)) if (len(rvals) > 1 and abs(r_mean) > 1e-9) else float("nan")

    print(f"\ntwo_phase_score aggregate ({N_SEEDS} seeds):")
    print(f"  mean = {score_mean:.4f}")
    print(f"  std  = {score_std:.4f}")
    print(f"  CV   = {score_cv:.4f}")
    print(f"\ndensity_speed_r aggregate ({len(rvals)} valid seeds):")
    print(f"  mean = {r_mean:.4f}")
    print(f"  std  = {r_std:.4f}")
    print(f"  CV   = {r_cv:.4f}")

    output: Dict[str, Any] = {
        "sprint": 60,
        "description": "P2 MIPS multi-seed two_phase_score extension",
        "params": {
            "n_particles": N_PARTICLES,
            "phi": PHI,
            "peclet": PECLET,
            "rho_star": RHO_STAR,
            "v0": V0,
            "dt": DT,
            "n_steps_total": N_STEPS_TOTAL,
            "n_burn": N_BURN,
            "record_every": RECORD_EVERY,
            "n_seeds": N_SEEDS,
            "seed_base": SEED_BASE,
        },
        "primary_observable": "two_phase_score = min(f_gas, f_liquid)",
        "per_seed": per_seed,
        "aggregate": {
            "score_mean": score_mean,
            "score_std": score_std,
            "score_cv": score_cv,
            "score_min": float(scores.min()),
            "score_max": float(scores.max()),
            "r_mean": r_mean,
            "r_std": r_std,
            "r_cv": r_cv,
        },
        "elapsed_seconds": time.time() - t0,
    }

    out_path = "analysis/outputs/p2_multiseed.json"
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(output, f, indent=2)
    print(f"\nSaved → {out_path}  (elapsed {output['elapsed_seconds']:.1f}s)")


if __name__ == "__main__":
    main()
