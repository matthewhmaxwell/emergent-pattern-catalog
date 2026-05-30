"""
P22 (SIR epidemic wavefront speed) multi-seed extension — Sprint 60 dim2 closure.

Runs N_SEEDS independent seeds of the Datta-Acharyya (2021) paper-exact CA
and reports per-seed wavefront speed plus aggregate (mean, std, CV).

Primary observable: wavefront speed (cells/step) via linear fit of R(t) over
steps 5–100. Published value: 0.4405 ± 0.0008 cells/step.

Parameters match the existing reproduction script (p22_dattaacharyya2005.py):
  L=200, p0=0.25, p1=0.97, p2=0.10, t_tau=4, Von Neumann neighbourhood,
  single-seed IC at centre, fit window [5, 100].

Output: analysis/outputs/p22_multiseed.json
"""

from __future__ import annotations

import json
import os
import sys
import time
from typing import Any, Dict, List

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))

# Paper parameters
PAPER_P0: float = 0.25
PAPER_P1: float = 0.97
PAPER_P2: float = 0.10
PAPER_T_TAU: int = 4
MODEL_L: int = 200
N_STEPS: int = 110
FIT_T_MIN: int = 5
FIT_T_MAX: int = 100

N_SEEDS = 20
SEED_BASE = 100


def _vn_infected_count(infected: np.ndarray) -> np.ndarray:
    """Count Von Neumann infected neighbours for every cell."""
    inf_int = infected.astype(np.int32)
    return (
        np.roll(inf_int, -1, axis=0)
        + np.roll(inf_int, +1, axis=0)
        + np.roll(inf_int, -1, axis=1)
        + np.roll(inf_int, +1, axis=1)
    )


def run_paper_ca(
    L: int, p0: float, p1: float, p2: float,
    t_tau: int, n_steps: int, seed: int,
) -> list[dict[str, Any]]:
    """Implement the Datta-Acharyya (2021) SIR CA exactly."""
    rng = np.random.default_rng(seed)
    center = L // 2
    grid = np.zeros((L, L), dtype=np.int32)
    grid[center, center] = 1
    R_STATE = t_tau + 1
    history: list[dict[str, Any]] = [{"grid": grid.copy(), "step": 0}]

    for t in range(1, n_steps + 1):
        infected = (grid >= 1) & (grid <= t_tau)
        susceptible = grid == 0
        recovered = grid == R_STATE
        n_inf_nb = _vn_infected_count(infected)
        new_grid = grid.copy()

        # S → I
        s_cand = susceptible & (n_inf_nb > 0)
        p_infect = 1.0 - (1.0 - p0) ** n_inf_nb.astype(float)
        si_mask = s_cand & (rng.random((L, L)) < p_infect)
        new_grid[si_mask] = 1

        # I aging
        aging_mask = infected & (grid < t_tau)
        new_grid[aging_mask] = grid[aging_mask] + 1

        # I → R or dead
        final_mask = infected & (grid == t_tau)
        rand3 = rng.random((L, L))
        new_grid[final_mask & (rand3 < p1)] = R_STATE
        new_grid[final_mask & (rand3 >= p1)] = -1

        # R → I (re-infection)
        r_cand = recovered & (n_inf_nb > 0)
        p_reinfect = 1.0 - (1.0 - p2) ** n_inf_nb.astype(float)
        ri_mask = r_cand & (rng.random((L, L)) < p_reinfect)
        new_grid[ri_mask] = 1

        grid = new_grid
        history.append({"grid": grid.copy(), "step": t})
        if not infected.any():
            break

    return history


def wavefront_radius(grid: np.ndarray, center_r: int, center_c: int) -> float:
    """Euclidean radius of outermost affected cell from centre."""
    affected = (grid >= 1) | (grid == -1)
    if not affected.any():
        return 0.0
    r_idx, c_idx = np.where(affected)
    dists = np.sqrt((r_idx - center_r) ** 2 + (c_idx - center_c) ** 2)
    return float(dists.max())


def run_seed(seed: int) -> Dict[str, Any]:
    """Run one SIR CA seed and return wavefront speed."""
    center = MODEL_L // 2
    history = run_paper_ca(
        L=MODEL_L, p0=PAPER_P0, p1=PAPER_P1, p2=PAPER_P2,
        t_tau=PAPER_T_TAU, n_steps=N_STEPS, seed=seed,
    )

    times: List[float] = []
    radii: List[float] = []
    for snap in history:
        t = snap["step"]
        if FIT_T_MIN <= t <= FIT_T_MAX:
            r = wavefront_radius(snap["grid"], center, center)
            if r > 0:
                times.append(float(t))
                radii.append(r)

    if len(times) < 5:
        return {
            "seed": seed, "speed": None, "r2": None,
            "n_fit_points": len(times), "skipped": True,
        }

    times_arr = np.array(times)
    radii_arr = np.array(radii)
    coeffs = np.polyfit(times_arr, radii_arr, 1)
    slope = float(coeffs[0])
    residuals = radii_arr - np.polyval(coeffs, times_arr)
    ss_res = float(np.sum(residuals ** 2))
    ss_tot = float(np.sum((radii_arr - radii_arr.mean()) ** 2))
    r2 = float(1.0 - ss_res / ss_tot) if ss_tot > 0 else 0.0

    return {
        "seed": seed, "speed": slope, "r2": r2,
        "n_fit_points": len(times), "skipped": False,
    }


def main() -> None:
    t0 = time.time()
    print(f"P22 SIR wavefront speed multi-seed: L={MODEL_L}, "
          f"p0={PAPER_P0}, p1={PAPER_P1}, p2={PAPER_P2}, "
          f"N_seeds={N_SEEDS}")

    per_seed: List[Dict[str, Any]] = []
    for i in range(N_SEEDS):
        seed = SEED_BASE + i
        t_seed = time.time()
        result = run_seed(seed)
        elapsed = time.time() - t_seed
        speed_str = f"{result['speed']:.4f}" if result['speed'] is not None else "SKIP"
        per_seed.append(result)
        print(f"  seed {seed:3d}: speed={speed_str}  ({elapsed:.1f}s)")

    valid_speeds = np.array([r["speed"] for r in per_seed
                             if r["speed"] is not None])
    n_valid = len(valid_speeds)

    if n_valid == 0:
        raise RuntimeError("All seeds died out — parameters are sub-critical.")

    speed_mean = float(np.mean(valid_speeds))
    speed_std = float(np.std(valid_speeds, ddof=1))
    speed_cv = float(speed_std / speed_mean) if speed_mean > 0 else float("nan")

    print(f"\nwavefront_speed aggregate ({n_valid} valid seeds):")
    print(f"  mean = {speed_mean:.4f}")
    print(f"  std  = {speed_std:.4f}")
    print(f"  CV   = {speed_cv:.4f}")

    output: Dict[str, Any] = {
        "sprint": 60,
        "description": "P22 SIR wavefront speed multi-seed extension",
        "params": {
            "model_L": MODEL_L,
            "paper_p0": PAPER_P0,
            "paper_p1": PAPER_P1,
            "paper_p2": PAPER_P2,
            "paper_t_tau": PAPER_T_TAU,
            "n_steps": N_STEPS,
            "fit_range": [FIT_T_MIN, FIT_T_MAX],
            "n_seeds": N_SEEDS,
            "seed_base": SEED_BASE,
        },
        "primary_observable": "wavefront speed (cells/step, linear fit R(t) over steps 5-100)",
        "per_seed": per_seed,
        "aggregate": {
            "speed_mean": speed_mean,
            "speed_std": speed_std,
            "speed_cv": speed_cv,
            "speed_min": float(valid_speeds.min()),
            "speed_max": float(valid_speeds.max()),
            "n_valid_seeds": n_valid,
        },
        "elapsed_seconds": time.time() - t0,
    }

    out_path = "analysis/outputs/p22_multiseed.json"
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(output, f, indent=2)
    print(f"\nSaved → {out_path}  (elapsed {output['elapsed_seconds']:.1f}s)")


if __name__ == "__main__":
    main()
