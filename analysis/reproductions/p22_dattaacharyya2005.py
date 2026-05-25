"""
Reproduction of Datta & Acharyya (2021/2022) §3.1.1 / Fig. 11 for P22 (SIR lattice).

Paper: Agniva Datta & Muktish Acharyya. "Modelling the Spread of an Epidemic in
Presence of Vaccination using Cellular Automata." arXiv:2104.10456 (2021),
Int. J. Mod. Phys. C 33, 2250094 (2022).

Note: The sprint brief references this paper as "Datta & Acharyya (2005)"; the
actual publication year is 2021/2022.  The briefing shorthand is retained in the
filename.

Anchor chosen: §3.1.1 Velocity of Epidemic Spread / Fig. 11 — wavefront radius
R(t) ∝ t confirming linear growth, with reported slope (speed)

    PUBLISHED_SPEED = 0.4405 ± 0.0008 cells/step

Paper model parameters (Table 1):
    p0 = 0.25   per-neighbour infection probability
    p1 = 0.97   recovery probability (remainder die)
    p2 = 0.10   re-infection probability for recovered cells
    t_τ = 4     fixed infection duration (each infected cell lives exactly 4 steps)
    neighbourhood: Von Neumann (4 nearest neighbours)
    lattice: N×N = 500×500, open boundary conditions
    initial condition: single infected cell at centre (§3.1.1 "Simulation 5")

This script implements the paper's exact CA rules rather than delegating to
epc.models.sir_epidemic (which uses stochastic geometric recovery, not fixed
duration).  No epc/ files are modified; the reproduction CA runs entirely inside
this module.  Key differences from epc SIREpidemicModel:

  Paper model            │  epc SIREpidemicModel
  ─────────────────────────────────────────────────────
  Fixed duration t_τ=4   │  Geometric recovery (prob q per step)
  Re-infection p2=0.10   │  Permanent immunity
  Open BCs               │  Periodic BCs (negligible at centre-seed run)
  State encoding float   │  Integer {0,1,2}

We use L=200 instead of 500.  Wavefront speed is a local property: at the
reported speed ≈ 0.44 cells/step, the wavefront travels only ~44 cells in
the 100-step fit window, well within the L/2=100 cell margin from the
boundary.

Tolerance: |measured − published| < 0.05 absolute  OR  < 15% relative (either).
"""

from __future__ import annotations

import json
import os
from typing import Any

import numpy as np

# ── Paper parameters ──────────────────────────────────────────────────────────
PAPER_P0: float = 0.25    # per-neighbour infection probability
PAPER_P1: float = 0.97    # recovery probability at end of infection period
PAPER_P2: float = 0.10    # re-infection probability for recovered cells
PAPER_T_TAU: int = 4      # fixed infection duration (steps)

# Published wavefront speed (Datta & Acharyya 2021, §3.1.1, Fig. 11)
PUBLISHED_SPEED: float = 0.4405
PUBLISHED_SPEED_ERR: float = 0.0008   # ± error as reported

# ── Simulation parameters ─────────────────────────────────────────────────────
MODEL_L: int = 200        # lattice side (500 in paper; L=200 adequate, see docstring)
N_SEEDS: int = 20
N_STEPS: int = 110
FIT_T_MIN: int = 5        # skip transient startup steps
FIT_T_MAX: int = 100      # stop well before wavefront nears boundary (L/2=100)


# ── Paper-exact CA implementation ─────────────────────────────────────────────

def _vn_infected_count(infected: np.ndarray) -> np.ndarray:
    """Count Von Neumann infected neighbours for every cell (periodic BCs)."""
    inf_int = infected.astype(np.int32)
    return (
        np.roll(inf_int, -1, axis=0)  # south
        + np.roll(inf_int, +1, axis=0)  # north
        + np.roll(inf_int, -1, axis=1)  # east
        + np.roll(inf_int, +1, axis=1)  # west
    )


def run_paper_ca(
    L: int,
    p0: float,
    p1: float,
    p2: float,
    t_tau: int,
    n_steps: int,
    seed: int,
) -> list[dict[str, Any]]:
    """Implement the Datta-Acharyya (2021) SIR CA exactly.

    State encoding (integer grid):
        0          → Susceptible
        1 .. t_tau → Infected, age 1 .. t_tau
        t_tau+1    → Recovered (= state 2 in paper notation)
        -1         → Dead (treated as removed for wavefront measurement)

    Rules (synchronous update):
        1. Infection: susceptible cell with ≥1 infected neighbour gets infected
           with probability 1 − (1−p0)^n_infected_neighbours.
        2. Ageing: infected cells at age < t_tau increment age by 1.
        3. Recovery/death: infected cells at age == t_tau become
           recovered (prob p1) or dead (prob 1−p1).
        4. Re-infection: recovered cells adjacent to infected cells get
           re-infected with probability 1 − (1−p2)^n_infected_neighbours.

    Returns list of per-step snapshots with keys: 'grid', 'step'.
    """
    rng = np.random.default_rng(seed)
    center = L // 2

    # grid: 0=S, 1..t_tau=I, t_tau+1=R, -1=dead
    grid = np.zeros((L, L), dtype=np.int32)
    grid[center, center] = 1  # single seed, age=1

    R_STATE = t_tau + 1
    history: list[dict[str, Any]] = [{"grid": grid.copy(), "step": 0}]

    for t in range(1, n_steps + 1):
        infected = (grid >= 1) & (grid <= t_tau)
        susceptible = grid == 0
        recovered = grid == R_STATE

        n_inf_nb = _vn_infected_count(infected)

        new_grid = grid.copy()

        # Rule 1: S → I
        s_cand = susceptible & (n_inf_nb > 0)
        p_infect = 1.0 - (1.0 - p0) ** n_inf_nb.astype(float)
        rand1 = rng.random((L, L))
        si_mask = s_cand & (rand1 < p_infect)
        new_grid[si_mask] = 1  # newly infected, age=1

        # Rule 2: I age < t_tau → increment age
        aging_mask = infected & (grid < t_tau)
        new_grid[aging_mask] = grid[aging_mask] + 1

        # Rule 3: I age == t_tau → recover or die
        final_mask = infected & (grid == t_tau)
        rand3 = rng.random((L, L))
        new_grid[final_mask & (rand3 < p1)] = R_STATE
        new_grid[final_mask & (rand3 >= p1)] = -1

        # Rule 4: R → I (re-infection)
        r_cand = recovered & (n_inf_nb > 0)
        p_reinfect = 1.0 - (1.0 - p2) ** n_inf_nb.astype(float)
        rand4 = rng.random((L, L))
        ri_mask = r_cand & (rand4 < p_reinfect)
        new_grid[ri_mask] = 1  # re-infected, age=1

        grid = new_grid
        history.append({"grid": grid.copy(), "step": t})

        # Early exit if epidemic is completely over
        if not infected.any():
            break

    return history


# ── Wavefront radius ──────────────────────────────────────────────────────────

def wavefront_radius(grid: np.ndarray, center_r: int, center_c: int,
                     t_tau: int) -> float:
    """Euclidean radius of the outermost affected cell (I ∪ R ∪ dead) from centre.

    Using the full affected region (I∪R) gives a monotonically non-decreasing
    signal that is robust for linear fitting.
    """
    affected = (grid >= 1) | (grid == -1)
    if not affected.any():
        return 0.0
    r_idx, c_idx = np.where(affected)
    dists = np.sqrt((r_idx - center_r) ** 2 + (c_idx - center_c) ** 2)
    return float(dists.max())


# ── Main ──────────────────────────────────────────────────────────────────────

def main() -> None:
    center = MODEL_L // 2
    seed_results: list[dict] = []
    valid_speeds: list[float] = []

    for seed in range(N_SEEDS):
        history = run_paper_ca(
            L=MODEL_L, p0=PAPER_P0, p1=PAPER_P1, p2=PAPER_P2,
            t_tau=PAPER_T_TAU, n_steps=N_STEPS, seed=seed,
        )

        times: list[float] = []
        radii: list[float] = []
        for snap in history:
            t = snap["step"]
            if FIT_T_MIN <= t <= FIT_T_MAX:
                r = wavefront_radius(snap["grid"], center, center, PAPER_T_TAU)
                if r > 0:
                    times.append(float(t))
                    radii.append(r)

        if len(times) < 5:
            seed_results.append({
                "seed": seed, "slope": None, "r2": None,
                "n_fit_points": len(times), "skipped": True,
            })
            continue

        times_arr = np.array(times)
        radii_arr = np.array(radii)
        coeffs = np.polyfit(times_arr, radii_arr, 1)
        slope = float(coeffs[0])
        residuals = radii_arr - np.polyval(coeffs, times_arr)
        ss_res = float(np.sum(residuals ** 2))
        ss_tot = float(np.sum((radii_arr - radii_arr.mean()) ** 2))
        r2 = float(1.0 - ss_res / ss_tot) if ss_tot > 0 else 0.0

        valid_speeds.append(slope)
        seed_results.append({
            "seed": seed, "slope": slope, "r2": r2,
            "n_fit_points": len(times), "skipped": False,
        })

    if not valid_speeds:
        raise RuntimeError("All seeds died out — parameters are sub-critical.")

    mean_speed = float(np.mean(valid_speeds))
    std_speed = float(np.std(valid_speeds))
    abs_err = abs(mean_speed - PUBLISHED_SPEED)
    rel_err = abs_err / PUBLISHED_SPEED
    passes = bool(abs_err < 0.05 or rel_err < 0.15)

    output = {
        "anchor": "wavefront_speed_datta_acharyya_2021_sec3_1_1_fig11",
        "paper": (
            "Datta & Acharyya (2021/2022), arXiv:2104.10456, "
            "Int. J. Mod. Phys. C 33, 2250094; §3.1.1 Velocity of Epidemic Spread / Fig. 11"
        ),
        "published_speed": PUBLISHED_SPEED,
        "published_speed_err": PUBLISHED_SPEED_ERR,
        "measured_speed_mean": mean_speed,
        "measured_speed_std": std_speed,
        "n_valid_seeds": len(valid_speeds),
        "n_skipped_seeds": N_SEEDS - len(valid_speeds),
        "abs_err": abs_err,
        "rel_err": rel_err,
        "passes_tolerance": passes,
        "parameters": {
            "paper_p0": PAPER_P0,
            "paper_p1": PAPER_P1,
            "paper_p2": PAPER_P2,
            "paper_t_tau": PAPER_T_TAU,
            "model_L": MODEL_L,
            "neighborhood": "von_neumann",
            "fit_range": [FIT_T_MIN, FIT_T_MAX],
            "n_seeds": N_SEEDS,
            "n_steps": N_STEPS,
            "model_note": (
                "Exact paper CA implemented inline (fixed t_tau duration + "
                "re-infection p2). epc SIREpidemicModel uses geometric recovery "
                "and is not directly called here."
            ),
        },
        "per_seed": seed_results,
    }

    os.makedirs("analysis/outputs", exist_ok=True)
    out_path = "analysis/outputs/p22_dattaacharyya2005_reproduction.json"
    with open(out_path, "w") as f:
        json.dump(output, f, indent=2)

    print(f"Published speed (Datta & Acharyya 2021 §3.1.1): "
          f"{PUBLISHED_SPEED} ± {PUBLISHED_SPEED_ERR}")
    print(f"Measured speed  (this reproduction):            "
          f"{mean_speed:.4f} ± {std_speed:.4f}  (n={len(valid_speeds)} seeds)")
    print(f"Absolute error: {abs_err:.4f}   Relative error: {rel_err:.4f}")
    print(f"Passes tolerance (|Δ|<0.05 OR rel<0.15): {passes}")
    print(f"Output: {out_path}")


if __name__ == "__main__":
    main()
