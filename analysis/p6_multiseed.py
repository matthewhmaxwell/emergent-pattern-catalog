"""
P6 (D'Orsogna milling) multi-seed extension — Sprint 56 dim2 closure.

Runs N_SEEDS independent random initializations of the D'Orsogna SPP model
at canonical milling parameters, computes mean normalised angular momentum
|L| per seed, and reports aggregate (mean, std, CV).

Parameters: N=100, C_a=0.5, C_r=1.0, l_a=3.0, l_r=0.5, α=1.0, β=0.5,
dt=0.05, init_mode="random" (seed-dependent IC), n_steps=3000.

init_mode="random" is used so that each seed genuinely provides a
different initial configuration — particles are drawn from uniform
random positions in [-R, R]² with random headings. This tests that the
milling attractor is reliably reached from diverse initial conditions.
With "ring" init the IC is deterministic (angle-spaced), so seeds would
all give identical dynamics.

Observable: |L| (normalised angular momentum, frame-averaged over the
steady-state window, steps 1501–3000).  |L| ≈ 1 → milling; |L| ≈ 0 →
no rotational order.

Output: analysis/outputs/p6_multiseed.json
"""

from __future__ import annotations

import json
import time
from typing import Any, Dict, List

import numpy as np

from epc.models.dorsogna_spp import DOrsognaSPPModel
from epc.metrics.collective_motion import AngularMomentumMetric

# --- Parameters ---
N_PARTICLES: int = 100
C_A: float = 0.5
C_R: float = 1.0
L_A: float = 3.0
L_R: float = 0.5
ALPHA: float = 1.0
BETA: float = 0.5
DT: float = 0.05
INIT_MODE: str = "random"
INIT_RADIUS: float = 5.0
WARMUP_STEPS: int = 2500       # discard initial transient (mill forms by ~2000 steps)
MEASURE_STEPS: int = 500       # measurement window after warmup
MEASURE_EVERY: int = 100       # snapshot stride → 5 snapshots per seed

N_SEEDS: int = 20
SEED_BASE: int = 100           # seeds 100–119; avoids panel seeds 0–4


def run_seed(seed: int) -> Dict[str, Any]:
    """Run one D'Orsogna seed; return mean |L| over steady-state window.

    Two-phase protocol:
      1. Warmup: WARMUP_STEPS steps with no recording (mill forms).
      2. Measure: MEASURE_STEPS steps, snapshot every MEASURE_EVERY steps.
    Uses model.setup() + manual step() loop to separate warmup from measurement
    (model.run() calls setup() internally, so split phases use the raw step API).
    """
    model = DOrsognaSPPModel(
        n_particles=N_PARTICLES,
        C_a=C_A,
        C_r=C_R,
        l_a=L_A,
        l_r=L_R,
        alpha=ALPHA,
        beta=BETA,
        dt=DT,
        init_mode=INIT_MODE,
        init_radius=INIT_RADIUS,
        seed=seed,
    )
    # Phase 1: warmup (no recording)
    model.setup()
    for _ in range(WARMUP_STEPS):
        model.step()
    # Phase 2: measurement window
    snapshots = []
    for t in range(1, MEASURE_STEPS + 1):
        model.step()
        if t % MEASURE_EVERY == 0:
            snapshots.append(model.get_state())
    L_series = np.array([
        AngularMomentumMetric.compute_instant(s, box_size=None)
        for s in snapshots
    ])
    L_abs_mean = float(np.mean(np.abs(L_series)))
    L_abs_std = float(np.std(np.abs(L_series), ddof=1)) if len(L_series) > 1 else 0.0
    return {
        "seed": seed,
        "L_abs_mean": L_abs_mean,
        "L_abs_std": L_abs_std,
        "n_frames": len(snapshots),
    }


def main() -> None:
    t0 = time.time()
    print(
        f"P6 D'Orsogna milling multi-seed: N={N_PARTICLES}, "
        f"C_a={C_A}, C_r={C_R}, l_a={L_A}, l_r={L_R}, "
        f"α={ALPHA}, β={BETA}, dt={DT}, init_mode={INIT_MODE!r}, "
        f"warmup={WARMUP_STEPS}, measure={MEASURE_STEPS}, stride={MEASURE_EVERY}, N_seeds={N_SEEDS}"
    )

    per_seed: List[Dict[str, Any]] = []
    for i in range(N_SEEDS):
        seed = SEED_BASE + i
        t_seed = time.time()
        result = run_seed(seed)
        elapsed = time.time() - t_seed
        per_seed.append(result)
        print(
            f"  seed {seed:3d}: |L|={result['L_abs_mean']:.4f} "
            f"± {result['L_abs_std']:.4f}  ({elapsed:.1f}s)"
        )

    L_vals = np.array([r["L_abs_mean"] for r in per_seed])
    L_mean = float(np.mean(L_vals))
    L_std = float(np.std(L_vals, ddof=1))
    L_cv = float(L_std / L_mean) if L_mean > 0 else float("nan")

    print(f"\n|L| aggregate ({N_SEEDS} seeds):")
    print(f"  mean = {L_mean:.4f}")
    print(f"  std  = {L_std:.4f}")
    print(f"  CV   = {L_cv:.4f}")

    output: Dict[str, Any] = {
        "sprint": 56,
        "description": "P6 D'Orsogna milling multi-seed |L| extension",
        "params": {
            "n_particles": N_PARTICLES,
            "C_a": C_A,
            "C_r": C_R,
            "l_a": L_A,
            "l_r": L_R,
            "alpha": ALPHA,
            "beta": BETA,
            "dt": DT,
            "init_mode": INIT_MODE,
            "init_radius": INIT_RADIUS,
            "warmup_steps": WARMUP_STEPS,
            "measure_steps": MEASURE_STEPS,
            "measure_every": MEASURE_EVERY,
            "n_seeds": N_SEEDS,
            "seed_base": SEED_BASE,
        },
        "primary_observable": "|L| (normalised angular momentum, frame-averaged over steady-state window)",
        "per_seed": per_seed,
        "aggregate": {
            "L_abs_mean": L_mean,
            "L_abs_std": L_std,
            "L_abs_cv": L_cv,
            "L_abs_min": float(L_vals.min()),
            "L_abs_max": float(L_vals.max()),
        },
        "elapsed_seconds": time.time() - t0,
    }

    out_path = "analysis/outputs/p6_multiseed.json"
    with open(out_path, "w") as f:
        json.dump(output, f, indent=2)
    print(f"\nSaved → {out_path}  (elapsed {output['elapsed_seconds']:.1f}s)")


if __name__ == "__main__":
    main()
