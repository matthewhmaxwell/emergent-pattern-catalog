"""
P7 dim2 — Multi-seed campaign at canonical regime.

20+ seeds at the canonical lane-forming regime (N=200, ρ=2.5/m², corridor 20×4).
Report: lane order parameter mean ± std ± CV.
Output: analysis/outputs/p7_multiseed.json
"""

import json
import sys
import os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

from epc.models.lane_formation import LaneFormationModel
from epc.detectors.p7_lane_formation import compute_lane_order_parameter


def run_multiseed():
    """Run 20-seed campaign at canonical P7 regime."""
    W, H = 20.0, 4.0
    N = 200
    n_seeds = 20
    seeds = list(range(n_seeds))
    n_steps = 2000
    record_interval = 20
    n_bins = 10

    phi_values = []

    for seed in seeds:
        model = LaneFormationModel(
            n_agents=N,
            corridor_width=W,
            corridor_height=H,
            desired_speed=1.0,
            repulsion_amplitude=5.0,
            repulsion_range=0.3,
            tau=0.5,
            dt=0.05,
            seed=seed,
        )
        history = model.run(n_steps=n_steps, record_interval=record_interval)

        # Steady state: last 25% of history
        late_start = len(history) * 3 // 4
        phi_late = []
        for state in history[late_start:]:
            phi = compute_lane_order_parameter(
                state["positions"], state["labels"], H, n_bins
            )
            phi_late.append(phi)
        phi_values.append(float(np.mean(phi_late)))

    phi_arr = np.array(phi_values)
    mean_phi = float(np.mean(phi_arr))
    std_phi = float(np.std(phi_arr))
    cv = std_phi / mean_phi if mean_phi > 0 else float("inf")

    results = {
        "campaign": "p7_multiseed",
        "pattern": "P7",
        "metric": "lane_order_parameter",
        "n_seeds": n_seeds,
        "parameters": {
            "n_agents": N,
            "corridor_width": W,
            "corridor_height": H,
            "density": N / (W * H),
            "desired_speed": 1.0,
            "repulsion_amplitude": 5.0,
            "repulsion_range": 0.3,
            "n_steps": n_steps,
        },
        "results": {
            "mean": mean_phi,
            "std": std_phi,
            "cv": cv,
            "min": float(np.min(phi_arr)),
            "max": float(np.max(phi_arr)),
            "n_valid": int(np.sum(phi_arr > 0.5)),
            "per_seed": phi_values,
        },
    }

    output_path = os.path.join(
        os.path.dirname(__file__), '..', 'outputs', 'p7_multiseed.json'
    )
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, 'w') as f:
        json.dump(results, f, indent=2)

    print(f"P7 multi-seed (n={n_seeds}): φ = {mean_phi:.4f} ± {std_phi:.4f} (CV={cv:.3f})")
    print(f"  Range: [{np.min(phi_arr):.3f}, {np.max(phi_arr):.3f}]")
    print(f"  Seeds with φ > 0.5: {np.sum(phi_arr > 0.5)}/{n_seeds}")

    return results


if __name__ == "__main__":
    run_multiseed()
