"""
P17 dim2 — Multi-seed campaign at canonical regime.

20+ seeds at the canonical collective-sensing regime (N=50, Berdahl parameters).
Report: chemotactic index mean ± std ± CV.
Output: analysis/outputs/p17_multiseed.json
"""

import json
import sys
import os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

from epc.models.collective_sensing import CollectiveSensingModel, compute_chemotactic_index


def run_multiseed():
    """Run 20-seed campaign at canonical P17 regime."""
    box_size = 20.0
    N = 50
    n_seeds = 20
    seeds = list(range(n_seeds))
    n_steps = 1000
    record_interval = 5

    ci_values = []

    for seed in seeds:
        model = CollectiveSensingModel(
            n_agents=N,
            box_size=box_size,
            v_max=0.4,
            turn_noise=0.3,
            sensing_noise=0.8,
            alpha=0.95,
            social_strength=0.2,
            field_sigma=5.0,
            field_amplitude=1.0,
            dt=1.0,
            init_mode="offset",
            seed=seed,
        )
        history = model.run(n_steps=n_steps, record_interval=record_interval)
        ci = compute_chemotactic_index(history, model.field_center, box_size)
        ci_values.append(float(ci))

    ci_arr = np.array(ci_values)
    mean_ci = float(np.mean(ci_arr))
    std_ci = float(np.std(ci_arr))
    cv = std_ci / abs(mean_ci) if abs(mean_ci) > 1e-10 else float("inf")

    results = {
        "reproduction": "p17_multiseed",
        "pattern": "P17",
        "n_agents": N,
        "n_seeds": n_seeds,
        "n_steps": n_steps,
        "parameters": {
            "box_size": box_size,
            "v_max": 0.4,
            "turn_noise": 0.3,
            "sensing_noise": 0.8,
            "alpha": 0.95,
            "social_strength": 0.2,
            "field_sigma": 5.0,
        },
        "ci_values": ci_values,
        "ci_mean": mean_ci,
        "ci_std": std_ci,
        "ci_cv": cv,
        "all_positive": bool(np.all(ci_arr > 0)),
        "fraction_positive": float(np.mean(ci_arr > 0)),
    }

    output_dir = os.path.join(os.path.dirname(__file__), '..', 'outputs')
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, 'p17_multiseed.json')
    with open(output_path, 'w') as f:
        json.dump(results, f, indent=2)

    print(f"P17 multi-seed (N={N}, {n_seeds} seeds):")
    print(f"  CI mean: {mean_ci:.4f}")
    print(f"  CI std:  {std_ci:.4f}")
    print(f"  CI CV:   {cv:.2%}")
    print(f"  Fraction positive: {results['fraction_positive']:.2%}")
    print(f"Output: {output_path}")
    return results


if __name__ == "__main__":
    run_multiseed()
