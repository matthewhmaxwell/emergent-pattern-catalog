"""
P19 dim2 — Multi-seed validation at canonical ρ=0.1.

Runs ≥20 seeds at the canonical informed fraction ρ=0.1 and reports
accuracy + influence-asymmetry statistics (mean ± std ± CV).

Output: analysis/outputs/p19_multiseed.json
"""

import json
import sys
import os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
from epc.models.informed_minority import InformedMinorityModel, group_accuracy
from epc.detectors.p19_emergent_leadership import detect_p19


def run_multiseed() -> dict:
    """Run 20-seed campaign at canonical ρ=0.1."""
    n_seeds = 20
    rho = 0.1
    n_particles = 200
    n_steps = 500

    accuracy_values = []
    te_asymmetry_values = []

    for seed in range(n_seeds):
        model = InformedMinorityModel(
            n_particles=n_particles,
            box_size=10.0,
            speed=0.03,
            noise=0.1,
            interaction_radius=1.0,
            informed_fraction=rho,
            bias_weight=0.3,
            preferred_direction=0.0,
            seed=seed,
        )
        history = model.run(n_steps)
        metadata = model.get_metadata()

        # Accuracy
        warmup = len(history) * 3 // 4
        accs = []
        for t in range(warmup, len(history)):
            accs.append(group_accuracy(history[t], 0.0))
        accuracy_values.append(float(np.mean(accs)))

        # Detector run for influence asymmetry
        result = detect_p19(history, metadata, n_permutations=99, seed=seed)
        pull = result.secondary_metrics.get("pull_mean", 0.0)
        te_asymmetry_values.append(float(pull))

        print(f"  seed {seed:2d}: accuracy={accuracy_values[-1]:.4f}, "
              f"pull={pull:.4f}, tier={result.tier}")

    acc_arr = np.array(accuracy_values)
    te_arr = np.array(te_asymmetry_values)

    output = {
        "reproduction": "p19_multiseed",
        "pattern": "P19",
        "informed_fraction": rho,
        "n_particles": n_particles,
        "n_seeds": n_seeds,
        "n_steps": n_steps,
        "parameters": {
            "box_size": 10.0,
            "speed": 0.03,
            "noise": 0.1,
            "interaction_radius": 1.0,
            "bias_weight": 0.3,
            "preferred_direction": 0.0,
        },
        "accuracy_values": [float(x) for x in accuracy_values],
        "accuracy_mean": float(np.mean(acc_arr)),
        "accuracy_std": float(np.std(acc_arr)),
        "accuracy_cv": float(np.std(acc_arr) / np.mean(acc_arr)) if np.mean(acc_arr) > 0 else 0.0,
        "influence_pull_values": [float(x) for x in te_asymmetry_values],
        "influence_pull_mean": float(np.mean(te_arr)),
        "influence_pull_std": float(np.std(te_arr)),
        "influence_pull_cv": float(np.std(te_arr) / abs(np.mean(te_arr))) if abs(np.mean(te_arr)) > 1e-8 else float('inf'),
        "fraction_detected": float(np.mean(np.array(accuracy_values) > 0.3)),
    }
    return output


if __name__ == "__main__":
    print("P19 multi-seed validation: ρ=0.1, 20 seeds")
    print("=" * 60)

    result = run_multiseed()

    outpath = os.path.join(
        os.path.dirname(__file__), '..', 'outputs',
        'p19_multiseed.json',
    )
    with open(outpath, 'w') as f:
        json.dump(result, f, indent=2)

    print(f"\nAccuracy: {result['accuracy_mean']:.4f} ± {result['accuracy_std']:.4f} "
          f"(CV={result['accuracy_cv']*100:.1f}%)")
    print(f"Influence pull: {result['influence_pull_mean']:.4f} ± {result['influence_pull_std']:.4f} "
          f"(CV={result['influence_pull_cv']*100:.1f}%)")
    print(f"Fraction detected: {result['fraction_detected']:.0%}")
    print(f"Output saved to {outpath}")
