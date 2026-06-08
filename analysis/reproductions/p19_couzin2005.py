"""
P19 dim1 reproduction — Emergent leadership (Couzin et al. 2005).

Target figure/anchor:
    Couzin, I. D., Krause, J., Franks, N. R., & Levin, S. A. (2005).
    Effective leadership and decision-making in animal groups on the move.
    Nature, 433(7025), 513–516. Fig. 2a: group directional accuracy
    increases with informed fraction ρ, with diminishing returns.

Reproduction targets:
    1. Accuracy rises monotonically with ρ.
       Published: accuracy(ρ=0.05) > chance, accuracy(ρ=0.5) ≈ 1.0.
       Tolerance: Spearman correlation ρ > 0.8 across sweep; monotone.

    2. Small informed fraction achieves high accuracy (diminishing returns).
       Published: even ρ=0.05 achieves moderate accuracy; ρ≈0.15 near maximum.
       Tolerance: accuracy(ρ=0.1) > 0.5; accuracy(ρ=0.05) > 0.2.

    3. ρ=0 produces chance-level alignment with preferred direction.
       Tolerance: |accuracy(ρ=0)| < 0.3 (random).

Parameters:
    N=200, box_size=10.0, speed=0.03, noise=0.1, interaction_radius=1.0,
    bias_weight=0.3, preferred_direction=0.0, n_steps=500
    Informed fractions: [0.0, 0.025, 0.05, 0.1, 0.15, 0.2, 0.3, 0.5]
    Seeds per ρ: 10
"""

import json
import sys
import os
import numpy as np
from scipy.stats import spearmanr

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))
from epc.models.informed_minority import InformedMinorityModel, group_accuracy


def run_reproduction() -> dict:
    """Run the Couzin 2005 accuracy-vs-ρ reproduction."""
    rho_values = [0.0, 0.025, 0.05, 0.1, 0.15, 0.2, 0.3, 0.5]
    n_seeds = 10
    n_steps = 500
    n_particles = 200

    results_by_rho = {}

    for rho in rho_values:
        accuracies = []
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
            # Compute accuracy over last 25% of history
            warmup = len(history) * 3 // 4
            accs = []
            for t in range(warmup, len(history)):
                accs.append(group_accuracy(history[t], 0.0))
            accuracies.append(float(np.mean(accs)))

        results_by_rho[str(rho)] = {
            "rho": rho,
            "accuracies": accuracies,
            "mean": float(np.mean(accuracies)),
            "std": float(np.std(accuracies)),
        }
        print(f"  ρ={rho:.3f}: accuracy = {np.mean(accuracies):.3f} ± {np.std(accuracies):.3f}")

    # Tolerance checks
    mean_values = [results_by_rho[str(r)]["mean"] for r in rho_values]
    rho_nonzero = [r for r in rho_values if r > 0]
    mean_nonzero = [results_by_rho[str(r)]["mean"] for r in rho_nonzero]

    # 1. Spearman correlation: accuracy should rise with ρ
    rho_corr, p_corr = spearmanr(rho_values, mean_values)
    check_monotone = rho_corr > 0.8

    # 2. Small ρ achieves moderate accuracy
    acc_01 = results_by_rho["0.1"]["mean"]
    check_small_rho = acc_01 > 0.5

    acc_005 = results_by_rho["0.05"]["mean"]
    check_tiny_rho = acc_005 > 0.2

    # 3. ρ=0 is chance-level
    acc_0 = results_by_rho["0.0"]["mean"]
    check_chance = abs(acc_0) < 0.3

    # 4. ρ=0.5 near maximum
    acc_05 = results_by_rho["0.5"]["mean"]
    check_high_rho = acc_05 > 0.8

    passes_tolerance = bool(all([
        check_monotone,
        check_small_rho,
        check_tiny_rho,
        check_chance,
        check_high_rho,
    ]))

    output = {
        "reproduction": "p19_couzin2005",
        "pattern": "P19",
        "reference": "Couzin et al. (2005) Nature 433:513–516, Fig. 2a",
        "n_particles": n_particles,
        "n_steps": n_steps,
        "n_seeds": n_seeds,
        "parameters": {
            "box_size": 10.0,
            "speed": 0.03,
            "noise": 0.1,
            "interaction_radius": 1.0,
            "bias_weight": 0.3,
            "preferred_direction": 0.0,
        },
        "rho_values": rho_values,
        "results_by_rho": results_by_rho,
        "tolerance_checks": {
            "monotone_spearman_rho": float(rho_corr),
            "monotone_spearman_p": float(p_corr),
            "check_monotone": bool(check_monotone),
            "accuracy_rho_0.1": float(acc_01),
            "check_small_rho": bool(check_small_rho),
            "accuracy_rho_0.05": float(acc_005),
            "check_tiny_rho": bool(check_tiny_rho),
            "accuracy_rho_0.0": float(acc_0),
            "check_chance": bool(check_chance),
            "accuracy_rho_0.5": float(acc_05),
            "check_high_rho": bool(check_high_rho),
        },
        "passes_tolerance": passes_tolerance,
    }
    return output


if __name__ == "__main__":
    print("P19 Couzin 2005 reproduction: accuracy vs informed fraction ρ")
    print("=" * 60)

    result = run_reproduction()

    outpath = os.path.join(
        os.path.dirname(__file__), '..', 'outputs',
        'p19_couzin2005_reproduction.json',
    )
    with open(outpath, 'w') as f:
        json.dump(result, f, indent=2)

    print(f"\nTolerance checks:")
    for k, v in result["tolerance_checks"].items():
        print(f"  {k}: {v}")

    print(f"\nOverall PASS: {result['passes_tolerance']}")
    print(f"Output saved to {outpath}")
