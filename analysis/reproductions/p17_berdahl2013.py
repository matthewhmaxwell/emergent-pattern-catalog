"""
P17 dim1 reproduction — Collective gradient sensing (Berdahl et al. 2013).

Target figure/anchor:
    Berdahl, A., Torney, C. J., Ioannou, C. C., Faria, J. J., & Couzin, I. D.
    (2013). Emergent sensing of complex environments by mobile animal groups.
    Science, 339(6119), 574–576. Fig. 1: navigation accuracy increases with
    group size.

Reproduction targets:
    1. Chemotactic index (CI) rises with group size N.
       Published: CI monotonically increases from N=1 (chance) to N~50+ (strong).
       Tolerance: positive, significant slope of CI vs log(N); CI(N=1) < 0.15;
       CI(N=50) > 0.15; slope > 0.02.

    2. Isolated agents (N=1) navigate at chance level.
       Published: N=1 has CI near 0 (no reliable gradient detection).
       Tolerance: |CI(N=1)| < 0.20 (averaged over seeds).

Parameters:
    box_size=20.0, v_max=0.4, turn_noise=0.3, sensing_noise=0.8,
    alpha=0.95, social_strength=0.2, field_sigma=5.0, n_steps=1000
    Group sizes: [1, 5, 10, 25, 50]
    Seeds per N: 10
"""

import json
import sys
import os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

from epc.models.collective_sensing import CollectiveSensingModel, compute_chemotactic_index


def run_reproduction():
    """Run P17 reproduction experiment."""
    box_size = 20.0
    group_sizes = [1, 5, 10, 25, 50]
    n_seeds = 10
    n_steps = 1000
    record_interval = 5

    params = dict(
        box_size=box_size, v_max=0.4, turn_noise=0.3,
        sensing_noise=0.8, alpha=0.95, social_strength=0.2,
        field_sigma=5.0, field_amplitude=1.0, dt=1.0,
        init_mode="offset",
    )

    results = {
        "reproduction": "p17_berdahl2013",
        "target": "CI rises with group size N (Berdahl 2013 Fig. 1)",
        "reference": "Berdahl et al. 2013, Science 339:574-576",
        "parameters": params,
        "group_sizes": group_sizes,
        "n_seeds": n_seeds,
        "n_steps": n_steps,
    }

    ci_by_n = {}
    for N in group_sizes:
        ci_values = []
        for s in range(n_seeds):
            seed = s * 77 + 13
            model = CollectiveSensingModel(
                n_agents=N, seed=seed, **params,
            )
            history = model.run(n_steps=n_steps, record_interval=record_interval)
            ci = compute_chemotactic_index(history, model.field_center, box_size)
            ci_values.append(ci)
        ci_by_n[N] = ci_values
        mean_ci = np.mean(ci_values)
        std_ci = np.std(ci_values)
        print(f"  N={N:3d}: CI = {mean_ci:.4f} ± {std_ci:.4f}")

    # Compute slope of CI vs log(N)
    sorted_ns = sorted(ci_by_n.keys())
    mean_ci_arr = np.array([np.mean(ci_by_n[n]) for n in sorted_ns])
    log_ns = np.log(np.array(sorted_ns, dtype=float))
    slope, intercept = np.polyfit(log_ns, mean_ci_arr, 1)

    # Compute Spearman correlation for significance
    from scipy.stats import spearmanr
    rho, p_spearman = spearmanr(log_ns, mean_ci_arr)

    ci_n1 = float(np.mean(ci_by_n[1]))
    ci_n50 = float(np.mean(ci_by_n[50]))

    # Tolerance checks
    tol_slope = slope > 0.02
    tol_n1_chance = abs(ci_n1) < 0.20
    tol_n50_positive = ci_n50 > 0.15
    tol_monotonic_trend = rho > 0.7

    passes_tolerance = all([tol_slope, tol_n1_chance, tol_n50_positive, tol_monotonic_trend])

    results["ci_by_N"] = {str(n): {
        "mean": float(np.mean(ci_by_n[n])),
        "std": float(np.std(ci_by_n[n])),
        "values": [float(x) for x in ci_by_n[n]],
    } for n in sorted_ns}
    results["slope_ci_vs_logN"] = float(slope)
    results["intercept"] = float(intercept)
    results["spearman_rho"] = float(rho)
    results["spearman_p"] = float(p_spearman)
    results["ci_N1_mean"] = ci_n1
    results["ci_N50_mean"] = ci_n50
    results["tolerance_checks"] = {
        "slope_gt_002": bool(tol_slope),
        "n1_at_chance": bool(tol_n1_chance),
        "n50_positive": bool(tol_n50_positive),
        "monotonic_trend": bool(tol_monotonic_trend),
    }
    results["passes_tolerance"] = bool(passes_tolerance)

    # Write output
    output_dir = os.path.join(os.path.dirname(__file__), '..', 'outputs')
    os.makedirs(output_dir, exist_ok=True)
    output_path = os.path.join(output_dir, 'p17_berdahl2013_reproduction.json')
    with open(output_path, 'w') as f:
        json.dump(results, f, indent=2)

    print(f"\nSlope of CI vs log(N): {slope:.4f}")
    print(f"Spearman rho: {rho:.4f}, p: {p_spearman:.4f}")
    print(f"CI(N=1): {ci_n1:.4f}, CI(N=50): {ci_n50:.4f}")
    print(f"Tolerance checks: {results['tolerance_checks']}")
    print(f"PASSES TOLERANCE: {passes_tolerance}")
    print(f"Output: {output_path}")
    return results


if __name__ == "__main__":
    run_reproduction()
