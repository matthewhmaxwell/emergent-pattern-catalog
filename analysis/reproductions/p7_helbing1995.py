"""
P7 dim1 reproduction — Lane formation in counterflow (Helbing & Molnár 1995).

Target figure/anchor:
    Helbing, D. & Molnár, P. (1995). Social force model for pedestrian dynamics.
    Physical Review E, 51(5), 4282–4286.

    Helbing, D., Farkas, I. J., & Vicsek, T. (2000). Simulating dynamical features
    of escape panic. Nature, 407, 487–490. (Counterflow efficiency gain.)

    Nowak, S. & Schadschneider, A. (2012). Quantitative analysis of pedestrian
    counterflow in a cellular automaton model. Phys. Rev. E, 85(6), 066128.

Reproduction targets:
    1. Lane formation from mixed initial: Lane order parameter φ rises from
       initial mixing (φ_init ≈ random baseline) to φ_final ≥ 0.6 at moderate
       density (canonical: ρ=2.5 agents/m²). Published: φ saturates in [0.6, 0.95].
       Tolerance: φ_final ∈ [0.5, 1.0] AND φ_final > φ_init + 0.2.

    2. Throughput gain: Mean speed in desired direction is higher in steady state
       (lanes formed) vs early transient (mixed). Published: counterflow with
       lanes achieves 50-80% of free-flow speed vs 30-50% when mixed.
       Tolerance: throughput_final / v0 > 0.4 AND throughput_final > throughput_early.

Parameters:
    Corridor: W=20, H=4 (80 m²)
    N = 200 (density ρ=2.5/m²)
    v0 = 1.0 m/s, A = 5.0 N, B = 0.3 m, tau = 0.5 s, dt = 0.05 s
    n_steps = 3000 (150 time units), record_interval = 20
    5 seeds
"""

import json
import sys
import os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

from epc.models.lane_formation import LaneFormationModel
from epc.detectors.p7_lane_formation import (
    compute_lane_order_parameter,
    compute_throughput,
)


def run_reproduction():
    """Run P7 reproduction experiment."""
    W, H = 20.0, 4.0
    N = 200
    density = N / (W * H)
    seeds = [42, 123, 456, 789, 1011]
    n_steps = 3000
    record_interval = 20
    n_bins = 10

    results = {
        "reproduction": "p7_helbing1995",
        "target": "Lane formation from mixed initial condition in counterflow",
        "reference": "Helbing & Molnar 1995; Nowak & Schadschneider 2012",
        "parameters": {
            "corridor_width": W,
            "corridor_height": H,
            "n_agents": N,
            "density": density,
            "desired_speed": 1.0,
            "repulsion_amplitude": 5.0,
            "repulsion_range": 0.3,
            "tau": 0.5,
            "dt": 0.05,
            "n_steps": n_steps,
            "n_bins": n_bins,
            "seeds": seeds,
        },
        "seed_results": [],
    }

    phi_finals = []
    phi_inits = []
    throughput_finals = []
    throughput_inits = []

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

        # Early window (first 10% of history)
        early_end = max(2, len(history) // 10)
        early_states = history[1:early_end]  # skip t=0 (trivially initialized)

        # Steady state window (last 25%)
        late_start = len(history) * 3 // 4
        late_states = history[late_start:]

        # Lane order: early vs late
        phi_early = []
        for state in early_states:
            phi = compute_lane_order_parameter(
                state["positions"], state["labels"], H, n_bins
            )
            phi_early.append(phi)

        phi_late = []
        for state in late_states:
            phi = compute_lane_order_parameter(
                state["positions"], state["labels"], H, n_bins
            )
            phi_late.append(phi)

        # Throughput: early vs late
        tp_early = [compute_throughput(s["velocities"], s["labels"]) for s in early_states]
        tp_late = [compute_throughput(s["velocities"], s["labels"]) for s in late_states]

        phi_init = float(np.mean(phi_early))
        phi_final = float(np.mean(phi_late))
        tp_init = float(np.mean(tp_early))
        tp_final = float(np.mean(tp_late))

        phi_inits.append(phi_init)
        phi_finals.append(phi_final)
        throughput_inits.append(tp_init)
        throughput_finals.append(tp_final)

        results["seed_results"].append({
            "seed": seed,
            "phi_init": phi_init,
            "phi_final": phi_final,
            "phi_gain": phi_final - phi_init,
            "throughput_init": tp_init,
            "throughput_final": tp_final,
            "throughput_fraction": tp_final / 1.0,  # divide by v0
        })

    # Aggregate
    mean_phi_init = float(np.mean(phi_inits))
    mean_phi_final = float(np.mean(phi_finals))
    std_phi_final = float(np.std(phi_finals))
    mean_tp_init = float(np.mean(throughput_inits))
    mean_tp_final = float(np.mean(throughput_finals))
    std_tp_final = float(np.std(throughput_finals))

    # Tolerance checks
    # 1. Lane order rises to ≥ 0.5 AND gains ≥ 0.2 over initial
    tol1_level = mean_phi_final >= 0.5
    tol1_gain = (mean_phi_final - mean_phi_init) >= 0.2
    tol1_pass = tol1_level and tol1_gain

    # 2. Throughput: final > 0.4 * v0 AND final > init
    tol2_level = mean_tp_final >= 0.4
    tol2_gain = mean_tp_final > mean_tp_init
    tol2_pass = tol2_level and tol2_gain

    passes_tolerance = tol1_pass and tol2_pass

    results["aggregate"] = {
        "phi_init_mean": mean_phi_init,
        "phi_final_mean": mean_phi_final,
        "phi_final_std": std_phi_final,
        "phi_gain": mean_phi_final - mean_phi_init,
        "throughput_init_mean": mean_tp_init,
        "throughput_final_mean": mean_tp_final,
        "throughput_final_std": std_tp_final,
        "throughput_fraction_of_v0": mean_tp_final / 1.0,
    }

    results["tolerances"] = {
        "phi_final_threshold": 0.5,
        "phi_final_pass": bool(tol1_level),
        "phi_gain_threshold": 0.2,
        "phi_gain_pass": bool(tol1_gain),
        "throughput_threshold": 0.4,
        "throughput_level_pass": bool(tol2_level),
        "throughput_gain_pass": bool(tol2_gain),
    }
    results["passes_tolerance"] = bool(passes_tolerance)

    # Write output
    output_path = os.path.join(
        os.path.dirname(__file__), '..', 'outputs',
        'p7_helbing1995_reproduction.json'
    )
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    with open(output_path, 'w') as f:
        json.dump(results, f, indent=2)

    print(f"Reproduction {'PASS' if passes_tolerance else 'FAIL'}")
    print(f"  φ_init = {mean_phi_init:.3f}, φ_final = {mean_phi_final:.3f} ± {std_phi_final:.3f}")
    print(f"    Level ≥ 0.5: {'PASS' if tol1_level else 'FAIL'}")
    print(f"    Gain ≥ 0.2: {'PASS' if tol1_gain else 'FAIL'} (gain={mean_phi_final - mean_phi_init:.3f})")
    print(f"  Throughput: init={mean_tp_init:.3f}, final={mean_tp_final:.3f} ± {std_tp_final:.3f}")
    print(f"    Level ≥ 0.4·v0: {'PASS' if tol2_level else 'FAIL'}")
    print(f"    Gain > init: {'PASS' if tol2_gain else 'FAIL'}")

    return results


if __name__ == "__main__":
    run_reproduction()
