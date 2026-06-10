"""
P32 dim1 reproduction — Bonabeau, Theraulaz & Deneubourg (1996) response-threshold
division of labor.

Anchor:
    Bonabeau, E., Theraulaz, G. & Deneubourg, J.-L. (1996). Quantitative study
    of the fixed threshold model for the regulation of division of labour in
    insect societies. Proceedings of the Royal Society B, 263(1376), 1565–1569.

    Canonical result: identical agents with reinforcement-based threshold
    adaptation spontaneously differentiate into specialized roles. Per-agent
    task entropy declines from near-maximal to low while population-level task
    coverage is maintained, and collective efficiency exceeds the non-specialized
    baseline.

Reproduction targets:
    1. Entropy decline: mean per-agent late entropy < 0.5 * log2(n_tasks).
       Published: agents converge to near-exclusive task specialization.
       Tolerance: late entropy < 50% of maximum.

    2. Population task coverage maintained: all tasks have at least one
       specialist agent in the late window.
       Published: division of labor ensures all tasks are covered.
       Tolerance: late coverage = 1.0.

    3. Efficiency gain: late-window task-coverage efficiency exceeds
       early-window efficiency.
       Published: specialization improves collective performance.
       Tolerance: late_efficiency > early_efficiency.

Parameters:
    n_agents = 20
    n_tasks = 3
    reinforcement_rate = 0.15
    forgetting_rate = 0.01
    n_steps = 2000
    seed = 42
"""

import json
import sys
import os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

from epc.models.division_of_labor import (
    ResponseThresholdModel,
    compute_windowed_entropy,
    compute_efficiency,
)
from epc.detectors.p32_specialization import (
    extract_observation_bundle,
    _per_agent_entropy_window,
    _task_coverage,
)


def run_reproduction() -> dict:
    """Run the Bonabeau 1996 division-of-labor reproduction."""
    n_agents = 20
    n_tasks = 3
    reinforcement_rate = 0.15
    forgetting_rate = 0.01
    n_steps = 1000
    seed = 42

    model = ResponseThresholdModel(
        n_agents=n_agents,
        n_tasks=n_tasks,
        reinforcement_rate=reinforcement_rate,
        forgetting_rate=forgetting_rate,
        seed=seed,
    )
    history = model.run(n_steps)

    # Compute metrics
    bundle = extract_observation_bundle(history)
    assignments = bundle["task_assignments"]
    T = assignments.shape[0]
    window = max(T // 4, 50)

    early_entropy = _per_agent_entropy_window(assignments, n_tasks, 0, window)
    late_entropy = _per_agent_entropy_window(assignments, n_tasks, T - window, T)
    max_entropy = np.log2(n_tasks)

    mean_early = float(np.mean(early_entropy))
    mean_late = float(np.mean(late_entropy))
    entropy_decline = mean_early - mean_late

    late_coverage = _task_coverage(assignments, n_tasks, T - window, T)

    early_efficiency = compute_efficiency(history[:window], n_tasks)
    late_efficiency = compute_efficiency(history[T - window:], n_tasks)
    efficiency_gain = late_efficiency - early_efficiency

    # Tolerances
    entropy_pass = mean_late < 0.5 * max_entropy
    coverage_pass = late_coverage >= 1.0
    efficiency_pass = late_efficiency > early_efficiency

    passes_all = bool(entropy_pass and coverage_pass and efficiency_pass)

    result = {
        "pattern_id": "P32",
        "reference": "Bonabeau, Theraulaz & Deneubourg 1996",
        "parameters": {
            "n_agents": n_agents,
            "n_tasks": n_tasks,
            "reinforcement_rate": reinforcement_rate,
            "forgetting_rate": forgetting_rate,
            "n_steps": n_steps,
            "seed": seed,
        },
        "metrics": {
            "mean_early_entropy": round(mean_early, 4),
            "mean_late_entropy": round(mean_late, 4),
            "max_entropy": round(max_entropy, 4),
            "entropy_decline": round(entropy_decline, 4),
            "late_coverage": round(late_coverage, 4),
            "early_efficiency": round(early_efficiency, 4),
            "late_efficiency": round(late_efficiency, 4),
            "efficiency_gain": round(efficiency_gain, 4),
        },
        "tolerances": {
            "entropy": f"late < 0.5 * max ({0.5 * max_entropy:.3f}): {bool(entropy_pass)}",
            "coverage": f"late >= 1.0: {bool(coverage_pass)}",
            "efficiency": f"late > early: {bool(efficiency_pass)}",
        },
        "passes_tolerance": passes_all,
    }

    return result


if __name__ == "__main__":
    result = run_reproduction()

    out_dir = os.path.join(os.path.dirname(__file__), '..', 'outputs')
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, 'p32_bonabeau1996_reproduction.json')

    with open(out_path, 'w') as f:
        json.dump(result, f, indent=2)

    print(f"Reproduction result: passes_tolerance = {result['passes_tolerance']}")
    for k, v in result['metrics'].items():
        print(f"  {k}: {v}")
    for k, v in result['tolerances'].items():
        print(f"  {k}: {v}")

    if not result['passes_tolerance']:
        sys.exit(1)
