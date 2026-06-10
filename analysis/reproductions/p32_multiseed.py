"""
P32 dim2 — multi-seed campaign for emergent specialization.

Runs ≥20 seeds of the canonical ResponseThresholdModel and reports
final per-agent entropy + efficiency-gain statistics.

Output: analysis/outputs/p32_multiseed.json
"""

import json
import sys
import os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

from epc.models.division_of_labor import ResponseThresholdModel
from epc.detectors.p32_specialization import (
    P32SpecializationDetector,
    extract_observation_bundle,
    _per_agent_entropy_window,
    _task_coverage,
)
from epc.detector_result import DetectionTier


def run_multiseed(n_seeds: int = 20) -> dict:
    """Run multi-seed campaign."""
    n_agents = 20
    n_tasks = 3
    reinforcement_rate = 0.05
    forgetting_rate = 0.01
    n_steps = 1000

    final_entropies = []
    efficiency_gains = []
    tiers = []
    confidences = []

    for seed in range(n_seeds):
        model = ResponseThresholdModel(
            n_agents=n_agents,
            n_tasks=n_tasks,
            reinforcement_rate=reinforcement_rate,
            forgetting_rate=forgetting_rate,
            seed=seed,
        )
        history = model.run(n_steps)

        bundle = extract_observation_bundle(history)
        assignments = bundle["task_assignments"]
        T = assignments.shape[0]
        window = max(T // 4, 50)

        late_entropy = _per_agent_entropy_window(assignments, n_tasks, T - window, T)
        mean_late_entropy = float(np.mean(late_entropy))
        final_entropies.append(mean_late_entropy)

        # Efficiency
        early_covered = 0
        for t in range(min(window, T)):
            tasks = set(assignments[t][assignments[t] >= 0])
            if len(tasks) >= n_tasks:
                early_covered += 1
        early_eff = early_covered / window

        late_covered = 0
        for t in range(T - window, T):
            tasks = set(assignments[t][assignments[t] >= 0])
            if len(tasks) >= n_tasks:
                late_covered += 1
        late_eff = late_covered / window

        efficiency_gains.append(late_eff - early_eff)

        # Run detector
        det = P32SpecializationDetector(n_permutations=199, seed=seed)
        r = det.detect(history, model.get_metadata())
        tiers.append(r.tier.value)
        confidences.append(r.confidence)

    final_entropies = np.array(final_entropies)
    efficiency_gains = np.array(efficiency_gains)

    result = {
        "pattern_id": "P32",
        "n_seeds": n_seeds,
        "parameters": {
            "n_agents": n_agents,
            "n_tasks": n_tasks,
            "reinforcement_rate": reinforcement_rate,
            "forgetting_rate": forgetting_rate,
            "n_steps": n_steps,
        },
        "final_entropy": {
            "mean": round(float(final_entropies.mean()), 4),
            "std": round(float(final_entropies.std()), 4),
            "cv": round(float(final_entropies.std() / max(final_entropies.mean(), 1e-6)), 4),
        },
        "efficiency_gain": {
            "mean": round(float(efficiency_gains.mean()), 4),
            "std": round(float(efficiency_gains.std()), 4),
            "cv": round(float(efficiency_gains.std() / max(abs(efficiency_gains.mean()), 1e-6)), 4),
        },
        "detection": {
            "tiers": tiers,
            "n_definitive": sum(1 for t in tiers if t == "definitive"),
            "n_confirmation": sum(1 for t in tiers if t == "confirmation"),
            "n_screening": sum(1 for t in tiers if t == "screening"),
            "n_detected": sum(1 for c in confidences if c > 0),
            "mean_confidence": round(float(np.mean(confidences)), 4),
        },
    }

    return result


if __name__ == "__main__":
    result = run_multiseed(n_seeds=20)

    out_dir = os.path.join(os.path.dirname(__file__), '..', 'outputs')
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, 'p32_multiseed.json')

    with open(out_path, 'w') as f:
        json.dump(result, f, indent=2)

    print(f"Multi-seed results ({result['n_seeds']} seeds):")
    print(f"  Final entropy: {result['final_entropy']['mean']:.4f} "
          f"± {result['final_entropy']['std']:.4f} "
          f"(CV={result['final_entropy']['cv']:.3f})")
    print(f"  Efficiency gain: {result['efficiency_gain']['mean']:.4f} "
          f"± {result['efficiency_gain']['std']:.4f} "
          f"(CV={result['efficiency_gain']['cv']:.3f})")
    print(f"  Detection tiers: {result['detection']['n_definitive']} definitive, "
          f"{result['detection']['n_confirmation']} confirmation, "
          f"{result['detection']['n_screening']} screening")
