"""
P16 dim2 — multi-seed stability at a fixed low-load α.

Runs ≥20 seeds at α = 0.05 (N=100, P=5) — the reliable retrieval
regime — and reports completion accuracy mean ± std ± CV.
Emits analysis/outputs/p16_multiseed.json.
"""

import json
import sys
import os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

from epc.models.hopfield import HopfieldNetwork, HopfieldParams
from epc.detectors.p16_associative_memory import P16AssociativeMemoryDetector


def run_multiseed() -> dict:
    """Run multi-seed campaign at canonical low-load α."""
    N = 100
    P = 5
    corruption = 0.2
    n_seeds = 20

    completion_accs = []
    recall_accs = []
    tiers = []
    confidences = []

    det = P16AssociativeMemoryDetector(n_permutations=199, seed=42)

    for s in range(n_seeds):
        seed = 100 + s
        params = HopfieldParams(N=N, P=P, corruption=corruption, seed=seed)
        model = HopfieldNetwork(params)
        history = model.simulate(n_cues=P, cue_pattern_indices=list(range(P)))

        result = det.detect(history, model.get_metadata())

        # Extract completion accuracy from result
        ca = result.primary_metric.get('mean_completion_accuracy', 0.0)
        ra = result.secondary_metrics.get('recall_accuracy', 0.0)

        completion_accs.append(ca)
        recall_accs.append(ra)
        tiers.append(result.tier.value)
        confidences.append(result.confidence)

    ca_arr = np.array(completion_accs)
    ra_arr = np.array(recall_accs)

    output = {
        'reproduction': 'p16_multiseed',
        'pattern': 'P16',
        'parameters': {
            'N': N,
            'P': P,
            'alpha': P / N,
            'corruption': corruption,
            'n_seeds': n_seeds,
        },
        'completion_accuracy': {
            'mean': float(np.mean(ca_arr)),
            'std': float(np.std(ca_arr)),
            'cv': float(np.std(ca_arr) / np.mean(ca_arr)) if np.mean(ca_arr) > 0 else 0.0,
            'min': float(np.min(ca_arr)),
            'max': float(np.max(ca_arr)),
            'values': [float(v) for v in ca_arr],
        },
        'recall_accuracy': {
            'mean': float(np.mean(ra_arr)),
            'std': float(np.std(ra_arr)),
            'values': [float(v) for v in ra_arr],
        },
        'detection': {
            'tiers': tiers,
            'confidences': [float(c) for c in confidences],
            'n_detected': sum(1 for t in tiers if t != 'screening' or confidences[tiers.index(t)] > 0),
            'n_definitive': sum(1 for t in tiers if t == 'definitive'),
            'n_confirmation': sum(1 for t in tiers if t == 'confirmation'),
        },
    }

    return output


if __name__ == '__main__':
    print("Running P16 multi-seed campaign...")
    output = run_multiseed()

    outpath = os.path.join(
        os.path.dirname(__file__), '..', 'outputs', 'p16_multiseed.json',
    )
    with open(outpath, 'w') as f:
        json.dump(output, f, indent=2)
    print(f"\nWritten to {outpath}")
    ca = output['completion_accuracy']
    print(f"Completion accuracy: {ca['mean']:.3f} ± {ca['std']:.3f} (CV={ca['cv']:.1%})")
    ra = output['recall_accuracy']
    print(f"Recall accuracy: {ra['mean']:.3f} ± {ra['std']:.3f}")
    d = output['detection']
    print(f"Detection: {d['n_definitive']} definitive, {d['n_confirmation']} confirmation")
