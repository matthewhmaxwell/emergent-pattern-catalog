"""
P30 dim2 multi-seed campaign: 20 seeds at canonical regime.

Output: analysis/outputs/p30_multiseed.json
"""

import sys
import os
import json
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

from epc.models.autopoiesis import AutopoiesisModel
from epc.detectors.p30_autopoiesis import P30AutopoiesisDetector


def run_multiseed():
    """Run 20-seed P30 campaign at canonical regime."""
    n_seeds = 20
    results = []

    for seed in range(n_seeds):
        m = AutopoiesisModel(
            n_substrate=100, n_catalyst=3, box_size=20.0,
            production_rate=0.15, decay_rate=0.01, link_attraction=0.3,
            membrane_equilibrium_radius=3.0, catalyst_link_attraction=0.5,
            production_radius=3.0, seed=seed,
        )
        h = m.run(150, record_every=5)
        d = P30AutopoiesisDetector(n_permutations=199, association_radius=4.0, seed=42)
        r = d.detect(h, m.get_metadata())

        results.append({
            'seed': seed,
            'detected': r.detected,
            'tier': r.tier.value,
            'confidence': float(r.confidence),
            'association_score': float(r.primary_metric['association_score']),
            'closure_fraction': float(r.primary_metric['closure_fraction']),
            'enrichment_ratio': float(r.primary_metric['enrichment_ratio']),
            'p_value': float(r.null_p_value),
        })
        print(f"  seed {seed}: {r.tier.value} (assoc={r.primary_metric['association_score']:.3f}, p={r.null_p_value:.4f})")

    # Summary
    assoc_vals = [r['association_score'] for r in results]
    closure_vals = [r['closure_fraction'] for r in results]
    enrich_vals = [r['enrichment_ratio'] for r in results]
    n_detected = sum(1 for r in results if r['detected'])
    tiers = [r['tier'] for r in results]
    n_confirmation = sum(1 for t in tiers if t in ('confirmation', 'definitive'))
    n_definitive = sum(1 for t in tiers if t == 'definitive')

    output = {
        'reproduction': 'p30_multiseed',
        'pattern': 'P30',
        'parameters': AutopoiesisModel(seed=0).get_metadata(),
        'n_seeds': n_seeds,
        'association_score': {
            'mean': float(np.mean(assoc_vals)),
            'std': float(np.std(assoc_vals)),
            'cv': float(np.std(assoc_vals) / np.mean(assoc_vals)) if np.mean(assoc_vals) > 0 else 0,
            'min': float(np.min(assoc_vals)),
            'max': float(np.max(assoc_vals)),
            'values': [float(v) for v in assoc_vals],
        },
        'closure_fraction': {
            'mean': float(np.mean(closure_vals)),
            'std': float(np.std(closure_vals)),
            'values': [float(v) for v in closure_vals],
        },
        'enrichment_ratio': {
            'mean': float(np.mean(enrich_vals)),
            'std': float(np.std(enrich_vals)),
            'values': [float(v) for v in enrich_vals],
        },
        'detection': {
            'n_detected': n_detected,
            'n_confirmation': n_confirmation,
            'n_definitive': n_definitive,
            'tiers': tiers,
            'confidences': [r['confidence'] for r in results],
        },
        'per_seed': results,
    }

    os.makedirs('analysis/outputs', exist_ok=True)
    with open('analysis/outputs/p30_multiseed.json', 'w') as f:
        json.dump(output, f, indent=2)

    print(f"\nP30 multi-seed (n={n_seeds}): assoc = {np.mean(assoc_vals):.3f} ± {np.std(assoc_vals):.3f} (CV={100*np.std(assoc_vals)/np.mean(assoc_vals):.1f}%)")
    print(f"  Detected: {n_detected}/{n_seeds}, Confirmation+: {n_confirmation}/{n_seeds}")
    return output


if __name__ == '__main__':
    run_multiseed()
