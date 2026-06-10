"""P29 dim2: multi-seed campaign for Physarum network efficiency.

Runs ≥20 seeds and computes length/MST ratio mean±std±CV.
Emits analysis/outputs/p29_multiseed.json.
"""

from __future__ import annotations

import json
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

import numpy as np

from epc.models.trail_network import PhysarumModel, PhysarumParams, _mst_weight
from epc.detectors.p29_trail_network import P29TrailNetworkDetector, _network_metrics


def main() -> None:
    """Run 20-seed campaign and emit JSON."""
    n_seeds = 20
    seeds = list(range(n_seeds))

    results = []
    for seed in seeds:
        params = PhysarumParams(
            n_nodes=5,
            grid_size=80,
            gamma=1.8,
            decay_rate=0.01,
            n_steps=2000,
            snapshot_interval=50,
            node_layout="random",
            seed=seed,
        )
        model = PhysarumModel(params)
        history = model.simulate()

        last = history[-1]
        node_pos = np.asarray(last['node_positions'])
        edge_weights = np.asarray(last['edge_weights'])
        metrics = _network_metrics(edge_weights, node_pos)

        detector = P29TrailNetworkDetector(n_permutations=199, seed=seed)
        det_result = detector.detect(history, model.get_metadata())

        results.append({
            'seed': seed,
            'length_ratio': round(metrics['length_ratio'], 4),
            'fault_tolerance': round(metrics['fault_tolerance'], 4),
            'weight_dist_corr': round(metrics['weight_dist_corr'], 4),
            'concentration': round(metrics['concentration'], 4),
            'connectivity': round(metrics['connectivity'], 4),
            'tier': det_result.tier.value,
            'confidence': round(det_result.confidence, 4),
            'p_value': round(det_result.null_p_value, 4),
            'detected': det_result.detected,
        })
        print(f"  seed {seed}: ratio={metrics['length_ratio']:.3f}, "
              f"corr={metrics['weight_dist_corr']:.3f}, "
              f"tier={det_result.tier.value}, "
              f"detected={det_result.detected}")

    # Aggregate statistics
    ratios = [r['length_ratio'] for r in results]
    corrs = [r['weight_dist_corr'] for r in results]
    detected_count = sum(1 for r in results if r['detected'])
    definitive_count = sum(1 for r in results if r['tier'] == 'definitive')

    ratio_mean = float(np.mean(ratios))
    ratio_std = float(np.std(ratios))
    ratio_cv = ratio_std / ratio_mean * 100 if ratio_mean > 0 else 0.0

    corr_mean = float(np.mean(corrs))
    corr_std = float(np.std(corrs))

    output = {
        'pattern': 'P29',
        'model': 'PhysarumModel',
        'n_seeds': n_seeds,
        'length_ratio_mean': round(ratio_mean, 4),
        'length_ratio_std': round(ratio_std, 4),
        'length_ratio_cv_pct': round(ratio_cv, 2),
        'weight_dist_corr_mean': round(corr_mean, 4),
        'weight_dist_corr_std': round(corr_std, 4),
        'detected_count': detected_count,
        'definitive_count': definitive_count,
        'detection_rate': round(detected_count / n_seeds, 4),
        'per_seed': results,
    }

    out_path = os.path.join(
        os.path.dirname(__file__), '..', 'outputs', 'p29_multiseed.json'
    )
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, 'w') as f:
        json.dump(output, f, indent=2)

    print(f"\nlength/MST: {ratio_mean:.4f} ± {ratio_std:.4f} (CV={ratio_cv:.1f}%)")
    print(f"weight-dist corr: {corr_mean:.4f} ± {corr_std:.4f}")
    print(f"Detected: {detected_count}/{n_seeds}, Definitive: {definitive_count}/{n_seeds}")
    print(f"Written to {out_path}")


if __name__ == '__main__':
    main()
