"""P29 dim1 reproduction: Physarum/ant network efficiency.

Reproduce the network-efficiency signature quantitatively:
emergent network length within a documented factor of the MST/Steiner
optimum (length/MST in ~[1.0, 1.5]) with nonzero fault tolerance.
Cf. Tero et al. (2010), Science 327(5964), 439-442.

Emits analysis/outputs/p29_physarum_reproduction.json with passes_tolerance.
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
    """Run dim1 reproduction and emit JSON."""
    # Canonical parameters: 5-node grid, Physarum flux-reinforcement
    params = PhysarumParams(
        n_nodes=5,
        grid_size=80,
        gamma=1.8,
        decay_rate=0.01,
        n_steps=2000,
        snapshot_interval=50,
        node_layout="grid",
        seed=42,
    )
    model = PhysarumModel(params)
    history = model.simulate()

    # Extract final-snapshot metrics
    last = history[-1]
    node_pos = np.asarray(last['node_positions'])
    edge_weights = np.asarray(last['edge_weights'])
    mst_w = _mst_weight(node_pos)
    metrics = _network_metrics(edge_weights, node_pos)

    # Run detector
    detector = P29TrailNetworkDetector(n_permutations=199, seed=42)
    result = detector.detect(history, model.get_metadata())

    # Tolerance checks (Tero et al. 2010 result: length/MST ∈ [1.0, 1.5])
    length_ratio = metrics['length_ratio']
    fault_tolerance = metrics['fault_tolerance']
    weight_dist_corr = metrics['weight_dist_corr']

    passes_ratio = 1.0 <= length_ratio <= 1.5
    passes_fault_tolerance = fault_tolerance > 0.0
    passes_correlation = weight_dist_corr > 0.5
    passes_tolerance = passes_ratio and passes_fault_tolerance and passes_correlation

    output = {
        'pattern': 'P29',
        'model': 'PhysarumModel',
        'reference': 'Tero et al. (2010), Science 327(5964), 439-442',
        'anchor': 'Network length within [1.0, 1.5] × MST with fault tolerance > 0',
        'params': {
            'n_nodes': params.n_nodes,
            'gamma': params.gamma,
            'decay_rate': params.decay_rate,
            'n_steps': params.n_steps,
            'node_layout': params.node_layout,
            'seed': params.seed,
        },
        'mst_weight': round(mst_w, 4),
        'network_length': round(metrics['network_length'], 4),
        'length_ratio': round(length_ratio, 4),
        'fault_tolerance': round(fault_tolerance, 4),
        'weight_dist_corr': round(weight_dist_corr, 4),
        'concentration_gini': round(metrics['concentration'], 4),
        'n_strong_edges': metrics['n_strong_edges'],
        'connectivity': round(metrics['connectivity'], 4),
        'detector_tier': result.tier.value,
        'detector_confidence': round(result.confidence, 4),
        'detector_p_value': round(result.null_p_value, 4),
        'cohens_d': round(result.effect_size['cohens_d'], 4),
        'tolerance_length_ratio': '[1.0, 1.5]',
        'passes_ratio': passes_ratio,
        'passes_fault_tolerance': passes_fault_tolerance,
        'passes_correlation': passes_correlation,
        'passes_tolerance': passes_tolerance,
    }

    out_path = os.path.join(
        os.path.dirname(__file__), '..', 'outputs', 'p29_physarum_reproduction.json'
    )
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, 'w') as f:
        json.dump(output, f, indent=2)

    print(f"length/MST = {length_ratio:.4f}  (tolerance [1.0, 1.5]: {'PASS' if passes_ratio else 'FAIL'})")
    print(f"fault_tolerance = {fault_tolerance:.4f}  (> 0: {'PASS' if passes_fault_tolerance else 'FAIL'})")
    print(f"weight_dist_corr = {weight_dist_corr:.4f}  (> 0.5: {'PASS' if passes_correlation else 'FAIL'})")
    print(f"detector: {result.tier.value} (confidence={result.confidence:.2f}, p={result.null_p_value:.4f})")
    print(f"passes_tolerance: {passes_tolerance}")
    print(f"Written to {out_path}")


if __name__ == '__main__':
    main()
