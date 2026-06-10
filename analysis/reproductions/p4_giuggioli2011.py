"""P4 dim1 reproduction: territoriality / exclusion boundaries.

Anchor: Giuggioli, Potts & Harris (2011) — territorial random walks
with scent marking produce stable, non-overlapping home ranges with
persistent boundaries.

Tolerance checks:
  1. Exclusivity index > 0.85 (strong single-agent cell dominance)
  2. Pairwise home-range overlap < 0.10 (exclusive territories)
  3. Boundary persistence > 0.70 (stable partition)
  4. Detector tier ∈ {confirmation, definitive}
"""

from __future__ import annotations

import json
from pathlib import Path

from epc.models.territoriality import ScentMarkingModel, ScentMarkingParams
from epc.detectors.p4_territoriality import (
    P4TerritorialityDetector,
    extract_observation_bundle,
    _exclusivity_index,
    _pairwise_overlap,
    _boundary_persistence,
)


def main() -> dict:
    """Run P4 dim1 reproduction and return results."""

    # ── Canonical regime ──
    params = ScentMarkingParams(
        n_agents=4, grid_size=48, n_steps=20000,
        snapshot_interval=100, seed=42,
    )
    model = ScentMarkingModel(params)
    history = model.simulate()
    metadata = model.get_metadata()

    # ── Extract metrics ──
    final_occ = history[-1]['occupancy']
    excl = _exclusivity_index(final_occ)
    overlap = _pairwise_overlap(final_occ)

    n_snap = len(history)
    quarter = max(1, n_snap // 4)
    early_occ = history[quarter]['occupancy']
    persistence = _boundary_persistence(early_occ, final_occ)

    # ── Run detector ──
    detector = P4TerritorialityDetector(n_permutations=199, seed=42)
    result = detector.detect(history, metadata)

    # ── Tolerance checks ──
    excl_pass = excl > 0.85
    overlap_pass = overlap < 0.10
    pers_pass = persistence > 0.70
    tier_pass = result.tier.value in ('confirmation', 'definitive')
    all_pass = excl_pass and overlap_pass and pers_pass and tier_pass

    results = {
        'pattern': 'P4',
        'anchor': 'Giuggioli, Potts & Harris (2011)',
        'model': 'ScentMarkingModel',
        'params': {
            'n_agents': params.n_agents,
            'grid_size': params.grid_size,
            'deposition_rate': params.deposition_rate,
            'decay_rate': params.decay_rate,
            'repulsion_strength': params.repulsion_strength,
            'home_attraction': params.home_attraction,
            'temperature': params.temperature,
            'n_steps': params.n_steps,
            'seed': params.seed,
        },
        'exclusivity_index': excl,
        'pairwise_overlap': overlap,
        'boundary_persistence': persistence,
        'detector_tier': result.tier.value,
        'detector_confidence': result.confidence,
        'null_p_value': result.null_p_value,
        'cohens_d': result.effect_size.get('cohens_d', 0.0),
        'occupancy_scent_correlation': result.primary_metric.get(
            'occupancy_scent_correlation', 0.0,
        ),
        'tolerances': {
            'exclusivity_gt_0.85': excl_pass,
            'overlap_lt_0.10': overlap_pass,
            'persistence_gt_0.70': pers_pass,
            'tier_confirmed_or_definitive': tier_pass,
        },
        'passes_tolerance': all_pass,
    }

    # ── Emit JSON ──
    output_dir = Path(__file__).resolve().parent.parent / 'outputs'
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / 'p4_giuggioli2011_reproduction.json'
    with open(output_path, 'w') as f:
        json.dump(results, f, indent=2)

    print(f"P4 dim1 reproduction:")
    print(f"  Exclusivity index = {excl:.4f} (> 0.85: {'PASS' if excl_pass else 'FAIL'})")
    print(f"  Pairwise overlap = {overlap:.4f} (< 0.10: {'PASS' if overlap_pass else 'FAIL'})")
    print(f"  Boundary persistence = {persistence:.4f} (> 0.70: {'PASS' if pers_pass else 'FAIL'})")
    print(f"  Detector tier = {result.tier.value} ({'PASS' if tier_pass else 'FAIL'})")
    print(f"  Overall: {'PASS' if all_pass else 'FAIL'}")
    print(f"  Written to {output_path}")

    return results


if __name__ == '__main__':
    main()
