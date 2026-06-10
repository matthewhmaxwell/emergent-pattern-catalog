"""P4 dim2 multi-seed campaign: territoriality.

20 seeds; report exclusivity index, pairwise overlap, boundary persistence
as mean ± std ± CV. Emit analysis/outputs/p4_multiseed.json.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from epc.models.territoriality import ScentMarkingModel, ScentMarkingParams
from epc.detectors.p4_territoriality import (
    P4TerritorialityDetector,
    _exclusivity_index,
    _pairwise_overlap,
    _boundary_persistence,
)


def main() -> dict:
    """Run 20-seed P4 campaign and return summary statistics."""
    n_seeds = 20
    seeds = list(range(42, 42 + n_seeds))

    excls = []
    overlaps = []
    persistences = []
    tiers = []
    confidences = []

    for seed in seeds:
        params = ScentMarkingParams(
            n_agents=4, grid_size=48, n_steps=20000,
            snapshot_interval=100, seed=seed,
        )
        model = ScentMarkingModel(params)
        history = model.simulate()
        metadata = model.get_metadata()

        final_occ = history[-1]['occupancy']
        excl = _exclusivity_index(final_occ)
        overlap = _pairwise_overlap(final_occ)

        n_snap = len(history)
        quarter = max(1, n_snap // 4)
        pers = _boundary_persistence(history[quarter]['occupancy'], final_occ)

        detector = P4TerritorialityDetector(n_permutations=199, seed=seed)
        result = detector.detect(history, metadata)

        excls.append(excl)
        overlaps.append(overlap)
        persistences.append(pers)
        tiers.append(result.tier.value)
        confidences.append(result.confidence)

        print(f"  seed={seed}: excl={excl:.3f}, overlap={overlap:.4f}, "
              f"pers={pers:.3f}, tier={result.tier.value}")

    excl_arr = np.array(excls)
    ov_arr = np.array(overlaps)
    pers_arr = np.array(persistences)

    summary = {
        'pattern': 'P4',
        'n_seeds': n_seeds,
        'seeds': seeds,
        'exclusivity_index': {
            'mean': float(np.mean(excl_arr)),
            'std': float(np.std(excl_arr)),
            'cv': float(np.std(excl_arr) / np.mean(excl_arr)) if np.mean(excl_arr) > 0 else 0.0,
            'values': excls,
        },
        'pairwise_overlap': {
            'mean': float(np.mean(ov_arr)),
            'std': float(np.std(ov_arr)),
            'cv': float(np.std(ov_arr) / np.mean(ov_arr)) if np.mean(ov_arr) > 0 else 0.0,
            'values': overlaps,
        },
        'boundary_persistence': {
            'mean': float(np.mean(pers_arr)),
            'std': float(np.std(pers_arr)),
            'cv': float(np.std(pers_arr) / np.mean(pers_arr)) if np.mean(pers_arr) > 0 else 0.0,
            'values': persistences,
        },
        'tiers': tiers,
        'confidences': confidences,
        'fraction_definitive': sum(1 for t in tiers if t == 'definitive') / n_seeds,
        'fraction_detected': sum(1 for t in tiers if t in ('confirmation', 'definitive')) / n_seeds,
    }

    output_dir = Path(__file__).resolve().parent.parent / 'outputs'
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / 'p4_multiseed.json'
    with open(output_path, 'w') as f:
        json.dump(summary, f, indent=2)

    print(f"\nP4 dim2 multi-seed summary:")
    print(f"  Exclusivity: {np.mean(excl_arr):.4f} ± {np.std(excl_arr):.4f} "
          f"(CV={np.std(excl_arr)/np.mean(excl_arr)*100:.1f}%)")
    print(f"  Overlap: {np.mean(ov_arr):.4f} ± {np.std(ov_arr):.4f}")
    print(f"  Persistence: {np.mean(pers_arr):.4f} ± {np.std(pers_arr):.4f}")
    print(f"  Fraction definitive: {summary['fraction_definitive']:.1%}")
    print(f"  Written to {output_path}")

    return summary


if __name__ == '__main__':
    main()
