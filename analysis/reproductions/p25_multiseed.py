"""
P25 dim2 multi-seed campaign — Canalized restoration / equifinality.

Runs 20 seeds of the canonical CanalizedLandscape at fixed parameters
and reports convergence-variance-ratio mean ± std ± CV.

Parameters:
    CanalizedLandscape: n_dims=10, basin_strength=2.0, ic_spread=5.0,
    n_ics=20, n_steps=200, dt=0.05, seeds=0..19.
"""

import json
import sys
import os

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

from epc.models.canalization import CanalizedLandscape, CanalizedLandscapeParams
from epc.detectors.p25_equifinality import (
    P25EquifinalityDetector,
    extract_observation_bundle,
    _convergence_variance_ratio,
    _basin_volume,
)


def run_multiseed(n_seeds: int = 20) -> dict:
    """Run multi-seed campaign."""
    ratios = []
    basins = []
    tiers = []
    confidences = []

    for seed in range(n_seeds):
        params = CanalizedLandscapeParams(
            n_dims=10,
            n_ics=20,
            n_steps=200,
            basin_strength=2.0,
            ic_spread=5.0,
            dt=0.05,
            seed=seed,
        )
        model = CanalizedLandscape(params)
        h = model.simulate()

        bundle = extract_observation_bundle(h)
        ratio = _convergence_variance_ratio(bundle['ics'], bundle['finals'])
        basin = _basin_volume(bundle['converged'])

        det = P25EquifinalityDetector(n_permutations=199, seed=seed)
        r = det.detect(h, model.get_metadata())

        ratios.append(ratio)
        basins.append(basin)
        tiers.append(r.tier.value)
        confidences.append(r.confidence)

    ratios_arr = np.array(ratios)
    basins_arr = np.array(basins)

    ratio_mean = float(np.mean(ratios_arr))
    ratio_std = float(np.std(ratios_arr))
    ratio_cv = ratio_std / ratio_mean if ratio_mean > 0 else 0.0

    basin_mean = float(np.mean(basins_arr))
    basin_std = float(np.std(basins_arr))

    n_detected = sum(1 for t in tiers if t in ("screening", "confirmation", "definitive"))
    n_definitive = sum(1 for t in tiers if t == "definitive")

    output = {
        "campaign": "p25_multiseed",
        "pattern": "P25",
        "n_seeds": n_seeds,
        "parameters": {
            "n_dims": 10,
            "basin_strength": 2.0,
            "ic_spread": 5.0,
            "n_ics": 20,
            "n_steps": 200,
            "dt": 0.05,
        },
        "convergence_variance_ratio": {
            "mean": round(ratio_mean, 8),
            "std": round(ratio_std, 8),
            "cv": round(ratio_cv, 4),
            "values": [round(v, 8) for v in ratios],
        },
        "basin_volume": {
            "mean": round(basin_mean, 4),
            "std": round(basin_std, 4),
            "values": [round(v, 4) for v in basins],
        },
        "detection": {
            "n_detected": n_detected,
            "n_definitive": n_definitive,
            "fraction_detected": round(n_detected / n_seeds, 4),
            "tiers": tiers,
            "confidences": [round(c, 4) for c in confidences],
        },
    }

    return output


if __name__ == "__main__":
    result = run_multiseed()
    out_path = os.path.join(
        os.path.dirname(__file__), '..', 'outputs',
        'p25_multiseed.json',
    )
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, 'w') as f:
        json.dump(result, f, indent=2)
    print(json.dumps(result, indent=2))
    print(f"\nWritten to {out_path}")
