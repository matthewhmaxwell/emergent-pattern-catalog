"""P20 dim2 multi-seed campaign: quorum sensing.

20 seeds; report critical density, hysteresis width, step-function R²
as mean ± std ± CV. Emit analysis/outputs/p20_multiseed.json.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from epc.models.quorum_sensing import AutoinducerQuorum, AutoinducerParams
from epc.detectors.p20_quorum_sensing import (
    P20QuorumSensingDetector,
    extract_observation_bundle,
    _compute_equilibrium_curve,
    _best_step_function_r2,
    _compute_hysteresis,
)


def main() -> dict:
    """Run 20-seed P20 campaign and return summary statistics."""
    n_seeds = 20
    seeds = list(range(42, 42 + n_seeds))

    r2s = []
    crit_dens = []
    hyst_widths = []
    tiers = []
    confidences = []

    for seed in seeds:
        params = AutoinducerParams(seed=seed)
        model = AutoinducerQuorum(params)
        history = model.simulate()
        metadata = model.get_metadata()

        bundle = extract_observation_bundle(history)
        up_d, up_f = _compute_equilibrium_curve(bundle, 'up')
        down_d, down_f = _compute_equilibrium_curve(bundle, 'down')

        step_fit = _best_step_function_r2(up_d, up_f)
        hyst = _compute_hysteresis(up_d, up_f, down_d, down_f)

        detector = P20QuorumSensingDetector(n_permutations=199, seed=seed)
        result = detector.detect(history, metadata)

        r2s.append(step_fit['r2'])
        crit_dens.append(step_fit['critical_density'])
        hyst_widths.append(hyst['hysteresis_width'])
        tiers.append(result.tier.value)
        confidences.append(result.confidence)

        print(f"  seed={seed}: R²={step_fit['r2']:.4f}, "
              f"d_c={step_fit['critical_density']:.3f}, "
              f"hyst={hyst['hysteresis_width']:.3f}, "
              f"tier={result.tier.value}")

    r2_arr = np.array(r2s)
    cd_arr = np.array(crit_dens)
    hw_arr = np.array(hyst_widths)

    summary = {
        'pattern': 'P20',
        'n_seeds': n_seeds,
        'seeds': seeds,
        'step_r2': {
            'mean': float(np.mean(r2_arr)),
            'std': float(np.std(r2_arr)),
            'cv': float(np.std(r2_arr) / np.mean(r2_arr)) if np.mean(r2_arr) > 0 else 0.0,
            'values': r2s,
        },
        'critical_density': {
            'mean': float(np.mean(cd_arr)),
            'std': float(np.std(cd_arr)),
            'cv': float(np.std(cd_arr) / np.mean(cd_arr)) if np.mean(cd_arr) > 0 else 0.0,
            'values': crit_dens,
        },
        'hysteresis_width': {
            'mean': float(np.mean(hw_arr)),
            'std': float(np.std(hw_arr)),
            'cv': float(np.std(hw_arr) / np.mean(hw_arr)) if np.mean(hw_arr) > 0 else 0.0,
            'values': hyst_widths,
        },
        'tiers': tiers,
        'confidences': confidences,
        'fraction_definitive': sum(1 for t in tiers if t == 'definitive') / n_seeds,
        'fraction_detected': sum(1 for t in tiers if t != 'screening' or confidences[i] > 0
                                  for i, _ in enumerate(tiers)) / n_seeds,
    }

    # Fix fraction_detected computation
    n_detected = sum(1 for i, t in enumerate(tiers) if t in ('confirmation', 'definitive'))
    summary['fraction_detected'] = n_detected / n_seeds

    output_dir = Path(__file__).resolve().parent.parent / 'outputs'
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / 'p20_multiseed.json'
    with open(output_path, 'w') as f:
        json.dump(summary, f, indent=2)

    print(f"\nP20 dim2 multi-seed summary:")
    print(f"  Step R²: {np.mean(r2_arr):.4f} ± {np.std(r2_arr):.4f} (CV={np.std(r2_arr)/np.mean(r2_arr)*100:.1f}%)")
    print(f"  Critical density: {np.mean(cd_arr):.3f} ± {np.std(cd_arr):.3f}")
    print(f"  Hysteresis width: {np.mean(hw_arr):.3f} ± {np.std(hw_arr):.3f}")
    print(f"  Fraction definitive: {summary['fraction_definitive']:.1%}")
    print(f"  Written to {output_path}")

    return summary


if __name__ == '__main__':
    main()
