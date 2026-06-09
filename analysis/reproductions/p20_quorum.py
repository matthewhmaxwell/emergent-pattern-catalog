"""P20 dim1 reproduction: quorum sensing / threshold-activated response.

Anchor: bacterial quorum sensing (Waters & Bassler 2005; Nealson et al. 1970).
Canonical signature: a SHARP collective-activation transition at a critical
density (transition width small relative to the density range) AND a
hysteresis loop (activation density meaningfully > deactivation density).

Tolerance checks:
  1. Step-function R² > 0.9 (sharp, switch-like transition)
  2. Hysteresis width > 0.1 (activation density > deactivation density)
  3. Detector tier ∈ {confirmation, definitive}
"""

from __future__ import annotations

import json
import os
from pathlib import Path

from epc.models.quorum_sensing import AutoinducerQuorum, AutoinducerParams
from epc.detectors.p20_quorum_sensing import (
    P20QuorumSensingDetector,
    extract_observation_bundle,
    _compute_equilibrium_curve,
    _best_step_function_r2,
    _compute_hysteresis,
)


def main() -> dict:
    """Run P20 dim1 reproduction and return results."""

    # ── Canonical regime ──
    params = AutoinducerParams(seed=42)
    model = AutoinducerQuorum(params)
    history = model.simulate()
    metadata = model.get_metadata()

    # ── Extract curves ──
    bundle = extract_observation_bundle(history)
    up_dens, up_frac = _compute_equilibrium_curve(bundle, 'up')
    down_dens, down_frac = _compute_equilibrium_curve(bundle, 'down')

    # ── Step-function fit ──
    step_fit = _best_step_function_r2(up_dens, up_frac)

    # ── Hysteresis ──
    hyst = _compute_hysteresis(up_dens, up_frac, down_dens, down_frac)

    # ── Run detector ──
    detector = P20QuorumSensingDetector(n_permutations=199, seed=42)
    result = detector.detect(history, metadata)

    # ── Tolerance checks ──
    r2_pass = step_fit['r2'] > 0.9
    hyst_pass = hyst['hysteresis_width'] > 0.1
    tier_pass = result.tier.value in ('confirmation', 'definitive')
    all_pass = r2_pass and hyst_pass and tier_pass

    results = {
        'pattern': 'P20',
        'anchor': 'Waters & Bassler (2005); Nealson et al. (1970)',
        'model': 'AutoinducerQuorum',
        'params': {
            'base_production': params.base_production,
            'enhanced_production': params.enhanced_production,
            'degradation_rate': params.degradation_rate,
            'activation_threshold': params.activation_threshold,
            'deactivation_threshold': params.deactivation_threshold,
            'n_density_levels': params.n_density_levels,
            'density_range': [params.density_min, params.density_max],
            'seed': params.seed,
        },
        'step_r2': step_fit['r2'],
        'critical_density': step_fit['critical_density'],
        'f_low': step_fit['f_low'],
        'f_high': step_fit['f_high'],
        'activation_density': hyst['activation_density'],
        'deactivation_density': hyst['deactivation_density'],
        'hysteresis_width': hyst['hysteresis_width'],
        'has_hysteresis': hyst['has_hysteresis'],
        'detector_tier': result.tier.value,
        'detector_confidence': result.confidence,
        'null_p_value': result.null_p_value,
        'cohens_d': result.effect_size.get('cohens_d', 0.0),
        'tolerances': {
            'r2_gt_0.9': r2_pass,
            'hysteresis_gt_0.1': hyst_pass,
            'tier_confirmed_or_definitive': tier_pass,
        },
        'passes_tolerance': all_pass,
    }

    # ── Emit JSON ──
    output_dir = Path(__file__).resolve().parent.parent / 'outputs'
    output_dir.mkdir(parents=True, exist_ok=True)
    output_path = output_dir / 'p20_quorum_reproduction.json'
    with open(output_path, 'w') as f:
        json.dump(results, f, indent=2)

    print(f"P20 dim1 reproduction:")
    print(f"  Step-function R² = {step_fit['r2']:.4f} (> 0.9: {'PASS' if r2_pass else 'FAIL'})")
    print(f"  Critical density = {step_fit['critical_density']:.3f}")
    print(f"  Hysteresis width = {hyst['hysteresis_width']:.3f} (> 0.1: {'PASS' if hyst_pass else 'FAIL'})")
    print(f"  Activation density = {hyst['activation_density']:.3f}")
    print(f"  Deactivation density = {hyst['deactivation_density']:.3f}")
    print(f"  Detector tier = {result.tier.value} (PASS: {'PASS' if tier_pass else 'FAIL'})")
    print(f"  Overall: {'PASS' if all_pass else 'FAIL'}")
    print(f"  Written to {output_path}")

    return results


if __name__ == '__main__':
    main()
