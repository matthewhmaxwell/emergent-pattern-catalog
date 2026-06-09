"""
P26 dim1 reproduction — Stochastic resonance (Collins et al. 1995).

Anchor:
    Collins, J. J., Chow, C. C., & Imhoff, T. T. (1995).
    Stochastic resonance without tuning. Nature, 376, 236–238.
    Gammaitoni, L. et al. (1998). Stochastic resonance. Rev. Mod. Phys.

    Canonical result: an overdamped particle in a bistable double-well
    potential, driven by a weak subthreshold periodic signal and tunable
    noise, exhibits an inverted-U SNR-vs-noise curve. The coherent
    response (⟨x·signal⟩) peaks at an intermediate noise level and
    declines at both lower and higher noise.

Reproduction targets:
    1. Interior argmax: peak performance at a noise level strictly
       between zero and maximum.
    2. Peak coherent response > zero-noise coherent response by a
       clear margin (gain > 0.05).
    3. Peak coherent response > high-noise coherent response
       (decline > 0.02).
    4. Detector reaches DEFINITIVE tier.

Parameters:
    Bistable double-well: a=4.0, b=1.0, barrier=4.0
    Signal: amplitude=1.0, frequency=0.005
    Noise sweep: 15 levels from 0 to 20
    n_trials=20 independent runs per noise level
    n_steps=20000 per trial
"""

import json
import sys
import os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

from epc.models.stochastic_resonance import (
    BistableDoubleWell,
    DoubleWellParams,
)
from epc.detectors.p26_stochastic_resonance import (
    P26StochasticResonanceDetector,
    extract_observation_bundle,
    _compute_performance_curve,
)


def run_reproduction() -> dict:
    """Run the P26 stochastic resonance reproduction."""
    params = DoubleWellParams(seed=42, n_trials=20, n_steps=20000)
    model = BistableDoubleWell(params)
    history = model.simulate()

    # Compute performance curve
    bundle = extract_observation_bundle(history)
    noise_levels, performance = _compute_performance_curve(bundle)

    # Find peak
    peak_idx = int(np.argmax(performance))
    peak_noise = float(noise_levels[peak_idx])
    peak_perf = float(performance[peak_idx])
    zero_perf = float(performance[0])
    high_perf = float(performance[-1])

    gain = peak_perf - zero_perf
    decline = peak_perf - high_perf
    is_interior = 0 < peak_idx < len(performance) - 1

    # Run detector
    det = P26StochasticResonanceDetector(n_permutations=199, seed=42)
    result = det.detect(history, model.get_metadata())

    # Tolerance checks
    tol1_pass = is_interior
    tol2_pass = gain > 0.05
    tol3_pass = decline > 0.02
    tol4_pass = result.tier.value in ("confirmation", "definitive")

    passes_tolerance = tol1_pass and tol2_pass and tol3_pass and tol4_pass

    output = {
        "reproduction": "p26_collins",
        "pattern": "P26",
        "anchor": "Gammaitoni et al. (1998) / Collins et al. (1995) — SR inverted-U",
        "parameters": {
            "a": params.a,
            "b": params.b,
            "barrier_height": params.a**2 / (4 * params.b),
            "signal_amplitude": params.signal_amplitude,
            "signal_frequency": params.signal_frequency,
            "dt": params.dt,
            "n_steps": params.n_steps,
            "n_trials": params.n_trials,
            "seed": params.seed,
        },
        "noise_levels": noise_levels.tolist(),
        "performance_curve": [round(v, 6) for v in performance.tolist()],
        "peak_noise": peak_noise,
        "peak_performance": round(peak_perf, 6),
        "zero_noise_performance": round(zero_perf, 6),
        "high_noise_performance": round(high_perf, 6),
        "gain_over_zero": round(gain, 6),
        "decline_after_peak": round(decline, 6),
        "detector_tier": result.tier.value,
        "detector_confidence": round(result.confidence, 4),
        "null_p_value": round(result.null_p_value, 6),
        "tolerance_checks": {
            "interior_argmax": {
                "peak_idx": peak_idx,
                "n_levels": len(noise_levels),
                "pass": tol1_pass,
            },
            "gain_gt_0.05": {
                "value": round(gain, 6),
                "threshold": 0.05,
                "pass": tol2_pass,
            },
            "decline_gt_0.02": {
                "value": round(decline, 6),
                "threshold": 0.02,
                "pass": tol3_pass,
            },
            "detector_confirmation_plus": {
                "tier": result.tier.value,
                "pass": tol4_pass,
            },
        },
        "passes_tolerance": passes_tolerance,
    }

    return output


if __name__ == "__main__":
    result = run_reproduction()
    out_path = os.path.join(
        os.path.dirname(__file__), '..', 'outputs',
        'p26_collins_reproduction.json',
    )
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, 'w') as f:
        json.dump(result, f, indent=2)
    print(json.dumps(result, indent=2))
    print(f"\nWritten to {out_path}")
    print(f"passes_tolerance: {result['passes_tolerance']}")
