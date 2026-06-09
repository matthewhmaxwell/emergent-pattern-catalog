"""
P26 dim2 multi-seed campaign — Stochastic resonance.

Run ≥20 seeds of the bistable double-well and report peak-noise
location and peak coherent response mean ± std ± CV.

Parameters:
    Bistable double-well: a=4.0, b=1.0, barrier=4.0
    Signal: amplitude=1.0, frequency=0.005
    n_trials=20 per noise level per seed, n_steps=20000
    n_seeds=20
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


def run_multiseed() -> dict:
    """Run 20-seed multi-seed campaign."""
    n_seeds = 20
    n_trials = 20
    n_steps = 20000

    peak_noises = []
    peak_perfs = []
    gains = []
    tiers = []

    for seed in range(n_seeds):
        params = DoubleWellParams(seed=seed, n_trials=n_trials, n_steps=n_steps)
        model = BistableDoubleWell(params)
        history = model.simulate()

        bundle = extract_observation_bundle(history)
        noise_levels, performance = _compute_performance_curve(bundle)

        peak_idx = int(np.argmax(performance))
        peak_noises.append(float(noise_levels[peak_idx]))
        peak_perfs.append(float(performance[peak_idx]))
        gains.append(float(performance[peak_idx] - performance[0]))

        det = P26StochasticResonanceDetector(n_permutations=199, seed=seed)
        result = det.detect(history, model.get_metadata())
        tiers.append(result.tier.value)

    pn_arr = np.array(peak_noises)
    pp_arr = np.array(peak_perfs)
    g_arr = np.array(gains)

    output = {
        "reproduction": "p26_multiseed",
        "pattern": "P26",
        "n_seeds": n_seeds,
        "n_trials": n_trials,
        "n_steps": n_steps,
        "parameters": {
            "a": 4.0,
            "b": 1.0,
            "signal_amplitude": 1.0,
            "signal_frequency": 0.005,
            "dt": 0.01,
        },
        "peak_noise_values": [round(v, 2) for v in peak_noises],
        "peak_noise_mean": round(float(np.mean(pn_arr)), 2),
        "peak_noise_std": round(float(np.std(pn_arr)), 4),
        "peak_performance_values": [round(v, 6) for v in peak_perfs],
        "peak_performance_mean": round(float(np.mean(pp_arr)), 6),
        "peak_performance_std": round(float(np.std(pp_arr)), 6),
        "peak_performance_cv": round(
            float(np.std(pp_arr) / np.mean(pp_arr) * 100), 1
        ) if np.mean(pp_arr) > 0 else 0.0,
        "gain_values": [round(v, 6) for v in gains],
        "gain_mean": round(float(np.mean(g_arr)), 6),
        "gain_std": round(float(np.std(g_arr)), 6),
        "gain_cv": round(
            float(np.std(g_arr) / np.mean(g_arr) * 100), 1
        ) if np.mean(g_arr) > 0 else 0.0,
        "tiers": tiers,
        "fraction_definitive": round(
            sum(1 for t in tiers if t == "definitive") / n_seeds, 2
        ),
    }

    return output


if __name__ == "__main__":
    result = run_multiseed()
    out_path = os.path.join(
        os.path.dirname(__file__), '..', 'outputs',
        'p26_multiseed.json',
    )
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, 'w') as f:
        json.dump(result, f, indent=2)
    print(json.dumps(result, indent=2))
    print(f"\nWritten to {out_path}")
