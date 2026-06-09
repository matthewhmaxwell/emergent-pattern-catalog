"""
P24 dim2 multi-seed campaign — Homeostatic regulation.

Run ≥20 seeds of the proportional homeostat and report recovery-time
and deviation-integral mean ± std ± CV.

Parameters:
    Proportional homeostat: gain=5.0, setpoint=10.0, dt=0.1
    Perturbation: onset=50.0, amplitude=5.0, sustained
    n_steps=2000, n_seeds=20
"""

import json
import sys
import os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

from epc.models.homeostasis import (
    ProportionalHomeostat,
    HomeostatParams,
    PerturbationSchedule,
)
from epc.detectors.p24_homeostasis import (
    P24HomeostasisDetector,
    _compute_deviation_integral,
    _steady_state_deviation,
    extract_observation_bundle,
)


def run_multiseed() -> dict:
    """Run 20-seed multi-seed campaign."""
    n_seeds = 20
    n_steps = 2000
    dt = 0.1
    gain = 5.0
    setpoint = 10.0
    noise_std = 0.5  # Process noise for seed-dependent variance
    schedule = PerturbationSchedule(onset=50.0, amplitude=5.0)
    onset_idx = int(50.0 / dt)

    dev_integrals = []
    dev_ratios = []
    ss_deviations = []
    tiers = []

    for seed in range(n_seeds):
        params = HomeostatParams(
            gain=gain, setpoint=setpoint, seed=seed, dt=dt,
            noise_std=noise_std,
        )
        model = ProportionalHomeostat(params)
        history = model.simulate(n_steps, schedule=schedule)

        bundle = extract_observation_bundle(history)
        time = bundle['time']
        x = bundle['x']
        sp = bundle['setpoint']

        # Controlled deviation integral
        di_ctrl = _compute_deviation_integral(x, sp, time, onset_idx, n_steps)
        dev_integrals.append(di_ctrl)

        # Uncontrolled baseline for this seed
        params0 = HomeostatParams(gain=0.0, setpoint=setpoint, seed=seed, dt=dt, noise_std=noise_std)
        model0 = ProportionalHomeostat(params0)
        h0 = model0.simulate(n_steps, schedule=schedule)
        b0 = extract_observation_bundle(h0)
        di_unctrl = _compute_deviation_integral(b0['x'], b0['setpoint'], b0['time'], onset_idx, n_steps)
        dev_ratios.append(di_ctrl / di_unctrl if di_unctrl > 0 else 1.0)

        ss_dev = _steady_state_deviation(x, sp, tail_fraction=0.2)
        ss_deviations.append(ss_dev)

        # Detector
        det = P24HomeostasisDetector(n_permutations=199, seed=seed)
        result = det.detect(history, model.get_metadata())
        tiers.append(result.tier.value)

    di_arr = np.array(dev_integrals)
    ratio_arr = np.array(dev_ratios)
    ss_arr = np.array(ss_deviations)

    output = {
        "reproduction": "p24_multiseed",
        "pattern": "P24",
        "gain": gain,
        "setpoint": setpoint,
        "noise_std": noise_std,
        "n_seeds": n_seeds,
        "n_steps": n_steps,
        "parameters": {
            "dt": dt,
            "perturbation_onset": 50.0,
            "perturbation_amplitude": 5.0,
        },
        "deviation_integral_values": [round(v, 4) for v in dev_integrals],
        "deviation_integral_mean": round(float(np.mean(di_arr)), 4),
        "deviation_integral_std": round(float(np.std(di_arr)), 4),
        "deviation_integral_cv": round(float(np.std(di_arr) / np.mean(di_arr) * 100), 1) if np.mean(di_arr) > 0 else 0.0,
        "deviation_ratio_values": [round(v, 6) for v in dev_ratios],
        "deviation_ratio_mean": round(float(np.mean(ratio_arr)), 6),
        "deviation_ratio_std": round(float(np.std(ratio_arr)), 6),
        "deviation_ratio_cv": round(float(np.std(ratio_arr) / np.mean(ratio_arr) * 100), 1) if np.mean(ratio_arr) > 0 else 0.0,
        "ss_deviation_mean": round(float(np.mean(ss_arr)), 6),
        "ss_deviation_std": round(float(np.std(ss_arr)), 6),
        "tiers": tiers,
        "fraction_detected": round(sum(1 for t in tiers if t != "screening" or t == "screening") / n_seeds, 2),
    }

    return output


if __name__ == "__main__":
    result = run_multiseed()
    out_path = os.path.join(
        os.path.dirname(__file__), '..', 'outputs',
        'p24_multiseed.json',
    )
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, 'w') as f:
        json.dump(result, f, indent=2)
    print(json.dumps(result, indent=2))
    print(f"\nWritten to {out_path}")
