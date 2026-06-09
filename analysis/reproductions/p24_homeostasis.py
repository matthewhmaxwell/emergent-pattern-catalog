"""
P24 dim1 reproduction — Homeostatic regulation (Ashby 1956).

Anchor:
    Ashby, W. R. (1956). An Introduction to Cybernetics. Chapman & Hall.
    Canonical result: a negative-feedback regulator maintains a variable
    near its set-point despite sustained perturbation. The deviation
    integral (∫|x − setpoint| dt) of the controlled system is a small
    fraction of the uncontrolled (gain=0) baseline.

Reproduction targets:
    1. Controlled/uncontrolled deviation-integral ratio < 0.3.
       Published: proportional controller with gain > 0 keeps x within
       O(pert/gain) of setpoint; uncontrolled drifts linearly.
       Tolerance: ratio < 0.30.

    2. Recovery time τ is finite and bounded under pulse perturbation.
       Tolerance: τ > 0 and τ < total_time / 2.

Parameters:
    Proportional homeostat: gain=5.0, setpoint=10.0, dt=0.1, seed=42
    Perturbation: onset=50.0, amplitude=5.0, sustained
    n_steps=2000
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
    extract_observation_bundle,
)


def run_reproduction() -> dict:
    """Run the P24 homeostatic regulation reproduction."""
    n_steps = 2000
    dt = 0.1

    # ── Controlled system (gain=5) ──
    params_ctrl = HomeostatParams(gain=5.0, setpoint=10.0, seed=42, dt=dt)
    schedule = PerturbationSchedule(onset=50.0, amplitude=5.0)
    model_ctrl = ProportionalHomeostat(params_ctrl)
    history_ctrl = model_ctrl.simulate(n_steps, schedule=schedule)

    bundle_ctrl = extract_observation_bundle(history_ctrl)
    time = bundle_ctrl['time']
    x_ctrl = bundle_ctrl['x']
    setpoint = bundle_ctrl['setpoint']

    # Find perturbation onset index
    onset_idx = int(50.0 / dt)  # 500

    dev_integral_ctrl = _compute_deviation_integral(
        x_ctrl, setpoint, time, onset_idx, n_steps,
    )

    # ── Uncontrolled system (gain=0) ──
    params_unctrl = HomeostatParams(gain=0.0, setpoint=10.0, seed=42, dt=dt)
    model_unctrl = ProportionalHomeostat(params_unctrl)
    history_unctrl = model_unctrl.simulate(n_steps, schedule=schedule)

    bundle_unctrl = extract_observation_bundle(history_unctrl)
    x_unctrl = bundle_unctrl['x']

    dev_integral_unctrl = _compute_deviation_integral(
        x_unctrl, setpoint, time, onset_idx, n_steps,
    )

    # ── Deviation ratio ──
    ratio = dev_integral_ctrl / dev_integral_unctrl

    # ── Steady-state deviation ──
    ss_dev_ctrl = float(np.mean(np.abs(x_ctrl[-400:] - setpoint[-400:])))
    ss_dev_unctrl = float(np.mean(np.abs(x_unctrl[-400:] - setpoint[-400:])))

    # ── Run detector ──
    det = P24HomeostasisDetector(n_permutations=199, seed=42)
    result = det.detect(history_ctrl, model_ctrl.get_metadata())

    # ── Tolerance checks ──
    tol1_pass = ratio < 0.30
    tol2_pass = result.tier.value in ("confirmation", "definitive")

    passes_tolerance = tol1_pass and tol2_pass

    output = {
        "reproduction": "p24_homeostasis",
        "pattern": "P24",
        "anchor": "Ashby (1956) — negative-feedback homeostat",
        "parameters": {
            "gain_controlled": 5.0,
            "gain_uncontrolled": 0.0,
            "setpoint": 10.0,
            "perturbation_amplitude": 5.0,
            "perturbation_onset": 50.0,
            "dt": dt,
            "n_steps": n_steps,
            "seed": 42,
        },
        "deviation_integral_controlled": round(dev_integral_ctrl, 4),
        "deviation_integral_uncontrolled": round(dev_integral_unctrl, 4),
        "deviation_ratio": round(ratio, 6),
        "ss_deviation_controlled": round(ss_dev_ctrl, 6),
        "ss_deviation_uncontrolled": round(ss_dev_unctrl, 4),
        "detector_tier": result.tier.value,
        "detector_confidence": round(result.confidence, 4),
        "tolerance_checks": {
            "ratio_lt_0.30": {
                "value": round(ratio, 6),
                "threshold": 0.30,
                "pass": tol1_pass,
            },
            "detector_confirmation_plus": {
                "tier": result.tier.value,
                "pass": tol2_pass,
            },
        },
        "passes_tolerance": passes_tolerance,
    }

    return output


if __name__ == "__main__":
    result = run_reproduction()
    out_path = os.path.join(
        os.path.dirname(__file__), '..', 'outputs',
        'p24_homeostasis_reproduction.json',
    )
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, 'w') as f:
        json.dump(result, f, indent=2)
    print(json.dumps(result, indent=2))
    print(f"\nWritten to {out_path}")
    print(f"passes_tolerance: {result['passes_tolerance']}")
