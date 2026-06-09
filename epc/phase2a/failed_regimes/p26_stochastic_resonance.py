"""P26 (stochastic resonance) — Class C failed regimes.

Two canonical failure modes per the Sprint 75 brief:
  1. suprathreshold_signal: Signal amplitude exceeds the barrier/threshold so
     detection works at ALL noise levels (monotone performance, no inverted-U).
     The detector should reject — no interior peak in the performance curve.
  2. extreme_noise_only: Only very strong noise levels are swept (all far above
     optimal). Performance is flat or monotonically declining → no inverted-U.
"""

from __future__ import annotations

from typing import Any, Dict, List

import numpy as np


CONFIG: Dict[str, Any] = {
    "substrate_id": "P26_stochastic_resonance_class_c",
    "format": "noise_sweep",
    "status": "populated",
    "regimes": [
        {
            "label": "suprathreshold_signal",
            "description": (
                "Signal amplitude far exceeds threshold (A=3.0 >> threshold gap). "
                "Output tracks signal at all noise levels → monotone high performance, "
                "no inverted-U."
            ),
            "params": {
                "model": "threshold_unit",
                "threshold": 1.0,
                "signal_amplitude": 3.0,
                "signal_frequency": 0.05,
                "dt": 0.01,
                "n_steps": 5000,
                "noise_levels": [0.0, 0.1, 0.3, 0.5, 1.0, 2.0, 4.0, 8.0, 12.0, 20.0],
            },
            "seed": 42,
        },
        {
            "label": "extreme_noise_only",
            "description": (
                "All noise levels are far above optimal (D >= 10). "
                "Performance is flat near zero (random hopping dominates). "
                "No interior peak."
            ),
            "params": {
                "model": "bistable_double_well",
                "a": 4.0,
                "b": 1.0,
                "signal_amplitude": 1.0,
                "signal_frequency": 0.005,
                "dt": 0.01,
                "n_steps": 5000,
                "n_trials": 5,
                "noise_levels": [10.0, 15.0, 20.0, 30.0, 40.0, 60.0, 80.0, 100.0, 150.0, 200.0],
            },
            "seed": 42,
        },
    ],
}


def build_substrate(regime: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Build a P26 failed-regime substrate."""
    p = regime["params"]
    seed = regime["seed"]

    if p["model"] == "threshold_unit":
        from epc.models.stochastic_resonance import ThresholdUnit, ThresholdUnitParams
        model = ThresholdUnit(ThresholdUnitParams(
            threshold=p["threshold"],
            signal_amplitude=p["signal_amplitude"],
            signal_frequency=p["signal_frequency"],
            dt=p["dt"],
            n_steps=p["n_steps"],
            noise_levels=p["noise_levels"],
            seed=seed,
        ))
        return model.simulate()

    elif p["model"] == "bistable_double_well":
        from epc.models.stochastic_resonance import BistableDoubleWell, DoubleWellParams
        model = BistableDoubleWell(DoubleWellParams(
            a=p["a"],
            b=p["b"],
            signal_amplitude=p["signal_amplitude"],
            signal_frequency=p["signal_frequency"],
            dt=p["dt"],
            n_steps=p["n_steps"],
            n_trials=p["n_trials"],
            noise_levels=p["noise_levels"],
            seed=seed,
        ))
        return model.simulate()

    else:
        raise ValueError(f"Unknown model type: {p['model']}")
