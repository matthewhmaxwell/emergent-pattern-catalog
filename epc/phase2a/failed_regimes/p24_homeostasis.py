"""P24 (homeostatic regulation) — Class C failed regimes.

Two canonical failure modes per the Sprint 73 brief:
  1. gain_zero: feedback gain = 0 → system drifts under perturbation (no
     regulation). The detector should reject at screening (growth_ratio >> 2.0).
  2. no_perturbation: perturbation amplitude = 0 → variable stays trivially at
     setpoint but regulation is NOT demonstrated. The detector should reject
     at the "No perturbation detected" prerequisite.
"""

from __future__ import annotations

from typing import Any, Dict, List


CONFIG: Dict[str, Any] = {
    "substrate_id": "P24_homeostasis_class_c",
    "format": "scalar_timeseries",
    "status": "populated",
    "regimes": [
        {
            "label": "gain_zero_drift",
            "description": "Feedback gain = 0; system drifts under sustained perturbation.",
            "params": {
                "gain": 0.0,
                "setpoint": 10.0,
                "dt": 0.1,
                "noise_std": 0.5,
                "n_steps": 1000,
                "pert_onset": 50.0,
                "pert_amplitude": 5.0,
            },
            "seed": 42,
        },
        {
            "label": "no_perturbation",
            "description": "No external perturbation; variable stays at setpoint trivially.",
            "params": {
                "gain": 5.0,
                "setpoint": 10.0,
                "dt": 0.1,
                "noise_std": 0.5,
                "n_steps": 1000,
                "pert_onset": None,
                "pert_amplitude": 0.0,
            },
            "seed": 42,
        },
    ],
}


def build_substrate(regime: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Build a P24 failed-regime substrate."""
    from epc.models.homeostasis import (
        ProportionalHomeostat,
        HomeostatParams,
        PerturbationSchedule,
    )

    p = regime["params"]
    model = ProportionalHomeostat(HomeostatParams(
        setpoint=p["setpoint"],
        gain=p["gain"],
        dt=p["dt"],
        noise_std=p["noise_std"],
        seed=regime["seed"],
    ))

    if p["pert_onset"] is None or p["pert_amplitude"] == 0.0:
        schedule = PerturbationSchedule(onset=float("inf"), amplitude=0.0)
    else:
        schedule = PerturbationSchedule(
            onset=p["pert_onset"],
            amplitude=p["pert_amplitude"],
        )

    return model.simulate(n_steps=p["n_steps"], schedule=schedule)
