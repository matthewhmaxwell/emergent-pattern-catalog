"""P17 (collective sensing) failed regimes — Class C.

Two failure modes that produce motion toward the field peak but do NOT exhibit
the P17 signature (group CI scaling with N):

1. **social_off**: social_strength=0 — agents sense independently. Each agent's
   noisy estimate is unaveraged. CI is low and does NOT improve with N because
   there is no group-level information pooling. True negative for P17.

2. **field_too_strong**: field_amplitude=10.0, sensing_noise=0.1 — the field is
   so strong relative to noise that even a single agent can reliably climb the
   gradient. CI is high at ALL N including N=1, so there is no group advantage
   (no CI-vs-N scaling). True negative for P17.

Both regimes violate the defining P17 signature: CI(N) >> CI(1) via collective
information pooling. The detector should reject because:
- social_off: slope(CI vs log(N)) ≈ 0, CI(N_max) ≈ CI(1) (both low)
- field_too_strong: slope ≈ 0, CI(1) already high (no *gain* from group)
"""

from __future__ import annotations

from typing import Any, Dict, List

import numpy as np


CONFIG: Dict[str, Any] = {
    "substrate_id": "P17_collective_sensing_failed",
    "format": "particles",
    "description": (
        "P17 failed regimes: (1) social_off — no social coupling, agents sense "
        "independently, no group CI scaling; (2) field_too_strong — individuals "
        "succeed alone, no group advantage. Both lack the N-scaling signature."
    ),
    "regimes": [
        # Class C-1: social coupling off (5 regimes at varying noise)
        {
            "label": "social_off_noise=0.8",
            "params": {
                "n_agents": 50,
                "box_size": 20.0,
                "v_max": 0.4,
                "turn_noise": 0.3,
                "sensing_noise": 0.8,
                "alpha": 0.95,
                "social_strength": 0.0,  # KEY: no social coupling
                "field_sigma": 5.0,
                "field_amplitude": 1.0,
                "n_steps": 1000,
            },
            "seed": 300,
        },
        {
            "label": "social_off_noise=1.0",
            "params": {
                "n_agents": 50,
                "box_size": 20.0,
                "v_max": 0.4,
                "turn_noise": 0.3,
                "sensing_noise": 1.0,
                "alpha": 0.95,
                "social_strength": 0.0,
                "field_sigma": 5.0,
                "field_amplitude": 1.0,
                "n_steps": 1000,
            },
            "seed": 301,
        },
        {
            "label": "social_off_noise=1.2",
            "params": {
                "n_agents": 50,
                "box_size": 20.0,
                "v_max": 0.4,
                "turn_noise": 0.3,
                "sensing_noise": 1.2,
                "alpha": 0.95,
                "social_strength": 0.0,
                "field_sigma": 5.0,
                "field_amplitude": 1.0,
                "n_steps": 1000,
            },
            "seed": 302,
        },
        {
            "label": "social_off_noise=0.5",
            "params": {
                "n_agents": 50,
                "box_size": 20.0,
                "v_max": 0.4,
                "turn_noise": 0.3,
                "sensing_noise": 0.5,
                "alpha": 0.95,
                "social_strength": 0.0,
                "field_sigma": 5.0,
                "field_amplitude": 1.0,
                "n_steps": 1000,
            },
            "seed": 303,
        },
        {
            "label": "social_off_noise=1.5",
            "params": {
                "n_agents": 50,
                "box_size": 20.0,
                "v_max": 0.4,
                "turn_noise": 0.3,
                "sensing_noise": 1.5,
                "alpha": 0.95,
                "social_strength": 0.0,
                "field_sigma": 5.0,
                "field_amplitude": 1.0,
                "n_steps": 1000,
            },
            "seed": 304,
        },
        # Class C-2: field too strong (individuals succeed alone)
        {
            "label": "field_strong_amp=10_noise=0.1",
            "params": {
                "n_agents": 50,
                "box_size": 20.0,
                "v_max": 0.4,
                "turn_noise": 0.3,
                "sensing_noise": 0.1,  # Very low noise
                "alpha": 0.95,
                "social_strength": 0.2,
                "field_sigma": 5.0,
                "field_amplitude": 10.0,  # KEY: very strong field
                "n_steps": 1000,
            },
            "seed": 310,
        },
        {
            "label": "field_strong_amp=5_noise=0.1",
            "params": {
                "n_agents": 50,
                "box_size": 20.0,
                "v_max": 0.4,
                "turn_noise": 0.3,
                "sensing_noise": 0.1,
                "alpha": 0.95,
                "social_strength": 0.2,
                "field_sigma": 5.0,
                "field_amplitude": 5.0,
                "n_steps": 1000,
            },
            "seed": 311,
        },
        {
            "label": "field_strong_amp=10_noise=0.05",
            "params": {
                "n_agents": 50,
                "box_size": 20.0,
                "v_max": 0.4,
                "turn_noise": 0.1,  # Also reduced turn noise
                "sensing_noise": 0.05,
                "alpha": 0.95,
                "social_strength": 0.2,
                "field_sigma": 5.0,
                "field_amplitude": 10.0,
                "n_steps": 1000,
            },
            "seed": 312,
        },
        {
            "label": "field_strong_amp=8_noise=0.15",
            "params": {
                "n_agents": 50,
                "box_size": 20.0,
                "v_max": 0.4,
                "turn_noise": 0.2,
                "sensing_noise": 0.15,
                "alpha": 0.95,
                "social_strength": 0.2,
                "field_sigma": 5.0,
                "field_amplitude": 8.0,
                "n_steps": 1000,
            },
            "seed": 313,
        },
        {
            "label": "field_strong_amp=15_noise=0.05",
            "params": {
                "n_agents": 50,
                "box_size": 20.0,
                "v_max": 0.4,
                "turn_noise": 0.1,
                "sensing_noise": 0.05,
                "alpha": 0.95,
                "social_strength": 0.2,
                "field_sigma": 5.0,
                "field_amplitude": 15.0,
                "n_steps": 1000,
            },
            "seed": 314,
        },
    ],
}


def build_substrate(regime: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Run a P17 failed regime and return particle-history."""
    from epc.models.collective_sensing import CollectiveSensingModel

    p = regime["params"]
    model = CollectiveSensingModel(
        n_agents=p["n_agents"],
        box_size=p["box_size"],
        v_max=p["v_max"],
        turn_noise=p["turn_noise"],
        sensing_noise=p["sensing_noise"],
        alpha=p["alpha"],
        social_strength=p["social_strength"],
        field_sigma=p["field_sigma"],
        field_amplitude=p["field_amplitude"],
        dt=1.0,
        init_mode="offset",
        seed=regime["seed"],
    )
    return model.run(n_steps=p["n_steps"], record_interval=5)
