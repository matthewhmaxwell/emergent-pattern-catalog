"""P19 (emergent leadership) failed regimes — Class C.

Two failure modes per Sprint 70 brief:

1. **rho_zero**: informed_fraction=0 — no informed agents at all. The flock
   has no preferred direction and cannot exhibit minority-guided leadership.
   Group accuracy is chance-level (~0). True negative for P19.

2. **bias_zero**: informed agents exist (ρ=0.1) but bias_weight=0 — informed
   agents have no preference and behave identically to naive agents. No
   directional pull, no asymmetric influence. True negative for P19.

Both regimes violate the defining P19 signature: a small informed minority
steers the group toward a preferred direction via asymmetric influence
(Couzin et al. 2005).
"""

from __future__ import annotations

from typing import Any, Dict, List


CONFIG: Dict[str, Any] = {
    "substrate_id": "P19_informed_minority_failed",
    "format": "particles",
    "description": (
        "P19 failed regimes: (1) rho_zero — no informed agents, group accuracy "
        "is chance-level; (2) bias_zero — informed agents have no preference, "
        "no directional pull. Both lack the P19 minority-guidance signature."
    ),
    "regimes": [
        # Class C-1: ρ=0 (no informed agents) — 5 seeds
        {
            "label": "rho_zero_seed=400",
            "params": {
                "n_particles": 200,
                "box_size": 10.0,
                "speed": 0.03,
                "noise": 0.1,
                "interaction_radius": 1.0,
                "informed_fraction": 0.0,
                "bias_weight": 0.3,
                "preferred_direction": 0.0,
                "n_steps": 500,
            },
            "seed": 400,
        },
        {
            "label": "rho_zero_seed=401",
            "params": {
                "n_particles": 200,
                "box_size": 10.0,
                "speed": 0.03,
                "noise": 0.1,
                "interaction_radius": 1.0,
                "informed_fraction": 0.0,
                "bias_weight": 0.3,
                "preferred_direction": 0.0,
                "n_steps": 500,
            },
            "seed": 401,
        },
        {
            "label": "rho_zero_seed=402",
            "params": {
                "n_particles": 200,
                "box_size": 10.0,
                "speed": 0.03,
                "noise": 0.1,
                "interaction_radius": 1.0,
                "informed_fraction": 0.0,
                "bias_weight": 0.3,
                "preferred_direction": 0.0,
                "n_steps": 500,
            },
            "seed": 402,
        },
        {
            "label": "rho_zero_seed=403",
            "params": {
                "n_particles": 200,
                "box_size": 10.0,
                "speed": 0.03,
                "noise": 0.1,
                "interaction_radius": 1.0,
                "informed_fraction": 0.0,
                "bias_weight": 0.3,
                "preferred_direction": 0.0,
                "n_steps": 500,
            },
            "seed": 403,
        },
        {
            "label": "rho_zero_seed=404",
            "params": {
                "n_particles": 200,
                "box_size": 10.0,
                "speed": 0.03,
                "noise": 0.1,
                "interaction_radius": 1.0,
                "informed_fraction": 0.0,
                "bias_weight": 0.3,
                "preferred_direction": 0.0,
                "n_steps": 500,
            },
            "seed": 404,
        },
        # Class C-2: bias_weight=0 (informed agents have no preference) — 5 seeds
        {
            "label": "bias_zero_seed=410",
            "params": {
                "n_particles": 200,
                "box_size": 10.0,
                "speed": 0.03,
                "noise": 0.1,
                "interaction_radius": 1.0,
                "informed_fraction": 0.1,
                "bias_weight": 0.0,
                "preferred_direction": 0.0,
                "n_steps": 500,
            },
            "seed": 410,
        },
        {
            "label": "bias_zero_seed=411",
            "params": {
                "n_particles": 200,
                "box_size": 10.0,
                "speed": 0.03,
                "noise": 0.1,
                "interaction_radius": 1.0,
                "informed_fraction": 0.1,
                "bias_weight": 0.0,
                "preferred_direction": 0.0,
                "n_steps": 500,
            },
            "seed": 411,
        },
        {
            "label": "bias_zero_seed=412",
            "params": {
                "n_particles": 200,
                "box_size": 10.0,
                "speed": 0.03,
                "noise": 0.1,
                "interaction_radius": 1.0,
                "informed_fraction": 0.1,
                "bias_weight": 0.0,
                "preferred_direction": 0.0,
                "n_steps": 500,
            },
            "seed": 412,
        },
        {
            "label": "bias_zero_seed=413",
            "params": {
                "n_particles": 200,
                "box_size": 10.0,
                "speed": 0.03,
                "noise": 0.1,
                "interaction_radius": 1.0,
                "informed_fraction": 0.1,
                "bias_weight": 0.0,
                "preferred_direction": 0.0,
                "n_steps": 500,
            },
            "seed": 413,
        },
        {
            "label": "bias_zero_seed=414",
            "params": {
                "n_particles": 200,
                "box_size": 10.0,
                "speed": 0.03,
                "noise": 0.1,
                "interaction_radius": 1.0,
                "informed_fraction": 0.1,
                "bias_weight": 0.0,
                "preferred_direction": 0.0,
                "n_steps": 500,
            },
            "seed": 414,
        },
    ],
}


def build_substrate(regime: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Run a P19 failed regime and return particle-history."""
    from epc.models.informed_minority import InformedMinorityModel

    p = regime["params"]
    model = InformedMinorityModel(
        n_particles=p["n_particles"],
        box_size=p["box_size"],
        speed=p["speed"],
        noise=p["noise"],
        interaction_radius=p["interaction_radius"],
        informed_fraction=p["informed_fraction"],
        bias_weight=p["bias_weight"],
        preferred_direction=p["preferred_direction"],
        seed=regime["seed"],
    )
    return model.run(n_steps=p["n_steps"])
