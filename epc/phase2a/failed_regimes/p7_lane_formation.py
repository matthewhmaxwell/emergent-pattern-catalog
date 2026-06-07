"""P7 (Lane Formation) failed regimes: 10 weak-repulsion counterflow regimes.

Lane formation requires sufficient social-force repulsion (A >> 1) for agents
to laterally segregate by direction. At weak repulsion (A ≤ 2.0), agents
interpenetrate freely and no directional lanes form: φ_lane ≈ 0.

Canonical lane-forming regime: A = 5.0, B = 0.3 (Helbing & Molnár 1995).

Also includes single-population regimes (no counterflow) at indices 8-9:
when all agents share the same desired direction, there is no opposing flow
to segregate from, so lane order parameter is meaningless.
"""

from __future__ import annotations

from typing import Any, Dict, List

import numpy as np


REPULSION_VALUES = list(np.linspace(0.1, 0.8, 8))

CONFIG: Dict[str, Any] = {
    "substrate_id": "P7_lane_formation_weak_repulsion",
    "format": "particles",
    "description": (
        "10 counterflow regimes that suppress lane formation. "
        "Regimes 0-7: repulsion_amplitude A ∈ linspace(0.1, 0.8, 8) with "
        "canonical corridor 20×4, N=200. A ≤ 0.8 produces insufficient "
        "repulsion for lane segregation; φ_lane ≈ 0. "
        "Regimes 8-9: single-population free-flow (no counterflow). "
        "Canonical lane regime uses A=5.0."
    ),
    "regimes": [
        *(
            {
                "label": f"A={a:.3f}",
                "params": {
                    "n_agents": 200,
                    "corridor_width": 20.0,
                    "corridor_height": 4.0,
                    "desired_speed": 1.0,
                    "repulsion_amplitude": float(a),
                    "repulsion_range": 0.3,
                    "tau": 0.5,
                    "dt": 0.05,
                    "n_steps": 400,
                    "single_population": False,
                },
                "seed": 300 + i,
            }
            for i, a in enumerate(REPULSION_VALUES)
        ),
        # Single-population: no counterflow → no lane formation possible.
        {
            "label": "single_pop_low_density",
            "params": {
                "n_agents": 100,
                "corridor_width": 20.0,
                "corridor_height": 4.0,
                "desired_speed": 1.0,
                "repulsion_amplitude": 5.0,
                "repulsion_range": 0.3,
                "tau": 0.5,
                "dt": 0.05,
                "n_steps": 400,
                "single_population": True,
            },
            "seed": 310,
        },
        {
            "label": "single_pop_high_density",
            "params": {
                "n_agents": 200,
                "corridor_width": 20.0,
                "corridor_height": 4.0,
                "desired_speed": 1.0,
                "repulsion_amplitude": 5.0,
                "repulsion_range": 0.3,
                "tau": 0.5,
                "dt": 0.05,
                "n_steps": 400,
                "single_population": True,
            },
            "seed": 311,
        },
    ],
}


def build_substrate(regime: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Run a weak-repulsion or single-population counterflow regime."""
    from epc.models.lane_formation import LaneFormationModel

    p = regime["params"]
    model = LaneFormationModel(
        n_agents=p["n_agents"],
        corridor_width=p["corridor_width"],
        corridor_height=p["corridor_height"],
        desired_speed=p["desired_speed"],
        repulsion_amplitude=p["repulsion_amplitude"],
        repulsion_range=p["repulsion_range"],
        tau=p["tau"],
        dt=p["dt"],
        seed=regime["seed"],
    )
    history = model.run(n_steps=p["n_steps"])

    if p.get("single_population", False):
        # Override labels to all-same direction → no counterflow.
        for state in history:
            state["labels"] = np.zeros_like(state["labels"])

    return history
