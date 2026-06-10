"""P32 (specialization / division of labor) failed regimes for Phase-2a panel.

Class C for P32 covers two failure modes:

1. No reinforcement (ξ=0, φ=0): thresholds never update → agents respond
   randomly based on fixed identical thresholds; entropy does not decline
   → P32 screening rejects.

2. High forgetting, low reinforcement: thresholds drift uniformly upward
   faster than reinforcement can lock agents into roles → agents switch
   tasks frequently; entropy remains high → P32 screening/confirmation
   rejects.

Both regimes use the existing division-of-labor models with modified parameters.
"""

from __future__ import annotations

from typing import Any, Dict, List

import numpy as np


CONFIG: Dict[str, Any] = {
    "substrate_id": "P32_specialization_failed",
    "format": "task_allocation_timeseries",
    "status": "populated",
    "description": (
        "10 P32 failed regimes: 5 with no reinforcement "
        "(NoReinforcementModel, ξ=0, φ=0 → no entropy decline) "
        "+ 5 with high forgetting / low reinforcement "
        "(ResponseThresholdModel, φ=0.2, ξ=0.005–0.02 → thresholds "
        "drift uniformly, agents never lock into roles). "
        "All regimes should be rejected by P32: no entropy decline."
    ),
    "regimes": [
        # Regimes 0-4: no reinforcement (NoReinforcementModel — identical
        # thresholds throughout, agents respond stochastically).
        *[
            {
                "label": f"no_reinforcement_seed={100 + i}",
                "params": {
                    "model": "no_reinforcement",
                    "n_agents": 20,
                    "n_tasks": 3,
                    "stimulus_rate": 0.1,
                    "initial_threshold": 0.5,
                    "n_steps": 500,
                },
                "seed": 100 + i,
            }
            for i in range(5)
        ],
        # Regimes 5-9: high forgetting, low reinforcement (ResponseThresholdModel
        # with φ=0.2 >> ξ → forgetting dominates, no stable roles emerge).
        *[
            {
                "label": f"high_forgetting_reinf={er:.3f}",
                "params": {
                    "model": "response_threshold",
                    "n_agents": 20,
                    "n_tasks": 3,
                    "reinforcement_rate": float(er),
                    "forgetting_rate": 0.2,
                    "stimulus_rate": 0.1,
                    "initial_threshold": 0.5,
                    "n_steps": 500,
                },
                "seed": 200 + i,
            }
            for i, er in enumerate(np.linspace(0.005, 0.02, 5))
        ],
    ],
}


def build_substrate(regime: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Run a failed P32 regime and return task-allocation state history."""
    p = regime["params"]
    model_type = p["model"]

    if model_type == "no_reinforcement":
        from epc.models.division_of_labor import NoReinforcementModel
        model = NoReinforcementModel(
            n_agents=p["n_agents"],
            n_tasks=p["n_tasks"],
            stimulus_rate=p["stimulus_rate"],
            initial_threshold=p["initial_threshold"],
            seed=regime["seed"],
        )
        model.setup()
        history: List[Dict[str, Any]] = []
        for _ in range(p["n_steps"]):
            history.append(model.step())
        return history

    if model_type == "response_threshold":
        from epc.models.division_of_labor import ResponseThresholdModel
        model = ResponseThresholdModel(
            n_agents=p["n_agents"],
            n_tasks=p["n_tasks"],
            reinforcement_rate=p["reinforcement_rate"],
            forgetting_rate=p["forgetting_rate"],
            stimulus_rate=p["stimulus_rate"],
            initial_threshold=p["initial_threshold"],
            seed=regime["seed"],
        )
        model.setup()
        history = []
        for _ in range(p["n_steps"]):
            history.append(model.step())
        return history

    raise ValueError(f"unknown model type: {model_type}")
