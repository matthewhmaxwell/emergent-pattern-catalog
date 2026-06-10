"""P29 (trail / network formation) failed regimes for Phase-2a panel.

Class C for P29 covers two failure modes:

1. No reinforcement (evaporation too fast): pheromone evaporates before
   agents can exploit it, so no persistent network forms. Weight-distance
   correlation ≈ 0 → P29 screening rejects.

2. No reinforcement (uniform random walk): agents move randomly with no
   pheromone-proportional choice (NoReinforcementModel). Edge weights
   accumulate uniformly, not preferentially on short edges → P29 rejects.

Both regimes use the existing trail-network models with modified parameters.
"""

from __future__ import annotations

from typing import Any, Dict, List

import numpy as np


CONFIG: Dict[str, Any] = {
    "substrate_id": "P29_trail_network_failed",
    "format": "trail_network",
    "status": "populated",
    "description": (
        "10 P29 failed regimes: 5 with high evaporation (evaporation_rate 0.80–0.99, "
        "pheromone vanishes before reinforcement can accumulate → no efficient network) "
        "+ 5 with no reinforcement (NoReinforcementModel, uniform random choice → "
        "uniform edge weights, no weight-distance correlation). "
        "All regimes should be rejected by P29: low weight-distance correlation."
    ),
    "regimes": [
        # Regimes 0-4: high evaporation (AntTrailModel with evaporation so high
        # that pheromone cannot accumulate — no persistent reinforced network).
        *[
            {
                "label": f"high_evaporation_rate={er:.2f}",
                "params": {
                    "model": "ant_trail",
                    "n_nodes": 7,
                    "n_agents": 40,
                    "grid_size": 100,
                    "alpha": 2.0,
                    "beta": 3.0,
                    "deposition_rate": 10.0,
                    "evaporation_rate": float(er),
                    "n_steps": 500,
                    "snapshot_interval": 10,
                },
                "seed": 200 + i,
            }
            for i, er in enumerate(np.linspace(0.80, 0.99, 5))
        ],
        # Regimes 5-9: no reinforcement (NoReinforcementModel — agents choose
        # next node uniformly at random, ignoring pheromone).
        *[
            {
                "label": f"no_reinforcement_seed={300 + i}",
                "params": {
                    "model": "no_reinforcement",
                    "n_nodes": 7,
                    "n_agents": 40,
                    "grid_size": 100,
                    "n_steps": 500,
                    "snapshot_interval": 10,
                },
                "seed": 300 + i,
            }
            for i in range(5)
        ],
    ],
}


def build_substrate(regime: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Run a failed P29 regime and return trail-network history."""
    p = regime["params"]
    model_type = p["model"]

    if model_type == "ant_trail":
        from epc.models.trail_network import AntTrailModel, AntTrailParams
        params = AntTrailParams(
            n_nodes=p["n_nodes"],
            n_agents=p["n_agents"],
            grid_size=p["grid_size"],
            alpha=p["alpha"],
            beta=p["beta"],
            deposition_rate=p["deposition_rate"],
            evaporation_rate=p["evaporation_rate"],
            n_steps=p["n_steps"],
            snapshot_interval=p["snapshot_interval"],
            seed=regime["seed"],
        )
        model = AntTrailModel(params)
        return model.simulate()

    if model_type == "no_reinforcement":
        from epc.models.trail_network import NoReinforcementModel, NoReinforcementParams
        params = NoReinforcementParams(
            n_nodes=p["n_nodes"],
            n_agents=p["n_agents"],
            grid_size=p["grid_size"],
            n_steps=p["n_steps"],
            snapshot_interval=p["snapshot_interval"],
            seed=regime["seed"],
        )
        model = NoReinforcementModel(params)
        return model.simulate()

    raise ValueError(f"unknown model type: {model_type}")
