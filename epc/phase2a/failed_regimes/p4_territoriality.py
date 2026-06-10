"""P4 (territoriality) failed regimes: overlapping-range + fast-decay regimes.

Class C for P4 covers two failure modes per the sprint 87 brief:

1. No scent response (overlapping ranges): agents deposit scent but do NOT
   avoid foreign scent, so home ranges overlap freely. Exclusivity ≈ 1/N
   for N agents → P4 screening rejects.

2. Scent decay too fast: agents respond to foreign scent but it decays so
   quickly that no persistent boundaries form. Boundary persistence drops
   below the confirmation threshold → P4 rejects at screening or fails
   confirmation.

Both regimes use the ScentMarkingModel with modified parameters. The model
itself is identical to the canonical positive; only the parameter values
differ (no code-level failure mode).
"""

from __future__ import annotations

from typing import Any, Dict, List

import numpy as np


# --- Regime 1: no scent response (repulsion_strength → 0) -------------------
# At repulsion_strength = 0, the foreign-scent threshold is 0 × deposition_rate = 0,
# so ALL cells with ANY foreign scent are blocked. This is too aggressive — agents
# can't move, producing trivially high exclusivity from frozen initial positions.
# Instead, we set repulsion_strength very HIGH (→ ∞), meaning the threshold
# foreign > repulsion_strength × deposition_rate is never met (the agent tolerates
# even massive foreign scent). Agents walk freely through foreign territory.
#
# Regimes 0-4: increasing repulsion_strength from 100 to 500 (all effectively
# infinite tolerance → overlapping ranges). home_attraction=0 disables own-scent
# bias, ensuring agents do genuine random walks (Sprint 87: with home_attraction>0,
# agents self-organize via own-scent attraction even without foreign avoidance,
# producing territories that are correctly detected).

# --- Regime 2: scent decay too fast -----------------------------------------
# At decay_rate → 1.0, scent disappears within 1 step. Foreign scent never
# accumulates enough to block movement, and boundaries can't persist.
#
# Regimes 5-9: increasing decay_rate from 0.5 to 0.95.

N_AGENTS = 4
GRID_SIZE = 32
N_STEPS = 3000
SNAPSHOT_INTERVAL = 30

CONFIG: Dict[str, Any] = {
    "substrate_id": "P4_territoriality_failed",
    "format": "territorial_agent_field",
    "status": "populated",
    "description": (
        "10 P4 failed regimes: 5 with infinite tolerance (repulsion_strength 100-500, "
        "agents ignore foreign scent → overlapping ranges) + 5 with fast decay "
        "(decay_rate 0.5-0.95, scent vanishes before boundaries form). "
        "All regimes should be rejected by P4: either low exclusivity (overlapping) "
        "or low boundary persistence (fast decay)."
    ),
    "regimes": [
        # Regimes 0-4: overlapping ranges (high tolerance, low scent response)
        *[
            {
                "label": f"high_tolerance_repulsion={rs:.0f}",
                "params": {
                    "n_agents": N_AGENTS,
                    "grid_size": GRID_SIZE,
                    "deposition_rate": 0.1,
                    "decay_rate": 0.03,
                    "repulsion_strength": float(rs),
                    "home_attraction": 0.0,  # No own-scent attraction either
                    "temperature": 0.5,
                    "n_steps": N_STEPS,
                    "snapshot_interval": SNAPSHOT_INTERVAL,
                },
                "seed": 200 + i,
            }
            for i, rs in enumerate(np.linspace(100, 500, 5))
        ],
        # Regimes 5-9: fast scent decay (no persistent boundaries)
        *[
            {
                "label": f"fast_decay_rate={dr:.2f}",
                "params": {
                    "n_agents": N_AGENTS,
                    "grid_size": GRID_SIZE,
                    "deposition_rate": 0.1,
                    "decay_rate": float(dr),
                    "repulsion_strength": 2.0,
                    "home_attraction": 2.0,
                    "temperature": 0.5,
                    "n_steps": N_STEPS,
                    "snapshot_interval": SNAPSHOT_INTERVAL,
                },
                "seed": 300 + i,
            }
            for i, dr in enumerate(np.linspace(0.5, 0.95, 5))
        ],
    ],
}


def build_substrate(regime: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Run a failed P4 regime and return territorial-agent-field history."""
    from epc.models.territoriality import ScentMarkingModel, ScentMarkingParams

    p = regime["params"]
    params = ScentMarkingParams(
        n_agents=p["n_agents"],
        grid_size=p["grid_size"],
        deposition_rate=p["deposition_rate"],
        decay_rate=p["decay_rate"],
        repulsion_strength=p["repulsion_strength"],
        home_attraction=p["home_attraction"],
        temperature=p["temperature"],
        n_steps=p["n_steps"],
        snapshot_interval=p["snapshot_interval"],
        seed=regime["seed"],
    )
    model = ScentMarkingModel(params)
    return model.simulate()
