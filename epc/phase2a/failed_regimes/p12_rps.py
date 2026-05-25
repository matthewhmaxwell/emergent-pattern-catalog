"""P12 (spatial RPS) — Class C failed regimes: 10 high-mobility extinction regimes.

Per Sprint 44 brief: mobility = linspace(5e-3, 5e-2, 10) — above the
extinction threshold. Reichenbach et al. (2007) report a critical mobility
M_c ≈ (4.5 ± 0.5) × 10⁻⁴ for σ = µ = 1 on an N=100² lattice. At
mobility ≫ M_c the three-species coexistence state becomes unstable: one
species goes extinct first, then the remaining two equilibrate, collapsing
the cyclic dominance structure to a monoculture (or a two-species limit
that is not cyclic). The smallest value in this sweep (5e-3) is already
~11× M_c, so every regime is firmly in the extinction-collapse phase.

P12's confirmation tier requires all three species to persist above 5%
for ≥ 80% of the second-half trajectory. In a monoculture collapse, one
or two species vanish entirely → persistence fails → P12 rejects.
"""

from __future__ import annotations

from typing import Any, Dict, List

import numpy as np


ROWS = 50
COLS = 50
N_STEPS = 200

MOBILITY_VALUES = list(np.linspace(5e-3, 5e-2, 10))


CONFIG: Dict[str, Any] = {
    "substrate_id": "P12_rps_high_mobility_extinction",
    "format": "grid",
    "description": (
        "10 spatial RPS regimes at mobility ∈ linspace(5e-3, 5e-2, 10) on 50×50 "
        "lattice, n_steps=200. All values >> M_c ≈ 4.5×10⁻⁴ (Reichenbach 2007). "
        "Above M_c the cyclic coexistence state collapses: one species goes extinct, "
        "leaving a monoculture or two-species state with no cyclic dominance. "
        "P12 rejects when three-species persistence fails."
    ),
    "regimes": [
        {
            "label": f"mobility={m:.4f}",
            "params": {
                "rows": ROWS,
                "cols": COLS,
                "mobility": float(m),
                "n_steps": N_STEPS,
            },
            "seed": 120 + i,
        }
        for i, m in enumerate(MOBILITY_VALUES)
    ],
}


def build_substrate(regime: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Run a high-mobility RPS extinction regime and return grid history."""
    from epc.models.rps_spatial import RPSSpatialModel

    p = regime["params"]
    model = RPSSpatialModel(
        rows=p["rows"],
        cols=p["cols"],
        mobility=p["mobility"],
        seed=regime["seed"],
    )
    return model.run(n_steps=p["n_steps"])
