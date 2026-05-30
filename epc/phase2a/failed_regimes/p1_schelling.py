"""P1 (Schelling) failed regimes: 10 sub-threshold segregation regimes.

Sprint 61 correction (was: linspace(0.05, 0.25, 10) on 32×32, density=0.9):
the empirical critical segregation threshold at density=0.9 is ≈0.13–0.15,
NOT 0.375 as documented for lower densities. Thresholds 0.161–0.250 at
density=0.9 genuinely segregate — they were true positives mislabeled as
negatives (brief-author error, same class as Sprint 40 P22 fix).

Corrected: threshold ∈ linspace(0.01, 0.10, 10), grid_size=50 (reduces
finite-size random-clustering noise that caused marginal FPs at 32×32),
density=0.9. At these thresholds, all agents are content in virtually
any mixed neighbourhood → no movement → system stays at random initial
state → Moran's I ≈ 0.01 (well below screening floor of 0.05).

Schelling (1971) established that segregation emerges above a critical
tolerance threshold that depends on density and neighbourhood size. At
density=0.9 on Moore neighbourhood, the critical threshold is empirically
≈0.13. Below that, agents are tolerant enough that the random initial
configuration is a stable equilibrium.
"""

from __future__ import annotations

from typing import Any, Dict, List

import numpy as np


GRID_SIZE = 50
DENSITY = 0.9
N_STEPS = 80

THRESHOLD_VALUES = list(np.linspace(0.01, 0.10, 10))


CONFIG: Dict[str, Any] = {
    "substrate_id": "P1_schelling_subthreshold",
    "format": "grid",
    "description": (
        "10 Schelling segregation regimes at threshold ∈ linspace(0.01, 0.10, 10) "
        f"on {GRID_SIZE}×{GRID_SIZE} grid, density={DENSITY}, n_steps={N_STEPS}. "
        "All values below empirical critical segregation threshold ≈0.13 at density=0.9; "
        "agents tolerate mixed neighbourhoods, system stays well-mixed, Moran's I stays low → P1 rejects."
    ),
    "regimes": [
        {
            "label": f"threshold={t:.3f}",
            "params": {
                "grid_size": GRID_SIZE,
                "density": DENSITY,
                "threshold": float(t),
                "n_steps": N_STEPS,
            },
            "seed": 100 + i,
        }
        for i, t in enumerate(THRESHOLD_VALUES)
    ],
}


def build_substrate(regime: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Run a sub-threshold Schelling regime and return grid history."""
    from epc.models.schelling import run_schelling

    p = regime["params"]
    return run_schelling(
        grid_size=p["grid_size"],
        density=p["density"],
        threshold=p["threshold"],
        n_steps=p["n_steps"],
        seed=regime["seed"],
    )
