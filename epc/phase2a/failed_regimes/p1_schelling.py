"""P1 (Schelling) failed regimes: 10 sub-threshold segregation regimes.

Per Sprint 42 brief: threshold = linspace(0.05, 0.25, 10) on a 32×32 grid.
All values are strictly below the canonical Schelling segregation threshold
of 0.375. At sub-threshold tolerance values, agents accept most mixed
neighbourhoods and the system remains well-mixed without forming clusters.
Moran's I stays low; P1 should reject.

Schelling (1971) established that segregation emerges above a critical
tolerance threshold (~0.3–0.4 for typical neighbourhood and density).
Below the threshold, tolerant agents have no incentive to move and the
initial (random) mixed configuration is stable.
"""

from __future__ import annotations

from typing import Any, Dict, List

import numpy as np


GRID_SIZE = 32
DENSITY = 0.9
N_STEPS = 80

THRESHOLD_VALUES = list(np.linspace(0.05, 0.25, 10))


CONFIG: Dict[str, Any] = {
    "substrate_id": "P1_schelling_subthreshold",
    "format": "grid",
    "description": (
        "10 Schelling segregation regimes at threshold ∈ linspace(0.05, 0.25, 10) "
        f"on {GRID_SIZE}×{GRID_SIZE} grid, density={DENSITY}, n_steps={N_STEPS}. "
        "All values below canonical segregation threshold 0.375; agents tolerate "
        "mixed neighbourhoods, system stays well-mixed, Moran's I stays low → P1 rejects."
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
