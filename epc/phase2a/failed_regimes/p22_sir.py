"""P22 (SIR) failed regimes: 10 sub-percolation epidemic regimes.

Per Sprint 39 brief: infection_prob = linspace(0.05, 0.18, 10) — these
are sub-percolation regimes where the epidemic fails to spread widely.
Below the effective cascade detection threshold (~0.2) the epidemic
dies out or remains too localised to produce a statistically detectable
information cascade, so the P22 detector should reject all regimes.

The SIR CA physical percolation threshold is p_c ≈ 0.038 (Moore, q=0.1).
At p ∈ [0.05, 0.18] the epidemic may technically percolate, but the
effective detection threshold for P22 (which requires ≥5% of susceptibles
infected AND significant spatial clustering Moran's I > random) is higher:
cascades in this range are too small, slow, or spatially diffuse to satisfy
both criteria simultaneously. In practice these regimes produce tiny halos
around the initial seed that die out well before reaching the detection floor.
"""

from __future__ import annotations

from typing import Any, Dict, List

import numpy as np


ROWS = 64
COLS = 64
RECOVERY_PROB = 0.1

INFECTION_PROB_VALUES = list(np.linspace(0.05, 0.18, 10))


CONFIG: Dict[str, Any] = {
    "substrate_id": "P22_sir_subpercolation",
    "format": "grid",
    "description": (
        "10 SIR epidemic regimes at infection_prob ∈ linspace(0.05, 0.18, 10) "
        "on 64×64 lattice, recovery_prob=0.1, single_seed init. "
        "Below effective cascade detection threshold ~0.2; "
        "epidemic fails to spread widely enough for P22 to detect."
    ),
    "regimes": [
        {
            "label": f"p={p:.3f}",
            "params": {
                "rows": ROWS,
                "cols": COLS,
                "infection_prob": float(p),
                "recovery_prob": RECOVERY_PROB,
                "init_mode": "single_seed",
                "n_steps": 200,
            },
            "seed": 400 + i,
        }
        for i, p in enumerate(INFECTION_PROB_VALUES)
    ],
}


def build_substrate(regime: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Run a sub-percolation SIR regime and return grid history."""
    from epc.models.sir_epidemic import SIREpidemicModel

    p = regime["params"]
    model = SIREpidemicModel(
        rows=p["rows"],
        cols=p["cols"],
        infection_prob=p["infection_prob"],
        recovery_prob=p["recovery_prob"],
        init_mode=p["init_mode"],
        seed=regime["seed"],
    )
    return model.run(n_steps=p["n_steps"])
