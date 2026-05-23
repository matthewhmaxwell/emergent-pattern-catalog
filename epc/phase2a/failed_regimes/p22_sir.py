"""P22 (SIR) failed regimes: 10 sub-percolation epidemic regimes.

Sprint 40 correction: infection_prob = linspace(0.005, 0.030, 10).
All values are strictly below the Moore-neighbourhood percolation threshold
p_c ≈ 0.038 (q=0.1, single_seed init). At sub-threshold infection
probabilities the epidemic from a single seed almost surely dies out
within a few steps; no global spread occurs.

Sprint 39 used linspace(0.05, 0.18, 10) — a brief-author error.
Those values are all above p_c ≈ 0.038, so every "failed" regime
was actually a genuine cascade, which caused Class C TNR = 0.000.
The fix is to use infection_prob ∈ [0.005, 0.030], entirely below the
percolation threshold, so the epidemic reliably fails to spread and
P22 correctly rejects all regimes.
"""

from __future__ import annotations

from typing import Any, Dict, List

import numpy as np


ROWS = 64
COLS = 64
RECOVERY_PROB = 0.1

INFECTION_PROB_VALUES = list(np.linspace(0.005, 0.030, 10))


CONFIG: Dict[str, Any] = {
    "substrate_id": "P22_sir_subpercolation",
    "format": "grid",
    "description": (
        "10 SIR epidemic regimes at infection_prob ∈ linspace(0.005, 0.030, 10) "
        "on 64×64 lattice, recovery_prob=0.1, single_seed init. "
        "All values below Moore-neighbourhood percolation threshold p_c ≈ 0.038 "
        "(q=0.1); epidemic from single seed dies out, no global spread → P22 rejects. "
        "Sprint 40 correction: Sprint 39 used [0.05, 0.18] which was above p_c, "
        "causing genuine cascades in the 'failed' class (brief-author error)."
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
