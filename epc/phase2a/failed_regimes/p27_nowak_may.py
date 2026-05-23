"""P27 (Nowak-May) failed regimes: 10 sub-cooperation extinction regimes.

Per Sprint 39 brief: b = linspace(2.0, 2.5, 10) — these are extinction
regimes where defector temptation is large enough that cooperation collapses
and cooperators go extinct.

Nowak & May (1992) establish b=2.0 as the cooperation extinction threshold:
for b ≥ 2, defectors dominate unconditionally on a finite lattice with Moore
neighbourhood + synchronous imitation update (T > 2R → always defect).
All regimes in this sweep should produce f_C → 0, making them strong negatives
for the P27 detector (which requires f_C > 0.02 after 100 generations).
"""

from __future__ import annotations

from typing import Any, Dict, List

import numpy as np


ROWS = 50
COLS = 50
N_STEPS = 200
INIT_COOP_FRACTION = 0.5

B_VALUES = list(np.linspace(2.0, 2.5, 10))


CONFIG: Dict[str, Any] = {
    "substrate_id": "P27_nowak_may_extinction",
    "format": "grid",
    "description": (
        "10 Nowak-May regimes at b ∈ linspace(2.0, 2.5, 10) on 50×50 lattice, "
        "n_steps=200, init_coop_fraction=0.5. "
        "b ≥ 2.0 is the cooperation extinction regime (T > 2R → defectors dominate). "
        "f_C → 0 in all regimes; P27 should reject."
    ),
    "regimes": [
        {
            "label": f"b={b:.3f}",
            "params": {
                "rows": ROWS,
                "cols": COLS,
                "b": float(b),
                "init_coop_fraction": INIT_COOP_FRACTION,
                "n_steps": N_STEPS,
            },
            "seed": 200 + i,
        }
        for i, b in enumerate(B_VALUES)
    ],
}


def build_substrate(regime: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Run a Nowak-May extinction regime and return grid history."""
    from epc.models.nowak_may import NowakMayModel

    p = regime["params"]
    model = NowakMayModel(
        rows=p["rows"],
        cols=p["cols"],
        b=p["b"],
        init_coop_fraction=p["init_coop_fraction"],
        seed=regime["seed"],
    )
    return model.run(n_steps=p["n_steps"])
