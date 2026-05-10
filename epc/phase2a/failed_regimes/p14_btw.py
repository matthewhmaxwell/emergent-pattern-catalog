"""P14 (BTW sandpile) failed regimes: 10 dissipative-sandpile sub-critical runs.

The brief calls these "sub-critical drive rate" regimes. In the original BTW
formulation drive rate doesn't gate criticality (the system self-organizes
regardless of slow-drive cadence). The canonical sub-critical null model
documented in the codebase is :func:`epc.models.btw_sandpile.run_dissipative_sandpile`,
which adds a per-toppling grain-loss probability ``p_diss`` that breaks the
power-law avalanche distribution. We use 10 dissipation rates evenly spaced
in [0.05, 0.5] — all within the documented "exponentially distributed
avalanches" regime where P14 should reject.

This is documented in
``tests/test_sandpile_p14_e2e.py::test_dissipative_negative`` as the
canonical P14 negative.
"""

from __future__ import annotations

from typing import Any, Dict, List

import numpy as np


L = 32
N_DRIVE = 10_000
N_BURN = 1_000
P_DISS_VALUES = list(np.linspace(0.05, 0.5, 10))


CONFIG: Dict[str, Any] = {
    "substrate_id": "P14_btw_dissipative",
    "format": "avalanches",
    "description": (
        "10 dissipative-sandpile regimes at L=32, n_drive=10_000, n_burn=1_000 "
        "with p_diss ∈ linspace(0.05, 0.5, 10). Bulk dissipation breaks SOC; "
        "avalanches become exponentially distributed."
    ),
    "regimes": [
        {
            "label": f"p_diss={p:.3f}",
            "params": {"L": L, "n_drive": N_DRIVE, "n_burn": N_BURN, "p_diss": float(p)},
            "seed": 300 + i,
        }
        for i, p in enumerate(P_DISS_VALUES)
    ],
}


def build_substrate(regime: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Run a dissipative sandpile and return its avalanche-size array."""
    from epc.models.btw_sandpile import run_dissipative_sandpile, BTWSandpileParams

    p = regime["params"]
    params = BTWSandpileParams(L=p["L"], n_drive=p["n_drive"], n_burn=p["n_burn"], seed=regime["seed"])
    result = run_dissipative_sandpile(params, p_diss=p["p_diss"])
    return [{"avalanche_sizes": result.avalanche_sizes, "step": 0}]
