"""P8 (Nagel-Schreckenberg) Class C failed regimes: 10 low-density / no-jamming.

Per Sprint 47 brief: density ∈ linspace(0.05, 0.20, 10), below / near jamming-
onset density.  At rho ≤ 0.10, stopped_fraction ≈ 0 → P8 rejects at screening
(stopped_fraction ≤ 0.05 floor).  At rho ≥ 0.12 (near / above onset at p=0.3),
some jams may form and some regimes may reach screening or confirmation tier.
Carry-forward C-p8-class-c-near-onset documents any FPs; see Sprint 47 return.

These are the SAME-MODEL strong negatives: density sweep on the canonical NS
model below the canonical jam regime (rho=0.30).  If the panel TNR < 1.0 on
this class, the detector is sensitive enough to flag near-onset behaviour.
"""

from __future__ import annotations

from typing import Any, Dict, List

import numpy as np


L = 1000          # ring length matching canonical positive
P_SLOW = 0.3      # randomization probability (NS 1992 canonical)
V_MAX = 5         # max velocity matching canonical positive
N_STEPS = 2500    # total steps: burn_in=1000 + 1500 measurement
SEED_OFFSET = 200

DENSITY_VALUES = list(np.linspace(0.05, 0.20, 10))

CONFIG: Dict[str, Any] = {
    "substrate_id": "P8_nagel_schreckenberg_low_density",
    "format": "sequence",
    "description": (
        "10 NS regimes at density ∈ linspace(0.05, 0.20, 10) with p_slow=0.3, "
        f"L={L}, v_max={V_MAX}. Below / near jamming-onset (~rho≈0.12 at p=0.3). "
        "Free-flow regimes (rho ≤ 0.10) have stopped_fraction ≈ 0 → P8 rejects "
        "at screening. Near-onset regimes may show partial jamming."
    ),
    "regimes": [
        {
            "label": f"rho={rho:.4f}",
            "params": {
                "L": L,
                "density": float(rho),
                "v_max": V_MAX,
                "p_slow": P_SLOW,
                "n_steps": N_STEPS,
            },
            "seed": SEED_OFFSET + i,
        }
        for i, rho in enumerate(DENSITY_VALUES)
    ],
}


def build_substrate(regime: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Run NS at low density and return the native model history.

    The history contains ``velocities`` (int8 array, per car per step) and
    ``positions`` keys — the native format that P8 requires for full detection.
    Low-density runs produce stopped_fraction ≈ 0 → P8 rejects at screening.
    """
    from epc.models.nagel_schreckenberg import NagelSchreckenberg

    p = regime["params"]
    model = NagelSchreckenberg(
        L=p["L"],
        density=p["density"],
        v_max=p["v_max"],
        p_slow=p["p_slow"],
        seed=regime["seed"],
    )
    return model.run(n_steps=p["n_steps"])
