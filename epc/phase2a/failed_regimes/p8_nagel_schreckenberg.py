"""P8 (Nagel-Schreckenberg) Class C failed regimes: 10 low-density / no-jamming.

Sprint 62 correction: density range narrowed from linspace(0.05, 0.20, 10) to
linspace(0.02, 0.07, 10).  The jamming onset for v_max=5, p_slow=0.3 at L=1000
is empirically ρ_c ≈ 0.10 (stopped_fraction jumps from ~0 at ρ=0.09 to ~0.045
at ρ=0.11, multi-seed validated).  The original range placed 6/10 regimes at
ρ ≥ 0.1167, well above onset — they genuinely jam and are mislabeled negatives
(brief-author error, same class as Sprint 40 P22 / Sprint 61 P1 corrections).

Corrected range [0.02, 0.07] sits entirely in the free-flow phase (max
stopped_fraction = 0.0 across all densities and seeds).  These are genuine
negatives: same-model runs that do NOT exhibit spontaneous jamming.
"""

from __future__ import annotations

from typing import Any, Dict, List

import numpy as np


L = 1000          # ring length matching canonical positive
P_SLOW = 0.3      # randomization probability (NS 1992 canonical)
V_MAX = 5         # max velocity matching canonical positive
N_STEPS = 2500    # total steps: burn_in=1000 + 1500 measurement
SEED_OFFSET = 200

DENSITY_VALUES = list(np.linspace(0.02, 0.07, 10))

CONFIG: Dict[str, Any] = {
    "substrate_id": "P8_nagel_schreckenberg_low_density",
    "format": "sequence",
    "description": (
        "10 NS regimes at density ∈ linspace(0.02, 0.07, 10) with p_slow=0.3, "
        f"L={L}, v_max={V_MAX}. All densities well below jamming onset "
        "(ρ_c ≈ 0.10 at v_max=5, p=0.3). Free-flow phase: stopped_fraction = 0 "
        "→ P8 rejects at screening. Sprint 62 correction from original [0.05, 0.20]."
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
