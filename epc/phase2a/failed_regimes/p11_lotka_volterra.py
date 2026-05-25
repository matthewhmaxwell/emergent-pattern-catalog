"""P11 (Lotka-Volterra) failed regimes: 10 predator-extinction-rate regimes.

Per Sprint 45 brief: predator_death_rate = linspace(2.0, 5.0, 10) — at and
above this range predators die out and the system collapses to prey-only.

Mobilia-Georgiev-Täuber (2007) identify the coexistence regime as
λ >> μ (predation rate >> predator death rate). The canonical positive
uses λ=2.0, μ=1.0 (Mobilia 2007 Sec. III). When μ ≥ λ, predators die
faster than they can reproduce by predation; the absorbing (prey-only)
state is reached within O(100) generations. The stochastic lattice DP
extinction threshold is at μ ≈ λ (mean-field crossing), so all 10
regimes here (μ ∈ [2.0, 5.0], λ=2.0) are at or above the extinction
boundary.

P11's prerequisites fail in these regimes:
  - n_species prerequisite: after predator extinction only prey occupies
    the lattice (states 0=empty, 1=prey) → only 1 non-zero species
    observed → prerequisite fails.
  - Alternatively (early extinction, short history): predator fraction
    std ≈ 0 after absorption → species_variance prerequisite fails.

Either failure causes P11 to return detected=False with the appropriate
prerequisite-failure reason.
"""

from __future__ import annotations

from typing import Any, Dict, List

import numpy as np


ROWS = 50
COLS = 50
PREDATION_RATE = 2.0          # λ: canonical Mobilia 2007 value
PREY_REPRODUCTION_RATE = 1.0  # σ
N_STEPS = 400                  # enough for predators to go extinct at μ ≥ λ

DEATH_RATE_VALUES = list(np.linspace(2.0, 5.0, 10))


CONFIG: Dict[str, Any] = {
    "substrate_id": "P11_lv_predator_extinction",
    "format": "grid",
    "description": (
        "10 lattice Lotka-Volterra regimes at predator_death_rate ∈ "
        "linspace(2.0, 5.0, 10) on 50×50 lattice, n_steps=400, "
        "predation_rate=2.0 (canonical Mobilia 2007). "
        "At μ ≥ λ=2.0 the extinction threshold is reached: predators die "
        "faster than they reproduce, collapsing the system to prey-only. "
        "P11 rejects when n_species=1 (only prey) or species variance ≈ 0."
    ),
    "regimes": [
        {
            "label": f"mu={mu:.4f}",
            "params": {
                "rows": ROWS,
                "cols": COLS,
                "predation_rate": PREDATION_RATE,
                "prey_reproduction_rate": PREY_REPRODUCTION_RATE,
                "predator_death_rate": float(mu),
                "n_steps": N_STEPS,
            },
            "seed": 300 + i,
        }
        for i, mu in enumerate(DEATH_RATE_VALUES)
    ],
}


def build_substrate(regime: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Run a predator-extinction LV regime and return grid history."""
    from epc.models.lotka_volterra_lattice import LotkaVolterraLattice

    p = regime["params"]
    model = LotkaVolterraLattice(
        rows=p["rows"],
        cols=p["cols"],
        predation_rate=p["predation_rate"],
        prey_reproduction_rate=p["prey_reproduction_rate"],
        predator_death_rate=p["predator_death_rate"],
        seed=regime["seed"],
    )
    return model.run(n_steps=p["n_steps"])
