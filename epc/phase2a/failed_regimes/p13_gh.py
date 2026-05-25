"""P13 (Greenberg-Hastings) — Class C failed regimes: 10 low-density init regimes.

Per Sprint 44 brief: init_density = linspace(0.01, 0.1, 10) — below the
excitation threshold needed for persistent spiral/target waves. In the
Greenberg-Hastings CA (Greenberg & Hastings 1978), spiral waves require a
sufficiently dense random initial seeding so that local clusters of excited
cells can trigger neighbors and produce self-sustaining propagating fronts.
At very low initial density, excited clusters are too sparse and isolated to
produce persistent wavefronts: each small cluster burns out before it can
nucleate a spiral, leaving the lattice quiescent after a few steps.

The canonical spiral-forming regime uses init_density ≈ 0.3 (30% of cells
initialized to random states). At densities in [0.01, 0.10] (1%–10%) the
excitable system fails to sustain wavefront propagation: wavefront_count
drops to ~0 within a few steps, and P13's screening prerequisite
(persistent wavefront for ≥ 5 × T_prop) is not satisfied → P13 rejects.
"""

from __future__ import annotations

from typing import Any, Dict, List

import numpy as np


ROWS = 50
COLS = 50
N_STEPS = 200
N_STATES = 8
THRESHOLD = 1

DENSITY_VALUES = list(np.linspace(0.01, 0.1, 10))


CONFIG: Dict[str, Any] = {
    "substrate_id": "P13_gh_low_density_no_waves",
    "format": "grid",
    "description": (
        "10 Greenberg-Hastings regimes at init_density ∈ linspace(0.01, 0.10, 10) "
        "on 50×50 lattice, n_states=8, threshold=1, moore neighbourhood, n_steps=200. "
        "Below the excitation threshold for persistent wavefronts: sparse isolated "
        "excited clusters burn out before sustaining spiral nucleation. "
        "P13 rejects when wavefront persistence fails at screening tier."
    ),
    "regimes": [
        {
            "label": f"density={d:.3f}",
            "params": {
                "rows": ROWS,
                "cols": COLS,
                "n_states": N_STATES,
                "threshold": THRESHOLD,
                "init_density": float(d),
                "n_steps": N_STEPS,
            },
            "seed": 130 + i,
        }
        for i, d in enumerate(DENSITY_VALUES)
    ],
}


def build_substrate(regime: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Run a low-density GH regime and return grid history."""
    from epc.models.greenberg_hastings import GreenbergHastings

    p = regime["params"]
    model = GreenbergHastings(
        rows=p["rows"],
        cols=p["cols"],
        n_states=p["n_states"],
        threshold=p["threshold"],
        neighborhood="moore",
        boundary="periodic",
        init_mode="random",
        init_density=p["init_density"],
        seed=regime["seed"],
    )
    return model.run(n_steps=p["n_steps"])
