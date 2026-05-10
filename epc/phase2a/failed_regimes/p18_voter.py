"""P18 (voter) failed regimes: 10 highly-biased-init voter runs.

The voter model on a finite lattice always reaches consensus eventually
— there is no parameter value that suppresses the consensus end-state.
However, the P18 detector's screening signature is the *coarsening
dynamics* (rapid early growth of Moran's I, persistent wall-density
decay, ≥30% Moran plateau). When the system is initialised with a
strong bias toward one opinion (init_fraction ≥ ~0.93), there are no
substantial domains to coarsen — minority cells flip immediately into
the majority and the trajectory shows no Moran growth.

Empirically the detector rejects (screening fail) for init_fraction
≥ 0.93 on L=64 (sub-threshold ``moran_spearman_early`` ≈ 0.50; below
the gate of 0.70). Below ~0.91 the detector still fires — these are
genuine canonical voter runs, just with biased starts.

The 10 regimes here use init_fraction ∈ linspace(0.93, 0.999, 10),
giving a clean failed-regime sweep along the model's only available
"no-coarsening" axis. See ``docs/sprint_returns/sprint_30_return.md``
for the design discussion.
"""

from __future__ import annotations

from typing import Any, Dict, List

import numpy as np


L = 64
N_STEPS = 200
INIT_FRACTIONS = list(np.linspace(0.93, 0.999, 10))


CONFIG: Dict[str, Any] = {
    "substrate_id": "P18_voter_biased_init",
    "format": "grid",
    "description": (
        "10 voter regimes at L=64, n_steps=200, init_mode=biased with "
        "init_fraction ∈ linspace(0.93, 0.999, 10). High bias suppresses "
        "the early Moran's I growth that P18 screening keys on; the "
        "system reaches trivial consensus without exhibiting coarsening "
        "dynamics. Voter has no parameter regime that suppresses "
        "consensus itself, so we vary the initial-condition regime."
    ),
    "regimes": [
        {
            "label": f"init_frac={frac:.3f}",
            "params": {"rows": L, "cols": L, "init_mode": "biased", "init_fraction": float(frac), "n_steps": N_STEPS},
            "seed": 200 + i,
        }
        for i, frac in enumerate(INIT_FRACTIONS)
    ],
}


def build_substrate(regime: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Run a biased-init voter regime and return grid-history."""
    from epc.models.voter import VoterModel

    p = regime["params"]
    model = VoterModel(
        rows=p["rows"], cols=p["cols"],
        init_mode=p["init_mode"], init_fraction=p["init_fraction"],
        seed=regime["seed"],
    )
    return model.run(n_steps=p["n_steps"])
