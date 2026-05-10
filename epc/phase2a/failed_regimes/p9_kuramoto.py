"""P9 (Kuramoto) failed regimes: 10 sub-critical K values.

Per Sprint 30 brief: K evenly spaced in [0.05·K_c, 0.5·K_c] where K_c = 2γ.
With γ=0.5 this is K_c = 1.0, so K ∈ [0.05, 0.50].

In this regime the Kuramoto model does not synchronize: order parameter
r ≈ 1/√N (incoherent), so a specific synchronization detector should
reject. These are the strongest negatives — same model, sub-critical
parameter — so a low TNR here would indicate the detector is firing on
something other than the synchronization signature.
"""

from __future__ import annotations

from typing import Any, Dict, List

import numpy as np


GAMMA = 0.5
K_C = 2 * GAMMA  # Kuramoto critical coupling for Lorentzian γ
N = 200
N_STEPS = 6000  # ≥10 oscillation periods at γ=0.5, dt=0.05 (T_osc ≈ 251 steps)
RECORD_EVERY = 10
DT = 0.05

K_VALUES = list(np.linspace(0.05 * K_C, 0.5 * K_C, 10))

CONFIG: Dict[str, Any] = {
    "substrate_id": "P9_kuramoto_subcritical",
    "format": "phases",
    "description": (
        "10 Kuramoto regimes at K ∈ linspace(0.05·K_c, 0.5·K_c) "
        "with γ=0.5, N=200, dt=0.05. K_c = 2γ = 1.0."
    ),
    "regimes": [
        {
            "label": f"K={K:.3f}",
            "params": {
                "N": N,
                "K": float(K),
                "gamma": GAMMA,
                "dt": DT,
                "n_steps": N_STEPS,
                "record_every": RECORD_EVERY,
            },
            "seed": 100 + i,
        }
        for i, K in enumerate(K_VALUES)
    ],
}


def build_substrate(regime: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Run a Kuramoto sub-critical regime and return phase-history."""
    from epc.models.kuramoto import KuramotoModel, KuramotoParams

    p = regime["params"]
    model = KuramotoModel(KuramotoParams(
        N=p["N"], K=p["K"], gamma=p["gamma"], dt=p["dt"], seed=regime["seed"],
    ))
    return model.run(n_steps=p["n_steps"], record_every=p["record_every"])
