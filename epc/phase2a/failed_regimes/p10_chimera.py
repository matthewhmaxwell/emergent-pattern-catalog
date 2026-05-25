"""P10 (chimera) Class C failed regimes: 10 ordinary all-to-all Kuramoto regimes.

Per Sprint 47 brief: K range above K_c on regular (all-to-all) Kuramoto, no
non-local coupling.  These produce uniform synchronisation (r→1 as K→∞), NOT
chimera states.

Discriminator per P10 detector docstring (ADR 50):
  - Chimera: pos_vel_ac[lag=4] ≈ 0.93 ± 0.01 (spatial coupling order)
  - Ordinary Kuramoto K=1.0: pos_vel_ac[lag=4] ≈ 0.31 ± 0.13

At K > K_c=1.0, all-to-all Kuramoto fully synchronises; pos_vel_ac[lag=4]
drops further below the screening floor of 0.55 → P10 rejects at screening.
Additionally, at full sync there are no incoherent windows → coexistence gate
(n_persistent_incoh ≥ 1) fails → screening rejection.

ADR 52: P10 DEFINITIVE additionally requires ``has_nonlocal_coupling=True``
in metadata, which ordinary Kuramoto does NOT provide.
"""

from __future__ import annotations

from typing import Any, Dict, List

import numpy as np


GAMMA = 0.5
K_C = 2 * GAMMA     # = 1.0 for Lorentzian distribution
N = 128             # same N as canonical chimera positive
N_STEPS = 6000      # total steps
RECORD_EVERY = 10   # → 600 frames; burn_fraction=0.30 in P10 → 420 post-burn
DT = 0.05
SEED_OFFSET = 300

# 10 K values above K_c — all produce synchronisation, NOT chimera.
K_VALUES = list(np.linspace(1.5 * K_C, 4.0 * K_C, 10))  # [1.5, 4.0]

CONFIG: Dict[str, Any] = {
    "substrate_id": "P10_ordinary_kuramoto_sync",
    "format": "phases",
    "description": (
        f"10 all-to-all Kuramoto regimes at K ∈ linspace(1.5·K_c, 4.0·K_c) "
        f"with γ={GAMMA}, N={N}, dt={DT}. K_c = 2γ = {K_C}. "
        "All produce full/partial synchronisation (r→1), NOT chimera coexistence. "
        "P10 rejects at screening via pos_vel_ac < 0.55 and/or coexistence gate."
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
            "seed": SEED_OFFSET + i,
        }
        for i, K in enumerate(K_VALUES)
    ],
}


def build_substrate(regime: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Run ordinary all-to-all Kuramoto above K_c and return phase history.

    Returns a list of dicts with ``theta`` (N,), ``r`` (float), ``step`` (int)
    keys — the native format that P10 accepts for processing.  Above K_c these
    produce synchronisation (r → 1), NOT chimera coexistence → P10 rejects.
    """
    from epc.models.kuramoto import KuramotoModel, KuramotoParams

    p = regime["params"]
    model = KuramotoModel(KuramotoParams(
        N=p["N"], K=p["K"], gamma=p["gamma"], dt=p["dt"], seed=regime["seed"],
    ))
    return model.run(n_steps=p["n_steps"], record_every=p["record_every"])
