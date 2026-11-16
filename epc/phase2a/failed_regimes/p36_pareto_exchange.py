"""P36 (heterogeneous kinetic exchange -> Pareto wealth tail) -- Class C failed regimes.

The target's own model, ``epc.models.kinetic_exchange`` (Chakraborti-Chakrabarti
kinetic exchange, N=3000, 5000 sweeps, 40 frames, the positive's protocol), with
the heterogeneity that produces the Pareto tail removed or bounded.

  homogeneous lambda=0.00   Boltzmann-Gibbs exponential (Dragulescu & Yakovenko
                            2000); standalone neg_boltzmann, reused
  homogeneous lambda=0.50   Gamma-like distribution (Chakraborti & Chakrabarti
                            2000; Patriarca, Chakraborti & Kaski 2004); standalone
                            neg_uniform_saving, reused
  homogeneous lambda=0.20 / 0.80   added; same literature result (Gamma, no
                            power-law tail) at a lower / higher saving propensity
  bounded heterogeneity lambda_i ~ U[0, 0.5)   added; the power-law tail needs a
                            non-vanishing density of saving propensities near 1
                            (Chatterjee, Chakrabarti & Manna 2004; Mohanty 2006,
                            PRE 74 011117), so bounding lambda well below 1 removes
                            it while keeping the agents heterogeneous

Two seeds per regime (200-209), distinct from the positive seeds (0-4).
"""
from __future__ import annotations

from typing import Any, Dict, List

import numpy as np

_KE = {"N": 3000, "n_sweeps": 5000, "n_frames": 40}

CONFIG: Dict[str, Any] = {
    "substrate_id": "P36_pareto_exchange_class_c",
    "format": "scalar_wealth",
    "status": "populated",
    "regimes": (
        [{"label": f"homogeneous_lambda={lam:.2f}",
          "params": dict(_KE, hetero=False, lam_fixed=lam), "seed": 200 + 2 * i + k}
         for i, lam in enumerate((0.0, 0.2, 0.5, 0.8)) for k in (0, 1)]
        + [{"label": "bounded_heterogeneity_lambda~U[0,0.5)",
            "params": dict(_KE, hetero=True, lam_max=0.5), "seed": 208 + k} for k in (0, 1)]
    ),
}


def build_substrate(regime: Dict[str, Any]) -> List[Dict[str, Any]]:
    from epc.models.kinetic_exchange import _cc_sweeps, kinetic_exchange
    p = regime["params"]
    seed = regime["seed"]
    if not p["hetero"]:
        return kinetic_exchange(seed, N=p["N"], hetero=False, lam_fixed=p["lam_fixed"],
                                n_sweeps=p["n_sweeps"], n_frames=p["n_frames"])[0]
    # Same construction as kinetic_exchange(hetero=True), with the range bounded.
    rng = np.random.default_rng(seed * 7919 + 1)
    lam = rng.uniform(0.0, p["lam_max"], size=p["N"])
    return _cc_sweeps(seed, p["N"], lam, p["n_sweeps"], p["n_frames"])
