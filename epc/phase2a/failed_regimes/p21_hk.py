"""P21 (HK polarization) — Class C failed regimes: high-ε consensus regimes.

Hegselmann-Krause (2002) converges to CONSENSUS at ε ≥ 0.5 (all agents interact,
yielding a single opinion cluster).  For ε ∈ [0.45, 0.60] the system either
reaches full consensus or produces very weak fragmentation with at most 1–2
short-lived clusters that immediately collapse.

At ε = 0.50 (the published threshold from HK 2002), all 100 agents converge to
a single cluster within ~10–20 steps from uniform IC.  The P21 detector should
reject these regimes because:
  - n_clusters = 1 at steady state → dip_p is high (unimodal)
  - persistence < 50 (convergence to consensus happens within ~20 steps)
  - OR dip_p does not drop below 0.05 once consensus is reached

The 10 regimes sweep ε linearly from 0.45 to 0.60 (inclusive), testing the
detector's ability to reject near-threshold and well-above-threshold consensus
regimes.  Seeds 400–409 (distinct from other failed-regime seed ranges).
"""

from __future__ import annotations

import numpy as np

from typing import Any, Dict, List


def build_substrate(regime: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Run HK at high ε and return native opinions-format history.

    P21's ``detect_p21`` reads ``history[-1]["opinions"]`` directly, so we
    return the native HK history (list of dicts with ``opinions`` key) without
    adaptation.  The convergence-to-consensus dynamics produce unimodal opinion
    distributions that P21 should reject.
    """
    from epc.models.hegselmann_krause import HegselmannKrauseModel

    p = regime["params"]
    model = HegselmannKrauseModel(
        n_agents=p["n_agents"],
        epsilon=p["epsilon"],
        init_mode=p.get("init_mode", "uniform"),
        seed=regime["seed"],
    )
    return model.run(n_steps=p["n_steps"])


# 10 high-ε regimes: ε ∈ linspace(0.45, 0.60, 10) — consensus territory per HK 2002.
_EPSILONS = list(np.linspace(0.45, 0.60, 10))

CONFIG: Dict[str, Any] = {
    "substrate_id": "P21_hk_high_epsilon",
    "format": "opinions",
    "status": "populated",
    "regimes": [
        {
            "label": f"hk_eps_{eps:.4f}",
            "params": {
                "n_agents": 100,
                "epsilon": eps,
                "init_mode": "uniform",
                "n_steps": 200,
            },
            "seed": 400 + i,
        }
        for i, eps in enumerate(_EPSILONS)
    ],
}
