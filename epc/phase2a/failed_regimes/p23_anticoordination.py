"""P23 (anti-coordination) — Class C failed regimes.

Two canonical failure modes per the Sprint 77 brief:
  1. random_agents: Random-choice agents with no strategy adaptation. Attendance
     follows Binomial(N, 0.5) → σ²/N ≈ 0.25 (random baseline). No anti-correlation.
     The detector should reject — variance is AT the random baseline, not below it.
  2. herding_regime: Majority-imitation agents that produce σ² ABOVE the random
     baseline (crowding, not balancing). Positive autocorrelation (persistence,
     not anti-persistence). The detector should reject — this is the opposite
     of anti-coordination.
"""

from __future__ import annotations

from typing import Any, Dict, List


CONFIG: Dict[str, Any] = {
    "substrate_id": "P23_anticoordination_class_c",
    "format": "choice_timeseries",
    "status": "populated",
    "regimes": [
        {
            "label": "random_agents",
            "description": (
                "Random-choice agents (no adaptation). Attendance ~ Binomial(N, 0.5). "
                "σ²/N ≈ 0.25 = random baseline. No temporal structure."
            ),
            "params": {
                "model": "random_choice",
                "n_agents": 101,
                "n_rounds": 3000,
            },
            "seed": 42,
        },
        {
            "label": "herding_regime",
            "description": (
                "Majority-imitation herding: agents copy majority with p=0.8. "
                "Attendance drifts toward consensus → σ²/N >> 0.25 (above random "
                "baseline). Positive autocorrelation (persistence). Opposite of "
                "anti-coordination."
            ),
            "params": {
                "model": "herding",
                "n_agents": 101,
                "n_rounds": 3000,
                "imitation_prob": 0.8,
            },
            "seed": 42,
        },
    ],
}


def build_substrate(regime: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Build a P23 failed-regime substrate."""
    import numpy as np

    p = regime["params"]
    seed = regime["seed"]
    rng = np.random.default_rng(seed)
    n_agents = p["n_agents"]
    n_rounds = p["n_rounds"]

    if p["model"] == "random_choice":
        # i.i.d. Binomial(N, 0.5) — no strategy adaptation.
        records: List[Dict[str, Any]] = []
        for t in range(n_rounds):
            att = int(rng.binomial(n_agents, 0.5))
            records.append({
                'round': t,
                'attendance': att,
                'n_agents': n_agents,
            })
        return records

    elif p["model"] == "herding":
        # Majority-imitation: agents copy majority with high probability.
        imitation_prob = p.get("imitation_prob", 0.8)
        att = n_agents // 2
        records = []
        for t in range(n_rounds):
            records.append({
                'round': t,
                'attendance': int(att),
                'n_agents': n_agents,
            })
            majority_side = 1 if att > n_agents // 2 else 0
            p_choose_1 = imitation_prob if majority_side == 1 else (1.0 - imitation_prob)
            att = int(rng.binomial(n_agents, p_choose_1))
        return records

    else:
        raise ValueError(f"Unknown model type: {p['model']}")
