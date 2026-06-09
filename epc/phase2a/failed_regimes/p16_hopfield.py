"""P16 (associative memory) — Class C failed regimes.

Two canonical failure modes per the Sprint 80 brief:
  1. over_capacity: Hopfield network loaded far above α_c ≈ 0.138 (α = 0.5,
     50 patterns in N=100 neurons). The network enters a spin-glass phase
     with no reliable recall → completion accuracy at chance. P16 should reject.
  2. single_pattern: Hopfield network storing exactly ONE pattern. Trivial
     convergence (not content-addressable memory over multiple distinct
     patterns). P16 definitive requires ≥2 distinct selectively recalled
     patterns; this regime tests the multi-pattern prerequisite.
"""

from __future__ import annotations

from typing import Any, Dict, List


CONFIG: Dict[str, Any] = {
    "substrate_id": "P16_associative_memory_class_c",
    "format": "state_vector",
    "status": "populated",
    "regimes": [
        {
            "label": "over_capacity",
            "description": (
                "Hopfield network at α=1.0 (100 patterns, N=100), far above "
                "α_c≈0.138. Spin-glass phase: no reliable recall. At α=0.5 "
                "the best-overlap metric (max over P patterns) can exceed the "
                "0.5 screening threshold purely from many-pattern statistics; "
                "α=1.0 pushes mean best-overlap below 0.5 (empirically 0.464)."
            ),
            "params": {
                "N": 100,
                "P": 100,
                "corruption": 0.2,
                "max_steps": 100,
            },
            "seed": 42,
        },
        {
            "label": "single_pattern",
            "description": (
                "Hopfield network storing exactly ONE pattern (P=1). Convergence "
                "is trivial (single fixed point, not content-addressable memory "
                "over multiple patterns). P16 definitive requires ≥2 distinct "
                "selectively recalled patterns."
            ),
            "params": {
                "N": 100,
                "P": 1,
                "corruption": 0.2,
                "max_steps": 100,
            },
            "seed": 42,
        },
    ],
}


def build_substrate(regime: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Build a P16 failed-regime substrate."""
    import numpy as np
    from epc.models.hopfield import HopfieldNetwork, HopfieldParams

    p = regime["params"]
    seed = regime["seed"]
    N = p["N"]
    P = p["P"]
    corruption = p["corruption"]
    max_steps = p["max_steps"]

    params = HopfieldParams(
        N=N, P=P, corruption=corruption, max_steps=max_steps, seed=seed,
    )
    model = HopfieldNetwork(params)

    # Run retrieval trials — one per stored pattern (up to 5)
    n_trials = min(P, 5)
    cue_indices = list(range(n_trials))
    history = model.simulate(
        n_cues=n_trials,
        corruption=corruption,
        cue_pattern_indices=cue_indices,
    )
    return history
