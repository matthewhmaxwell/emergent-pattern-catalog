"""Class C failed-regime substrates for P30 (autopoiesis).

Three regimes where the autopoiesis model fails to produce a closed
self-maintaining membrane:

1. Non-bonding: production_rate=0 — no links ever form.
2. High-decay: decay_rate=0.5 — links decay faster than produced, no stable membrane.
3. No-attraction: catalyst_link_attraction=0 — links drift away, no closure.

All should be true negatives: P30 detector rejects at screening
(association_score ≈ 1.0 or closure_fraction ≈ 0).
"""

from __future__ import annotations

from typing import Any, Dict, List

CONFIG: Dict[str, Any] = {
    "status": "populated",
    "regimes": [
        {
            "label": "non_bonding_production_rate=0.00",
            "params": {
                "model": "autopoiesis",
                "n_substrate": 100,
                "n_catalyst": 3,
                "box_size": 20.0,
                "production_rate": 0.0,
                "decay_rate": 0.01,
                "membrane_equilibrium_radius": 3.0,
                "catalyst_link_attraction": 0.5,
                "n_steps": 200,
            },
            "seed": 200,
        },
        {
            "label": "high_decay_rate=0.50",
            "params": {
                "model": "autopoiesis",
                "n_substrate": 100,
                "n_catalyst": 3,
                "box_size": 20.0,
                "production_rate": 0.15,
                "decay_rate": 0.50,
                "membrane_equilibrium_radius": 3.0,
                "catalyst_link_attraction": 0.5,
                "n_steps": 200,
            },
            "seed": 201,
        },
        {
            "label": "no_attraction_high_decay",
            "params": {
                "model": "autopoiesis",
                "n_substrate": 100,
                "n_catalyst": 3,
                "box_size": 20.0,
                "production_rate": 0.02,
                "decay_rate": 0.20,
                "membrane_equilibrium_radius": 3.0,
                "catalyst_link_attraction": 0.0,
                "link_attraction": 0.0,
                "substrate_diffusion": 0.5,
                "n_steps": 200,
            },
            "seed": 202,
        },
        {
            "label": "weak_production_rate=0.01",
            "params": {
                "model": "autopoiesis",
                "n_substrate": 100,
                "n_catalyst": 3,
                "box_size": 20.0,
                "production_rate": 0.01,
                "decay_rate": 0.05,
                "membrane_equilibrium_radius": 3.0,
                "catalyst_link_attraction": 0.5,
                "n_steps": 200,
            },
            "seed": 203,
        },
        {
            "label": "large_box_low_production_size=80.0",
            "params": {
                "model": "autopoiesis",
                "n_substrate": 100,
                "n_catalyst": 3,
                "box_size": 80.0,
                "production_rate": 0.02,
                "decay_rate": 0.05,
                "membrane_equilibrium_radius": 3.0,
                "catalyst_link_attraction": 0.5,
                "n_steps": 200,
            },
            "seed": 204,
        },
    ],
}


def build_substrate(regime: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Build a failed-regime autopoiesis substrate."""
    from epc.models.autopoiesis import AutopoiesisModel

    p = regime["params"]
    seed = regime["seed"]
    n_steps = p.get("n_steps", 200)

    model = AutopoiesisModel(
        n_substrate=p.get("n_substrate", 100),
        n_catalyst=p.get("n_catalyst", 3),
        box_size=p.get("box_size", 20.0),
        production_rate=p.get("production_rate", 0.15),
        decay_rate=p.get("decay_rate", 0.01),
        production_radius=p.get("production_radius", 3.0),
        membrane_equilibrium_radius=p.get("membrane_equilibrium_radius", 3.0),
        catalyst_link_attraction=p.get("catalyst_link_attraction", 0.5),
        link_attraction=p.get("link_attraction", 0.3),
        substrate_diffusion=p.get("substrate_diffusion", 0.3),
        link_diffusion=p.get("link_diffusion", 0.01),
        catalyst_diffusion=p.get("catalyst_diffusion", 0.002),
        seed=seed,
    )
    model.setup()
    history: List[Dict[str, Any]] = []
    for _ in range(n_steps):
        history.append(model.step())
    return history
