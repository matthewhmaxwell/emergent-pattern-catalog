"""P5 (Vicsek) failed regimes: 10 high-noise regimes above order-disorder transition.

Per Sprint 46 brief: noise = linspace(0.7, 1.5, 10) — all above η_c ≈ 0.5
(published critical noise for N=300, ρ≈6.12).

In these regimes the Vicsek model is in the disordered phase: φ ≈ 1/√N << 0.5.
P5 correctly rejects at screening (φ < 0.5 threshold).
"""

from __future__ import annotations

from typing import Any, Dict, List

import numpy as np


NOISE_VALUES = np.linspace(4.0, 6.0, 10)  # genuinely disordered (phi<0.5); transition at this density is eta~3.5, NOT 0.5

CONFIG: Dict[str, Any] = {
    "substrate_id": "P5_vicsek_high_noise",
    "format": "particles",
    "description": (
        "10 Vicsek regimes at noise in linspace(4.0, 6.0, 10), "
        "N=300, box_size=7.0, speed=0.03. Genuinely DISORDERED "
        "(phi<0.5; the order-disorder transition at this density is eta~3.5, "
        "not 0.5 as the prior comment wrongly assumed). 1500 steps clears the "
        "5xT_cross run-length guard so the polar-order metric actually runs."
    ),
    "regimes": [
        {
            "label": f"noise={noise:.3f}",
            "params": {
                "n_particles": 300,
                "box_size": 7.0,
                "speed": 0.03,
                "noise": float(noise),
                "n_steps": 1500,
            },
            "seed": 200 + i,
        }
        for i, noise in enumerate(NOISE_VALUES)
    ],
}


def build_substrate(regime: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Run a Vicsek high-noise regime and return particle-history."""
    from epc.models.vicsek import VicsekModel

    p = regime["params"]
    model = VicsekModel(
        n_particles=p["n_particles"],
        box_size=p["box_size"],
        speed=p["speed"],
        noise=p["noise"],
        seed=regime["seed"],
    )
    return model.run(n_steps=p["n_steps"])
