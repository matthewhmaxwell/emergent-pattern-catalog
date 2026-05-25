"""P6 (D'Orsogna SPP) failed regimes: 10 mismatched-radii regimes.

Per Sprint 46 brief: 10 high-noise / mismatched-radii regimes that suppress milling.

Milling requires attraction range (l_a) >> repulsion range (l_r).
Canonical regime: l_a=3.0, l_r=0.5 (l_a/l_r=6).

Mismatched-radii failed regimes: l_a ∈ linspace(0.10, 0.49, 10) with l_r=0.5 fixed.
At l_a ≤ l_r: short-range attraction cannot compete with repulsion across the ring;
the mill does not form. Angular momentum |L| ≈ 0 → P6 correctly rejects.
"""

from __future__ import annotations

from typing import Any, Dict, List

import numpy as np


LA_VALUES = list(np.linspace(0.10, 0.49, 10))

CONFIG: Dict[str, Any] = {
    "substrate_id": "P6_dorsogna_mismatched_radii",
    "format": "particles",
    "description": (
        "10 D'Orsogna SPP regimes with l_a ∈ linspace(0.10, 0.49, 10), l_r=0.5 fixed. "
        "Attraction range ≤ repulsion range suppresses milling (|L|≈0). "
        "Canonical milling regime uses l_a=3.0 >> l_r=0.5."
    ),
    "regimes": [
        {
            "label": f"l_a={la:.3f}",
            "params": {
                "n_particles": 100,
                "C_a": 0.5,
                "C_r": 1.0,
                "l_a": float(la),
                "l_r": 0.5,
                "alpha": 1.0,
                "beta": 0.5,
                "dt": 0.05,
                "init_mode": "random",
                "init_radius": 5.0,
                "n_steps": 500,
            },
            "seed": 400 + i,
        }
        for i, la in enumerate(LA_VALUES)
    ],
}


def build_substrate(regime: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Run a mismatched-radii D'Orsogna regime and return particle-history."""
    from epc.models.dorsogna_spp import DOrsognaSPPModel

    p = regime["params"]
    model = DOrsognaSPPModel(
        n_particles=p["n_particles"],
        C_a=p["C_a"],
        C_r=p["C_r"],
        l_a=p["l_a"],
        l_r=p["l_r"],
        alpha=p["alpha"],
        beta=p["beta"],
        dt=p["dt"],
        init_mode=p["init_mode"],
        init_radius=p["init_radius"],
        seed=regime["seed"],
    )
    return model.run(n_steps=p["n_steps"])
