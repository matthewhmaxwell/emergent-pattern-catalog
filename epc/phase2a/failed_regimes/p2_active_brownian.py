"""P2 (Active Brownian Particles) failed regimes: 10 low-Péclet regimes.

Per Sprint 46 brief: 10 low-activity / high-noise regimes that suppress MIPS.

Low Pe = v0 / D_r (rotational correlation time × speed). Below Pe ≈ 50 the
density-dependent slowdown does not produce macroscopic phase separation.

Pe values: linspace(0.5, 10.0, 10) — all well below the MIPS onset (Pe ≈ 50).
"""

from __future__ import annotations

from typing import Any, Dict, List

import numpy as np


PE_VALUES = list(np.linspace(0.5, 10.0, 10))

# Fixed base parameters (same phi as canonical positive).
_V0 = 1.0
_N = 400
_PHI = 0.5
_RHO_STAR = 4.0
_L = float(np.sqrt(_N * np.pi / 4.0 / _PHI))   # box_size ≈ 25.1

CONFIG: Dict[str, Any] = {
    "substrate_id": "P2_abp_low_peclet",
    "format": "particles",
    "description": (
        f"10 ABP regimes at Pe ∈ linspace(0.5, 10.0, 10), "
        f"N={_N}, phi={_PHI}, v0={_V0}. "
        "Below MIPS onset (Pe≈50): no phase separation, two_phase_score≈0."
    ),
    "regimes": [
        {
            "label": f"Pe={pe:.2f}",
            "params": {
                "n_particles": _N,
                "box_size": _L,
                "v0": _V0,
                "D_r": float(_V0 / pe),
                "rho_star": _RHO_STAR,
                "r_cg": 1.0,
                "dt": 0.05,
                "n_steps": 600,
            },
            "seed": 300 + i,
        }
        for i, pe in enumerate(PE_VALUES)
    ],
}


def build_substrate(regime: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Run a low-Pe ABP regime and return particle-history."""
    from epc.models.active_brownian_particles import ActiveBrownianParticles

    p = regime["params"]
    model = ActiveBrownianParticles(
        n_particles=p["n_particles"],
        box_size=p["box_size"],
        v0=p["v0"],
        D_r=p["D_r"],
        rho_star=p["rho_star"],
        r_cg=p["r_cg"],
        dt=p["dt"],
        seed=regime["seed"],
    )
    return model.run(n_steps=p["n_steps"])
