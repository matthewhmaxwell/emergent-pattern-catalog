"""P25 (equifinality) — Class C failed regimes.

Two canonical failure modes per the Sprint 82 brief:
  1. narrow_basin: Only ICs very close to the target converge. Basin
     volume < 0.8 (most of IC space does NOT reach the target).
     The detector should reject at confirmation (basin_volume < 0.8).
  2. divergent_dynamics: ICs actively diverge from the target (repulsive
     dynamics). Convergence ratio >> 0.1. The detector should reject at
     screening.
"""

from __future__ import annotations

from typing import Any, Dict, List

import numpy as np


CONFIG: Dict[str, Any] = {
    "substrate_id": "P25_equifinality_class_c",
    "format": "canalization_bundle",
    "status": "populated",
    "regimes": [
        {
            "label": "narrow_basin",
            "description": (
                "Only ICs within 0.5 of the target converge; ICs beyond "
                "that radius diverge. Basin volume << 0.8."
            ),
            "params": {
                "n_dims": 10,
                "n_ics": 20,
                "n_steps": 200,
                "basin_strength": 2.0,
                "ic_spread": 0.3,       # narrow IC range (near target only)
                "diverge_radius": 0.5,   # ICs beyond this diverge
            },
            "seed": 42,
        },
        {
            "label": "divergent_dynamics",
            "description": (
                "Repulsive dynamics: all ICs move away from the target. "
                "Convergence ratio >> 0.1."
            ),
            "params": {
                "n_dims": 10,
                "n_ics": 20,
                "n_steps": 200,
                "repulsion_strength": 1.0,
            },
            "seed": 42,
        },
    ],
}


def build_substrate(regime: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Build a P25 failed-regime substrate."""
    p = regime["params"]
    seed = regime["seed"]
    rng = np.random.default_rng(seed)
    n_dims = p["n_dims"]
    n_ics = p["n_ics"]
    n_steps = p["n_steps"]
    target = rng.standard_normal(n_dims)

    label = regime["label"]
    history: List[Dict[str, Any]] = []

    if label == "narrow_basin":
        # Near-target ICs converge; far ICs don't. But ic_spread is so
        # narrow that ALL ICs start near the target — giving the appearance
        # of convergence. However, the *basin volume* metric tests whether
        # a WIDE range of ICs converge, so we generate ICs from a wide
        # spread and only let near ones converge.
        wide_spread = 5.0
        diverge_r = p["diverge_radius"]
        basin_strength = p["basin_strength"]
        dt = 0.05

        for trial in range(n_ics):
            direction = rng.standard_normal(n_dims)
            direction /= np.linalg.norm(direction) + 1e-12
            radius = wide_spread * rng.random() ** (1.0 / n_dims)
            x = target + direction * radius
            ic = x.copy()

            for step in range(n_steps):
                dist = float(np.linalg.norm(x - target))
                converged = dist < 0.1

                history.append({
                    'state': x.copy(),
                    'step': step,
                    'trial': trial,
                    'ic': ic.copy(),
                    'target': target.copy(),
                    'distance_to_target': dist,
                    'converged': converged,
                })

                # Near ICs converge; far ICs diffuse
                if dist < diverge_r:
                    force = -basin_strength * (x - target)
                    x = x + force * dt
                else:
                    x = x + rng.standard_normal(n_dims) * 0.1

    elif label == "divergent_dynamics":
        repulsion = p["repulsion_strength"]
        dt = 0.05

        for trial in range(n_ics):
            x = target + rng.standard_normal(n_dims) * 2.0
            ic = x.copy()

            for step in range(n_steps):
                dist = float(np.linalg.norm(x - target))
                history.append({
                    'state': x.copy(),
                    'step': step,
                    'trial': trial,
                    'ic': ic.copy(),
                    'target': target.copy(),
                    'distance_to_target': dist,
                    'converged': dist < 0.1,
                })

                # Repulsive dynamics: push away from target
                delta = x - target
                norm = np.linalg.norm(delta) + 1e-12
                x = x + repulsion * delta / norm * dt + rng.standard_normal(n_dims) * 0.05

    return history
