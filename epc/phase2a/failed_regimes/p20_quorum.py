"""P20 (quorum sensing) — Class C failed regimes.

Two canonical failure modes per the Sprint 84 brief:
  1. sub_threshold: Agent density never crosses the activation threshold.
     System stays permanently OFF. The detector should reject at screening
     (no OFF→ON transition: f_low ≈ f_high ≈ 0).
  2. graded_response: System response increases smoothly/linearly with
     density — no sharp switch, no hysteresis. This is the key P18
     (consensus) discrimination case. The detector should reject at
     confirmation (step-function R² too low, transition too gradual).
"""

from __future__ import annotations

from typing import Any, Dict, List

import numpy as np


CONFIG: Dict[str, Any] = {
    "substrate_id": "P20_quorum_sensing_class_c",
    "format": "density_sweep",
    "status": "populated",
    "regimes": [
        {
            "label": "sub_threshold",
            "description": (
                "Density never reaches the activation threshold. System "
                "remains permanently OFF (fraction_on ≈ 0 at all densities)."
            ),
            "params": {
                "n_density_levels": 40,
                "n_steps_per_density": 100,
                "density_min": 0.01,
                "density_max": 0.5,  # well below activation threshold
                "activation_threshold": 5.0,  # unreachably high
            },
            "seed": 42,
        },
        {
            "label": "graded_response",
            "description": (
                "Smooth sigmoid response to density with LARGE amplitude "
                "(fraction_on spans ~0.05→0.95, jump>0.5 so it passes the "
                "amplitude screen) but a GENTLE slope — the transition spans "
                "most of the density range (relative_width ~0.74). No sharp "
                "switch, no hysteresis. The detector must reject this at "
                "screening BY SHARPNESS (relative_width >= 0.4), not by "
                "amplitude. This is the key P18 (consensus) discrimination "
                "case and exercises the sharpness gate as a live discriminator."
            ),
            "params": {
                "n_density_levels": 40,
                "n_steps_per_density": 100,
                "density_min": 0.1,
                "density_max": 6.0,
                "sigmoid_center": 3.0,
                "sigmoid_slope": 1.0,  # gentle slope, wide span → graded, large amplitude
                "noise_std": 0.03,
            },
            "seed": 42,
        },
    ],
}


def build_substrate(regime: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Build a P20 failed-regime substrate."""
    p = regime["params"]
    seed = regime["seed"]
    rng = np.random.default_rng(seed)
    n_levels = p["n_density_levels"]
    n_steps = p["n_steps_per_density"]
    d_min = p["density_min"]
    d_max = p["density_max"]

    label = regime["label"]
    history: List[Dict[str, Any]] = []

    densities_up = np.linspace(d_min, d_max, n_levels)
    densities_down = np.linspace(d_max, d_min, n_levels)
    all_densities = np.concatenate([densities_up, densities_down])

    if label == "sub_threshold":
        # System stays permanently OFF at all density levels
        activation_threshold = p["activation_threshold"]
        for d_idx, density in enumerate(all_densities):
            direction = 'up' if d_idx < n_levels else 'down'
            for step_i in range(n_steps):
                # Fraction_on stays near 0 (small noise)
                fraction_on = float(np.clip(
                    rng.normal(0.02, 0.01), 0.0, 0.1
                ))
                concentration = density * fraction_on
                history.append({
                    'density': float(density),
                    'concentration': float(concentration),
                    'collective_state': 0,
                    'fraction_on': fraction_on,
                    'step': step_i,
                    'density_idx': d_idx,
                    'sweep_direction': direction,
                })

    elif label == "graded_response":
        # Smooth sigmoid response: fraction_on = sigmoid(density)
        # No hysteresis — same curve up and down
        center = p["sigmoid_center"]
        slope = p["sigmoid_slope"]
        noise_std = p["noise_std"]
        for d_idx, density in enumerate(all_densities):
            direction = 'up' if d_idx < n_levels else 'down'
            # Gentle logistic curve
            base_frac = 1.0 / (1.0 + np.exp(-slope * (density - center)))
            for step_i in range(n_steps):
                fraction_on = float(np.clip(
                    base_frac + rng.normal(0.0, noise_std), 0.0, 1.0
                ))
                concentration = density * fraction_on
                collective = 1 if fraction_on > 0.5 else 0
                history.append({
                    'density': float(density),
                    'concentration': float(concentration),
                    'collective_state': collective,
                    'fraction_on': fraction_on,
                    'step': step_i,
                    'density_idx': d_idx,
                    'sweep_direction': direction,
                })

    return history
