"""P3 (Turing wavelength) — Class C adversarial negatives.

REBUILT (validation-rebuild, 2026-06-14). Previously this module declared
Class C "N/A", on the argument that non-Turing Gray-Scott regimes simply
produce homogeneous fields rejected by the field_std prerequisite. That left
the detector's hardest failure mode untested: the original peak-to-mean /
peak_k metric answers only "is there a dominant spatial wavelength?", which is
satisfied by ANY periodic field — an imposed sinusoid, oriented stripes, or a
travelling / rotating wave — none of which are Turing patterns. With those
adversaries removed, a static imposed sinusoid reached DEFINITIVE (a false
positive), and the green panel was dishonest.

This Class C restores genuine near-pattern adversaries: fields that DO have a
strong characteristic wavelength (they clear peak-to-mean and the field_std /
n_unique prerequisites) but are NOT Turing patterns. The rebuilt detector must
reject each via its Turing-specificity gates:

  * imposed_static_sinusoid / imposed_static_stripes_diag — stationary but
    ANISOTROPIC (a single wave vector): rejected by the angular-isotropy gate.
  * rotating_spiral / target_waves — (near-)isotropic ring but NON-STATIONARY
    (the pattern advects): rejected by the stationarity gate.
  * travelling_plane_wave — anisotropic AND non-stationary: rejected by both.

A genuine Turing pattern is a STATIONARY, ISOTROPIC standing wave selected by a
diffusion-driven instability; these five impostors each violate exactly one (or
both) of those defining properties, so they are the discriminating test.
"""

from __future__ import annotations

from typing import Any, Dict, List

import numpy as np


N_GRID = 64
N_SNAP = 100
K_CYCLES = 5  # ~12.8 px wavelength on a 64-grid, matching the labyrinth scale


CONFIG: Dict[str, Any] = {
    "substrate_id": "P3_non_turing_periodic",
    "format": "field",
    "description": (
        "5 non-Turing periodic-field adversaries on a 64x64 grid (100 "
        "snapshots each), all with a strong ~12.8 px characteristic "
        "wavelength so they clear peak-to-mean and the field_std / n_unique "
        "prerequisites: imposed_static_sinusoid and imposed_static_stripes_diag "
        "(stationary but anisotropic single wavevector -> isotropy gate), "
        "rotating_spiral and target_waves (isotropic ring but non-stationary "
        "-> stationarity gate), and travelling_plane_wave (both). A true Turing "
        "pattern is a stationary ISOTROPIC standing wave; each impostor "
        "violates one or both properties."
    ),
    "regimes": [
        {"label": "imposed_static_sinusoid",
         "params": {"kind": "static_sinusoid"}, "seed": 701},
        {"label": "imposed_static_stripes_diag",
         "params": {"kind": "static_stripes_diag"}, "seed": 702},
        {"label": "rotating_spiral",
         "params": {"kind": "rotating_spiral"}, "seed": 703},
        {"label": "target_waves",
         "params": {"kind": "target_waves"}, "seed": 704},
        {"label": "travelling_plane_wave",
         "params": {"kind": "travelling_plane_wave"}, "seed": 705},
    ],
}


def build_substrate(regime: Dict[str, Any]) -> List[Dict[str, Any]]:
    """Generate a non-Turing periodic field history for one regime."""
    kind = regime["params"]["kind"]
    rng = np.random.default_rng(regime.get("seed", 0))
    N = N_GRID
    x = np.arange(N)
    X, Y = np.meshgrid(x, x)
    cx = cy = N / 2.0
    R = np.sqrt((X - cx) ** 2 + (Y - cy) ** 2)
    TH = np.arctan2(Y - cy, X - cx)
    k = 2.0 * np.pi * K_CYCLES / N

    history: List[Dict[str, Any]] = []
    for t in range(N_SNAP):
        if kind == "static_sinusoid":
            # Stationary single-direction sinusoid + small per-frame noise.
            field = np.sin(k * X) + 0.01 * rng.standard_normal((N, N))
        elif kind == "static_stripes_diag":
            # Stationary diagonal stripes (single oblique wavevector).
            field = np.sin(k * (X + Y) / np.sqrt(2.0)) \
                + 0.01 * rng.standard_normal((N, N))
        elif kind == "rotating_spiral":
            # Single-arm spiral rotating in time (isotropic ring, non-stationary).
            field = np.sin(TH + k * R - 0.3 * t)
        elif kind == "target_waves":
            # Concentric waves expanding in time (isotropic, non-stationary).
            field = np.sin(k * R - 0.3 * t)
        elif kind == "travelling_plane_wave":
            # Plane wave advecting in x (anisotropic AND non-stationary).
            field = np.sin(k * (X - 0.3 * t))
        else:
            raise ValueError(f"unknown P3 Class C kind: {kind!r}")
        history.append({"field": field.astype(np.float32), "step": t})
    return history
