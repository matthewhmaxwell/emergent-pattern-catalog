"""P35 (entropy-driven crystallization) -- Class C failed regimes.

The target's own model, ``epc.models.entropy_crystal.hard_disk_crystallization``
(N=196 soft disks, L=14, compression with annealed noise over 2500 steps, the
positive's protocol), compressed only to packing fractions BELOW the freezing
range. For hard disks the liquid is stable up to eta ~ 0.700, with liquid-hexatic
coexistence to 0.716 and the solid above ~0.72 (Bernard & Krauth 2011, PRL 107
155704; the Alder transition). The positive ends at eta = 0.74.

  dilute_gas eta_end=0.20   standalone negative neg_dilute_gas, reused (same call)
  fluid eta_end=0.40, 0.55, 0.65   added; all below the literature freezing range

Ground truth from the model (analysis/validation_rebuild/p33_p37_class_c_groundtruth.py):
crystalline-particle fraction (local |psi6| > 0.7) of the final configuration is
0.02-0.05 at eta_end 0.20-0.55 and 0.13 at 0.65.

Two seeds per regime (200-207), distinct from the positive seeds (0-4). Frames
carry no box_size; the panel scores them with the canonical metadata box (L=14),
which is their own box.
"""
from __future__ import annotations

from typing import Any, Dict, List

_HD = {"N": 196, "L": 14.0, "n_steps": 2500, "n_frames": 80, "noise0": 0.16}

CONFIG: Dict[str, Any] = {
    "substrate_id": "P35_entropy_crystal_class_c",
    "format": "positions",
    "status": "populated",
    "regimes": [
        {"label": f"{'dilute_gas' if e == 0.20 else 'fluid'}_eta_end={e:.2f}",
         "params": dict(_HD, eta_end=e), "seed": 200 + 2 * i + k}
        for i, e in enumerate((0.20, 0.40, 0.55, 0.65)) for k in (0, 1)
    ],
}


def build_substrate(regime: Dict[str, Any]) -> List[Dict[str, Any]]:
    from epc.models.entropy_crystal import hard_disk_crystallization
    p = regime["params"]
    return hard_disk_crystallization(regime["seed"], N=p["N"], L=p["L"], eta_end=p["eta_end"],
                                     n_steps=p["n_steps"], n_frames=p["n_frames"],
                                     noise0=p["noise0"])[0]
