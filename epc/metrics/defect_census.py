"""Topological-defect census lens — winding-number singularities in a continuous
orientation/phase field (a documented missing axis: defect-driven order).

Global order parameters AVERAGE defects out: a balanced +-1/2 defect sea (active-nematic
turbulence, the BKT regime) has ~zero net order yet is strongly organized. Vorticity sees
advected angular momentum, not field singularities. This lens counts the SIGNED winding of
the angle field around each plaquette -> defect density + net topological charge.

Denoise by COARSE-GRAINING: smooth cos(m*theta)/sin(m*theta) before counting. Real
topological defects are protected and survive smoothing; noise spawns +/- pairs that
annihilate. m = 2 for a nematic (head-tail, pi-periodic 'theta_field' -> +-1/2 defects),
m = 1 for a polar phase field ('phases'/'theta', 2pi-periodic -> +-1 defects).

Scalars (late-window mean):
  defect_density — fraction of plaquettes carrying a (post-smoothing) defect. THE
                   discriminator: a defect SEA scores high; uniform/aligned ~0; a single
                   vortex ~tiny; iid noise ~0 (annihilates under smoothing).
  net_charge     — summed signed winding (≈0 for a balanced sea; nonzero if charge imposed).
Reads a 2D 'theta_field' (nematic) or 'phases'/'theta' (polar). Returns None otherwise.
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional, Tuple

import numpy as np


def _angle_field(frame: Dict[str, Any]) -> Tuple[Optional[np.ndarray], int]:
    if isinstance(frame, dict) and "theta_field" in frame:
        a = np.asarray(frame["theta_field"], dtype=float)
        return (a, 2) if a.ndim == 2 else (None, 0)          # nematic
    for k in ("phases", "theta"):
        if isinstance(frame, dict) and k in frame:
            a = np.asarray(frame[k], dtype=float)
            return (a, 1) if a.ndim == 2 else (None, 0)      # polar
    return None, 0


def _defect_stats(theta: np.ndarray, m: int, smooth: float = 1.5,
                  min_charge: float = 0.4, coh_thr: float = 0.35) -> Tuple[float, float, float]:
    from scipy.ndimage import gaussian_filter
    c = gaussian_filter(np.cos(m * theta), smooth, mode="wrap")
    s = gaussian_filter(np.sin(m * theta), smooth, mode="wrap")
    coh = float(np.mean(np.sqrt(c ** 2 + s ** 2)))           # mean coarse-grained coherence
    ang = np.arctan2(s, c)                                    # = smoothed m*theta in (-pi,pi]

    def wrap(d):
        return (d + np.pi) % (2 * np.pi) - np.pi

    a00 = ang[:-1, :-1]; a01 = ang[:-1, 1:]; a11 = ang[1:, 1:]; a10 = ang[1:, :-1]
    loop = wrap(a01 - a00) + wrap(a11 - a01) + wrap(a10 - a11) + wrap(a00 - a10)
    charge = loop / (2 * np.pi) / m                           # winding of theta (nematic -> +-1/2)
    n_def = int((np.abs(charge) > min_charge).sum())
    # COHERENCE GATE: in an incoherent (noise) field, windings are spurious -> not defects.
    # A real defect is a singularity in an otherwise locally-ordered field.
    density = (n_def / charge.size) if coh >= coh_thr else 0.0
    return density, float(charge.sum()), coh


def defect_census(history: List[Dict[str, Any]]) -> Optional[Dict[str, float]]:
    frames = [f for f in history if isinstance(f, dict)]
    late = frames[len(frames) // 2:] or frames
    dens, nets, cohs = [], [], []
    for f in late:
        th, m = _angle_field(f)
        if th is None or th.ndim != 2 or min(th.shape) < 6:
            continue
        d, net, coh = _defect_stats(th, m)
        dens.append(d); nets.append(net); cohs.append(coh)
    if not dens:
        return None
    return {"defect_density": round(float(np.mean(dens)), 5),
            "net_charge": round(float(np.mean(nets)), 3),
            "coherence": round(float(np.mean(cohs)), 3)}
