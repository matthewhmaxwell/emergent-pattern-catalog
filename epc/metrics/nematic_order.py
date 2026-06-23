"""Nematic-order metrics for P33 (active nematic with ±1/2 topological defects).

Substrate-faithful estimators on an orientation field:
  local_nematic_order        block-scale |<e^{2i theta}>|  (apolar local order)
  polar_order                global |<e^{i ang}>|          (polarity; ~0 for nematic)
  half_integer_defect_density density of ±1/2 topological defects (the nematic
                             fingerprint: charge ±0.5 plaquette winding) + integer
  angular_momentum           |mean r_hat x v_hat|          (milling discriminator)
  director_field             a frame's director field theta∈[0,pi): native
                             theta_field, else nematic-binned velocity headings.
"""
from __future__ import annotations

from typing import Any, Dict, Optional, Tuple

import numpy as np
from scipy.ndimage import uniform_filter

PI = np.pi


def wrap_nematic(d: np.ndarray) -> np.ndarray:
    """Wrap an angle difference into (-pi/2, pi/2] (nematic: theta ~ theta+pi)."""
    return d - PI * np.round(d / PI)


def local_nematic_order(theta: np.ndarray, block: int = 6) -> float:
    """Mean magnitude of the local (block-scale) nematic mean field |<e^{2i theta}>|.
    ~1 for locally aligned director, ~0 for isotropic."""
    c = np.cos(2.0 * theta); s = np.sin(2.0 * theta)
    cf = uniform_filter(c, size=block, mode="wrap")
    sf = uniform_filter(s, size=block, mode="wrap")
    return float(np.mean(np.sqrt(cf * cf + sf * sf)))


def polar_order(angles: np.ndarray) -> float:
    """Global polar order |<e^{i ang}>| of velocity headings. High=flock; ~0=apolar."""
    a = np.asarray(angles, dtype=float).ravel()
    return float(np.hypot(np.mean(np.cos(a)), np.mean(np.sin(a))))


def half_integer_defect_density(theta: np.ndarray) -> Tuple[float, float]:
    """Density of ±1/2 defects in a periodic director field theta∈[0,pi). Winding of
    the nematic angle around each plaquette; charge = winding/2pi; a ±1/2 defect has
    |charge|≈0.5. Returns (half_density, integer_density)."""
    er = wrap_nematic(np.roll(theta, -1, axis=0) - theta)
    eu = wrap_nematic(np.roll(theta, -1, axis=1) - theta)
    winding = (er + np.roll(eu, -1, axis=0) - np.roll(er, -1, axis=1) - eu)
    charge = winding / (2.0 * PI)
    half = np.abs(np.abs(charge) - 0.5) < 0.15
    integ = np.abs(np.abs(charge) - 1.0) < 0.15
    n = charge.size
    return float(half.sum() / n), float(integ.sum() / n)


def angular_momentum(positions: np.ndarray, headings: np.ndarray, box: float) -> float:
    """|mean (r_hat x v_hat)| about the (periodic) centre of mass. High = milling."""
    p = np.asarray(positions, dtype=float); h = np.asarray(headings, dtype=float)
    tx = 2 * PI * p[:, 0] / box; ty = 2 * PI * p[:, 1] / box
    cx = box / (2 * PI) * np.arctan2(np.mean(np.sin(tx)), np.mean(np.cos(tx))) % box
    cy = box / (2 * PI) * np.arctan2(np.mean(np.sin(ty)), np.mean(np.cos(ty))) % box
    dx = p[:, 0] - cx; dy = p[:, 1] - cy
    dx -= box * np.round(dx / box); dy -= box * np.round(dy / box)
    dist = np.maximum(np.hypot(dx, dy), 1e-12)
    cross = (dx / dist) * np.sin(h) - (dy / dist) * np.cos(h)
    return float(abs(np.mean(cross)))


def director_field(f: Dict[str, Any], G: int = 48) -> Optional[np.ndarray]:
    """Director field theta∈[0,pi) for a frame: native theta_field, else nematic-bin
    velocity headings onto a GxG grid (empty cells filled by nematic neighbour mean)."""
    if "theta_field" in f:
        return np.asarray(f["theta_field"], dtype=float)
    if "velocities" in f and "positions" in f:
        v = np.asarray(f["velocities"], dtype=float)
        p = np.asarray(f["positions"], dtype=float)
        box = float(f.get("box_size", p.max() + 1e-9))
        ang = np.arctan2(v[:, 1], v[:, 0])
        gx = np.clip((p[:, 0] / box * G).astype(int), 0, G - 1)
        gy = np.clip((p[:, 1] / box * G).astype(int), 0, G - 1)
        cs = np.zeros((G, G)); sn = np.zeros((G, G)); cnt = np.zeros((G, G))
        np.add.at(cs, (gx, gy), np.cos(2 * ang))
        np.add.at(sn, (gx, gy), np.sin(2 * ang))
        np.add.at(cnt, (gx, gy), 1.0)
        for _ in range(6):
            empty = cnt == 0
            if not empty.any():
                break
            csf = uniform_filter(cs, 3, mode="wrap"); snf = uniform_filter(sn, 3, mode="wrap")
            cs[empty] = csf[empty]; sn[empty] = snf[empty]; cnt[empty] = 1.0
        return (np.arctan2(sn, cs) / 2.0) % PI
    return None
