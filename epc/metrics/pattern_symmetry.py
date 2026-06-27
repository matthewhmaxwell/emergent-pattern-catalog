"""Pattern-symmetry / lattice-class lens — the AZIMUTHAL structure of the power spectrum
(a documented missing axis). structure_factor reports only the principal |k| (the scale);
stripes, squares, hexagons and a labyrinth can all share that |k| yet differ in symmetry.

At the dominant radial ring |k|=k*, take the azimuthal power profile a(phi) and its angular
Fourier harmonics. The DOMINANT harmonic m is the lattice symmetry: m=2 stripes/lamellar,
m=4 square, m=6 hexagonal. Its relative amplitude (symmetry_strength) separates a discrete
lattice (one harmonic dominates) from an isotropic labyrinth (flat ring, power spread over all
harmonics). ring_contrast (radial peak / median) gates out a structureless gas (no ring at all).

Discriminator = the angular HARMONIC ORDER + its concentration, NOT the radial |k| (which the
structure-factor lens already owns) and NOT total spectral power.

Scalars (late-window mean):
  lattice_symmetry   — dominant angular harmonic m (2/4/6 = stripes/square/hex).
  symmetry_strength  — amplitude of harmonic m / total azimuthal AC power (~1 pure lattice, ~0 labyrinth).
  ring_contrast      — radial-peak / radial-median (>~2 a real ring exists; ~1 = gas).
Reads a 2D field ('field'/'u'/'theta_field'/'concentration') or a density grid from 'positions'.
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional

import numpy as np


def _field_2d(frame: Dict[str, Any], G: int = 96) -> Optional[np.ndarray]:
    if not isinstance(frame, dict):
        return None
    for k in ("field", "u", "concentration", "theta_field"):
        if k in frame:
            a = np.asarray(frame[k], dtype=float)
            if a.ndim == 2 and min(a.shape) >= 12:
                return np.cos(a) if k == "theta_field" else a   # director -> a smooth field
    pos = frame.get("positions")
    if pos is not None:
        pos = np.asarray(pos, dtype=float)
        if pos.ndim == 2 and pos.shape[1] == 2 and len(pos) >= 16:
            lo = pos.min(0); span = np.ptp(pos, axis=0); span[span == 0] = 1.0
            ij = np.clip(((pos - lo) / span * (G - 1)).astype(int), 0, G - 1)
            grid = np.zeros((G, G))
            np.add.at(grid, (ij[:, 1], ij[:, 0]), 1.0)
            return grid
    return None


def _symmetry(field: np.ndarray) -> Optional[Dict[str, float]]:
    f = field - field.mean()
    if not np.any(f):
        return None
    P = np.abs(np.fft.fftshift(np.fft.fft2(f))) ** 2
    ny, nx = P.shape; cy, cx = ny // 2, nx // 2
    Y, X = np.indices(P.shape)
    R = np.hypot(X - cx, Y - cy); TH = np.arctan2(Y - cy, X - cx)
    rmax = min(cy, cx) - 1
    if rmax < 4:
        return None
    rbin = R.astype(int)
    radial = np.array([P[rbin == rr].mean() for rr in range(2, rmax)])
    if radial.size < 3 or not np.any(radial):
        return None
    kstar = int(np.argmax(radial)) + 2
    ring_contrast = float(radial.max() / (np.median(radial) + 1e-12))

    msk = (R >= kstar - 1.0) & (R <= kstar + 1.0)
    th = TH[msk] % (2 * np.pi); pw = P[msk]
    nb = 180
    a = np.zeros(nb); cnt = np.zeros(nb)
    ii = (th / (2 * np.pi) * nb).astype(int) % nb
    np.add.at(a, ii, pw); np.add.at(cnt, ii, 1)
    a = a / np.maximum(cnt, 1)                                # raw ring azimuthal profile (>=0)
    mean_a = float(a.mean())
    # CONCENTRATION = how spiky the ring is: discrete lattice (sharp peaks) high, labyrinth (flat) ~0
    concentration = float(a.std() / (mean_a + 1e-12))
    A = np.abs(np.fft.rfft(a - mean_a))
    if A.size < 9 or A[1:].sum() < 1e-12:
        return {"lattice_symmetry": 0, "azimuthal_concentration": 0.0, "ring_contrast": round(ring_contrast, 2)}
    m = int(np.argmax(A[1:9]) + 1)
    return {"lattice_symmetry": m, "azimuthal_concentration": round(concentration, 3),
            "ring_contrast": round(ring_contrast, 2)}


def pattern_symmetry(history: List[Dict[str, Any]]) -> Optional[Dict[str, float]]:
    frames = [f for f in history if isinstance(f, dict)]
    late = frames[len(frames) // 2:] or frames
    syms, conc, rings = [], [], []
    for f in late:
        fld = _field_2d(f)
        if fld is None:
            continue
        r = _symmetry(fld)
        if r is None:
            continue
        syms.append(r["lattice_symmetry"]); conc.append(r["azimuthal_concentration"]); rings.append(r["ring_contrast"])
    if not conc:
        return None
    # report a lattice symmetry only where a ring is present AND the ring is spiky (a real lattice)
    voted = [s for s, rc, cc in zip(syms, rings, conc) if rc >= 2.0 and cc >= 3.5]
    sym = int(np.bincount(voted).argmax()) if voted else 0
    return {"lattice_symmetry": sym,
            "azimuthal_concentration": round(float(np.mean(conc)), 3),
            "ring_contrast": round(float(np.mean(rings)), 2)}
