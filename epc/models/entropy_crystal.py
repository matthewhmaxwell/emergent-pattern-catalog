"""Entropy-driven crystallization model for P35 — the none_entropy substrate.

Repulsive disks (excluded volume only, NO attraction) under slow compression with
annealed noise self-organize hexatic order from a disordered fluid (the Alder
transition). Frames carry {"positions": Nx2}; metadata carries box_size. The neg_*
generators are the Phase-2a lookalikes.
"""
from __future__ import annotations

from typing import Any, Dict, List, Tuple

import numpy as np

Built = Tuple[List[Dict[str, Any]], Dict[str, Any]]


def _step_forces(pos: np.ndarray, d: float, L: float) -> np.ndarray:
    rij = pos[:, None, :] - pos[None, :, :]
    rij -= L * np.round(rij / L)
    dist = np.sqrt((rij ** 2).sum(-1))
    np.fill_diagonal(dist, 1e9)
    overlap = np.maximum(d - dist, 0.0)
    fhat = rij / dist[..., None]
    return (overlap[..., None] * fhat).sum(1)


def hard_disk_crystallization(seed: int = 0, N: int = 196, L: float = 14.0,
                              eta_end: float = 0.74, n_steps: int = 2500,
                              n_frames: int = 80, dt: float = 0.05,
                              noise0: float = 0.16) -> Built:
    """Entropy/packing-driven crystallization (the none_entropy positive)."""
    rng = np.random.default_rng(seed)
    pos = rng.uniform(0, L, size=(N, 2))
    d_end = 2.0 * np.sqrt(eta_end * L * L / (N * np.pi))
    d0 = 0.30 * d_end
    frames: List[Dict[str, Any]] = []
    snap = max(1, n_steps // n_frames)
    for step in range(n_steps):
        frac = step / n_steps
        d = d0 + (d_end - d0) * frac
        F = _step_forces(pos, d, L)
        pos = (pos + dt * F + rng.normal(0, noise0 * (1.0 - frac), size=(N, 2))) % L
        if step % snap == 0:
            frames.append({"positions": pos.copy()})
    frames.append({"positions": pos.copy()})
    return frames, {"model": "hard_disk_crystallization", "eta_end": eta_end,
                    "N": N, "box_size": L}


def neg_dilute_gas(seed: int = 0) -> Built:
    """Same dynamics at LOW density -> stays a fluid (no hexatic order)."""
    return hard_disk_crystallization(seed, N=196, L=14.0, eta_end=0.20,
                                     n_steps=2500, n_frames=80)


def neg_static_crystal(seed: int = 0) -> Built:
    """A pre-formed crystal held FIXED: high psi6 but ZERO gain (no emergence)."""
    h, _ = hard_disk_crystallization(seed, N=196, L=14.0, eta_end=0.74,
                                     n_steps=2500, n_frames=80, noise0=0.16)
    last = h[-1]["positions"]
    return ([{"positions": last.copy()} for _ in range(60)],
            {"model": "static_crystal", "box_size": 14.0})
