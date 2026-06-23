"""Ring 1 — entropy-driven ordering gap (interact = none_entropy, empty in catalog).

Ordering with NO agent-agent interaction FORCE — only excluded volume / packing
entropy. The Alder transition (Alder & Wainwright 1957): purely repulsive disks
spontaneously CRYSTALLIZE into a hexagonal lattice at high density, driven by entropy
(free-volume maximization), not by any attraction. Distinct from P1 (type
segregation), P2 (active motility): identical, passive, repulsion-only.

Minimal model: soft repulsive disks under slow compression (growing diameter) with
annealing noise; from a random fluid they self-organize into a hexatic/crystalline
arrangement. Frames carry {"positions": Nx2}. Null: held at low density (stays fluid).
"""
from __future__ import annotations

from typing import Any, Dict, List, Tuple

import numpy as np

Built = Tuple[List[Dict[str, Any]], Dict[str, Any]]


def _step_forces(pos: np.ndarray, d: float, L: float) -> np.ndarray:
    rij = pos[:, None, :] - pos[None, :, :]
    rij -= L * np.round(rij / L)                       # periodic min-image
    dist = np.sqrt((rij ** 2).sum(-1))
    np.fill_diagonal(dist, 1e9)
    overlap = np.maximum(d - dist, 0.0)                # >0 only when overlapping
    fhat = rij / dist[..., None]
    return (overlap[..., None] * fhat).sum(1)          # harmonic push-apart


def hard_disk_crystallization(seed: int = 0, N: int = 196, L: float = 14.0,
                              eta_end: float = 0.80, n_steps: int = 1600,
                              n_frames: int = 80, dt: float = 0.05,
                              noise0: float = 0.06) -> Built:
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
        noise = rng.normal(0, noise0 * (1.0 - frac), size=(N, 2))   # annealed
        pos = (pos + dt * F + noise) % L
        if step % snap == 0:
            frames.append({"positions": pos.copy()})
    frames.append({"positions": pos.copy()})
    return frames, {"model": "hard_disk_crystallization", "eta_end": eta_end,
                    "N": N, "box_size": L}


def null_dilute_gas(seed: int = 0, N: int = 196, L: float = 14.0,
                    eta_end: float = 0.20, n_steps: int = 1600, n_frames: int = 80) -> Built:
    """Same dynamics held at LOW density: stays a disordered fluid (no crystal).
    none_entropy NULL -> must read NO-EMERGENCE."""
    return hard_disk_crystallization(seed, N=N, L=L, eta_end=eta_end,
                                     n_steps=n_steps, n_frames=n_frames)


def local_psi6(pos: np.ndarray, L=None, n_neighbors: int = 6) -> float:
    """Mean per-particle |psi6_i| — local hexatic bond-orientational order over each
    particle's 6 nearest neighbours. ~0.9 single crystal, ~0.5 polycrystal/RCP, ~0.3
    fluid/random (finite-N background). L=None -> no periodic wrap (raw distances)."""
    pos = np.asarray(pos, dtype=float)
    N = pos.shape[0]
    if N < n_neighbors + 1:
        return 0.0
    rij = pos[:, None, :] - pos[None, :, :]
    if L:
        rij -= L * np.round(rij / L)
    dist = np.sqrt((rij ** 2).sum(-1))
    np.fill_diagonal(dist, 1e9)
    vals = np.empty(N)
    for i in range(N):
        nn = np.argsort(dist[i])[:n_neighbors]
        ang = np.arctan2(rij[i, nn, 1], rij[i, nn, 0])
        vals[i] = abs(np.mean(np.exp(6j * ang)))
    return float(np.mean(vals))


def psi6(pos: np.ndarray, L: float, n_neighbors: int = 6) -> float:
    """Global hexatic bond-orientational order |<(1/6) sum exp(6 i theta_ij)>| over
    each particle's 6 nearest neighbours. ~1 for a triangular crystal; ~0 for a fluid."""
    N = pos.shape[0]
    rij = pos[:, None, :] - pos[None, :, :]
    rij -= L * np.round(rij / L)
    dist = np.sqrt((rij ** 2).sum(-1))
    np.fill_diagonal(dist, 1e9)
    psi = np.zeros(N, dtype=complex)
    for i in range(N):
        nn = np.argsort(dist[i])[:n_neighbors]
        ang = np.arctan2(rij[i, nn, 1], rij[i, nn, 0])
        psi[i] = np.mean(np.exp(6j * ang))
    return float(abs(np.mean(psi)))
