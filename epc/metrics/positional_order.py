"""Positional-order metrics for P35 (entropy-driven crystallization).

Local hexatic bond-orientational order psi6 on a 2-D point set.
"""
from __future__ import annotations

import numpy as np


def local_psi6(pos: np.ndarray, L=None, n_neighbors: int = 6) -> float:
    """Mean per-particle |psi6_i| — local hexatic bond-orientational order over each
    particle's 6 nearest neighbours. ~0.9 single crystal, ~0.5 polycrystal/RCP,
    ~0.3 fluid/random (finite-N background). L=None -> raw distances (no periodic wrap)."""
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
