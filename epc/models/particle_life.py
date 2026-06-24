"""Particle-life ("Clusters") substrate for the Ring-2 sourcing pilot.

K particle types drift in a periodic box; each ordered pair of types (a,b) has a
signed interaction strength F[a,b] (attraction>0 / repulsion<0) acting in a mid-range
band, plus a universal short-range repulsion that prevents overlap. The PARAMETER SET
varied by the search is the K x K matrix F (and K). Same fixed rules, different F ->
drifting gases, clusters, chasers, oscillating cells, membranes, etc.

Emits frames with positions / velocities / types, the observables the instrument's
morphology + orientation lenses read.

Refs: Ventrella "Clusters"; the widely-reproduced "particle life" force law.
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional, Tuple

import numpy as np

Built = Tuple[List[Dict[str, Any]], Dict[str, Any]]


def particle_life(seed: int = 0, K: int = 4, F: Optional[np.ndarray] = None,
                  N: int = 250, L: float = 10.0, r_max: float = 2.0,
                  r_min_frac: float = 0.3, friction: float = 0.85, dt: float = 0.1,
                  steps: int = 300, n_frames: int = 36, v_max: float = 3.0) -> Built:
    rng = np.random.default_rng(seed)
    if F is None:
        F = rng.uniform(-1.0, 1.0, size=(K, K))
    F = np.asarray(F, dtype=float)
    K = F.shape[0]
    types = rng.integers(0, K, size=N)
    pos = rng.uniform(0, L, size=(N, 2))
    vel = np.zeros((N, 2))
    r_min = r_min_frac * r_max
    snap = max(1, steps // n_frames)
    frames: List[Dict[str, Any]] = [{"positions": pos.copy(), "velocities": vel.copy(),
                                     "types": types.copy(), "step": 0}]
    Fab = F[types[:, None], types[None, :]]                 # (N,N) coeff per ordered pair

    for s in range(1, steps + 1):
        d = pos[None, :, :] - pos[:, None, :]              # d[i,j] = pos_j - pos_i (i->j)
        d -= L * np.round(d / L)                            # periodic minimum image
        r = np.sqrt((d ** 2).sum(-1))
        np.fill_diagonal(r, np.inf)
        f = np.zeros_like(r)
        close = r < r_min
        f[close] = -(1.0 - r[close] / r_min)               # universal short-range repulsion
        band = (r >= r_min) & (r < r_max)
        tent = 1.0 - np.abs(2.0 * r - r_min - r_max) / (r_max - r_min)   # peak mid-band
        f[band] = Fab[band] * tent[band]                   # F>0 attract, F<0 repel
        inv = np.where(r > 0, 1.0 / r, 0.0)
        ax = (f * d[:, :, 0] * inv).sum(1)
        ay = (f * d[:, :, 1] * inv).sum(1)
        vel = vel * friction + np.column_stack([ax, ay]) * dt
        np.clip(vel, -v_max, v_max, out=vel)
        pos = (pos + vel * dt) % L
        if s % snap == 0:
            frames.append({"positions": pos.copy(), "velocities": vel.copy(),
                           "types": types.copy(), "step": s})
    return frames, {"model": "particle_life", "K": int(K), "N": N, "L": L,
                    "r_max": r_max, "F": F.round(4).tolist()}
