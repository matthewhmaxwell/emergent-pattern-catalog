"""Ring 0 — NEW out-of-catalog emergent systems (minimal canonical models).

Each is a textbook emergent phenomenon NOT among the 32 catalog patterns, in the
repo's observation-frame format (`field`/`grid` 2-D arrays, `positions` Nx2,
`phases`, `vector`). Distinctness from the nearest catalog pattern is noted in each
docstring so a false MATCH can be read as an instrument finding, not a bug.

References are in `ring0_systems.py` (the registry).
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional, Tuple

import numpy as np

Frame = Dict[str, Any]
History = List[Frame]
Built = Tuple[History, Optional[Dict[str, Any]]]


def _lap(f: np.ndarray) -> np.ndarray:
    """5-point periodic Laplacian."""
    return (np.roll(f, 1, 0) + np.roll(f, -1, 0) +
            np.roll(f, 1, 1) + np.roll(f, -1, 1) - 4.0 * f)


def cahn_hilliard_spinodal(seed: int = 0, G: int = 64, n_steps: int = 30000,
                           n_frames: int = 80, dt: float = 0.01,
                           kappa: float = 1.0, M: float = 1.0) -> Built:
    """Cahn-Hilliard spinodal decomposition: a CONSERVED order parameter phase-
    separates into bicontinuous domains that coarsen as L(t)~t^(1/3).
    ∂c/∂t = M ∇²(c³ - c - κ∇²c). Distinct from P2 MIPS (active, motility-induced,
    non-conserved free energy), P18 voter coarsening (non-conserved), and P3 Turing
    (Turing selects an INTRINSIC wavelength; spinodal domains coarsen without one)."""
    rng = np.random.default_rng(seed)
    c = rng.normal(0.0, 0.02, size=(G, G))
    frames: History = []
    snap = max(1, n_steps // n_frames)
    for t in range(n_steps):
        mu = c ** 3 - c - kappa * _lap(c)
        c = c + dt * M * _lap(mu)
        if t % snap == 0:
            frames.append({"field": c.copy()})
    frames.append({"field": c.copy()})
    return frames, None


def gray_scott_selfrep(seed: int = 0, G: int = 96, n_steps: int = 9000,
                       n_frames: int = 80, Du: float = 0.16, Dv: float = 0.08,
                       F: float = 0.0367, k: float = 0.0649, dt: float = 1.0) -> Built:
    """Gray-Scott reaction-diffusion in the self-replicating ("λ"/spot-mitosis)
    regime (Pearson 1993): localized spots GROW and DIVIDE, filling the domain by
    repeated replication. Distinct from P3 Turing (stationary patterns from a
    diffusion-driven instability) — here the spots are dynamic and self-replicate."""
    rng = np.random.default_rng(seed)
    u = np.ones((G, G)); v = np.zeros((G, G))
    r = 10; c = G // 2
    u[c - r:c + r, c - r:c + r] = 0.50
    v[c - r:c + r, c - r:c + r] = 0.25
    u = np.clip(u + 0.02 * rng.standard_normal((G, G)), 0, 1)
    v = np.clip(v + 0.02 * rng.standard_normal((G, G)), 0, 1)
    frames: History = []
    snap = max(1, n_steps // n_frames)
    for t in range(n_steps):
        uvv = u * v * v
        u = u + dt * (Du * _lap(u) - uvv + F * (1.0 - u))
        v = v + dt * (Dv * _lap(v) + uvv - (F + k) * v)
        if t % snap == 0:
            frames.append({"field": v.copy()})
    frames.append({"field": v.copy()})
    return frames, None


def swarmalators(seed: int = 0, N: int = 120, J: float = 1.0, K: float = 1.0,
                 dt: float = 0.1, n_steps: int = 1200, n_frames: int = 80) -> Built:
    """Swarmalators (O'Keeffe-Hong-Strogatz 2017): agents whose SPATIAL motion and
    internal PHASE are bidirectionally coupled — they swarm in space and sync in
    phase simultaneously, and the coupling correlates the two. Distinct from P9
    Kuramoto (phase only, no space) and P5 flocking (space/velocity only, no phase).
    Here (J=K=1) → static synchrony: a compact disk with phase coherence."""
    rng = np.random.default_rng(seed)
    x = rng.uniform(-1.0, 1.0, size=(N, 2))
    th = rng.uniform(-np.pi, np.pi, size=N)
    frames: History = []
    snap = max(1, n_steps // n_frames)
    for t in range(n_steps):
        dx = x[None, :, :] - x[:, None, :]        # x_j - x_i : (N,N,2)
        d = np.sqrt((dx ** 2).sum(-1))
        np.fill_diagonal(d, np.inf)               # self terms -> 0
        dth = th[None, :] - th[:, None]           # θ_j - θ_i : (N,N)
        attract = dx / d[..., None] * (1.0 + J * np.cos(dth)[..., None])
        repel = dx / (d[..., None] ** 2)
        vx = (attract - repel).sum(1) / N
        dthi = (np.sin(dth) / d).sum(1) * (K / N)
        x = x + dt * vx
        th = (th + dt * dthi) % (2 * np.pi)
        if t % snap == 0:
            frames.append({"positions": x.copy(), "phases": th.copy()})
    frames.append({"positions": x.copy(), "phases": th.copy()})
    return frames, None


def langton_ant(seed: int = 0, G: int = 128, n_steps: int = 12000,
                n_frames: int = 80) -> Built:
    """Langton's ant (1986): a single agent on a binary grid (white→turn right+flip+
    advance; black→turn left+flip+advance). After ~10⁴ chaotic-looking steps an
    ordered, recurrent diagonal "highway" emerges spontaneously — order from
    deterministic disorder. Deterministic (seed ignored). Distinct from every
    catalog pattern (single agent; emergent transport structure)."""
    grid = np.zeros((G, G), dtype=np.float32)
    x = y = G // 2
    dx, dy = 0, 1
    frames: History = []
    snap = max(1, n_steps // n_frames)
    for t in range(n_steps):
        if grid[x, y] == 0:
            dx, dy = dy, -dx                       # turn right
            grid[x, y] = 1.0
        else:
            dx, dy = -dy, dx                       # turn left
            grid[x, y] = 0.0
        x = (x + dx) % G; y = (y + dy) % G
        if t % snap == 0:
            frames.append({"grid": grid.copy()})
    frames.append({"grid": grid.copy()})
    return frames, None


def kauffman_rbn_critical(seed: int = 0, N: int = 100, K: int = 2,
                          n_steps: int = 200, n_frames: int = 80) -> Built:
    """Kauffman random Boolean network at the critical connectivity K=2: each node
    has K random inputs and a random Boolean function; synchronous update. K=2 sits
    at the order/chaos boundary — a frozen core coexists with a fluctuating periphery,
    settling onto short attractor cycles. Distinct from P16 Hopfield (designed
    attractors / associative memory). Emergence is dynamical, not spatial."""
    rng = np.random.default_rng(seed)
    inputs = np.array([rng.choice(N, size=K, replace=False) for _ in range(N)])
    tables = rng.integers(0, 2, size=(N, 1 << K))
    state = rng.integers(0, 2, size=N)
    powers = (1 << np.arange(K))[::-1]
    frames: History = []
    for t in range(n_steps):
        idx = (state[inputs] * powers).sum(1)      # per-node lookup index
        state = tables[np.arange(N), idx].astype(np.int64)
        frames.append({"state": state.copy()})
    return frames, None
