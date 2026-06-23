"""Adaptive (coevolving) network models for P34 — the evolving-network substrate.

The canonical positive is the adaptive voter (Holme & Newman 2006; Vazquez et al.
2008): N nodes carry a binary opinion on a network; a discordant edge (i,j) is, with
prob phi, REWIRED (i drops j, links to a random same-opinion node), else i ADOPTS
j's opinion. Above the fragmentation transition (phi~0.99 here for mean degree 4) the
network spontaneously FRAGMENTS into disconnected same-opinion communities — the
emergence is in the TOPOLOGY, which coevolves with the node states.

Frames carry {"adjacency": NxN, "opinions": N}. The `neg_*` generators are the
Phase-2a lookalikes (negative panel): topology-change-without-community, static
modular, fixed-topology voter, static ER, pre-fragmented.
"""
from __future__ import annotations

from typing import Any, Dict, List, Tuple

import numpy as np

Built = Tuple[List[Dict[str, Any]], Dict[str, Any]]


def adaptive_voter(seed: int = 0, N: int = 200, mean_deg: int = 4, phi: float = 0.99,
                   max_steps: int = 400000, n_frames: int = 80) -> Built:
    """Adaptive voter -> fragmentation transition (the evolving-network positive)."""
    rng = np.random.default_rng(seed)
    op = rng.integers(0, 2, N)
    p = mean_deg / (N - 1)
    A = (np.triu(rng.random((N, N)), 1) < p).astype(np.int8)
    A = A + A.T
    frames: List[Dict[str, Any]] = []
    snap = max(1, max_steps // n_frames)
    for step in range(max_steps):
        i = int(rng.integers(N))
        nbrs = np.nonzero(A[i])[0]
        if nbrs.size == 0:
            if step % snap == 0:
                frames.append({"adjacency": A.astype(float).copy(), "opinions": op.copy()})
            continue
        j = int(nbrs[rng.integers(nbrs.size)])
        if op[i] != op[j]:
            if rng.random() < phi:
                same = np.nonzero(op == op[i])[0]
                same = same[(same != i)]
                if same.size:
                    l = int(same[rng.integers(same.size)])
                    A[i, j] = A[j, i] = 0
                    A[i, l] = A[l, i] = 1
            else:
                op[i] = op[j]
        if step % snap == 0:
            frames.append({"adjacency": A.astype(float).copy(), "opinions": op.copy()})
    frames.append({"adjacency": A.astype(float).copy(), "opinions": op.copy()})
    return frames, {"model": "adaptive_voter", "phi": phi, "N": N, "mean_deg": mean_deg}


def null_random_rewire(seed: int = 0, N: int = 200, mean_deg: int = 4,
                       max_steps: int = 400000, n_frames: int = 80) -> Built:
    """Opinion-blind random rewiring: topology drifts but NO community emerges
    (stays Erdos-Renyi-like). Evolving-network NULL -> NO-EMERGENCE."""
    rng = np.random.default_rng(seed)
    p = mean_deg / (N - 1)
    A = (np.triu(rng.random((N, N)), 1) < p).astype(np.int8)
    A = A + A.T
    frames: List[Dict[str, Any]] = []
    snap = max(1, max_steps // n_frames)
    for step in range(max_steps):
        i = int(rng.integers(N))
        nbrs = np.nonzero(A[i])[0]
        if nbrs.size:
            j = int(nbrs[rng.integers(nbrs.size)])
            l = int(rng.integers(N))
            if l != i and A[i, l] == 0:
                A[i, j] = A[j, i] = 0
                A[i, l] = A[l, i] = 1
        if step % snap == 0:
            frames.append({"adjacency": A.astype(float).copy()})
    frames.append({"adjacency": A.astype(float).copy()})
    return frames, {"model": "null_random_rewire", "N": N}


def neg_static_modular(seed: int = 0, N: int = 200, n_frames: int = 40) -> Built:
    """Pre-built 2-community SBM held FIXED: high Q but NO emergence (no gain)."""
    rng = np.random.default_rng(seed)
    op = (np.arange(N) >= N // 2).astype(int)
    same = op[:, None] == op[None, :]
    P = np.where(same, 0.10, 0.005)
    A = (np.triu(rng.random((N, N)), 1) < np.triu(P, 1)).astype(float)
    A = A + A.T
    return [{"adjacency": A.copy(), "opinions": op.copy()} for _ in range(n_frames)], \
        {"model": "static_modular", "N": N}


def neg_fixed_network_voter(seed: int = 0, N: int = 200, mean_deg: int = 4,
                            n_frames: int = 60) -> Built:
    """Voter dynamics on a FIXED ER network (no rewiring): states change, topology
    constant -> Q_gain ~ 0 (P18-style on a network)."""
    rng = np.random.default_rng(seed)
    p = mean_deg / (N - 1)
    A = (np.triu(rng.random((N, N)), 1) < p).astype(float)
    A = A + A.T
    op = rng.integers(0, 2, N)
    frames, steps = [], 60000
    snap = max(1, steps // n_frames)
    for t in range(steps):
        i = int(rng.integers(N))
        nbrs = np.nonzero(A[i])[0]
        if nbrs.size:
            op[i] = op[int(nbrs[rng.integers(nbrs.size)])]
        if t % snap == 0:
            frames.append({"adjacency": A.copy(), "opinions": op.copy()})
    frames.append({"adjacency": A.copy(), "opinions": op.copy()})
    return frames, {"model": "fixed_network_voter", "N": N}


def neg_er_static(seed: int = 0, N: int = 200, mean_deg: int = 4, n_frames: int = 40) -> Built:
    """Static Erdos-Renyi graph held fixed: low Q, no structure, no emergence."""
    rng = np.random.default_rng(seed)
    p = mean_deg / (N - 1)
    A = (np.triu(rng.random((N, N)), 1) < p).astype(float)
    A = A + A.T
    return [{"adjacency": A.copy()} for _ in range(n_frames)], {"model": "er_static", "N": N}


def neg_pre_fragmented(seed: int = 0, N: int = 200, mean_deg: int = 4, n_frames: int = 40) -> Built:
    """Two disconnected ER components from the START, held fixed: high Q from t=0
    -> Q_gain ~ 0 (structure present but not emergent)."""
    rng = np.random.default_rng(seed)
    h = N // 2
    p = mean_deg / (h - 1)
    A = np.zeros((N, N))
    for off in (0, h):
        b = (np.triu(rng.random((h, h)), 1) < p).astype(float)
        A[off:off + h, off:off + h] = b + b.T
    return [{"adjacency": A.copy()} for _ in range(n_frames)], {"model": "pre_fragmented", "N": N}
