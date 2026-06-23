"""Ring 1 — gap-targeted systems (constructed from empty ontology cells).

Each model is built to occupy a dimensional-coverage gap that the 32-pattern
catalog leaves EMPTY (see coverage_map.py): the periodic-table logic — search
where the ontology predicts a pattern *could* live but the catalog doesn't cover.
Run through the instrument; the target read is EMERGENT-UNCLASSIFIED (emergence
detected, no false MATCH) -> candidate P34+.
"""
from __future__ import annotations

from typing import Any, Dict, List, Tuple

import numpy as np

Built = Tuple[List[Dict[str, Any]], Dict[str, Any]]


def adaptive_voter(seed: int = 0, N: int = 200, mean_deg: int = 8, phi: float = 0.7,
                   max_steps: int = 60000, n_frames: int = 80) -> Built:
    """GAP: substrate = EVOLVING NETWORK (empty in the catalog).

    Adaptive voter model (Holme & Newman 2006; Vazquez et al. 2008): N nodes carry
    a binary opinion on a network. Repeatedly pick a discordant edge (i,j); with
    prob phi REWIRE it (i drops j, links to a random same-opinion node), else i
    ADOPTS j's opinion. Above a critical phi the network spontaneously FRAGMENTS
    into disconnected same-opinion communities — the emergence is in the TOPOLOGY,
    which coevolves with the node states. Distinct from P18 voter (fixed lattice),
    P21 polarization (fixed network), P22 cascade (fixed network)."""
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
                frames.append({"adjacency": A.astype(float).copy()})
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
            frames.append({"adjacency": A.astype(float).copy()})
    frames.append({"adjacency": A.astype(float).copy()})
    return frames, {"model": "adaptive_voter", "phi": phi, "N": N, "mean_deg": mean_deg}


# Null control for the evolving-network arm: a NON-fragmenting rewiring (random
# rewiring, opinion-blind) keeps the graph Erdos-Renyi-like — no community emerges.
def null_random_rewire(seed: int = 0, N: int = 200, mean_deg: int = 8,
                       max_steps: int = 60000, n_frames: int = 80) -> Built:
    """Opinion-blind random rewiring: pick a random edge, move one endpoint to a
    random node. The degree sequence drifts but NO community structure emerges
    (stays ER-like). Evolving-network NULL -> must read NO-EMERGENCE."""
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
