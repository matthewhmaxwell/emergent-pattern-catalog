"""Network-emergence metrics for P34 (adaptive-network fragmentation).

Operates on an 'adjacency' observable (NxN). Reuses the spectral modularity from
the round-2 network channel; adds connected-component fragmentation and
state-topology assortativity (homophily) for the coevolution signature.
"""
from __future__ import annotations

from typing import Any, Dict, Optional

import numpy as np
from scipy.sparse import csr_matrix
from scipy.sparse.csgraph import connected_components

from epc.phase2a.info_channels import _spectral_modularity


def binarize(f: Dict[str, Any]) -> Optional[np.ndarray]:
    """Symmetric binary adjacency from a frame's 'adjacency' observable."""
    if not isinstance(f, dict) or "adjacency" not in f:
        return None
    A = np.asarray(f["adjacency"], dtype=float)
    if A.ndim != 2 or A.shape[0] != A.shape[1]:
        return None
    A = ((A + A.T) > 0).astype(float)
    np.fill_diagonal(A, 0.0)
    return A


def modularity(A: np.ndarray) -> float:
    """Newman leading-eigenvector 2-partition modularity Q (>0 = community structure)."""
    return _spectral_modularity(A)


def giant_fraction(A: np.ndarray) -> float:
    """Largest connected component as a fraction of N (1.0 = connected; ~0.5 = split)."""
    n = A.shape[0]
    _, lab = connected_components(csr_matrix(A > 0), directed=False)
    return float(np.bincount(lab).max() / n) if n else 1.0


def state_assortativity(A: np.ndarray, state: np.ndarray) -> float:
    """Fraction of edges joining same-state nodes (homophily). ~0.5 = chance for
    balanced binary; ->1 = state-segregated topology (opinion-driven rewiring)."""
    i, j = np.nonzero(np.triu(A, 1) > 0)
    if i.size == 0:
        return 0.5
    return float(np.mean(state[i] == state[j]))
