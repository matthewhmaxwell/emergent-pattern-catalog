"""Graph-structure lens — a Tier-2 network-topology lens (networkx).

Covers the substrate class every positions/field lens is blind to: explicit
interaction GRAPHS. Reads frames carrying an 'adjacency' matrix and reports the two
structural axes that separate an organized network from a random one:

  degree_cv  — coefficient of variation of the degree sequence. HIGH when a few hubs
               dominate (scale-free / preferential attachment); low for an
               Erdős–Rényi random graph at the same density (Poisson degrees).
  modularity — best-partition modularity (greedy). HIGH for community (block) AND for
               ring-lattice / small-world structure. NOTE: greedy modularity over-fits
               to a finite-size floor (~0.35) even on a random graph, so the structured
               threshold is ~0.46, not 0. clustering disambiguates the two structured
               families: small-world high (~0.45), block low (~0.13).

A graph is "structured" if it is high on EITHER axis; the Erdős–Rényi null is low on
degree_cv (~0.4), modularity (~0.35 floor) AND clustering (~0.06). clustering,
assortativity and mean_degree are reported as secondary context. Validated gaps:
degree_cv hubs 0.89-1.00 vs random 0.39-0.45 (+0.45); modularity communities 0.54-0.58
vs random 0.35-0.39 (+0.14). degree_cv catches hubs, modularity catches communities;
neither alone covers both families, which is why both are reported.

Requires networkx (project venv). Returns None for non-graph substrates.
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional

import numpy as np

try:
    import networkx as nx
    _HAVE_NX = True
except Exception:
    _HAVE_NX = False


def _gini(x: np.ndarray) -> float:
    x = np.sort(np.asarray(x, dtype=float))
    n = len(x); s = x.sum()
    if n == 0 or s == 0:
        return 0.0
    idx = np.arange(1, n + 1)
    return float((2 * (idx * x).sum() - (n + 1) * s) / (n * s))


def graph_structure(history: List[Dict[str, Any]]) -> Optional[Dict[str, float]]:
    if not _HAVE_NX:
        return None
    frames = [f for f in history if isinstance(f, dict) and "adjacency" in f]
    if not frames:
        return None
    A = np.asarray(frames[-1]["adjacency"], dtype=float)        # mature graph
    if A.ndim != 2 or A.shape[0] != A.shape[1] or A.shape[0] < 4:
        return None
    A = ((A + A.T) > 0).astype(int)                             # symmetric, binary
    np.fill_diagonal(A, 0)
    deg = A.sum(1)
    if deg.sum() / 2.0 < 1:
        return None
    mean_deg = float(deg.mean())
    degree_cv = float(deg.std() / mean_deg) if mean_deg > 0 else 0.0
    G = nx.from_numpy_array(A)
    try:
        comms = nx.community.greedy_modularity_communities(G)
        modularity = float(nx.community.modularity(G, comms))
    except Exception:
        modularity = 0.0
    clustering = float(nx.average_clustering(G))
    try:
        assort = float(nx.degree_assortativity_coefficient(G))
        if not np.isfinite(assort):
            assort = 0.0
    except Exception:
        assort = 0.0
    return {"degree_cv": round(degree_cv, 4),
            "degree_gini": round(_gini(deg), 4),
            "modularity": round(modularity, 4),
            "clustering": round(clustering, 4),
            "assortativity": round(assort, 4),
            "mean_degree": round(mean_deg, 2)}
