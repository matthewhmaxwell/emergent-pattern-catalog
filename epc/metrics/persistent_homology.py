"""Persistent homology — a Tier-2 topology lens (Vietoris-Rips via ripser).

Captures the one structural family every other lens misses: TOPOLOGY. H1 persistence
counts loops/voids/holes (rings, membranes, porous/cellular structures); H0 counts
connected components. The density (clustering), scale (structure-factor), and
orientation lenses are all blind to a hole in the configuration.

Scalars (averaged over the late window):
  h1_max   — persistence of the single most prominent loop; scale-NORMALIZED. THE
             discriminator: a genuine ring/void scores high (validated: noisy-ring
             control 0.58, vortex 0.43) while blobs / gas / nulls / lattices stay
             <= 0.12 — a clean gap. This is the admitted topology axis.
  h1_total — SUM of all 1-cycle persistences. NOTE: noise-confounded — a dense random
             cloud has many spurious tiny holes whose persistences sum high and swamp
             the one real loop; kept only as a secondary, NOT the discriminator.
  n_components — persistent H0 components (topological clusters).

Positions are normalized to a unit bounding box per frame, so persistence is in
scale-invariant units. Requires `ripser` (project venv).
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional

import numpy as np

try:
    from ripser import ripser
    _HAVE_RIPSER = True
except Exception:
    _HAVE_RIPSER = False


def persistent_homology(history: List[Dict[str, Any]], max_points: int = 200,
                        seed: int = 0) -> Optional[Dict[str, float]]:
    if not _HAVE_RIPSER:
        return None
    frames = [f for f in history if isinstance(f, dict) and "positions" in f]
    if len(frames) < 2:
        return None
    rng = np.random.default_rng(seed)
    late = frames[len(frames) // 2:]
    h1t, h1m, h0n = [], [], []
    for f in late:
        p = np.asarray(f["positions"], dtype=float)
        if p.ndim != 2 or p.shape[1] < 2:
            return None
        p = p[:, :2]
        if len(p) < 5:
            continue
        if len(p) > max_points:
            p = p[rng.choice(len(p), max_points, replace=False)]
        ext = float(np.ptp(p, axis=0).max()) + 1e-9
        pn = (p - p.min(0)) / ext                       # unit box -> scale-invariant
        dgms = ripser(pn, maxdim=1)["dgms"]
        h1 = dgms[1]
        if len(h1):
            pe = h1[:, 1] - h1[:, 0]; pe = pe[np.isfinite(pe)]
            h1t.append(float(pe.sum())); h1m.append(float(pe.max()) if len(pe) else 0.0)
        else:
            h1t.append(0.0); h1m.append(0.0)
        h0 = dgms[0]; pe0 = h0[:, 1] - h0[:, 0]; pe0 = pe0[np.isfinite(pe0)]
        h0n.append(int((pe0 > 0.05).sum()))
    if not h1t:
        return None
    return {"h1_total": round(float(np.mean(h1t)), 4),
            "h1_max": round(float(np.mean(h1m)), 4),
            "n_components": round(float(np.mean(h0n)), 2)}
