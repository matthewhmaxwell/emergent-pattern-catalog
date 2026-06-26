"""Persistent homology — a Tier-2 topology lens (two substrate paths).

Captures the one structural family every other lens misses: TOPOLOGY (loops/voids/holes).
The density (clustering), scale (structure-factor), and orientation lenses are all blind
to a hole in the configuration.

PATH 1 — POSITIONS (Vietoris-Rips via ripser). Scalars (late-window mean):
  h1_max   — persistence of the single most prominent loop; scale-NORMALIZED. THE
             discriminator: a genuine ring/void scores high (noisy-ring control 0.58,
             vortex 0.43) while blobs / gas / nulls / lattices stay <= 0.12 — a clean gap.
  h1_total — SUM of 1-cycle persistences; noise-confounded, secondary only.
  n_components — persistent H0 components.

PATH 2 — FIELD / GRID (superlevel-set hole count via scipy.ndimage). Closes the field
topology gap the positions path could not reach (RD labyrinths, Lenia creatures): a 2D
field has no point cloud, but its superlevel sets DO have holes. Sweeping the threshold,
count enclosed background components that are PERSISTENT — enclosed across a large fraction
of the threshold sweep (a structural void stays enclosed across many levels; a noise pocket
is born and dies within one or two). This is the field analogue of the h1_max lesson: a raw
hole COUNT is noise-confounded (iid noise yields hundreds of transient pockets), so we
filter by PERSISTENCE, not count. Scalars (late-window mean):
  field_loops     — number of persistent enclosed holes (b1 of the foreground, denoised).
  field_loop_area — fraction of grid area inside those holes (a persistence-style weight).
A labyrinth / Swiss-cheese scores high; spots / blobs / smooth fields / noise score ~0.

Requires `ripser` (positions path) and scipy.ndimage (field path); both in the venv.
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional

import numpy as np

try:
    from ripser import ripser
    _HAVE_RIPSER = True
except Exception:
    _HAVE_RIPSER = False

try:
    from scipy import ndimage as _ndi
    _HAVE_NDI = True
except Exception:
    _HAVE_NDI = False


def _field_loops(field: np.ndarray, n_levels: int = 16, min_hole: int = 9,
                 persist_frac: float = 0.35) -> tuple:
    """PERSISTENT enclosed-hole count of the superlevel set. A pixel is 'enclosed' at a
    threshold if it is in a background component not touching the border; we accumulate how
    many levels each pixel is enclosed, keep pixels enclosed across >= persist_frac of the
    sweep (structural voids, not transient noise pockets), and count the resulting size-
    filtered components. Returns (persistent_holes, hole_area_fraction)."""
    a = np.asarray(field, dtype=float)
    if a.ndim != 2 or a.size == 0:
        return None
    a = _ndi.gaussian_filter(a, 2.0)                     # denoise: kills spatially-uncorrelated
    lo, hi = float(a.min()), float(a.max())              # (iid) pits; smooth fields survive
    if hi - lo < 1e-9:
        return 0, 0.0
    levels = np.linspace(lo, hi, n_levels)[1:-1]
    enclosed = np.zeros(a.shape, dtype=int)
    for tau in levels:
        bg = a <= tau
        lab, n = _ndi.label(bg)
        if n == 0:
            continue
        border = np.unique(np.concatenate([lab[0, :], lab[-1, :], lab[:, 0], lab[:, -1]]))
        border = border[border > 0]
        enclosed += (bg & (lab > 0) & ~np.isin(lab, border))
    K = max(1, int(persist_frac * len(levels)))
    lab2, n2 = _ndi.label(enclosed >= K)
    if n2 == 0:
        return 0, 0.0
    sizes = _ndi.sum(np.ones_like(lab2), lab2, index=range(1, n2 + 1))
    keep = sizes >= min_hole
    return int(keep.sum()), float(sizes[keep].sum()) / a.size


def persistent_homology(history: List[Dict[str, Any]], max_points: int = 200,
                        seed: int = 0) -> Optional[Dict[str, float]]:
    frames = [f for f in history if isinstance(f, dict)]
    pos_frames = [f for f in frames if "positions" in f]
    field_frames = [f for f in frames if ("field" in f or "grid" in f) and "positions" not in f]

    # PATH 1 — positions (Rips)
    if _HAVE_RIPSER and len(pos_frames) >= 2:
        rng = np.random.default_rng(seed)
        late = pos_frames[len(pos_frames) // 2:]
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
            pn = (p - p.min(0)) / ext
            dgms = ripser(pn, maxdim=1)["dgms"]
            h1 = dgms[1]
            if len(h1):
                pe = h1[:, 1] - h1[:, 0]; pe = pe[np.isfinite(pe)]
                h1t.append(float(pe.sum())); h1m.append(float(pe.max()) if len(pe) else 0.0)
            else:
                h1t.append(0.0); h1m.append(0.0)
            h0 = dgms[0]; pe0 = h0[:, 1] - h0[:, 0]; pe0 = pe0[np.isfinite(pe0)]
            h0n.append(int((pe0 > 0.05).sum()))
        if h1t:
            return {"kind": "rips", "h1_total": round(float(np.mean(h1t)), 4),
                    "h1_max": round(float(np.mean(h1m)), 4),
                    "n_components": round(float(np.mean(h0n)), 2)}

    # PATH 2 — field / grid (superlevel-set holes)
    if _HAVE_NDI and len(field_frames) >= 2:
        late = field_frames[len(field_frames) // 2:]
        loops, areas = [], []
        for f in late:
            g = f.get("field", f.get("grid"))
            res = _field_loops(np.asarray(g, dtype=float))
            if res is None:
                continue
            loops.append(res[0]); areas.append(res[1])
        if loops:
            return {"kind": "field", "field_loops": round(float(np.mean(loops)), 2),
                    "field_loop_area": round(float(np.mean(areas)), 4)}
    return None
