"""Fractal-dimension lens — scale-FREE spatial structure (box-counting).

STATUS: ADMITTED, SCOPED to sparse self-similar AGGREGATES (DLA-type). Box-counting D
ALONE is confounded — it conflates fractality with density/boundary (random gas D~1.57 ~=
Sierpinski D~1.585; percolation D~1.83 ~= uniform disk D~1.90). The fix is LACUNARITY
(gliding-box mass heterogeneity), which separates structured fractal gaps from random
texture: validated dla 5.69 / Sierpinski 4.09 vs random_gas 1.96 / uniform 1.54 (gap +2.1
gated by D>1.2 to exclude lines like a filament, lacun 41.9 / D 0.94). The discriminating
signal is therefore `lacunarity` gated by `fractal_dim > 1.2`, NOT D alone.
SCOPE LIMIT (honest): this catches SPARSE fractals (open aggregates). DENSE near-critical
fractals (percolation, fill ~0.5) are near-space-filling and homogeneous (lacunarity ~1.1),
so they fall OUTSIDE this lens — they need the full multifractal spectrum (f(alpha) / D_q),
a larger build deferred. Such systems are still DETECTED by the emergence indicator; they
just carry no dedicated fractal label here.

The complement of structure_factor: where S(k) keys on a CHARACTERISTIC length scale
(a peak in reciprocal space), this lens keys on the ABSENCE of one — self-similar
structure that repeats across scales, so a fractal slips straight through the
structure-factor lens. Box-counting dimension D of the occupied set:

    N(eps) ~ eps^(-D)            =>   D = slope of log N(eps) vs log(1/eps)

A space-filling region has D -> 2, a filament D -> 1, and a fractal a non-integer D
in between WITH a clean power-law (high R^2). Reads a binary/threshold 'grid' or
'field', or rasterizes 'positions' to an occupancy grid.

Scalars (averaged over the late window):
  fractal_dim — box-counting slope D.
  fit_r2      — R^2 of the log-log regression: how clean the self-similar scaling is.
                A clean fractal scales over many octaves (R^2 ~ 0.99); random/uniform
                fills are noisier or pinned at D ~ 2.
  fill_frac   — occupied fraction (context; D is only meaningful away from empty/full).

A structure reads "fractal/scale-free" when fractal_dim is sub-space-filling with a
high fit_r2. Returns None for substrates with no spatial occupancy.
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional

import numpy as np


def _occupancy(frame: Dict[str, Any], grid_n: int) -> Optional[np.ndarray]:
    if "grid" in frame:
        g = np.asarray(frame["grid"], dtype=float)
        if g.ndim != 2:
            return None
        return (g > 0).astype(int) if g.max() <= 1.0 else (g > g.mean()).astype(int)
    if "field" in frame:
        f = np.asarray(frame["field"], dtype=float)
        if f.ndim != 2:
            return None
        return (f > f.mean() + f.std()).astype(int)
    if "positions" in frame:
        p = np.asarray(frame["positions"], dtype=float)
        if p.ndim != 2 or p.shape[1] < 2 or len(p) < 5:
            return None
        p = p[:, :2]
        lo, span = p.min(0), (p.max(0) - p.min(0)); span[span == 0] = 1.0
        ij = np.clip(np.floor((p - lo) / span * (grid_n - 1)).astype(int), 0, grid_n - 1)
        occ = np.zeros((grid_n, grid_n), int); occ[ij[:, 0], ij[:, 1]] = 1
        return occ
    return None


def _box_count_dim(occ: np.ndarray) -> Optional[tuple]:
    n = min(occ.shape)
    occ = occ[:n, :n]
    sizes = []
    s = 1
    while s <= n // 2:
        sizes.append(s); s *= 2
    counts = []
    for s in sizes:
        nb = n // s
        if nb < 2:
            continue
        sub = occ[:nb * s, :nb * s].reshape(nb, s, nb, s)
        counts.append((s, int((sub.sum(axis=(1, 3)) > 0).sum())))
    counts = [(s, c) for s, c in counts if c > 0]
    if len(counts) < 3:
        return None
    s_arr = np.array([c[0] for c in counts], float)
    n_arr = np.array([c[1] for c in counts], float)
    x = np.log(1.0 / s_arr); y = np.log(n_arr)
    A = np.vstack([x, np.ones_like(x)]).T
    slope, intercept = np.linalg.lstsq(A, y, rcond=None)[0]
    yhat = A @ np.array([slope, intercept])
    ss_res = float(((y - yhat) ** 2).sum()); ss_tot = float(((y - y.mean()) ** 2).sum())
    r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0
    return float(slope), float(r2)


def _lacunarity(occ: np.ndarray, boxes=(2, 4, 8, 16)) -> float:
    """Gliding-box lacunarity Lambda(r) = E[S^2]/E[S]^2 (S = box mass), averaged over box
    sizes. Measures mass HETEROGENEITY at scale: a structured fractal (Sierpinski, DLA) is
    very uneven (high Lambda) while a uniform-random set of the SAME density is homogeneous
    (Lambda ~ 1). This is the axis that separates fractal from random texture where the
    box-counting dimension cannot (they can share D)."""
    from scipy.ndimage import uniform_filter
    a = occ.astype(float)
    vals = []
    for b in boxes:
        if b >= min(a.shape):
            continue
        S = uniform_filter(a, size=b, mode="constant") * float(b * b)   # box mass per position
        m1 = float(S.mean())
        if m1 <= 0:
            continue
        vals.append(float((S ** 2).mean() / (m1 ** 2)))
    return float(np.mean(vals)) if vals else 0.0


def fractal_dimension(history: List[Dict[str, Any]], grid_n: int = 128) -> Optional[Dict[str, float]]:
    frames = [f for f in history if isinstance(f, dict)]
    if not frames:
        return None
    late = frames[len(frames) // 2:] or frames
    dims, r2s, fills, lacs = [], [], [], []
    for f in late:
        occ = _occupancy(f, grid_n)
        if occ is None:
            continue
        fill = float(occ.mean())
        if fill < 0.002 or fill > 0.98:          # near-empty / near-full: D not meaningful
            continue
        res = _box_count_dim(occ)
        if res is None:
            continue
        d, r2 = res
        dims.append(d); r2s.append(r2); fills.append(fill); lacs.append(_lacunarity(occ))
    if not dims:
        return None
    return {"fractal_dim": round(float(np.mean(dims)), 4),
            "fit_r2": round(float(np.mean(r2s)), 4),
            "lacunarity": round(float(np.mean(lacs)), 4),
            "fill_frac": round(float(np.mean(fills)), 4)}
