"""Cross-substrate coupling lens — joint structure provably ABSENT from either marginal (a
documented missing axis: phenomena that live across two substrates at once). Single-substrate
lenses read each field alone; this reads their coupling. Two estimators on a pair of co-evolving
fields A, B:

  cross_coherence  — max radial spectral coherence |<S_AB(k)>|^2 / (<S_AA(k)><S_BB(k)>) over k>0,
                     azimuthally averaged. Coherence normalizes OUT the marginal power, so it is
                     ~1 where A,B are linearly related at a wavelength even if each auto-spectrum
                     is FLAT (reciprocal cross-diffusion patterning), and ~0 for independent
                     fields -- including two fields that are each individually patterned but
                     uncorrelated. THE discriminator (a joint quantity, zero in each marginal).
  loop_coupling    — min(TE(A->B), TE(B->A)) on the field means: a CLOSED feedback loop (stigmergy
                     / sensorimotor) has BOTH directions nonzero, distinct from one-way drive
                     (which directed_info_flow owns) and from independence.

Reads two 2D fields from a frame: ('field_a','field_b'), ('u','v'), ('field','field2'), or the
first two distinct 2D arrays present. Returns None if a pair isn't available.
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional, Tuple

import numpy as np

_PAIRS = [("field_a", "field_b"), ("u", "v"), ("field", "field2"), ("concentration", "field2")]


def _two_fields(frame: Dict[str, Any]) -> Optional[Tuple[np.ndarray, np.ndarray]]:
    if not isinstance(frame, dict):
        return None
    for ka, kb in _PAIRS:
        if ka in frame and kb in frame:
            a = np.asarray(frame[ka], dtype=float); b = np.asarray(frame[kb], dtype=float)
            if a.ndim == 2 and b.ndim == 2 and a.shape == b.shape:
                return a, b
    grids = [np.asarray(v, dtype=float) for v in frame.values()
             if isinstance(v, (list, np.ndarray)) and np.asarray(v).ndim == 2]
    grids = [g for g in grids if min(g.shape) >= 12]
    if len(grids) >= 2 and grids[0].shape == grids[1].shape:
        return grids[0], grids[1]
    return None


def _radial_coherence(A: np.ndarray, B: np.ndarray) -> float:
    A = A - A.mean(); B = B - B.mean()
    if not np.any(A) or not np.any(B):
        return 0.0
    FA = np.fft.fftshift(np.fft.fft2(A)); FB = np.fft.fftshift(np.fft.fft2(B))
    Sab = FA * np.conj(FB); Saa = np.abs(FA) ** 2; Sbb = np.abs(FB) ** 2
    ny, nx = A.shape; cy, cx = ny // 2, nx // 2
    Y, X = np.indices(A.shape); R = np.hypot(X - cx, Y - cy).astype(int)
    rmax = min(cy, cx) - 1
    best = 0.0
    for r in range(1, rmax):
        m = R == r
        if m.sum() < 4:
            continue
        num = np.abs(Sab[m].sum()) ** 2
        den = Saa[m].sum() * Sbb[m].sum()
        if den > 0:
            best = max(best, float(num / den))
    return best


def _te(x: np.ndarray, y: np.ndarray, bins: int = 4) -> float:
    # transfer entropy x->y, order 1, on binned series
    if len(y) < 12:
        return 0.0
    xb = np.clip((np.digitize(x, np.quantile(x, np.linspace(0, 1, bins + 1)[1:-1]))), 0, bins - 1)
    yb = np.clip((np.digitize(y, np.quantile(y, np.linspace(0, 1, bins + 1)[1:-1]))), 0, bins - 1)
    yf, yp, xp = yb[1:], yb[:-1], xb[:-1]
    from collections import Counter
    c_yp = Counter(zip(yp)); c_yfyp = Counter(zip(yf, yp)); c_ypxp = Counter(zip(yp, xp))
    c_all = Counter(zip(yf, yp, xp)); n = len(yf)
    te = 0.0
    for (a, b, cc), cnt in c_all.items():
        p_all = cnt / n
        p1 = c_yfyp[(a, b)] / n; p2 = c_ypxp[(b, cc)] / n; p3 = c_yp[(b,)] / n
        if p1 > 0 and p2 > 0 and p3 > 0:
            te += p_all * np.log2((p_all * p3) / (p1 * p2))
    return max(0.0, float(te))


def cross_substrate(history: List[Dict[str, Any]]) -> Optional[Dict[str, float]]:
    frames = [f for f in history if _two_fields(f) is not None]
    if not frames:
        return None
    late = frames[len(frames) // 2:] or frames
    cohs, ma, mb = [], [], []
    for f in late:
        A, B = _two_fields(f)
        cohs.append(_radial_coherence(A, B))
        ma.append(float(A.mean())); mb.append(float(B.mean()))
    out = {"cross_coherence": round(float(np.mean(cohs)), 4)}
    if len(ma) >= 12:
        a = np.asarray(ma); b = np.asarray(mb)
        if a.std() > 1e-9 and b.std() > 1e-9:
            out["loop_coupling"] = round(min(_te(a, b), _te(b, a)), 4)
    return out
