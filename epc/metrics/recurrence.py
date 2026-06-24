"""Recurrence quantification — a Tier-2 lens imported from nonlinear dynamics.

Recurrence Quantification Analysis (RQA; Marwan et al. 2007) measures whether a
system's state TRAJECTORY revisits earlier states and evolves deterministically
from them. The two headline measures:
  DET (determinism) — fraction of recurrence points lying on diagonal lines;
      high for deterministic/structured dynamics (periodic, chaotic-attractor,
      patterned), low for stochastic churn.
  LAM (laminarity) — fraction on vertical lines; high for intermittent / trapped
      (laminar) states.

This adds an axis the catalog lacks: deterministic structure vs. stochastic churn.
It is computed on a TRANSLATION + PERMUTATION-INVARIANT structure trajectory (each
frame's component values are sorted and mean-centered), so global drift and
component relabeling do not fool it (the failure mode that limited the 1-D MPR
macro). eps is auto-set to a target recurrence rate (standard RQA practice), so
DET/LAM are comparable across systems.

Refs: Marwan, Romano, Thiel, Kurths (2007), "Recurrence plots for the analysis of
complex systems," Physics Reports 438.
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional

import numpy as np

from epc.phase2a.info_channels import micro_macro


def _invariant_trajectory(history: List[Dict[str, Any]]) -> Optional[np.ndarray]:
    """(T, d) trajectory of sorted+centered configuration vectors from any history."""
    micro, _ = micro_macro(history)
    if micro is None:
        return None
    M = np.asarray(micro, dtype=float)            # (T, n)
    if M.shape[0] < 8 or M.shape[1] < 2:
        return None
    S = np.sort(M, axis=1)                          # permutation-invariant
    S = S - S.mean(axis=1, keepdims=True)           # translation-invariant
    return S


def _line_fraction(R: np.ndarray, lmin: int, vertical: bool) -> float:
    """Fraction of recurrence points on diagonal (or vertical) lines of length >= lmin."""
    M = R.T if vertical else R
    total = float(R.sum())
    if total <= 0:
        return 0.0
    on_lines = 0
    T = M.shape[0]
    # diagonals of M (offsets), skip the main diagonal (k=0) for DET
    offsets = range(-(T - 1), T) if not vertical else [0]
    if vertical:
        # vertical lines = consecutive 1s down each column
        for c in range(M.shape[1]):
            run = 0
            for r in range(M.shape[0]):
                if M[r, c]:
                    run += 1
                else:
                    if run >= lmin:
                        on_lines += run
                    run = 0
            if run >= lmin:
                on_lines += run
        return on_lines / total
    for k in offsets:
        if k == 0:
            continue
        diag = np.diagonal(M, offset=k)
        run = 0
        for v in diag:
            if v:
                run += 1
            else:
                if run >= lmin:
                    on_lines += run
                run = 0
        if run >= lmin:
            on_lines += run
    return on_lines / total


def rqa(traj: np.ndarray, target_rr: float = 0.10, lmin: int = 2,
        theiler: Optional[int] = None) -> Dict[str, float]:
    """RQA on a (T, d) trajectory. eps auto-set to hit target recurrence rate.

    A THEILER window excludes temporally-adjacent pairs (|i-j| <= theiler) so DET
    counts only NON-TRIVIAL recurrences (revisiting distant states). Without it a
    smooth stochastic trajectory (random walk) scores high DET purely from
    continuity; the window is standard practice for continuous-time data."""
    T = traj.shape[0]
    if theiler is None:
        theiler = max(2, T // 12)
    diff = traj[:, None, :] - traj[None, :, :]
    D = np.sqrt((diff ** 2).sum(-1))               # (T, T)
    off = D[~np.eye(T, dtype=bool)]
    if off.size == 0 or off.max() <= 0:
        return {"DET": 0.0, "LAM": 0.0, "RR": 0.0}
    eps = float(np.quantile(off, target_rr))
    R = (D <= eps).astype(np.int8)
    ii, jj = np.indices((T, T))
    R[np.abs(ii - jj) <= theiler] = 0              # Theiler window: drop trivial near-diagonal recurrences
    denom = float(R.size - np.count_nonzero(np.abs(ii - jj) <= theiler))
    rr = float(R.sum() / denom) if denom > 0 else 0.0
    return {"DET": round(_line_fraction(R, lmin, vertical=False), 4),
            "LAM": round(_line_fraction(R, lmin, vertical=True), 4),
            "RR": round(rr, 4), "eps": round(eps, 4)}


def recurrence_determinism(history: List[Dict[str, Any]]) -> Optional[Dict[str, float]]:
    """RQA from a raw observation history (None if no usable trajectory)."""
    traj = _invariant_trajectory(history)
    if traj is None:
        return None
    return rqa(traj)
