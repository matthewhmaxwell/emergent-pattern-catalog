"""Nonequilibrium-current / time-irreversibility lens — broken detailed balance (a documented
missing axis; NOTHING else in the battery measures time's arrow). A nonequilibrium steady state
carries a circulating probability current even when nothing structural changes and the spectrum
is flat; an equilibrium / detailed-balanced system is time-reversible (forward == reverse stats).

  circulation — normalized state-space current J = <A_t(B_{t+1}-B_t) - B_t(A_{t+1}-A_t)> /
                (std A * std B), the discrete line integral of (A dB - B dA): net rotation /
                enclosed-area rate in the 2-channel plane. |J| large = DIRECTED cycle
                (nonequilibrium flux); ~0 = reversible. THE discriminator -- it separates a
                directed cycle from a STANDING oscillation (both have a spectral peak; only the
                directed one circulates) and from equilibrium noise.

Two channels from: explicit 'chan_a'/'chan_b' (or 'channels' Tx2); else (COM_x, COM_y) from
'positions', (mean cos, mean sin) from 'phases', or (mean, std) of a 'field'. Needs >=12 frames.
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional, Tuple

import numpy as np


def _channels(history: List[Dict[str, Any]]) -> Optional[Tuple[np.ndarray, np.ndarray]]:
    A, B = [], []
    for f in history:
        if not isinstance(f, dict):
            continue
        if "chan_a" in f and "chan_b" in f:
            A.append(float(f["chan_a"])); B.append(float(f["chan_b"])); continue
        if "channels" in f:
            c = np.asarray(f["channels"], dtype=float).ravel()
            if c.size >= 2:
                A.append(float(c[0])); B.append(float(c[1])); continue
        if f.get("positions") is not None:
            p = np.asarray(f["positions"], dtype=float)
            if p.ndim == 2 and p.shape[1] == 2:
                com = p.mean(0); A.append(float(com[0])); B.append(float(com[1])); continue
        if "phases" in f:
            th = np.asarray(f["phases"], dtype=float).ravel()
            A.append(float(np.cos(th).mean())); B.append(float(np.sin(th).mean())); continue
        done = False
        for k in ("field", "u", "concentration"):
            if k in f:
                a = np.asarray(f[k], dtype=float)
                A.append(float(a.mean())); B.append(float(a.std())); done = True; break
        if done:
            continue
    if len(A) >= 12:
        return np.asarray(A), np.asarray(B)
    return None


def nonequilibrium_current(history: List[Dict[str, Any]]) -> Optional[Dict[str, float]]:
    ch = _channels(history)
    if ch is None:
        return None
    A, B = ch
    sa, sb = A.std(), B.std()
    if sa < 1e-9 or sb < 1e-9:
        return {"circulation": 0.0}
    A = A - A.mean(); B = B - B.mean()
    dA = np.diff(A); dB = np.diff(B)
    J = float(np.mean(A[:-1] * dB - B[:-1] * dA)) / (sa * sb)
    return {"circulation": round(abs(J), 4)}
