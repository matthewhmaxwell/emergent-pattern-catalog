"""Algorithmic-compressibility lens — Lempel-Ziv redundancy (a documented missing axis). The
covered MPR / permutation complexity uses short local windows and reads a deterministic-but-
aperiodic sequence (Thue-Morse, rule-90 Sierpinski) as high-entropy noise. LZ parsing captures
LONG-RANGE algorithmic redundancy across scales.

  lz_ratio — normalized LZ76 complexity of the symbol sequence divided by the mean LZ76 of its
             SHUFFLES. <1 = algorithmically structured (compressible beyond its own histogram);
             ~1 = no structure beyond the symbol frequencies (random). THE discriminator (the
             shuffle baseline removes the trivial frequency contribution).
  lz_norm  — LZ76 complexity normalized by n/log2(n) (1 ~ random, <1 structured); context.

Symbol sequence comes from an explicit 'symbols' (1D int/bool), else a binarized (median-split)
flattened field raster, else the binarized macrostate scalar series. Needs >=256 symbols (LZ is
length-hungry); returns None below. (Maximally random-looking deterministic systems -- rule 30,
logistic -- need an epsilon-machine grammar to recover; LZ catches detectable redundancy, the
honest scope.)
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional

import numpy as np


def _symbols(history: List[Dict[str, Any]]) -> Optional[np.ndarray]:
    for f in history:
        if isinstance(f, dict) and "symbols" in f:
            a = np.asarray(f["symbols"]).ravel()
            if a.size >= 256:
                return (a > np.median(a)).astype(np.int8) if a.dtype != np.int8 else a.astype(np.int8)
    # field raster (late frame), binarized
    fields = []
    for f in history:
        if isinstance(f, dict):
            for k in ("field", "u", "concentration"):
                if k in f:
                    a = np.asarray(f[k], dtype=float)
                    if a.ndim == 2:
                        fields.append(a.ravel())
                    break
    if fields:
        a = fields[-1]
        if a.size >= 256:
            return (a > np.median(a)).astype(np.int8)
    # macrostate scalar series
    s = []
    for f in history:
        if isinstance(f, dict):
            if "scalar" in f:
                s.append(float(f["scalar"]))
            elif "order_parameter" in f:
                s.append(float(f["order_parameter"]))
    s = np.asarray(s, dtype=float)
    if s.size >= 256 and np.ptp(s) > 0:
        return (s > np.median(s)).astype(np.int8)
    return None


def _lz76(s: np.ndarray) -> int:
    n = len(s); i = 0; c = 1; l = 1; k = 1; kmax = 1
    while l + k <= n:
        if s[i + k - 1] == s[l + k - 1]:
            k += 1
            if l + k > n:
                c += 1; break
        else:
            if k > kmax:
                kmax = k
            i += 1
            if i == l:
                c += 1; l += kmax; i = 0; k = 1; kmax = 1
            else:
                k = 1
    return c


def compressibility(history: List[Dict[str, Any]]) -> Optional[Dict[str, float]]:
    s = _symbols(history)
    if s is None:
        return None
    n = len(s)
    c = _lz76(s)
    lz_norm = c / (n / np.log2(n))
    rng = np.random.default_rng(0)
    sh = []
    for _ in range(8):
        p = s.copy(); rng.shuffle(p)
        sh.append(_lz76(p))
    base = float(np.mean(sh))
    lz_ratio = float(c / base) if base > 0 else 1.0
    return {"lz_ratio": round(lz_ratio, 3), "lz_norm": round(float(lz_norm), 3)}
