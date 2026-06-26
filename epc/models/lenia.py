"""Lenia — continuous cellular automaton (Bert Wang-Chak Chan 2019), a multi-family
substrate of autonomous MOVING creatures. The canonical Orbium glider is a localized
coherent structure that self-propels — the qualitative regime RD (static Turing) and
Kuramoto (phase waves) lack. It is also ASAL's own search space (Kumar et al. 2024), so
running our interpretable lens battery on it is directly relevant to the differentiation
narrative (named/validated axes vs one opaque foundation-model descriptor).

    U = K * A            (ring-kernel convolution, FFT, periodic)
    A <- clip(A + (1/T) * (bell(U, mu, sigma)*2 - 1), 0, 1)

Emits the 'field' A per frame. Orbium params: R=13, T=10, mu=0.15, sigma=0.015,
single Gaussian-ring kernel bell(r, 0.5, 0.15). Seed = the published 20x20 Orbium.
"""
from __future__ import annotations

from typing import Any, Dict, List

import numpy as np

_ORBIUM = np.array([
    [0,0,0,0,0,0,0.1,0.14,0.1,0,0,0.03,0.03,0,0,0.3,0,0,0,0],
    [0,0,0,0,0,0.08,0.24,0.3,0.3,0.18,0.14,0.15,0.16,0.15,0.09,0.2,0,0,0,0],
    [0,0,0,0,0,0.15,0.34,0.44,0.46,0.38,0.18,0.14,0.11,0.13,0.19,0.18,0.45,0,0,0],
    [0,0,0,0,0.06,0.13,0.39,0.5,0.5,0.37,0.06,0,0,0,0.02,0.16,0.68,0,0,0],
    [0,0,0,0,0.11,0.17,0.17,0.33,0.4,0.38,0.28,0.14,0,0,0,0,0.18,0.42,0,0],
    [0,0,0.09,0.18,0.13,0.06,0.08,0.26,0.32,0.32,0.27,0,0,0,0,0,0,0.82,0,0],
    [0.27,0,0.16,0.12,0,0,0,0.25,0.38,0.44,0.45,0.34,0,0,0,0,0,0.22,0.17,0],
    [0,0.07,0.2,0.02,0,0,0,0.31,0.48,0.57,0.6,0.57,0,0,0,0,0,0,0.49,0],
    [0,0.59,0.19,0,0,0,0,0.2,0.57,0.69,0.76,0.76,0.49,0,0,0,0,0,0.36,0],
    [0,0.58,0.19,0,0,0,0,0,0.67,0.83,0.9,0.92,0.87,0.12,0,0,0,0,0.22,0.07],
    [0,0,0.46,0,0,0,0,0,0.7,0.93,1,1,1,0.61,0,0,0,0,0.18,0.11],
    [0,0,0.82,0,0,0,0,0,0.47,1,1,0.98,1,0.96,0.27,0,0,0,0.19,0.1],
    [0,0,0.46,0,0,0,0,0,0.25,1,1,0.84,0.92,0.97,0.54,0.14,0.04,0.1,0.21,0.05],
    [0,0,0,0.4,0,0,0,0,0.09,0.8,1,0.82,0.8,0.85,0.63,0.31,0.18,0.19,0.2,0.01],
    [0,0,0,0.36,0.1,0,0,0,0.05,0.54,0.86,0.79,0.74,0.72,0.6,0.39,0.28,0.24,0.13,0],
    [0,0,0,0.01,0.3,0.07,0,0,0.08,0.36,0.64,0.7,0.64,0.6,0.51,0.39,0.29,0.19,0.04,0],
    [0,0,0,0,0.1,0.24,0.14,0.1,0.15,0.29,0.45,0.53,0.52,0.46,0.4,0.31,0.21,0.08,0,0],
    [0,0,0,0,0,0.08,0.21,0.21,0.22,0.29,0.36,0.39,0.37,0.33,0.26,0.18,0.09,0,0,0],
    [0,0,0,0,0,0,0.03,0.13,0.19,0.22,0.24,0.24,0.23,0.18,0.13,0.05,0,0,0,0],
    [0,0,0,0,0,0,0,0,0.02,0.06,0.08,0.09,0.07,0.05,0.01,0,0,0,0,0]], dtype=float)


def _bell(x, m, s):
    return np.exp(-((x - m) / s) ** 2 / 2.0)


def _kernel_fft(N: int, R: int) -> np.ndarray:
    x = np.arange(N)
    X, Y = np.meshgrid(x, x)
    D = np.sqrt((X - N // 2) ** 2 + (Y - N // 2) ** 2) / R
    K = (D < 1) * _bell(D, 0.5, 0.15)
    K = K / K.sum()
    return np.fft.fft2(np.fft.fftshift(K))


def lenia(seed: int = 0, N: int = 64, R: int = 13, mu: float = 0.15, sigma: float = 0.015,
          T: int = 10, steps: int = 400, record: int = 80, n_creatures: int = 1,
          soup: bool = False) -> List[Dict[str, Any]]:
    rng = np.random.default_rng(seed)
    fK = _kernel_fft(N, R)
    A = np.zeros((N, N))
    if soup:
        # smooth random initialisation: low-pass filtered noise
        A = rng.random((N, N))
        A = np.real(np.fft.ifft2(np.fft.fft2(A) * fK)); A = (A - A.min()) / (np.ptp(A) + 1e-9)
    else:
        for _ in range(n_creatures):
            o = np.rot90(_ORBIUM, int(rng.integers(0, 4)))
            oy, ox = o.shape
            py, px = int(rng.integers(0, N - oy)), int(rng.integers(0, N - ox))
            A[py:py + oy, px:px + ox] = np.maximum(A[py:py + oy, px:px + ox], o)
    hist: List[Dict[str, Any]] = []
    rec_every = max(1, steps // record)
    for t in range(steps):
        U = np.real(np.fft.ifft2(fK * np.fft.fft2(A)))
        A = np.clip(A + (1.0 / T) * (_bell(U, mu, sigma) * 2 - 1), 0, 1)
        if t % rec_every == 0:
            hist.append({"field": A.copy(), "step": t})
    return hist
