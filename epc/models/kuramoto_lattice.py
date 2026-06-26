"""Kuramoto oscillators on a 2D lattice — a multi-family substrate for exercising the
lens battery out-of-distribution, especially directed_info_flow (causal direction) and
the spectral/temporal channel.

    dtheta_i/dt = omega_i + K * sum_{j in 4-neighbourhood} sin(theta_j - theta_i)

Regimes (set by coupling K, natural-frequency spread, and initial condition):
  - incoherent : weak K + heterogeneous omega -> no coupling structure (null-ish).
  - sync       : strong K + identical omega -> global synchrony (symmetric coupling).
  - spiral     : identical omega + a topological spiral seed -> a sustained rotating
                 wave (the winding is topologically protected, so it persists).
  - plane      : identical omega + a phase-gradient seed -> a sustained travelling wave
                 with a clear propagation DIRECTION.

Emits per-frame 'phases' (the N*N grid), so micro_macro / the emergence indicator read
it directly; build sin(phases) component series for directed_info_flow.
"""
from __future__ import annotations

from typing import Any, Dict, List

import numpy as np


def kuramoto_lattice(seed: int = 0, N: int = 40, K: float = 1.5, omega0: float = 1.0,
                     omega_spread: float = 0.0, steps: int = 3000, dt: float = 0.05,
                     record: int = 120, init: str = "random") -> List[Dict[str, Any]]:
    rng = np.random.default_rng(seed)
    if init == "spiral":
        y, x = np.indices((N, N))
        theta = np.arctan2(y - N / 2.0, x - N / 2.0) % (2 * np.pi)
    elif init == "plane":
        y, x = np.indices((N, N))
        theta = (2 * np.pi * x / (N / 3.0)) % (2 * np.pi)
    else:
        theta = rng.uniform(0, 2 * np.pi, (N, N))
    omega = omega0 + omega_spread * rng.standard_normal((N, N))
    hist: List[Dict[str, Any]] = []
    rec_every = max(1, steps // record)
    for t in range(steps):
        s = (np.sin(np.roll(theta, 1, 0) - theta) + np.sin(np.roll(theta, -1, 0) - theta)
             + np.sin(np.roll(theta, 1, 1) - theta) + np.sin(np.roll(theta, -1, 1) - theta))
        theta = np.mod(theta + dt * (omega + K * s), 2 * np.pi)
        if t % rec_every == 0:
            hist.append({"phases": theta.copy(), "step": t})
    return hist


def order_parameter(history: List[Dict[str, Any]]) -> float:
    """Kuramoto global coherence r in [0,1] on the final frame (sanity check)."""
    th = np.asarray(history[-1]["phases"], dtype=float).ravel()
    return float(np.abs(np.exp(1j * th).mean()))
