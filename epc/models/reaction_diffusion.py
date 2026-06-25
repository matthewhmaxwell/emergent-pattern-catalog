"""Gray-Scott reaction-diffusion — a multi-family substrate for exercising the lens
battery out-of-distribution.

    du/dt = Du ∇²u - u v² + F(1-u)
    dv/dt = Dv ∇²v + u v² - (F+k) v

A continuous FIELD substrate (not agents, not grid-CA), so it stresses the field lenses
(structure_factor, fractal_dimension) and the generic emergence indicator on inputs
unlike anything in the catalog's agent/network corpus. Sweeping (F, k) walks the classic
Pearson atlas: solitons, self-replicating "mitosis" chaos, stripes/labyrinths, spots,
worms, and (at extreme params) uniform death. Records the v FIELD only — RD is a field
substrate and should advertise ONE primary representation. (Hardening lesson: when a
frame carries both 'field' and 'positions', field lenses dispatch on 'positions' first
and read the weaker view; positions lenses on RD want the field's sublevel-set loops,
not spot centroids — that field-PH path is deferred, so centroids are not emitted.)
"""
from __future__ import annotations

from typing import Any, Dict, List

import numpy as np


def _lap(a: np.ndarray) -> np.ndarray:
    return (-4.0 * a + np.roll(a, 1, 0) + np.roll(a, -1, 0)
            + np.roll(a, 1, 1) + np.roll(a, -1, 1))


def reaction_diffusion(seed: int = 0, F: float = 0.030, k: float = 0.062, N: int = 96,
                       steps: int = 8000, Du: float = 0.16, Dv: float = 0.08,
                       record: int = 24, dt: float = 1.0) -> List[Dict[str, Any]]:
    rng = np.random.default_rng(seed)
    u = np.ones((N, N)); v = np.zeros((N, N))
    r, c = N // 8, N // 2
    u[c - r:c + r, c - r:c + r] = 0.50
    v[c - r:c + r, c - r:c + r] = 0.25
    u = np.clip(u + 0.02 * rng.standard_normal((N, N)), 0, 1)
    v = np.clip(v + 0.02 * rng.standard_normal((N, N)), 0, 1)
    hist: List[Dict[str, Any]] = []
    rec_every = max(1, steps // record)
    for t in range(steps):
        uvv = u * v * v
        u += dt * (Du * _lap(u) - uvv + F * (1.0 - u))
        v += dt * (Dv * _lap(v) + uvv - (F + k) * v)
        if t % rec_every == 0:
            hist.append({"field": v.copy(), "step": t})
    return hist
