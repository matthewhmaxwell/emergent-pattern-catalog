"""Active nematic model (P33) — lattice director field + lookalike generators.

Active nematics (Sanchez et al. 2012, Nature 491:431; Doostmohammadi et al. 2018,
Nat Commun 9:3246): apolar (head-tail-symmetric) rods build LOCAL nematic order
while continually nucleating/annihilating HALF-INTEGER (±1/2) topological defects —
"active turbulence". The ±1/2 defect is the topological fingerprint absent from any
polar (flocking) or integer-vortex (milling) system.

`active_nematic_field` is the canonical minimal model: a director theta∈[0,pi)
relaxes to the local nematic mean field with angular "activity" noise eta; above the
defect-unbinding point a steady-state ±1/2 defect gas persists. Self-propulsion is
apolar (random ± along the director) -> polar order φ≈0, local nematic order high.
The `neg_*` generators are the Phase-2a lookalikes for the negative panel.
"""
from __future__ import annotations

from typing import Any, Dict, List, Tuple

import numpy as np

PI = np.pi
Built = Tuple[List[Dict[str, Any]], Dict[str, Any]]


def active_nematic_field(seed: int = 0, G: int = 64, eta: float = 0.18,
                         n_steps: int = 400, n_frames: int = 60,
                         speed: float = 1.0, init_mode: str = "random") -> Built:
    """Canonical active nematic. eta≈0.18 = active-turbulent regime (locally ordered
    + sustained ±1/2 defect gas). init_mode='aligned' starts defect-free (the
    ordered-nematic negative)."""
    rng = np.random.default_rng(seed)
    if init_mode == "aligned":
        theta = (np.full((G, G), PI / 3.0) + 0.02 * rng.standard_normal((G, G))) % PI
    else:
        theta = rng.uniform(0, PI, size=(G, G))
    xs, ys = np.meshgrid(np.arange(G), np.arange(G), indexing="ij")
    pos = np.column_stack([xs.ravel(), ys.ravel()]).astype(float)
    frames: List[Dict[str, Any]] = []
    snap = max(1, n_steps // n_frames)
    for t in range(n_steps):
        m = (np.exp(2j * np.roll(theta, 1, 0)) + np.exp(2j * np.roll(theta, -1, 0))
             + np.exp(2j * np.roll(theta, 1, 1)) + np.exp(2j * np.roll(theta, -1, 1)))
        theta = (np.angle(m) / 2.0 + eta * rng.standard_normal((G, G))) % PI
        if t % snap == 0:
            ang = (theta.ravel() + rng.integers(0, 2, size=(G, G)).ravel() * PI)
            vel = speed * np.column_stack([np.cos(ang), np.sin(ang)])
            frames.append({"theta_field": theta.copy(), "velocities": vel,
                           "positions": pos.copy(), "headings": ang, "box_size": float(G)})
    return frames, {"model": "active_nematic_field", "eta": eta, "box_size": float(G),
                    "init_mode": init_mode}


def _emit(theta: np.ndarray, ang: np.ndarray, speed: float = 1.0) -> Dict[str, Any]:
    G = theta.shape[0]
    xs, ys = np.meshgrid(np.arange(G), np.arange(G), indexing="ij")
    pos = np.column_stack([xs.ravel(), ys.ravel()]).astype(float)
    vel = speed * np.column_stack([np.cos(ang.ravel()), np.sin(ang.ravel())])
    return {"theta_field": theta.copy(), "velocities": vel, "positions": pos,
            "headings": ang.ravel(), "box_size": float(G)}


def neg_polar_flock(seed: int = 0, G: int = 64, n_frames: int = 40) -> Built:
    """Polar flock: one global heading (small noise) -> φ high (rejects: polar)."""
    rng = np.random.default_rng(seed); base = rng.uniform(0, 2 * PI); frames = []
    for _ in range(n_frames):
        ang = (base + 0.15 * rng.standard_normal((G, G))).ravel() % (2 * PI)
        frames.append(_emit(ang.reshape(G, G) % PI, ang))
    return frames, {"model": "polar_flock", "box_size": float(G)}


def neg_milling(seed: int = 0, G: int = 64, n_frames: int = 40) -> Built:
    """Milling: a single +1 vortex director (tangential) -> |L| high, INTEGER charge."""
    rng = np.random.default_rng(seed)
    xs, ys = np.meshgrid(np.arange(G), np.arange(G), indexing="ij")
    c = (G - 1) / 2.0; phi = np.arctan2(ys - c, xs - c); frames = []
    for _ in range(n_frames):
        ang = (phi + PI / 2 + 0.12 * rng.standard_normal((G, G))).ravel() % (2 * PI)
        frames.append(_emit(ang.reshape(G, G) % PI, ang))
    return frames, {"model": "milling_vortex", "box_size": float(G)}


def neg_isotropic(seed: int = 0, G: int = 64, n_frames: int = 40) -> Built:
    """Isotropic noise: independent random orientations -> S_loc low (rejects)."""
    rng = np.random.default_rng(seed); frames = []
    for _ in range(n_frames):
        ang = rng.uniform(0, 2 * PI, size=(G, G)).ravel()
        frames.append(_emit(ang.reshape(G, G) % PI, ang))
    return frames, {"model": "isotropic_noise", "box_size": float(G)}


def neg_uniform_nematic(seed: int = 0, G: int = 64, n_frames: int = 40) -> Built:
    """Uniform nematic: aligned-init low-noise field -> ordered, ~0 defects (rejects:
    nematic order WITHOUT the active-turbulence topological signature)."""
    return active_nematic_field(seed, G=G, eta=0.05, n_steps=300, n_frames=n_frames,
                                init_mode="aligned")
