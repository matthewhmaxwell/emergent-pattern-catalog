"""Probe systems for the blind-spot / detection-recall audit.

Each probe is an unambiguously-emergent (or null) system spanning the six
emergence families, DELIBERATELY including morphologies the current generic
indicator is suspected to miss (rotation, connectivity, transient waves,
synergy, chaos, temporal oscillation). Each returns:

  history : list[frame]          for the current generic-emergence indicator
  micro   : ndarray [T, n]       per-part time series (for Psi_CE)
  macro   : ndarray [T]          a global summary series (for Psi_CE + MPR-C)
  truth   : "emergent" | "null"  ground truth
  family  : morphology/dynamics family label

Generators are minimal by design — the point is to test which detection CHANNEL
fires on which morphology, not to be canonical models.
"""
from __future__ import annotations

import math
from typing import Any, Dict, List

import numpy as np

from analysis import t2c_systems as T

Probe = Dict[str, Any]
TWO_PI = 2 * math.pi


def _torus_d(a, b, L):
    d = np.abs(a - b); return np.minimum(d, L - d)


# ---------------- positive controls (current instrument SHOULD catch) --------
def flocking(seed=0, N=80, L=8.0, steps=90, R=1.2, v=0.3, noise=0.25) -> Probe:
    rng = np.random.default_rng(seed)
    pos = rng.uniform(0, L, (N, 2)); th = rng.uniform(-math.pi, math.pi, N)
    H, heads, order = [], [], []
    for _ in range(steps):
        nth = th.copy()
        for i in range(N):
            dx = _torus_d(pos[i, 0], pos[:, 0], L); dy = _torus_d(pos[i, 1], pos[:, 1], L)
            nb = (dx * dx + dy * dy) <= R * R
            if nb.sum() > 1:
                mean = math.atan2(np.sin(th[nb]).mean(), np.cos(th[nb]).mean())
                nth[i] = mean + noise * rng.uniform(-math.pi, math.pi)
        th = (nth + math.pi) % TWO_PI - math.pi
        pos = (pos + v * np.c_[np.cos(th), np.sin(th)]) % L
        H.append({"positions": pos.copy(), "velocities": v * np.c_[np.cos(th), np.sin(th)]})
        heads.append(th.copy()); order.append(abs(np.exp(1j * th).mean()))
    return dict(history=H, micro=np.array(heads), macro=np.array(order),
                truth="emergent", family="alignment(control)")


def aggregation(seed=0, N=120, L=12.0, steps=90, v=0.25) -> Probe:
    rng = np.random.default_rng(seed)
    pos = rng.uniform(0, L, (N, 2)); H, dists, mnn = [], [], []
    for _ in range(steps):
        c = pos.mean(0)
        step = (c - pos)
        step /= (np.linalg.norm(step, axis=1, keepdims=True) + 1e-9)
        pos = pos + v * step + rng.normal(0, 0.05, pos.shape)
        H.append({"positions": pos.copy()})
        d = np.sqrt(((pos[:, None] - pos[None]) ** 2).sum(-1)); np.fill_diagonal(d, np.inf)
        dists.append(np.linalg.norm(pos - c, axis=1)); mnn.append(d.min(1).mean())
    return dict(history=H, micro=np.array(dists), macro=np.array(mnn),
                truth="emergent", family="clustering(control)")


def active_nematic(seed=0) -> Probe:
    h, _ = T.nov_active_nematic(seed)
    vel = [np.asarray(f["velocities"]) for f in h]
    th = [np.arctan2(v[:, 1], v[:, 0]) for v in vel]
    nem = np.array([abs(np.exp(2j * t).mean()) for t in th])
    return dict(history=h, micro=np.array(th), macro=nem,
                truth="emergent", family="orientational(control)")


# ---------------- suspected blind spots --------------------------------------
def vortex_milling(seed=0, N=100, L=20.0, steps=90, v=0.3) -> Probe:
    """Self-organizing mill: cohesion to centroid + tangential (curl) bias.
    Rotational order grows; clustering and polar/nematic order stay low."""
    rng = np.random.default_rng(seed)
    pos = rng.uniform(L * 0.3, L * 0.7, (N, 2)); th = rng.uniform(-math.pi, math.pi, N)
    H, ang, Lmom = [], [], []
    for _ in range(steps):
        c = pos.mean(0); rel = pos - c
        rad = np.arctan2(rel[:, 1], rel[:, 0])
        tang = rad + math.pi / 2                      # tangential heading
        coh = np.arctan2((c - pos)[:, 1], (c - pos)[:, 0])
        target = np.arctan2(0.7 * np.sin(tang) + 0.3 * np.sin(coh),
                            0.7 * np.cos(tang) + 0.3 * np.cos(coh))
        th = th + 0.5 * ((target - th + math.pi) % TWO_PI - math.pi) + rng.normal(0, 0.1, N)
        pos = pos + v * np.c_[np.cos(th), np.sin(th)]
        H.append({"positions": pos.copy(), "velocities": v * np.c_[np.cos(th), np.sin(th)]})
        rel = pos - pos.mean(0); r = np.linalg.norm(rel, axis=1) + 1e-9
        Lz = (rel[:, 0] * np.sin(th) - rel[:, 1] * np.cos(th)) / r   # per-agent ang. mom.
        ang.append(np.arctan2(rel[:, 1], rel[:, 0])); Lmom.append(abs(Lz.mean()))
    return dict(history=H, micro=np.array(ang), macro=np.array(Lmom),
                truth="emergent", family="rotation")


def lane_banding(seed=0, N=120, L=16.0, steps=90, v=0.3) -> Probe:
    """Counterflow: two groups moving +x / -x, repelling laterally → lanes."""
    rng = np.random.default_rng(seed)
    pos = rng.uniform(0, L, (N, 2)); grp = rng.integers(0, 2, N) * 2 - 1   # ±1
    H, lat, order = [], [], []
    for _ in range(steps):
        # lateral repulsion from opposite-group near neighbors
        dy = pos[:, 1][:, None] - pos[:, 1][None]
        dx = _torus_d(pos[:, 0][:, None], pos[:, 0][None], L)
        opp = (grp[:, None] != grp[None]) & (dx < 1.0) & (np.abs(dy) < 1.5)
        push = np.where(opp, np.sign(dy + 1e-9), 0).sum(1)
        pos[:, 1] = (pos[:, 1] + 0.15 * np.sign(push)) % L
        pos[:, 0] = (pos[:, 0] + v * grp) % L
        H.append({"positions": pos.copy(), "velocities": np.c_[v * grp, np.zeros(N)]})
        lat.append(pos[:, 1].copy())
        # lane order: variance of group-mean lateral positions (segregation)
        order.append(abs(pos[grp == 1, 1].mean() - pos[grp == -1, 1].mean()))
    return dict(history=H, micro=np.array(lat), macro=np.array(order),
                truth="emergent", family="banding")


def dla_fractal(seed=0) -> Probe:
    h, _ = T.nov_dla(seed)
    grids = [np.asarray(f["grid"], float) for f in h]
    flat = np.array([g.ravel() for g in grids])
    idx = np.linspace(0, flat.shape[1] - 1, 24).astype(int)
    macro = np.array([g.sum() for g in grids])           # cluster mass grows
    return dict(history=h, micro=flat[:, idx], macro=macro,
                truth="emergent", family="fractal")


def percolation(seed=0) -> Probe:
    from scipy import ndimage
    h, _ = T.probe_percolation(seed)
    grids = [np.asarray(f["grid"], float) for f in h]
    flat = np.array([g.ravel() for g in grids]); idx = np.linspace(0, flat.shape[1] - 1, 24).astype(int)
    giant = []
    for g in grids:
        lab, n = ndimage.label(g > 0.5)
        giant.append(int(np.bincount(lab.ravel())[1:].max()) if n else 0)
    return dict(history=h, micro=flat[:, idx], macro=np.array(giant, float),
                truth="emergent", family="connectivity")


def traveling_wave(seed=0, L=60, steps=110) -> Probe:
    """SIR front: susceptible→infected→recovered wave crosses, final state uniform."""
    rng = np.random.default_rng(seed)
    S, I, Rc = 0, 1, 2
    g = np.zeros((L, L), int); g[:, :2] = I               # front on the left
    H, micro_cells, active = [], [], []
    idx = rng.integers(0, L, (24, 2))
    for _ in range(steps):
        inf = (g == I)
        nb = (np.roll(inf, 1, 0) | np.roll(inf, -1, 0) | np.roll(inf, 1, 1) | np.roll(inf, -1, 1))
        newg = g.copy()
        newg[(g == S) & nb] = I
        newg[g == I] = Rc
        g = newg
        H.append({"grid": g.copy()})
        micro_cells.append(g[idx[:, 0], idx[:, 1]].astype(float))
        active.append((g == I).mean())
    return dict(history=H, micro=np.array(micro_cells), macro=np.array(active),
                truth="emergent", family="transient-wave")


def limit_cycle(seed=0, G=6, steps=240, D=0.12) -> Probe:
    """COUPLED Lotka-Volterra patches (nearest-neighbor diffusion) → the patches
    SYNCHRONIZE into a collective limit cycle. The coupling is what makes the
    global oscillation emergent (independent patches would average out); a
    coordination-gated detector must see this and reject uncoupled oscillators."""
    rng = np.random.default_rng(seed)
    prey = rng.uniform(0.3, 0.7, (G, G)); pred = rng.uniform(0.3, 0.7, (G, G))
    a, b, c, d, dt = 1.0, 0.5, 0.5, 0.5, 0.05

    def lap(f):
        return (np.roll(f, 1, 0) + np.roll(f, -1, 0) +
                np.roll(f, 1, 1) + np.roll(f, -1, 1) - 4 * f)

    H, micro, tot = [], [], []
    for _ in range(steps):
        prey = np.clip(prey + dt * (a * prey - b * prey * pred) + D * lap(prey), 0, 8)
        pred = np.clip(pred + dt * (-c * pred + d * prey * pred) + D * lap(pred), 0, 8)
        H.append({"field": prey.copy()})
        micro.append(prey.ravel().copy()); tot.append(float(prey.sum()))
    return dict(history=H, micro=np.array(micro), macro=np.array(tot),
                truth="emergent", family="temporal-oscillation")


def xor_synergy(seed=0, n=8, steps=600) -> Probe:
    """Parity flips deterministically; individual bits carry no info about it.
    Pure synergy — invisible to clustering/correlation/order-parameter channels."""
    rng = np.random.default_rng(seed)
    parity = np.zeros(steps, int); bits = np.zeros((steps, n), int)
    for t in range(steps):
        parity[t] = (1 - parity[t - 1]) if t > 0 else 0
        b = rng.integers(0, 2, n); b[-1] = parity[t] ^ (int(b[:-1].sum()) % 2)
        bits[t] = b
    H = [{"state": bits[t].astype(float)} for t in range(steps)]
    return dict(history=H, micro=bits.astype(float), macro=parity.astype(float),
                truth="emergent", family="synergy")


def spatiotemporal_chaos(seed=0, L=64, steps=200, eps=0.1, r=3.9) -> Probe:
    """Coupled logistic map lattice → turbulent chaos (low spatial autocorr,
    high dynamical complexity)."""
    rng = np.random.default_rng(seed)
    x = rng.uniform(0, 1, L)
    H, micro, mean = [], [], []
    idx = np.linspace(0, L - 1, 24).astype(int)
    for _ in range(steps):
        fx = r * x * (1 - x)
        x = (1 - eps) * fx + eps * 0.5 * (np.roll(fx, 1) + np.roll(fx, -1))
        H.append({"field": np.tile(x, (8, 1)).astype(float)})
        micro.append(x.copy()); mean.append(x.mean())
    return dict(history=H, micro=np.array(micro)[:, idx], macro=np.array(mean),
                truth="emergent", family="chaos")


# ---------------- nulls (instrument should NOT catch) ------------------------
def null_noise(seed=0) -> Probe:
    h, _ = T.null_spatial_noise(seed)
    pos = [np.asarray(f["positions"]) for f in h]
    micro = np.array([p[:24, 0] for p in pos]); macro = np.array([p[:, 0].mean() for p in pos])
    return dict(history=h, micro=micro, macro=macro, truth="null", family="null-noise")


def null_walk(seed=0) -> Probe:
    h, _ = T.null_random_walk(seed)
    pos = [np.asarray(f["positions"]) for f in h]
    micro = np.array([p[:24, 0] for p in pos]); macro = np.array([p[:, 0].mean() for p in pos])
    return dict(history=h, micro=micro, macro=macro, truth="null", family="null-walk")


PROBES = [
    flocking, aggregation, active_nematic,           # positive controls
    vortex_milling, lane_banding, dla_fractal, percolation,
    traveling_wave, limit_cycle, xor_synergy, spatiotemporal_chaos,
    null_noise, null_walk,                            # nulls
]
