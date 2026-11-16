"""Class A synthetic substrates for the discovery-phase detector formats (P33-P37).

The ten Class A generators in :mod:`epc.phase2a.synthetic` dispatch here for the
four detector formats introduced with the P33-P37 full panels (validation
rebuild, panels v12). P36 needs no new format: it reads the existing
``scalar_wealth`` format (P28's), whose Class A branches are reused unchanged.

Every substrate below is STRUCTURELESS by construction: i.i.d. draws, independent
per-unit random walks, or fixed deterministic textures (constant / gradient /
checkerboard) with no dynamics. None contains the mechanism of any catalog
pattern. They are written in exactly the frame schema the canonical positive of
the corresponding pattern uses, so the detector reads them through the same code
path as its positives.

Formats
-------
``orientation_field`` (P33)
    Frames ``{theta_field (G,G) in [0,pi), velocities (G*G,2), positions (G*G,2)
    lattice sites, headings (G*G,), box_size}`` -- the schema of
    ``epc.models.active_nematic.active_nematic_field`` (G=64). For the random
    kinds the per-site heading is the director plus a random head/tail flip
    (apolar emission, as in the positive); for the deterministic kinds the
    heading equals the director.
``adjacency_timeseries`` (P34)
    Frames ``{adjacency (N,N) int8 symmetric, opinions (N,) in {0,1}}`` -- the
    schema of ``epc.models.adaptive_network.adaptive_voter`` (N=200, mean degree
    4 where a degree is chosen).
``positions`` (P35)
    Frames ``{positions (N,2) in [0,L)^2, box_size L}`` -- the schema of
    ``epc.models.entropy_crystal.hard_disk_crystallization`` (N=196, L=14), plus
    the per-frame ``box_size`` so the P35 panel scores every substrate with its
    own periodic box.
``abundance_resource_timeseries`` (P37)
    Frames ``{abundances (5,), resources (5,), step}`` -- the schema of
    ``epc.models.resource_competition.resource_competition`` (5 species, 5
    resources).

Kinds (one per Class A generator id)
------------------------------------
random_uniform, random_gaussian, random_binary, spatial_white_noise,
temporal_white_noise, constant, linear_gradient, checkerboard. The two
positive-derived generators (permutation_shuffled, time_shuffled) are handled in
:mod:`epc.phase2a.synthetic`, see :func:`permutation_shuffled_discovery`.
"""
from __future__ import annotations

from typing import Any, Dict, List

import numpy as np

PI = np.pi

FORMATS = ("orientation_field", "adjacency_timeseries", "positions",
           "abundance_resource_timeseries")

# Geometry / size defaults = the canonical positives' sizes.
OF_G = 64                 # P33 lattice side
ADJ_N = 200               # P34 nodes
ADJ_MEAN_DEG = 4.0        # P34 mean degree
POS_N = 196               # P35 particles
POS_L = 14.0              # P35 periodic box
AR_SPECIES = 5            # P37 species
AR_RESOURCES = 5          # P37 resources
AR_SUPPLY = np.array([6.0, 10.0, 14.0, 4.0, 9.0])   # P37 positive's supply levels (scale only)


# --- orientation_field (P33) -------------------------------------------------

def _of_lattice(G: int) -> np.ndarray:
    xs, ys = np.meshgrid(np.arange(G), np.arange(G), indexing="ij")
    return np.column_stack([xs.ravel(), ys.ravel()]).astype(float)


def _of_frame(theta: np.ndarray, ang: np.ndarray, pos: np.ndarray, t: int) -> Dict[str, Any]:
    a = np.asarray(ang, dtype=float).ravel()
    return {"theta_field": np.asarray(theta, dtype=float).copy(),
            "velocities": np.column_stack([np.cos(a), np.sin(a)]),
            "positions": pos, "headings": a, "box_size": float(theta.shape[0]), "step": t}


def _orientation_field(rng: np.random.Generator, kind: str, n_steps: int,
                       G: int = OF_G) -> List[Dict[str, Any]]:
    pos = _of_lattice(G)
    out: List[Dict[str, Any]] = []

    def apolar(theta):  # heading = director +- (random head/tail), as the positive emits
        return theta.ravel() + PI * rng.integers(0, 2, size=G * G)

    if kind in ("random_uniform", "spatial_white_noise"):
        for t in range(n_steps):
            th = rng.uniform(0.0, PI, size=(G, G))
            out.append(_of_frame(th, apolar(th), pos, t))
    elif kind == "random_gaussian":
        # i.i.d. Gaussian angle, std = half the nematic range (pi/2), wrapped to [0,pi)
        # -- the orientation analog of the phases-format branch (std pi on [0,2pi)).
        for t in range(n_steps):
            th = (rng.standard_normal((G, G)) * (PI / 2.0)) % PI
            out.append(_of_frame(th, apolar(th), pos, t))
    elif kind == "random_binary":
        # Bernoulli(0.5) between two orthogonal directors {0, pi/2}.
        for t in range(n_steps):
            th = rng.integers(0, 2, size=(G, G)) * (PI / 2.0)
            out.append(_of_frame(th, apolar(th), pos, t))
    elif kind == "temporal_white_noise":
        # Each site's director performs an independent random walk (no coupling).
        th = rng.uniform(0.0, PI, size=(G, G))
        flip = PI * rng.integers(0, 2, size=G * G)
        for t in range(n_steps):
            th = (th + rng.normal(0.0, 0.3, size=(G, G))) % PI
            out.append(_of_frame(th, th.ravel() + flip, pos, t))
    elif kind == "constant":
        th = np.zeros((G, G))
        for t in range(n_steps):
            out.append(_of_frame(th, th, pos, t))
    elif kind == "linear_gradient":
        # Director rotates smoothly by pi across the box (defect-free texture), static.
        th = np.broadcast_to((PI * np.arange(G) / G)[:, None], (G, G)).copy()
        for t in range(n_steps):
            out.append(_of_frame(th, th, pos, t))
    elif kind == "checkerboard":
        i, j = np.meshgrid(np.arange(G), np.arange(G), indexing="ij")
        th = ((i + j) % 2) * (PI / 2.0)
        for t in range(n_steps):
            out.append(_of_frame(th, th, pos, t))
    else:
        raise ValueError(f"unknown kind: {kind}")
    return out


# --- adjacency_timeseries (P34) ----------------------------------------------

def _er(rng: np.random.Generator, n: int, p: float) -> np.ndarray:
    u = np.triu(rng.random((n, n)) < p, 1)
    return (u | u.T).astype(np.int8)


def _ring(n: int, offsets) -> np.ndarray:
    A = np.zeros((n, n), dtype=np.int8)
    idx = np.arange(n)
    for o in offsets:
        A[idx, (idx + o) % n] = 1
        A[(idx + o) % n, idx] = 1
    np.fill_diagonal(A, 0)
    return A


def _adjacency(rng: np.random.Generator, kind: str, n_steps: int,
               N: int = ADJ_N) -> List[Dict[str, Any]]:
    p = ADJ_MEAN_DEG / (N - 1)
    out: List[Dict[str, Any]] = []
    if kind in ("random_uniform", "spatial_white_noise"):
        # Fresh Erdos-Renyi graph (mean degree 4) + fresh random opinions each frame.
        for t in range(n_steps):
            out.append({"adjacency": _er(rng, N, p), "opinions": rng.integers(0, 2, N), "step": t})
    elif kind == "random_gaussian":
        # Fresh Chung-Lu graph with Gaussian expected degrees |N(4, 2)|: degree
        # heterogeneity but no community structure.
        for t in range(n_steps):
            w = np.abs(rng.normal(ADJ_MEAN_DEG, 2.0, size=N)) + 1e-3
            P = np.minimum(np.outer(w, w) / w.sum(), 1.0)
            u = np.triu(rng.random((N, N)) < P, 1)
            out.append({"adjacency": (u | u.T).astype(np.int8),
                        "opinions": rng.integers(0, 2, N), "step": t})
    elif kind == "random_binary":
        # Every node pair is an edge with probability 0.5 (dense Bernoulli(0.5) graph).
        for t in range(n_steps):
            out.append({"adjacency": _er(rng, N, 0.5), "opinions": rng.integers(0, 2, N), "step": t})
    elif kind == "temporal_white_noise":
        # Edge-flip random walk at stationary density p: every present edge is
        # removed w.p. q and every absent pair added w.p. q*p/(1-p) per frame;
        # every opinion flips w.p. q. No coupling between topology and opinions.
        q = 0.05
        A = _er(rng, N, p).astype(bool)
        op = rng.integers(0, 2, N)
        for t in range(n_steps):
            r = np.triu(rng.random((N, N)), 1)
            up = np.triu(A, 1)
            new = np.where(up, r >= q, r < q * p / (1.0 - p))
            new = np.triu(new, 1)
            A = new | new.T
            op = np.where(rng.random(N) < q, 1 - op, op)
            out.append({"adjacency": A.astype(np.int8), "opinions": op.copy(), "step": t})
    elif kind == "constant":
        A = np.ones((N, N), dtype=np.int8)
        np.fill_diagonal(A, 0)              # complete graph, static, all one opinion
        op = np.zeros(N, dtype=int)
        out = [{"adjacency": A, "opinions": op, "step": t} for t in range(n_steps)]
    elif kind == "linear_gradient":
        # Static ring lattice (each node linked to its 2 nearest on each side,
        # degree 4) with an opinion gradient along the ring (first half 0, rest 1).
        A = _ring(N, (1, 2))
        op = (np.arange(N) >= N // 2).astype(int)
        out = [{"adjacency": A, "opinions": op, "step": t} for t in range(n_steps)]
    elif kind == "checkerboard":
        # Static bipartite ring (odd offsets 1 and 3, degree 4) with alternating opinions.
        A = _ring(N, (1, 3))
        op = (np.arange(N) % 2).astype(int)
        out = [{"adjacency": A, "opinions": op, "step": t} for t in range(n_steps)]
    else:
        raise ValueError(f"unknown kind: {kind}")
    return out


# --- positions (P35) -----------------------------------------------------------

def _positions(rng: np.random.Generator, kind: str, n_steps: int,
               N: int = POS_N, L: float = POS_L) -> List[Dict[str, Any]]:
    out: List[Dict[str, Any]] = []
    if kind in ("random_uniform", "spatial_white_noise"):
        for t in range(n_steps):        # ideal gas re-drawn every frame
            out.append({"positions": rng.uniform(0.0, L, size=(N, 2)), "box_size": L, "step": t})
    elif kind == "random_gaussian":
        for t in range(n_steps):        # Gaussian blob (std L/4) re-drawn every frame
            out.append({"positions": (L / 2.0 + rng.normal(0.0, L / 4.0, size=(N, 2))) % L,
                        "box_size": L, "step": t})
    elif kind == "random_binary":
        # Random half-occupation of a square lattice with 2N-ish sites (20x20),
        # re-drawn every frame.
        m = 20
        a = L / m
        xs, ys = np.meshgrid(np.arange(m), np.arange(m), indexing="ij")
        sites = (np.column_stack([xs.ravel(), ys.ravel()]) + 0.5) * a
        for t in range(n_steps):
            idx = rng.choice(len(sites), size=N, replace=False)
            out.append({"positions": sites[idx].copy(), "box_size": L, "step": t})
    elif kind == "temporal_white_noise":
        pos = rng.uniform(0.0, L, size=(N, 2))     # independent Brownian walkers
        for t in range(n_steps):
            pos = (pos + rng.normal(0.0, 0.1, size=(N, 2))) % L
            out.append({"positions": pos.copy(), "box_size": L, "step": t})
    elif kind == "constant":
        pos = rng.uniform(0.0, L, size=(N, 2))     # one random configuration, frozen
        out = [{"positions": pos, "box_size": L, "step": t} for t in range(n_steps)]
    elif kind == "linear_gradient":
        x = L * np.sqrt(rng.random(N))             # density rising linearly in x, frozen
        pos = np.column_stack([x, rng.uniform(0.0, L, N)]) % L
        out = [{"positions": pos, "box_size": L, "step": t} for t in range(n_steps)]
    elif kind == "checkerboard":
        m = int(round(np.sqrt(N)))                 # perfect square lattice 14x14, frozen
        a = L / m
        xs, ys = np.meshgrid(np.arange(m), np.arange(m), indexing="ij")
        pos = (np.column_stack([xs.ravel(), ys.ravel()]) + 0.5) * a
        out = [{"positions": pos, "box_size": L, "step": t} for t in range(n_steps)]
    else:
        raise ValueError(f"unknown kind: {kind}")
    return out


# --- abundance_resource_timeseries (P37) --------------------------------------

def _abundance_resource(rng: np.random.Generator, kind: str, n_steps: int,
                        n: int = AR_SPECIES, k: int = AR_RESOURCES) -> List[Dict[str, Any]]:
    S = AR_SUPPLY[:k] if k <= len(AR_SUPPLY) else np.full(k, 8.0)
    out: List[Dict[str, Any]] = []

    def fr(a, r, t):
        return {"abundances": np.asarray(a, dtype=float), "resources": np.asarray(r, dtype=float), "step": t}

    if kind == "random_uniform":
        for t in range(n_steps):
            out.append(fr(rng.uniform(0.0, 1.0, n), rng.uniform(0.0, 1.0, k) * S, t))
    elif kind == "random_gaussian":
        for t in range(n_steps):
            out.append(fr(np.clip(rng.normal(1.0, 0.3, n), 0.0, None),
                          np.clip(rng.normal(1.0, 0.3, k), 0.0, None) * S, t))
    elif kind == "random_binary":
        for t in range(n_steps):
            out.append(fr(rng.choice([0.5, 1.5], size=n), rng.choice([0.5, 1.5], size=k) * S, t))
    elif kind == "spatial_white_noise":
        for t in range(n_steps):        # fresh i.i.d. exponential draws every frame
            out.append(fr(rng.exponential(1.0, n), rng.exponential(1.0, k) * S, t))
    elif kind == "temporal_white_noise":
        la = np.log(rng.uniform(0.5, 1.5, n))     # independent log-random walks
        lr = np.log(rng.uniform(0.5, 1.5, k))
        for t in range(n_steps):
            la = la + rng.normal(0.0, 0.05, n)
            lr = lr + rng.normal(0.0, 0.05, k)
            out.append(fr(np.exp(la), np.exp(lr) * S, t))
    elif kind == "constant":
        for t in range(n_steps):
            out.append(fr(np.ones(n), S.copy(), t))
    elif kind == "linear_gradient":
        # Monotonic drift, no oscillation: each abundance grows linearly in time
        # (species-specific scale), each resource is drawn down linearly.
        scale = 0.5 + 0.25 * np.arange(n)
        for t in range(n_steps):
            f = (t + 1) / n_steps
            out.append(fr(scale * f, S * (1.0 - 0.5 * f), t))
    elif kind == "checkerboard":
        # Deterministic period-2 alternation (neighbouring species in antiphase).
        for t in range(n_steps):
            sgn_a = (-1.0) ** (t + np.arange(n))
            sgn_r = (-1.0) ** (t + np.arange(k))
            out.append(fr(1.0 + 0.5 * sgn_a, S * (1.0 + 0.5 * sgn_r), t))
    else:
        raise ValueError(f"unknown kind: {kind}")
    return out


_DISPATCH = {
    "orientation_field": _orientation_field,
    "adjacency_timeseries": _adjacency,
    "positions": _positions,
    "abundance_resource_timeseries": _abundance_resource,
}


def null_history(format: str, kind: str, seed: int, n_steps: int = 200) -> List[Dict[str, Any]]:
    """Class A substrate ``kind`` in discovery ``format`` (seeded, deterministic)."""
    rng = np.random.default_rng(seed)
    return _DISPATCH[format](rng, kind, int(n_steps))


def permutation_shuffled_discovery(format: str, seed: int,
                                   positive: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Class A ``permutation_shuffled`` for the discovery formats.

    orientation_field     the positive's FINAL frame with its lattice sites
                          permuted (director, velocity and heading permuted
                          together; site coordinates kept), held static for the
                          positive's length -- destroys spatial structure,
                          preserves the marginal orientation distribution (the
                          grid-format convention).
    adjacency_timeseries  the positive's FINAL graph with its edges re-placed at
                          random (upper-triangle entries permuted: same edge
                          count, no community structure), held static -- the
                          trail_network convention (weights permuted across edges).
    positions             consistent relabelling of particle indices across all
                          frames (psi6 is invariant to it: degenerate, skipped
                          when the detector declares permutation invariance).
    abundance_resource_timeseries
                          consistent relabelling of species (and of resources)
                          across all frames (per-species statistics are invariant
                          to it: degenerate, skipped likewise).
    """
    rng = np.random.default_rng(seed)
    T = len(positive)
    if format == "orientation_field":
        last = positive[-1]
        th = np.asarray(last["theta_field"], dtype=float)
        G = th.shape[0]
        perm = rng.permutation(th.size)
        th_s = th.ravel()[perm].reshape(th.shape)
        hd = np.asarray(last["headings"], dtype=float).ravel()[perm]
        pos = np.asarray(last["positions"], dtype=float)
        frame = {"theta_field": th_s, "velocities": np.column_stack([np.cos(hd), np.sin(hd)]),
                 "positions": pos, "headings": hd, "box_size": float(last.get("box_size", G))}
        return [dict(frame, step=t) for t in range(T)]
    if format == "adjacency_timeseries":
        last = positive[-1]
        A = (np.asarray(last["adjacency"]) > 0)
        n = A.shape[0]
        iu = np.triu_indices(n, 1)
        vals = A[iu].copy()
        rng.shuffle(vals)
        B = np.zeros((n, n), dtype=np.int8)
        B[iu] = vals
        B = B + B.T
        frame = {"adjacency": B}
        if "opinions" in last:
            frame["opinions"] = np.asarray(last["opinions"]).copy()
        return [dict(frame, step=t) for t in range(T)]
    if format == "positions":
        n = len(positive[0]["positions"])
        perm = rng.permutation(n)
        return [dict(f, positions=np.asarray(f["positions"])[perm], step=t)
                for t, f in enumerate(positive)]
    if format == "abundance_resource_timeseries":
        n = len(positive[0]["abundances"])
        k = len(positive[0]["resources"])
        ps, pr = rng.permutation(n), rng.permutation(k)
        return [dict(f, abundances=np.asarray(f["abundances"])[ps],
                     resources=np.asarray(f["resources"])[pr], step=t)
                for t, f in enumerate(positive)]
    raise ValueError(f"unknown format: {format}")
