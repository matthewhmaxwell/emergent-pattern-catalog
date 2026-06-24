"""Evolutionary games on networks — network reciprocity (Pattern P38 candidate).

Ring-1 combination cell: substrate = fixed_net AND update = reproduction.

A weak Prisoner's Dilemma is played on a FIXED network; nodes reproduce their
strategy by IMITATING more-successful neighbours (asynchronous proportional
imitation, Santos-Pacheco 2005). Each node accumulates payoff over ALL neighbours,
so high-degree hubs that cooperate become very successful and are widely copied.
On a heterogeneous (scale-free) network this makes cooperation survive far above
the well-mixed limit AND concentrate on the hubs: cooperation frequency
correlates with node degree. That degree-cooperation correlation is the signature
of NETWORK reciprocity, structurally impossible on a uniform-degree lattice
(P27 spatial reciprocity) — the key separation from the catalog.

Game (weak PD): payoff to i = sum over neighbours of M[s_i, s_j], with
    R(C,C)=1, S(C,D)=0, T(D,C)=b, P(D,D)=0,   1 < b < 2 (temptation).
Update (async proportional imitation): focal x picks random neighbour y; if
    pi_y > pi_x, x copies s_y with probability (pi_y - pi_x)/(b*max(k_x,k_y)).

References
----------
Santos, F.C. & Pacheco, J.M. (2005). Scale-free networks provide a unifying
    framework for the emergence of cooperation. Phys. Rev. Lett. 95, 098104.
Ohtsuki, H., Hauert, C., Lieberman, E. & Nowak, M.A. (2006). A simple rule for
    the evolution of cooperation on graphs and social networks. Nature 441, 502.
"""
from __future__ import annotations

from typing import Any, Dict, List, Tuple

import numpy as np

Built = Tuple[List[Dict[str, Any]], Dict[str, Any]]


def _ba_network(N: int, m: int, rng: np.random.Generator):
    """Barabasi-Albert preferential-attachment scale-free graph -> CSR."""
    adj = [set(j for j in range(m) if j != i) for i in range(m)]
    repeated = list(range(m)) * m
    for new in range(m, N):
        adj.append(set())
        chosen = set()
        while len(chosen) < m:
            t = repeated[rng.integers(len(repeated))]
            if t != new:
                chosen.add(t)
        for t in chosen:
            adj[new].add(t); adj[t].add(new)
        repeated.extend(chosen); repeated.extend([new] * m)
    return _csr(adj, N)


def _lattice_network(L: int):
    """Periodic 2D von-Neumann lattice (uniform degree 4) -> CSR."""
    N = L * L
    adj = [[] for _ in range(N)]
    for x in range(L):
        for y in range(L):
            i = x * L + y
            for dx, dy in ((1, 0), (-1, 0), (0, 1), (0, -1)):
                adj[i].append(((x + dx) % L) * L + (y + dy) % L)
    return _csr([set(a) for a in adj], N)


def _csr(adj, N):
    cnt = np.array([len(a) for a in adj])
    start = np.zeros(N, dtype=int); start[1:] = np.cumsum(cnt)[:-1]
    flat = np.concatenate([np.fromiter(a, dtype=int, count=len(a)) for a in adj])
    return flat, start, cnt, cnt.astype(float)


def _async_game(flat, start, cnt, deg, s, b, generations, n_frames, rng):
    """Asynchronous proportional-imitation PD. s is a bool array (True=cooperate)."""
    N = len(s)
    coop_nbr = np.empty(N)
    for i in range(N):
        coop_nbr[i] = s[flat[start[i]:start[i] + cnt[i]]].sum()
    snap = max(1, generations // n_frames)
    frames = [{"strategies": s.astype(float).copy(), "degree": deg.copy(), "step": 0}]
    for g in range(1, generations + 1):
        xs = rng.integers(0, N, size=N)
        rn = rng.random(N)
        rp = rng.random(N)
        for t in range(N):
            x = int(xs[t]); kx = int(cnt[x])
            y = int(flat[start[x] + int(rn[t] * kx)])
            px = coop_nbr[x] * (1.0 if s[x] else b)
            py = coop_nbr[y] * (1.0 if s[y] else b)
            if py > px and s[x] != s[y]:
                if rp[t] < (py - px) / (b * max(deg[x], deg[y])):
                    s[x] = s[y]
                    nb = flat[start[x]:start[x] + kx]
                    coop_nbr[nb] += (1.0 if s[x] else -1.0)
        if g % snap == 0:
            frames.append({"strategies": s.astype(float).copy(), "degree": deg.copy(), "step": g})
    return frames


def network_reciprocity(seed: int = 0, N: int = 500, m: int = 2, b: float = 1.7,
                        generations: int = 2000, n_frames: int = 40) -> Built:
    """Positive: weak PD on a scale-free network -> cooperation survives and
    concentrates on hubs (degree-cooperation correlation)."""
    rng = np.random.default_rng(seed * 7919 + 1)
    flat, start, cnt, deg = _ba_network(N, m, rng)
    s = rng.random(N) < 0.5
    frames = _async_game(flat, start, cnt, deg, s, b, generations, n_frames, rng)
    return frames, {"model": "network_reciprocity", "N": N, "mean_degree": float(deg.mean()),
                    "scale_free": True, "temptation_b": b}


def neg_lattice_reciprocity(seed: int = 0, L: int = 22, b: float = 1.7,
                            generations: int = 2000, n_frames: int = 40) -> Built:
    """Same game on a uniform-degree LATTICE (P27 spatial reciprocity): cooperation
    may persist in clusters, but degree is constant -> NO degree-cooperation correlation."""
    rng = np.random.default_rng(seed * 104729 + 3)
    flat, start, cnt, deg = _lattice_network(L)
    s = rng.random(L * L) < 0.5
    frames = _async_game(flat, start, cnt, deg, s, b, generations, n_frames, rng)
    return frames, {"model": "lattice_reciprocity", "N": L * L, "mean_degree": 4.0, "scale_free": False}


def neg_wellmixed_defection(seed: int = 0, N: int = 500, k: int = 40, b: float = 1.7,
                            generations: int = 2000, n_frames: int = 40) -> Built:
    """High-degree random-regular (well-mixed) game -> tragedy of the commons: defection."""
    rng = np.random.default_rng(seed * 13 + 5)
    adj = [set() for _ in range(N)]
    for i in range(N):
        while len(adj[i]) < k:
            j = int(rng.integers(N))
            if j != i:
                adj[i].add(j); adj[j].add(i)
    flat, start, cnt, deg = _csr(adj, N)
    s = rng.random(N) < 0.5
    frames = _async_game(flat, start, cnt, deg, s, b, generations, n_frames, rng)
    return frames, {"model": "wellmixed_defection", "N": N, "mean_degree": float(deg.mean()), "scale_free": False}


def neg_scalefree_random(seed: int = 0, N: int = 500, m: int = 2, n_frames: int = 40) -> Built:
    """Scale-free network, strategies RE-RANDOMISED each frame (no imitation) ->
    no degree-cooperation correlation builds up."""
    rng = np.random.default_rng(seed * 31 + 7)
    flat, start, cnt, deg = _ba_network(N, m, rng)
    frames = [{"strategies": (rng.random(N) < 0.5).astype(float), "degree": deg.copy(), "step": fr}
              for fr in range(n_frames)]
    return frames, {"model": "scalefree_random_strategies", "N": N,
                    "mean_degree": float(deg.mean()), "scale_free": True}
