"""Voter model — canonical model for P18 (consensus).

Reference:
    Clifford, P. & Sudbury, A. (1973). "A model for spatial conflict."
    Biometrika, 60(3), 581-588. DOI: 10.1093/biomet/60.3.581

    Holley, R. & Liggett, T.M. (1975). "Ergodic theorems for weakly
    interacting infinite systems and the voter model."
    The Annals of Probability, 3(4), 643-663.

    Dornic, I., Chaté, H., Chave, J., & Hinrichsen, H. (2001). "Critical
    coarsening without surface tension: the universality class of the
    voter model." Phys. Rev. Lett. 87(4), 045701.
    DOI: 10.1103/PhysRevLett.87.045701

Update rule (asynchronous):
    One "step" = N site-updates (one Monte Carlo sweep), where N = rows * cols.
    Per site-update: pick a random site uniformly, pick a random Moore-neighbor
    uniformly, copy the neighbor's state onto the site.

    This is the canonical asynchronous Glauber-like dynamics. Synchronous
    updating introduces stripe artifacts and is not used here.

Key known results (2D, Moore neighborhood, binary states):
    - Coarsening without surface tension: domain-wall density ρ_w(t) ∝ ln(t)/t
      (logarithmic correction; Dornic et al. 2001). Often empirically fit as
      a power law with exponent ~-1/2 over accessible simulation times.
    - On the infinite 2D lattice: no true consensus (recurrent but does not
      absorb). On finite L×L torus: consensus reached in mean time
      τ_c(L) ∝ L² ln(L) for 2D.
    - No surface tension → no Lifshitz-Slyozov L(t) ∝ t^{1/2} scaling as in
      Ising; voter model is in a distinct universality class.

State history keys:
    grid       : np.ndarray (rows, cols), dtype int8 — opinion {0, 1}
    grid_dims  : tuple (rows, cols)
    step       : int — sweep number (1 sweep = N site-updates)
    magnetization : float — m = (2 * frac_1) - 1 ∈ [-1, 1]
    abs_magnetization : float — |m|
    wall_density : float — fraction of Moore-neighbor pairs with differing states
    moran_i    : float — Moran's I of opinion indicator
    consensus_reached : bool — True iff all sites share one opinion
"""

from __future__ import annotations

from typing import Any

import numpy as np


class VoterModel:
    """2D voter model on a torus (Moore neighborhood, asynchronous updating).

    Parameters
    ----------
    rows, cols : int
        Grid dimensions. Default 64×64.
    init_mode : str
        'random' (50/50), 'biased' (fraction p of opinion 1), 'half_and_half'
        (left half 0, right half 1).
    init_fraction : float
        For 'biased' mode: initial fraction of opinion 1.
    neighborhood : str
        'moore' (8 neighbors) or 'von_neumann' (4 neighbors). Default 'moore'.
    boundary : str
        'periodic' (torus) only supported for now.
    seed : int
        Random seed.
    """

    def __init__(
        self,
        rows: int = 64,
        cols: int = 64,
        init_mode: str = "random",
        init_fraction: float = 0.5,
        neighborhood: str = "moore",
        boundary: str = "periodic",
        seed: int = 42,
    ):
        self.rows = rows
        self.cols = cols
        self.init_mode = init_mode
        self.init_fraction = init_fraction
        self.neighborhood = neighborhood
        self.boundary = boundary
        self.seed = seed

        if neighborhood not in ("moore", "von_neumann"):
            raise ValueError(f"Unknown neighborhood: {neighborhood}")
        if boundary != "periodic":
            raise ValueError(
                f"Only 'periodic' boundary supported, got {boundary}"
            )

        # State
        self.grid: np.ndarray = np.empty(0, dtype=np.int8)
        self.rng: np.random.Generator = np.random.default_rng(seed)
        self._step_count: int = 0

        # Precompute neighbor offsets
        if neighborhood == "moore":
            self._n_offsets = [
                (-1, -1), (-1, 0), (-1, 1),
                (0, -1),           (0, 1),
                (1, -1),  (1, 0),  (1, 1),
            ]
        else:  # von_neumann
            self._n_offsets = [(-1, 0), (1, 0), (0, -1), (0, 1)]
        self._n_neighbors = len(self._n_offsets)

    def setup(self) -> dict[str, Any]:
        """Initialize opinions."""
        self.rng = np.random.default_rng(self.seed)
        self._step_count = 0
        R, C = self.rows, self.cols

        if self.init_mode == "random":
            self.grid = (self.rng.random((R, C)) < 0.5).astype(np.int8)
        elif self.init_mode == "biased":
            self.grid = (self.rng.random((R, C)) < self.init_fraction).astype(np.int8)
        elif self.init_mode == "half_and_half":
            self.grid = np.zeros((R, C), dtype=np.int8)
            self.grid[:, C // 2:] = 1
        else:
            raise ValueError(f"Unknown init_mode: {self.init_mode}")

        return self._state_dict()

    def step(self) -> dict[str, Any]:
        """One Monte Carlo sweep = R*C site-updates.

        Per site-update: pick uniform random site, pick uniform random neighbor,
        copy neighbor's opinion onto site. This is the canonical asynchronous
        voter model dynamics (Clifford & Sudbury 1973; Liggett 1985).

        A parallelized "checkerboard" step was trialed but removed: while it
        yields statistically compatible late-time coarsening exponents, its
        early-time trajectories differ by several standard deviations from
        the canonical dynamics. The Sprint 20 characterization therefore
        uses this canonical step exclusively.
        """
        R, C = self.rows, self.cols
        N = R * C

        # Pre-generate randomness for the full sweep
        site_rs = self.rng.integers(0, R, size=N)
        site_cs = self.rng.integers(0, C, size=N)
        nb_idx = self.rng.integers(0, self._n_neighbors, size=N)

        for k in range(N):
            r = site_rs[k]
            c = site_cs[k]
            dr, dc = self._n_offsets[nb_idx[k]]
            nr = (r + dr) % R
            nc = (c + dc) % C
            self.grid[r, c] = self.grid[nr, nc]

        self._step_count += 1
        return self._state_dict()

    def run(self, n_steps: int = 500) -> list[dict[str, Any]]:
        """Run for n_steps Monte Carlo sweeps (each = R*C site-updates)."""
        history = [self.setup()]
        for _ in range(n_steps):
            history.append(self.step())
        return history

    def run_until_consensus(
        self, max_steps: int = 100_000, record_every: int = 1
    ) -> list[dict[str, Any]]:
        """Run until consensus reached or max_steps exceeded.

        Records every `record_every` sweeps to limit history length for large L.
        """
        history = [self.setup()]
        for t in range(1, max_steps + 1):
            state = self.step()
            if t % record_every == 0 or state["consensus_reached"]:
                history.append(state)
            if state["consensus_reached"]:
                break
        return history

    def _state_dict(self) -> dict[str, Any]:
        """Build state dictionary."""
        frac_1 = float(np.mean(self.grid == 1))
        magnetization = 2.0 * frac_1 - 1.0
        abs_m = abs(magnetization)
        wall_density = _wall_density(self.grid, self._n_offsets)
        moran = _moran_i_fast(self.grid.astype(float), self._n_offsets)
        consensus = bool(frac_1 == 0.0 or frac_1 == 1.0)

        return {
            "grid": self.grid.copy(),
            "grid_dims": (self.rows, self.cols),
            "step": self._step_count,
            "magnetization": magnetization,
            "abs_magnetization": abs_m,
            "wall_density": wall_density,
            "moran_i": moran,
            "consensus_reached": consensus,
            "n_states": 2,
        }

    def get_metadata(self) -> dict[str, Any]:
        """Return model metadata."""
        return {
            "model": "voter",
            "rows": self.rows,
            "cols": self.cols,
            "neighborhood": self.neighborhood,
            "boundary": self.boundary,
            "substrate": "lattice_2d",
            "update": "asynchronous_copy_neighbor",
            "has_movement": False,
        }


def _wall_density(grid: np.ndarray, offsets: list[tuple[int, int]]) -> float:
    """Fraction of neighbor pairs with differing states (periodic torus).

    Each unordered pair is counted once: we only sum over a canonical half of
    offsets. For Moore: (−1,−1), (−1,0), (−1,1), (0,1) covers all 4 unique
    directions. For Von Neumann: (−1,0), (0,1).
    """
    # Use half the offsets to avoid double-counting
    # Half-space: dr < 0, or (dr == 0 and dc > 0)
    half_offsets = [(dr, dc) for (dr, dc) in offsets
                    if dr < 0 or (dr == 0 and dc > 0)]
    total = 0
    diff = 0
    N = grid.size
    for dr, dc in half_offsets:
        shifted = np.roll(np.roll(grid, -dr, axis=0), -dc, axis=1)
        diff += int(np.sum(grid != shifted))
        total += N
    return diff / total if total > 0 else 0.0


def _moran_i_fast(grid: np.ndarray, offsets: list[tuple[int, int]]) -> float:
    """Moran's I with the same neighborhood (periodic torus)."""
    N = grid.size
    x = grid - grid.mean()
    var = float(np.mean(x ** 2))
    if var < 1e-12:
        # Consensus: no spatial variance; convention return 0
        return 0.0

    cross = 0.0
    W = 0
    for dr, dc in offsets:
        shifted = np.roll(np.roll(x, -dr, axis=0), -dc, axis=1)
        cross += float(np.sum(x * shifted))
        W += N

    return (N * cross) / (W * N * var)
