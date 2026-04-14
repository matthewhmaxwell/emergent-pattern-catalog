"""Greenberg-Hastings excitable cellular automaton.

Reference: Greenberg, J.M. & Hastings, S.P. (1978). "Spatial patterns for
discrete models of diffusion in excitable media." SIAM J. Applied Mathematics,
34(3), 515-523.

Also: Fisch, R., Gravner, J. & Griffeath, D. (1991). "Threshold-range scaling
of excitable cellular automata." Statistics and Computing, 1, 23-39.

The GH automaton is parameterized by (κ, θ, ρ):
  κ (n_states): number of states arranged on a color wheel (0, 1, ..., κ-1).
    State 0 = resting, state 1 = excited, states 2..κ-1 = refractory.
  θ (threshold): minimum number of excited neighbors to trigger excitation.
  ρ (neighborhood): 'von_neumann' (4-connected) or 'moore' (8-connected).

Update rules (synchronous):
  - Resting (0): if ≥ θ neighbors are excited (state 1) → become excited (1),
    otherwise remain resting (0).
  - Excited (1): automatically advance to first refractory state (2).
  - Refractory (s, 2 ≤ s < κ-1): advance to s+1.
  - Final refractory (κ-1): return to resting (0).

Boundary conditions: periodic (torus) or fixed (dead border).

Initial conditions:
  - 'random': fraction p of cells randomly assigned to each state
  - 'single_seed': one excited cell at center
  - 'broken_wave': horizontal wavefront with a gap to seed spiral
  - 'custom': user-supplied grid

State history keys:
    grid       : np.ndarray (rows, cols), dtype int — cell states [0, κ)
    grid_dims  : tuple (rows, cols)
    n_states   : int — κ
    step       : int — timestep
    excited_count   : int — number of cells in state 1
    resting_count   : int — number of cells in state 0
    refractory_count: int — number of cells in states 2..κ-1
    wavefront_count : int — resting→excited transitions this step
    activity_density: float — fraction of excited cells
"""

from __future__ import annotations

from typing import Any

import numpy as np

from epc.base_model import BaseModel


class GreenbergHastings(BaseModel):
    """Greenberg-Hastings excitable cellular automaton.

    Parameters
    ----------
    rows : int
        Grid height (default 100).
    cols : int
        Grid width (default 100).
    n_states : int
        Number of states κ ≥ 3 (default 3: rest/excited/refractory).
    threshold : int
        Excitation threshold θ ≥ 1 (default 1).
    neighborhood : str
        'von_neumann' (4-connected) or 'moore' (8-connected).
    boundary : str
        'periodic' (torus) or 'fixed' (dead border).
    init_mode : str
        'random', 'single_seed', 'broken_wave', or 'custom'.
    init_density : float
        For 'random' mode: fraction of cells initially non-resting.
        Distributed uniformly among states 1..κ-1.
    init_grid : np.ndarray, optional
        For 'custom' mode: user-supplied initial grid.
    seed : int
        Random seed.
    """

    def __init__(
        self,
        rows: int = 100,
        cols: int = 100,
        n_states: int = 3,
        threshold: int = 1,
        neighborhood: str = "moore",
        boundary: str = "periodic",
        init_mode: str = "random",
        init_density: float = 0.3,
        init_grid: np.ndarray | None = None,
        seed: int = 42,
    ) -> None:
        super().__init__(seed=seed)
        if n_states < 3:
            raise ValueError(f"n_states must be ≥ 3, got {n_states}")
        if threshold < 1:
            raise ValueError(f"threshold must be ≥ 1, got {threshold}")
        if neighborhood not in ("von_neumann", "moore"):
            raise ValueError(f"neighborhood must be 'von_neumann' or 'moore'")
        if boundary not in ("periodic", "fixed"):
            raise ValueError(f"boundary must be 'periodic' or 'fixed'")

        self.rows = rows
        self.cols = cols
        self.n_states = n_states
        self.threshold = threshold
        self.neighborhood = neighborhood
        self.boundary = boundary
        self.init_mode = init_mode
        self.init_density = init_density
        self.init_grid = init_grid

        self._grid: np.ndarray | None = None
        self._wavefront_count: int = 0

        # Precompute neighbor offsets
        if neighborhood == "von_neumann":
            self._offsets = [(-1, 0), (1, 0), (0, -1), (0, 1)]
        else:  # moore
            self._offsets = [
                (-1, -1), (-1, 0), (-1, 1),
                (0, -1),           (0, 1),
                (1, -1),  (1, 0),  (1, 1),
            ]

    def setup(self) -> dict[str, Any]:
        if self.init_mode == "random":
            self._grid = self._init_random()
        elif self.init_mode == "single_seed":
            self._grid = self._init_single_seed()
        elif self.init_mode == "broken_wave":
            self._grid = self._init_broken_wave()
        elif self.init_mode == "custom":
            if self.init_grid is None:
                raise ValueError("init_grid required for 'custom' mode")
            self._grid = self.init_grid.copy().astype(int)
        else:
            raise ValueError(f"Unknown init_mode: {self.init_mode}")

        self._step_count = 0
        self._wavefront_count = 0
        self._is_setup = True

        return self._snapshot()

    def step(self) -> dict[str, Any]:
        """Synchronous update of all cells."""
        old = self._grid
        new = np.zeros_like(old)
        wavefront = 0

        for r in range(self.rows):
            for c in range(self.cols):
                state = old[r, c]
                if state == 0:
                    # Resting: check if enough neighbors are excited
                    excited_neighbors = self._count_excited_neighbors(old, r, c)
                    if excited_neighbors >= self.threshold:
                        new[r, c] = 1
                        wavefront += 1
                    else:
                        new[r, c] = 0
                elif state == 1:
                    # Excited → first refractory
                    new[r, c] = 2
                else:
                    # Refractory: advance, or return to resting
                    new[r, c] = (state + 1) % self.n_states
                    # Note: if n_states=3, state 2 → 0 (resting)

        self._grid = new
        self._wavefront_count = wavefront
        self._step_count += 1

        return self._snapshot()

    def step_vectorized(self) -> dict[str, Any]:
        """Vectorized synchronous update using numpy — much faster for large grids."""
        old = self._grid
        new = np.zeros_like(old)

        # Count excited neighbors for each cell
        excited_map = (old == 1).astype(int)
        neighbor_count = np.zeros((self.rows, self.cols), dtype=int)

        for dr, dc in self._offsets:
            if self.boundary == "periodic":
                shifted = np.roll(np.roll(excited_map, -dr, axis=0), -dc, axis=1)
            else:
                shifted = np.zeros_like(excited_map)
                # Source and destination slices
                sr = slice(max(0, dr), min(self.rows, self.rows + dr))
                sc = slice(max(0, dc), min(self.cols, self.cols + dc))
                dr_slice = slice(max(0, -dr), min(self.rows, self.rows - dr))
                dc_slice = slice(max(0, -dc), min(self.cols, self.cols - dc))
                shifted[dr_slice, dc_slice] = excited_map[sr, sc]
            neighbor_count += shifted

        # Apply rules
        resting = old == 0
        excited = old == 1
        refractory = (old >= 2)

        # Resting cells with enough excited neighbors → excited
        newly_excited = resting & (neighbor_count >= self.threshold)
        new[newly_excited] = 1

        # Resting cells without enough neighbors → stay resting
        new[resting & ~newly_excited] = 0

        # Excited → first refractory (state 2)
        new[excited] = 2

        # Refractory → advance (or return to resting if at last state)
        new[refractory] = (old[refractory] + 1) % self.n_states

        self._wavefront_count = int(newly_excited.sum())
        self._grid = new
        self._step_count += 1

        return self._snapshot()

    def run(self, n_steps: int, record_every: int = 1, vectorized: bool = True) -> list[dict[str, Any]]:
        """Execute simulation with optional vectorized stepping.

        Parameters
        ----------
        n_steps : int
            Number of timesteps.
        record_every : int
            Record state every N steps.
        vectorized : bool
            Use numpy-vectorized update (much faster). Default True.
        """
        if not self._is_setup:
            initial_state = self.setup()
            self._state_history = [initial_state]

        step_fn = self.step_vectorized if vectorized else self.step

        for i in range(n_steps):
            state = step_fn()
            if (i + 1) % record_every == 0:
                self._state_history.append(state)

        return self._state_history

    def get_metadata(self) -> dict[str, Any]:
        return {
            "model_name": "greenberg_hastings",
            "model_class": "excitable_ca",
            "rows": self.rows,
            "cols": self.cols,
            "n_states": self.n_states,
            "threshold": self.threshold,
            "neighborhood": self.neighborhood,
            "boundary": self.boundary,
            "init_mode": self.init_mode,
            "init_density": self.init_density,
            "interaction_type": "local_state_transition",
            "update_mode": "synchronous",
            "seed": self.seed,
            "reference": "Greenberg & Hastings (1978), Fisch, Gravner & Griffeath (1991)",
        }

    def get_timescale(self) -> float:
        """System-intrinsic timescale: wave propagation time across grid.

        T_prop = max(rows, cols) for periodic boundary (wave crosses grid).
        This is the natural timescale for P13 detection.
        """
        return float(max(self.rows, self.cols))

    def is_converged(self) -> bool:
        """GH converges when all cells are resting (activity dies out)."""
        if self._grid is None:
            return False
        return int((self._grid > 0).sum()) == 0

    # --- Initialization methods ---

    def _init_random(self) -> np.ndarray:
        """Random initial state with specified density of non-resting cells."""
        grid = np.zeros((self.rows, self.cols), dtype=int)
        n_active = int(self.rows * self.cols * self.init_density)
        # Randomly choose positions
        positions = self.rng.choice(
            self.rows * self.cols, size=n_active, replace=False
        )
        # Distribute among states 1..n_states-1
        states = self.rng.integers(1, self.n_states, size=n_active)
        for pos, state in zip(positions, states):
            r, c = divmod(pos, self.cols)
            grid[r, c] = state
        return grid

    def _init_single_seed(self) -> np.ndarray:
        """Single excited cell at center, everything else resting."""
        grid = np.zeros((self.rows, self.cols), dtype=int)
        grid[self.rows // 2, self.cols // 2] = 1
        return grid

    def _init_broken_wave(self) -> np.ndarray:
        """Horizontal wavefront with a gap — seeds spiral formation.

        Creates a vertical line of excited cells at the center column,
        with refractory cells to the left (simulating a wavefront that
        already passed), and a gap in the middle to create a broken
        wave tip that will curl into a spiral.
        """
        grid = np.zeros((self.rows, self.cols), dtype=int)
        mid_c = self.cols // 2
        mid_r = self.rows // 2
        gap_half = max(2, self.rows // 10)

        # Excited front (vertical line)
        for r in range(self.rows):
            # Gap in the middle
            if mid_r - gap_half <= r <= mid_r + gap_half:
                continue
            grid[r, mid_c] = 1

        # Refractory wake behind the front (to the left)
        for s in range(2, self.n_states):
            col = mid_c - (s - 1)
            if col < 0:
                if self.boundary == "periodic":
                    col = col % self.cols
                else:
                    continue
            for r in range(self.rows):
                if mid_r - gap_half <= r <= mid_r + gap_half:
                    continue
                grid[r, col] = s

        return grid

    # --- Neighbor counting ---

    def _count_excited_neighbors(self, grid: np.ndarray, r: int, c: int) -> int:
        """Count neighbors in state 1 (excited)."""
        count = 0
        for dr, dc in self._offsets:
            nr, nc = r + dr, c + dc
            if self.boundary == "periodic":
                nr = nr % self.rows
                nc = nc % self.cols
            elif not (0 <= nr < self.rows and 0 <= nc < self.cols):
                continue
            if grid[nr, nc] == 1:
                count += 1
        return count

    # --- State snapshot ---

    def _snapshot(self) -> dict[str, Any]:
        grid = self._grid.copy()
        excited = int((grid == 1).sum())
        resting = int((grid == 0).sum())
        refractory = int(((grid >= 2)).sum())
        total = self.rows * self.cols

        return {
            "grid": grid,
            "grid_dims": (self.rows, self.cols),
            "n_states": self.n_states,
            "step": self._step_count,
            "excited_count": excited,
            "resting_count": resting,
            "refractory_count": refractory,
            "wavefront_count": self._wavefront_count,
            "activity_density": excited / total,
        }
