"""Conway's Game of Life — canonical model for P15 (persistent computation).

Reference: Gardner, M. (1970). "The fantastic combinations of John Conway's
new solitaire game 'life'." Scientific American, 223(4), 120-123.

Rule B3/S23: A dead cell with exactly 3 live neighbors becomes alive.
A live cell with 2 or 3 live neighbors survives; otherwise it dies.

Deterministic, synchronous, binary-state, totalistic 2D CA on Moore
neighborhood. The simplest known system supporting universal computation
(Rendell 2002, via glider-based logic gates).

For P15 detection, GoL produces:
- Persistent propagating structures (gliders, spaceships)
- Information-bearing collisions (input-dependent outputs)
- Transfer Entropy > 0 at collision sites (Lizier et al. 2007, 2012)

State history keys:
    grid       : np.ndarray (rows, cols), dtype int — cell states {0, 1}
    grid_dims  : tuple (rows, cols)
    n_states   : int — always 2
    step       : int — timestep
    alive_count     : int — number of live cells
    activity_density: float — fraction of live cells
"""

from __future__ import annotations

from typing import Any

import numpy as np

from epc.base_model import BaseModel


class GameOfLife(BaseModel):
    """Conway's Game of Life (B3/S23).

    Parameters
    ----------
    rows : int
        Grid height (default 100).
    cols : int
        Grid width (default 100).
    boundary : str
        'periodic' (torus) or 'fixed' (dead border).
    init_mode : str
        'random', 'glider_collision', 'r_pentomino', or 'custom'.
    init_density : float
        For 'random' mode: initial fraction of live cells.
    init_grid : np.ndarray, optional
        For 'custom' mode.
    seed : int
        Random seed.
    """

    def __init__(
        self,
        rows: int = 100,
        cols: int = 100,
        boundary: str = "periodic",
        init_mode: str = "random",
        init_density: float = 0.3,
        init_grid: np.ndarray | None = None,
        seed: int = 42,
    ) -> None:
        super().__init__(seed=seed)
        self.rows = rows
        self.cols = cols
        self.boundary = boundary
        self.init_mode = init_mode
        self.init_density = init_density
        self.init_grid = init_grid

        self._grid: np.ndarray | None = None

    def setup(self) -> dict[str, Any]:
        if self.init_mode == "random":
            self._grid = (self.rng.random((self.rows, self.cols)) < self.init_density).astype(int)
        elif self.init_mode == "glider_collision":
            self._grid = self._init_glider_collision()
        elif self.init_mode == "r_pentomino":
            self._grid = self._init_r_pentomino()
        elif self.init_mode == "lwss":
            self._grid = self._init_lwss()
        elif self.init_mode == "custom":
            if self.init_grid is None:
                raise ValueError("init_grid required for 'custom' mode")
            self._grid = self.init_grid.copy().astype(int)
        else:
            raise ValueError(f"Unknown init_mode: {self.init_mode}")

        self._step_count = 0
        self._is_setup = True
        return self._snapshot()

    def step(self) -> dict[str, Any]:
        """Synchronous B3/S23 update using numpy convolution."""
        # Count live neighbors via convolution
        kernel = np.array([[1, 1, 1], [1, 0, 1], [1, 1, 1]])
        neighbor_count = self._convolve2d(self._grid, kernel)

        # Apply B3/S23 rules
        birth = (self._grid == 0) & (neighbor_count == 3)
        survive = (self._grid == 1) & ((neighbor_count == 2) | (neighbor_count == 3))
        self._grid = (birth | survive).astype(int)

        self._step_count += 1
        return self._snapshot()

    def _convolve2d(self, grid: np.ndarray, kernel: np.ndarray) -> np.ndarray:
        """2D convolution with boundary handling."""
        rows, cols = grid.shape
        kh, kw = kernel.shape
        ph, pw = kh // 2, kw // 2
        result = np.zeros_like(grid)

        if self.boundary == "periodic":
            padded = np.pad(grid, ((ph, ph), (pw, pw)), mode='wrap')
        else:
            padded = np.pad(grid, ((ph, ph), (pw, pw)), mode='constant', constant_values=0)

        for dr in range(kh):
            for dc in range(kw):
                result += kernel[dr, dc] * padded[dr:dr + rows, dc:dc + cols]

        return result

    def run(self, n_steps: int, record_every: int = 1) -> list[dict[str, Any]]:
        if not self._is_setup:
            initial_state = self.setup()
            self._state_history = [initial_state]

        for i in range(n_steps):
            state = self.step()
            if (i + 1) % record_every == 0:
                self._state_history.append(state)

        return self._state_history

    def get_metadata(self) -> dict[str, Any]:
        return {
            "model_name": "game_of_life",
            "model_class": "binary_ca",
            "rule": "B3/S23",
            "rows": self.rows,
            "cols": self.cols,
            "boundary": self.boundary,
            "init_mode": self.init_mode,
            "init_density": self.init_density,
            "interaction_type": "local_state_transition",
            "update_mode": "synchronous",
            "seed": self.seed,
            "reference": "Gardner (1970), Rendell (2002)",
        }

    def get_timescale(self) -> float:
        """T_prop: time for a glider to cross the grid.

        Glider speed = c/4 (1 cell per 4 steps diagonally).
        T_prop = 4 × max(rows, cols).
        """
        return 4.0 * max(self.rows, self.cols)

    def is_converged(self) -> bool:
        """GoL converges when all cells are dead."""
        if self._grid is None:
            return False
        return int(self._grid.sum()) == 0

    # --- Initialization ---

    def _init_glider_collision(self) -> np.ndarray:
        """Two gliders heading toward each other — will collide.

        Standard SE-moving glider and NW-moving glider placed to
        collide near the center.
        """
        grid = np.zeros((self.rows, self.cols), dtype=int)
        mid_r, mid_c = self.rows // 2, self.cols // 2

        # SE-moving glider (top-left region)
        # Pattern:
        #  .#.
        #  ..#
        #  ###
        r0, c0 = mid_r - 20, mid_c - 20
        glider_se = [(0, 1), (1, 2), (2, 0), (2, 1), (2, 2)]
        for dr, dc in glider_se:
            grid[(r0 + dr) % self.rows, (c0 + dc) % self.cols] = 1

        # NW-moving glider (bottom-right region) — rotated 180°
        # Pattern:
        #  ###
        #  #..
        #  .#.
        r1, c1 = mid_r + 20, mid_c + 20
        glider_nw = [(0, 0), (0, 1), (0, 2), (1, 0), (2, 1)]
        for dr, dc in glider_nw:
            grid[(r1 + dr) % self.rows, (c1 + dc) % self.cols] = 1

        return grid

    def _init_r_pentomino(self) -> np.ndarray:
        """R-pentomino: small pattern that produces complex, long-lived dynamics.

        The R-pentomino produces gliders, blocks, blinkers, and other
        structures over ~1100 generations before stabilizing. Good for
        testing P15 detection because it creates multiple structure types
        and collision events.

        Pattern:
         .##
         ##.
         .#.
        """
        grid = np.zeros((self.rows, self.cols), dtype=int)
        r0, c0 = self.rows // 2, self.cols // 2
        pattern = [(0, 1), (0, 2), (1, 0), (1, 1), (2, 1)]
        for dr, dc in pattern:
            grid[r0 + dr, c0 + dc] = 1
        return grid

    def _init_lwss(self) -> np.ndarray:
        """Lightweight spaceship (LWSS): travels at c/2 rightward.

        Period 4, displaces 2 cells right per period → speed c/2.
        9 cells. Rightward-moving pattern (from LifeWiki / secretGeek):
          .####.
          #...#.
          ....#.
          #..#..
        """
        grid = np.zeros((self.rows, self.cols), dtype=int)
        r0, c0 = self.rows // 2, self.cols // 4
        pattern = [
            (0, 1), (0, 2), (0, 3), (0, 4),
            (1, 0),                  (1, 4),
                                     (2, 4),
            (3, 0),         (3, 3),
        ]
        for dr, dc in pattern:
            grid[(r0 + dr) % self.rows, (c0 + dc) % self.cols] = 1
        return grid

    def _snapshot(self) -> dict[str, Any]:
        grid = self._grid.copy()
        alive = int(grid.sum())
        total = self.rows * self.cols

        return {
            "grid": grid,
            "grid_dims": (self.rows, self.cols),
            "n_states": 2,
            "step": self._step_count,
            "alive_count": alive,
            "activity_density": alive / total,
        }
