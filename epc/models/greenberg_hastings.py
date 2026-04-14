"""Greenberg-Hastings excitable cellular automaton.

2D grid with periodic boundary conditions and synchronous update.
States: 0 = rest, 1 = excited, 2..n_refractory+1 = refractory stages.

Rules:
  - Rest cell → excited if ≥ threshold Moore neighbors are excited
  - Excited cell → highest refractory state (n_refractory + 1)
  - Refractory cell → decrement; state 2 → rest (0)

Canonical model for P13 (excitable spiral and target waves).

Reference: Greenberg, J.M. & Hastings, S.P. (1978). Spatial patterns
for discrete models of diffusion in excitable media. SIAM J. Appl. Math.
"""

from __future__ import annotations

from typing import Any

import numpy as np
from scipy.signal import convolve2d

from epc.base_model import BaseModel

# Moore neighborhood kernel (8 neighbors, not counting center)
_MOORE_KERNEL = np.array([[1, 1, 1],
                          [1, 0, 1],
                          [1, 1, 1]], dtype=np.int32)


class GreenbergHastingsCA(BaseModel):
    """Greenberg-Hastings excitable CA on a 2D periodic grid.

    Parameters
    ----------
    grid_size : int
        Side length of the square grid.
    n_refractory : int
        Number of distinct refractory stages. Total state count = n_refractory + 2.
    threshold : int
        Minimum excited Moore neighbors to excite a resting cell.
    density : float
        Initial fraction of cells in excited state (for 'random' init).
    init_mode : str
        'random' — sparse random excitation.
        'broken_wave' — planar wave with a cut, promotes spiral formation.
    """

    def __init__(
        self,
        grid_size: int = 100,
        n_refractory: int = 3,
        threshold: int = 1,
        density: float = 0.05,
        init_mode: str = "random",
        seed: int = 42,
    ) -> None:
        super().__init__(seed=seed)
        self.grid_size = grid_size
        self.n_refractory = n_refractory
        self.threshold = threshold
        self.density = density
        self.init_mode = init_mode
        self._grid: np.ndarray = np.zeros((grid_size, grid_size), dtype=np.int32)

    def setup(self) -> dict[str, Any]:
        g = self.grid_size
        self._grid = np.zeros((g, g), dtype=np.int32)

        if self.init_mode == "broken_wave":
            # Planar wave in left half, cut in upper-left quadrant
            # Creates a free end that curls into a spiral
            mid = g // 2
            quarter = g // 4
            # Excited wavefront
            self._grid[:, mid] = 1
            # Refractory trail behind wavefront
            for k in range(1, min(self.n_refractory + 1, mid)):
                self._grid[:, mid - k] = min(k + 1, self.n_refractory + 1)
            # Cut: remove upper half of the wave to create a free tip
            self._grid[:quarter, mid - self.n_refractory:mid + 1] = 0
        else:
            # Random sparse excitation
            n_excited = int(self.density * g * g)
            indices = self.rng.choice(g * g, size=n_excited, replace=False)
            rows, cols = np.unravel_index(indices, (g, g))
            self._grid[rows, cols] = 1

        self._step_count = 0
        self._is_setup = True
        return self._snapshot()

    def step(self) -> dict[str, Any]:
        grid = self._grid
        new_grid = np.zeros_like(grid)

        # Count excited neighbors (state == 1) using convolution with wrap
        excited_mask = (grid == 1).astype(np.int32)
        neighbor_count = convolve2d(
            excited_mask, _MOORE_KERNEL, mode="same", boundary="wrap"
        )

        # Rule 1: Rest (0) → Excited (1) if enough excited neighbors
        rest_mask = grid == 0
        new_grid[rest_mask & (neighbor_count >= self.threshold)] = 1

        # Rule 2: Excited (1) → highest refractory state
        new_grid[grid == 1] = self.n_refractory + 1

        # Rule 3: Refractory states decrement toward rest
        # State 2 → 0 (rest), state k>2 → k-1
        # Must handle state 2 specially: decrement to 0 (rest), not 1 (excited)
        refrac_high = grid > 2
        new_grid[refrac_high] = grid[refrac_high] - 1
        refrac_last = grid == 2
        new_grid[refrac_last] = 0

        self._grid = new_grid
        self._step_count += 1
        return self._snapshot()

    def _snapshot(self) -> dict[str, Any]:
        return {
            "grid": self._grid.copy(),
            "step": self._step_count,
            "n_excited": int((self._grid == 1).sum()),
            "n_refractory": int((self._grid > 1).sum()),
            "n_rest": int((self._grid == 0).sum()),
            "grid_dims": (self.grid_size, self.grid_size),
        }

    def get_metadata(self) -> dict[str, Any]:
        return {
            "model_name": "greenberg_hastings_ca",
            "model_class": "excitable_medium",
            "grid_size": self.grid_size,
            "n_refractory": self.n_refractory,
            "threshold": self.threshold,
            "density": self.density,
            "init_mode": self.init_mode,
            "interaction_type": "local_ca",
            "update_mode": "synchronous",
            "params": {
                "grid_size": self.grid_size,
                "n_refractory": self.n_refractory,
                "threshold": self.threshold,
                "density": self.density,
            },
            "seed": self.seed,
        }

    def get_timescale(self) -> float:
        """T_prop = grid_size / wavefront_speed.

        For threshold=1 on Moore neighborhood, speed ≈ 1 cell/step.
        """
        return float(self.grid_size)

    def is_converged(self) -> bool:
        """Quiescent if no excited or refractory cells remain."""
        return int((self._grid > 0).sum()) == 0
