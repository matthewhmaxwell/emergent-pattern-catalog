"""Nowak-May Spatial Prisoner's Dilemma — canonical model for P27.

Reference:
    Nowak, M.A. & May, R.M. (1992). "Evolutionary games and spatial chaos."
    Nature, 359(6398), 826-829. DOI: 10.1038/359826a0

Simplified payoff matrix (Nowak-May parameterization):
    Both cooperate: R = 1
    Both defect:    P = 0
    C vs D:         S = 0 (sucker)
    D vs C:         T = b (temptation)

    b > 1: defection temptation exceeds mutual cooperation reward.
    In well-mixed populations, defection dominates for any b > 1.
    On a lattice, cooperator clusters survive for b < 2 via spatial reciprocity.

Update rule (synchronous, deterministic):
    Each cell plays PD with all Moore neighbors + self (9 games).
    Total payoff = sum of individual game payoffs.
    Copy strategy of neighbor (including self) with highest total payoff.
    Ties broken: keep current strategy (Nowak-May convention).

Key results (Nowak & May 1992):
    b = 1.0:  all cooperate (no temptation)
    b = 1.5:  ~60% cooperators, dynamic clusters
    b = 1.8:  ~30% cooperators, kaleidoscope patterns
    b = 2.0:  cooperators extinct (too much temptation)
    b = 1.8 from 50% random init: fractal-like boundary dynamics

State history keys:
    grid       : np.ndarray (rows, cols), dtype int — strategy {0=C, 1=D}
    grid_dims  : tuple (rows, cols)
    step       : int — generation number
    coop_fraction : float — fraction of cooperators
    moran_i    : float — Moran's I of cooperator indicator
"""

from __future__ import annotations

from typing import Any, Optional

import numpy as np


class NowakMayModel:
    """Nowak-May Spatial Prisoner's Dilemma (1992).

    Parameters
    ----------
    rows : int
        Grid height.
    cols : int
        Grid width.
    b : float
        Temptation payoff (T = b). Must be > 1 for PD structure.
    init_mode : str
        'random' (fraction p cooperators), 'single_defector' (all C except center),
        'checkerboard' (alternating C/D).
    init_coop_fraction : float
        For 'random' mode: initial fraction of cooperators.
    boundary : str
        'periodic' (torus) or 'fixed' (dead border, all-C neighbors outside).
    seed : int
        Random seed.
    """

    def __init__(
        self,
        rows: int = 100,
        cols: int = 100,
        b: float = 1.8,
        init_mode: str = "random",
        init_coop_fraction: float = 0.5,
        boundary: str = "periodic",
        seed: int = 42,
    ):
        self.rows = rows
        self.cols = cols
        self.b = b
        self.init_mode = init_mode
        self.init_coop_fraction = init_coop_fraction
        self.boundary = boundary
        self.seed = seed

        # State
        self.grid: np.ndarray = np.empty(0)
        self.rng: np.random.Generator = np.random.default_rng(seed)
        self._step_count: int = 0

    def setup(self) -> dict[str, Any]:
        """Initialize strategies."""
        self.rng = np.random.default_rng(self.seed)
        self._step_count = 0
        R, C = self.rows, self.cols

        if self.init_mode == "random":
            self.grid = (self.rng.random((R, C)) >= self.init_coop_fraction).astype(np.int8)
        elif self.init_mode == "single_defector":
            self.grid = np.zeros((R, C), dtype=np.int8)  # all cooperate
            self.grid[R // 2, C // 2] = 1  # one defector
        elif self.init_mode == "checkerboard":
            self.grid = np.indices((R, C)).sum(axis=0).astype(np.int8) % 2
        else:
            raise ValueError(f"Unknown init_mode: {self.init_mode}")

        return self._state_dict()

    def step(self) -> dict[str, Any]:
        """One synchronous generation update."""
        R, C = self.rows, self.cols
        b = self.b

        # Compute total payoff for each cell
        # Strategy: 0 = C, 1 = D
        # Payoff from game between i and j:
        #   C vs C: 1, C vs D: 0, D vs C: b, D vs D: 0
        # So payoff_i from playing j = (1 - s_i) * (1 - s_j) * 1 + s_i * (1 - s_j) * b
        #                              = (1-s_j) * [(1-s_i) + s_i * b]
        #                              = (1-s_j) * [1 + s_i * (b-1)]
        # Equivalently: each C neighbor gives 1 to C, b to D. Each D neighbor gives 0 to both.

        # Count cooperator neighbors (Moore + self)
        if self.boundary == "periodic":
            padded = np.pad(self.grid, 1, mode="wrap")
        else:
            padded = np.pad(self.grid, 1, mode="constant", constant_values=0)

        # Number of cooperators in Moore neighborhood + self
        coop_neighbors = np.zeros((R, C), dtype=np.float64)
        for dr in (-1, 0, 1):
            for dc in (-1, 0, 1):
                coop_neighbors += (1 - padded[1 + dr:R + 1 + dr, 1 + dc:C + 1 + dc])

        # Payoff: cooperators get 1 per C-neighbor, defectors get b per C-neighbor
        payoff = np.where(self.grid == 0, coop_neighbors * 1.0, coop_neighbors * b)

        # Find best strategy in neighborhood: copy highest-payoff neighbor
        # (ties: keep current strategy)
        new_grid = self.grid.copy()

        if self.boundary == "periodic":
            padded_payoff = np.pad(payoff, 1, mode="wrap")
            padded_strat = np.pad(self.grid, 1, mode="wrap")
        else:
            padded_payoff = np.pad(payoff, 1, mode="constant", constant_values=-1)
            padded_strat = np.pad(self.grid, 1, mode="constant", constant_values=0)

        # For each cell, find neighbor with max payoff
        best_payoff = payoff.copy()  # self payoff
        best_strat = self.grid.copy()

        for dr in (-1, 0, 1):
            for dc in (-1, 0, 1):
                if dr == 0 and dc == 0:
                    continue
                nb_pay = padded_payoff[1 + dr:R + 1 + dr, 1 + dc:C + 1 + dc]
                nb_str = padded_strat[1 + dr:R + 1 + dr, 1 + dc:C + 1 + dc]

                # Strict > for update (ties keep current)
                better = nb_pay > best_payoff
                best_payoff = np.where(better, nb_pay, best_payoff)
                best_strat = np.where(better, nb_str, best_strat)

        self.grid = best_strat.astype(np.int8)
        self._step_count += 1
        return self._state_dict()

    def run(self, n_steps: int = 500) -> list[dict[str, Any]]:
        """Run for n_steps generations."""
        history = [self.setup()]
        for _ in range(n_steps):
            history.append(self.step())
        return history

    def _state_dict(self) -> dict[str, Any]:
        """Build state dictionary."""
        coop_frac = float(np.mean(self.grid == 0))

        # Moran's I for cooperator indicator
        indicator = (self.grid == 0).astype(float)
        moran_i = _moran_i_fast(indicator, self.boundary == "periodic")

        return {
            "grid": self.grid.copy(),
            "grid_dims": (self.rows, self.cols),
            "step": self._step_count,
            "coop_fraction": coop_frac,
            "moran_i": moran_i,
            "n_states": 2,
        }

    def get_metadata(self) -> dict[str, Any]:
        """Return model metadata for detector card checks."""
        return {
            "model": "nowak_may",
            "b": self.b,
            "payoff_T": self.b,
            "payoff_R": 1.0,
            "payoff_P": 0.0,
            "payoff_S": 0.0,
            "pd_structure": self.b > 1.0,  # T > R > P = S
            "substrate": "lattice_2d",
            "update": "synchronous_imitation",
            "has_movement": False,
        }


def _moran_i_fast(grid: np.ndarray, periodic: bool = True) -> float:
    """Compute Moran's I for a 2D binary grid (Moore neighborhood)."""
    N = grid.size
    x = grid - grid.mean()
    var = float(np.mean(x ** 2))
    if var < 1e-12:
        return 0.0

    # Sum of x_i * x_j for all adjacent pairs
    cross = 0.0
    W = 0
    if periodic:
        for dr, dc in [(-1, 0), (1, 0), (0, -1), (0, 1),
                        (-1, -1), (-1, 1), (1, -1), (1, 1)]:
            shifted = np.roll(np.roll(x, -dr, axis=0), -dc, axis=1)
            cross += np.sum(x * shifted)
            W += N
    else:
        rows, cols = grid.shape
        for dr, dc in [(-1, 0), (1, 0), (0, -1), (0, 1),
                        (-1, -1), (-1, 1), (1, -1), (1, 1)]:
            if dr == -1:
                s1, s2 = slice(1, None), slice(None, -1)
            elif dr == 1:
                s1, s2 = slice(None, -1), slice(1, None)
            else:
                s1 = s2 = slice(None)
            if dc == -1:
                t1, t2 = slice(1, None), slice(None, -1)
            elif dc == 1:
                t1, t2 = slice(None, -1), slice(1, None)
            else:
                t1 = t2 = slice(None)
            cross += np.sum(x[s1, t1] * x[s2, t2])
            W += x[s1, t1].size

    return float(cross / (W * var))
