"""SIR epidemic cellular automaton on a 2D lattice.

Reference: Fuks, H. & Lawniczak, A.T. (2002). "Individual-based lattice model
for spatial spread of epidemics." Discrete Dynamics in Nature and Society,
6(3), 191-200.

Also: Datta, A. & Acharyya, M. (2021). "Modelling the Spread of an Epidemic in
Presence of Vaccination using Cellular Automata." arXiv:2104.10456.

Also: Rousseau, G., Giorgini, G., Livi, R. & Chaté, H. (1997). "Dynamical
phases in a cellular automaton model for epidemic propagation." Physica D,
103, 554-563.

The SIR CA is parameterized by (p, q, neighborhood):
  p (infection_prob): probability that a susceptible cell becomes infected
    given ≥1 infected neighbor.
  q (recovery_prob): probability that an infected cell recovers each timestep.
  neighborhood: 'moore' (8-connected) or 'von_neumann' (4-connected).

Update rules (synchronous):
  - Susceptible (0): if ≥1 neighbor is infected, become infected with prob p.
    For multiple infected neighbors, infection is attempted independently per
    neighbor: P(infection) = 1 - (1-p)^n_infected_neighbors (standard model).
  - Infected (1): recover (→ state 2) with probability q.
  - Recovered (2): permanent immunity, no state change.

Boundary conditions: periodic (torus) or fixed (dead border).

Initial conditions:
  - 'single_seed': one infected cell at center — produces circular wavefront
  - 'random_fraction': fraction f of cells randomly set to infected
  - 'custom': user-supplied grid

State history keys:
    grid              : np.ndarray (rows, cols), dtype int — cell states {0,1,2}
    grid_dims         : tuple (rows, cols)
    n_states          : int — always 3
    step              : int — timestep
    s_count           : int — susceptible cells
    i_count           : int — infected cells
    r_count           : int — recovered cells
    newly_infected    : int — S→I transitions this step
    newly_recovered   : int — I→R transitions this step
    activity_density  : float — fraction of infected cells

Key published results for replication:
  1. Epidemic curve: I(t) rises then falls (bell shape), matching Kermack-McKendrick.
  2. Single-seed wavefront: circular expansion with linear radius growth.
  3. Phase transition: critical infection probability p_c depends on q and
     neighborhood; below p_c, epidemic dies out; above, percolates.
  4. Final epidemic size: R_∞/N is a function of the effective R0 = p*n_neighbors/q.

P13 boundary test:
  SIR has 3 states (S/I/R) → passes the P13 n_states≥3 guard.
  However, SIR wavefronts are SINGLE-PASS: susceptible cells can only be
  infected once, then recover permanently. Unlike GH excitable media where
  refractory cells return to resting and can be re-excited, SIR recovered
  cells never become susceptible again. This means:
  - No spiral formation (wavefront has no medium to re-enter)
  - Activity dies out once the epidemic burns through
  - P13 should reject on the persistence check (activity dies out)
  The discrimination is through epidemic transience, not wavefront speed.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from epc.base_model import BaseModel


class SIREpidemicModel(BaseModel):
    """SIR epidemic cellular automaton.

    Parameters
    ----------
    rows : int
        Grid height (default 100).
    cols : int
        Grid width (default 100).
    infection_prob : float
        Per-neighbor infection probability p ∈ (0, 1). Default 0.3.
    recovery_prob : float
        Per-timestep recovery probability q ∈ (0, 1). Default 0.1.
    neighborhood : str
        'moore' (8-connected) or 'von_neumann' (4-connected).
    boundary : str
        'periodic' (torus) or 'fixed' (dead border).
    init_mode : str
        'single_seed', 'random_fraction', or 'custom'.
    init_fraction : float
        For 'random_fraction': fraction of cells initially infected.
    init_grid : np.ndarray, optional
        For 'custom' mode: user-supplied initial grid.
    independent_neighbors : bool
        If True (default), each infected neighbor independently attempts
        infection: P(infect) = 1 - (1-p)^n_infected. This is the standard
        model in the literature. If False, use threshold model: P(infect) = p
        if ≥1 neighbor is infected.
    seed : int
        Random seed.
    """

    # State constants
    SUSCEPTIBLE = 0
    INFECTED = 1
    RECOVERED = 2

    def __init__(
        self,
        rows: int = 100,
        cols: int = 100,
        infection_prob: float = 0.3,
        recovery_prob: float = 0.1,
        neighborhood: str = "moore",
        boundary: str = "periodic",
        init_mode: str = "single_seed",
        init_fraction: float = 0.01,
        init_grid: np.ndarray | None = None,
        independent_neighbors: bool = True,
        seed: int = 42,
    ) -> None:
        super().__init__(seed=seed)
        if not 0 < infection_prob <= 1:
            raise ValueError(f"infection_prob must be in (0, 1], got {infection_prob}")
        if not 0 < recovery_prob <= 1:
            raise ValueError(f"recovery_prob must be in (0, 1], got {recovery_prob}")
        if neighborhood not in ("von_neumann", "moore"):
            raise ValueError(f"neighborhood must be 'von_neumann' or 'moore'")
        if boundary not in ("periodic", "fixed"):
            raise ValueError(f"boundary must be 'periodic' or 'fixed'")

        self.rows = rows
        self.cols = cols
        self.infection_prob = infection_prob
        self.recovery_prob = recovery_prob
        self.neighborhood = neighborhood
        self.boundary = boundary
        self.init_mode = init_mode
        self.init_fraction = init_fraction
        self.init_grid = init_grid
        self.independent_neighbors = independent_neighbors

        self._grid: np.ndarray | None = None
        self._newly_infected: int = 0
        self._newly_recovered: int = 0

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
        if self.init_mode == "single_seed":
            self._grid = self._init_single_seed()
        elif self.init_mode == "random_fraction":
            self._grid = self._init_random_fraction()
        elif self.init_mode == "custom":
            if self.init_grid is None:
                raise ValueError("init_grid required for 'custom' mode")
            self._grid = self.init_grid.copy().astype(int)
        else:
            raise ValueError(f"Unknown init_mode: {self.init_mode}")

        self._step_count = 0
        self._newly_infected = 0
        self._newly_recovered = 0
        self._is_setup = True

        return self._snapshot()

    def step(self) -> dict[str, Any]:
        """Vectorized synchronous update."""
        old = self._grid
        new = old.copy()

        # --- Infection: susceptible cells check infected neighbors ---
        infected_map = (old == self.INFECTED).astype(np.float64)

        # Count infected neighbors for each cell
        neighbor_count = np.zeros((self.rows, self.cols), dtype=np.float64)
        for dr, dc in self._offsets:
            if self.boundary == "periodic":
                shifted = np.roll(np.roll(infected_map, -dr, axis=0), -dc, axis=1)
            else:
                shifted = np.zeros_like(infected_map)
                sr = slice(max(0, dr), min(self.rows, self.rows + dr))
                sc = slice(max(0, dc), min(self.cols, self.cols + dc))
                dr_slice = slice(max(0, -dr), min(self.rows, self.rows - dr))
                dc_slice = slice(max(0, -dc), min(self.cols, self.cols - dc))
                shifted[dr_slice, dc_slice] = infected_map[sr, sc]
            neighbor_count += shifted

        susceptible_mask = old == self.SUSCEPTIBLE
        has_infected_neighbor = neighbor_count > 0

        # Compute infection probability
        if self.independent_neighbors:
            # P(infection) = 1 - (1 - p)^n_infected_neighbors
            # Only compute where susceptible AND has neighbors
            candidates = susceptible_mask & has_infected_neighbor
            infection_probs = np.zeros((self.rows, self.cols))
            infection_probs[candidates] = 1.0 - (
                (1.0 - self.infection_prob) ** neighbor_count[candidates]
            )
        else:
            # Threshold model: P = p if any neighbor infected
            infection_probs = np.where(
                susceptible_mask & has_infected_neighbor,
                self.infection_prob,
                0.0,
            )

        # Draw random numbers and infect
        rand_infect = self.rng.random((self.rows, self.cols))
        newly_infected_mask = susceptible_mask & (rand_infect < infection_probs)
        new[newly_infected_mask] = self.INFECTED

        # --- Recovery: infected cells recover with probability q ---
        infected_mask = old == self.INFECTED
        rand_recover = self.rng.random((self.rows, self.cols))
        newly_recovered_mask = infected_mask & (rand_recover < self.recovery_prob)
        new[newly_recovered_mask] = self.RECOVERED

        self._newly_infected = int(newly_infected_mask.sum())
        self._newly_recovered = int(newly_recovered_mask.sum())
        self._grid = new
        self._step_count += 1

        return self._snapshot()

    def run(
        self, n_steps: int, record_every: int = 1
    ) -> list[dict[str, Any]]:
        """Execute simulation.

        Parameters
        ----------
        n_steps : int
            Number of timesteps.
        record_every : int
            Record state every N steps.
        """
        if not self._is_setup:
            initial_state = self.setup()
            self._state_history = [initial_state]

        for i in range(n_steps):
            state = self.step()
            if (i + 1) % record_every == 0:
                self._state_history.append(state)

            # Early termination: no infected cells left
            if state["i_count"] == 0:
                break

        return self._state_history

    def get_metadata(self) -> dict[str, Any]:
        n_neighbors = 4 if self.neighborhood == "von_neumann" else 8
        # Effective R0 approximation for lattice SIR
        # R0 ≈ (1 - (1-p)^n_neighbors) / q for independent model
        if self.independent_neighbors:
            effective_infection = 1.0 - (1.0 - self.infection_prob) ** n_neighbors
        else:
            effective_infection = self.infection_prob
        r0_approx = effective_infection / self.recovery_prob

        return {
            "model_name": "sir_epidemic",
            "model_class": "epidemic_ca",
            "rows": self.rows,
            "cols": self.cols,
            "n_states": 3,
            "infection_prob": self.infection_prob,
            "recovery_prob": self.recovery_prob,
            "neighborhood": self.neighborhood,
            "n_neighbors": n_neighbors,
            "boundary": self.boundary,
            "init_mode": self.init_mode,
            "init_fraction": self.init_fraction,
            "independent_neighbors": self.independent_neighbors,
            "r0_approx": r0_approx,
            "interaction_type": "local_state_transition",
            "update_mode": "synchronous",
            "seed": self.seed,
            "reference": (
                "Fuks & Lawniczak (2002), Datta & Acharyya (2021), "
                "Rousseau et al. (1997)"
            ),
        }

    def get_timescale(self) -> float:
        """System-intrinsic timescale: wave propagation time across grid.

        For epidemic wavefronts, T_prop ≈ max(rows, cols) / v_front.
        Wavefront speed v ≈ 1 cell/step for typical parameters, so
        T_prop ≈ max(rows, cols). Same convention as GH.
        """
        return float(max(self.rows, self.cols))

    def is_converged(self) -> bool:
        """SIR converges when no infected cells remain (epidemic over)."""
        if self._grid is None:
            return False
        return int((self._grid == self.INFECTED).sum()) == 0

    # --- Initialization methods ---

    def _init_single_seed(self) -> np.ndarray:
        """Single infected cell at center, all others susceptible."""
        grid = np.zeros((self.rows, self.cols), dtype=int)
        grid[self.rows // 2, self.cols // 2] = self.INFECTED
        return grid

    def _init_random_fraction(self) -> np.ndarray:
        """Random fraction of cells initially infected."""
        grid = np.zeros((self.rows, self.cols), dtype=int)
        n_infected = max(1, int(self.rows * self.cols * self.init_fraction))
        positions = self.rng.choice(
            self.rows * self.cols, size=n_infected, replace=False
        )
        for pos in positions:
            r, c = divmod(pos, self.cols)
            grid[r, c] = self.INFECTED
        return grid

    # --- State snapshot ---

    def _snapshot(self) -> dict[str, Any]:
        grid = self._grid.copy()
        s_count = int((grid == self.SUSCEPTIBLE).sum())
        i_count = int((grid == self.INFECTED).sum())
        r_count = int((grid == self.RECOVERED).sum())
        total = self.rows * self.cols

        return {
            "grid": grid,
            "grid_dims": (self.rows, self.cols),
            "n_states": 3,
            "step": self._step_count,
            "s_count": s_count,
            "i_count": i_count,
            "r_count": r_count,
            "newly_infected": self._newly_infected,
            "newly_recovered": self._newly_recovered,
            "activity_density": i_count / total,
            # Epidemic-specific observables
            "s_fraction": s_count / total,
            "i_fraction": i_count / total,
            "r_fraction": r_count / total,
        }
