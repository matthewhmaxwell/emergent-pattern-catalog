"""Nagel-Schreckenberg cellular-automaton traffic model.

Primary reference:
  Nagel, K. & Schreckenberg, M. (1992). "A cellular automaton model for
  freeway traffic." J. Phys. I France 2, 2221-2229.
  https://doi.org/10.1051/jp1:1992277

Related references:
  Schadschneider, A. & Schreckenberg, M. (1997). "Cellular automata for
  traffic flow: analytical results." arXiv:cond-mat/9511037. — Mean-field
  and n-cluster analytic results; exact v_max = 1 case (particle-hole
  symmetric, peak flow at rho = 1/2).

  Bette, H. M., Habel, L., Emig, T. & Schreckenberg, M. (2017).
  "Mechanisms of jamming in the Nagel-Schreckenberg model for traffic
  flow." Phys. Rev. E 95, 012311.
  https://doi.org/10.1103/PhysRevE.95.012311
  — Analyses P(v=0) (stopped-car fraction) as the order parameter; this
  is EPC's P8 primary metric.

  Sasvari, M. & Kertesz, J. (1997). "Cellular automata models of single
  lane traffic." arXiv:cond-mat/9708114. — Jamming-transition analysis
  for v_max > 1.

MODEL RULES (faithful to Nagel-Schreckenberg 1992):

One-dimensional ring of L cells with periodic boundary conditions. Each
cell is either empty or holds exactly one car. Each car has an integer
velocity v in {0, 1, ..., v_max}. Cars move rightward.

Each time step is a PARALLEL update applying these four rules IN ORDER
to every car simultaneously:

  (1) Acceleration:   v_i <- min(v_i + 1, v_max)
  (2) Slowing down:   v_i <- min(v_i, d_i)
                      where d_i is the number of EMPTY cells between car
                      i and the car directly in front (so d_i = gap size
                      measured in empty cells). This is Schadschneider-
                      Schreckenberg's convention `v -> d - 1` with
                      d = gap_including_own_cell, which is equivalent.
  (3) Randomization:  v_i <- max(v_i - 1, 0)  with probability p_slow
                      (conditional on v_i > 0; cars at v=0 stay at 0).
  (4) Motion:         x_i <- (x_i + v_i) mod L

Notation note: in the EPC codebase, the "gap" d_i is defined as the
number of empty cells strictly between car i and the car in front. At
v_max = 5 a car with an empty cell ahead has d = 1 and can move at
v = 1; a car adjacent to the next has d = 0 and must set v = 0.

State history keys produced per snapshot:
    positions   : np.ndarray (n_cars,) int32  — cell indices, ascending on ring
    velocities  : np.ndarray (n_cars,) int8   — per-car velocities in 0..v_max
    gaps        : np.ndarray (n_cars,) int32  — gap to car in front (empty cells)
    n_cars      : int
    density     : float
    mean_velocity : float
    stopped_fraction : float
    step        : int

Canonical regimes (Sprint 15 characterization at L=1000, v_max=5,
3-seed average, 1000 burn-in, 2000 measurement steps):

  Free flow       (rho=0.05, p=0.3): stopped=0.000, flow=0.234
  Near-transition (rho=0.10, p=0.3): stopped=0.003, flow=0.459
  Onset           (rho=0.12, p=0.3): stopped=0.082, flow=0.463 (peak)
  Canonical jam   (rho=0.15, p=0.3): stopped=0.181, flow=0.457
  Deep jam        (rho=0.30, p=0.3): stopped=0.431, flow=0.393
  Deterministic   (rho=0.15, p=0.0): stopped=0.000, flow=0.750
  Density sat.    (rho=0.80, p=0.0): stopped=0.750, flow=0.200
                                     (trivial pigeonhole, NOT jamming)

Published quantitative anchors:
  - Free-flow mean velocity in the dilute limit: <v> -> v_max - p_slow
    (our code: at rho=0.02, p=0.3, v_max=5 -> <v> = 4.70 ✓).
  - Deterministic p=0 sharp transition at rho_c = 1/(v_max+1) = 1/6
    (Wikipedia, analytic).
  - Peak flow at p=0.3: rho* in [0.10, 0.15] with flow ~ 0.46-0.47
    (Nagel-Schreckenberg 1992; Wikipedia).
  - Spontaneous jams (free-flow -> jammed) appear around rho ~ 0.12 at
    p=0.3 (Bette et al. 2017).

CONSERVATION: total number of cars n_cars is exactly conserved. No births,
deaths, lane changes, overtaking.

model_class: "traffic_1d". Deliberately does NOT contain "ca" or
"excitable" substrings, so P13's placeholder exclusion does not
spuriously fire (same pattern as LV's "predator_prey" and Gray-Scott's
"reaction_diffusion").
"""

from __future__ import annotations

from typing import Any

import numpy as np

from epc.base_model import BaseModel


class NagelSchreckenberg(BaseModel):
    """Nagel-Schreckenberg traffic CA on a one-dimensional ring.

    Parameters
    ----------
    L : int
        Number of cells in the ring (default 1000, matching NS 1992 and
        Bette et al. 2017).
    n_cars : int | None
        Number of cars. If None, computed from ``density * L``.
    density : float
        Car density rho (default 0.15, canonical jam regime). Ignored if
        ``n_cars`` is given.
    v_max : int
        Maximum velocity in cells per step (default 5, matching NS 1992).
    p_slow : float
        Randomization probability (default 0.3, matching NS 1992's
        illustrative parameter choice). Must be in [0, 1].
    init_mode : str
        'uniform' : cars spaced as evenly as possible, all starting at
                    v = v_max.
        'random'  : cars placed at n_cars distinct random cells, each
                    starting with v drawn uniformly from 0..v_max.
        'jammed'  : cars packed at positions 0, 1, ..., n_cars-1 with
                    all v = 0 (worst-case initial condition for measuring
                    jam dissipation).
    seed : int
        RNG seed.

    Notes
    -----
    One call to ``step()`` advances one NS time step (parallel update
    over all cars). Our ``get_timescale`` returns ``L / v_max`` (one
    free-flow traversal of the ring), and the canonical P8 measurement
    uses 1000 burn-in steps + at least 1500 measurement steps.

    Examples
    --------
    >>> m = NagelSchreckenberg(L=1000, density=0.15, p_slow=0.3,
    ...                         v_max=5, seed=42)
    >>> m.run(2000)
    >>> m.state_history[-1]["stopped_fraction"]  # doctest: +SKIP
    0.18...
    """

    def __init__(
        self,
        L: int = 1000,
        n_cars: int | None = None,
        density: float = 0.15,
        v_max: int = 5,
        p_slow: float = 0.3,
        init_mode: str = "uniform",
        seed: int = 42,
    ) -> None:
        super().__init__(seed=seed)

        if L <= 0:
            raise ValueError(f"L must be positive, got {L}")
        if v_max < 1:
            raise ValueError(f"v_max must be >= 1, got {v_max}")
        if not 0.0 <= p_slow <= 1.0:
            raise ValueError(f"p_slow must be in [0, 1], got {p_slow}")
        if init_mode not in ("uniform", "random", "jammed"):
            raise ValueError(
                f"init_mode must be 'uniform', 'random', or 'jammed'; got {init_mode!r}"
            )
        if n_cars is None:
            if not 0.0 < density < 1.0:
                raise ValueError(
                    f"density must be in (0, 1) when n_cars is None, got {density}"
                )
            n_cars = max(1, int(round(density * L)))
        if n_cars < 1 or n_cars > L:
            raise ValueError(
                f"n_cars must be in [1, L]; got n_cars={n_cars}, L={L}"
            )

        self.L = int(L)
        self.n_cars = int(n_cars)
        self.density = self.n_cars / self.L
        self.v_max = int(v_max)
        self.p_slow = float(p_slow)
        self.init_mode = init_mode

        self._positions: np.ndarray | None = None
        self._velocities: np.ndarray | None = None

    # ------------------------------------------------------------------
    # Setup / step
    # ------------------------------------------------------------------

    def setup(self) -> dict[str, Any]:
        N, L = self.n_cars, self.L
        if self.init_mode == "uniform":
            positions = np.linspace(0, L - 1, N, dtype=np.int64) % L
            positions = np.sort(positions).astype(np.int32)
            velocities = np.full(N, self.v_max, dtype=np.int8)
        elif self.init_mode == "random":
            # Choose N distinct random cells, sorted
            chosen = self.rng.choice(L, size=N, replace=False)
            positions = np.sort(chosen).astype(np.int32)
            velocities = self.rng.integers(
                0, self.v_max + 1, size=N, dtype=np.int8
            )
        elif self.init_mode == "jammed":
            positions = np.arange(N, dtype=np.int32)
            velocities = np.zeros(N, dtype=np.int8)
        else:  # pragma: no cover — guarded above
            raise ValueError(self.init_mode)

        self._positions = positions
        self._velocities = velocities

        self._step_count = 0
        self._is_setup = True
        return self._snapshot()

    def step(self) -> dict[str, Any]:
        """One parallel-update NS step."""
        assert self._positions is not None and self._velocities is not None, \
            "call setup() first"

        L = self.L
        N = self.n_cars
        v_max = self.v_max
        positions = self._positions
        velocities = self._velocities

        # Compute gaps: empty cells between car i and car (i+1) on the ring.
        if N == 1:
            gaps = np.array([L - 1], dtype=np.int32)
        else:
            next_pos = np.roll(positions, -1)
            gaps = ((next_pos - positions) % L - 1).astype(np.int32)
            # On the ring, gaps always in [0, L - N].

        # Rule 1: accelerate
        v = np.minimum(velocities.astype(np.int32) + 1, v_max)
        # Rule 2: slow to gap (cars limited by the empty cells ahead)
        v = np.minimum(v, gaps)
        # Rule 3: randomize (only cars with v > 0 can be slowed)
        if self.p_slow > 0.0:
            rands = self.rng.random(N)
            mask = (rands < self.p_slow) & (v > 0)
            v = np.where(mask, v - 1, v)
        # Rule 4: motion
        new_positions = (positions + v) % L
        new_velocities = v.astype(np.int8)

        # After a parallel update, car-order on the ring is preserved
        # (no car can overtake because rule 2 enforces v <= gap). But the
        # index-0 car may wrap past position L-1, rotating the array.
        if N > 1:
            diffs = np.diff(new_positions)
            wrap_idx = np.where(diffs < 0)[0]
            if len(wrap_idx) > 0:
                # Roll so positions[0] is minimum
                k = int(wrap_idx[0]) + 1
                new_positions = np.roll(new_positions, -k)
                new_velocities = np.roll(new_velocities, -k)

        self._positions = new_positions
        self._velocities = new_velocities
        self._step_count += 1
        return self._snapshot()

    # ------------------------------------------------------------------
    # Metadata + timescale
    # ------------------------------------------------------------------

    def get_metadata(self) -> dict[str, Any]:
        return {
            "model_name": "nagel_schreckenberg",
            "model_class": "traffic_1d",  # deliberately not "ca" / "excitable"
            "L": self.L,
            "n_cars": self.n_cars,
            "density": self.density,
            "v_max": self.v_max,
            "p_slow": self.p_slow,
            "init_mode": self.init_mode,
            "boundary": "periodic",
            "interaction_type": "nearest_neighbor_gap",
            "update_mode": "parallel",
            "seed": self.seed,
            "reference": (
                "Nagel & Schreckenberg (1992), J. Phys. I France 2, 2221; "
                "Bette et al. (2017), Phys. Rev. E 95, 012311."
            ),
        }

    def get_timescale(self) -> float:
        """System-intrinsic timescale tau: one ring traversal at v_max.

        tau = L / v_max (in timesteps). At L=1000, v_max=5 this is 200.
        The canonical measurement regime is >= 5 * tau = 1000 burn-in
        steps + >= 7.5 * tau = 1500 measurement steps.
        """
        return float(self.L) / float(self.v_max)

    # ------------------------------------------------------------------
    # Snapshot
    # ------------------------------------------------------------------

    def _snapshot(self) -> dict[str, Any]:
        positions = self._positions
        velocities = self._velocities
        N = self.n_cars
        L = self.L
        if N == 1:
            gaps = np.array([L - 1], dtype=np.int32)
        else:
            next_pos = np.roll(positions, -1)
            gaps = ((next_pos - positions) % L - 1).astype(np.int32)
        mean_v = float(velocities.mean()) if N > 0 else 0.0
        stopped = float((velocities == 0).mean()) if N > 0 else 0.0
        return {
            "positions": positions.copy(),
            "velocities": velocities.copy(),
            "gaps": gaps.copy(),
            "n_cars": N,
            "density": self.density,
            "L": L,
            "v_max": self.v_max,
            "mean_velocity": mean_v,
            "flow": mean_v * self.density,
            "stopped_fraction": stopped,
            "step": self._step_count,
        }
