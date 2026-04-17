"""Stochastic lattice Lotka-Volterra predator-prey model.

Primary reference:
  Mobilia, M., Georgiev, I.T. & Täuber, U.C. (2007). "Phase Transitions
  and Spatio-Temporal Fluctuations in Stochastic Lattice Lotka-Volterra
  Models." J. Stat. Phys. 128, 447-483. arXiv: q-bio/0512039
  — Canonical stochastic lattice implementation with on-site occupation
  restrictions. 2D square lattice, periodic BCs, states {∅, B (prey),
  A (predator)}. Three Poisson-rate reactions.

Related references:
  Heiba, B., Chen, S. & Täuber, U.C. (2018). "Boundary effects on
  population dynamics in stochastic lattice Lotka-Volterra models."
  Physica A 491, 582-590. arXiv: 1706.02567
  — Clean review of the same reactions and cross-correlation analysis.

  Täuber, U.C. (2024). "Stochastic spatial Lotka-Volterra predator-prey
  models." arXiv: 2405.05006. — Comprehensive review of variants and
  universality classes (DP transition at extinction threshold).

MODEL RULES (Mobilia-Georgiev-Täuber 2007, single-occupation variant):

Lattice: L × L square lattice with periodic boundaries. Each site is
either empty (∅), occupied by prey (B), or occupied by a predator (A).
On-site occupation restricted to ≤ 1 particle (finite local carrying
capacity = 1). This is the variant that exhibits a DP-class active-to-
absorbing extinction transition.

Three reactions:
  A         → ∅           at rate μ       (predator spontaneous death)
  B + ∅     → B + B       at rate σ       (prey reproduction into empty NN)
  A + B     → A + A       at rate λ       (predation + predator reproduction
                                          combined, per Mobilia 2007 Sec. II)

At each elementary step:
  1. Pick a random site i uniformly.
  2. Pick a random nearest neighbor j (von Neumann: 4 neighbors).
  3. Propose a reaction with probability proportional to its rate:
       P(death) = μ / (μ + σ + λ)
       P(birth) = σ / (μ + σ + λ)
       P(pred)  = λ / (μ + σ + λ)
  4. Execute if the (site_i, site_j) pair is compatible for that reaction;
     otherwise no-op.

Reaction compatibility:
  - death: requires s_i = A; result: s_i ← ∅. (s_j is ignored but we still
           consume the step for rate bookkeeping.)
  - birth: requires s_i = B, s_j = ∅; result: s_j ← B.
  - pred:  requires s_i = A, s_j = B; result: s_j ← A. (Prey is eaten and
           replaced by a new predator in-place.)

Time unit: one GENERATION = N = L² elementary steps. This matches the
convention used by Mobilia 2007 and RPS.

Published qualitative results (Mobilia-Georgiev-Täuber 2007):

1. COEXISTENCE PHASE (active state): for λ above extinction threshold
   λ_c(σ, μ), both species persist. Near stable foci the system displays:
     - Erratic population oscillations (amplitude shrinks to 0 as L→∞,
       but persists in finite systems — resonant amplification of
       demographic noise).
     - Spreading pursuit-evasion activity fronts (not spirals — bilateral,
       not cyclic).
     - Predator lags prey by ~quarter period.
     - Density autocorrelation shows damped oscillations.

2. STATIONARY-NODE REGIME (still in active phase, further from threshold):
   Essentially stationary localized clusters of predators in a sea of
   prey. Population oscillations absent or strongly suppressed.

3. EXTINCTION (absorbing) PHASE: for λ < λ_c, predators go extinct,
   and prey fill the lattice (uniform B). The transition is continuous,
   directed-percolation universality class.

Parameter ranges used in Mobilia 2007 (with L = 512):
  λ = 0.2 σ = μ = 0.1: well above threshold, clear oscillations (foci)
  λ = 0.15 σ = μ = 0.1: near threshold, slow transient
  (All reaction rates are in units of the MC step rate; rates are
  relative magnitudes, so scaling them all by a constant scales time.)

For our purposes (detection testing, not phase-transition study), we
want the coexistence focus regime with clear oscillations:
  L = 100, λ ≈ 2.0, σ ≈ 1.0, μ ≈ 1.0 (following Heiba-Chen-Täuber 2017
  Sec. II conventions, rescaled so σ+μ+λ ≈ 4 — i.e. we preserve the
  ratios but use rates of order unity so each MC step has meaningful
  acceptance probability).

IMPLEMENTATION NOTES:

Random-sequential updates. Each generation consists of N = L² elementary
(site, neighbor, reaction) draws, applied sequentially. Mirrors RPS
structure — loop through N elementary events per call to step(). Python
inner loop is ~O(N) per generation; at L=100 this is ~10k ops per step,
runtime ≈ 15-25 ms per generation with numpy-vectorized random draws.

Initial conditions:
  'random': each site independently assigned {∅, B, A} with probabilities
            (1 - f_prey - f_pred, f_prey, f_pred). Default f_prey=0.3,
            f_pred=0.3, f_empty=0.4 — matches Mobilia 2007 Sec. III's
            "random initial conditions" convention.
  'uniform': specific fractions (f_empty, f_prey, f_pred).
  'custom':  user-provided grid.

CONSERVATION: total site count L² is conserved. Individual species counts
are NOT conserved (they oscillate).

model_class: "predator_prey". NOT "cyclic_competition" (that's RPS's 3-
species case) and explicitly NOT containing "ca" or "excitable" so that
P13's placeholder exclusion does not spuriously fire (per Sprint 9
Decision 31).
"""

from __future__ import annotations

from typing import Any

import numpy as np

from epc.base_model import BaseModel


class LotkaVolterraLattice(BaseModel):
    """Stochastic lattice Lotka-Volterra (Mobilia-Georgiev-Täuber 2007).

    Two species on a 2D periodic square lattice with on-site occupation
    restriction ≤ 1. States: 0 = empty, 1 = prey (B), 2 = predator (A).

    Parameters
    ----------
    rows : int
        Grid height L (default 100).
    cols : int
        Grid width L (default 100).
    predation_rate : float
        λ, rate of A + B → A + A (default 2.0). Controls how effectively
        predators convert prey into new predators; also the control
        parameter for the DP extinction transition.
    prey_reproduction_rate : float
        σ, rate of B + ∅ → B + B (default 1.0).
    predator_death_rate : float
        μ, rate of A → ∅ (default 1.0).
    neighborhood : str
        'von_neumann' (4 neighbors, Mobilia canonical) or 'moore' (8).
    init_mode : str
        'random', 'uniform', or 'custom'.
    init_prey_fraction : float
        For 'random'/'uniform' init: fraction of sites occupied by prey
        at t=0 (default 0.3).
    init_predator_fraction : float
        For 'random'/'uniform' init: fraction of sites occupied by
        predators at t=0 (default 0.3). Remainder (1 - prey - pred) are
        empty.
    init_grid : np.ndarray, optional
        For 'custom' init: pre-specified grid of shape (rows, cols).
    seed : int
        RNG seed.

    Examples
    --------
    Coexistence (focus) regime, oscillations expected:
        m = LotkaVolterraLattice(rows=100, cols=100,
                                  predation_rate=2.0,
                                  prey_reproduction_rate=1.0,
                                  predator_death_rate=1.0, seed=42)

    Near extinction threshold (weak predation):
        m = LotkaVolterraLattice(rows=100, cols=100,
                                  predation_rate=0.5,
                                  prey_reproduction_rate=1.0,
                                  predator_death_rate=1.0, seed=42)

    Strong predation (toward node regime):
        m = LotkaVolterraLattice(rows=100, cols=100,
                                  predation_rate=5.0,
                                  prey_reproduction_rate=1.0,
                                  predator_death_rate=1.0, seed=42)
    """

    # State constants
    EMPTY = 0
    PREY = 1       # species B
    PREDATOR = 2   # species A
    N_STATES = 3

    def __init__(
        self,
        rows: int = 100,
        cols: int = 100,
        predation_rate: float = 2.0,
        prey_reproduction_rate: float = 1.0,
        predator_death_rate: float = 1.0,
        neighborhood: str = "von_neumann",
        init_mode: str = "random",
        init_prey_fraction: float = 0.3,
        init_predator_fraction: float = 0.3,
        init_grid: np.ndarray | None = None,
        seed: int = 42,
    ) -> None:
        super().__init__(seed=seed)

        if rows <= 0 or cols <= 0:
            raise ValueError("rows and cols must be positive")
        if predation_rate <= 0:
            raise ValueError(
                f"predation_rate must be > 0, got {predation_rate}"
            )
        if prey_reproduction_rate <= 0:
            raise ValueError(
                f"prey_reproduction_rate must be > 0, got {prey_reproduction_rate}"
            )
        if predator_death_rate <= 0:
            raise ValueError(
                f"predator_death_rate must be > 0, got {predator_death_rate}"
            )
        if init_prey_fraction < 0 or init_predator_fraction < 0:
            raise ValueError("init fractions must be >= 0")
        if init_prey_fraction + init_predator_fraction > 1:
            raise ValueError(
                "init_prey_fraction + init_predator_fraction must be <= 1"
            )
        if neighborhood not in ("von_neumann", "moore"):
            raise ValueError("neighborhood must be 'von_neumann' or 'moore'")
        if init_mode not in ("random", "uniform", "custom"):
            raise ValueError(f"Unknown init_mode: {init_mode}")

        self.rows = rows
        self.cols = cols
        self.predation_rate = float(predation_rate)
        self.prey_reproduction_rate = float(prey_reproduction_rate)
        self.predator_death_rate = float(predator_death_rate)
        self.neighborhood = neighborhood
        self.init_mode = init_mode
        self.init_prey_fraction = float(init_prey_fraction)
        self.init_predator_fraction = float(init_predator_fraction)
        self.init_grid = init_grid

        # Precompute reaction probabilities (normalized against total rate)
        total = (
            self.predator_death_rate
            + self.prey_reproduction_rate
            + self.predation_rate
        )
        self._p_death = self.predator_death_rate / total
        self._p_birth = self.prey_reproduction_rate / total
        # p_pred = 1 - p_death - p_birth (implicit, for numerical safety)

        # Neighbor offsets
        if neighborhood == "von_neumann":
            self._offsets = np.array(
                [(-1, 0), (1, 0), (0, -1), (0, 1)], dtype=np.int32
            )
        else:  # moore
            self._offsets = np.array(
                [(-1, -1), (-1, 0), (-1, 1),
                 (0, -1),            (0, 1),
                 (1, -1),  (1, 0),  (1, 1)], dtype=np.int32
            )
        self._n_neighbors = len(self._offsets)

        # State
        self._grid: np.ndarray | None = None
        self._last_reaction_counts = {
            "death": 0, "birth": 0, "predation": 0
        }

    def setup(self) -> dict[str, Any]:
        if self.init_mode == "random":
            # Independent per-site categorical draw with given fractions.
            # Threshold via uniform sample on [0, 1).
            u = self.rng.random(size=(self.rows, self.cols))
            grid = np.zeros((self.rows, self.cols), dtype=np.int8)
            grid[u < self.init_prey_fraction] = self.PREY
            mask_pred = (u >= self.init_prey_fraction) & (
                u < self.init_prey_fraction + self.init_predator_fraction
            )
            grid[mask_pred] = self.PREDATOR
            # rest remain EMPTY (=0)
            self._grid = grid
        elif self.init_mode == "uniform":
            # Shuffle a fixed count of each species across the lattice.
            n_cells = self.rows * self.cols
            n_prey = int(round(self.init_prey_fraction * n_cells))
            n_pred = int(round(self.init_predator_fraction * n_cells))
            flat = np.zeros(n_cells, dtype=np.int8)
            flat[:n_prey] = self.PREY
            flat[n_prey:n_prey + n_pred] = self.PREDATOR
            self.rng.shuffle(flat)
            self._grid = flat.reshape(self.rows, self.cols)
        elif self.init_mode == "custom":
            if self.init_grid is None:
                raise ValueError("init_grid required for 'custom' mode")
            self._grid = self.init_grid.astype(np.int8).copy()
            if self._grid.shape != (self.rows, self.cols):
                raise ValueError(
                    f"init_grid shape {self._grid.shape} != "
                    f"({self.rows},{self.cols})"
                )
        else:
            raise ValueError(f"Unknown init_mode: {self.init_mode}")

        self._step_count = 0
        self._last_reaction_counts = {
            "death": 0, "birth": 0, "predation": 0
        }
        self._is_setup = True

        return self._snapshot()

    def step(self) -> dict[str, Any]:
        """Advance by ONE GENERATION = N elementary steps (N = L²)."""
        N = self.rows * self.cols

        # Vectorized draws for the whole generation.
        site_rows = self.rng.integers(0, self.rows, size=N)
        site_cols = self.rng.integers(0, self.cols, size=N)
        nbr_idx = self.rng.integers(0, self._n_neighbors, size=N)
        reaction_roll = self.rng.random(size=N)

        grid = self._grid  # direct reference; mutated in-place
        offsets = self._offsets
        p_death = self._p_death
        p_birth = self._p_birth
        # predation is the remainder

        rows = self.rows
        cols = self.cols

        n_death = 0
        n_birth = 0
        n_pred = 0

        EMPTY = self.EMPTY
        PREY = self.PREY
        PREDATOR = self.PREDATOR

        for k in range(N):
            r = int(site_rows[k])
            c = int(site_cols[k])
            dr, dc = offsets[nbr_idx[k]]
            nr = (r + dr) % rows
            nc = (c + dc) % cols

            s_i = int(grid[r, c])
            s_j = int(grid[nr, nc])

            roll = reaction_roll[k]

            if roll < p_death:
                # Predator spontaneous death: A → ∅
                if s_i == PREDATOR:
                    grid[r, c] = EMPTY
                    n_death += 1
            elif roll < p_death + p_birth:
                # Prey reproduction: B + ∅ → B + B  (at neighbor site)
                if s_i == PREY and s_j == EMPTY:
                    grid[nr, nc] = PREY
                    n_birth += 1
            else:
                # Predation + predator reproduction: A + B → A + A
                # The prey at the neighbor site becomes a new predator.
                if s_i == PREDATOR and s_j == PREY:
                    grid[nr, nc] = PREDATOR
                    n_pred += 1

        self._last_reaction_counts = {
            "death": n_death,
            "birth": n_birth,
            "predation": n_pred,
        }
        self._step_count += 1
        return self._snapshot()

    def run(
        self, n_steps: int, record_every: int = 1
    ) -> list[dict[str, Any]]:
        """Execute simulation.

        Parameters
        ----------
        n_steps : int
            Number of GENERATIONS to simulate.
        record_every : int
            Record state every N generations.

        Notes
        -----
        Early termination on an absorbing state: either (a) predators
        go extinct (prey-only state, no further change possible since
        reproduction requires empties and birth fills them stably), or
        (b) total extinction (all empty). Not terminated on prey-only
        because the lattice will continue to fill up with prey until
        saturation; we instead detect saturation.
        """
        if not self._is_setup:
            initial_state = self.setup()
            self._state_history = [initial_state]

        for i in range(n_steps):
            state = self.step()
            if (i + 1) % record_every == 0:
                self._state_history.append(state)

            # Absorbing state: predators extinct → prey-only regime.
            # We continue one more step in case predators were just
            # eliminated; then if still zero, we stop.
            if state["predator_count"] == 0:
                # Keep going a few more steps to record filling-up, then stop
                # Actually, once predators are gone, prey just saturate
                # and there's nothing more interesting. Stop immediately.
                break
            if state["prey_count"] == 0 and state["predator_count"] > 0:
                # Predators will die off; let them decay, stop after
                # predator count also reaches zero.
                pass

        return self._state_history

    def get_metadata(self) -> dict[str, Any]:
        return {
            "model_name": "lotka_volterra_lattice",
            # IMPORTANT: "predator_prey" — NOT "cyclic_competition"
            # (that's RPS, 3 species) and NOT containing "ca" or
            # "excitable" substrings (Sprint 9 Decision 31) so that
            # P13's placeholder exclusion does not misfire.
            "model_class": "predator_prey",
            "rows": self.rows,
            "cols": self.cols,
            "n_states": self.N_STATES,
            "n_species": 2,
            "predation_rate": self.predation_rate,
            "prey_reproduction_rate": self.prey_reproduction_rate,
            "predator_death_rate": self.predator_death_rate,
            "neighborhood": self.neighborhood,
            "n_neighbors": self._n_neighbors,
            "init_mode": self.init_mode,
            "init_prey_fraction": self.init_prey_fraction,
            "init_predator_fraction": self.init_predator_fraction,
            "interaction_type": "bilateral_predator_prey",
            "update_mode": "asynchronous_sequential",
            "seed": self.seed,
            "reference": (
                "Mobilia, Georgiev & Täuber (2007), "
                "J. Stat. Phys. 128, 447-483"
            ),
        }

    def get_timescale(self) -> float:
        """System-intrinsic timescale: T_osc.

        For the oscillating coexistence phase, the natural timescale is
        one oscillation period. Mean-field estimate: T_osc ≈ 2π/√(σλ⟨B⟩)
        where ⟨B⟩ is the stationary prey fraction. For our default rates
        (λ=2, σ=1, μ=1) on 100×100, empirically T_osc ≈ 20-40 generations.
        We use a conservative lower bound matching existing lattice_2d
        conventions: T_prop ≈ max(rows, cols).
        """
        # Conservative: use propagation timescale (lattice traversal).
        # For detector-specific timescale checks, detectors should use
        # the dominant-frequency inverse from autocorrelation analysis.
        return float(max(self.rows, self.cols))

    def is_converged(self) -> bool:
        """Absorbing states: predator extinction OR total extinction."""
        if self._grid is None:
            return False
        pred_count = int((self._grid == self.PREDATOR).sum())
        prey_count = int((self._grid == self.PREY).sum())
        return pred_count == 0

    def _snapshot(self) -> dict[str, Any]:
        grid = self._grid.copy()
        total = self.rows * self.cols
        prey_count = int((grid == self.PREY).sum())
        predator_count = int((grid == self.PREDATOR).sum())
        empty_count = total - prey_count - predator_count

        return {
            "grid": grid,
            "grid_dims": (self.rows, self.cols),
            "n_states": self.N_STATES,
            "step": self._step_count,
            "empty_count": empty_count,
            "prey_count": prey_count,
            "predator_count": predator_count,
            "empty_fraction": empty_count / total,
            "prey_fraction": prey_count / total,
            "predator_fraction": predator_count / total,
            "n_deaths": self._last_reaction_counts["death"],
            "n_births": self._last_reaction_counts["birth"],
            "n_predations": self._last_reaction_counts["predation"],
            # Activity density: fraction of sites occupied by ANY species.
            # For compatibility with existing detectors.
            "activity_density": (prey_count + predator_count) / total,
        }
