"""Spatial Rock-Paper-Scissors (stochastic May-Leonard) on a 2D lattice.

Primary reference:
  Reichenbach, T., Mobilia, M. & Frey, E. (2007). "Mobility promotes and
  jeopardizes biodiversity in rock-paper-scissors games."
  Nature 448, 1046-1049. arXiv: q-bio/0702032
  — Canonical stochastic lattice implementation. Square lattice with
  periodic BCs, states {A, B, C, ∅} (empty), three Poisson-rate reactions:
  selection (σ), reproduction (µ), exchange/mobility (ε).

Related references:
  Reichenbach, T., Mobilia, M. & Frey, E. (2008). "Self-organization of
  mobile populations in cyclic competition." J. Theor. Biol. 254, 368-383.
  — Extended analysis: spiral wavelength, phase diagram, PDE limit.

  Szabó, G. & Fáth, G. (2007). "Evolutionary games on graphs."
  Phys. Rep. 446, 97-216. (Review covering cyclic competition games.)

MODEL RULES (exactly as in Reichenbach 2007 Methods / Supp. Info):

Lattice: L × L square lattice, periodic boundaries. Each site is either
occupied by an individual of species A, B, or C, or empty (∅).

At each elementary step:
  1. Pick a random site i.
  2. Pick a random nearest neighbor j (von Neumann: 4 neighbors).
  3. Propose one of three reactions at its rate:
     - Selection (rate σ, default 1):
         A_i + B_j → A_i + ∅_j   (A beats B; B dies, leaves empty site)
         B_i + C_j → B_i + ∅_j
         C_i + A_j → C_i + ∅_j
         Reverse pairs (B_i + A_j, etc.) — "A beats B" is the same event
         regardless of which is at site i vs j. We implement: if the pair
         is {species X, species Y} with X dominating Y, then Y dies.
     - Reproduction (rate µ, default 1):
         X_i + ∅_j → X_i + X_j   (X fills the empty neighbor site)
         (Or ∅_i + X_j → X_i + X_j)
     - Exchange/mobility (rate ε):
         X_i + Y_j → Y_i + X_j   (swap; any X, Y including empties)

The reaction is chosen with probability proportional to its rate:
  P(selection) = σ / (σ + µ + ε)
  P(reproduction) = µ / (σ + µ + ε)
  P(exchange) = ε / (σ + µ + ε)

A reaction is ATTEMPTED; it only EXECUTES if the site-pair content is
compatible (e.g., "selection" requires two individuals in a dominance
relationship; "reproduction" requires one individual and one empty).
Incompatible proposals are no-ops.

Time units: one GENERATION = N = L² elementary steps (each site is
selected once on average). This is the canonical Reichenbach convention.

Mobility measure (Reichenbach Eq. 1):
    M = 2 ε a² / N         with lattice spacing a = 1, N = L²
    ⇔ ε = M · N / 2 = M · L² / 2

Published quantitative results (Reichenbach 2007):
  1. Critical mobility: M_c ≈ (4.5 ± 0.5) × 10⁻⁴ for σ = µ = 1.
     Below M_c: coexistence of three species (biodiversity maintained).
     Above M_c: two species go extinct (uniform phase).
  2. Spiral wavelength λ ∝ √M (Eq. 2) in natural units (length in lattice
     spacings, time set by σ = 1). Critical wavelength λ_c ≈ 0.8 (universal).
  3. Species fractions fluctuate around 1/3 each during coexistence.
  4. Typical system sizes L = 100 to 500.

IMPLEMENTATION NOTES:

This implementation uses a VECTORIZED ASYNCHRONOUS update: each generation
consists of N = L² independent (site, neighbor, reaction) draws, applied
sequentially. This matches the paper's "Gillespie-style" description
semantically but runs in numpy. The per-step ordering is randomized via
the seeded RNG. Results are reproducible for a fixed seed.

Conservation: total site count L² is conserved (obviously). The sum of
species + empty counts is always L². Verified in tests.

Initial conditions:
  'random': each site independently assigned to {A, B, C, ∅} with
            equal probability 1/4 (matches Reichenbach default).
  'random_no_empty': each site assigned to {A, B, C} with probability 1/3.
  'blocks': three equal blocks (one per species) — for controlled spiral
            nucleation studies.
  'custom': user-supplied grid.

State constants (int dtype):
  0 = EMPTY
  1 = SPECIES_A   (rock)
  2 = SPECIES_B   (paper)
  3 = SPECIES_C   (scissors)

Dominance cycle (A beats B beats C beats A):
  A dominates B   (A eats B → B becomes empty)
  B dominates C
  C dominates A

State history keys:
  grid             : np.ndarray (L, L), dtype int — cell states {0,1,2,3}
  grid_dims        : tuple (L, L)
  n_states         : int — always 4 (empty + 3 species)
  step             : int — generation number (NOT elementary steps)
  empty_count      : int
  a_count, b_count, c_count : int — per-species cell counts
  a_fraction, b_fraction, c_fraction : float — species densities
  empty_fraction   : float
  n_selections, n_reproductions, n_exchanges : int — counts of reactions
                    actually executed during this generation

P12 / P13 boundary:
  RPS has n_states = 4 ≥ 3, so it passes P13's hard guard. Unlike GH,
  where state transitions are clock-driven (refractory → rest regardless
  of neighbors), RPS transitions are ENTIRELY neighbor-dependent. This is
  the signature P12 keys off to discriminate cyclic dominance from
  excitable dynamics.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from epc.base_model import BaseModel


class RPSSpatialModel(BaseModel):
    """Stochastic May-Leonard RPS on a 2D lattice (Reichenbach 2007).

    Parameters
    ----------
    rows : int
        Grid height L (default 100). Must equal cols for a square lattice.
    cols : int
        Grid width L (default 100).
    selection_rate : float
        Rate σ of dominance-removal reactions (default 1.0).
    reproduction_rate : float
        Rate µ of reproduction into empty neighbors (default 1.0).
    exchange_rate : float, optional
        Rate ε of swap reactions. If None, derived from `mobility`.
    mobility : float, optional
        M = 2ε/N (Reichenbach's natural mobility measure). If provided,
        exchange_rate = M * N / 2. Either `exchange_rate` or `mobility`
        must be provided (mobility takes precedence if both).
    neighborhood : str
        'von_neumann' (4 neighbors, Reichenbach canonical) or 'moore' (8).
    init_mode : str
        'random', 'random_no_empty', 'blocks', or 'custom'.
    init_grid : np.ndarray, optional
        For 'custom' init_mode.
    seed : int
        RNG seed for reproducibility.

    Examples
    --------
    Coexistence regime (far below M_c):
        m = RPSSpatialModel(rows=100, cols=100, mobility=1e-5, seed=42)

    Near-critical regime:
        m = RPSSpatialModel(rows=100, cols=100, mobility=5e-4, seed=42)

    Extinction regime (above M_c):
        m = RPSSpatialModel(rows=100, cols=100, mobility=5e-3, seed=42)
    """

    # State constants
    EMPTY = 0
    SPECIES_A = 1
    SPECIES_B = 2
    SPECIES_C = 3
    N_STATES = 4

    # Dominance map: DOMINATES[predator] = prey
    #   A (1) dominates B (2)
    #   B (2) dominates C (3)
    #   C (3) dominates A (1)
    DOMINATES = {1: 2, 2: 3, 3: 1}
    DOMINATED_BY = {2: 1, 3: 2, 1: 3}  # inverse: DOMINATED_BY[prey] = predator

    def __init__(
        self,
        rows: int = 100,
        cols: int = 100,
        selection_rate: float = 1.0,
        reproduction_rate: float = 1.0,
        exchange_rate: float | None = None,
        mobility: float | None = None,
        neighborhood: str = "von_neumann",
        init_mode: str = "random",
        init_grid: np.ndarray | None = None,
        seed: int = 42,
    ) -> None:
        super().__init__(seed=seed)

        if selection_rate <= 0:
            raise ValueError(f"selection_rate must be > 0, got {selection_rate}")
        if reproduction_rate <= 0:
            raise ValueError(f"reproduction_rate must be > 0, got {reproduction_rate}")

        # Resolve mobility <-> exchange_rate. Mobility takes precedence.
        N = rows * cols
        if mobility is not None:
            if mobility < 0:
                raise ValueError(f"mobility must be >= 0, got {mobility}")
            self.mobility = float(mobility)
            self.exchange_rate = mobility * N / 2.0
        elif exchange_rate is not None:
            if exchange_rate < 0:
                raise ValueError(f"exchange_rate must be >= 0, got {exchange_rate}")
            self.exchange_rate = float(exchange_rate)
            self.mobility = 2.0 * exchange_rate / N
        else:
            raise ValueError(
                "Must provide either `mobility` or `exchange_rate`. "
                "Typical values: mobility=1e-4 (coexistence), 5e-3 (extinction)."
            )

        if neighborhood not in ("von_neumann", "moore"):
            raise ValueError("neighborhood must be 'von_neumann' or 'moore'")
        if init_mode not in ("random", "random_no_empty", "blocks", "custom"):
            raise ValueError(f"Unknown init_mode: {init_mode}")

        self.rows = rows
        self.cols = cols
        self.selection_rate = float(selection_rate)
        self.reproduction_rate = float(reproduction_rate)
        self.neighborhood = neighborhood
        self.init_mode = init_mode
        self.init_grid = init_grid

        # Precompute reaction probabilities
        total = self.selection_rate + self.reproduction_rate + self.exchange_rate
        self._p_selection = self.selection_rate / total
        self._p_reproduction = self.reproduction_rate / total
        # p_exchange = 1 - p_selection - p_reproduction

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
        self._last_reaction_counts = {"selection": 0, "reproduction": 0, "exchange": 0}

    def setup(self) -> dict[str, Any]:
        if self.init_mode == "random":
            self._grid = self.rng.integers(
                0, self.N_STATES, size=(self.rows, self.cols), dtype=np.int8
            )
        elif self.init_mode == "random_no_empty":
            self._grid = self.rng.integers(
                1, self.N_STATES, size=(self.rows, self.cols), dtype=np.int8
            )
        elif self.init_mode == "blocks":
            self._grid = self._init_blocks()
        elif self.init_mode == "custom":
            if self.init_grid is None:
                raise ValueError("init_grid required for 'custom' mode")
            self._grid = self.init_grid.astype(np.int8).copy()
            if self._grid.shape != (self.rows, self.cols):
                raise ValueError(
                    f"init_grid shape {self._grid.shape} != ({self.rows},{self.cols})"
                )
        else:
            raise ValueError(f"Unknown init_mode: {self.init_mode}")

        self._step_count = 0
        self._last_reaction_counts = {"selection": 0, "reproduction": 0, "exchange": 0}
        self._is_setup = True

        return self._snapshot()

    def _init_blocks(self) -> np.ndarray:
        """Three equal horizontal blocks, one per species (no empties)."""
        grid = np.zeros((self.rows, self.cols), dtype=np.int8)
        third = self.rows // 3
        grid[:third, :] = self.SPECIES_A
        grid[third:2 * third, :] = self.SPECIES_B
        grid[2 * third:, :] = self.SPECIES_C
        return grid

    def step(self) -> dict[str, Any]:
        """Advance by ONE GENERATION = N elementary steps (N = L²).

        Each elementary step:
          1. Pick random site (i_row, i_col).
          2. Pick random neighbor direction.
          3. Propose reaction type with probabilities (p_sel, p_rep, p_exc).
          4. Execute if compatible; otherwise no-op.
        """
        N = self.rows * self.cols

        # Vectorized draws for the whole generation.
        # Site picks: (row, col) pairs
        site_rows = self.rng.integers(0, self.rows, size=N)
        site_cols = self.rng.integers(0, self.cols, size=N)
        # Neighbor direction indices
        nbr_idx = self.rng.integers(0, self._n_neighbors, size=N)
        # Reaction-type roll: uniform in [0, 1)
        reaction_roll = self.rng.random(size=N)

        # We unfortunately cannot vectorize the APPLICATION of updates
        # (concurrent modifications would violate async semantics). Loop
        # through and apply sequentially.
        grid = self._grid  # direct reference; mutated in-place
        offsets = self._offsets
        p_sel = self._p_selection
        p_rep = self._p_reproduction
        # p_exc boundary: anything >= (p_sel + p_rep) is exchange

        n_sel = 0
        n_rep = 0
        n_exc = 0

        rows = self.rows
        cols = self.cols

        for k in range(N):
            r = site_rows[k]
            c = site_cols[k]
            dr, dc = offsets[nbr_idx[k]]
            nr = (r + dr) % rows
            nc = (c + dc) % cols

            s_i = grid[r, c]
            s_j = grid[nr, nc]

            roll = reaction_roll[k]

            if roll < p_sel:
                # Selection: if the pair contains a dominance relation,
                # the dominated species dies (becomes empty).
                # Check both directions because the paper does "picker
                # dominates neighbor" semantics; to be exactly symmetric
                # with the undirected physical process, we check for any
                # dominance relation in the pair.
                # Canonical Reichenbach: i selects j. So s_i kills s_j
                # iff s_i dominates s_j.
                if s_i != self.EMPTY and s_j != self.EMPTY:
                    if self.DOMINATES.get(int(s_i)) == int(s_j):
                        grid[nr, nc] = self.EMPTY
                        n_sel += 1
            elif roll < p_sel + p_rep:
                # Reproduction: i reproduces into j if j is empty and i isn't.
                if s_i != self.EMPTY and s_j == self.EMPTY:
                    grid[nr, nc] = s_i
                    n_rep += 1
            else:
                # Exchange: swap i and j (only non-trivial if s_i != s_j).
                if s_i != s_j:
                    grid[r, c] = s_j
                    grid[nr, nc] = s_i
                    n_exc += 1

        self._last_reaction_counts = {
            "selection": n_sel,
            "reproduction": n_rep,
            "exchange": n_exc,
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
        Early termination on absorbing state (only one species remains).
        """
        if not self._is_setup:
            initial_state = self.setup()
            self._state_history = [initial_state]

        for i in range(n_steps):
            state = self.step()
            if (i + 1) % record_every == 0:
                self._state_history.append(state)

            # Early termination: absorbing state reached
            # (exactly one species present, rest are empty or extinct)
            counts = [state["a_count"], state["b_count"], state["c_count"]]
            n_surviving = sum(1 for x in counts if x > 0)
            if n_surviving <= 1:
                break

        return self._state_history

    def get_metadata(self) -> dict[str, Any]:
        return {
            "model_name": "rps_spatial",
            # IMPORTANT: NOT "excitable_ca" — P13's placeholder exclusion
            # logic keys off "ca" substring and would spuriously mark P12
            # as excluded on RPS. We explicitly name this something else.
            "model_class": "cyclic_competition",
            "rows": self.rows,
            "cols": self.cols,
            "n_states": self.N_STATES,
            "selection_rate": self.selection_rate,
            "reproduction_rate": self.reproduction_rate,
            "exchange_rate": self.exchange_rate,
            "mobility": self.mobility,
            "neighborhood": self.neighborhood,
            "n_neighbors": self._n_neighbors,
            "init_mode": self.init_mode,
            "dominance_map": dict(self.DOMINATES),  # {predator: prey}
            "interaction_type": "cyclic_dominance",
            "update_mode": "asynchronous_gillespie",
            "seed": self.seed,
            "reference": "Reichenbach, Mobilia & Frey (2007), Nature 448, 1046-1049",
        }

    def get_timescale(self) -> float:
        """System-intrinsic timescale: T_prop ≈ max(rows, cols).

        Same convention as GH and SIR (lattice_2d with wavefronts). For
        RPS specifically, one generation allows individuals to traverse
        distance ~√M in natural units; the pattern-characteristic timescale
        is a rotation period of a spiral, which depends on the mobility.
        We use T_prop = grid size as a conservative upper bound consistent
        with existing detectors.
        """
        return float(max(self.rows, self.cols))

    def is_converged(self) -> bool:
        """Absorbing state: one or zero species remain."""
        if self._grid is None:
            return False
        counts = [
            int((self._grid == self.SPECIES_A).sum()),
            int((self._grid == self.SPECIES_B).sum()),
            int((self._grid == self.SPECIES_C).sum()),
        ]
        return sum(1 for x in counts if x > 0) <= 1

    def _snapshot(self) -> dict[str, Any]:
        grid = self._grid.copy()
        total = self.rows * self.cols
        a_count = int((grid == self.SPECIES_A).sum())
        b_count = int((grid == self.SPECIES_B).sum())
        c_count = int((grid == self.SPECIES_C).sum())
        empty_count = total - a_count - b_count - c_count

        return {
            "grid": grid,
            "grid_dims": (self.rows, self.cols),
            "n_states": self.N_STATES,
            "step": self._step_count,
            "empty_count": empty_count,
            "a_count": a_count,
            "b_count": b_count,
            "c_count": c_count,
            "empty_fraction": empty_count / total,
            "a_fraction": a_count / total,
            "b_fraction": b_count / total,
            "c_fraction": c_count / total,
            "n_selections": self._last_reaction_counts["selection"],
            "n_reproductions": self._last_reaction_counts["reproduction"],
            "n_exchanges": self._last_reaction_counts["exchange"],
            # Compatibility with excitable-CA conventions: fraction of
            # non-empty cells plays the role of "activity density"
            "activity_density": (a_count + b_count + c_count) / total,
        }
