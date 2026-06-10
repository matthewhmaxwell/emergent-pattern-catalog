"""
Territorial random walks with scent marking and decay.

Two independent implementations for T1b cross-model generalization:

1. ScentMarkingModel — Agents on a 2D torus deposit their own-ID scent
   (a decaying field); movement avoids cells bearing strong FOREIGN scent
   and is attracted to own-scent cells. Emergent: stable, non-overlapping
   home ranges with persistent boundaries. Canonical model following
   Giuggioli, Potts, & Harris (2011).

2. PheromoneRepulsionModel — Simplified pheromone-based territory model
   with exponential decay and hard avoidance threshold. Provides an
   independent implementation for T1b cross-model generalization.

Both models return agent positions + per-agent scent fields over time
via simulate(), deterministic from seed.

Reference: Giuggioli, L., Potts, J. R., & Harris, S. (2011).
Animal interactions and the emergence of territoriality.
PLoS Computational Biology, 7(3), e1002008.
"""

from __future__ import annotations

from dataclasses import dataclass, field as dc_field
from typing import Any, Dict, List, Optional

import numpy as np


@dataclass
class ScentMarkingParams:
    """Parameters for the scent-marking territorial model.

    Agents deposit scent on their current cell at each step; scent decays
    exponentially. Movement is biased: foreign scent repels, own scent
    weakly attracts (home fidelity).
    """
    n_agents: int = 4
    grid_size: int = 48          # L x L torus
    deposition_rate: float = 0.1   # scent deposited per step
    decay_rate: float = 0.03       # exponential decay per step (fraction removed)
    repulsion_strength: float = 2.0  # foreign scent rejection ratio threshold
    home_attraction: float = 2.0    # weight of own-scent attraction
    temperature: float = 0.5       # softmax temperature (higher = more random)
    n_steps: int = 20000
    snapshot_interval: int = 100   # record state every N steps
    seed: int = 42


@dataclass
class PheromoneRepulsionParams:
    """Parameters for the pheromone-repulsion territorial model (T1b)."""
    n_agents: int = 4
    grid_size: int = 48
    deposition_rate: float = 5.0
    decay_rate: float = 0.005
    avoidance_threshold: float = 0.3  # hard threshold: don't enter if foreign > this
    n_steps: int = 5000
    snapshot_interval: int = 50
    seed: int = 42


def _neighbors_torus(r: int, c: int, L: int) -> List[tuple]:
    """Return von Neumann neighbors on an L x L torus (+ stay)."""
    return [
        (r, c),  # stay
        ((r - 1) % L, c),  # up
        ((r + 1) % L, c),  # down
        (r, (c - 1) % L),  # left
        (r, (c + 1) % L),  # right
    ]


class ScentMarkingModel:
    """Territorial random walks with scent marking (Giuggioli et al., 2011).

    N agents on an L x L torus. Each agent deposits own-ID scent that decays
    over time. Movement uses a two-stage rule inspired by Giuggioli's
    boundary-reflection mechanism:

    Stage 1 — Avoidance: Among the von Neumann neighbors (+ stay), exclude
    cells where foreign scent exceeds repulsion_strength × own scent (the
    foreign-dominance threshold). This creates hard territorial boundaries.

    Stage 2 — Among allowed cells, move with probability proportional to
    (1 + home_attraction × own_scent), biasing toward familiar territory.

    If all neighbor cells are blocked (rare, boundary deadlock), agent stays.

    Emergent outcome: stable, mutually exclusive home ranges with persistent
    boundaries that partition the torus.

    Parameters
    ----------
    params : ScentMarkingParams
        Model parameters.
    """

    def __init__(self, params: Optional[ScentMarkingParams] = None) -> None:
        if params is None:
            params = ScentMarkingParams()
        self.params = params
        self.rng = np.random.default_rng(params.seed)

    def simulate(self) -> List[Dict[str, Any]]:
        """Run territorial simulation.

        Returns
        -------
        list[dict]
            One dict per snapshot with keys:
            'positions': ndarray(n_agents, 2) — (row, col) per agent
            'scent_fields': ndarray(n_agents, L, L) — per-agent scent field
            'step': int
            'occupancy': ndarray(n_agents, L, L) — cumulative visit counts
        """
        p = self.params
        L = p.grid_size
        N = p.n_agents

        scent = np.zeros((N, L, L), dtype=np.float64)
        occupancy = np.zeros((N, L, L), dtype=np.float64)

        # Spread agents evenly across grid
        positions = np.zeros((N, 2), dtype=int)
        side = int(np.ceil(np.sqrt(N)))
        spacing = L // side
        for i in range(N):
            positions[i, 0] = ((i // side) * spacing + spacing // 2) % L
            positions[i, 1] = ((i % side) * spacing + spacing // 2) % L

        history: List[Dict[str, Any]] = []

        for step in range(p.n_steps):
            # Decay all scent fields
            scent *= (1.0 - p.decay_rate)

            # Each agent deposits scent
            for i in range(N):
                r, c = positions[i]
                scent[i, r, c] += p.deposition_rate
                occupancy[i, r, c] += 1.0

            # Move agents via two-stage avoidance + attraction
            for i in range(N):
                r, c = positions[i]
                nbrs = _neighbors_torus(r, c, L)

                # Stage 1: Filter — exclude cells where foreign scent
                # exceeds a detection threshold. The threshold is set by
                # repulsion_strength × deposition_rate; cells with active
                # foreign presence are blocked, but stale marks (decayed
                # below threshold) become available, creating dynamic
                # territory expansion at the frontier.
                threshold = p.repulsion_strength * p.deposition_rate
                allowed = []
                weights = []
                for nr, nc in nbrs:
                    own = scent[i, nr, nc]
                    foreign = sum(
                        scent[k, nr, nc] for k in range(N) if k != i
                    )
                    if foreign <= threshold:
                        allowed.append((nr, nc))
                        # Stage 2: weight by own-scent familiarity
                        w = 1.0 + p.home_attraction * own
                        weights.append(w)

                if not allowed:
                    # Boundary deadlock — stay in place
                    continue

                # Normalize weights and choose
                warr = np.array(weights)
                # Add temperature-based noise
                warr = warr ** (1.0 / max(p.temperature, 0.01))
                warr /= warr.sum()

                choice = self.rng.choice(len(allowed), p=warr)
                positions[i] = allowed[choice]

            # Record snapshot
            if step % p.snapshot_interval == 0 or step == p.n_steps - 1:
                history.append({
                    'positions': positions.copy(),
                    'scent_fields': scent.copy(),
                    'occupancy': occupancy.copy(),
                    'step': step,
                    'n_agents': N,
                    'grid_size': L,
                })

        return history

    def get_metadata(self) -> Dict[str, Any]:
        """Return model metadata."""
        p = self.params
        return {
            'model_name': 'scent_marking_territory',
            'model_class': 'territorial_random_walk',
            'n_agents': p.n_agents,
            'grid_size': p.grid_size,
            'deposition_rate': p.deposition_rate,
            'decay_rate': p.decay_rate,
            'repulsion_strength': p.repulsion_strength,
            'home_attraction': p.home_attraction,
            'temperature': p.temperature,
            'n_steps': p.n_steps,
            'seed': p.seed,
            'substrate': 'territorial_agent_field',
            'has_scent_marking': True,
            'has_foreign_avoidance': True,
            'has_home_attraction': True,
            'reference': 'Giuggioli, Potts & Harris (2011)',
        }


class PheromoneRepulsionModel:
    """Pheromone-repulsion territorial model (T1b cross-model).

    Simplified territorial model: agents deposit pheromone and avoid cells
    where foreign pheromone exceeds a hard threshold. When all reachable
    neighbors are blocked, the agent stays in place. Provides an independent
    implementation mechanism (hard threshold vs softmax) for T1b.

    Parameters
    ----------
    params : PheromoneRepulsionParams
        Model parameters.
    """

    def __init__(self, params: Optional[PheromoneRepulsionParams] = None) -> None:
        if params is None:
            params = PheromoneRepulsionParams()
        self.params = params
        self.rng = np.random.default_rng(params.seed)

    def simulate(self) -> List[Dict[str, Any]]:
        """Run pheromone-repulsion simulation.

        Returns
        -------
        list[dict]
            Same schema as ScentMarkingModel.
        """
        p = self.params
        L = p.grid_size
        N = p.n_agents

        scent = np.zeros((N, L, L), dtype=np.float64)
        occupancy = np.zeros((N, L, L), dtype=np.float64)

        positions = np.zeros((N, 2), dtype=int)
        spacing = L // int(np.ceil(np.sqrt(N)))
        for i in range(N):
            row_idx = i // int(np.ceil(np.sqrt(N)))
            col_idx = i % int(np.ceil(np.sqrt(N)))
            positions[i, 0] = (row_idx * spacing + spacing // 2) % L
            positions[i, 1] = (col_idx * spacing + spacing // 2) % L

        history: List[Dict[str, Any]] = []

        for step in range(p.n_steps):
            scent *= (1.0 - p.decay_rate)

            for i in range(N):
                r, c = positions[i]
                scent[i, r, c] += p.deposition_rate
                occupancy[i, r, c] += 1.0

            for i in range(N):
                r, c = positions[i]
                nbrs = _neighbors_torus(r, c, L)

                # Filter: exclude cells with strong foreign scent
                allowed = []
                for nr, nc in nbrs:
                    foreign = sum(scent[k, nr, nc] for k in range(N) if k != i)
                    if foreign < p.avoidance_threshold:
                        allowed.append((nr, nc))

                if not allowed:
                    # Blocked — stay in place
                    continue

                # Uniform random among allowed cells
                choice = self.rng.integers(len(allowed))
                positions[i] = allowed[choice]

            if step % p.snapshot_interval == 0 or step == p.n_steps - 1:
                history.append({
                    'positions': positions.copy(),
                    'scent_fields': scent.copy(),
                    'occupancy': occupancy.copy(),
                    'step': step,
                    'n_agents': N,
                    'grid_size': L,
                })

        return history

    def get_metadata(self) -> Dict[str, Any]:
        """Return model metadata."""
        p = self.params
        return {
            'model_name': 'pheromone_repulsion_territory',
            'model_class': 'territorial_random_walk',
            'n_agents': p.n_agents,
            'grid_size': p.grid_size,
            'deposition_rate': p.deposition_rate,
            'decay_rate': p.decay_rate,
            'avoidance_threshold': p.avoidance_threshold,
            'n_steps': p.n_steps,
            'seed': p.seed,
            'substrate': 'territorial_agent_field',
            'has_scent_marking': True,
            'has_foreign_avoidance': True,
            'has_home_attraction': False,
            'reference': 'pheromone-threshold territorial model',
        }


class PlainRandomWalkModel:
    """Plain random walk with scent but NO avoidance (negative control).

    Agents deposit scent identically but do not respond to it. Used as the
    non-territorial baseline: home ranges overlap freely.

    Parameters
    ----------
    n_agents, grid_size, deposition_rate, decay_rate, n_steps,
    snapshot_interval, seed: same as ScentMarkingModel.
    """

    def __init__(
        self,
        n_agents: int = 4,
        grid_size: int = 64,
        deposition_rate: float = 1.0,
        decay_rate: float = 0.02,
        n_steps: int = 3000,
        snapshot_interval: int = 10,
        seed: int = 42,
    ) -> None:
        self.n_agents = n_agents
        self.grid_size = grid_size
        self.deposition_rate = deposition_rate
        self.decay_rate = decay_rate
        self.n_steps = n_steps
        self.snapshot_interval = snapshot_interval
        self.seed = seed
        self.rng = np.random.default_rng(seed)

    def simulate(self) -> List[Dict[str, Any]]:
        """Run plain random walk (no scent response)."""
        L = self.grid_size
        N = self.n_agents

        scent = np.zeros((N, L, L), dtype=np.float64)
        occupancy = np.zeros((N, L, L), dtype=np.float64)

        positions = np.zeros((N, 2), dtype=int)
        spacing = L // int(np.ceil(np.sqrt(N)))
        for i in range(N):
            row_idx = i // int(np.ceil(np.sqrt(N)))
            col_idx = i % int(np.ceil(np.sqrt(N)))
            positions[i, 0] = (row_idx * spacing + spacing // 2) % L
            positions[i, 1] = (col_idx * spacing + spacing // 2) % L

        history: List[Dict[str, Any]] = []

        for step in range(self.n_steps):
            scent *= (1.0 - self.decay_rate)

            for i in range(N):
                r, c = positions[i]
                scent[i, r, c] += self.deposition_rate
                occupancy[i, r, c] += 1.0

            # Move randomly (no scent response)
            for i in range(N):
                r, c = positions[i]
                nbrs = _neighbors_torus(r, c, L)
                choice = self.rng.integers(len(nbrs))
                positions[i] = nbrs[choice]

            if step % self.snapshot_interval == 0 or step == self.n_steps - 1:
                history.append({
                    'positions': positions.copy(),
                    'scent_fields': scent.copy(),
                    'occupancy': occupancy.copy(),
                    'step': step,
                    'n_agents': N,
                    'grid_size': L,
                })

        return history

    def get_metadata(self) -> Dict[str, Any]:
        """Return model metadata."""
        return {
            'model_name': 'plain_random_walk',
            'model_class': 'random_walk',
            'n_agents': self.n_agents,
            'grid_size': self.grid_size,
            'deposition_rate': self.deposition_rate,
            'decay_rate': self.decay_rate,
            'n_steps': self.n_steps,
            'seed': self.seed,
            'substrate': 'territorial_agent_field',
            'has_scent_marking': True,
            'has_foreign_avoidance': False,
            'has_home_attraction': False,
            'reference': 'plain random walk (negative control)',
        }
