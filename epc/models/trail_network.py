"""
Trail / network formation models for P29 detection.

Two independent implementations for T1b cross-model generalization:

1. AntTrailModel — Ant colony optimization on a complete graph connecting
   food/source nodes. Edge pheromone reinforces frequently-used routes and
   evaporates elsewhere → a sparse efficient trail network emerges.
   Canonical: Deneubourg et al. (1990); Dorigo & Stützle (2004).

2. PhysarumModel — Flux-reinforcement model on a complete graph.
   Edge conductance grows with flow (positive feedback) and decays otherwise.
   Emergent: a transport network near-optimal relative to the MST.
   Canonical: Tero et al. (2010). Independent T1b cross-model.

3. NoReinforcementModel — Agents choose edges uniformly at random with no
   pheromone response (negative control). No efficient network forms.

All models return node positions + edge weights over time via simulate(),
deterministic from seed.

References:
  Tero, A. et al. (2010). Rules for biologically inspired adaptive network
    design. Science, 327(5964), 439-442.
  Deneubourg, J.-L. et al. (1990). The self-organizing exploratory pattern
    of the Argentine ant. J. Insect Behavior, 3(2), 159-168.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, List, Optional, Tuple

import numpy as np


# ── Utilities ────────────────────────────────────────────────────────

def _mst_weight(positions: np.ndarray) -> float:
    """Compute MST total weight for a set of 2D points using Prim's algorithm.

    Parameters
    ----------
    positions : ndarray(N, 2)
        Node coordinates.

    Returns
    -------
    float
        Total MST edge weight (Euclidean distance).
    """
    n = len(positions)
    if n <= 1:
        return 0.0
    diff = positions[:, None, :] - positions[None, :, :]
    dist = np.sqrt((diff ** 2).sum(axis=-1))

    in_tree = np.zeros(n, dtype=bool)
    min_edge = np.full(n, np.inf)
    min_edge[0] = 0.0
    total = 0.0

    for _ in range(n):
        candidates = np.where(~in_tree)[0]
        u = candidates[np.argmin(min_edge[candidates])]
        total += min_edge[u]
        in_tree[u] = True
        for v in range(n):
            if not in_tree[v] and dist[u, v] < min_edge[v]:
                min_edge[v] = dist[u, v]

    return total


def _place_nodes(n_nodes: int, grid_size: int, rng: np.random.Generator,
                 layout: str = "random") -> np.ndarray:
    """Place food/source nodes on a 2D domain.

    Parameters
    ----------
    n_nodes : int
    grid_size : int
    rng : Generator
    layout : str
        'random' or 'grid'.

    Returns
    -------
    ndarray(n_nodes, 2) of float
    """
    if layout == "grid":
        side = int(np.ceil(np.sqrt(n_nodes)))
        spacing = grid_size / (side + 1)
        coords = []
        for i in range(n_nodes):
            r = ((i // side) + 1) * spacing
            c = ((i % side) + 1) * spacing
            coords.append([r, c])
        return np.array(coords[:n_nodes], dtype=float)
    else:
        margin = grid_size * 0.1
        return rng.uniform(margin, grid_size - margin, size=(n_nodes, 2))


def _pheromone_field_from_edges(
    nodes: np.ndarray,
    edge_weights: np.ndarray,
    grid_size: int,
    n_samples: int = 30,
) -> np.ndarray:
    """Paint edge weights onto a 2D field for visualization/observation."""
    n = len(nodes)
    field = np.zeros((grid_size, grid_size), dtype=np.float64)
    for i in range(n):
        for j in range(i + 1, n):
            w = edge_weights[i, j]
            if w > 0.01:
                for t in np.linspace(0, 1, n_samples):
                    x = nodes[i, 0] * (1 - t) + nodes[j, 0] * t
                    y = nodes[i, 1] * (1 - t) + nodes[j, 1] * t
                    gx = int(np.clip(x, 0, grid_size - 1))
                    gy = int(np.clip(y, 0, grid_size - 1))
                    field[gx, gy] += w
    return field


# ── Model 1: Ant Trail (edge-based ACO) ─────────────────────────────

@dataclass
class AntTrailParams:
    """Parameters for ant trail optimization model.

    n_nodes food/source nodes on a grid_size × grid_size domain.
    n_agents ants travel between nodes, choosing next-edge with probability
    proportional to pheromone^alpha / distance^beta. Pheromone deposited
    on traversed edges; evaporation each step.
    """
    n_nodes: int = 7
    n_agents: int = 40
    grid_size: int = 100
    alpha: float = 2.0              # pheromone importance exponent
    beta: float = 3.0               # distance importance exponent (1/dist^beta)
    deposition_rate: float = 10.0   # pheromone deposited = deposition_rate / edge_dist
    evaporation_rate: float = 0.05  # fraction removed per step
    n_steps: int = 500
    snapshot_interval: int = 10
    node_layout: str = "random"
    seed: int = 42


class AntTrailModel:
    """Ant colony optimization on a complete graph of food/source nodes.

    Agents shuttle between nodes, choosing edges probabilistically based on
    pheromone concentration and inverse distance. After each step, pheromone
    is deposited on traversed edges (proportional to 1/distance) and all
    edges evaporate. Over many iterations, frequently-reinforced short edges
    consolidate into a sparse, efficient network near the MST.

    Parameters
    ----------
    params : AntTrailParams
        Model parameters.
    """

    def __init__(self, params: Optional[AntTrailParams] = None) -> None:
        if params is None:
            params = AntTrailParams()
        self.params = params
        self.rng = np.random.default_rng(params.seed)

    def simulate(self) -> List[Dict[str, Any]]:
        """Run ant trail simulation.

        Returns
        -------
        list[dict]
            One dict per snapshot with keys:
            'node_positions': ndarray(n_nodes, 2)
            'pheromone_field': ndarray(grid_size, grid_size) — visualization
            'edge_weights': ndarray(n_nodes, n_nodes) — pheromone on each edge
            'step': int
            'n_nodes': int
            'grid_size': int
        """
        p = self.params
        rng = self.rng
        n = p.n_nodes

        nodes = _place_nodes(n, p.grid_size, rng, p.node_layout)

        # Distance matrix
        diff = nodes[:, None, :] - nodes[None, :, :]
        dist = np.sqrt((diff ** 2).sum(axis=-1))
        np.fill_diagonal(dist, 1.0)

        # Initial pheromone: small uniform value
        pheromone = np.ones((n, n), dtype=np.float64) * 0.1
        np.fill_diagonal(pheromone, 0.0)

        # Heuristic: 1/distance^beta
        heuristic = (1.0 / dist) ** p.beta
        np.fill_diagonal(heuristic, 0.0)

        # Agent positions: current node index
        agent_node = rng.integers(0, n, size=p.n_agents)

        history: List[Dict[str, Any]] = []

        for step in range(p.n_steps):
            # Evaporate
            pheromone *= (1.0 - p.evaporation_rate)

            # Move each agent: choose next node
            deposit = np.zeros((n, n), dtype=np.float64)

            for a in range(p.n_agents):
                current = agent_node[a]
                # Probability: pheromone^alpha * heuristic
                prob = (pheromone[current] ** p.alpha) * heuristic[current]
                prob[current] = 0.0  # don't stay
                total = prob.sum()
                if total <= 0:
                    # Fallback: uniform random
                    others = [j for j in range(n) if j != current]
                    next_node = rng.choice(others)
                else:
                    prob /= total
                    next_node = rng.choice(n, p=prob)

                # Deposit pheromone proportional to 1/distance
                d = dist[current, next_node]
                dep = p.deposition_rate / d if d > 0 else p.deposition_rate
                deposit[current, next_node] += dep
                deposit[next_node, current] += dep

                agent_node[a] = next_node

            pheromone += deposit
            np.fill_diagonal(pheromone, 0.0)

            # Snapshot
            if step % p.snapshot_interval == 0 or step == p.n_steps - 1:
                field = _pheromone_field_from_edges(
                    nodes, pheromone, p.grid_size
                )
                history.append({
                    'node_positions': nodes.copy(),
                    'pheromone_field': field,
                    'edge_weights': pheromone.copy(),
                    'step': step,
                    'n_nodes': n,
                    'grid_size': p.grid_size,
                })

        return history

    def get_metadata(self) -> Dict[str, Any]:
        """Return model metadata."""
        p = self.params
        return {
            'model_name': 'ant_trail_network',
            'model_class': 'trail_network',
            'n_nodes': p.n_nodes,
            'n_agents': p.n_agents,
            'grid_size': p.grid_size,
            'alpha': p.alpha,
            'beta': p.beta,
            'deposition_rate': p.deposition_rate,
            'evaporation_rate': p.evaporation_rate,
            'n_steps': p.n_steps,
            'seed': p.seed,
            'substrate': 'trail_network',
            'has_pheromone_reinforcement': True,
            'has_evaporation': True,
            'has_multiple_nodes': True,
            'reference': 'Deneubourg et al. (1990); Dorigo & Stützle (2004)',
        }


# ── Model 2: Physarum tube-reinforcement ────────────────────────────

@dataclass
class PhysarumParams:
    """Parameters for Physarum polycephalum flux-reinforcement model.

    Tero et al. (2010) adaptive network: edge conductance grows with flow
    (Murray's law feedback) and decays otherwise. Converges to an efficient
    transport network connecting source/food nodes.
    """
    n_nodes: int = 7
    grid_size: int = 100
    gamma: float = 1.8           # reinforcement exponent (>1 = positive feedback)
    decay_rate: float = 0.01     # conductance decay per step
    n_steps: int = 2000
    snapshot_interval: int = 40
    node_layout: str = "random"
    seed: int = 42


class PhysarumModel:
    """Physarum polycephalum flux-reinforcement network model.

    Implements Tero et al. (2010): a complete graph connecting n_nodes food
    sources. Each edge has a conductance D_ij that evolves by:
        dD_ij/dt = |Q_ij|^gamma - decay_rate * D_ij
    where Q_ij is the flux through the edge. Over time, high-flow edges
    strengthen and low-flow edges decay → efficient network near the MST.

    Parameters
    ----------
    params : PhysarumParams
        Model parameters.
    """

    def __init__(self, params: Optional[PhysarumParams] = None) -> None:
        if params is None:
            params = PhysarumParams()
        self.params = params
        self.rng = np.random.default_rng(params.seed)

    def simulate(self) -> List[Dict[str, Any]]:
        """Run Physarum flux-reinforcement simulation.

        Returns
        -------
        list[dict]
            Same schema as AntTrailModel.
        """
        p = self.params
        rng = self.rng
        n = p.n_nodes

        nodes = _place_nodes(n, p.grid_size, rng, p.node_layout)

        # Initial conductance: proportional to 1/distance
        diff = nodes[:, None, :] - nodes[None, :, :]
        dist = np.sqrt((diff ** 2).sum(axis=-1))
        np.fill_diagonal(dist, 1.0)
        conductance = 1.0 / dist
        np.fill_diagonal(conductance, 0.0)

        history: List[Dict[str, Any]] = []

        for step in range(p.n_steps):
            # Compute total flux across all source-sink pairs
            total_flux = np.zeros((n, n), dtype=np.float64)

            for src in range(n):
                for sink in range(src + 1, n):
                    flux = self._solve_flow(conductance, src, sink, n)
                    total_flux += np.abs(flux)

            # Normalize flux
            flux_max = total_flux.max()
            if flux_max > 0:
                total_flux /= flux_max

            # Update conductance: dD/dt = |Q|^gamma - decay * D
            dD = total_flux ** p.gamma - p.decay_rate * conductance
            conductance += dD * 0.1  # Euler step with dt=0.1
            conductance = np.maximum(conductance, 1e-8)
            np.fill_diagonal(conductance, 0.0)

            # Snapshot
            if step % p.snapshot_interval == 0 or step == p.n_steps - 1:
                field = _pheromone_field_from_edges(
                    nodes, conductance, p.grid_size
                )
                history.append({
                    'node_positions': nodes.copy(),
                    'pheromone_field': field,
                    'edge_weights': conductance.copy(),
                    'step': step,
                    'n_nodes': n,
                    'grid_size': p.grid_size,
                })

        return history

    @staticmethod
    def _solve_flow(
        conductance: np.ndarray, src: int, sink: int, n: int,
    ) -> np.ndarray:
        """Solve network flow for a single source-sink pair.

        Uses Kirchhoff's current law: L * p = b where L is the
        conductance-weighted Laplacian.
        """
        L = -conductance.copy()
        np.fill_diagonal(L, 0.0)
        np.fill_diagonal(L, -L.sum(axis=1))

        b = np.zeros(n, dtype=np.float64)
        b[src] = 1.0
        b[sink] = -1.0

        mask = np.ones(n, dtype=bool)
        mask[sink] = False
        L_red = L[np.ix_(mask, mask)]
        b_red = b[mask]

        try:
            p_red = np.linalg.solve(L_red, b_red)
        except np.linalg.LinAlgError:
            return np.zeros((n, n), dtype=np.float64)

        pressure = np.zeros(n, dtype=np.float64)
        pressure[mask] = p_red

        flux = np.zeros((n, n), dtype=np.float64)
        for i in range(n):
            for j in range(n):
                flux[i, j] = conductance[i, j] * (pressure[i] - pressure[j])

        return flux

    def get_metadata(self) -> Dict[str, Any]:
        """Return model metadata."""
        p = self.params
        return {
            'model_name': 'physarum_network',
            'model_class': 'trail_network',
            'n_nodes': p.n_nodes,
            'grid_size': p.grid_size,
            'gamma': p.gamma,
            'decay_rate': p.decay_rate,
            'n_steps': p.n_steps,
            'seed': p.seed,
            'substrate': 'trail_network',
            'has_pheromone_reinforcement': True,
            'has_evaporation': True,
            'has_multiple_nodes': True,
            'reference': 'Tero et al. (2010)',
        }


# ── Model 3: No-reinforcement negative control ─────────────────────

@dataclass
class NoReinforcementParams:
    """Parameters for no-reinforcement random walk (negative control)."""
    n_nodes: int = 7
    n_agents: int = 40
    grid_size: int = 100
    n_steps: int = 500
    snapshot_interval: int = 10
    node_layout: str = "random"
    seed: int = 42


class NoReinforcementModel:
    """Random walk on a complete graph with no pheromone reinforcement.

    Agents choose the next node uniformly at random (no pheromone bias,
    no distance bias). Deposited pheromone is proportional to 1/distance
    (same as AntTrailModel) but does not influence choices. The resulting
    edge weights reflect uniform random traffic → no efficient network.

    Parameters
    ----------
    params : NoReinforcementParams
        Model parameters.
    """

    def __init__(self, params: Optional[NoReinforcementParams] = None) -> None:
        if params is None:
            params = NoReinforcementParams()
        self.params = params
        self.rng = np.random.default_rng(params.seed)

    def simulate(self) -> List[Dict[str, Any]]:
        """Run random-walk simulation (no reinforcement).

        Returns same schema as AntTrailModel.
        """
        p = self.params
        rng = self.rng
        n = p.n_nodes

        nodes = _place_nodes(n, p.grid_size, rng, p.node_layout)

        diff = nodes[:, None, :] - nodes[None, :, :]
        dist = np.sqrt((diff ** 2).sum(axis=-1))
        np.fill_diagonal(dist, 1.0)

        # Edge weights: accumulate traffic (no reinforcement)
        edge_weights = np.zeros((n, n), dtype=np.float64)
        agent_node = rng.integers(0, n, size=p.n_agents)

        history: List[Dict[str, Any]] = []

        for step in range(p.n_steps):
            # Mild decay so weights reflect recent traffic
            edge_weights *= 0.98

            for a in range(p.n_agents):
                current = agent_node[a]
                # Uniform random choice (no pheromone, no distance bias)
                others = [j for j in range(n) if j != current]
                next_node = rng.choice(others)

                # Uniform deposit — no distance bias (negative control)
                edge_weights[current, next_node] += 1.0
                edge_weights[next_node, current] += 1.0

                agent_node[a] = next_node

            if step % p.snapshot_interval == 0 or step == p.n_steps - 1:
                field = _pheromone_field_from_edges(
                    nodes, edge_weights, p.grid_size
                )
                history.append({
                    'node_positions': nodes.copy(),
                    'pheromone_field': field,
                    'edge_weights': edge_weights.copy(),
                    'step': step,
                    'n_nodes': n,
                    'grid_size': p.grid_size,
                })

        return history

    def get_metadata(self) -> Dict[str, Any]:
        """Return model metadata."""
        p = self.params
        return {
            'model_name': 'no_reinforcement_walk',
            'model_class': 'random_walk',
            'n_nodes': p.n_nodes,
            'n_agents': p.n_agents,
            'grid_size': p.grid_size,
            'n_steps': p.n_steps,
            'seed': p.seed,
            'substrate': 'trail_network',
            'has_pheromone_reinforcement': False,
            'has_evaporation': False,
            'has_multiple_nodes': True,
            'reference': 'random walk negative control',
        }
