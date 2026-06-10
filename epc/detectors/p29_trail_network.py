"""P29 Trail / network formation detector.

Detects emergent efficient transport networks: agents collectively build
connective pathways linking source/food nodes that approach the efficiency
of the minimum spanning tree (MST) while retaining fault tolerance.

Distinct from P4 (territoriality): P29 builds CONNECTIVE transport networks;
P4 makes EXCLUSIVE domains. From P7 (lane formation): P29 trails are
persistent built paths, not dynamic flow lanes.

Primary metrics:
  - weight_distance_correlation: Spearman rank correlation between edge
    weight and 1/distance. For reinforced networks, short edges accumulate
    more weight (positive correlation). For random traffic, no correlation.
  - network_efficiency: total trail length of the strong-edge network
    relative to the MST. Efficient networks have ratio in [1.0, ~2.0].
  - fault_tolerance: fraction of node pairs remaining connected after
    removing the strongest edge.

Detection tiers:
  Screening:  Positive weight-distance correlation + connected strong-edge
              network + significantly better than null.
  Confirmation: length/MST ratio < 2.0, correlation > 0.3, null p < 0.05.
  Definitive: Correlation > 0.5, length/MST < 2.0, fault_tolerance > 0,
              metadata confirms reinforcement mechanism.

Null model: edge-weight shuffle — randomly permute edge weights across
the complete graph, destroying the weight-distance correlation while
preserving the marginal weight distribution.

T1a observation contract:
  Input: trail-network observation bundle — time series of
  (node_positions, edge_weights, step, n_nodes, grid_size)
  read via extract_observation_bundle().

References:
  Tero, A. et al. (2010). Science, 327(5964), 439-442.
  Deneubourg, J.-L. et al. (1990). J. Insect Behavior, 3(2), 159-168.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional

import numpy as np
from scipy import stats as sp_stats

from epc.detector_result import (
    DetectionTier,
    DetectorResult,
    NullType,
    compute_confidence,
)


# ── T1a: observation-bundle adapter ──────────────────────────────────

def extract_observation_bundle(
    history: List[Dict[str, Any]],
) -> Dict[str, Any]:
    """Extract the trail-network observation bundle from history.

    This adapter is the T1a observation contract for P29. Any system that
    emits a time series with keys 'node_positions', 'edge_weights',
    'step', 'n_nodes', 'grid_size' can be consumed.

    Parameters
    ----------
    history : list[dict]
        State trajectory.

    Returns
    -------
    dict with keys:
        'node_positions': list of ndarray(N, 2)
        'edge_weights': list of ndarray(N, N)
        'pheromone_fields': list of ndarray(G, G)
        'steps': ndarray(T,) int
        'n_nodes': int
        'grid_size': int
    """
    node_positions = [np.asarray(h['node_positions']) for h in history]
    edge_weights = [
        np.asarray(h['edge_weights']) if 'edge_weights' in h
        else np.zeros((h['n_nodes'], h['n_nodes']), dtype=np.float64)
        for h in history
    ]
    pheromone_fields = [
        np.asarray(h['pheromone_field']) if 'pheromone_field' in h
        else np.zeros((h['grid_size'], h['grid_size']), dtype=np.float64)
        for h in history
    ]
    steps = np.array([h['step'] for h in history], dtype=int)
    n_nodes = int(history[0]['n_nodes'])
    grid_size = int(history[0]['grid_size'])
    return {
        'node_positions': node_positions,
        'edge_weights': edge_weights,
        'pheromone_fields': pheromone_fields,
        'steps': steps,
        'n_nodes': n_nodes,
        'grid_size': grid_size,
    }


# ── Helpers ──────────────────────────────────────────────────────────

def _mst_weight(positions: np.ndarray) -> float:
    """Compute MST total weight using Prim's algorithm."""
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


def _connectivity_fraction(adj: np.ndarray, n: int) -> float:
    """Fraction of node pairs connected in the graph."""
    if n <= 1:
        return 1.0
    connected = 0
    total = n * (n - 1) // 2
    for start in range(n):
        visited = {start}
        queue = [start]
        while queue:
            node = queue.pop(0)
            for nb in range(n):
                if adj[node, nb] and nb not in visited:
                    visited.add(nb)
                    queue.append(nb)
        for j in range(start + 1, n):
            if j in visited:
                connected += 1
    return connected / total if total > 0 else 1.0


def _network_metrics(
    edge_weights: np.ndarray,
    node_positions: np.ndarray,
) -> Dict[str, float]:
    """Compute network metrics from edge weights.

    Returns dict with: weight_dist_corr, length_ratio, connectivity,
    fault_tolerance, concentration (Gini), n_strong_edges.
    """
    n = len(node_positions)
    diff = node_positions[:, None, :] - node_positions[None, :, :]
    dist = np.sqrt((diff ** 2).sum(axis=-1))

    # Upper-triangle edge indices
    idx = np.triu_indices(n, k=1)
    weights = edge_weights[idx]
    distances = dist[idx]

    mst_w = _mst_weight(node_positions)

    # Weight-distance correlation
    inv_dist = 1.0 / np.maximum(distances, 1e-10)
    if weights.std() > 0 and inv_dist.std() > 0:
        corr, _ = sp_stats.spearmanr(weights, inv_dist)
    else:
        corr = 0.0

    # Gini coefficient of edge weights (concentration)
    sorted_w = np.sort(weights)
    n_edges = len(sorted_w)
    cum = np.cumsum(sorted_w)
    total = sorted_w.sum()
    if total > 0 and n_edges > 1:
        gini = (2 * np.sum((np.arange(1, n_edges + 1) * sorted_w)) /
                (n_edges * total)) - (n_edges + 1) / n_edges
    else:
        gini = 0.0

    # Strong-edge network: edges above median weight (if > 0)
    positive = weights[weights > 0]
    if len(positive) >= 2:
        threshold = np.median(positive)
    elif len(positive) == 1:
        threshold = positive[0] * 0.5
    else:
        return {
            'weight_dist_corr': 0.0,
            'length_ratio': float('inf'),
            'connectivity': 0.0,
            'fault_tolerance': 0.0,
            'concentration': 0.0,
            'n_strong_edges': 0,
            'mst_weight': float(mst_w),
            'network_length': 0.0,
        }

    adj = (edge_weights > threshold) & (edge_weights > 0)
    np.fill_diagonal(adj, False)

    # Network length
    net_length = 0.0
    n_strong = 0
    for i in range(n):
        for j in range(i + 1, n):
            if adj[i, j]:
                net_length += dist[i, j]
                n_strong += 1

    length_ratio = net_length / mst_w if mst_w > 0 else float('inf')
    connectivity = _connectivity_fraction(adj, n)

    # Fault tolerance: remove strongest edge
    ft = 0.0
    if n_strong > 0:
        upper_w = edge_weights.copy()
        np.fill_diagonal(upper_w, 0)
        max_idx = np.unravel_index(upper_w.argmax(), upper_w.shape)
        adj_ft = adj.copy()
        adj_ft[max_idx[0], max_idx[1]] = False
        adj_ft[max_idx[1], max_idx[0]] = False
        ft = _connectivity_fraction(adj_ft, n)

    return {
        'weight_dist_corr': float(corr),
        'length_ratio': float(length_ratio),
        'connectivity': float(connectivity),
        'fault_tolerance': float(ft),
        'concentration': float(gini),
        'n_strong_edges': n_strong,
        'mst_weight': float(mst_w),
        'network_length': float(net_length),
    }


# ── Detector ─────────────────────────────────────────────────────────

class P29TrailNetworkDetector:
    """P29 Trail / network formation detector.

    Detects emergent efficient transport networks via weight-distance
    correlation and network efficiency relative to the MST.

    Parameters
    ----------
    n_permutations : int
        Number of edge-weight shuffle permutations for null model.
    seed : int or None
        RNG seed for reproducibility.
    """

    def __init__(
        self,
        n_permutations: int = 199,
        seed: int = 42,
    ) -> None:
        self.n_permutations = n_permutations
        self.seed = seed

    def detect(
        self,
        history: List[Dict[str, Any]],
        metadata: Optional[Dict[str, Any]] = None,
    ) -> DetectorResult:
        """Run P29 detection on a trail-network history.

        Parameters
        ----------
        history : list[dict]
            State trajectory with keys per the T1a contract.
        metadata : dict, optional
            Model metadata (enables DEFINITIVE tier).

        Returns
        -------
        DetectorResult
        """
        bundle = extract_observation_bundle(history)
        rng = np.random.default_rng(self.seed)

        n_nodes = bundle['n_nodes']
        node_pos = bundle['node_positions'][-1]
        edge_weights = bundle['edge_weights'][-1]

        has_metadata = metadata is not None

        # Content prerequisite: need positive edge weights
        if edge_weights.max() <= 0:
            return self._no_detection(
                edge_weights, node_pos, has_metadata,
                ['no positive edge weights — no network formed'],
            )

        # Compute observed metrics
        obs = _network_metrics(edge_weights, node_pos)

        # ── Null model: edge-weight shuffle ──────────────────────────
        idx = np.triu_indices(n_nodes, k=1)
        observed_corrs = []
        null_corrs = []

        for _ in range(self.n_permutations):
            # Shuffle weights across edges
            shuffled_ew = np.zeros_like(edge_weights)
            w_vals = edge_weights[idx].copy()
            rng.shuffle(w_vals)
            shuffled_ew[idx] = w_vals
            shuffled_ew = shuffled_ew + shuffled_ew.T

            null_m = _network_metrics(shuffled_ew, node_pos)
            null_corrs.append(null_m['weight_dist_corr'])

        null_corrs = np.array(null_corrs)

        # P-value: fraction of null correlations ≥ observed
        p_value = float(
            (np.sum(null_corrs >= obs['weight_dist_corr']) + 1)
            / (self.n_permutations + 1)
        )

        null_mean = float(np.mean(null_corrs))
        null_std = float(np.std(null_corrs)) if len(null_corrs) > 1 else 1.0
        cohens_d = (
            (obs['weight_dist_corr'] - null_mean) / null_std
            if null_std > 0 else 0.0
        )

        # ── Tier assignment ──────────────────────────────────────────
        warnings: List[str] = []

        # Screening: positive correlation + connected network + p < 0.10
        screening_pass = (
            obs['weight_dist_corr'] > 0.1
            and obs['connectivity'] >= 0.6
            and p_value < 0.10
        )

        # Confirmation: stronger correlation + efficient + p < 0.05
        confirmation_pass = (
            screening_pass
            and obs['weight_dist_corr'] > 0.3
            and obs['length_ratio'] < 2.0
            and p_value < 0.05
            and cohens_d > 1.0
        )

        # Definitive: strong correlation + efficient + metadata
        definitive_pass = (
            confirmation_pass
            and obs['weight_dist_corr'] > 0.5
            and obs['length_ratio'] <= 2.0
            and obs['fault_tolerance'] > 0.0
            and p_value <= 0.005
            and has_metadata
            and metadata.get('has_pheromone_reinforcement', False)
        )

        if not screening_pass:
            return self._no_detection(
                edge_weights, node_pos, has_metadata, warnings,
                obs=obs, p_value=p_value, cohens_d=cohens_d,
                null_mean=null_mean, null_std=null_std,
            )

        if definitive_pass:
            tier = DetectionTier.DEFINITIVE
            bonuses = {
                'all_exclusions_cleared': True,
                'both_null_types_rejected': p_value <= 0.001,
                'finite_size_robustness': obs['fault_tolerance'] >= 0.5,
            }
        elif confirmation_pass:
            tier = DetectionTier.CONFIRMATION
            bonuses = {
                'null_p_0001': p_value < 0.001,
                'effect_size_gt_1': cohens_d > 1.0,
                'all_secondaries': obs['fault_tolerance'] > 0.0,
            }
        else:
            tier = DetectionTier.SCREENING
            bonuses = {
                'secondaries_pass': obs['fault_tolerance'] > 0.0,
                'shuffle_null_p_001': p_value < 0.01,
            }

        confidence = compute_confidence(tier, bonuses)

        return DetectorResult(
            pattern_id='P29',
            detected=True,
            tier=tier,
            confidence=confidence,
            primary_metric={
                'weight_dist_corr': obs['weight_dist_corr'],
                'length_ratio': obs['length_ratio'],
                'connectivity': obs['connectivity'],
                'fault_tolerance': obs['fault_tolerance'],
            },
            secondary_metrics={
                'mst_weight': obs['mst_weight'],
                'network_length': obs['network_length'],
                'n_strong_edges': obs['n_strong_edges'],
                'concentration': obs['concentration'],
                'n_nodes': n_nodes,
            },
            effect_size={
                'cohens_d': float(cohens_d),
                'observed_corr': obs['weight_dist_corr'],
                'null_mean_corr': float(null_mean),
                'null_std_corr': float(null_std),
            },
            null_p_value=float(p_value),
            null_type=NullType.SHUFFLE,
            metadata_available=has_metadata,
            warnings=warnings,
        )

    def _no_detection(
        self,
        edge_weights: np.ndarray,
        node_pos: np.ndarray,
        has_metadata: bool,
        warnings: List[str],
        obs: Optional[Dict[str, float]] = None,
        p_value: float = 1.0,
        cohens_d: float = 0.0,
        null_mean: float = 0.0,
        null_std: float = 0.0,
    ) -> DetectorResult:
        """Return a not-detected result."""
        if obs is None:
            obs = _network_metrics(edge_weights, node_pos)
        return DetectorResult(
            pattern_id='P29',
            detected=False,
            tier=DetectionTier.SCREENING,
            confidence=0.0,
            primary_metric={
                'weight_dist_corr': obs.get('weight_dist_corr', 0.0),
                'length_ratio': obs.get('length_ratio', float('inf')),
                'connectivity': obs.get('connectivity', 0.0),
                'fault_tolerance': obs.get('fault_tolerance', 0.0),
            },
            secondary_metrics={
                'mst_weight': obs.get('mst_weight', 0.0),
                'network_length': obs.get('network_length', 0.0),
                'n_strong_edges': obs.get('n_strong_edges', 0),
                'concentration': obs.get('concentration', 0.0),
                'n_nodes': len(node_pos),
            },
            effect_size={
                'cohens_d': float(cohens_d),
                'observed_corr': obs.get('weight_dist_corr', 0.0),
                'null_mean_corr': float(null_mean),
                'null_std_corr': float(null_std),
            },
            null_p_value=float(p_value),
            null_type=NullType.SHUFFLE,
            metadata_available=has_metadata,
            warnings=warnings,
            notes='P29 not detected.',
        )
