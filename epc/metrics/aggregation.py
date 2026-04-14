"""P1 — Similarity-driven aggregation metrics.

Computes spatial autocorrelation and clustering statistics from type-labeled
spatial data. Works for both 1D arrays and 2D grids.

Supported state history formats:
  1D: 'cell_types' (list[str]) + 'n' (int)
  2D: 'type_labels_at_pos' (N,) + 'grid_dims' (rows, cols)
      or 'grid' (rows, cols) + 'grid_dims' (rows, cols)
"""

from __future__ import annotations

from typing import Any

import numpy as np
from scipy import ndimage

from epc.base_metric import BaseMetric


def _extract_labels_and_adjacency(
    state: dict[str, Any],
) -> tuple[np.ndarray, list[list[int]], int]:
    """Extract type labels and build adjacency list from state.

    Returns (labels, adjacency_list, n) where:
      labels: 1D array of type labels (int-encoded)
      adjacency_list: adj[i] = [j1, j2, ...] for neighbors of position i
      n: number of positions
    """
    # --- 1D case: cell_types list ---
    if "cell_types" in state and "n" in state:
        types = state["cell_types"]
        n = state["n"]
        # Encode string types as integers
        unique = sorted(set(types))
        type_map = {t: i for i, t in enumerate(unique)}
        labels = np.array([type_map[t] for t in types], dtype=int)
        # 1D adjacency: left and right neighbors
        adj = [[] for _ in range(n)]
        for i in range(n):
            if i > 0:
                adj[i].append(i - 1)
            if i < n - 1:
                adj[i].append(i + 1)
        return labels, adj, n

    # --- 2D case: grid with grid_dims ---
    if "grid_dims" in state:
        rows, cols = state["grid_dims"]
        n = rows * cols
        if "type_labels_at_pos" in state:
            labels = np.asarray(state["type_labels_at_pos"], dtype=int).ravel()
        elif "grid" in state:
            labels = np.asarray(state["grid"], dtype=int).ravel()
        else:
            raise KeyError("Need 'type_labels_at_pos' or 'grid' for 2D")
        # 8-connected adjacency
        adj = [[] for _ in range(n)]
        for r in range(rows):
            for c in range(cols):
                i = r * cols + c
                for dr in (-1, 0, 1):
                    for dc in (-1, 0, 1):
                        if dr == 0 and dc == 0:
                            continue
                        nr, nc = r + dr, c + dc
                        if 0 <= nr < rows and 0 <= nc < cols:
                            adj[i].append(nr * cols + nc)
        return labels, adj, n

    raise KeyError("State must contain ('cell_types' + 'n') or ('grid_dims' + labels)")


class MoransI(BaseMetric):
    """Moran's I spatial autocorrelation for categorical data.

    Works for both 1D and 2D spatial structures.
    """

    def __init__(self) -> None:
        super().__init__(name="morans_i")

    def required_keys(self) -> list[str]:
        return []  # Flexible — checked at runtime

    def compute(
        self,
        state_history: list[dict[str, Any]],
        **kwargs: Any,
    ) -> dict[str, Any]:
        t = kwargs.get("timestep", -1)
        state = state_history[t]
        labels, adj, n = _extract_labels_and_adjacency(state)

        unique_types = np.unique(labels)
        w_sum = sum(len(neighbors) for neighbors in adj)

        per_type = {}
        weights = {}

        for t_label in unique_types:
            x = (labels == t_label).astype(float)
            x_bar = x.mean()
            dev = x - x_bar
            denom = np.sum(dev ** 2)

            if denom == 0 or w_sum == 0:
                per_type[int(t_label)] = 0.0
                weights[int(t_label)] = float(x.sum())
                continue

            numer = 0.0
            for i in range(n):
                for j in adj[i]:
                    numer += dev[i] * dev[j]

            I = (n / w_sum) * numer / denom
            per_type[int(t_label)] = float(I)
            weights[int(t_label)] = float(x.sum())

        total_weight = sum(weights.values())
        if total_weight > 0:
            morans_i = sum(per_type[t] * weights[t] / total_weight for t in per_type)
        else:
            morans_i = 0.0

        e_i = -1.0 / (n - 1)

        return {
            "morans_i": float(morans_i),
            "z_score": float(morans_i - e_i),
            "expected_i": float(e_i),
            "per_type": per_type,
        }


class SegregationIndex(BaseMetric):
    """Segregation index: mean same-type neighbor fraction."""

    def __init__(self) -> None:
        super().__init__(name="segregation_index")

    def required_keys(self) -> list[str]:
        return []

    def compute(
        self,
        state_history: list[dict[str, Any]],
        **kwargs: Any,
    ) -> dict[str, Any]:
        t = kwargs.get("timestep", -1)
        state = state_history[t]
        labels, adj, n = _extract_labels_and_adjacency(state)

        fractions = []
        for i in range(n):
            neighbors = adj[i]
            if not neighbors:
                continue
            same = sum(1 for j in neighbors if labels[j] == labels[i])
            fractions.append(same / len(neighbors))

        seg_idx = float(np.mean(fractions)) if fractions else 0.0
        return {
            "segregation_index": seg_idx,
            "segregation_std": float(np.std(fractions)) if fractions else 0.0,
        }


class ClusterStats(BaseMetric):
    """Cluster count and mean size via connected components per type."""

    def __init__(self) -> None:
        super().__init__(name="cluster_stats")

    def required_keys(self) -> list[str]:
        return []

    def compute(
        self,
        state_history: list[dict[str, Any]],
        **kwargs: Any,
    ) -> dict[str, Any]:
        t = kwargs.get("timestep", -1)
        state = state_history[t]
        labels, adj, n = _extract_labels_and_adjacency(state)

        unique_types = np.unique(labels)
        total_clusters = 0
        all_sizes = []

        for t_label in unique_types:
            mask = (labels == t_label)
            # BFS connected components using adjacency
            visited = np.zeros(n, dtype=bool)
            for start in range(n):
                if not mask[start] or visited[start]:
                    continue
                # BFS from start
                queue = [start]
                visited[start] = True
                size = 0
                while queue:
                    node = queue.pop(0)
                    size += 1
                    for nb in adj[node]:
                        if mask[nb] and not visited[nb]:
                            visited[nb] = True
                            queue.append(nb)
                total_clusters += 1
                all_sizes.append(size)

        return {
            "cluster_count": total_clusters,
            "mean_cluster_size": float(np.mean(all_sizes)) if all_sizes else 0.0,
            "max_cluster_size": int(np.max(all_sizes)) if all_sizes else 0,
            "cluster_sizes": all_sizes,
        }


# --- Null model helpers ---

def label_shuffle_null(
    state: dict[str, Any],
    n_permutations: int = 999,
    rng: np.random.Generator | None = None,
) -> np.ndarray:
    """Generate Moran's I distribution under label-shuffle null.

    Permutes type labels across positions, preserving type proportions.
    Works for both 1D and 2D.
    """
    if rng is None:
        rng = np.random.default_rng(0)

    labels, adj, n = _extract_labels_and_adjacency(state)
    w_sum = sum(len(neighbors) for neighbors in adj)

    unique_types = np.unique(labels)
    null_values = np.empty(n_permutations)

    for perm in range(n_permutations):
        shuffled = rng.permutation(labels)

        per_type_i = {}
        type_weights = {}
        for t_label in unique_types:
            x = (shuffled == t_label).astype(float)
            x_bar = x.mean()
            dev = x - x_bar
            denom = np.sum(dev ** 2)
            if denom == 0 or w_sum == 0:
                per_type_i[int(t_label)] = 0.0
                type_weights[int(t_label)] = float(x.sum())
                continue
            numer = sum(dev[i] * dev[j] for i in range(n) for j in adj[i])
            I = (n / w_sum) * numer / denom
            per_type_i[int(t_label)] = float(I)
            type_weights[int(t_label)] = float(x.sum())

        tw = sum(type_weights.values())
        if tw > 0:
            null_values[perm] = sum(
                per_type_i[t] * type_weights[t] / tw for t in per_type_i
            )
        else:
            null_values[perm] = 0.0

    return null_values
