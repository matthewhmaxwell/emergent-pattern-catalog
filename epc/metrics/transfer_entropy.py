"""Transfer Entropy for spatiotemporal systems.

Implements the local and average Transfer Entropy (TE) for 2D cellular
automata, following Lizier, Prokopenko & Zomaya (2007, 2012).

Reference: Lizier, J.T., Prokopenko, M. & Zomaya, A.Y. (2007).
"Information transfer by particles in cellular automata."
Progress in Artificial Life, LNCS 4828, 49-60.

Also: Schreiber, T. (2000). "Measuring information transfer."
Physical Review Letters, 85(2), 461-464.

Transfer Entropy from source X to target Y with history embeddings:
    TE(X → Y) = H(Y_{t+1} | Y_t^k) - H(Y_{t+1} | Y_t^k, X_t^l)

For discrete CAs with small state alphabets, we use the plug-in
estimator (direct frequency counting). This is exact for deterministic
CAs and well-suited for the small alphabets (2-14 states) used in
our models.

For the P13/P15 discriminator:
    - GH (excitable waves): TE ≈ 0 (waves annihilate on collision,
      output independent of input configuration)
    - GoL (persistent computation): TE > 0 (glider collisions produce
      input-dependent outputs — computation)
"""

from __future__ import annotations

from typing import Any

import numpy as np

from epc.base_metric import BaseMetric


def _conditional_entropy(
    joint_counts: dict[tuple, int],
    condition_counts: dict[tuple, int],
    total: int,
) -> float:
    """H(A | B) from joint and marginal counts.

    H(A|B) = -Σ p(a,b) log p(a|b) = -Σ p(a,b) log [p(a,b) / p(b)]
    """
    if total == 0:
        return 0.0

    h = 0.0
    for key, n_ab in joint_counts.items():
        if n_ab == 0:
            continue
        # Extract the conditioning part (all but the first element)
        b_key = key[1:]
        n_b = condition_counts.get(b_key, 0)
        if n_b == 0:
            continue
        p_ab = n_ab / total
        p_a_given_b = n_ab / n_b
        h -= p_ab * np.log2(p_a_given_b)

    return h


class TransferEntropy(BaseMetric):
    """Average Transfer Entropy between neighboring cells in a 2D grid CA.

    Computes TE(neighbor → target) averaged over all cells and all
    neighbor directions. Uses plug-in (frequency counting) estimator
    with history embedding k=1 for target and l=1 for source.

    For the P13/P15 discriminator:
        TE ≈ 0 → waves annihilate (P13)
        TE > 0 → information-bearing collisions (P15 candidate)
    """

    def __init__(self, k: int = 1, l: int = 1) -> None:
        super().__init__(name="transfer_entropy")
        self.k = k  # target history embedding
        self.l = l  # source history embedding

    def required_keys(self) -> list[str]:
        return ["grid", "grid_dims"]

    def compute(
        self,
        state_history: list[dict[str, Any]],
        **kwargs: Any,
    ) -> dict[str, Any]:
        """Compute average TE across all cells and neighbor directions.

        Parameters
        ----------
        state_history : list[dict]
            Sequence of grid states.
        n_sample_cells : int, optional
            Number of cells to sample (default: min(500, total)).
            Sampling reduces compute cost on large grids.
        directions : list[tuple], optional
            Neighbor directions to measure TE from.
            Default: [(0,1), (1,0), (0,-1), (-1,0)] (VN neighbors).

        Returns
        -------
        dict with:
            transfer_entropy_mean : float — mean TE across cells and directions
            transfer_entropy_max  : float — max TE over directions
            transfer_entropy_per_direction : dict — TE per direction
            n_transitions : int — number of transitions counted
        """
        if len(state_history) < self.k + 2:
            return {
                "transfer_entropy_mean": 0.0,
                "transfer_entropy_max": 0.0,
                "transfer_entropy_per_direction": {},
                "n_transitions": 0,
            }

        rows, cols = state_history[0]["grid_dims"]
        n_states = state_history[0].get("n_states", 2)
        n_steps = len(state_history)

        n_sample = kwargs.get("n_sample_cells", min(500, rows * cols))
        directions = kwargs.get("directions", [(0, 1), (1, 0), (0, -1), (-1, 0)])

        # Sample cells
        rng = np.random.default_rng(0)
        all_indices = np.arange(rows * cols)
        if n_sample < rows * cols:
            sample_idx = rng.choice(all_indices, size=n_sample, replace=False)
        else:
            sample_idx = all_indices
        sample_rc = [(idx // cols, idx % cols) for idx in sample_idx]

        # Preload grids as arrays for fast access
        grids = [state["grid"] for state in state_history]

        te_per_dir = {}

        for dr, dc in directions:
            dir_key = f"({dr},{dc})"

            # Count joint and conditional frequencies
            # Joint: (y_{t+1}, y_t, x_t)
            # Conditional on y_t only: (y_{t+1}, y_t)
            # Conditional on (y_t, x_t): (y_{t+1}, y_t, x_t) [same as joint]
            joint_counts: dict[tuple, int] = {}  # (y_next, y_curr, x_curr)
            yt_counts: dict[tuple, int] = {}     # (y_curr,)
            yt_xt_counts: dict[tuple, int] = {}  # (y_curr, x_curr)
            ynext_yt_counts: dict[tuple, int] = {}  # (y_next, y_curr)

            total = 0

            for t in range(n_steps - 1):
                g_curr = grids[t]
                g_next = grids[t + 1]

                for r, c in sample_rc:
                    y_curr = int(g_curr[r, c])
                    y_next = int(g_next[r, c])
                    nr = (r + dr) % rows
                    nc = (c + dc) % cols
                    x_curr = int(g_curr[nr, nc])

                    # Joint: (y_next, y_curr, x_curr)
                    jk = (y_next, y_curr, x_curr)
                    joint_counts[jk] = joint_counts.get(jk, 0) + 1

                    # y_t marginal
                    yk = (y_curr,)
                    yt_counts[yk] = yt_counts.get(yk, 0) + 1

                    # (y_t, x_t) marginal
                    yxk = (y_curr, x_curr)
                    yt_xt_counts[yxk] = yt_xt_counts.get(yxk, 0) + 1

                    # (y_next, y_t) marginal
                    ynk = (y_next, y_curr)
                    ynext_yt_counts[ynk] = ynext_yt_counts.get(ynk, 0) + 1

                    total += 1

            # TE = H(Y_{t+1} | Y_t) - H(Y_{t+1} | Y_t, X_t)
            # H(Y_{t+1} | Y_t) using (y_next, y_t) and (y_t)
            h_y_given_yt = _conditional_entropy(ynext_yt_counts, yt_counts, total)
            # H(Y_{t+1} | Y_t, X_t) using (y_next, y_t, x_t) and (y_t, x_t)
            h_y_given_yt_xt = _conditional_entropy(joint_counts, yt_xt_counts, total)

            te = h_y_given_yt - h_y_given_yt_xt
            te_per_dir[dir_key] = max(0.0, te)  # TE is non-negative in expectation

        te_values = list(te_per_dir.values())
        te_mean = float(np.mean(te_values)) if te_values else 0.0
        te_max = float(np.max(te_values)) if te_values else 0.0

        return {
            "transfer_entropy_mean": te_mean,
            "transfer_entropy_max": te_max,
            "transfer_entropy_per_direction": te_per_dir,
            "n_transitions": total,
        }


class LocalTransferEntropy(BaseMetric):
    """Local Transfer Entropy at each space-time point.

    Computes the local TE value at individual (cell, timestep) locations.
    Useful for identifying WHERE information transfer is happening
    (e.g., at glider positions and collision sites in GoL).

    Returns a spatial map of time-averaged local TE, which can be
    compared against known structure locations.
    """

    def __init__(self) -> None:
        super().__init__(name="local_transfer_entropy")

    def required_keys(self) -> list[str]:
        return ["grid", "grid_dims"]

    def compute(
        self,
        state_history: list[dict[str, Any]],
        **kwargs: Any,
    ) -> dict[str, Any]:
        """Compute time-averaged local TE map.

        Returns
        -------
        dict with:
            te_map : np.ndarray (rows, cols) — time-averaged local TE
            te_map_max : float — maximum local TE value
            te_map_mean : float — mean local TE across grid
            active_fraction : float — fraction of cells with TE > threshold
        """
        if len(state_history) < 3:
            rows = state_history[0]["grid_dims"][0] if state_history else 1
            cols = state_history[0]["grid_dims"][1] if state_history else 1
            return {
                "te_map": np.zeros((rows, cols)),
                "te_map_max": 0.0,
                "te_map_mean": 0.0,
                "active_fraction": 0.0,
            }

        rows, cols = state_history[0]["grid_dims"]
        n_steps = len(state_history)
        grids = [state["grid"] for state in state_history]

        # Compute conditional probabilities from global frequencies first
        # Then compute local TE at each (r, c, t)

        # Build global frequency tables (same as TransferEntropy)
        directions = [(0, 1), (1, 0), (0, -1), (-1, 0)]

        # For efficiency: compute for all 4 VN neighbors combined
        # Count (y_next, y_curr, n_excited_neighbors)
        joint_counts: dict[tuple, int] = {}
        yt_counts: dict[tuple, int] = {}
        total = 0

        for t in range(n_steps - 1):
            g_curr = grids[t]
            g_next = grids[t + 1]
            for r in range(rows):
                for c in range(cols):
                    y_curr = int(g_curr[r, c])
                    y_next = int(g_next[r, c])
                    # Count excited neighbors (simplified: total neighbor sum)
                    n_exc = 0
                    for dr, dc in directions:
                        nr, nc = (r + dr) % rows, (c + dc) % cols
                        n_exc += int(g_curr[nr, nc])

                    jk = (y_next, y_curr, n_exc)
                    joint_counts[jk] = joint_counts.get(jk, 0) + 1
                    yk = (y_curr,)
                    yt_counts[yk] = yt_counts.get(yk, 0) + 1
                    total += 1

        # Build conditional probability tables
        # p(y_next | y_curr) and p(y_next | y_curr, n_exc)
        p_y_given_yt: dict[tuple, float] = {}
        p_y_given_yt_x: dict[tuple, float] = {}

        yt_x_counts: dict[tuple, int] = {}
        ynext_yt_counts: dict[tuple, int] = {}

        for (y_next, y_curr, n_exc), count in joint_counts.items():
            yx_key = (y_curr, n_exc)
            yt_x_counts[yx_key] = yt_x_counts.get(yx_key, 0) + count
            yn_key = (y_next, y_curr)
            ynext_yt_counts[yn_key] = ynext_yt_counts.get(yn_key, 0) + count

        for (y_next, y_curr, n_exc), count in joint_counts.items():
            yx_key = (y_curr, n_exc)
            p_y_given_yt_x[(y_next, y_curr, n_exc)] = count / yt_x_counts[yx_key]

        for (y_next, y_curr), count in ynext_yt_counts.items():
            yt_key = (y_curr,)
            p_y_given_yt[(y_next, y_curr)] = count / yt_counts[yt_key]

        # Compute local TE at each cell, averaged over time
        te_map = np.zeros((rows, cols))
        te_count = np.zeros((rows, cols))

        for t in range(n_steps - 1):
            g_curr = grids[t]
            g_next = grids[t + 1]
            for r in range(rows):
                for c in range(cols):
                    y_curr = int(g_curr[r, c])
                    y_next = int(g_next[r, c])
                    n_exc = 0
                    for dr, dc in directions:
                        nr, nc = (r + dr) % rows, (c + dc) % cols
                        n_exc += int(g_curr[nr, nc])

                    p_cond = p_y_given_yt_x.get((y_next, y_curr, n_exc), 0)
                    p_marg = p_y_given_yt.get((y_next, y_curr), 0)

                    if p_cond > 0 and p_marg > 0:
                        local_te = np.log2(p_cond / p_marg)
                        te_map[r, c] += local_te
                        te_count[r, c] += 1

        # Average over time
        mask = te_count > 0
        te_map[mask] /= te_count[mask]

        threshold = 0.01  # bits
        active = float((te_map > threshold).sum()) / (rows * cols)

        return {
            "te_map": te_map,
            "te_map_max": float(te_map.max()),
            "te_map_mean": float(te_map.mean()),
            "active_fraction": active,
        }
