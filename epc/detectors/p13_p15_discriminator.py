"""P13/P15 Transfer Entropy discriminator.

Two-stage test from detector_cards.md:
  Stage 1 (TE test): Compute boundary-conditioned Transfer Entropy.
    TE ≈ 0 at state boundaries → P13 (annihilative collisions).
    TE > 0 at state boundaries → P15 candidate.
  Stage 2 (Functional test): Do collision outcomes vary with input?
    (Deferred — requires collision detection and outcome cataloguing.)

The discriminator uses boundary-conditioned TE: TE is computed only at
cells with heterogeneous neighborhoods (where different states meet).
This isolates the signal at wavefronts and structure boundaries.

Key result: GH excitable waves have boundary TE ≈ 0.001 bits (resting
cells become excited deterministically when any neighbor is excited —
minimal per-neighbor information). GoL structures have boundary TE ≈
0.010 bits (birth/survival depends on neighbor count, creating richer
per-neighbor predictive structure).

Reference: Lizier, J.T., Prokopenko, M. & Zomaya, A.Y. (2007, 2012).
"""

from __future__ import annotations

from typing import Any

import numpy as np


class P13P15Discriminator:
    """Transfer Entropy discriminator for the P13/P15 boundary.

    Computes boundary-conditioned TE for a system and compares against
    a Greenberg-Hastings control (the canonical P13 null).

    Usage:
        disc = P13P15Discriminator()
        result = disc.discriminate(state_history, gh_control_history)
    """

    def __init__(self, n_permutations: int = 99) -> None:
        self.n_permutations = n_permutations

    def discriminate(
        self,
        state_history: list[dict[str, Any]],
        gh_control_history: list[dict[str, Any]] | None = None,
        model_metadata: dict[str, Any] | None = None,
    ) -> dict[str, Any]:
        """Run the P13/P15 TE discrimination test.

        Parameters
        ----------
        state_history : list[dict]
            State trajectory of the system under test.
        gh_control_history : list[dict], optional
            State trajectory from a GH control run (same grid size).
            If not provided, uses a default GH run.
        model_metadata : dict, optional
            Model metadata for context.

        Returns
        -------
        dict with:
            classification : str — 'P13', 'P15_candidate', or 'inconclusive'
            observed_te : float — boundary-conditioned TE of test system
            control_te : float — boundary-conditioned TE of GH control
            te_ratio : float — observed / control
            p_value : float — permutation test p-value
            stage1_pass : bool — TE significantly above control
            note : str
        """
        if not state_history or "grid" not in state_history[0]:
            return {
                "classification": "inconclusive",
                "observed_te": 0.0,
                "control_te": 0.0,
                "te_ratio": 0.0,
                "p_value": 1.0,
                "stage1_pass": False,
                "note": "No grid data available",
            }

        rows, cols = state_history[0]["grid_dims"]

        # Compute observed boundary TE
        obs_te, obs_te_max, obs_n = self._boundary_te(state_history)

        # Compute GH control boundary TE
        if gh_control_history is not None:
            ctrl_te, ctrl_te_max, ctrl_n = self._boundary_te(gh_control_history)
        else:
            ctrl_te = self._run_gh_control(rows, cols, len(state_history))

        # Permutation test: shuffle grids spatially, recompute boundary TE
        rng = np.random.default_rng(0)
        null_tes = []
        # Subsample trajectory for null runs
        subsample = max(1, len(state_history) // 80)
        subsampled = state_history[::subsample]

        for _ in range(self.n_permutations):
            shuffled = []
            for state in subsampled:
                grid = state["grid"].copy()
                flat = grid.ravel()
                rng.shuffle(flat)
                new_state = dict(state)
                new_state["grid"] = flat.reshape(rows, cols)
                shuffled.append(new_state)
            null_te, _, _ = self._boundary_te(shuffled)
            null_tes.append(null_te)

        null_tes = np.array(null_tes)
        null_mean = float(null_tes.mean())
        null_std = float(null_tes.std())

        # P-value: fraction of null TE values >= observed TE
        p_value = float(np.mean(null_tes >= obs_te))
        if p_value == 0:
            p_value = 1.0 / (self.n_permutations + 1)

        # TE ratio vs GH control
        te_ratio = obs_te / ctrl_te if ctrl_te > 0 else float("inf")

        # Classification
        # Stage 1: TE significantly above permutation null AND above GH control
        # Note: permutation p-value uses ≤ threshold (standard convention)
        stage1_pass = (p_value <= 0.05) and (te_ratio > 3.0)

        if not stage1_pass:
            classification = "P13"
            note = (f"TE ratio {te_ratio:.1f}× vs GH control. "
                    f"Below threshold for P15 (need >3× and p<0.05).")
        else:
            classification = "P15_candidate"
            note = (f"TE ratio {te_ratio:.1f}× vs GH control (p={p_value:.4f}). "
                    f"Stage 2 functional test needed for definitive P15.")

        return {
            "classification": classification,
            "observed_te": obs_te,
            "control_te": ctrl_te,
            "te_ratio": te_ratio,
            "p_value": p_value,
            "stage1_pass": stage1_pass,
            "null_te_mean": null_mean,
            "null_te_std": null_std,
            "note": note,
        }

    def _boundary_te(
        self,
        state_history: list[dict[str, Any]],
    ) -> tuple[float, float, int]:
        """Compute TE conditioned on boundary cells (heterogeneous neighborhoods).

        Uses numpy vectorization for boundary detection (10-50× faster than
        the naive Python loop implementation).

        Returns (mean_te, max_te, n_boundary_observations).
        """
        grids = [s["grid"] for s in state_history]
        rows, cols = state_history[0]["grid_dims"]
        directions_te = [(0, 1), (1, 0), (0, -1), (-1, 0)]

        # Frequency tables
        joint_per_dir: dict[int, dict[tuple, int]] = {i: {} for i in range(4)}
        yt_xt_per_dir: dict[int, dict[tuple, int]] = {i: {} for i in range(4)}
        ynext_yt: dict[tuple, int] = {}
        yt_counts: dict[tuple, int] = {}
        total = 0

        for t in range(len(grids) - 1):
            g = grids[t]
            g_next = grids[t + 1]

            # Vectorized boundary detection: cell differs from ANY Moore neighbor
            is_boundary = np.zeros((rows, cols), dtype=bool)
            for dr in (-1, 0, 1):
                for dc in (-1, 0, 1):
                    if dr == 0 and dc == 0:
                        continue
                    shifted = np.roll(np.roll(g, -dr, axis=0), -dc, axis=1)
                    is_boundary |= (g != shifted)

            # Get boundary cell indices
            boundary_rc = np.argwhere(is_boundary)
            if len(boundary_rc) == 0:
                continue

            # Precompute VN neighbor grids
            vn_grids = []
            for dr, dc in directions_te:
                vn_grids.append(np.roll(np.roll(g, -dr, axis=0), -dc, axis=1))

            # Collect frequency counts from boundary cells
            for r, c in boundary_rc:
                y = int(g[r, c])
                y_next = int(g_next[r, c])

                ynyt_key = (y_next, y)
                ynext_yt[ynyt_key] = ynext_yt.get(ynyt_key, 0) + 1
                yt_key = (y,)
                yt_counts[yt_key] = yt_counts.get(yt_key, 0) + 1
                total += 1

                for d_idx in range(4):
                    x = int(vn_grids[d_idx][r, c])
                    jk = (y_next, y, x)
                    joint_per_dir[d_idx][jk] = joint_per_dir[d_idx].get(jk, 0) + 1
                    yxk = (y, x)
                    yt_xt_per_dir[d_idx][yxk] = yt_xt_per_dir[d_idx].get(yxk, 0) + 1

        if total == 0:
            return 0.0, 0.0, 0

        # H(Y_next | Y_curr) for boundary cells
        h_y_yt = 0.0
        for (y_next, y), n in ynext_yt.items():
            n_y = yt_counts[(y,)]
            p = n / total
            p_cond = n / n_y
            if p_cond > 0:
                h_y_yt -= p * np.log2(p_cond)

        # TE per direction
        te_dirs = []
        for d_idx in range(4):
            h_y_yt_xt = 0.0
            for (y_next, y, x), n in joint_per_dir[d_idx].items():
                n_yx = yt_xt_per_dir[d_idx].get((y, x), 0)
                if n_yx == 0:
                    continue
                p = n / total
                p_cond = n / n_yx
                if p_cond > 0:
                    h_y_yt_xt -= p * np.log2(p_cond)

            te = max(0.0, h_y_yt - h_y_yt_xt)
            te_dirs.append(te)

        return float(np.mean(te_dirs)), float(np.max(te_dirs)), total

    def _run_gh_control(self, rows: int, cols: int, n_steps: int) -> float:
        """Run a GH control and return its boundary TE."""
        from epc.models.greenberg_hastings import GreenbergHastings

        gh = GreenbergHastings(
            rows=rows, cols=cols, n_states=5, threshold=1,
            init_mode="broken_wave", neighborhood="moore",
            boundary="periodic", seed=0,
        )
        h = gh.run(n_steps=min(n_steps, 300), vectorized=True)
        te, _, _ = self._boundary_te(h)
        return te
