"""P12 — Cyclic (intransitive) dominance detector.

Detects intransitive species dominance: a cyclic A→B→C→A structure where
no species globally dominates the others, and the local neighbor structure
DRIVES transitions between species. This is the canonical mechanism in
spatial rock-paper-scissors games (Reichenbach, Mobilia & Frey 2007).

Observable scope: state-history only.
Works for 2D grid models with ≥ 3 distinguishable species (plus optional
empty cells).

DISCRIMINATION FROM NEIGHBORS:

The interesting boundary is P13 (excitable waves). Both RPS and GH
Greenberg-Hastings can produce persistent spiral-like wavefronts on
a 2D lattice with ≥ 3 states. The mechanisms, however, are different:
  - GH: transitions are CLOCK-DRIVEN. Once a cell is in state k > 0, it
    MUST transition to state (k+1) mod n_states at the next step,
    regardless of neighbors.
  - RPS: transitions are NEIGHBOR-DRIVEN. A cell only changes state when
    a neighbor of the specific species that dominates it (predator) is
    present.

This detector keys off the neighbor-conditional replacement ratio:

    ρ(X, Y) = P(cell becomes Y | cell was X AND had a Y-neighbor) /
              P(cell becomes Y | cell was X AND had NO Y-neighbor)

For RPS dominance edges (prey → predator), ρ is large (empirically
≈ 70–200 at coexistence mobility). For GH clock edges, ρ = 1.0 exactly
because the transition is independent of neighbors. This gives a clean
2+ order-of-magnitude separation.

Detection tiers:
  Screening:    Model has ≥ 3 non-empty species; primary score
                log10(min forward ρ) > 1.0 for at least one cyclic triple.
  Confirmation: All three species persist (each > 5% for ≥ 80% of the
                second-half trajectory); primary score > 1.3; null
                p-value < 0.01.
  Definitive:   Confirmation + primary score > 2.0 + null p-value < 0.005
                + cycle direction is stable across the trajectory halves.

Null model: spatial shuffle of the grid at each timestep. This destroys
the neighbor–transition correlation while preserving species marginals.
Under this null, ρ values collapse to ≈ 1.0 (no conditioning effect).
This is a shuffle-type null (moderate strength); a mechanistic null
(e.g., removing dominance from the rule set) requires model modification
and is not attempted here.

Nearest-neighbor exclusions:
  P13: Excitable waves have clock-driven transitions (ρ ≈ 1.0). P12 fires
       only when ρ >> 1. The two are mechanistically distinct.
  P22: Information cascade is single-pass: cells become recovered and
       stay recovered. RPS cells cycle through species repeatedly.
"""

from __future__ import annotations

from itertools import permutations
from typing import Any

import numpy as np

from epc.base_detector import BaseDetector
from epc.detector_result import DetectionTier, NullType


class P12CyclicDominanceDetector(BaseDetector):
    """Detector for P12 — intransitive (cyclic) dominance.

    The primary signal is a neighbor-conditional replacement ratio ρ(X,Y)
    that is large for exactly three edges forming a cyclic triangle
    (A→B→C→A) and small in the reverse direction.
    """

    # Smoothing constant for ρ computation (avoids div-by-zero and
    # stabilizes tiny denominators).
    _EPS = 1e-6

    def __init__(
        self,
        n_permutations: int = 199,
        neighborhood: str = "von_neumann",
        seed: int = 42,
    ) -> None:
        super().__init__(
            pattern_id="P12",
            excluded_patterns=["P13", "P22"],
            allowed_co_occurrences=["P1"],
            observable_scope="state_history_only",
        )
        if neighborhood not in ("von_neumann", "moore"):
            raise ValueError("neighborhood must be 'von_neumann' or 'moore'")
        self.n_permutations = n_permutations
        self.neighborhood = neighborhood
        self._seed = seed

        if neighborhood == "von_neumann":
            self._offsets = [(-1, 0), (1, 0), (0, -1), (0, 1)]
        else:
            self._offsets = [
                (-1, -1), (-1, 0), (-1, 1),
                (0, -1),            (0, 1),
                (1, -1),  (1, 0),  (1, 1),
            ]

        # Cache identified cyclic triple and direction for use by secondaries.
        self._identified_triple: tuple[int, int, int] | None = None
        self._all_species: list[int] = []

    def _estimate_timescale(
        self,
        state_history: list[dict[str, Any]],
        model_metadata: dict[str, Any] | None,
    ) -> float:
        """T_prop ≈ max(rows, cols): same convention as P13/P22 on lattice_2d."""
        if state_history and "grid_dims" in state_history[0]:
            rows, cols = state_history[0]["grid_dims"]
            return float(max(rows, cols))
        return max(1.0, len(state_history) / 10.0)

    def _validate_prerequisites(
        self,
        state_history: list[dict[str, Any]],
        timescale: float,
    ) -> list[str]:
        warnings = super()._validate_prerequisites(state_history, timescale)
        if not state_history:
            return warnings

        s0 = state_history[0]
        if "grid" not in s0:
            warnings.append("no 'grid' key in state history — P12 needs 2D grid")
            return warnings
        if "grid_dims" not in s0:
            warnings.append("no 'grid_dims' — cannot determine grid size")

        # Count distinct non-empty species across the trajectory.
        # We consider state 0 as "empty" (absent) by convention; state 0 is
        # counted as a species only if state 0 is the ONLY non-species state.
        # We actually just need >=3 distinct states observed, to apply the
        # cyclic-triple logic. Empty cells (state 0) in RPS don't participate
        # in the cycle; we filter state 0 out ONLY if we have >=4 states.
        states_present = set()
        for s in state_history:
            if "grid" in s:
                states_present.update(np.unique(s["grid"]).tolist())
        self._all_species = sorted(states_present)

        if len(states_present) < 3:
            warnings.append(
                f"only {len(states_present)} distinct states observed "
                f"(need ≥ 3 for cyclic triple)"
            )

        return warnings

    # -----------------------------------------------------------------
    # Core metric: neighbor-conditional replacement ratio ρ(X, Y)
    # -----------------------------------------------------------------

    def _count_neighbor_transitions(
        self,
        history: list[dict[str, Any]],
        X: int,
        Y: int,
    ) -> tuple[int, int, int, int]:
        """Accumulate counts for ρ(X, Y) across the trajectory.

        Returns
        -------
        (n_X_with_Y_nbr, n_X_without_Y_nbr,
         n_became_Y_with, n_became_Y_without)
        """
        n_with = 0
        n_without = 0
        n_became_with = 0
        n_became_without = 0

        for t in range(len(history) - 1):
            old = history[t]["grid"]
            new = history[t + 1]["grid"]
            rows, cols = old.shape

            # Count Y-neighbors at each cell (at time t).
            nbr_Y = np.zeros((rows, cols), dtype=np.int32)
            for dr, dc in self._offsets:
                shifted = np.roll(
                    np.roll((old == Y).astype(np.int32), -dr, axis=0),
                    -dc,
                    axis=1,
                )
                nbr_Y += shifted
            has_Y_nbr = nbr_Y > 0

            X_mask = (old == X)
            became_Y = (new == Y)

            X_with = X_mask & has_Y_nbr
            X_without = X_mask & ~has_Y_nbr

            n_with += int(X_with.sum())
            n_without += int(X_without.sum())
            n_became_with += int((X_with & became_Y).sum())
            n_became_without += int((X_without & became_Y).sum())

        return n_with, n_without, n_became_with, n_became_without

    def _rho(
        self,
        history: list[dict[str, Any]],
        X: int,
        Y: int,
    ) -> float:
        """Neighbor-conditional replacement ratio for ordered pair (X, Y)."""
        n_w, n_wo, bw, bwo = self._count_neighbor_transitions(history, X, Y)
        p_with = bw / max(n_w, 1)
        p_without = bwo / max(n_wo, 1)
        # Laplace-smoothed ratio
        return (p_with + self._EPS) / (p_without + self._EPS)

    def _candidate_species(self) -> list[int]:
        """Species IDs to consider for cyclic triples.

        Heuristic: if state 0 is present along with ≥ 3 other states, treat
        state 0 as "empty" (excluded). Otherwise treat all observed states.
        """
        all_s = list(self._all_species)
        if 0 in all_s and len(all_s) >= 4:
            return [s for s in all_s if s != 0]
        return all_s

    def _best_intransitive_triple(
        self,
        history: list[dict[str, Any]],
    ) -> tuple[float, tuple[int, int, int] | None, dict[str, float]]:
        """Search species triples for the one with max min-forward-ρ.

        Returns
        -------
        score : float
            log10(min ρ) over the best cyclic triple direction, or -inf
            if fewer than 3 candidate species.
        best_triple : tuple or None
            (a, b, c) meaning b replaces a (prey a → predator b), c replaces b,
            a replaces c. This is the cyclic direction.
        diagnostics : dict
            Includes best_rho_min, best_rho_max, all three ρ values for
            the winning triple.

        Note: candidates are recomputed locally from the PASSED history.
        This lets the helper be called on sub-trajectories (quarters of
        the run) without depending on the full-trajectory species set
        cached in _validate_prerequisites.
        """
        # Local candidate computation: same rule as _candidate_species() but
        # based on only the states observed in the sub-history.
        states_local: set[int] = set()
        for h in history:
            if "grid" in h:
                states_local.update(np.unique(h["grid"]).tolist())
        all_s = sorted(states_local)
        if 0 in all_s and len(all_s) >= 4:
            candidates = [s for s in all_s if s != 0]
        else:
            candidates = all_s

        if len(candidates) < 3:
            return float("-inf"), None, {
                "n_candidate_species": len(candidates),
            }

        best_min_rho = -np.inf
        best_triple: tuple[int, int, int] | None = None
        best_rhos: tuple[float, float, float] = (0.0, 0.0, 0.0)

        # Try all ordered triples from the candidates (not just 3-combinations
        # in case there are 4+ species — e.g., RPS has states {0,1,2,3} and we
        # look at {1,2,3}; but in general-purpose use the detector should find
        # the cyclic triple within any 3+ species set).
        # For efficiency: first find all 3-subsets, then enumerate 6 cyclic
        # orderings per subset.
        from itertools import combinations

        for triple in combinations(candidates, 3):
            # 6 permutations, but each cyclic direction has 3 equivalent
            # rotations (e.g., (a,b,c) ≡ (b,c,a) ≡ (c,a,b)), and there are
            # 2 cyclic directions per triple (forward and reverse). To avoid
            # redundant work, enumerate 2 directions × 1 rotation = 2 per
            # subset.
            a, b, c = triple
            for perm in [(a, b, c), (a, c, b)]:
                x, y, z = perm
                r1 = self._rho(history, x, y)  # y replaces x
                r2 = self._rho(history, y, z)
                r3 = self._rho(history, z, x)
                min_r = min(r1, r2, r3)
                if min_r > best_min_rho:
                    best_min_rho = min_r
                    best_triple = perm
                    best_rhos = (r1, r2, r3)

        score = float(np.log10(max(best_min_rho, 1e-9)))
        diag = {
            "best_rho_min": float(best_min_rho),
            "best_rho_max": float(max(best_rhos)),
            "rho_1": float(best_rhos[0]),
            "rho_2": float(best_rhos[1]),
            "rho_3": float(best_rhos[2]),
        }
        return score, best_triple, diag

    # -----------------------------------------------------------------
    # BaseDetector contract
    # -----------------------------------------------------------------

    def _compute_primary(
        self,
        state_history: list[dict[str, Any]],
        timescale: float,
    ) -> dict[str, float]:
        """Primary metric: log10 of min forward-cycle ρ over best triple."""
        # Use second half of trajectory for steady-state characterization.
        half = max(1, len(state_history) // 2)
        late = state_history[half:]

        if len(late) < 2:
            return {
                "intransitivity_score": 0.0,
                "best_rho_min": 0.0,
                "n_candidate_species": len(self._all_species),
            }

        score, best_triple, diag = self._best_intransitive_triple(late)
        self._identified_triple = best_triple

        return {
            "intransitivity_score": float(score) if np.isfinite(score) else 0.0,
            **diag,
            "n_candidate_species": len(self._candidate_species()),
        }

    def _check_screening(
        self,
        primary_result: dict[str, float],
        timescale: float,
    ) -> bool:
        """Screening: primary score > 1.0 (at least 10× ρ over no-neighbor)."""
        if primary_result.get("n_candidate_species", 0) < 3:
            return False
        return primary_result.get("intransitivity_score", 0.0) > 1.0

    def _compute_secondaries(
        self,
        state_history: list[dict[str, Any]],
        timescale: float,
    ) -> dict[str, Any]:
        """Secondary metrics:
        - Coexistence stability: fraction of second-half snapshots where
          all three triple-members are present at ≥ 5% density.
        - Cycle direction stability: does the best triple/direction found
          on first-half late-transient match the one found on second-half?
        """
        half = max(1, len(state_history) // 2)
        late = state_history[half:]

        if self._identified_triple is None or len(late) < 2:
            return {
                "coexistence_fraction": 0.0,
                "min_species_density": 0.0,
                "direction_stable": False,
            }

        a, b, c = self._identified_triple
        total_cells = (
            state_history[0]["grid"].size
            if state_history and "grid" in state_history[0]
            else 1
        )

        # Fraction of (second-half) snapshots where all three are ≥ 5%.
        n_all_present = 0
        min_frac_overall = 1.0
        for h in late:
            grid = h["grid"]
            fa = float((grid == a).sum()) / total_cells
            fb = float((grid == b).sum()) / total_cells
            fc = float((grid == c).sum()) / total_cells
            mf = min(fa, fb, fc)
            min_frac_overall = min(min_frac_overall, mf)
            if mf >= 0.05:
                n_all_present += 1
        coexistence_fraction = n_all_present / max(len(late), 1)

        # Direction stability: split the second half into two quarters and
        # compare identified triples.
        direction_stable = False
        if len(late) >= 4:
            q3 = late[: len(late) // 2]
            q4 = late[len(late) // 2:]
            _, trip_q3, _ = self._best_intransitive_triple(q3)
            _, trip_q4, _ = self._best_intransitive_triple(q4)
            if trip_q3 is not None and trip_q4 is not None:
                # Consider direction stable if the cyclic *ordering* matches
                # up to rotation (not just equal tuples).
                direction_stable = self._cyclic_equal(trip_q3, trip_q4)

        return {
            "coexistence_fraction": coexistence_fraction,
            "min_species_density": min_frac_overall,
            "direction_stable": bool(direction_stable),
            "identified_triple": list(self._identified_triple)
            if self._identified_triple
            else None,
        }

    @staticmethod
    def _cyclic_equal(
        t1: tuple[int, int, int],
        t2: tuple[int, int, int],
    ) -> bool:
        """True iff t1 and t2 represent the same cyclic ordering."""
        rotations = [t1, (t1[1], t1[2], t1[0]), (t1[2], t1[0], t1[1])]
        return t2 in rotations

    def _run_null_model(
        self,
        state_history: list[dict[str, Any]],
        primary_result: dict[str, float],
        timescale: float,
    ) -> tuple[float, NullType, dict[str, float]]:
        """Null model: spatial shuffle of the grid at each timestep.

        This destroys the neighbor-transition correlation while preserving
        species marginals. Under the null, ρ → 1.0 and the intransitivity
        score → 0.
        """
        if not state_history or "grid" not in state_history[0]:
            return 1.0, NullType.SHUFFLE, {"mean": 0.0, "std": 0.0}

        observed_score = primary_result.get("intransitivity_score", 0.0)
        if primary_result.get("n_candidate_species", 0) < 3:
            return 1.0, NullType.SHUFFLE, {"mean": 0.0, "std": 0.0}

        half = max(1, len(state_history) // 2)
        late = state_history[half:]

        # Subsample for null speed if trajectory is long.
        target_len = min(len(late), 60)
        if len(late) > target_len:
            stride = max(1, len(late) // target_len)
            subsampled = late[::stride]
        else:
            subsampled = late

        rng = np.random.default_rng(self._seed)
        null_scores: list[float] = []

        for _ in range(self.n_permutations):
            shuffled = []
            for h in subsampled:
                g = h["grid"].copy()
                flat = g.ravel()
                rng.shuffle(flat)
                sh = dict(h)
                sh["grid"] = flat.reshape(g.shape)
                shuffled.append(sh)

            score, _, _ = self._best_intransitive_triple(shuffled)
            if not np.isfinite(score):
                score = 0.0
            null_scores.append(score)

        null_arr = np.asarray(null_scores)
        null_mean = float(null_arr.mean())
        null_std = float(null_arr.std())

        # p-value: fraction of null scores ≥ observed (one-sided; we want
        # to detect EXCESS intransitivity over random spatial arrangement).
        ge = int(np.sum(null_arr >= observed_score))
        p = (ge + 1) / (len(null_arr) + 1)  # add-one smoothing

        return float(p), NullType.SHUFFLE, {"mean": null_mean, "std": null_std}

    def _compute_effect_size(
        self,
        primary_result: dict[str, float],
        null_dist_stats: dict[str, float],
    ) -> dict[str, float]:
        """Cohen's d on intransitivity score vs null."""
        if not null_dist_stats or null_dist_stats.get("std", 0.0) == 0.0:
            return {
                "raw_value": primary_result.get("intransitivity_score", 0.0),
                "null_mean": null_dist_stats.get("mean", 0.0),
                "null_std": null_dist_stats.get("std", 0.0),
            }

        observed = primary_result.get("intransitivity_score", 0.0)
        nm = null_dist_stats.get("mean", 0.0)
        ns = null_dist_stats.get("std", 1.0)
        cohens_d = (observed - nm) / ns if ns > 0 else 0.0
        return {
            "cohens_d": cohens_d,
            "raw_value": observed,
            "null_mean": nm,
            "null_std": ns,
        }

    def _check_confirmation(
        self,
        primary_result: dict[str, float],
        secondary_result: dict[str, Any],
        null_p: float,
        timescale: float,
    ) -> bool:
        """Confirmation: score > 1.3 AND coexistence > 80% AND null p < 0.01."""
        if null_p >= 0.01:
            return False
        if primary_result.get("intransitivity_score", 0.0) <= 1.3:
            return False
        if secondary_result.get("coexistence_fraction", 0.0) < 0.8:
            return False
        return True

    def _check_definitive(
        self,
        primary_result: dict[str, float],
        secondary_result: dict[str, Any],
        null_p: float,
        null_type: NullType,
        state_history: list[dict[str, Any]],
        model_metadata: dict[str, Any] | None,
        timescale: float,
    ) -> bool:
        """Definitive: confirmation + score > 2.0 + direction stable + null p < 0.005."""
        if null_p >= 0.005:
            return False
        if primary_result.get("intransitivity_score", 0.0) <= 2.0:
            return False
        if not secondary_result.get("direction_stable", False):
            return False
        return True

    def _all_secondaries_pass(self, secondary_result: dict[str, Any]) -> bool:
        return (
            secondary_result.get("coexistence_fraction", 0.0) >= 0.8
            and secondary_result.get("direction_stable", False)
        )

    def _check_exclusions(
        self,
        state_history: list[dict[str, Any]],
        model_metadata: dict[str, Any] | None,
        timescale: float,
    ) -> tuple[list[str], dict[str, str]]:
        """Check P13 and P22 exclusions.

        P13: Excitable waves have clock-driven transitions (ρ ≈ 1.0). If our
             intransitivity score is large, that implies neighbor-conditional
             transitions dominate — which is incompatible with P13's
             clock-driven mechanism. We exclude P13 if our best_rho_min >> 1.

        P22: Information cascade is single-pass. If our run shows stable
             coexistence (all three species persist in the long run), that
             rules out single-pass cascade dynamics. We exclude P22 based
             on model_metadata if available; otherwise we mark inconclusive.
        """
        checked = ["P13", "P22"]
        results: dict[str, str] = {}

        # P13 exclusion: large neighbor-conditional ratio implies
        # NOT clock-driven. Use the same trajectory-derived metric.
        # We recompute min-ρ for efficiency check.
        half = max(1, len(state_history) // 2)
        late = state_history[half:]
        score, _, diag = self._best_intransitive_triple(late)
        if diag.get("best_rho_min", 1.0) > 10.0:
            results["P13"] = "excluded"
        elif diag.get("best_rho_min", 1.0) < 2.0:
            results["P13"] = "not_excluded"
        else:
            results["P13"] = "inconclusive"

        # P22 exclusion: metadata-guided.
        if model_metadata:
            model_class = model_metadata.get("model_class", "")
            model_name = model_metadata.get("model_name", "")
            if "cyclic_competition" in model_class or "rps" in model_name:
                # Cyclic dynamics → re-entrant, not single-pass.
                results["P22"] = "excluded"
            elif "epidemic" in model_class or "sir" in model_name:
                results["P22"] = "not_excluded"
            else:
                results["P22"] = "inconclusive"
        else:
            results["P22"] = "not_checked"

        return checked, results
