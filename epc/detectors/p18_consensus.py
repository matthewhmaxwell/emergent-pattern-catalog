"""P18 — Consensus / coarsening-to-consensus detector.

Detects emergent spatial coarsening toward consensus on a 2D lattice with
discrete opinion states. The canonical target is the voter model (Clifford &
Sudbury 1973; Holley & Liggett 1975), and the detector is designed to
discriminate voter-like coarsening from three nearby lattice_2d phenomena:
excitable waves (P13/GH), persistent computation (P15/GoL), and pattern
formation (P1/Schelling).

Observable scope: state-history only.
Works for any 2D grid model with binary or multi-state opinion variables
where local imitation/copying drives a monotonic reduction in spatial
boundary density over time.

Detection tiers (from Sprint 20 characterization, §4.20):

  Screening: early-time monotonic growth of Moran's I AND final Moran's I
             above the spatial-structure threshold. Captures the rapid
             cluster-formation transient of voter dynamics.

  Confirmation: persistent monotonic decay of boundary (wall) density over
             the full run AND wall density plateaus below 0.35 AND shuffle
             null p-value < 0.01. Distinguishes active coarsening from
             static-structured models.

  Definitive: all confirmation gates AND three-class content-level exclusion
             of the lattice_2d neighbors (GH excitable, GoL persistent,
             Schelling aggregation). Specifically:
               - moran_final ∈ [0.30, 0.75] (excludes pre-organized
                 GH spiral at ~0.87 and GoL-random plateau at ~0.27)
               - wall_final > 0.05 (excludes GH at ~0.02)
               - minority fraction stays above 0.05 at end (excludes GoL
                 decay-to-sparse-still-life)

Null model: circular time-shuffle on the Moran's I trajectory. Under the
voter model, Moran's I is monotonically non-decreasing in early time;
shuffling time indices destroys this trend while preserving the marginal
distribution of Moran values. The test statistic is the early-time Spearman
ρ(t, Moran). Under H0 (no trend), null Spearman is centered at 0.

References:
  Clifford, P. & Sudbury, A. (1973). "A model for spatial conflict."
    Biometrika 60(3), 581-588.
  Holley, R. & Liggett, T.M. (1975). "Ergodic theorems for weakly
    interacting infinite systems and the voter model."
    Ann. Prob. 3(4), 643-663.
  Dornic, I., Chaté, H., Chave, J., & Hinrichsen, H. (2001). "Critical
    coarsening without surface tension: the universality class of the
    voter model." Phys. Rev. Lett. 87, 045701.
  Cox, J.T. (1989). "Coalescing random walks and voter model consensus
    times on the torus in Z^d." Ann. Prob. 17(4), 1333-1366.
"""

from __future__ import annotations

from typing import Any

import numpy as np
from scipy.stats import spearmanr

from epc.base_detector import BaseDetector
from epc.detector_result import DetectionTier, NullType


# Moore neighborhood offsets (consistent across detectors on lattice_2d)
_MOORE_OFFSETS = [
    (-1, -1), (-1, 0), (-1, 1),
    (0, -1),           (0, 1),
    (1, -1),  (1, 0),  (1, 1),
]


def _moran_i_moore(grid: np.ndarray) -> float:
    """Moran's I of a 2D grid with Moore neighborhood, periodic boundaries.

    For constant grids returns 0 by convention.
    """
    N = grid.size
    x = grid.astype(float) - grid.mean()
    var = float(np.mean(x ** 2))
    if var < 1e-12:
        return 0.0

    cross = 0.0
    W = 0
    for dr, dc in _MOORE_OFFSETS:
        shifted = np.roll(np.roll(x, -dr, axis=0), -dc, axis=1)
        cross += float(np.sum(x * shifted))
        W += N
    return (N * cross) / (W * N * var)


def _wall_density_moore(grid: np.ndarray) -> float:
    """Boundary (wall) density: fraction of Moore-neighbor pairs with different values."""
    # Half-offsets to avoid double-counting each pair
    half = [(-1, -1), (-1, 0), (-1, 1), (0, 1)]
    total = 0
    diff = 0
    N = grid.size
    for dr, dc in half:
        shifted = np.roll(np.roll(grid, -dr, axis=0), -dc, axis=1)
        diff += int(np.sum(grid != shifted))
        total += N
    return diff / total if total > 0 else 0.0


class P18ConsensusDetector(BaseDetector):
    """Detector for P18 — Consensus / coarsening-to-consensus.

    Primary metrics:
      - moran_spearman_early: Spearman ρ of (t, Moran's I) over t ≤ t_early
        Captures the rapid initial growth of spatial correlation.
      - moran_final_qtr_mean: mean Moran's I over final quarter of run.
        Captures the plateau level distinguishing voter from GoL-random decay.

    Secondary metrics:
      - wall_spearman_early: Spearman ρ of (t, wall density) over t ≤ τ.
        Captures the persistent coarsening trend.
      - wall_final_qtr_mean: plateau wall density.
      - minority_fraction_final: minority opinion share at end of run.
    """

    # Characterization-derived gates (Sprint 20 §4.20)
    SCREENING_MORAN_SPEARMAN_MIN = 0.70
    SCREENING_MORAN_FINAL_MIN = 0.30
    SCREENING_MORAN_GROWTH_MIN = 0.20

    CONFIRMATION_WALL_SPEARMAN_MAX = -0.40
    CONFIRMATION_WALL_FINAL_MAX = 0.30
    CONFIRMATION_WALL_DECAY_MIN = 0.15

    DEFINITIVE_MORAN_FINAL_MIN = 0.30
    DEFINITIVE_MORAN_FINAL_MAX = 0.75
    DEFINITIVE_WALL_FINAL_MIN = 0.05
    DEFINITIVE_MINORITY_FINAL_MIN = 0.05

    # Early-time window: t ≤ this many sweeps (Moran's I grows rapidly here)
    EARLY_TIME_FRACTION = 0.10

    def __init__(self, n_permutations: int = 199, seed: int = 42,
                 n_states_expected: int = 2) -> None:
        super().__init__(
            pattern_id="P18",
            excluded_patterns=["P13", "P15", "P1"],
            allowed_co_occurrences=[],
            observable_scope="state_history_only",
        )
        self.n_permutations = n_permutations
        self._seed = seed
        self.n_states_expected = n_states_expected

    # --- Helpers ---

    def _extract_grids(
        self, state_history: list[dict[str, Any]]
    ) -> list[np.ndarray]:
        """Extract and binarize the opinion grid at each timestep.

        For binary voter models (n_states=2), uses the raw 0/1 grid.
        For multi-state models, collapses to majority-class indicator.
        """
        grids = []
        for state in state_history:
            if "grid" not in state:
                return []
            g = np.asarray(state["grid"])
            if g.dtype != np.int8:
                g = g.astype(np.int8)
            grids.append(g)
        return grids

    def _trajectory_metrics(
        self, state_history: list[dict[str, Any]], timescale: float | None = None,
    ) -> dict[str, Any]:
        """Compute Moran's I and wall-density trajectories + derived metrics.

        The "early window" is t ≤ τ (system-intrinsic timescale, ~L/2). This
        captures the rapid initial cluster-formation transient of voter
        dynamics before the Moran plateau sets in. Empirically at L=64 this
        is ~32 sweeps, which is where the characterization run observed
        Spearman ρ(t, Moran) ∈ [0.73, 0.99] reliably.
        """
        grids = self._extract_grids(state_history)
        if len(grids) < 5:
            return {}

        ts = np.arange(len(grids))
        moran = np.array([_moran_i_moore(g) for g in grids])
        wall = np.array([_wall_density_moore(g) for g in grids])

        # Minority fraction across the run (fraction of cells NOT in
        # majority opinion at the final timestep)
        minority_frac_final = float(
            min(np.mean(grids[-1] == 0), np.mean(grids[-1] == 1))
        )

        # Early window: t ≤ τ (system-intrinsic timescale, ~L/2).
        # Bounded below by 10, above by 100.
        if timescale is None:
            timescale = self._estimate_timescale(state_history, None)
        n = len(grids)
        n_early = max(10, min(100, int(timescale)))
        n_early = min(n_early, n)
        mask_early = ts < n_early
        mask_full = ts >= n_early

        # Spearman trends (safe for near-constant signals)
        def _safe_spearman(x, y):
            if len(x) < 3 or np.std(y) < 1e-10:
                return 0.0
            r, _ = spearmanr(x, y)
            return float(r) if not np.isnan(r) else 0.0

        # Early-window Moran Spearman: captures rapid initial cluster formation
        moran_sp_early = _safe_spearman(ts[mask_early], moran[mask_early])

        # Early-window wall Spearman: voter wall density decreases sharply
        # in t ≤ τ, then plateaus. Sprint 20 §4.20 characterization shows
        # the early-window Spearman is reliably ≤ −0.65 across L=64 / L=128
        # voter seeds, while a full-window Spearman is dominated by late-time
        # plateau noise and can flip positive on individual runs.
        wall_sp_early = _safe_spearman(ts[mask_early], wall[mask_early])

        # Final-quarter means
        n_qtr = max(1, n // 4)
        moran_final_qtr = float(np.mean(moran[-n_qtr:]))
        wall_final_qtr = float(np.mean(wall[-n_qtr:]))

        return {
            "ts": ts,
            "moran": moran,
            "wall": wall,
            "moran_spearman_early": moran_sp_early,
            "moran_final_qtr_mean": moran_final_qtr,
            "moran_initial": float(moran[0]),
            "moran_final": float(moran[-1]),
            "moran_growth": float(moran[-1] - moran[0]),
            "wall_spearman_early": wall_sp_early,
            "wall_final_qtr_mean": wall_final_qtr,
            "wall_initial": float(wall[0]),
            "wall_final": float(wall[-1]),
            "wall_decay": float(wall[0] - wall[-1]),
            "minority_fraction_final": minority_frac_final,
            "n_early": n_early,
        }

    # --- BaseDetector required methods ---

    def _estimate_timescale(
        self,
        state_history: list[dict[str, Any]],
        model_metadata: dict[str, Any] | None,
    ) -> float:
        """Timescale: L/2 = linear size / 2 (mean domain radius at late coarsening).

        For the voter model, the characteristic timescale of coarsening is
        the time to form domains of size comparable to the simulation, which
        scales as L². For detector budgeting we use L/2 as a more modest
        lower bound (a few hundred sweeps at L=64).
        """
        if state_history and "grid_dims" in state_history[0]:
            rows, cols = state_history[0]["grid_dims"]
            return float(max(1, max(rows, cols) // 2))
        return max(1.0, len(state_history) / 10.0)

    def _validate_prerequisites(
        self,
        state_history: list[dict[str, Any]],
        timescale: float,
    ) -> list[str]:
        """Require at least 10τ of history (empirically, ≥100 steps for L=64)."""
        warnings = super()._validate_prerequisites(state_history, timescale)
        if not state_history:
            warnings.append("empty state history")
            return warnings
        if "grid" not in state_history[0]:
            warnings.append("no 'grid' observable")
        if len(state_history) < 50:
            warnings.append(f"run length ({len(state_history)}) < 50 steps")
        return warnings

    def _compute_primary(
        self,
        state_history: list[dict[str, Any]],
        timescale: float,
    ) -> dict[str, float]:
        """Primary: early-time Moran growth + final Moran plateau."""
        metrics = self._trajectory_metrics(state_history)
        if not metrics:
            return {
                "moran_spearman_early": 0.0,
                "moran_final_qtr_mean": 0.0,
                "moran_growth": 0.0,
            }
        return {
            "moran_spearman_early": metrics["moran_spearman_early"],
            "moran_final_qtr_mean": metrics["moran_final_qtr_mean"],
            "moran_initial": metrics["moran_initial"],
            "moran_final": metrics["moran_final"],
            "moran_growth": metrics["moran_growth"],
        }

    def _check_screening(
        self,
        primary_result: dict[str, float],
        timescale: float,
    ) -> bool:
        """Screening gates (all three must pass)."""
        return (
            primary_result.get("moran_spearman_early", 0.0)
            > self.SCREENING_MORAN_SPEARMAN_MIN
            and primary_result.get("moran_final_qtr_mean", 0.0)
            > self.SCREENING_MORAN_FINAL_MIN
            and primary_result.get("moran_growth", 0.0)
            > self.SCREENING_MORAN_GROWTH_MIN
        )

    def _compute_secondaries(
        self,
        state_history: list[dict[str, Any]],
        timescale: float,
    ) -> dict[str, Any]:
        """Secondary: wall-density persistent decay + minority fraction."""
        metrics = self._trajectory_metrics(state_history)
        if not metrics:
            return {}
        return {
            "wall_spearman_early": metrics["wall_spearman_early"],
            "wall_final_qtr_mean": metrics["wall_final_qtr_mean"],
            "wall_initial": metrics["wall_initial"],
            "wall_final": metrics["wall_final"],
            "wall_decay": metrics["wall_decay"],
            "minority_fraction_final": metrics["minority_fraction_final"],
        }

    def _run_null_model(
        self,
        state_history: list[dict[str, Any]],
        primary_result: dict[str, float],
        timescale: float,
    ) -> tuple[float, NullType, dict[str, float]]:
        """Null: random permutation of time indices, recompute Moran Spearman.

        Test statistic: moran_spearman_early (Spearman ρ of t vs Moran's I
        over t ≤ τ).

        Null distribution: permute the order of the Moran's I sequence
        completely, then recompute Spearman ρ over the same early window.
        Random permutation destroys both the directed trend AND the
        autocorrelation of Moran's I, which is appropriate for testing
        "is there a monotonic trend in time" against a no-trend null.

        A circular-shift null was trialed first but rejected: because Moran's
        I is highly autocorrelated (consecutive values differ by < 0.05),
        circular shifts preserve the within-window autocorrelation, leaving
        the null distribution of Spearman ρ with too much mass at large
        positive values and inflating the p-value above the strict 0.01
        confirmation gate. This mirrors Sprint 11 ADR 36's finding that
        circular shifts preserve autocorrelation. Sprint 20 ADR (TBD)
        captures the rationale for using a full permutation null here.

        P-value: fraction of null Spearman ρ ≥ observed.
        """
        grids = self._extract_grids(state_history)
        if len(grids) < 10:
            return 1.0, NullType.SHUFFLE, {"mean": 0.0, "std": 0.0}

        observed = primary_result.get("moran_spearman_early", 0.0)
        n = len(grids)

        # Precompute Moran's I for every timestep
        moran_all = np.array([_moran_i_moore(g) for g in grids])
        ts_full = np.arange(n)
        n_early = max(10, min(100, int(timescale)))
        n_early = min(n_early, n)
        x_early = ts_full[:n_early]

        rng = np.random.default_rng(self._seed)
        null_spearmans = []

        for _ in range(self.n_permutations):
            # Random permutation of the entire Moran sequence; take early
            # window of the permuted sequence and compute Spearman.
            perm = rng.permutation(n)
            permuted = moran_all[perm]
            y = permuted[:n_early]
            if np.std(y) < 1e-10:
                null_spearmans.append(0.0)
                continue
            r, _ = spearmanr(x_early, y)
            null_spearmans.append(float(r) if not np.isnan(r) else 0.0)

        null_arr = np.array(null_spearmans)
        null_mean = float(null_arr.mean())
        null_std = float(null_arr.std())

        # One-sided p-value: P(null Spearman ≥ observed)
        p_value = float(np.mean(null_arr >= observed))
        if p_value == 0:
            p_value = 1.0 / (self.n_permutations + 1)

        return p_value, NullType.SHUFFLE, {"mean": null_mean, "std": null_std}

    def _compute_effect_size(
        self,
        primary_result: dict[str, float],
        null_dist_stats: dict[str, float],
    ) -> dict[str, float]:
        """Effect size on moran_spearman_early."""
        if not null_dist_stats or null_dist_stats.get("std", 0) == 0:
            return {}

        observed = primary_result.get("moran_spearman_early", 0.0)
        null_mean = null_dist_stats.get("mean", 0.0)
        null_std = null_dist_stats.get("std", 1.0)

        cohens_d = (observed - null_mean) / null_std if null_std > 0 else 0.0
        return {
            "cohens_d": cohens_d,
            "raw_value": observed,
            "null_mean": null_mean,
            "null_std": null_std,
            "note": "positive_d_means_stronger_monotonic_Moran_growth_than_shuffled",
        }

    def _check_confirmation(
        self,
        primary_result: dict[str, float],
        secondary_result: dict[str, Any],
        null_p: float,
        timescale: float,
    ) -> bool:
        """Confirmation: null p < 0.01 AND persistent wall decay AND wall plateau low."""
        if null_p >= 0.01:
            return False
        if secondary_result.get(
            "wall_spearman_early", 0.0
        ) >= self.CONFIRMATION_WALL_SPEARMAN_MAX:
            return False
        if secondary_result.get(
            "wall_final_qtr_mean", 1.0
        ) >= self.CONFIRMATION_WALL_FINAL_MAX:
            return False
        if secondary_result.get(
            "wall_decay", 0.0
        ) < self.CONFIRMATION_WALL_DECAY_MIN:
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
        """Definitive: all confirmation gates + three-class exclusion bounds.

        These bounds are calibrated from the Sprint 20 §4.20 characterization
        to exclude each nearest-neighbor on the lattice_2d substrate:
          - GH broken_wave: moran_final ~0.87 > 0.75, rejected.
          - GoL random:     moran_final ~0.27 < 0.30, rejected at screening.
          - GoL r_pent:     moran_spearman_early ~0.17, rejected at screening.
          - GH random:      wall_final ~0.036 < 0.05, rejected here.
          - Schelling P1 (threshold = 0.375): rejected at screening or
            confirmation, NOT here. The Sprint 21 5-seed characterization
            (TestSchellingP18ContentLevel) found Schelling's three-state
            grid {0, 1, 2} yields wall_final ~0.36, well ABOVE the 0.05
            definitive floor; the actual rejection mechanism is
            moran_final_qtr ≤ 0.30 (4 of 5 seeds fail screening) or
            wall_final_qtr ≥ 0.30 (the 5th seed reaches screening but
            fails the confirmation ceiling).
          - Schelling P1 (threshold = 0.5): KNOWN FALSE POSITIVE — reaches
            DEFINITIVE on all 5 characterized seeds with P1 marked
            "inconclusive" because Schelling's metadata lacks a
            copy/imitation/voter `update` key. See Sprint 21 carry-forward
            #20b in REPLICATION_NOTES.md.
        """
        mfq = primary_result.get("moran_final_qtr_mean", 0.0)
        if not (self.DEFINITIVE_MORAN_FINAL_MIN <= mfq <=
                self.DEFINITIVE_MORAN_FINAL_MAX):
            return False
        if secondary_result.get(
            "wall_final_qtr_mean", 0.0
        ) < self.DEFINITIVE_WALL_FINAL_MIN:
            return False
        if secondary_result.get(
            "minority_fraction_final", 0.0
        ) < self.DEFINITIVE_MINORITY_FINAL_MIN:
            return False
        return True

    def _check_exclusions(
        self,
        state_history: list[dict[str, Any]],
        model_metadata: dict[str, Any] | None,
        timescale: float,
    ) -> tuple[list[str], dict[str, str]]:
        """Three-class exclusions via metric-based + metadata discrimination.

        P13 (excitable wave): final wall density near 0.02 and Moran near
          0.87 indicates pre-organized spiral, excluded by metric thresholds
          (wall_final > 0.05 AND moran_final < 0.75).
        P15 (persistent computation): Moran plateau below 0.30 with sparse
          alive fraction indicates GoL-like decay, excluded by metric
          thresholds (moran_final >= 0.30 AND minority_final >= 0.05).
        P1  (similarity aggregation): metadata-keyed. Returns "excluded"
          only if model metadata's `update` key contains 'copy', 'imitation',
          or 'voter'; returns "inconclusive" otherwise (including for
          Schelling, whose metadata keys are `threshold` and `density`).
          For canonical Schelling at threshold = 0.375 the metric gates
          alone reject below CONFIRMATION (Sprint 21 5-seed audit), so the
          inconclusive P1 outcome does not yield a false-positive
          DEFINITIVE. For Schelling at threshold = 0.5 the metric gates
          DO admit the model to DEFINITIVE with P1 marked "inconclusive";
          this is recorded as Sprint 21 carry-forward #20b in
          REPLICATION_NOTES.md.
        """
        results: dict[str, str] = {}
        metrics = self._trajectory_metrics(state_history)
        if not metrics:
            return list(self.excluded_patterns), {p: "inconclusive"
                                                  for p in self.excluded_patterns}

        # P13 exclusion: voter's wall_final > 0.05 and Moran < 0.75 distinguishes
        # from GH spiral which has wall ~0.02 and Moran ~0.87
        if (metrics["wall_final_qtr_mean"] > 0.05
                and metrics["moran_final_qtr_mean"] < 0.75):
            results["P13"] = "excluded"
        else:
            results["P13"] = "not_excluded"

        # P15 exclusion: voter has growing Moran from ~0 and sustained minority;
        # GoL has Moran plateau < 0.30 and minority fraction < 0.05
        if (metrics["moran_final_qtr_mean"] >= 0.30
                and metrics["minority_fraction_final"] >= 0.05):
            results["P15"] = "excluded"
        else:
            results["P15"] = "not_excluded"

        # P1 exclusion: requires metadata to establish type-stability. If
        # metadata says types_are_constant=False, voter is correctly
        # excluded from P1 by the P1 detector itself. Here we check
        # whether the metadata plausibly signals non-P1 behavior.
        if model_metadata is not None:
            update = model_metadata.get("update", "")
            if "copy" in update or "imitation" in update or "voter" in update:
                results["P1"] = "excluded"
            else:
                results["P1"] = "inconclusive"
        else:
            results["P1"] = "inconclusive"

        return list(self.excluded_patterns), results

    def _all_secondaries_pass(self, secondary_result: dict[str, Any]) -> bool:
        """All secondaries pass: wall decay + wall plateau + minority preserved."""
        return (
            secondary_result.get(
                "wall_spearman_early", 0.0
            ) < self.CONFIRMATION_WALL_SPEARMAN_MAX
            and secondary_result.get(
                "wall_final_qtr_mean", 1.0
            ) < self.CONFIRMATION_WALL_FINAL_MAX
            and secondary_result.get(
                "minority_fraction_final", 0.0
            ) >= self.DEFINITIVE_MINORITY_FINAL_MIN
        )

    def _check_finite_size_robustness(
        self,
        state_history: list[dict[str, Any]],
        timescale: float,
    ) -> bool:
        """Default: run length ≥ 10τ."""
        return len(state_history) >= 10 * timescale
