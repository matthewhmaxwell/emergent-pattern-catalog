"""P22 — Information cascade / contagion detector.

Detects spatially-structured epidemic-like spreading where a signal (infection,
information, behavior) propagates through local contact, producing a wave-like
cascade that is statistically distinguishable from random transmission.

Observable scope: state-history only.
Works for 2D grid models with discrete states including at least one
"susceptible" and one "infected/active" class.

Detection tiers (from detector design):
  Screening:    Cascade reaches ≥ 5% of susceptible population AND
                spatial clustering of affected cells significantly exceeds
                random expectation (Moran's I test).
  Confirmation: Epidemic curve is unimodal (single peak) AND estimated R0 > 1
                AND cascade p-value < 0.01 against random-timing null.
  Definitive:   Confirmation + cascade reaches ≥ 30% AND R0 estimate
                consistent across estimation methods.

Null models:
  Random-timing null: preserve total number of new infections per step but
  assign them to random susceptible cells (breaks spatial correlation).

Nearest-neighbor exclusions:
  P13: SIR waves are single-pass (transient); P13 excitable waves are
       persistent (re-entrant). Temporal persistence discriminates.
  P1:  Cascade produces spatial clustering of recovered cells, but this is
       a consequence of spreading, not preference-based aggregation. The
       cascade detector looks at the PROCESS (spreading dynamics) while
       P1 looks at the OUTCOME (final spatial pattern).
"""

from __future__ import annotations

from typing import Any

import numpy as np

from epc.base_detector import BaseDetector
from epc.detector_result import DetectionTier, NullType


class P22CascadeDetector(BaseDetector):
    """Detector for P22 — Information cascade / contagion.

    Detects spatially-structured epidemic spreading through local contact.
    The primary signal is spatial correlation in cascade timing — cells near
    an infected cell become infected sooner than distant cells, producing
    a measurable wavefront.

    The key discrimination from P13 (excitable waves) is that cascades are
    transient and single-pass, while excitable waves are persistent and
    re-entrant. The key discrimination from P1 (aggregation) is that the
    cascade detector measures the SPREADING PROCESS, not just the final
    spatial pattern.
    """

    def __init__(self, n_permutations: int = 199, seed: int = 42) -> None:
        super().__init__(
            pattern_id="P22",
            excluded_patterns=["P13"],
            allowed_co_occurrences=["P1"],
            observable_scope="state_history_only",
        )
        self.n_permutations = n_permutations
        self._seed = seed

    def _estimate_timescale(
        self,
        state_history: list[dict[str, Any]],
        model_metadata: dict[str, Any] | None,
    ) -> float:
        """T_spread = max(rows, cols): time for cascade to cross the grid."""
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
        if state_history:
            s = state_history[0]
            if "grid" not in s:
                warnings.append("no 'grid' key in state history — P22 needs 2D grid")
            if "grid_dims" not in s:
                warnings.append("no 'grid_dims' — cannot determine grid size")
        if len(state_history) < 10:
            warnings.append("state history too short for epidemic analysis")
        return warnings

    def _compute_primary(
        self,
        state_history: list[dict[str, Any]],
        timescale: float,
    ) -> dict[str, float]:
        """Primary metrics: cascade reach, R0 estimate, epidemic curve shape.

        Cascade reach: final fraction of cells that transitioned out of the
        susceptible state (i.e., were infected at some point).

        R0 estimate: from early-epidemic exponential growth rate.

        Moran's I on infection time map: the key spatial metric. The binary
        recovered pattern becomes nearly uniform when the epidemic infects
        most cells, destroying the clustering signal. But the infection TIME
        map (when each cell first became infected) retains strong spatial
        structure — cells near the seed are infected early, distant cells
        late. Moran's I on this continuous variable captures the wavefront
        spreading process regardless of final epidemic size.

        Epidemic curve: unimodality check via peak count in infected fraction.
        """
        if not state_history or "grid" not in state_history[0]:
            return {
                "cascade_reach": 0.0,
                "r0_estimate": 0.0,
                "peak_infected_fraction": 0.0,
                "is_unimodal": 0.0,
                "epidemic_died_out": 1.0,
                "moran_i_time": 0.0,
            }

        rows, cols = state_history[0]["grid_dims"]
        total = rows * cols

        # Extract time series
        i_series = []
        r_series = []
        new_infected_series = []

        for t, state in enumerate(state_history):
            grid = state["grid"]
            i_count = int((grid == 1).sum())
            r_count = int((grid == 2).sum())
            i_series.append(i_count)
            r_series.append(r_count)
            new_inf = state.get("newly_infected", 0)
            new_infected_series.append(new_inf)

        i_series = np.array(i_series, dtype=float)
        r_series = np.array(r_series, dtype=float)
        new_infected_series = np.array(new_infected_series, dtype=float)

        # Cascade reach: max R fraction achieved
        cascade_reach = float(r_series[-1] / total) if total > 0 else 0.0
        # Also add the current infected (they will eventually recover)
        cascade_reach_total = float((r_series[-1] + i_series[-1]) / total)

        # Peak infected fraction
        peak_infected_fraction = float(i_series.max() / total) if total > 0 else 0.0

        # Epidemic died out? (no more infected at end)
        epidemic_died_out = 1.0 if i_series[-1] == 0 else 0.0

        # Unimodality check: count peaks in smoothed I(t)
        is_unimodal = self._check_unimodal(i_series)

        # R0 estimation from early growth
        r0_estimate = self._estimate_r0(new_infected_series, i_series)

        # Spatial clustering: Moran's I on infection TIME map (continuous)
        # This captures the wavefront structure even when the final binary
        # pattern is nearly uniform (high cascade reach).
        infection_time = self._build_infection_time_map(state_history)
        moran_i_time = self._moran_i_infection_time(infection_time, rows, cols)

        return {
            "cascade_reach": cascade_reach,
            "cascade_reach_total": cascade_reach_total,
            "r0_estimate": r0_estimate,
            "peak_infected_fraction": peak_infected_fraction,
            "is_unimodal": 1.0 if is_unimodal else 0.0,
            "epidemic_died_out": epidemic_died_out,
            "moran_i_time": moran_i_time,
        }

    def _check_screening(
        self,
        primary_result: dict[str, float],
        timescale: float,
    ) -> bool:
        """Screening: cascade reaches ≥5% AND shows spatial structure in timing."""
        reach = primary_result.get("cascade_reach_total", 0.0)
        if reach < 0.05:
            return False

        # Must show spatial structure in infection timing
        moran = primary_result.get("moran_i_time", 0.0)
        if moran < 0.1:
            return False

        return True

    def _compute_secondaries(
        self,
        state_history: list[dict[str, Any]],
        timescale: float,
    ) -> dict[str, Any]:
        """Secondary metrics: wavefront velocity, cascade duration."""
        if not state_history or "grid" not in state_history[0]:
            return {"wavefront_velocity": 0.0, "cascade_duration": 0}

        rows, cols = state_history[0]["grid_dims"]

        # Measure wavefront velocity from infection time map
        infection_time = self._build_infection_time_map(state_history)
        wavefront_velocity = self._measure_wavefront_velocity(
            infection_time, rows, cols
        )

        # Cascade duration: time from first infection to last recovery
        i_series = [int((s["grid"] == 1).sum()) for s in state_history]
        active_steps = [t for t, i in enumerate(i_series) if i > 0]
        cascade_duration = (
            (active_steps[-1] - active_steps[0]) if len(active_steps) >= 2 else 0
        )

        return {
            "wavefront_velocity": wavefront_velocity,
            "cascade_duration": cascade_duration,
        }

    def _run_null_model(
        self,
        state_history: list[dict[str, Any]],
        primary_result: dict[str, float],
        timescale: float,
    ) -> tuple[float, NullType, dict[str, float]]:
        """Null model: random-timing shuffle on infection time map.

        Preserves the NUMBER of new infections per timestep but assigns
        them to random cells instead of spatially correlated ones. This
        destroys the spatial wavefront while preserving the epidemic curve
        envelope. Computes Moran's I on the resulting infection-time map.
        """
        if not state_history or "grid" not in state_history[0]:
            return 1.0, NullType.SHUFFLE, {"mean": 0.0, "std": 0.0}

        observed_moran = primary_result.get("moran_i_time", 0.0)
        rows, cols = state_history[0]["grid_dims"]
        total = rows * cols

        # Extract new infection counts per step
        new_infections_per_step = []
        for state in state_history:
            new_infections_per_step.append(state.get("newly_infected", 0))

        rng = np.random.default_rng(self._seed)
        null_morans = []

        for _ in range(self.n_permutations):
            # Build null infection-time map: same counts per step,
            # random spatial assignment
            null_time = np.full((rows, cols), -1, dtype=int)
            available = list(range(total))
            rng.shuffle(available)
            cursor = 0

            for t, n_new in enumerate(new_infections_per_step):
                if n_new == 0:
                    continue
                n_to_assign = min(n_new, len(available) - cursor)
                if n_to_assign <= 0:
                    break
                for idx in available[cursor : cursor + n_to_assign]:
                    r, c = divmod(idx, cols)
                    null_time[r, c] = t
                cursor += n_to_assign

            null_moran = self._moran_i_infection_time(null_time, rows, cols)
            null_morans.append(null_moran)

        null_morans = np.array(null_morans)
        null_mean = float(null_morans.mean())
        null_std = float(null_morans.std())

        # P-value: fraction of null Moran's I ≥ observed
        p_value = float(np.mean(null_morans >= observed_moran))
        if p_value == 0:
            p_value = 1.0 / (self.n_permutations + 1)

        return p_value, NullType.SHUFFLE, {"mean": null_mean, "std": null_std}

    def _compute_effect_size(
        self,
        primary_result: dict[str, float],
        null_dist_stats: dict[str, float],
    ) -> dict[str, float]:
        """Effect size: how much more spatially structured is cascade timing than random."""
        if not null_dist_stats or null_dist_stats.get("std", 0) == 0:
            return {}

        observed = primary_result.get("moran_i_time", 0.0)
        null_mean = null_dist_stats.get("mean", 0.0)
        null_std = null_dist_stats.get("std", 1.0)

        cohens_d = (observed - null_mean) / null_std if null_std > 0 else 0.0

        return {
            "cohens_d": cohens_d,
            "raw_value": observed,
            "null_mean": null_mean,
            "null_std": null_std,
            "note": "positive_d_means_more_spatially_structured_than_random",
        }

    def _check_confirmation(
        self,
        primary_result: dict[str, float],
        secondary_result: dict[str, Any],
        null_p: float,
        timescale: float,
    ) -> bool:
        """Confirmation: unimodal curve AND R0 > 1 AND p < 0.01."""
        if null_p >= 0.01:
            return False
        if primary_result.get("is_unimodal", 0.0) < 0.5:
            return False
        if primary_result.get("r0_estimate", 0.0) <= 1.0:
            return False
        if primary_result.get("cascade_reach_total", 0.0) < 0.05:
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
        """Definitive: confirmation + cascade ≥ 30% + epidemic died out + strong signal."""
        if primary_result.get("cascade_reach_total", 0.0) < 0.30:
            return False
        if null_p >= 0.01:
            return False
        if primary_result.get("r0_estimate", 0.0) <= 1.0:
            return False
        # Definitive requires the epidemic to have completed (died out)
        if primary_result.get("epidemic_died_out", 0.0) < 0.5:
            return False
        return True

    def _check_exclusions(
        self,
        state_history: list[dict[str, Any]],
        model_metadata: dict[str, Any] | None,
        timescale: float,
    ) -> tuple[list[str], dict[str, str]]:
        """Check P13 exclusion (persistent waves vs. transient cascade).

        A cascade (P22) is single-pass: activity spreads then dies out.
        An excitable wave (P13) is re-entrant: activity persists indefinitely.

        The discrimination logic has three signals, any of which is
        sufficient to exclude P13:

        1. FINAL STATE: if the last recorded step has zero activity, the
           epidemic has definitively died out. Persistent waves cannot
           produce this outcome — they never stop.

        2. TRAILING-WINDOW DECAY: in a window sized to ~2×T_spread at the
           end of the series, activity must be < 10% of peak. This catches
           long simulations where activity hasn't quite hit zero yet but
           has clearly decayed.

        3. FINAL-OVER-PEAK RATIO: end-of-series activity < 1% of peak
           indicates an essentially finished epidemic.

        The previous implementation used a fixed "last quarter" window,
        which fails for late-peaking epidemics where the peak is inside
        that window (window max is then large even when the epidemic has
        finished by the end of the series).
        """
        checked = ["P13"]
        results: dict[str, str] = {}

        if not state_history:
            results["P13"] = "not_checked"
            return checked, results

        i_series = [int((s["grid"] == 1).sum()) for s in state_history]
        peak = max(i_series)

        if peak == 0:
            # No infection ever happened; P13 can't apply
            results["P13"] = "excluded"
            return checked, results

        final = i_series[-1]

        # Signal 1: final state is zero → epidemic completed, clearly excluded
        if final == 0:
            results["P13"] = "excluded"
            return checked, results

        # Signal 3: final-over-peak ratio is tiny → practically done
        if final < 0.01 * peak:
            results["P13"] = "excluded"
            return checked, results

        # Signal 2: trailing window of ~2×T_spread (or 10% of series, whichever
        # is larger) has decayed below 10% of peak
        trailing_len = max(int(2 * timescale), len(i_series) // 10, 5)
        trailing_len = min(trailing_len, len(i_series))
        trailing = i_series[-trailing_len:]
        if max(trailing) < 0.1 * peak:
            results["P13"] = "excluded"
            return checked, results

        # Otherwise activity is still substantial at end — truly inconclusive
        results["P13"] = "inconclusive"
        return checked, results

    def _all_secondaries_pass(self, secondary_result: dict[str, Any]) -> bool:
        """Secondary checks: wavefront velocity should be positive."""
        velocity = secondary_result.get("wavefront_velocity", 0.0)
        return velocity > 0

    # --- Helper methods ---

    def _estimate_r0(
        self,
        new_infected_series: np.ndarray,
        i_series: np.ndarray,
    ) -> float:
        """Estimate R0 from the ratio of secondary to primary infections.

        Uses the growth phase (before peak) where I(t) is increasing.
        R0 ≈ mean(new_infections[t] / current_infected[t-1]) during growth,
        which estimates the per-infected-per-step reproduction.
        """
        if len(new_infected_series) < 5:
            return 0.0

        # Find peak of epidemic
        peak_t = int(np.argmax(i_series))
        if peak_t < 3:
            return 0.0

        # Use growth phase: steps 1 to peak
        ratios = []
        for t in range(2, min(peak_t, len(new_infected_series))):
            if i_series[t - 1] > 0 and new_infected_series[t] > 0:
                # Generation ratio: new infections per current infected
                ratio = new_infected_series[t] / i_series[t - 1]
                ratios.append(ratio)

        if not ratios:
            return 0.0

        # R0 = mean reproduction ratio + 1 (since we're counting new per existing)
        # Actually, for SIR CA: each infected generates new_infections/n_infected
        # new cases, and persists for ~1/q steps, so R0 ≈ (new/existing) / q
        # But simpler: just use the average generation ratio as R_effective
        r_eff = float(np.median(ratios))

        # R_eff > 0 means epidemic is growing; convert to R0-like quantity
        # R0 ≈ 1 + r_eff for discrete-time SIR with these dynamics
        return 1.0 + r_eff

    def _check_unimodal(self, i_series: np.ndarray) -> bool:
        """Check if the infected time series is unimodal (single peak).

        Uses a simple algorithm: smooth with a window proportional to series
        length, then count significant local maxima. A single significant
        peak indicates a classic epidemic curve.
        """
        if len(i_series) < 5:
            return False

        max_val = i_series.max()
        if max_val == 0:
            return False

        # Smooth with moving average — wider window for robustness
        window = max(5, len(i_series) // 10)
        if window % 2 == 0:
            window += 1
        window = min(window, len(i_series))
        kernel = np.ones(window) / window
        smoothed = np.convolve(i_series, kernel, mode="same")

        # Find peaks: local maxima above 20% of global max
        # Use interior points only (skip boundary artifacts)
        margin = window // 2
        peaks = []
        smooth_max = smoothed[margin:-margin].max() if len(smoothed) > 2 * margin else smoothed.max()

        for t in range(max(1, margin), min(len(smoothed) - 1, len(smoothed) - margin)):
            if smoothed[t] > smoothed[t - 1] and smoothed[t] > smoothed[t + 1]:
                if smoothed[t] >= 0.2 * smooth_max:
                    peaks.append(t)

        return len(peaks) == 1

    def _moran_i_infection_time(
        self,
        infection_time: np.ndarray,
        rows: int,
        cols: int,
    ) -> float:
        """Compute Moran's I on the infection time map (continuous variable).

        Only includes cells that were actually infected (time >= 0).
        Uses rook contiguity (4-connected) with periodic boundaries.

        I = (N/W) * (Σ_ij w_ij (x_i - x̄)(x_j - x̄)) / (Σ_i (x_i - x̄)²)
        """
        valid = infection_time >= 0
        n_valid = int(valid.sum())
        if n_valid < 10:
            return 0.0

        x = infection_time.astype(float)
        x_bar = float(x[valid].mean())

        # Zero out non-infected cells so they don't contribute
        z = np.where(valid, x - x_bar, 0.0)

        denominator = float((z ** 2).sum())
        if denominator == 0:
            return 0.0

        # Sum of cross-products for rook neighbors (4-connected, periodic)
        numerator = 0.0
        for dr, dc in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
            shifted = np.roll(np.roll(z, -dr, axis=0), -dc, axis=1)
            # Only count pairs where both cells are valid
            valid_shifted = np.roll(np.roll(valid, -dr, axis=0), -dc, axis=1)
            both_valid = valid & valid_shifted
            numerator += float((z * shifted * both_valid).sum())

        # W = total number of valid neighbor pairs
        w = 0
        for dr, dc in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
            valid_shifted = np.roll(np.roll(valid, -dr, axis=0), -dc, axis=1)
            w += int((valid & valid_shifted).sum())

        if w == 0:
            return 0.0

        moran_i = (n_valid / w) * (numerator / denominator)
        return float(moran_i)

    def _build_infection_time_map(
        self,
        state_history: list[dict[str, Any]],
    ) -> np.ndarray:
        """Build a map of when each cell first became infected.

        Returns array of shape (rows, cols) with the timestep of first
        infection, or -1 for cells never infected.
        """
        rows, cols = state_history[0]["grid_dims"]
        infection_time = np.full((rows, cols), -1, dtype=int)

        prev_grid = state_history[0]["grid"]
        for t in range(1, len(state_history)):
            grid = state_history[t]["grid"]
            # Cells that just became infected: were 0, now 1
            newly_infected = (prev_grid == 0) & (grid == 1)
            infection_time[newly_infected & (infection_time == -1)] = t
            prev_grid = grid

        # Also mark initially infected cells
        if state_history:
            initially_infected = state_history[0]["grid"] == 1
            infection_time[initially_infected & (infection_time == -1)] = 0

        return infection_time

    def _measure_wavefront_velocity(
        self,
        infection_time: np.ndarray,
        rows: int,
        cols: int,
    ) -> float:
        """Estimate wavefront velocity from the infection time map.

        Compute average distance from seed vs. infection time. Wavefront
        velocity = slope of distance vs. time regression.
        """
        # Find seed location (first infected cell)
        infected_cells = np.argwhere(infection_time >= 0)
        if len(infected_cells) < 5:
            return 0.0

        seed_cells = np.argwhere(infection_time == 0)
        if len(seed_cells) == 0:
            return 0.0

        # Use centroid of initial infection as seed
        seed_r = float(seed_cells[:, 0].mean())
        seed_c = float(seed_cells[:, 1].mean())

        # Compute distances and times for all infected cells
        distances = []
        times = []
        for r, c in infected_cells:
            t = infection_time[r, c]
            if t <= 0:
                continue
            # Periodic distance
            dr = min(abs(r - seed_r), rows - abs(r - seed_r))
            dc = min(abs(c - seed_c), cols - abs(c - seed_c))
            dist = np.sqrt(dr**2 + dc**2)
            distances.append(dist)
            times.append(t)

        if len(distances) < 5:
            return 0.0

        distances = np.array(distances)
        times = np.array(times)

        # Linear regression: distance = velocity * time + intercept
        # Use least squares
        A = np.column_stack([times, np.ones_like(times)])
        result = np.linalg.lstsq(A, distances, rcond=None)
        velocity = float(result[0][0])

        return max(0.0, velocity)
