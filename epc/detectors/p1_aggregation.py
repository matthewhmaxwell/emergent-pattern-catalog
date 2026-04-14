"""P1 — Similarity-driven aggregation detector.

Updated to handle both 1D arrays (chimeric sorting) and 2D grids (Schelling etc.).
Implements the full P1 detector card from detector_cards.md.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from epc.base_detector import BaseDetector
from epc.detector_result import DetectionTier, NullType
from epc.metrics.aggregation import (
    ClusterStats,
    MoransI,
    SegregationIndex,
    label_shuffle_null,
)


class P1AggregationDetector(BaseDetector):
    """Detector for P1 — Similarity-driven aggregation.

    Observable scope: state-history only.
    Works for 1D (chimeric arrays) and 2D (grid models).
    """

    def __init__(self, n_permutations: int = 999) -> None:
        super().__init__(
            pattern_id="P1",
            excluded_patterns=["P2", "P3", "P30"],
            allowed_co_occurrences=["P31", "P27"],
            observable_scope="state_history_only",
        )
        self.n_permutations = n_permutations
        self._morans_i = MoransI()
        self._segregation = SegregationIndex()
        self._clusters = ClusterStats()

    def _estimate_timescale(self, state_history, model_metadata):
        if model_metadata and "algorithm" in model_metadata:
            return float(len(state_history))
        return max(1.0, len(state_history) / 10.0)

    def _validate_prerequisites(self, state_history, timescale):
        warnings = super()._validate_prerequisites(state_history, timescale)
        if state_history:
            s = state_history[0]
            n = s.get("n", 0)
            if "grid_dims" in s:
                n = s["grid_dims"][0] * s["grid_dims"][1]
            if n < 50:
                warnings.append(f"N = {n} < 50 (finite-size risk)")
        return warnings

    def _compute_primary(self, state_history, timescale):
        """Compute Moran's I — uses peak over trajectory for transient patterns.

        Chimeric array aggregation peaks during sorting then may decrease.
        We check both final and peak values.
        """
        result_final = self._morans_i.compute(state_history, timestep=-1)

        # Sample trajectory for peak detection (every ~1% of run)
        n_states = len(state_history)
        step = max(1, n_states // 100)
        peak_i = result_final["morans_i"]
        for t in range(0, n_states, step):
            mi = self._morans_i.compute(state_history, timestep=t)
            peak_i = max(peak_i, mi["morans_i"])

        return {
            "morans_i": max(result_final["morans_i"], peak_i),
            "morans_i_final": result_final["morans_i"],
            "morans_i_peak": peak_i,
            "expected_i": result_final["expected_i"],
        }

    def _check_screening(self, primary_result, timescale):
        return primary_result["morans_i"] > primary_result["expected_i"]

    def _compute_secondaries(self, state_history, timescale):
        seg = self._segregation.compute(state_history, timestep=-1)
        clust = self._clusters.compute(state_history, timestep=-1)

        n_states = len(state_history)
        last_20_start = max(0, int(n_states * 0.8))
        sustained_i = []
        step = max(1, (n_states - last_20_start) // 10)
        for t in range(last_20_start, n_states, step):
            mi = self._morans_i.compute(state_history, timestep=t)
            sustained_i.append(mi["morans_i"])

        mean_i = float(np.mean(sustained_i)) if sustained_i else 0.0
        cv_i = (float(np.std(sustained_i) / mean_i)
                if sustained_i and mean_i > 0 else float("inf"))

        return {
            "segregation_index": seg["segregation_index"],
            "segregation_std": seg["segregation_std"],
            "cluster_count": clust["cluster_count"],
            "mean_cluster_size": clust["mean_cluster_size"],
            "max_cluster_size": clust["max_cluster_size"],
            "sustained_i_mean": mean_i,
            "sustained_i_cv": cv_i,
        }

    def _run_null_model(self, state_history, primary_result, timescale):
        final_state = state_history[-1]
        rng = np.random.default_rng(0)
        null_values = label_shuffle_null(final_state, self.n_permutations, rng)

        observed_i = primary_result["morans_i"]
        null_mean = float(np.mean(null_values))
        null_std = float(np.std(null_values))

        p_value = float(np.mean(null_values >= observed_i))
        if p_value == 0:
            p_value = 1.0 / (self.n_permutations + 1)

        return (
            p_value,
            NullType.SHUFFLE,
            {"mean": null_mean, "std": null_std},
        )

    def _check_confirmation(self, primary_result, secondary_result, null_p, timescale):
        if null_p >= 0.01:
            return False
        seg = secondary_result.get("segregation_index", 0)
        if seg < 0.4:
            return False
        cv = secondary_result.get("sustained_i_cv", float("inf"))
        if cv > 0.3:
            return False
        return True

    def _check_definitive(self, primary_result, secondary_result, null_p, null_type,
                          state_history, model_metadata, timescale):
        return null_p < 0.001

    def _check_exclusions(self, state_history, model_metadata, timescale):
        checked = ["P2", "P3", "P30"]
        results = {"P2": "not_checked", "P3": "not_checked", "P30": "not_checked"}

        # P3 exclusion: FFT peak check (only meaningful for 2D)
        final = state_history[-1]
        if "grid_dims" in final:
            try:
                rows, cols = final["grid_dims"]
                if "type_labels_at_pos" in final:
                    grid = np.array(final["type_labels_at_pos"]).reshape(rows, cols).astype(float)
                else:
                    grid = np.array(final["grid"]).astype(float)
                fft2 = np.fft.fft2(grid - grid.mean())
                power = np.abs(fft2) ** 2
                power[0, 0] = 0
                peak_to_mean = power.max() / power.mean() if power.mean() > 0 else 0
                results["P3"] = "not_excluded" if peak_to_mean > 5.0 else "excluded"
            except Exception:
                results["P3"] = "inconclusive"

        return checked, results

    def _all_secondaries_pass(self, secondary_result):
        seg = secondary_result.get("segregation_index", 0)
        cv = secondary_result.get("sustained_i_cv", float("inf"))
        return seg > 0.4 and cv < 0.3
