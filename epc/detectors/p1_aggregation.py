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
            # Sprint 14 B.1: substrate-mismatch guard. P1 needs categorical
            # integer labels (1D 'cell_types' or 2D 'grid'/'type_labels_at_pos').
            # Continuous-valued fields (e.g. Gray-Scott 'field') are not valid
            # inputs — graceful reject with warning, matching P11/P13 pattern.
            if ("grid" not in s and "type_labels_at_pos" not in s
                    and "cell_types" not in s):
                warnings.append(
                    "no 'grid', 'type_labels_at_pos', or 'cell_types' — "
                    "P1 needs integer-labeled spatial data"
                )
                return warnings
            n = s.get("n", 0)
            if "grid_dims" in s:
                n = s["grid_dims"][0] * s["grid_dims"][1]
            if n < 50:
                warnings.append(f"N = {n} < 50 (finite-size risk)")
        return warnings

    def _compute_primary(self, state_history, timescale):
        """Compute Moran's I — uses FINAL-STATE I for aggregation.

        Sprint 10 change (was: peak over trajectory; final is more conservative
        and distinguishes transient wavefronts from genuine aggregation).

        RATIONALE: Models with transient wavefronts (SIR epidemic, traveling
        excitation waves) produce high peak Moran's I during the wavefront
        but near-zero final Moran's I. The peak-based primary metric caused
        SIR × P1 to falsely pass screening as a "transient aggregation"
        co-occurrence when in fact the final state is uniform recovered
        cells (I ≈ 0.02). A final-state primary correctly distinguishes:
          - Schelling, Nowak-May: final I ≈ 0.41–0.53 — genuine aggregation.
          - RPS spiral domains: final I ≈ 0.55 — rotating spirals maintain
            persistent clustering at every timestep.
          - SIR epidemic: final I ≈ 0.02 — no aggregation at steady state.

        The sustained Moran (mean over last 20%) is reported alongside as
        a diagnostic; for robust aggregators, sustained ≈ final, and for
        transient wavefronts like SIR, sustained >> final.

        See PROJECT_STATUS.md decision 32 and REPLICATION_NOTES.md Sprint 10
        entry for the empirical characterization across 6 canonical models.

        Moran's I is degenerate (0/0) with only one type, so we still
        require n_unique_types ≥ 2.
        """
        # Sprint 14 B.1: substrate-mismatch guard. Return benign zeros when
        # state_history lacks integer-labeled spatial data; _check_screening
        # will reject cleanly on n_unique_types=0 / observed<=expected=0.
        if not state_history:
            return {
                "morans_i": 0.0, "morans_i_final": 0.0,
                "morans_i_sustained": 0.0, "morans_i_peak": 0.0,
                "expected_i": 0.0, "n_unique_types": 0,
                "screening_rejection_reason": "empty_state_history",
            }
        s0 = state_history[0]
        if ("grid" not in s0 and "type_labels_at_pos" not in s0
                and "cell_types" not in s0):
            return {
                "morans_i": 0.0, "morans_i_final": 0.0,
                "morans_i_sustained": 0.0, "morans_i_peak": 0.0,
                "expected_i": 0.0, "n_unique_types": 0,
                "screening_rejection_reason": "substrate_mismatch",
            }

        result_final = self._morans_i.compute(state_history, timestep=-1)

        # Count unique types from the final state
        final_state = state_history[-1]
        n_unique_types = len(result_final.get("per_type", {}))
        if n_unique_types == 0:
            # Fallback: count from state directly
            if "cell_types" in final_state:
                n_unique_types = len(set(final_state["cell_types"]))
            elif "grid" in final_state:
                n_unique_types = len(np.unique(final_state["grid"]))

        # Sample trajectory to compute peak and sustained (diagnostics).
        n_states = len(state_history)
        step = max(1, n_states // 100)
        sampled_idxs = list(range(0, n_states, step))
        if sampled_idxs[-1] != n_states - 1:
            sampled_idxs.append(n_states - 1)

        sampled_i = []
        for t in sampled_idxs:
            mi = self._morans_i.compute(state_history, timestep=t)
            sampled_i.append(float(mi["morans_i"]))
        sampled_i = np.array(sampled_i)

        peak_i = float(np.max(sampled_i))
        final_i = float(result_final["morans_i"])

        # Sustained: mean over samples from the last 20% of the trajectory.
        last_20_start = max(0, int(n_states * 0.8))
        sustained_samples = [sampled_i[i] for i, t in enumerate(sampled_idxs)
                             if t >= last_20_start]
        if not sustained_samples:
            sustained_samples = [final_i]
        sustained_i = float(np.mean(sustained_samples))

        # Sprint 14.5 D.3: attach a screening-rejection-reason diagnostic
        # so users can see *why* a screening-reject happened on the
        # primary_metric alone (since secondary_metrics is empty on early
        # exit). The three rejection modes are:
        #   - uniform_state: n_unique_types < 2 (SIR after wavefront collapse)
        #   - below_expected: observed <= expected_i (random noise)
        #   - below_magnitude_floor: observed < 0.05 (transient dissipated)
        expected_i = float(result_final["expected_i"])
        if n_unique_types < 2:
            rejection_reason = "uniform_state"
        elif final_i <= expected_i:
            rejection_reason = "below_expected"
        elif final_i < 0.05:
            rejection_reason = "below_magnitude_floor"
        else:
            rejection_reason = "none"  # screening would pass

        return {
            # PRIMARY METRIC (Sprint 10): final-state Moran's I.
            # Kept under the key "morans_i" for compatibility with
            # existing secondaries/null machinery that reads this key.
            "morans_i": final_i,
            "morans_i_final": final_i,
            "morans_i_sustained": sustained_i,
            "morans_i_peak": peak_i,
            "expected_i": expected_i,
            "n_unique_types": n_unique_types,
            # Sprint 14.5 D.3: screening-rejection diagnostic. When
            # detected=False and tier=SCREENING, this names the gate
            # that rejected.
            "screening_rejection_reason": rejection_reason,
        }

    def _check_screening(self, primary_result, timescale):
        """Screening: final-state Moran's I meaningfully above expected.

        Sprint 10 change: added a minimum-magnitude floor of 0.05 on top
        of the original 'I > expected_I' check. The trivial 'I > expected_I'
        check is satisfied by pure random noise on any finite grid (because
        sampling variance gives occasional tiny positive values above
        E[I] = -1/(N-1) ≈ 0). The 0.05 floor excludes this:

          - Schelling final I = 0.41, Nowak-May = 0.53, RPS = 0.55: pass.
          - SIR final I = 0.02: correctly rejected (transient wavefront
            has collapsed, grid is near-uniform).
          - Random grid final I ≈ 0.01–0.03: correctly rejected as noise.

        Moran's I is degenerate (0/0) with only one type, so we still
        require n_unique_types ≥ 2.

        The null-model p-value is a separate, stronger test applied at
        the confirmation tier; screening is a prerequisite gate.
        """
        # Require at least 2 distinct types
        if primary_result.get("n_unique_types", 0) < 2:
            return False

        observed = primary_result["morans_i"]
        expected = primary_result["expected_i"]

        # Must be meaningfully above expected, not just numerically above
        if observed <= expected:
            return False

        # Minimum magnitude floor to reject noise and transient-dissipated states
        SCREENING_FLOOR = 0.05
        if observed < SCREENING_FLOOR:
            return False

        return True

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
        # Sprint 14.5 D.3: CV is undefined when mean is near-zero (e.g.,
        # a uniform post-collapse grid like SIR at t=T_end). Flag this
        # explicitly; keep cv_i=inf for backwards-compat with the
        # numeric `cv > 0.3` confirmation check, but add a boolean so
        # downstream consumers can distinguish "huge CV" from "mean
        # too small to define CV".
        CV_UNDEFINED_MEAN_FLOOR = 0.01
        cv_undefined = bool(sustained_i) and abs(mean_i) < CV_UNDEFINED_MEAN_FLOOR
        if sustained_i and abs(mean_i) >= CV_UNDEFINED_MEAN_FLOOR:
            cv_i = float(np.std(sustained_i) / abs(mean_i))
        else:
            cv_i = float("inf")

        return {
            "segregation_index": seg["segregation_index"],
            "segregation_std": seg["segregation_std"],
            "cluster_count": clust["cluster_count"],
            "mean_cluster_size": clust["mean_cluster_size"],
            "max_cluster_size": clust["max_cluster_size"],
            "sustained_i_mean": mean_i,
            "sustained_i_cv": cv_i,
            "sustained_i_cv_undefined": cv_undefined,
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

        # Sprint 14 B.1: skip P3 FFT exclusion on substrates without grid.
        if not state_history:
            return checked, results
        final = state_history[-1]
        if "grid" not in final and "type_labels_at_pos" not in final:
            return checked, results

        # P3 exclusion: FFT peak check (only meaningful for 2D)
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
