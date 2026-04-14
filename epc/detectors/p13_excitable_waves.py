"""P13 — Excitable spiral and target wave detector.

Detects sustained excitable waves with constant propagation speed,
spiral tips, and/or persistent wave sources in cellular automata
with rest/excited/refractory state dynamics.

See docs/detector_cards.md § P13 for full specification.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from epc.base_detector import BaseDetector
from epc.detector_result import NullType
from epc.metrics.excitable_waves import (
    SpiralTipDetector,
    WavefrontSpeed,
    WaveSourceCount,
)


class P13ExcitableWaveDetector(BaseDetector):
    """Detector for P13: excitable spiral and target waves.

    Primary metric: wavefront speed persistence (CV < 0.2 for screening).
    Null models: frame-shuffle (always available) or mechanistic
    no-refractory / high-threshold (if model metadata available).
    """

    def __init__(self, n_permutations: int = 199) -> None:
        super().__init__(
            pattern_id="P13",
            excluded_patterns=["P15", "P12"],
            allowed_co_occurrences=["P9"],
            observable_scope="state_history_only",
        )
        self.n_permutations = n_permutations
        self._speed_metric = WavefrontSpeed()
        self._tip_metric = SpiralTipDetector()
        self._source_metric = WaveSourceCount()

    def _estimate_timescale(
        self,
        state_history: list[dict[str, Any]],
        model_metadata: dict[str, Any] | None,
    ) -> float:
        """T_prop = grid_size / wavefront_speed."""
        if model_metadata and "grid_size" in model_metadata:
            return float(model_metadata["grid_size"])

        if state_history and "grid_dims" in state_history[0]:
            return float(state_history[0]["grid_dims"][0])

        return max(1.0, len(state_history) / 10.0)

    def _validate_prerequisites(
        self,
        state_history: list[dict[str, Any]],
        timescale: float,
    ) -> list[str]:
        warnings = super()._validate_prerequisites(state_history, timescale)
        if not state_history:
            warnings.append("empty state history")
            return warnings
        if "grid" not in state_history[0]:
            warnings.append("missing 'grid' key — not a grid-based model")
        return warnings

    def _compute_primary(
        self,
        state_history: list[dict[str, Any]],
        timescale: float,
    ) -> dict[str, float]:
        if not state_history or "grid" not in state_history[0]:
            return {
                "wavefront_speed_mean": 0.0,
                "wavefront_speed_cv": float("inf"),
                "activity_persistence": 0.0,
                "activity_cv": float("inf"),
                "n_wavefront_frames": 0,
            }

        speed_result = self._speed_metric.compute(state_history)

        # Activity persistence: fraction of frames with excited cells
        excited_counts = np.array([
            float((s["grid"] == 1).sum()) for s in state_history
        ])
        active_frames = float((excited_counts > 0).sum())
        persistence = active_frames / len(state_history)

        # Activity consistency: CV of excited cell count (over active frames)
        active_vals = excited_counts[excited_counts > 0]
        if len(active_vals) > 1:
            activity_cv = float(active_vals.std() / active_vals.mean())
        else:
            activity_cv = float("inf")

        return {
            "wavefront_speed_mean": speed_result["wavefront_speed_mean"],
            "wavefront_speed_cv": speed_result["wavefront_speed_cv"],
            "activity_persistence": persistence,
            "activity_cv": activity_cv,
            "n_wavefront_frames": float(speed_result["n_wavefront_frames"]),
        }

    def _check_screening(
        self,
        primary_result: dict[str, float],
        timescale: float,
    ) -> bool:
        """Screening: persistent activity ≥ 80% of frames AND activity CV < 0.5.

        Spiral waves maintain consistent excited-cell counts across frames.
        Transient activity dies out (low persistence) or fluctuates wildly.
        """
        persistence = primary_result.get("activity_persistence", 0)
        activity_cv = primary_result.get("activity_cv", float("inf"))
        return persistence >= 0.8 and activity_cv < 0.5

    def _compute_secondaries(
        self,
        state_history: list[dict[str, Any]],
        timescale: float,
    ) -> dict[str, Any]:
        if not state_history or "grid" not in state_history[0]:
            return {}

        n_states = int(state_history[0]["grid"].max()) + 1
        tip_result = self._tip_metric.compute(state_history, n_states=n_states)
        source_result = self._source_metric.compute(state_history)

        # Refractory consistency: std(n_refractory) / mean(n_refractory) over last 50%
        half = max(1, len(state_history) // 2)
        refrac_counts = [s.get("n_refractory", 0) for s in state_history[half:]]
        if refrac_counts and np.mean(refrac_counts) > 0:
            refrac_cv = float(np.std(refrac_counts) / np.mean(refrac_counts))
        else:
            refrac_cv = float("inf")

        return {
            "n_spiral_tips_max": tip_result["n_spiral_tips_max"],
            "max_rotations": tip_result["max_rotations"],
            "tip_persistence_frames": tip_result["tip_persistence_frames"],
            "topological_charge_net": tip_result["topological_charge_net"],
            "wave_source_count": source_result["wave_source_count"],
            "persistent_sources": source_result["persistent_sources"],
            "refractory_cv": refrac_cv,
        }

    def _run_null_model(
        self,
        state_history: list[dict[str, Any]],
        primary_result: dict[str, float],
        timescale: float,
    ) -> tuple[float, NullType, dict[str, float]]:
        """Frame-shuffle null: permute temporal order, recompute speed CV.

        If frames are shuffled, wavefront tracking breaks → CV → ∞.
        """
        rng = np.random.default_rng(0)
        observed_cv = primary_result.get("wavefront_speed_cv", float("inf"))

        if not state_history or "grid" not in state_history[0]:
            return 1.0, NullType.SHUFFLE, {"mean": float("inf"), "std": 0.0}

        null_cvs: list[float] = []
        indices = np.arange(len(state_history))

        for _ in range(self.n_permutations):
            perm = rng.permutation(indices)
            shuffled = [state_history[int(i)] for i in perm]
            result = self._speed_metric.compute(shuffled)
            null_cvs.append(result["wavefront_speed_cv"])

        null_arr = np.array(null_cvs)
        # P-value: fraction of null CVs <= observed CV (lower is better)
        p_value = float(np.mean(null_arr <= observed_cv))
        if p_value == 0:
            p_value = 1.0 / (self.n_permutations + 1)

        return (
            p_value,
            NullType.SHUFFLE,
            {"mean": float(null_arr.mean()), "std": float(null_arr.std())},
        )

    def _check_confirmation(
        self,
        primary_result: dict[str, float],
        secondary_result: dict[str, Any],
        null_p: float,
        timescale: float,
    ) -> bool:
        """Activity CV < 0.3 AND (spiral tips OR persistent sources)."""
        activity_cv = primary_result.get("activity_cv", float("inf"))
        if activity_cv >= 0.3:
            return False
        if null_p >= 0.01:
            return False

        tips = secondary_result.get("n_spiral_tips_max", 0)
        sources = secondary_result.get("persistent_sources", 0)
        rotations = secondary_result.get("max_rotations", 0)
        return tips > 0 or sources > 0 or rotations > 0

    def _check_exclusions(
        self,
        state_history: list[dict[str, Any]],
        model_metadata: dict[str, Any] | None,
        timescale: float,
    ) -> tuple[list[str], dict[str, str]]:
        results: dict[str, str] = {}

        # P15: Transfer entropy test — stub for Sprint 2
        results["P15"] = "not_checked"

        # P12: Excitable states (P13) vs nontransitive species (P12)
        if model_metadata:
            model_class = model_metadata.get("model_class", "")
            if model_class == "excitable_medium":
                results["P12"] = "excluded"
            elif "nontransitive" in model_class or "rps" in model_class:
                results["P12"] = "not_excluded"
            else:
                results["P12"] = "inconclusive"
        else:
            results["P12"] = "inconclusive"

        return (["P15", "P12"], results)
