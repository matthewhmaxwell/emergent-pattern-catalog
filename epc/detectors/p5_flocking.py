"""P5 — Translational alignment (flocking) detector.

Detects global directional alignment achieved through local
velocity-matching interactions. Primary metric: polarization order
parameter φ = |(1/N) Σ v̂_i|.

See docs/detector_cards.md § P5 for full specification.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from epc.base_detector import BaseDetector
from epc.detector_result import NullType
from epc.metrics.collective_motion import (
    AngularMomentum,
    GroupSpeedRatio,
    Polarization,
)


class P5FlockingDetector(BaseDetector):
    """Detector for P5: translational alignment (flocking).

    Primary metric: polarization φ (mean over sustained window).
    Null model: heading-shuffle (preserves speeds, destroys alignment).
    """

    def __init__(self, n_permutations: int = 999) -> None:
        super().__init__(
            pattern_id="P5",
            excluded_patterns=["P6", "P7", "P8"],
            allowed_co_occurrences=["P17", "P19", "P9"],
            observable_scope="state_history_only",
        )
        self.n_permutations = n_permutations
        self._pol_metric = Polarization()
        self._am_metric = AngularMomentum()
        self._gsr_metric = GroupSpeedRatio()

    def _estimate_timescale(
        self,
        state_history: list[dict[str, Any]],
        model_metadata: dict[str, Any] | None,
    ) -> float:
        """T_cross = box_size / speed."""
        if model_metadata:
            box = model_metadata.get("box_size", None)
            speed = model_metadata.get("speed", None)
            if box and speed and speed > 0:
                return box / speed

        if state_history and "box_size" in state_history[0]:
            box = state_history[0]["box_size"]
            v = state_history[0].get("velocities", None)
            if v is not None:
                speed = float(np.linalg.norm(v, axis=1).mean())
                if speed > 0:
                    return box / speed

        return max(1.0, len(state_history) / 10.0)

    def _validate_prerequisites(
        self,
        state_history: list[dict[str, Any]],
        timescale: float,
    ) -> list[str]:
        warnings = super()._validate_prerequisites(state_history, timescale)
        if state_history and "velocities" not in state_history[0]:
            warnings.append("missing 'velocities' key — not an active-matter model")
        return warnings

    def _compute_primary(
        self,
        state_history: list[dict[str, Any]],
        timescale: float,
    ) -> dict[str, float]:
        if not state_history or "velocities" not in state_history[0]:
            return {"polarization_mean": 0.0, "polarization_std": 1.0}

        result = self._pol_metric.compute(state_history)
        return {
            "polarization_mean": result.get("polarization_mean", result.get("polarization", 0.0)),
            "polarization_std": result.get("polarization_std", 0.0),
        }

    def _check_screening(
        self,
        primary_result: dict[str, float],
        timescale: float,
    ) -> bool:
        """φ_mean > 0.5."""
        return primary_result.get("polarization_mean", 0) > 0.5

    def _compute_secondaries(
        self,
        state_history: list[dict[str, Any]],
        timescale: float,
    ) -> dict[str, Any]:
        if not state_history or "velocities" not in state_history[0]:
            return {}

        gsr = self._gsr_metric.compute(state_history)
        am = self._am_metric.compute(state_history)

        # Directional persistence: heading autocorrelation at lag ~T_cross
        lag = max(1, min(int(timescale), len(state_history) - 1))
        headings_series = []
        for s in state_history:
            if "headings" in s:
                headings_series.append(s["headings"])

        persistence = 0.0
        if len(headings_series) > lag:
            # Circular correlation at lag
            cos_diffs = []
            for t in range(len(headings_series) - lag):
                diff = headings_series[t + lag] - headings_series[t]
                cos_diffs.append(float(np.mean(np.cos(diff))))
            persistence = float(np.mean(cos_diffs))

        return {
            "group_speed_ratio": gsr.get("group_speed_ratio", gsr.get("group_speed_ratio_mean", 0.0)),
            "angular_momentum_abs": am.get("angular_momentum_abs_mean", abs(am.get("angular_momentum", 0.0))),
            "directional_persistence": persistence,
        }

    def _run_null_model(
        self,
        state_history: list[dict[str, Any]],
        primary_result: dict[str, float],
        timescale: float,
    ) -> tuple[float, NullType, dict[str, float]]:
        """Heading-shuffle null: permute headings across agents.

        Expected φ_null ≈ 1/√N for random headings.
        Sample 10 evenly-spaced frames, run 999 permutations each.
        """
        rng = np.random.default_rng(0)
        observed_phi = primary_result.get("polarization_mean", 0)

        if not state_history or "velocities" not in state_history[0]:
            return 1.0, NullType.SHUFFLE, {"mean": 0.0, "std": 0.0}

        # Sample frames for null
        n_frames = len(state_history)
        n_sample = min(10, n_frames)
        sample_idx = np.linspace(0, n_frames - 1, n_sample, dtype=int)

        null_phis: list[float] = []
        for _ in range(self.n_permutations):
            frame_phis = []
            for idx in sample_idx:
                v = state_history[int(idx)]["velocities"].copy()
                # Shuffle heading assignment across agents
                perm = rng.permutation(len(v))
                v_shuffled = v[perm]
                speeds = np.linalg.norm(v_shuffled, axis=1, keepdims=True)
                speeds = np.where(speeds == 0, 1.0, speeds)
                v_hat = v_shuffled / speeds
                phi = float(np.linalg.norm(v_hat.mean(axis=0)))
                frame_phis.append(phi)
            null_phis.append(float(np.mean(frame_phis)))

        null_arr = np.array(null_phis)
        p_value = float(np.mean(null_arr >= observed_phi))
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
        """φ_mean > 0.7 AND R > 0.5 AND shuffle p < 0.01."""
        phi = primary_result.get("polarization_mean", 0)
        R = secondary_result.get("group_speed_ratio", 0)
        return phi > 0.7 and R > 0.5 and null_p < 0.01

    def _check_exclusions(
        self,
        state_history: list[dict[str, Any]],
        model_metadata: dict[str, Any] | None,
        timescale: float,
    ) -> tuple[list[str], dict[str, str]]:
        results: dict[str, str] = {}

        if not state_history or "velocities" not in state_history[0]:
            return (["P6", "P7", "P8"], {"P6": "inconclusive", "P7": "inconclusive", "P8": "not_checked"})

        # P6: milling — |L| > 0.5 → not_excluded (milling present, exclude P5)
        am = self._am_metric.compute(state_history)
        L_abs = am.get("angular_momentum_abs_mean", abs(am.get("angular_momentum", 0)))
        results["P6"] = "not_excluded" if L_abs > 0.5 else "excluded"

        # P7: lane formation — bimodal antiparallel heading distribution
        # Simple test: check if heading histogram has two peaks ~π apart
        all_headings = []
        n_sample = min(10, len(state_history))
        sample_idx = np.linspace(0, len(state_history) - 1, n_sample, dtype=int)
        for idx in sample_idx:
            if "headings" in state_history[int(idx)]:
                all_headings.extend(state_history[int(idx)]["headings"].tolist())

        if all_headings:
            h = np.array(all_headings)
            # Check bimodality via circular variance of doubled angles
            # If bimodal antiparallel, 2*theta should be unimodal
            doubled = 2.0 * h
            r2 = abs(np.mean(np.exp(1j * doubled)))
            # High r2 means antiparallel bimodal
            results["P7"] = "not_excluded" if r2 > 0.7 else "excluded"
        else:
            results["P7"] = "inconclusive"

        # P8: jamming — requires spacetime FFT, stub for now
        results["P8"] = "not_checked"

        return (["P6", "P7", "P8"], results)
