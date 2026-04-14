"""P6 — Milling / vortex formation detector.

Detects coherent rotational motion where agents orbit a common center,
producing high angular momentum and low group speed ratio.

See docs/detector_cards.md § P6 for full specification.
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


class P6MillingDetector(BaseDetector):
    """Detector for P6: milling / vortex formation.

    Primary metric: angular momentum |L|.
    Null model: heading-shuffle (|L|_null ≈ 0).
    """

    def __init__(self, n_permutations: int = 999) -> None:
        super().__init__(
            pattern_id="P6",
            excluded_patterns=["P5", "P12"],
            allowed_co_occurrences=["P9"],
            observable_scope="state_history_only",
        )
        self.n_permutations = n_permutations
        self._am_metric = AngularMomentum()
        self._gsr_metric = GroupSpeedRatio()
        self._pol_metric = Polarization()

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
        if state_history and "positions" not in state_history[0]:
            warnings.append("missing 'positions' key")
        if state_history and "velocities" not in state_history[0]:
            warnings.append("missing 'velocities' key")
        return warnings

    def _compute_primary(
        self,
        state_history: list[dict[str, Any]],
        timescale: float,
    ) -> dict[str, float]:
        if not state_history or "velocities" not in state_history[0]:
            return {"angular_momentum_abs_mean": 0.0, "angular_momentum_abs_std": 1.0}

        result = self._am_metric.compute(state_history)
        return {
            "angular_momentum_abs_mean": result.get("angular_momentum_abs_mean", abs(result.get("angular_momentum", 0.0))),
            "angular_momentum_abs_std": result.get("angular_momentum_abs_std", 0.0),
        }

    def _check_screening(
        self,
        primary_result: dict[str, float],
        timescale: float,
    ) -> bool:
        """|L| > 0.3."""
        return primary_result.get("angular_momentum_abs_mean", 0) > 0.3

    def _compute_secondaries(
        self,
        state_history: list[dict[str, Any]],
        timescale: float,
    ) -> dict[str, Any]:
        if not state_history or "velocities" not in state_history[0]:
            return {}

        gsr = self._gsr_metric.compute(state_history)
        pol = self._pol_metric.compute(state_history)

        # Ring density: radial density profile from COM
        # Check hollowness = ρ(0) / ρ(r_peak)
        hollowness = self._compute_hollowness(state_history)

        # Rotation coherence: std of per-agent angular velocities
        rot_coherence = self._compute_rotation_coherence(state_history)

        return {
            "group_speed_ratio": gsr.get("group_speed_ratio", gsr.get("group_speed_ratio_mean", 0.0)),
            "polarization": pol.get("polarization", pol.get("polarization_mean", 0.0)),
            "hollowness": hollowness,
            "rotation_coherence": rot_coherence,
        }

    @staticmethod
    def _compute_hollowness(state_history: list[dict[str, Any]]) -> float:
        """Compute ring hollowness from radial density profile."""
        from epc.metrics.collective_motion import _periodic_com

        densities_center = []
        densities_peak = []

        n_sample = min(10, len(state_history))
        sample_idx = np.linspace(0, len(state_history) - 1, n_sample, dtype=int)

        for idx in sample_idx:
            state = state_history[int(idx)]
            pos = state["positions"]
            box = state.get("box_size", None)
            if box:
                com = _periodic_com(pos, box)
                dr = pos - com
                dr = dr - box * np.round(dr / box)
            else:
                com = pos.mean(axis=0)
                dr = pos - com

            radii = np.linalg.norm(dr, axis=1)
            if len(radii) < 5:
                continue

            # Bin into 10 radial bins
            r_max = radii.max()
            if r_max == 0:
                continue
            bins = np.linspace(0, r_max, 11)
            counts, _ = np.histogram(radii, bins=bins)
            # Normalize by ring area
            areas = np.pi * (bins[1:] ** 2 - bins[:-1] ** 2)
            areas = np.where(areas == 0, 1, areas)
            density = counts / areas

            if density.max() > 0:
                densities_center.append(density[0])
                densities_peak.append(density.max())

        if not densities_peak or np.mean(densities_peak) == 0:
            return 1.0

        return float(np.mean(densities_center) / np.mean(densities_peak))

    @staticmethod
    def _compute_rotation_coherence(state_history: list[dict[str, Any]]) -> float:
        """Low variance of per-agent angular velocities → high coherence."""
        from epc.metrics.collective_motion import _periodic_com

        if len(state_history) < 3:
            return 0.0

        angular_vel_stds = []
        for t in range(1, len(state_history)):
            s0 = state_history[t - 1]
            s1 = state_history[t]
            pos = s1["positions"]
            box = s1.get("box_size", None)

            if box:
                com = _periodic_com(pos, box)
                dr = pos - com
                dr = dr - box * np.round(dr / box)
            else:
                com = pos.mean(axis=0)
                dr = pos - com

            angles = np.arctan2(dr[:, 1], dr[:, 0])

            pos0 = s0["positions"]
            if box:
                com0 = _periodic_com(pos0, box)
                dr0 = pos0 - com0
                dr0 = dr0 - box * np.round(dr0 / box)
            else:
                com0 = pos0.mean(axis=0)
                dr0 = pos0 - com0
            angles0 = np.arctan2(dr0[:, 1], dr0[:, 0])

            d_angle = angles - angles0
            d_angle = (d_angle + np.pi) % (2 * np.pi) - np.pi
            angular_vel_stds.append(float(np.std(d_angle)))

        if not angular_vel_stds:
            return 0.0
        # Coherence: inverse of mean std (normalized to [0,1] range)
        mean_std = float(np.mean(angular_vel_stds))
        return float(1.0 / (1.0 + mean_std))

    def _run_null_model(
        self,
        state_history: list[dict[str, Any]],
        primary_result: dict[str, float],
        timescale: float,
    ) -> tuple[float, NullType, dict[str, float]]:
        """Heading-shuffle null: permute headings, recompute |L|."""
        rng = np.random.default_rng(0)
        observed_L = primary_result.get("angular_momentum_abs_mean", 0)

        if not state_history or "velocities" not in state_history[0]:
            return 1.0, NullType.SHUFFLE, {"mean": 0.0, "std": 0.0}

        n_frames = len(state_history)
        n_sample = min(10, n_frames)
        sample_idx = np.linspace(0, n_frames - 1, n_sample, dtype=int)

        null_Ls: list[float] = []
        for _ in range(self.n_permutations):
            frame_Ls = []
            for idx in sample_idx:
                state = state_history[int(idx)]
                # Shuffle velocity assignment
                v = state["velocities"]
                perm = rng.permutation(len(v))
                shuffled_state = dict(state)
                shuffled_state["velocities"] = v[perm]
                L = abs(AngularMomentum._compute_L(shuffled_state))
                frame_Ls.append(L)
            null_Ls.append(float(np.mean(frame_Ls)))

        null_arr = np.array(null_Ls)
        p_value = float(np.mean(null_arr >= observed_L))
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
        """|L| > 0.5 AND R < 0.3."""
        L = primary_result.get("angular_momentum_abs_mean", 0)
        R = secondary_result.get("group_speed_ratio", 1.0)
        return L > 0.5 and R < 0.3 and null_p < 0.01

    def _check_exclusions(
        self,
        state_history: list[dict[str, Any]],
        model_metadata: dict[str, Any] | None,
        timescale: float,
    ) -> tuple[list[str], dict[str, str]]:
        results: dict[str, str] = {}

        # P5: R > 0.5 → translational motion present → not_excluded (exclude P6)
        if state_history and "velocities" in state_history[0]:
            gsr = self._gsr_metric.compute(state_history)
            R = gsr.get("group_speed_ratio", gsr.get("group_speed_ratio_mean", 0))
            results["P5"] = "not_excluded" if R > 0.5 else "excluded"
        else:
            results["P5"] = "inconclusive"

        # P12: nontransitive species spirals vs physical rotation
        if model_metadata:
            model_class = model_metadata.get("model_class", "")
            if "nontransitive" in model_class or "rps" in model_class:
                results["P12"] = "not_excluded"
            else:
                results["P12"] = "excluded"
        else:
            results["P12"] = "inconclusive"

        return (["P5", "P12"], results)
