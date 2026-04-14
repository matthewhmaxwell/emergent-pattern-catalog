"""Collective motion metrics for P5 (flocking) and P6 (milling) detection.

- Polarization: φ = |(1/N) Σ v̂_i| — global alignment order parameter
- AngularMomentum: L = (1/N) Σ (r̂_i × v̂_i) — rotational order parameter
- GroupSpeedRatio: R = |V_cm| / ⟨|v_i|⟩ — translational efficiency
"""

from __future__ import annotations

from typing import Any

import numpy as np

from epc.base_metric import BaseMetric


def _periodic_com(positions: np.ndarray, box_size: float) -> np.ndarray:
    """Center of mass in periodic domain via circular mean."""
    com = np.empty(2)
    for d in range(2):
        theta = 2.0 * np.pi * positions[:, d] / box_size
        com[d] = box_size * np.arctan2(
            np.mean(np.sin(theta)), np.mean(np.cos(theta))
        ) / (2.0 * np.pi)
        com[d] = com[d] % box_size
    return com


def _unit_vectors(velocities: np.ndarray) -> np.ndarray:
    """Normalize velocity vectors, skipping zero-magnitude agents."""
    speeds = np.linalg.norm(velocities, axis=1, keepdims=True)
    speeds = np.where(speeds == 0, 1.0, speeds)  # avoid div by zero
    return velocities / speeds


class Polarization(BaseMetric):
    """Polarization order parameter φ = |(1/N) Σ v̂_i|.

    φ ≈ 1 → aligned (flocking). φ ≈ 1/√N → disordered.
    """

    def __init__(self) -> None:
        super().__init__(name="polarization")

    def required_keys(self) -> list[str]:
        return ["velocities"]

    def compute(
        self, state_history: list[dict[str, Any]], **kwargs: Any
    ) -> dict[str, Any]:
        timestep: int = kwargs.get("timestep", -1)

        if timestep != -1 or len(state_history) == 1:
            # Single-frame computation
            v = state_history[timestep]["velocities"]
            v_hat = _unit_vectors(v)
            phi = float(np.linalg.norm(v_hat.mean(axis=0)))
            return {self.name: phi}

        # Full timeseries
        ts = self.compute_timeseries(state_history)
        return {
            self.name: float(ts.mean()),
            "polarization_mean": float(ts.mean()),
            "polarization_std": float(ts.std()),
            "polarization_min": float(ts.min()),
            "polarization_max": float(ts.max()),
            "polarization_timeseries": ts,
        }

    def compute_timeseries(
        self, state_history: list[dict[str, Any]], **kwargs: Any
    ) -> np.ndarray:
        values = np.empty(len(state_history))
        for i, state in enumerate(state_history):
            v_hat = _unit_vectors(state["velocities"])
            values[i] = np.linalg.norm(v_hat.mean(axis=0))
        return values


class AngularMomentum(BaseMetric):
    """Angular momentum L = (1/N) Σ (r̂_i × v̂_i) relative to COM.

    |L| ≈ 1 → milling/vortex. |L| ≈ 0 → no rotational order.
    Uses periodic COM computation for systems with periodic boundaries.
    """

    def __init__(self) -> None:
        super().__init__(name="angular_momentum")

    def required_keys(self) -> list[str]:
        return ["positions", "velocities"]

    def compute(
        self, state_history: list[dict[str, Any]], **kwargs: Any
    ) -> dict[str, Any]:
        timestep: int = kwargs.get("timestep", -1)

        if timestep != -1 or len(state_history) == 1:
            state = state_history[timestep]
            L = self._compute_L(state)
            return {self.name: L, "angular_momentum_abs": abs(L)}

        ts = self.compute_timeseries(state_history)
        return {
            self.name: float(ts.mean()),
            "angular_momentum_abs_mean": float(np.abs(ts).mean()),
            "angular_momentum_abs_std": float(np.abs(ts).std()),
            "angular_momentum_timeseries": ts,
        }

    def compute_timeseries(
        self, state_history: list[dict[str, Any]], **kwargs: Any
    ) -> np.ndarray:
        values = np.empty(len(state_history))
        for i, state in enumerate(state_history):
            values[i] = self._compute_L(state)
        return values

    @staticmethod
    def _compute_L(state: dict[str, Any]) -> float:
        """Compute angular momentum for one frame."""
        positions = state["positions"]
        velocities = state["velocities"]
        box_size = state.get("box_size", None)
        n = len(positions)

        # Compute COM
        if box_size is not None:
            com = _periodic_com(positions, box_size)
        else:
            com = positions.mean(axis=0)

        # Displacement from COM with periodic minimum image
        dr = positions - com
        if box_size is not None:
            dr = dr - box_size * np.round(dr / box_size)

        # Radial unit vectors
        r_mag = np.linalg.norm(dr, axis=1, keepdims=True)
        r_mag = np.where(r_mag == 0, 1.0, r_mag)
        r_hat = dr / r_mag

        # Velocity unit vectors
        v_hat = _unit_vectors(velocities)

        # 2D cross product: r̂_x·v̂_y - r̂_y·v̂_x
        cross = r_hat[:, 0] * v_hat[:, 1] - r_hat[:, 1] * v_hat[:, 0]

        return float(cross.mean())


class GroupSpeedRatio(BaseMetric):
    """Group speed ratio R = |V_cm| / ⟨|v_i|⟩.

    R ≈ 1 → coherent translational motion (flocking).
    R ≈ 0 → no net displacement despite individual motion (milling or disorder).
    """

    def __init__(self) -> None:
        super().__init__(name="group_speed_ratio")

    def required_keys(self) -> list[str]:
        return ["velocities"]

    def compute(
        self, state_history: list[dict[str, Any]], **kwargs: Any
    ) -> dict[str, Any]:
        timestep: int = kwargs.get("timestep", -1)

        if timestep != -1 or len(state_history) == 1:
            state = state_history[timestep]
            R = self._compute_R(state)
            return {self.name: R}

        values = np.array([self._compute_R(s) for s in state_history])
        return {
            self.name: float(values.mean()),
            "group_speed_ratio_mean": float(values.mean()),
            "group_speed_ratio_std": float(values.std()),
        }

    @staticmethod
    def _compute_R(state: dict[str, Any]) -> float:
        v = state["velocities"]
        v_cm = v.mean(axis=0)
        v_cm_mag = float(np.linalg.norm(v_cm))
        mean_speed = float(np.linalg.norm(v, axis=1).mean())
        if mean_speed == 0:
            return 0.0
        return v_cm_mag / mean_speed
