"""
Collective motion metrics for Cluster B patterns (P5-P8).

Implements the metrics specified in the P5 and P6 detector cards:
- Polarization order parameter φ (P5 primary)
- Group speed ratio R (P5 secondary, P6 secondary)
- Angular momentum L (P6 primary, P5 exclusion check)
- Heading autocorrelation / directional persistence (P5 secondary)

All metrics operate on state histories (list of dicts with 'positions',
'velocities', and optionally 'headings' keys).
"""

from __future__ import annotations

import numpy as np
from numpy.typing import NDArray
from typing import Any, Optional


class PolarizationMetric:
    """Polarization order parameter φ — primary metric for P5 flocking.

    φ(t) = |(1/N) Σ v̂_i(t)|

    where v̂_i = v_i / |v_i| is the unit heading vector.

    φ = 1: perfect alignment (all moving in same direction).
    φ ≈ 1/√N: disordered (random headings).

    Returns time series and summary statistics.
    """

    @staticmethod
    def compute_instant(state: dict[str, Any]) -> float:
        """Compute φ for a single state snapshot."""
        velocities = state["velocities"]
        speeds = np.linalg.norm(velocities, axis=1)
        mask = speeds > 1e-12  # exclude stationary agents
        if mask.sum() == 0:
            return 0.0
        unit_v = velocities[mask] / speeds[mask, np.newaxis]
        mean_v = np.mean(unit_v, axis=0)
        return float(np.linalg.norm(mean_v))

    @staticmethod
    def compute(history: list[dict[str, Any]]) -> dict[str, Any]:
        """Compute polarization time series and statistics.

        Args:
            history: List of state dicts with 'velocities' key.

        Returns:
            Dict with keys:
                phi_series: (T,) time series of φ
                phi_mean: mean φ over trajectory
                phi_std: std of φ over trajectory
                phi_final: φ at last timestep
                n_steps: number of timesteps
        """
        phi_series = np.array([
            PolarizationMetric.compute_instant(s) for s in history
        ])
        return {
            "phi_series": phi_series,
            "phi_mean": float(np.mean(phi_series)),
            "phi_std": float(np.std(phi_series)),
            "phi_final": float(phi_series[-1]),
            "n_steps": len(phi_series),
        }


class GroupSpeedRatioMetric:
    """Group speed ratio R = |V_cm| / ⟨|v_i|⟩ — secondary metric for P5.

    R = 1: all moving in same direction at same speed.
    R ≈ 0: random motion, COM velocity averages to zero.
    R = φ exactly for constant-speed models.

    For P5 confirmation: R > 0.5 required.
    For P6 exclusion: R < 0.3 indicates milling.
    """

    @staticmethod
    def compute_instant(state: dict[str, Any]) -> float:
        """Compute R for a single state snapshot."""
        velocities = state["velocities"]
        v_cm = np.mean(velocities, axis=0)
        v_cm_mag = np.linalg.norm(v_cm)
        mean_speed = np.mean(np.linalg.norm(velocities, axis=1))
        if mean_speed < 1e-12:
            return 0.0
        return float(v_cm_mag / mean_speed)

    @staticmethod
    def compute(history: list[dict[str, Any]]) -> dict[str, Any]:
        """Compute R time series and statistics."""
        R_series = np.array([
            GroupSpeedRatioMetric.compute_instant(s) for s in history
        ])
        return {
            "R_series": R_series,
            "R_mean": float(np.mean(R_series)),
            "R_std": float(np.std(R_series)),
            "n_steps": len(R_series),
        }


class AngularMomentumMetric:
    """Angular momentum L — primary metric for P6 milling, exclusion for P5.

    L(t) = (1/N) Σ (r̂_i × v̂_i)

    where r̂_i = (r_i - r_cm) / |r_i - r_cm| and v̂_i = v_i / |v_i|.
    Cross product in 2D gives a scalar.

    |L| ≈ 1: milling (coherent rotation).
    |L| ≈ 0: no rotational order.

    For P6 screening: |L| > 0.3.
    For P5 exclusion of P6: |L| > 0.5 → milling, not flocking.
    """

    @staticmethod
    def compute_instant(
        state: dict[str, Any],
        box_size: Optional[float] = None,
    ) -> float:
        """Compute angular momentum L for a single state snapshot.

        Args:
            state: State dict with 'positions' and 'velocities'.
            box_size: If provided, use periodic-aware COM calculation.

        Returns:
            Angular momentum L ∈ [-1, 1].
        """
        positions = state["positions"]
        velocities = state["velocities"]
        N = len(positions)

        # --- Center of mass ---
        if box_size is not None:
            # Periodic-aware COM via circular mean
            theta_x = 2 * np.pi * positions[:, 0] / box_size
            theta_y = 2 * np.pi * positions[:, 1] / box_size
            com_x = box_size / (2 * np.pi) * np.arctan2(
                np.mean(np.sin(theta_x)), np.mean(np.cos(theta_x))
            ) % box_size
            com_y = box_size / (2 * np.pi) * np.arctan2(
                np.mean(np.sin(theta_y)), np.mean(np.cos(theta_y))
            ) % box_size
            dx = positions[:, 0] - com_x
            dy = positions[:, 1] - com_y
            dx -= box_size * np.round(dx / box_size)
            dy -= box_size * np.round(dy / box_size)
        else:
            com = np.mean(positions, axis=0)
            dx = positions[:, 0] - com[0]
            dy = positions[:, 1] - com[1]

        # Unit radial vectors
        dist = np.sqrt(dx**2 + dy**2)
        dist = np.maximum(dist, 1e-12)
        rx_hat = dx / dist
        ry_hat = dy / dist

        # Unit heading vectors
        speeds = np.linalg.norm(velocities, axis=1)
        mask = speeds > 1e-12
        if mask.sum() == 0:
            return 0.0

        vx_hat = np.where(mask, velocities[:, 0] / np.maximum(speeds, 1e-12), 0.0)
        vy_hat = np.where(mask, velocities[:, 1] / np.maximum(speeds, 1e-12), 0.0)

        # 2D cross product: r̂ × v̂ = rx*vy - ry*vx
        cross = rx_hat * vy_hat - ry_hat * vx_hat
        return float(np.mean(cross))

    @staticmethod
    def compute(
        history: list[dict[str, Any]],
        box_size: Optional[float] = None,
    ) -> dict[str, Any]:
        """Compute angular momentum time series and statistics."""
        L_series = np.array([
            AngularMomentumMetric.compute_instant(s, box_size)
            for s in history
        ])
        return {
            "L_series": L_series,
            "L_mean": float(np.mean(L_series)),
            "L_abs_mean": float(np.mean(np.abs(L_series))),
            "L_std": float(np.std(L_series)),
            "n_steps": len(L_series),
        }


class HeadingAutocorrelationMetric:
    """Heading autocorrelation — directional persistence metric for P5.

    C(τ) = ⟨(1/N) Σ v̂_i(t) · v̂_i(t+τ)⟩_t

    Measures how long individual agents maintain their direction.
    High C(τ) at large τ indicates persistent directed motion.
    """

    @staticmethod
    def compute(
        history: list[dict[str, Any]],
        max_lag: int = 50,
    ) -> dict[str, Any]:
        """Compute heading autocorrelation function.

        Args:
            history: State history with 'velocities' key.
            max_lag: Maximum lag in timesteps.

        Returns:
            Dict with autocorrelation function and decay time.
        """
        T = len(history)
        max_lag = min(max_lag, T // 2)

        # Extract unit headings time series: (T, N, 2)
        unit_headings = []
        for state in history:
            v = state["velocities"]
            speeds = np.linalg.norm(v, axis=1, keepdims=True)
            speeds = np.maximum(speeds, 1e-12)
            unit_headings.append(v / speeds)
        unit_headings = np.array(unit_headings)  # (T, N, 2)

        # Compute autocorrelation at each lag
        C = np.zeros(max_lag + 1)
        for tau in range(max_lag + 1):
            # Dot product of unit headings at t and t+τ, averaged over agents and time
            dots = np.sum(
                unit_headings[:T - tau] * unit_headings[tau:],
                axis=2,  # dot product over (x,y)
            )  # (T-tau, N)
            C[tau] = np.mean(dots)

        # Find decorrelation time (first crossing below 1/e)
        decay_threshold = C[0] / np.e
        decay_time = max_lag  # default: never decorrelates
        for tau in range(1, max_lag + 1):
            if C[tau] < decay_threshold:
                # Linear interpolation
                if C[tau - 1] > decay_threshold:
                    frac = (C[tau - 1] - decay_threshold) / (C[tau - 1] - C[tau])
                    decay_time = tau - 1 + frac
                else:
                    decay_time = float(tau)
                break

        return {
            "autocorrelation": C,
            "lags": np.arange(max_lag + 1),
            "decay_time": float(decay_time),
            "C_0": float(C[0]),
        }


class HeadingDistributionMetric:
    """Heading distribution analysis — used for P7 (lane) exclusion in P5.

    Checks whether the heading distribution is unimodal (flocking)
    or bimodal/antiparallel (lanes).

    Uses the nematic order parameter S = |(1/N) Σ exp(2iθ_j)| to detect
    antiparallel alignment. S ≈ 1 when headings cluster at θ and θ+π
    (antiparallel). When combined with low polar order φ, this indicates
    counterflow / lanes rather than flocking.
    """

    @staticmethod
    def compute(
        history: list[dict[str, Any]],
        window_fraction: float = 0.5,
    ) -> dict[str, Any]:
        """Analyze heading distribution over late trajectory.

        Args:
            history: State history with 'velocities' key.
            window_fraction: Fraction of trajectory to analyze (from end).

        Returns:
            Dict with heading statistics and bimodality test.
        """
        # Collect headings from late trajectory
        start_idx = max(0, int(len(history) * (1 - window_fraction)))
        all_headings = []
        for state in history[start_idx:]:
            v = state["velocities"]
            speeds = np.linalg.norm(v, axis=1)
            mask = speeds > 1e-12
            headings = np.arctan2(v[mask, 1], v[mask, 0])
            all_headings.extend(headings)

        all_headings = np.array(all_headings)

        # Polar order (standard polarization resultant length)
        mean_sin = np.mean(np.sin(all_headings))
        mean_cos = np.mean(np.cos(all_headings))
        mean_heading = np.arctan2(mean_sin, mean_cos)
        resultant_length = np.sqrt(mean_sin**2 + mean_cos**2)

        # Nematic order parameter S = |(1/N) Σ exp(2iθ)|
        # S ≈ 1 for antiparallel alignment (θ and θ+π)
        # S ≈ 1 for polar alignment too, so must combine with polar order
        doubled = 2 * all_headings
        nematic_sin = np.mean(np.sin(doubled))
        nematic_cos = np.mean(np.cos(doubled))
        nematic_order = np.sqrt(nematic_sin**2 + nematic_cos**2)

        # Antiparallel detection: high nematic order + low polar order
        # Thresholds: nematic S > 0.5 AND polar φ < 0.3
        is_antiparallel = (nematic_order > 0.5) and (resultant_length < 0.3)

        # Circular variance (1 - resultant_length)
        circular_variance = 1 - resultant_length

        return {
            "mean_heading": float(mean_heading),
            "resultant_length": float(resultant_length),
            "nematic_order": float(nematic_order),
            "circular_variance": float(circular_variance),
            "antiparallel_fraction": float(nematic_order if is_antiparallel else 0.0),
            "is_bimodal_antiparallel": is_antiparallel,
            "n_samples": len(all_headings),
        }
