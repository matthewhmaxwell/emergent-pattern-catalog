"""Unit tests for metric implementations."""

import numpy as np
import pytest


# --- Excitable wave metrics ---

class TestWavefrontSpeed:
    def test_constant_speed_wavefront(self):
        """Synthetic wavefront moving at 1 cell/step should report speed ~1."""
        from epc.metrics.excitable_waves import WavefrontSpeed
        metric = WavefrontSpeed()

        # Create a 20x20 grid with a wavefront moving rightward
        history = []
        for t in range(20):
            grid = np.zeros((20, 20), dtype=int)
            col = t % 20
            grid[:, col] = 1  # excited wavefront column
            if t > 0:
                prev_col = (t - 1) % 20
                grid[:, prev_col] = 2  # refractory
            history.append({"grid": grid, "grid_dims": (20, 20), "n_states": 3})

        result = metric.compute(history)
        assert result["wavefront_count_mean"] > 5
        assert abs(result["wavefront_speed_mean"] - 1.0) < 0.2

    def test_no_wavefront(self):
        """All-rest grid should report zero speed."""
        from epc.metrics.excitable_waves import WavefrontSpeed
        metric = WavefrontSpeed()
        history = [
            {"grid": np.zeros((10, 10), dtype=int), "grid_dims": (10, 10), "n_states": 3}
            for _ in range(10)
        ]
        result = metric.compute(history)
        assert result["wavefront_count_mean"] == 0


# --- Collective motion metrics ---

class TestPolarization:
    def test_perfect_alignment(self):
        """All agents heading right -> phi = 1."""
        from epc.metrics.collective_motion import PolarizationMetric
        N = 50
        pos = np.random.default_rng(0).random((N, 2)) * 10
        vel = np.tile([1.0, 0.0], (N, 1))
        history = [{"positions": pos, "velocities": vel}] * 10
        result = PolarizationMetric.compute(history)
        assert abs(result["phi_mean"] - 1.0) < 0.01

    def test_random_headings(self):
        """Random headings -> phi ~ 1/sqrt(N)."""
        from epc.metrics.collective_motion import PolarizationMetric
        N = 200
        rng = np.random.default_rng(42)
        pos = rng.random((N, 2)) * 10
        theta = rng.uniform(0, 2 * np.pi, N)
        vel = np.column_stack([np.cos(theta), np.sin(theta)])
        history = [{"positions": pos, "velocities": vel}] * 10
        result = PolarizationMetric.compute(history)
        expected = 1 / np.sqrt(N)
        assert result["phi_mean"] < 0.3, \
            f"Expected phi ~ {expected:.3f}, got {result['phi_mean']:.3f}"


class TestAngularMomentum:
    def test_perfect_mill(self):
        """Agents arranged in circle, velocities tangent -> |L| ~ 1."""
        from epc.metrics.collective_motion import AngularMomentumMetric
        N = 50
        theta = np.linspace(0, 2 * np.pi, N, endpoint=False)
        r = 5.0
        pos = np.column_stack([r * np.cos(theta), r * np.sin(theta)]) + 10
        vel = np.column_stack([-np.sin(theta), np.cos(theta)])
        history = [{"positions": pos, "velocities": vel}] * 10
        result = AngularMomentumMetric.compute(history)
        assert result["L_abs_mean"] > 0.9, \
            f"Expected |L| ~ 1, got {result['L_abs_mean']:.3f}"

    def test_zero_angular_momentum(self):
        """Random velocities -> L ~ 0."""
        from epc.metrics.collective_motion import AngularMomentumMetric
        N = 200
        rng = np.random.default_rng(42)
        pos = rng.random((N, 2)) * 10
        theta = rng.uniform(0, 2 * np.pi, N)
        vel = np.column_stack([np.cos(theta), np.sin(theta)])
        history = [{"positions": pos, "velocities": vel}] * 10
        result = AngularMomentumMetric.compute(history)
        assert result["L_abs_mean"] < 0.3, \
            f"Expected |L| ~ 0, got {result['L_abs_mean']:.3f}"


class TestGroupSpeedRatio:
    def test_aligned_group(self):
        """All moving same direction -> R = 1."""
        from epc.metrics.collective_motion import GroupSpeedRatioMetric
        N = 50
        pos = np.random.default_rng(0).random((N, 2)) * 10
        vel = np.tile([1.0, 0.0], (N, 1))
        history = [{"positions": pos, "velocities": vel}] * 10
        result = GroupSpeedRatioMetric.compute(history)
        assert abs(result["R_mean"] - 1.0) < 0.01

    def test_opposing_groups(self):
        """Half left, half right -> R ~ 0."""
        from epc.metrics.collective_motion import GroupSpeedRatioMetric
        N = 100
        pos = np.random.default_rng(0).random((N, 2)) * 10
        vel = np.zeros((N, 2))
        vel[:N // 2, 0] = 1.0
        vel[N // 2:, 0] = -1.0
        history = [{"positions": pos, "velocities": vel}] * 10
        result = GroupSpeedRatioMetric.compute(history)
        assert result["R_mean"] < 0.1, \
            f"Expected R ~ 0, got {result['R_mean']:.3f}"


class TestSpiralTipDetector:
    def test_no_tips_uniform(self):
        """Uniform grid should have no spiral tips."""
        from epc.metrics.excitable_waves import SpiralTipDetector
        metric = SpiralTipDetector()
        grid = np.zeros((20, 20), dtype=int)
        # n_states must be in the state dict
        history = [{"grid": grid, "grid_dims": (20, 20), "n_states": 5}]
        result = metric.compute(history)
        assert result["tip_count_mean"] == 0


class TestWavePersistence:
    def test_no_persistence_quiescent(self):
        """All-rest grid should have zero persistence."""
        from epc.metrics.excitable_waves import WavePersistence
        metric = WavePersistence()
        history = [
            {"grid": np.zeros((10, 10), dtype=int), "grid_dims": (10, 10), "n_states": 3}
            for _ in range(20)
        ]
        result = metric.compute(history)
        assert result["active_fraction"] == 0.0
