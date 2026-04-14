"""Unit tests for metrics."""

import numpy as np
import pytest


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
            history.append({"grid": grid, "grid_dims": (20, 20)})

        result = metric.compute(history)
        assert result["n_wavefront_frames"] > 5
        assert result["wavefront_speed_cv"] < 0.5

    def test_no_wavefront(self):
        """All-rest grid should report zero speed."""
        from epc.metrics.excitable_waves import WavefrontSpeed
        metric = WavefrontSpeed()
        history = [{"grid": np.zeros((10, 10), dtype=int), "grid_dims": (10, 10)} for _ in range(10)]
        result = metric.compute(history)
        assert result["n_wavefront_frames"] == 0


class TestPolarization:
    def test_perfect_alignment(self):
        """All agents heading right → φ = 1."""
        from epc.metrics.collective_motion import Polarization
        metric = Polarization()
        v = np.column_stack([np.ones(50), np.zeros(50)])
        history = [{"velocities": v}]
        result = metric.compute(history, timestep=0)
        assert abs(result["polarization"] - 1.0) < 1e-10

    def test_random_headings(self):
        """Random headings → φ ≈ 1/√N."""
        from epc.metrics.collective_motion import Polarization
        metric = Polarization()
        rng = np.random.default_rng(42)
        n = 1000
        angles = rng.uniform(0, 2 * np.pi, n)
        v = np.column_stack([np.cos(angles), np.sin(angles)])
        history = [{"velocities": v}]
        result = metric.compute(history, timestep=0)
        expected = 1.0 / np.sqrt(n)
        assert result["polarization"] < 3 * expected  # should be near 1/sqrt(N)


class TestAngularMomentum:
    def test_perfect_mill(self):
        """Agents arranged in circle, velocities tangent → |L| ≈ 1."""
        from epc.metrics.collective_motion import AngularMomentum
        metric = AngularMomentum()
        n = 100
        angles = np.linspace(0, 2 * np.pi, n, endpoint=False)
        r = 5.0
        positions = np.column_stack([r * np.cos(angles) + 15, r * np.sin(angles) + 15])
        # Tangent velocities (counterclockwise)
        velocities = np.column_stack([-np.sin(angles), np.cos(angles)])
        history = [{"positions": positions, "velocities": velocities, "box_size": 30.0}]
        result = metric.compute(history, timestep=0)
        assert abs(result["angular_momentum"]) > 0.9

    def test_zero_angular_momentum(self):
        """Random velocities → L ≈ 0."""
        from epc.metrics.collective_motion import AngularMomentum
        metric = AngularMomentum()
        rng = np.random.default_rng(42)
        n = 500
        positions = rng.uniform(0, 20, (n, 2))
        angles = rng.uniform(0, 2 * np.pi, n)
        velocities = np.column_stack([np.cos(angles), np.sin(angles)])
        history = [{"positions": positions, "velocities": velocities, "box_size": 20.0}]
        result = metric.compute(history, timestep=0)
        assert abs(result["angular_momentum"]) < 0.2


class TestGroupSpeedRatio:
    def test_aligned_group(self):
        """All moving same direction → R = 1."""
        from epc.metrics.collective_motion import GroupSpeedRatio
        metric = GroupSpeedRatio()
        v = np.column_stack([np.ones(50), np.zeros(50)])
        history = [{"velocities": v}]
        result = metric.compute(history, timestep=0)
        assert abs(result["group_speed_ratio"] - 1.0) < 1e-10

    def test_opposing_groups(self):
        """Half left, half right → R ≈ 0."""
        from epc.metrics.collective_motion import GroupSpeedRatio
        metric = GroupSpeedRatio()
        v = np.zeros((100, 2))
        v[:50, 0] = 1.0
        v[50:, 0] = -1.0
        history = [{"velocities": v}]
        result = metric.compute(history, timestep=0)
        assert abs(result["group_speed_ratio"]) < 0.01


class TestSpiralTipDetector:
    def test_no_tips_uniform(self):
        """Uniform grid should have no spiral tips."""
        from epc.metrics.excitable_waves import SpiralTipDetector
        metric = SpiralTipDetector()
        grid = np.zeros((20, 20), dtype=int)
        history = [{"grid": grid, "grid_dims": (20, 20)}]
        result = metric.compute(history, n_states=5)
        assert result["n_spiral_tips_max"] == 0


class TestWaveSourceCount:
    def test_no_sources_quiescent(self):
        """All-rest grid should have no wave sources."""
        from epc.metrics.excitable_waves import WaveSourceCount
        metric = WaveSourceCount()
        history = [{"grid": np.zeros((10, 10), dtype=int), "grid_dims": (10, 10)} for _ in range(10)]
        result = metric.compute(history)
        assert result["wave_source_count"] == 0
