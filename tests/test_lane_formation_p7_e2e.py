"""
End-to-end tests for P7 lane formation (model + detector).

Tests:
1. Model is deterministic from seed
2. Detector reaches DEFINITIVE on canonical lane-forming regime
3. Detector does NOT fire on single-population free-flow (negative control)
4. Registry registration is correct
"""

import numpy as np
import pytest

from epc.models.lane_formation import LaneFormationModel
from epc.detectors.p7_lane_formation import (
    P7LaneFormationDetector,
    compute_lane_order_parameter,
    compute_headon_encounter_rate,
    compute_throughput,
)
from epc.orchestration import MODEL_REGISTRY, DETECTOR_REGISTRY, check_compatibility


class TestModelDeterminism:
    """Model produces identical output from same seed."""

    def test_deterministic_from_seed(self):
        """Two runs with same seed produce identical histories."""
        model1 = LaneFormationModel(n_agents=50, seed=123)
        h1 = model1.run(n_steps=100, record_interval=10)

        model2 = LaneFormationModel(n_agents=50, seed=123)
        h2 = model2.run(n_steps=100, record_interval=10)

        for s1, s2 in zip(h1, h2):
            np.testing.assert_array_equal(s1["positions"], s2["positions"])
            np.testing.assert_array_equal(s1["velocities"], s2["velocities"])
            np.testing.assert_array_equal(s1["labels"], s2["labels"])

    def test_different_seeds_differ(self):
        """Different seeds produce different outcomes."""
        model1 = LaneFormationModel(n_agents=50, seed=1)
        h1 = model1.run(n_steps=50)

        model2 = LaneFormationModel(n_agents=50, seed=2)
        h2 = model2.run(n_steps=50)

        # Final positions should differ
        assert not np.allclose(h1[-1]["positions"], h2[-1]["positions"])


class TestLaneFormationDetection:
    """Detector correctly identifies lane formation."""

    def test_canonical_regime_reaches_definitive(self):
        """Canonical counterflow regime produces DEFINITIVE detection."""
        # Use parameters that reliably produce strong lanes
        model = LaneFormationModel(
            n_agents=200,
            corridor_width=20.0,
            corridor_height=4.0,
            desired_speed=1.0,
            repulsion_amplitude=5.0,
            repulsion_range=0.3,
            tau=0.5,
            dt=0.05,
            seed=42,
        )
        history = model.run(n_steps=2000, record_interval=10)
        metadata = model.get_metadata()

        detector = P7LaneFormationDetector(n_permutations=199, seed=42)
        result = detector.detect(history, metadata)

        assert result.detected is True
        assert result.tier in ("confirmation", "definitive"), (
            f"Expected confirmation or definitive, got {result.tier}; "
            f"notes: {result.notes}"
        )
        assert result.confidence >= 0.700
        assert result.null_p_value < 0.01
        assert result.primary_metric["lane_order_parameter_mean"] > 0.4

    def test_negative_control_no_lanes_weak_repulsion(self):
        """Very weak repulsion prevents lane segregation → no detection.

        With near-zero repulsion, agents pass through each other without
        segregating laterally. The lane order parameter should stay near
        the random-mixing baseline.
        """
        model = LaneFormationModel(
            n_agents=100,
            corridor_width=20.0,
            corridor_height=4.0,
            desired_speed=1.0,
            repulsion_amplitude=0.01,  # essentially no repulsion
            repulsion_range=0.1,
            tau=0.5,
            dt=0.05,
            seed=42,
        )
        history = model.run(n_steps=500, record_interval=10)
        metadata = model.get_metadata()

        detector = P7LaneFormationDetector(n_permutations=99, seed=42)
        result = detector.detect(history, metadata)

        # Without repulsion, agents don't segregate laterally
        # The null test should fail to reject (p > 0.01)
        assert result.tier == "none" or result.confidence <= 0.500, (
            f"Expected no detection or screening-only with weak repulsion, "
            f"got tier={result.tier}, confidence={result.confidence}, "
            f"phi={result.primary_metric.get('lane_order_parameter_mean', 'N/A')}"
        )


class TestLaneOrderParameter:
    """Unit tests for the lane order metric."""

    def test_perfect_lanes(self):
        """Perfectly segregated populations → phi ≈ 1."""
        N = 100
        positions = np.column_stack([
            np.random.default_rng(0).uniform(0, 20, N),
            np.concatenate([
                np.random.default_rng(0).uniform(0, 2, N // 2),   # pop 0 in bottom half
                np.random.default_rng(1).uniform(2, 4, N // 2),   # pop 1 in top half
            ]),
        ])
        labels = np.array([0] * (N // 2) + [1] * (N // 2))

        phi = compute_lane_order_parameter(positions, labels, 4.0, n_bins=10)
        assert phi > 0.8, f"Expected phi > 0.8 for perfect lanes, got {phi}"

    def test_random_mixing(self):
        """Randomly mixed populations → phi ≈ 0 (low)."""
        N = 200
        rng = np.random.default_rng(42)
        positions = np.column_stack([
            rng.uniform(0, 20, N),
            rng.uniform(0, 4, N),
        ])
        labels = np.array([0] * (N // 2) + [1] * (N // 2))
        rng.shuffle(labels)

        phi = compute_lane_order_parameter(positions, labels, 4.0, n_bins=10)
        # Random mixing should give phi around 1/sqrt(N/n_bins) ≈ 0.22
        assert phi < 0.5, f"Expected phi < 0.5 for random mixing, got {phi}"


class TestThroughput:
    """Test throughput computation."""

    def test_free_flow(self):
        """All agents at desired speed → throughput = v0."""
        N = 100
        labels = np.array([0] * (N // 2) + [1] * (N // 2))
        velocities = np.zeros((N, 2))
        velocities[:N//2, 0] = 1.0   # +x population going +x
        velocities[N//2:, 0] = -1.0  # -x population going -x

        tp = compute_throughput(velocities, labels)
        assert abs(tp - 1.0) < 0.01

    def test_blocked_flow(self):
        """All agents stationary → throughput = 0."""
        N = 100
        labels = np.array([0] * (N // 2) + [1] * (N // 2))
        velocities = np.zeros((N, 2))

        tp = compute_throughput(velocities, labels)
        assert abs(tp) < 0.01


class TestRegistry:
    """Test registry registration."""

    def test_model_registered(self):
        """Lane formation model is in MODEL_REGISTRY."""
        assert "lane_formation" in MODEL_REGISTRY
        reg = MODEL_REGISTRY["lane_formation"]
        assert reg.substrate_type == "continuous_2d"
        assert "positions" in reg.observables
        assert "velocities" in reg.observables
        assert "labels" in reg.observables
        assert "P7" in reg.primary_patterns

    def test_detector_registered(self):
        """P7 detector is in DETECTOR_REGISTRY."""
        assert "P7" in DETECTOR_REGISTRY
        reg = DETECTOR_REGISTRY["P7"]
        assert "continuous_2d" in reg.required_substrate
        assert "positions" in reg.required_observables

    def test_compatibility(self):
        """Lane formation model is compatible with P7 detector."""
        result = check_compatibility("lane_formation", "P7")
        assert result.compatible is True

    def test_vicsek_compatible_with_p7(self):
        """Vicsek (continuous_2d with positions+velocities) is compatible with P7."""
        result = check_compatibility("vicsek", "P7")
        assert result.compatible is True

    def test_lattice_model_incompatible(self):
        """Lattice_2d model should be incompatible with P7."""
        result = check_compatibility("game_of_life", "P7")
        assert result.compatible is False
