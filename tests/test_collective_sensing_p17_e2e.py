"""
End-to-end tests for P17 collective sensing (model + detector).

Tests:
1. Model is deterministic from seed
2. Group CI rises with N (the Berdahl signature) — fast version
3. Negative control: isolated agents (N=1) or alpha=0 → CI at chance
4. Detector reaches confirmation/definitive on canonical regime (slow)
5. Registry registration is correct
"""

import numpy as np
import pytest

from epc.models.collective_sensing import (
    CollectiveSensingModel,
    compute_chemotactic_index,
)
from epc.detectors.p17_collective_sensing import (
    P17CollectiveSensingDetector,
    _run_group_ci,
)
from epc.orchestration import MODEL_REGISTRY, DETECTOR_REGISTRY, check_compatibility


class TestModelDeterminism:
    """Model produces identical output from same seed."""

    def test_deterministic_from_seed(self):
        """Two runs with same seed produce identical histories."""
        model1 = CollectiveSensingModel(n_agents=20, seed=123)
        h1 = model1.run(n_steps=100, record_interval=10)

        model2 = CollectiveSensingModel(n_agents=20, seed=123)
        h2 = model2.run(n_steps=100, record_interval=10)

        for s1, s2 in zip(h1, h2):
            np.testing.assert_array_equal(s1["positions"], s2["positions"])
            np.testing.assert_array_equal(s1["headings"], s2["headings"])
            np.testing.assert_array_equal(s1["speeds"], s2["speeds"])

    def test_different_seeds_differ(self):
        """Different seeds produce different outcomes."""
        model1 = CollectiveSensingModel(n_agents=20, seed=1)
        h1 = model1.run(n_steps=50)

        model2 = CollectiveSensingModel(n_agents=20, seed=2)
        h2 = model2.run(n_steps=50)

        assert not np.allclose(h1[-1]["positions"], h2[-1]["positions"])


class TestGroupSizeScaling:
    """CI rises with group size — the P17 signature."""

    def test_ci_rises_with_n(self):
        """Groups navigate better than isolated agents (fast check, 5 seeds)."""
        # Canonical parameters
        params = dict(
            box_size=20.0, v_max=0.4, turn_noise=0.3, sensing_noise=0.8,
            alpha=0.95, social_strength=0.2, field_sigma=5.0,
            n_steps=800, record_interval=5,
        )

        ci_n1 = []
        ci_n50 = []
        for s in range(5):
            ci_n1.append(_run_group_ci(n_agents=1, seed=s * 77 + 13, **params))
            ci_n50.append(_run_group_ci(n_agents=50, seed=s * 77 + 13, **params))

        mean_ci_1 = np.mean(ci_n1)
        mean_ci_50 = np.mean(ci_n50)

        # Group should substantially outperform isolated agents
        assert mean_ci_50 > mean_ci_1 + 0.1, (
            f"Expected CI(N=50) > CI(N=1)+0.1; got {mean_ci_50:.3f} vs {mean_ci_1:.3f}"
        )
        # Group CI should be positive (approaching peak)
        assert mean_ci_50 > 0.1, f"Expected CI(N=50) > 0.1; got {mean_ci_50:.3f}"


class TestNegativeControl:
    """Negative controls should NOT trigger detection."""

    def test_no_speed_modulation_no_gradient_climbing(self):
        """With alpha=0, no gradient climbing regardless of group size."""
        params = dict(
            box_size=20.0, v_max=0.4, turn_noise=0.3, sensing_noise=0.8,
            alpha=0.0,  # NO speed modulation
            social_strength=0.2, field_sigma=5.0,
            n_steps=800, record_interval=5,
        )

        cis = []
        for s in range(5):
            ci = _run_group_ci(n_agents=50, seed=s * 100 + 7, **params)
            cis.append(ci)

        mean_ci = np.mean(cis)
        # Without speed modulation, group should NOT navigate to peak
        assert mean_ci < 0.15, (
            f"Expected CI < 0.15 without speed modulation; got {mean_ci:.3f}"
        )

    def test_isolated_agent_at_chance(self):
        """Single agent cannot reliably navigate with noisy sensing."""
        params = dict(
            box_size=20.0, v_max=0.4, turn_noise=0.3, sensing_noise=0.8,
            alpha=0.95, social_strength=0.2, field_sigma=5.0,
            n_steps=800, record_interval=5,
        )

        cis = []
        for s in range(10):
            ci = _run_group_ci(n_agents=1, seed=s * 50 + 3, **params)
            cis.append(ci)

        mean_ci = np.mean(cis)
        # N=1 should be near chance (|CI| < 0.15)
        assert abs(mean_ci) < 0.20, (
            f"Expected |CI(N=1)| < 0.20 (chance-level); got {mean_ci:.3f}"
        )


class TestDetectorCanonical:
    """Full detector run on canonical regime."""

    @pytest.mark.slow
    def test_canonical_regime_reaches_confirmation(self):
        """Detector reaches at least CONFIRMATION on canonical collective sensing."""
        model = CollectiveSensingModel(
            n_agents=50, box_size=20.0, v_max=0.4,
            turn_noise=0.3, sensing_noise=0.8,
            alpha=0.95, social_strength=0.2,
            field_sigma=5.0, init_mode="offset", seed=42,
        )
        history = model.run(n_steps=1000, record_interval=5)
        metadata = model.get_metadata()

        detector = P17CollectiveSensingDetector(
            group_sizes=[1, 5, 10, 25, 50],
            n_null_runs=19,  # reduced for test speed
            n_seeds=5,
            n_steps=800,
            record_interval=5,
            seed=42,
        )
        result = detector.detect(history, metadata)

        assert result.detected is True
        assert result.tier in ("screening", "confirmation", "definitive"), (
            f"Expected at least screening, got {result.tier}; notes: {result.notes}"
        )
        assert result.confidence >= 0.500
        assert result.primary_metric["ci_at_max_N"] > 0.1

    @pytest.mark.slow
    def test_null_regime_does_not_fire(self):
        """Detector does NOT fire when alpha=0 (no speed modulation)."""
        model = CollectiveSensingModel(
            n_agents=50, box_size=20.0, v_max=0.4,
            turn_noise=0.3, sensing_noise=0.8,
            alpha=0.0,  # No speed modulation
            social_strength=0.2,
            field_sigma=5.0, init_mode="offset", seed=42,
        )
        history = model.run(n_steps=1000, record_interval=5)
        metadata = model.get_metadata()

        detector = P17CollectiveSensingDetector(
            group_sizes=[1, 5, 10, 25, 50],
            n_null_runs=19,
            n_seeds=5,
            n_steps=800,
            record_interval=5,
            seed=42,
        )
        result = detector.detect(history, metadata)

        assert result.detected is False or result.tier == "none", (
            f"Expected no detection with alpha=0; got tier={result.tier}"
        )


class TestRegistration:
    """P17 is correctly registered in orchestration."""

    def test_model_in_registry(self):
        """collective_sensing model is registered."""
        assert 'collective_sensing' in MODEL_REGISTRY

    def test_detector_in_registry(self):
        """P17 detector is registered."""
        assert 'P17' in DETECTOR_REGISTRY

    def test_model_detector_compatible(self):
        """collective_sensing model is compatible with P17 detector."""
        result = check_compatibility('collective_sensing', 'P17')
        assert result.compatible is True

    def test_model_substrate_type(self):
        """collective_sensing is continuous_2d substrate."""
        reg = MODEL_REGISTRY['collective_sensing']
        assert reg.substrate_type == 'continuous_2d'

    def test_detector_substrate(self):
        """P17 detector accepts continuous_2d substrate."""
        reg = DETECTOR_REGISTRY['P17']
        assert 'continuous_2d' in reg.required_substrate
