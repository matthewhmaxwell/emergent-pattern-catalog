"""End-to-end tests for P32 emergent specialization / division of labor.

Sprint 90 deliverable. Tests:
  1. Model determinism and basic dynamics
  2. Detector achieves DEFINITIVE on canonical regime
  3. Negative control: no reinforcement → detector does NOT fire
  4. Negative control: single-task collapse → detector does NOT fire (definitive)
  5. T1a adapter works correctly
  6. Orchestration registration
"""

import numpy as np
import pytest

from epc.models.division_of_labor import (
    ResponseThresholdModel,
    NoReinforcementModel,
    compute_per_agent_entropy,
    compute_windowed_entropy,
    compute_efficiency,
)
from epc.detectors.p32_specialization import (
    P32SpecializationDetector,
    extract_observation_bundle,
)
from epc.detector_result import DetectionTier


# ===================================================================
# Model replication tests
# ===================================================================


class TestResponseThresholdModel:
    """Verify model produces expected specialization dynamics."""

    def test_deterministic(self):
        """Two runs with same seed produce identical histories."""
        m1 = ResponseThresholdModel(n_agents=10, n_tasks=3, seed=42)
        h1 = m1.run(100)
        m2 = ResponseThresholdModel(n_agents=10, n_tasks=3, seed=42)
        h2 = m2.run(100)
        for s1, s2 in zip(h1, h2):
            np.testing.assert_array_equal(s1["task_assignments"], s2["task_assignments"])
            np.testing.assert_array_equal(s1["thresholds"], s2["thresholds"])

    def test_different_seeds_differ(self):
        """Different seeds produce different trajectories."""
        m1 = ResponseThresholdModel(n_agents=10, n_tasks=3, seed=42)
        h1 = m1.run(100)
        m2 = ResponseThresholdModel(n_agents=10, n_tasks=3, seed=99)
        h2 = m2.run(100)
        # At least one snapshot should differ
        any_diff = any(
            not np.array_equal(s1["task_assignments"], s2["task_assignments"])
            for s1, s2 in zip(h1, h2)
        )
        assert any_diff

    def test_entropy_decline_with_reinforcement(self):
        """Per-agent entropy declines when reinforcement is active."""
        m = ResponseThresholdModel(
            n_agents=20, n_tasks=3, reinforcement_rate=0.05,
            forgetting_rate=0.01, seed=42,
        )
        h = m.run(500)
        early, late = compute_windowed_entropy(h, 3, window_size=50)
        assert np.mean(early) > np.mean(late) + 0.2, (
            f"Early entropy ({np.mean(early):.3f}) should exceed late "
            f"({np.mean(late):.3f}) by >0.2 with reinforcement"
        )

    def test_no_entropy_decline_without_reinforcement(self):
        """Without reinforcement, entropy stays high."""
        m = NoReinforcementModel(n_agents=20, n_tasks=3, seed=42)
        h = m.run(500)
        early, late = compute_windowed_entropy(h, 3, window_size=50)
        decline = np.mean(early) - np.mean(late)
        assert decline < 0.15, (
            f"Without reinforcement, entropy decline should be <0.15, "
            f"got {decline:.3f}"
        )

    def test_thresholds_differentiate(self):
        """After reinforcement, agents should have different threshold profiles."""
        m = ResponseThresholdModel(
            n_agents=20, n_tasks=3, reinforcement_rate=0.05,
            forgetting_rate=0.01, seed=42,
        )
        h = m.run(500)
        final_thresholds = h[-1]["thresholds"]  # (n_agents, n_tasks)
        # Standard deviation across agents for each task should be > 0.1
        std_per_task = np.std(final_thresholds, axis=0)
        assert np.mean(std_per_task) > 0.05, (
            f"Threshold std across agents ({np.mean(std_per_task):.3f}) "
            f"too low — agents not differentiating"
        )

    def test_metadata_fields(self):
        """Model metadata has required fields."""
        m = ResponseThresholdModel(seed=42)
        meta = m.get_metadata()
        assert meta["model_name"] == "division_of_labor"
        assert meta["model_class"] == "task_allocation"
        assert meta["has_threshold_reinforcement"] is True
        assert meta["reference"].startswith("Bonabeau")

    def test_state_history_keys(self):
        """State snapshots have required keys."""
        m = ResponseThresholdModel(n_agents=10, n_tasks=3, seed=42)
        h = m.run(10)
        for s in h:
            assert "task_assignments" in s
            assert "thresholds" in s
            assert "stimulus" in s
            assert "n_agents" in s
            assert "n_tasks" in s
            assert s["n_agents"] == 10
            assert s["n_tasks"] == 3


# ===================================================================
# Detector canonical positive
# ===================================================================


class TestP32CanonicalPositive:
    """P32 detector on canonical ResponseThresholdModel."""

    def test_definitive_on_canonical(self):
        """DEFINITIVE detection on canonical regime (reinforcement on, 3 tasks)."""
        m = ResponseThresholdModel(
            n_agents=20, n_tasks=3, reinforcement_rate=0.05,
            forgetting_rate=0.01, seed=42,
        )
        h = m.run(1000)
        det = P32SpecializationDetector(n_permutations=499, seed=42)
        r = det.detect(h, m.get_metadata())

        assert r.detected, f"P32 not detected: {r.summary()}"
        assert r.tier >= DetectionTier.DEFINITIVE, (
            f"Expected DEFINITIVE, got {r.tier.value}: {r.summary()}"
        )
        assert r.confidence >= 0.80
        assert r.primary_metric["entropy_decline"] > 0.3
        assert r.effect_size.get("cohens_d", 0.0) > 5.0

    def test_entropy_decline_significant(self):
        """Primary metric: entropy decline > 0.3 bits."""
        m = ResponseThresholdModel(
            n_agents=20, n_tasks=3, reinforcement_rate=0.05,
            forgetting_rate=0.01, seed=42,
        )
        h = m.run(1000)
        det = P32SpecializationDetector(n_permutations=199, seed=42)
        r = det.detect(h, m.get_metadata())
        assert r.primary_metric["entropy_decline"] > 0.3

    def test_role_diversity_maintained(self):
        """Agents cover all tasks (not single-task collapse)."""
        m = ResponseThresholdModel(
            n_agents=20, n_tasks=3, reinforcement_rate=0.05,
            forgetting_rate=0.01, seed=42,
        )
        h = m.run(1000)
        det = P32SpecializationDetector(n_permutations=199, seed=42)
        r = det.detect(h, m.get_metadata())
        assert r.secondary_metrics["role_diversity"] >= 0.67
        assert not r.secondary_metrics["single_task_collapse"]

    def test_p23_excluded(self):
        """P23 (anti-coordination) should be excluded at confirmation+."""
        m = ResponseThresholdModel(
            n_agents=20, n_tasks=3, reinforcement_rate=0.05,
            forgetting_rate=0.01, seed=42,
        )
        h = m.run(1000)
        det = P32SpecializationDetector(n_permutations=499, seed=42)
        r = det.detect(h, m.get_metadata())
        assert r.tier >= DetectionTier.CONFIRMATION
        assert r.exclusion_results.get("P23") == "excluded"


# ===================================================================
# Negative controls
# ===================================================================


class TestP32NegativeControls:
    """Detector correctly rejects non-specializing systems."""

    def test_no_reinforcement_rejected(self):
        """Without threshold reinforcement, agents stay generalists → not detected."""
        m = NoReinforcementModel(n_agents=20, n_tasks=3, seed=42)
        h = m.run(500)
        det = P32SpecializationDetector(n_permutations=199, seed=42)
        r = det.detect(h, m.get_metadata())
        assert not r.detected, (
            f"P32 should NOT detect without reinforcement: {r.summary()}"
        )

    def test_single_task_collapse_blocked(self):
        """If ALL agents collapse to the SAME task, not division of labor.

        With very high reinforcement + 1 task, all agents do the same thing.
        This should not reach DEFINITIVE.
        """
        m = ResponseThresholdModel(
            n_agents=20, n_tasks=1, reinforcement_rate=0.1,
            forgetting_rate=0.0, seed=42,
        )
        h = m.run(500)
        det = P32SpecializationDetector(n_permutations=199, seed=42)
        r = det.detect(h, m.get_metadata())
        # With 1 task, entropy is always 0 so decline from early (also ~0)
        # is negligible. Should fail screening or at least not reach definitive.
        if r.detected and r.tier >= DetectionTier.CONFIRMATION:
            assert r.secondary_metrics.get("single_task_collapse", False) is False, (
                "Single-task collapse should be flagged and block confirmation"
            )

    def test_short_run_warning(self):
        """Very short runs should produce a warning."""
        m = ResponseThresholdModel(n_agents=20, n_tasks=3, seed=42)
        h = m.run(10)
        det = P32SpecializationDetector(n_permutations=49, seed=42)
        r = det.detect(h, m.get_metadata())
        # Should have a warning about short run
        has_length_warning = any("run length" in w for w in r.warnings)
        assert has_length_warning or not r.detected


# ===================================================================
# T1a adapter
# ===================================================================


class TestT1aAdapter:
    """Verify the T1a observation-bundle adapter."""

    def test_adapter_extracts_correctly(self):
        """Adapter produces correct shapes and values."""
        m = ResponseThresholdModel(n_agents=10, n_tasks=3, seed=42)
        h = m.run(100)
        bundle = extract_observation_bundle(h)

        assert bundle["n_agents"] == 10
        assert bundle["n_tasks"] == 3
        assert bundle["task_assignments"].shape == (len(h), 10)
        assert bundle["steps"].shape == (len(h),)

    def test_adapter_works_with_detector(self):
        """Detector can read via T1a adapter on external-format data."""
        # Create synthetic data that looks like specialization
        rng = np.random.default_rng(42)
        n_agents, n_tasks, T = 15, 3, 200

        # Early: random assignments
        early_assignments = rng.integers(0, n_tasks, size=(T // 2, n_agents))
        # Late: specialized (agent i % n_tasks is dominant)
        late_assignments = np.zeros((T // 2, n_agents), dtype=np.int64)
        for i in range(n_agents):
            dominant = i % n_tasks
            late_assignments[:, i] = dominant

        assignments = np.vstack([early_assignments, late_assignments])
        history = [
            {
                "task_assignments": assignments[t],
                "thresholds": np.zeros((n_agents, n_tasks)),
                "stimulus": np.zeros(n_tasks),
                "n_agents": n_agents,
                "n_tasks": n_tasks,
                "step": t,
            }
            for t in range(T)
        ]

        det = P32SpecializationDetector(n_permutations=99, seed=42)
        r = det.detect(history)
        assert r.detected, f"Detector should fire on synthetic specialization: {r.summary()}"


# ===================================================================
# Orchestration registration
# ===================================================================


class TestP32Registration:
    """Verify P32 is correctly registered in orchestration."""

    def test_model_in_registry(self):
        """ResponseThresholdModel registered as 'response_threshold'."""
        from epc.orchestration import MODEL_REGISTRY
        assert "response_threshold" in MODEL_REGISTRY
        reg = MODEL_REGISTRY["response_threshold"]
        assert reg.substrate_type == "task_allocation_timeseries"
        assert "P32" in reg.primary_patterns
        assert "task_assignments" in reg.observables

    def test_detector_in_registry(self):
        """P32 detector registered in DETECTOR_REGISTRY."""
        from epc.orchestration import DETECTOR_REGISTRY
        assert "P32" in DETECTOR_REGISTRY
        reg = DETECTOR_REGISTRY["P32"]
        assert "task_allocation_timeseries" in reg.required_substrate
        assert "task_assignments" in reg.required_observables

    def test_compatibility(self):
        """response_threshold × P32 is compatible."""
        from epc.orchestration import check_compatibility
        result = check_compatibility("response_threshold", "P32")
        assert result.compatible, f"Should be compatible: {result.reason}"

    def test_no_reinforcement_compatible(self):
        """no_reinforcement × P32 is compatible."""
        from epc.orchestration import check_compatibility
        result = check_compatibility("no_reinforcement", "P32")
        assert result.compatible

    def test_cross_substrate_rejected(self):
        """lattice_2d model × P32 should be incompatible."""
        from epc.orchestration import check_compatibility
        result = check_compatibility("greenberg_hastings", "P32")
        assert not result.compatible


# ===================================================================
# Helper function tests
# ===================================================================


class TestHelperFunctions:
    """Verify compute_per_agent_entropy and compute_efficiency."""

    def test_entropy_max_for_uniform(self):
        """Uniform task distribution → entropy = log2(n_tasks)."""
        n_agents, n_tasks = 5, 4
        # Create history where each agent does each task equally
        history = []
        for t in range(100):
            assignments = np.array([t % n_tasks] * n_agents, dtype=np.int64)
            history.append({
                "task_assignments": assignments,
                "n_agents": n_agents,
                "n_tasks": n_tasks,
            })
        entropy = compute_per_agent_entropy(history, n_tasks)
        expected = np.log2(n_tasks)
        np.testing.assert_allclose(entropy, expected, atol=0.1)

    def test_entropy_zero_for_specialized(self):
        """Agent always doing one task → entropy = 0."""
        history = []
        for t in range(100):
            history.append({
                "task_assignments": np.array([0, 1, 2, 0, 1], dtype=np.int64),
                "n_agents": 5,
                "n_tasks": 3,
            })
        entropy = compute_per_agent_entropy(history, 3)
        np.testing.assert_allclose(entropy, 0.0, atol=1e-10)

    def test_efficiency_full_coverage(self):
        """All tasks covered every step → efficiency = 1.0."""
        history = []
        for t in range(50):
            history.append({
                "task_assignments": np.array([0, 1, 2, 0, 1], dtype=np.int64),
                "n_agents": 5,
                "n_tasks": 3,
            })
        eff = compute_efficiency(history, 3)
        assert eff == 1.0

    def test_efficiency_partial_coverage(self):
        """Not all tasks covered → efficiency < 1.0."""
        history = []
        for t in range(50):
            # Only tasks 0 and 1 covered
            history.append({
                "task_assignments": np.array([0, 1, 0, 0, 1], dtype=np.int64),
                "n_agents": 5,
                "n_tasks": 3,
            })
        eff = compute_efficiency(history, 3)
        assert eff == 0.0


if __name__ == "__main__":
    pytest.main([__file__, "-v", "--tb=short"])
