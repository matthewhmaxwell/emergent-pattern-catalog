"""
End-to-end tests for P19 (emergent leadership / minority guidance).

Tests:
1. Canonical positive: ρ=0.1, moderate bias → DEFINITIVE detection
2. Negative control: ρ=0.0 → no asymmetric influence → NOT detected
3. Model produces deterministic output from seed
4. Detector uses existing TE infrastructure (KSG estimator)
5. Transfer matrix: P19 is registered, +1 model, +1 detector

Sprint 69 deliverable.
"""

from __future__ import annotations

import numpy as np
import pytest

from epc.models.informed_minority import (
    InformedMinorityModel,
    group_accuracy,
    polarization,
)
from epc.detectors.p19_emergent_leadership import detect_p19
from epc.orchestration import (
    MODEL_REGISTRY,
    DETECTOR_REGISTRY,
    check_compatibility,
)


# === Model tests ===

class TestInformedMinorityModel:
    """Basic model validation."""

    def test_deterministic_from_seed(self):
        """Two runs with same seed produce identical histories."""
        kwargs = dict(
            n_particles=50, box_size=7.0, speed=0.03, noise=0.1,
            informed_fraction=0.1, bias_weight=0.3, seed=42,
        )
        h1 = InformedMinorityModel(**kwargs).run(50)
        h2 = InformedMinorityModel(**kwargs).run(50)
        np.testing.assert_array_equal(h1[-1]["headings"], h2[-1]["headings"])
        np.testing.assert_array_equal(h1[-1]["positions"], h2[-1]["positions"])

    def test_informed_mask_correct_count(self):
        """Informed mask has the right number of True entries."""
        m = InformedMinorityModel(n_particles=100, informed_fraction=0.1, seed=42)
        s = m.setup()
        assert np.sum(s["informed_mask"]) == 10

    def test_zero_informed_fraction(self):
        """ρ=0 produces no informed agents."""
        m = InformedMinorityModel(n_particles=100, informed_fraction=0.0, seed=42)
        s = m.setup()
        assert np.sum(s["informed_mask"]) == 0

    def test_state_keys(self):
        """State dict has all expected keys."""
        m = InformedMinorityModel(n_particles=50, seed=42)
        s = m.setup()
        for key in ["positions", "velocities", "headings", "informed_mask", "step"]:
            assert key in s, f"Missing key: {key}"

    def test_metadata_keys(self):
        """Metadata has required keys."""
        m = InformedMinorityModel(n_particles=50, seed=42)
        md = m.get_metadata()
        for key in ["model_family", "informed_fraction", "bias_weight",
                     "preferred_direction", "has_informed_minority"]:
            assert key in md, f"Missing metadata key: {key}"

    def test_polarization_with_informed(self):
        """With informed agents, group polarizes after sufficient steps."""
        m = InformedMinorityModel(
            n_particles=100, box_size=7.0, speed=0.03, noise=0.1,
            informed_fraction=0.1, bias_weight=0.3, seed=42,
        )
        h = m.run(300)
        phi = polarization(h[-1])
        assert phi > 0.5, f"Expected polarization > 0.5, got {phi:.3f}"

    def test_accuracy_with_informed(self):
        """With informed agents, group aligns with preferred direction."""
        m = InformedMinorityModel(
            n_particles=100, box_size=7.0, speed=0.03, noise=0.1,
            informed_fraction=0.1, bias_weight=0.3,
            preferred_direction=0.0, seed=42,
        )
        h = m.run(300)
        acc = group_accuracy(h[-1], 0.0)
        assert acc > 0.5, f"Expected accuracy > 0.5, got {acc:.3f}"


# === Detector tests ===

class TestP19DetectorCanonical:
    """Canonical positive: ρ=0.1 informed minority steers group."""

    @pytest.fixture(scope="class")
    def canonical_result(self):
        """Run detector on canonical regime."""
        m = InformedMinorityModel(
            n_particles=100, box_size=7.0, speed=0.03, noise=0.1,
            informed_fraction=0.1, bias_weight=0.3,
            preferred_direction=0.0, seed=42,
        )
        history = m.run(500)
        metadata = m.get_metadata()
        # 199 permutations for floor p=0.005, needed for definitive tier
        return detect_p19(history, metadata, n_permutations=199, seed=42)

    def test_detected(self, canonical_result):
        assert canonical_result.detected is True

    def test_pattern_id(self, canonical_result):
        assert canonical_result.pattern_id == "P19"

    def test_tier_definitive(self, canonical_result):
        assert canonical_result.tier == "definitive", (
            f"Expected definitive, got {canonical_result.tier}"
        )

    def test_accuracy_positive(self, canonical_result):
        acc = canonical_result.primary_metric["group_accuracy_mean"]
        assert acc > 0.3, f"Expected accuracy > 0.3, got {acc:.3f}"

    def test_pull_positive(self, canonical_result):
        """Influence pull from informed agents should be positive."""
        pull = canonical_result.secondary_metrics.get("pull_mean", 0.0)
        assert pull > 0.0, f"Expected positive pull, got {pull:.4f}"


class TestP19NegativeControl:
    """Negative control: ρ=0 → no directional guidance → not detected."""

    @pytest.fixture(scope="class")
    def negative_result(self):
        """Run detector with no informed agents."""
        m = InformedMinorityModel(
            n_particles=100, box_size=7.0, speed=0.03, noise=0.1,
            informed_fraction=0.0, bias_weight=0.3,
            preferred_direction=0.0, seed=42,
        )
        history = m.run(500)
        metadata = m.get_metadata()
        # ρ=0 means informed_mask is all False
        # The group will polarize (Vicsek alignment) but NOT toward
        # the preferred direction reliably
        return detect_p19(history, metadata, n_permutations=99, seed=42)

    def test_not_detected_or_low_tier(self, negative_result):
        """ρ=0 should not produce definitive detection."""
        # With ρ=0, the group may randomly align with preferred direction
        # in some seeds, but the detector should not reach definitive
        # because there's no informed mask to compute TE from
        assert negative_result.tier in ("none", "screening"), (
            f"Expected none or screening for ρ=0, got {negative_result.tier}"
        )


# === Registry tests ===

class TestP19Registry:
    """P19 is properly registered in orchestration."""

    def test_informed_minority_in_model_registry(self):
        assert "informed_minority" in MODEL_REGISTRY

    def test_p19_in_detector_registry(self):
        assert "P19" in DETECTOR_REGISTRY

    def test_informed_minority_p19_compatible(self):
        result = check_compatibility("informed_minority", "P19")
        assert result.compatible, f"Expected compatible, got: {result.reason}"

    def test_informed_minority_substrate(self):
        assert MODEL_REGISTRY["informed_minority"].substrate_type == "continuous_2d"

    def test_p19_required_substrate(self):
        assert "continuous_2d" in DETECTOR_REGISTRY["P19"].required_substrate

    def test_vicsek_p19_compatible(self):
        """Vicsek should be substrate-compatible with P19 (both continuous_2d)."""
        result = check_compatibility("vicsek", "P19")
        assert result.compatible
