"""
End-to-end tests for P30 autopoiesis (spontaneous boundary formation).

Tests:
- Canonical model produces detectable autopoiesis at CONFIRMATION tier
- NonBondingParticleModel negative control: no closure, no detection
- DenseClusterModel negative control: clustering but no closed membrane
- Detector rejects non-autopoietic histories
- Orchestration registration
- T1a observation bundle adapter
"""

import numpy as np
import pytest

from epc.models.autopoiesis import (
    AutopoiesisModel,
    DenseClusterModel,
    NonBondingParticleModel,
)
from epc.detectors.p30_autopoiesis import (
    P30AutopoiesisDetector,
    extract_observation_bundle,
)
from epc.orchestration import (
    DETECTOR_REGISTRY,
    MODEL_REGISTRY,
    check_compatibility,
)


# === Canonical positive detection ===

class TestP30CanonicalDetection:
    """P30 detected on canonical autopoiesis model."""

    @pytest.fixture(scope="class")
    def canonical_result(self):
        m = AutopoiesisModel(
            n_substrate=100, n_catalyst=3, box_size=20.0,
            production_rate=0.15, decay_rate=0.01, link_attraction=0.3,
            membrane_equilibrium_radius=3.0, catalyst_link_attraction=0.5,
            production_radius=3.0, seed=42,
        )
        h = m.run(150, record_every=5)
        d = P30AutopoiesisDetector(n_permutations=199, association_radius=4.0, seed=42)
        return d.detect(h, m.get_metadata())

    def test_detected(self, canonical_result):
        assert canonical_result.detected, (
            f"P30 not detected; assoc={canonical_result.primary_metric.get('association_score', 0):.3f}"
        )

    def test_at_least_confirmation(self, canonical_result):
        from epc.detector_result import DetectionTier
        assert canonical_result.tier >= DetectionTier.CONFIRMATION, (
            f"Only {canonical_result.tier.value}, p={canonical_result.null_p_value:.4f}"
        )

    def test_closure_fraction(self, canonical_result):
        assert canonical_result.primary_metric['closure_fraction'] > 0.5, (
            f"Low closure: {canonical_result.primary_metric['closure_fraction']:.3f}"
        )

    def test_association_above_csr(self, canonical_result):
        assert canonical_result.primary_metric['association_score'] > 1.5, (
            f"Association {canonical_result.primary_metric['association_score']:.3f} not above 1.5"
        )

    def test_enrichment_above_1(self, canonical_result):
        assert canonical_result.primary_metric['enrichment_ratio'] > 1.0, (
            f"Enrichment {canonical_result.primary_metric['enrichment_ratio']:.3f} not above 1.0"
        )

    def test_null_p_below_001(self, canonical_result):
        assert canonical_result.null_p_value < 0.01, (
            f"p={canonical_result.null_p_value:.4f} not below 0.01"
        )

    def test_confidence_above_055(self, canonical_result):
        assert canonical_result.confidence >= 0.55, (
            f"Confidence {canonical_result.confidence:.3f} too low"
        )

    def test_p1_excluded(self, canonical_result):
        assert canonical_result.exclusion_results.get('P1') == 'excluded', (
            f"P1 not excluded: {canonical_result.exclusion_results}"
        )


# === Negative controls ===

class TestP30NegativeControls:
    """P30 correctly rejects non-autopoietic models."""

    def test_non_bonding_rejected(self):
        """Non-bonding particles: no types, no closure → screening failure."""
        m = NonBondingParticleModel(n_particles=100, box_size=20.0, seed=42)
        h = m.run(100, record_every=5)
        d = P30AutopoiesisDetector(n_permutations=49, association_radius=4.0, seed=42)
        r = d.detect(h, m.get_metadata())
        assert not r.detected, (
            f"P30 should not fire on non-bonding particles, got assoc={r.primary_metric['association_score']:.3f}"
        )

    def test_dense_cluster_rejected(self):
        """Dense cluster (P1-like): aggregation but no closed membrane."""
        m = DenseClusterModel(n_particles=50, box_size=20.0, seed=42)
        h = m.run(80, record_every=5)
        d = P30AutopoiesisDetector(n_permutations=49, association_radius=4.0, seed=42)
        r = d.detect(h, m.get_metadata())
        assert not r.detected, (
            f"P30 should not fire on dense cluster, got assoc={r.primary_metric['association_score']:.3f}"
        )


# === Orchestration registration ===

class TestP30Registration:
    """P30 correctly registered in orchestration."""

    def test_model_registered(self):
        assert 'autopoiesis' in MODEL_REGISTRY

    def test_detector_registered(self):
        assert 'P30' in DETECTOR_REGISTRY

    def test_compatible(self):
        result = check_compatibility('autopoiesis', 'P30')
        assert result.compatible, f"Not compatible: {result.reason}"

    def test_substrate_type(self):
        assert MODEL_REGISTRY['autopoiesis'].substrate_type == 'particle_membrane'

    def test_primary_pattern(self):
        assert 'P30' in MODEL_REGISTRY['autopoiesis'].primary_patterns


# === T1a observation bundle ===

class TestP30ObservationBundle:
    """T1a observation bundle adapter works correctly."""

    def test_bundle_keys(self):
        m = AutopoiesisModel(n_substrate=30, n_catalyst=2, box_size=10.0, seed=42)
        h = m.run(10, record_every=5)
        bundle = extract_observation_bundle(h)
        assert 'positions' in bundle
        assert 'types' in bundle
        assert 'bonds' in bundle
        assert 'steps' in bundle
        assert 'box_size' in bundle
        assert 'n_particles' in bundle

    def test_bundle_shapes(self):
        m = AutopoiesisModel(n_substrate=30, n_catalyst=2, box_size=10.0, seed=42)
        h = m.run(10, record_every=5)
        bundle = extract_observation_bundle(h)
        assert len(bundle['positions']) == len(h)
        assert bundle['positions'][0].shape == (32, 2)
        assert bundle['types'][0].shape == (32,)


# === Determinism ===

class TestP30Determinism:
    """Model and detector are deterministic from seed."""

    def test_model_deterministic(self):
        h1 = AutopoiesisModel(n_substrate=30, n_catalyst=2, box_size=10.0, seed=42).run(20)
        h2 = AutopoiesisModel(n_substrate=30, n_catalyst=2, box_size=10.0, seed=42).run(20)
        np.testing.assert_array_equal(h1[-1]['positions'], h2[-1]['positions'])
        np.testing.assert_array_equal(h1[-1]['types'], h2[-1]['types'])

    def test_detector_deterministic(self):
        m = AutopoiesisModel(n_substrate=50, n_catalyst=2, box_size=15.0, seed=42)
        h = m.run(50, record_every=5)
        d1 = P30AutopoiesisDetector(n_permutations=49, association_radius=4.0, seed=42)
        d2 = P30AutopoiesisDetector(n_permutations=49, association_radius=4.0, seed=42)
        r1 = d1.detect(h, m.get_metadata())
        r2 = d2.detect(h, m.get_metadata())
        assert r1.confidence == r2.confidence
        assert r1.null_p_value == r2.null_p_value


# === Model state consistency ===

class TestP30ModelState:
    """Model produces valid state throughout simulation."""

    def test_types_valid(self):
        m = AutopoiesisModel(n_substrate=50, n_catalyst=3, box_size=15.0, seed=42)
        h = m.run(30)
        for snap in h:
            types = snap['types']
            assert set(np.unique(types)).issubset({0, 1, 2}), f"Invalid types: {np.unique(types)}"

    def test_catalyst_count_constant(self):
        m = AutopoiesisModel(n_substrate=50, n_catalyst=3, box_size=15.0, seed=42)
        h = m.run(30)
        for snap in h:
            n_cat = np.sum(snap['types'] == 1)
            assert n_cat == 3, f"Catalyst count changed: {n_cat}"

    def test_total_particle_count_constant(self):
        m = AutopoiesisModel(n_substrate=50, n_catalyst=3, box_size=15.0, seed=42)
        h = m.run(30)
        for snap in h:
            assert len(snap['types']) == 53

    def test_positions_in_bounds(self):
        m = AutopoiesisModel(n_substrate=50, n_catalyst=3, box_size=15.0, seed=42)
        h = m.run(30)
        for snap in h:
            pos = snap['positions']
            assert np.all(pos >= 0) and np.all(pos < 15.0), "Positions out of bounds"

    def test_links_produced(self):
        m = AutopoiesisModel(n_substrate=50, n_catalyst=3, box_size=15.0,
                             production_rate=0.15, seed=42)
        h = m.run(50)
        final_links = np.sum(h[-1]['types'] == 2)
        assert final_links > 0, "No links produced"

    def test_metadata(self):
        m = AutopoiesisModel(seed=42)
        meta = m.get_metadata()
        assert meta['model_class'] == 'autopoiesis'
        assert meta['has_production'] is True
        assert meta['has_decay'] is True
        assert meta['has_bonding'] is True
