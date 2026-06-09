"""End-to-end tests for P20 quorum sensing / threshold-activated response.

Tests:
- Canonical positive: AutoinducerQuorum reaches DEFINITIVE with hysteresis.
- Negative control: GradedResponseModel (linear, no threshold) NOT detected.
- T1a observation bundle: manual history passes through adapter correctly.
- Determinism: same seed → identical results.
"""

import pytest
import numpy as np

from epc.models.quorum_sensing import (
    AutoinducerQuorum,
    AutoinducerParams,
    GradedResponseModel,
)
from epc.detectors.p20_quorum_sensing import (
    P20QuorumSensingDetector,
    extract_observation_bundle,
)
from epc.detector_result import DetectionTier


# ── Fixtures ──────────────────────────────────────────────────────────

@pytest.fixture(scope="module")
def canonical_result():
    """Run canonical positive and return (result, metadata)."""
    model = AutoinducerQuorum(AutoinducerParams(seed=42))
    history = model.simulate()
    metadata = model.get_metadata()
    detector = P20QuorumSensingDetector(n_permutations=199, seed=42)
    result = detector.detect(history, metadata)
    return result, metadata


# ── Canonical positive ────────────────────────────────────────────────

class TestP20CanonicalPositive:

    def test_detected(self, canonical_result):
        result, _ = canonical_result
        assert result.detected

    def test_definitive_tier(self, canonical_result):
        result, _ = canonical_result
        assert result.tier == DetectionTier.DEFINITIVE

    def test_confidence(self, canonical_result):
        result, _ = canonical_result
        assert result.confidence >= 0.75

    def test_step_r2_high(self, canonical_result):
        result, _ = canonical_result
        assert result.primary_metric['step_r2'] > 0.9

    def test_hysteresis_present(self, canonical_result):
        result, _ = canonical_result
        assert result.secondary_metrics['has_hysteresis']

    def test_hysteresis_width_positive(self, canonical_result):
        result, _ = canonical_result
        assert result.primary_metric['hysteresis_width'] > 0.1

    def test_activation_gt_deactivation(self, canonical_result):
        result, _ = canonical_result
        assert (result.primary_metric['activation_density']
                > result.primary_metric['deactivation_density'])

    def test_null_rejected(self, canonical_result):
        result, _ = canonical_result
        assert result.null_p_value <= 0.005

    def test_effect_size(self, canonical_result):
        result, _ = canonical_result
        assert result.effect_size['cohens_d'] > 1.0


# ── Negative control ─────────────────────────────────────────────────

class TestP20NegativeControl:

    def test_graded_not_detected(self):
        """Graded (linear) response should NOT be detected as quorum."""
        model = GradedResponseModel(seed=42)
        history = model.simulate()
        metadata = model.get_metadata()
        detector = P20QuorumSensingDetector(n_permutations=99, seed=42)
        result = detector.detect(history, metadata)
        assert not result.detected

    def test_graded_metadata_no_threshold(self):
        """Graded model metadata confirms no threshold activation."""
        model = GradedResponseModel(seed=42)
        metadata = model.get_metadata()
        assert metadata['has_threshold_activation'] is False
        assert metadata['has_hysteresis'] is False


# ── T1a observation bundle ───────────────────────────────────────────

class TestT1aObservationBundle:

    def test_bundle_keys(self):
        """Bundle has required keys."""
        history = [
            {'density': 0.5, 'concentration': 0.2, 'collective_state': 0,
             'fraction_on': 0.0, 'step': 0, 'density_idx': 0,
             'sweep_direction': 'up'},
            {'density': 1.0, 'concentration': 1.5, 'collective_state': 1,
             'fraction_on': 0.8, 'step': 0, 'density_idx': 1,
             'sweep_direction': 'up'},
        ]
        bundle = extract_observation_bundle(history)
        for key in ['density', 'concentration', 'collective_state',
                     'fraction_on', 'density_idx', 'sweep_direction']:
            assert key in bundle

    def test_bundle_values_match(self):
        """Bundle values match input history."""
        history = [
            {'density': 0.5, 'concentration': 0.2, 'collective_state': 0,
             'fraction_on': 0.1, 'step': 0, 'density_idx': 0,
             'sweep_direction': 'up'},
        ]
        bundle = extract_observation_bundle(history)
        assert bundle['density'][0] == pytest.approx(0.5)
        assert bundle['concentration'][0] == pytest.approx(0.2)
        assert bundle['collective_state'][0] == 0
        assert bundle['fraction_on'][0] == pytest.approx(0.1)

    def test_detector_uses_bundle(self):
        """Detector operates on manual history (not model object)."""
        # Construct a sharp step in manual history
        history = []
        for d_idx in range(20):
            density = 0.1 + d_idx * 0.15
            frac = 0.0 if density < 1.5 else 1.0
            for step in range(50):
                history.append({
                    'density': density,
                    'concentration': density * frac,
                    'collective_state': int(frac > 0.5),
                    'fraction_on': frac,
                    'step': step,
                    'density_idx': d_idx,
                    'sweep_direction': 'up',
                })
        # Add down-sweep with hysteresis
        for d_idx in range(20):
            density = 3.0 - d_idx * 0.15
            frac = 0.0 if density < 0.5 else 1.0
            for step in range(50):
                history.append({
                    'density': density,
                    'concentration': density * frac,
                    'collective_state': int(frac > 0.5),
                    'fraction_on': frac,
                    'step': step,
                    'density_idx': 20 + d_idx,
                    'sweep_direction': 'down',
                })

        detector = P20QuorumSensingDetector(n_permutations=99, seed=42)
        result = detector.detect(history)
        assert result.detected
        assert result.tier >= DetectionTier.CONFIRMATION


# ── Determinism ──────────────────────────────────────────────────────

class TestP20Determinism:

    def test_same_seed_same_result(self):
        """Same seed produces identical results."""
        model1 = AutoinducerQuorum(AutoinducerParams(seed=42))
        model2 = AutoinducerQuorum(AutoinducerParams(seed=42))
        h1 = model1.simulate()
        h2 = model2.simulate()
        assert len(h1) == len(h2)
        for i in range(0, len(h1), 100):
            assert h1[i]['fraction_on'] == h2[i]['fraction_on']
            assert h1[i]['concentration'] == h2[i]['concentration']
