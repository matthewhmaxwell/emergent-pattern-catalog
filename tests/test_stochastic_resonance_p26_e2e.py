"""End-to-end tests for P26 stochastic resonance.

Tests cover:
- Deterministic bistable double-well → DEFINITIVE detection
- Suprathreshold signal negative control → NOT detected
- Observation bundle adapter (T1a contract)
- Threshold unit cross-model test (T1b)
"""

import numpy as np
import pytest

from epc.models.stochastic_resonance import (
    BistableDoubleWell,
    ThresholdUnit,
    DoubleWellParams,
    ThresholdUnitParams,
)
from epc.detectors.p26_stochastic_resonance import (
    P26StochasticResonanceDetector,
    extract_observation_bundle,
)


# ── Canonical positive: bistable double-well ──

class TestP26DoubleWellPositive:
    """P26 on canonical bistable double-well with subthreshold signal."""

    @pytest.fixture(scope='class')
    def sr_result(self):
        params = DoubleWellParams(seed=42, n_trials=20, n_steps=20000)
        model = BistableDoubleWell(params)
        history = model.simulate()
        det = P26StochasticResonanceDetector(n_permutations=199, seed=42)
        return det.detect(history, model.get_metadata())

    def test_detected(self, sr_result):
        assert sr_result.detected, f"P26 not detected: {sr_result.notes}"

    def test_definitive_tier(self, sr_result):
        assert sr_result.tier.value == "definitive", (
            f"Expected definitive, got {sr_result.tier.value}"
        )

    def test_interior_peak(self, sr_result):
        assert sr_result.secondary_metrics['is_interior_peak'], (
            "Peak must be at interior noise level (not endpoints)"
        )

    def test_inverted_u_shape(self, sr_result):
        assert sr_result.secondary_metrics['has_rise'], "Missing rise"
        assert sr_result.secondary_metrics['has_fall'], "Missing fall"

    def test_gain_over_zero(self, sr_result):
        gain = sr_result.primary_metric['gain_over_zero']
        assert gain > 0.05, f"Gain {gain} not > 0.05"

    def test_deterministic_from_seed(self):
        """Same seed → same result."""
        params = DoubleWellParams(seed=42, n_trials=10, n_steps=10000)
        confidences = []
        for _ in range(3):
            model = BistableDoubleWell(params)
            history = model.simulate()
            det = P26StochasticResonanceDetector(n_permutations=99, seed=42)
            result = det.detect(history, model.get_metadata())
            confidences.append(result.confidence)
        assert confidences[0] == confidences[1] == confidences[2]


# ── Negative control: suprathreshold signal ──

class TestP26NegativeControl:
    """Suprathreshold signal → no inverted-U → NOT detected."""

    def test_suprathreshold_not_detected(self):
        """Signal above threshold → detected at all noise → monotone, no SR."""
        params = ThresholdUnitParams(
            threshold=1.0, signal_amplitude=1.5, seed=42,
        )
        model = ThresholdUnit(params)
        history = model.simulate()
        det = P26StochasticResonanceDetector(n_permutations=99, seed=42)
        result = det.detect(history, model.get_metadata())
        assert not result.detected, (
            f"P26 false positive on suprathreshold signal"
        )

    def test_zero_noise_only_not_detected(self):
        """Only zero noise → no sweep → not detected."""
        params = ThresholdUnitParams(
            seed=42, noise_levels=[0.0],
        )
        model = ThresholdUnit(params)
        history = model.simulate()
        det = P26StochasticResonanceDetector(n_permutations=99, seed=42)
        result = det.detect(history, model.get_metadata())
        assert not result.detected


# ── T1a observation bundle ──

class TestT1aObservationBundle:
    """Observation bundle adapter contract."""

    def test_bundle_keys(self):
        params = DoubleWellParams(
            seed=0, noise_levels=[0.0, 0.5], n_steps=100, n_trials=1,
        )
        model = BistableDoubleWell(params)
        history = model.simulate()
        bundle = extract_observation_bundle(history)
        expected_keys = {'time', 'x', 'signal', 'noise_level', 'noise_level_idx'}
        assert set(bundle.keys()) == expected_keys
        for k, v in bundle.items():
            assert isinstance(v, np.ndarray), f"{k} is not ndarray"
            assert v.shape == (200,), f"{k} shape is {v.shape}"

    def test_bundle_values_match_history(self):
        params = DoubleWellParams(
            seed=0, noise_levels=[0.0, 0.5], n_steps=50, n_trials=1,
        )
        model = BistableDoubleWell(params)
        history = model.simulate()
        bundle = extract_observation_bundle(history)
        for i in [0, 10, 50, 99]:
            assert bundle['x'][i] == pytest.approx(history[i]['x'])
            assert bundle['signal'][i] == pytest.approx(history[i]['signal'])
            assert bundle['noise_level'][i] == pytest.approx(history[i]['noise_level'])

    def test_detector_uses_bundle(self):
        """Detector reads through bundle adapter, not model object."""
        # Construct history dicts manually (no model object needed)
        history = []
        # Two noise levels: D=0 (no response) and D=0.5 (some response)
        for nl_idx, D in enumerate([0.0, 0.5, 1.0]):
            for step_i in range(1000):
                t = step_i * 0.01
                signal = 0.5 * np.sin(2.0 * np.pi * 0.05 * t)
                # At D=0.5, x tracks signal; at D=0, x=0; at D=1, x=noise
                if D == 0.0:
                    x = 0.0
                elif D == 0.5:
                    x = np.sign(signal) if abs(signal) > 0.3 else 0.0
                else:
                    x = 0.0
                history.append({
                    'time': t, 'x': float(x), 'signal': float(signal),
                    'noise_level': float(D), 'noise_level_idx': nl_idx,
                    'step': step_i,
                })
        det = P26StochasticResonanceDetector(n_permutations=99, seed=42)
        result = det.detect(history)
        # Should detect — performance peaks at D=0.5
        assert result.detected, f"Not detected on manual bundle"


# ── Threshold unit positive (T1b cross-model) ──

class TestP26ThresholdUnitPositive:
    """Threshold unit (independent SR implementation) also reaches DEFINITIVE."""

    @pytest.fixture(scope='class')
    def threshold_result(self):
        params = ThresholdUnitParams(seed=42)
        model = ThresholdUnit(params)
        history = model.simulate()
        det = P26StochasticResonanceDetector(n_permutations=199, seed=42)
        return det.detect(history, model.get_metadata())

    def test_detected(self, threshold_result):
        assert threshold_result.detected

    def test_definitive_tier(self, threshold_result):
        assert threshold_result.tier.value == "definitive", (
            f"Expected definitive, got {threshold_result.tier.value}"
        )

    def test_interior_peak(self, threshold_result):
        assert threshold_result.secondary_metrics['is_interior_peak']
