"""End-to-end tests for P24 homeostatic regulation.

Tests cover:
- Deterministic proportional homeostat → DEFINITIVE detection
- gain=0 negative control → NOT detected
- Observation bundle adapter (T1a contract)
- Integral controller cross-model test (T1b)
"""

import numpy as np
import pytest

from epc.models.homeostasis import (
    ProportionalHomeostat,
    IntegralHomeostat,
    HomeostatParams,
    PerturbationSchedule,
)
from epc.detectors.p24_homeostasis import (
    P24HomeostasisDetector,
    extract_observation_bundle,
)


# ── Canonical positive: proportional homeostat ──

class TestP24ProportionalPositive:
    """P24 on canonical proportional homeostat with sustained perturbation."""

    @pytest.fixture
    def regulated_history(self):
        params = HomeostatParams(gain=5.0, setpoint=10.0, seed=42, dt=0.1)
        schedule = PerturbationSchedule(onset=50.0, amplitude=5.0)
        model = ProportionalHomeostat(params)
        return model.simulate(2000, schedule=schedule), model.get_metadata()

    def test_detected(self, regulated_history):
        history, metadata = regulated_history
        det = P24HomeostasisDetector(n_permutations=199, seed=42)
        result = det.detect(history, metadata)
        assert result.detected, f"P24 not detected: {result.summary()}"

    def test_definitive_tier(self, regulated_history):
        history, metadata = regulated_history
        det = P24HomeostasisDetector(n_permutations=199, seed=42)
        result = det.detect(history, metadata)
        assert result.tier.value == "definitive", (
            f"Expected definitive, got {result.tier.value}: {result.summary()}"
        )

    def test_deviation_ratio_small(self, regulated_history):
        history, metadata = regulated_history
        det = P24HomeostasisDetector(n_permutations=199, seed=42)
        result = det.detect(history, metadata)
        ratio = result.secondary_metrics.get('deviation_ratio', 1.0)
        assert ratio < 0.3, f"deviation_ratio {ratio} not < 0.3"

    def test_bounded_deviation(self, regulated_history):
        """Regulated system has bounded deviation (growth ratio ≈ 1.0)."""
        history, metadata = regulated_history
        det = P24HomeostasisDetector(n_permutations=199, seed=42)
        result = det.detect(history, metadata)
        growth = result.primary_metric.get('deviation_growth_ratio', float('inf'))
        assert growth <= 2.0, f"growth_ratio {growth} not bounded"

    def test_deterministic_from_seed(self):
        """Same seed → same result."""
        params = HomeostatParams(gain=5.0, setpoint=10.0, seed=42, dt=0.1)
        schedule = PerturbationSchedule(onset=50.0, amplitude=5.0)

        for _ in range(3):
            model = ProportionalHomeostat(params)
            history = model.simulate(2000, schedule=schedule)
            det = P24HomeostasisDetector(n_permutations=99, seed=42)
            result = det.detect(history, model.get_metadata())
            assert result.detected
            assert result.confidence == pytest.approx(result.confidence)


# ── Negative control: gain=0 ──

class TestP24NegativeControl:
    """gain=0 → no feedback → variable drifts → NOT detected."""

    def test_gain_zero_not_detected(self):
        params = HomeostatParams(gain=0.0, setpoint=10.0, seed=42, dt=0.1)
        schedule = PerturbationSchedule(onset=50.0, amplitude=5.0)
        model = ProportionalHomeostat(params)
        history = model.simulate(2000, schedule=schedule)
        det = P24HomeostasisDetector(n_permutations=199, seed=42)
        result = det.detect(history, model.get_metadata())
        assert not result.detected, (
            f"P24 false positive on gain=0: {result.summary()}"
        )

    def test_gain_zero_growth_ratio_large(self):
        """Uncontrolled system has linearly growing deviation."""
        params = HomeostatParams(gain=0.0, setpoint=10.0, seed=42, dt=0.1)
        schedule = PerturbationSchedule(onset=50.0, amplitude=5.0)
        model = ProportionalHomeostat(params)
        history = model.simulate(2000, schedule=schedule)

        # Verify x drifts away from setpoint
        xs = np.array([h['x'] for h in history])
        assert xs[-1] > 100, f"Expected drift, got final x={xs[-1]}"


# ── T1a observation bundle ──

class TestT1aObservationBundle:
    """Observation bundle adapter contract."""

    def test_bundle_keys(self):
        params = HomeostatParams(gain=1.0, setpoint=10.0, seed=0, dt=0.1)
        schedule = PerturbationSchedule(onset=5.0, amplitude=3.0)
        model = ProportionalHomeostat(params)
        history = model.simulate(200, schedule=schedule)
        bundle = extract_observation_bundle(history)
        assert set(bundle.keys()) == {'time', 'x', 'setpoint', 'perturbation'}
        for k, v in bundle.items():
            assert isinstance(v, np.ndarray), f"{k} is not ndarray"
            assert v.shape == (200,), f"{k} shape is {v.shape}"

    def test_bundle_values_match_history(self):
        params = HomeostatParams(gain=1.0, setpoint=10.0, seed=0, dt=0.1)
        schedule = PerturbationSchedule(onset=5.0, amplitude=3.0)
        model = ProportionalHomeostat(params)
        history = model.simulate(100, schedule=schedule)
        bundle = extract_observation_bundle(history)
        for i in [0, 10, 50, 99]:
            assert bundle['x'][i] == pytest.approx(history[i]['x'])
            assert bundle['setpoint'][i] == pytest.approx(history[i]['setpoint'])
            assert bundle['perturbation'][i] == pytest.approx(history[i]['perturbation'])

    def test_detector_uses_bundle(self):
        """Detector reads through bundle adapter, not model object."""
        # Construct history dicts manually (no model object needed)
        history = []
        for i in range(500):
            t = i * 0.1
            pert = 5.0 if t >= 10.0 else 0.0
            # Simulate regulated x: deviation = pert/gain = 5/5 = 1
            x = 10.0 + (1.0 if t >= 10.0 else 0.0)
            history.append({
                'time': t, 'x': x, 'setpoint': 10.0,
                'perturbation': pert, 'step': i,
            })
        det = P24HomeostasisDetector(n_permutations=99, seed=42)
        result = det.detect(history)
        # Should detect — bounded deviation under perturbation
        assert result.detected, f"Not detected on manual bundle: {result.summary()}"


# ── Pulse perturbation ──

class TestP24PulsePerturbation:
    """Pulse perturbation (onset + offset) → recovery."""

    def test_pulse_detected(self):
        params = HomeostatParams(gain=5.0, setpoint=10.0, seed=42, dt=0.1)
        schedule = PerturbationSchedule(onset=50.0, offset=70.0, amplitude=5.0)
        model = ProportionalHomeostat(params)
        history = model.simulate(2000, schedule=schedule)
        det = P24HomeostasisDetector(n_permutations=199, seed=42)
        result = det.detect(history, model.get_metadata())
        assert result.detected, f"P24 not detected on pulse: {result.summary()}"


# ── Model metadata interaction ──

class TestP24MetadataInteraction:
    """Metadata no longer gates any tier; measured feedback does (f2d4044)."""

    def test_no_metadata_still_detects(self):
        """Detection works without metadata (at confirmation or lower)."""
        params = HomeostatParams(gain=5.0, setpoint=10.0, seed=42, dt=0.1)
        schedule = PerturbationSchedule(onset=50.0, amplitude=5.0)
        model = ProportionalHomeostat(params)
        history = model.simulate(2000, schedule=schedule)
        det = P24HomeostasisDetector(n_permutations=199, seed=42)
        result = det.detect(history, metadata=None)
        assert result.detected

    # Test maintenance: test_false_feedback_metadata_blocks_definitive was
    # replaced by the two tests below. f2d4044 removed the has_active_feedback
    # metadata gate on the definitive tier and made the confirmation tier
    # require MEASURED negative feedback (restoring_slope < 0, from regressing
    # dx/dt on x - setpoint over the perturbed window); S2 Text, P24 row:
    # "measured restoring feedback, no metadata flag".

    def test_feedback_metadata_flag_does_not_gate_tier(self):
        """The has_active_feedback flag no longer changes the tier (f2d4044)."""
        params = HomeostatParams(gain=5.0, setpoint=10.0, seed=42, dt=0.1)
        schedule = PerturbationSchedule(onset=50.0, amplitude=5.0)
        model = ProportionalHomeostat(params)
        history = model.simulate(2000, schedule=schedule)
        det = P24HomeostasisDetector(n_permutations=199, seed=42)
        tiers = {
            label: det.detect(history, metadata=md).tier.value
            for label, md in (("model", model.get_metadata()),
                              ("flag_false", {'has_active_feedback': False}),
                              ("none", None))
        }
        assert set(tiers.values()) == {"definitive"}, tiers
        result = det.detect(history, metadata={'has_active_feedback': False})
        assert result.primary_metric['restoring_slope'] < 0.0

    def test_no_negative_feedback_held_at_screening(self):
        """What now does the job: a trajectory without negative feedback is held
        at SCREENING by the measured restoring slope, even though it is bounded,
        stays near the setpoint and beats the uncontrolled-drift null.

        Synthetic slow runaway: deviation 0.3*exp(t/90) during a constant
        perturbation of 5 (positive feedback, restoring slope = +1/90). Its
        deviation ratio and null p would satisfy the definitive gates; only
        restoring_slope < 0 (required at confirmation since f2d4044) fails.
        """
        dt, n, onset_i, tau, sp = 0.1, 1000, 100, 90.0, 10.0
        t = np.arange(n) * dt
        x = sp + 0.3 * np.exp((t - t[onset_i]) / tau)
        pert = np.where(np.arange(n) >= onset_i, 5.0, 0.0)
        history = [{'time': float(t[i]), 'x': float(x[i]), 'setpoint': sp,
                    'perturbation': float(pert[i])} for i in range(n)]
        det = P24HomeostasisDetector(n_permutations=199, seed=42)
        result = det.detect(history, metadata=None)
        assert result.detected
        assert result.tier.value == "screening", result.summary()
        assert result.primary_metric['restoring_slope'] > 0.0
        # The other confirmation/definitive conditions are met:
        assert result.null_p_value <= 0.005
        assert result.secondary_metrics['deviation_ratio'] < 0.3
