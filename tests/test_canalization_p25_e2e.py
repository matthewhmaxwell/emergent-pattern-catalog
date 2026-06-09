"""End-to-end tests for P25 canalized restoration / equifinality.

Tests:
  1. Canonical positive: CanalizedLandscape reaches DEFINITIVE.
  2. Negative: DiffusiveDynamics (ICs do NOT converge) → NOT detected.
  3. Negative: TrivialCollapse (instant map, no dynamics) → NOT detected.
  4. Perturbation recovery: CanalizedLandscape with mid-trajectory perturbation
     restores to target → DEFINITIVE with restoration_accuracy ≥ 0.8.
  5. Determinism: two runs with same seed produce identical results.
"""

import numpy as np
import pytest

from epc.models.canalization import (
    CanalizedLandscape,
    CanalizedLandscapeParams,
    DiffusiveDynamics,
    MultiBasinGRN,
    TrivialCollapse,
)
from epc.detectors.p25_equifinality import (
    P25EquifinalityDetector,
    extract_observation_bundle,
)


class TestP25CanonicalPositive:
    """P25 on canonical CanalizedLandscape: DEFINITIVE detection."""

    def test_canalized_landscape_detected(self):
        """P25 detects canalized convergence at definitive tier."""
        params = CanalizedLandscapeParams(
            n_dims=10, n_ics=20, n_steps=200, seed=42,
            basin_strength=2.0, ic_spread=5.0,
        )
        model = CanalizedLandscape(params)
        h = model.simulate()
        det = P25EquifinalityDetector(n_permutations=199, seed=42)
        r = det.detect(h, model.get_metadata())

        assert r.detected, f"P25 not detected: {r.summary()}"
        assert r.tier.value == "definitive", f"Expected definitive, got {r.tier.value}"
        assert r.confidence >= 0.75, f"Confidence too low: {r.confidence}"

    def test_convergence_variance_ratio_small(self):
        """Final-state variance ≪ IC variance (ratio < 0.01)."""
        params = CanalizedLandscapeParams(
            n_dims=10, n_ics=20, n_steps=200, seed=42,
        )
        model = CanalizedLandscape(params)
        h = model.simulate()
        det = P25EquifinalityDetector(n_permutations=99, seed=42)
        r = det.detect(h, model.get_metadata())

        ratio = r.primary_metric['convergence_variance_ratio']
        assert ratio < 0.01, f"Convergence ratio {ratio} not ≪ 1"

    def test_basin_volume_high(self):
        """Basin volume ≥ 0.9 (most ICs converge)."""
        params = CanalizedLandscapeParams(
            n_dims=10, n_ics=20, n_steps=200, seed=42,
        )
        model = CanalizedLandscape(params)
        h = model.simulate()
        det = P25EquifinalityDetector(n_permutations=99, seed=42)
        r = det.detect(h, model.get_metadata())

        basin = r.secondary_metrics['basin_volume']
        assert basin >= 0.9, f"Basin volume {basin} too low"

    def test_relaxation_time_nontrivial(self):
        """Relaxation time > 1 step (not trivial collapse)."""
        params = CanalizedLandscapeParams(
            n_dims=10, n_ics=20, n_steps=200, seed=42,
        )
        model = CanalizedLandscape(params)
        h = model.simulate()
        det = P25EquifinalityDetector(n_permutations=99, seed=42)
        r = det.detect(h, model.get_metadata())

        relax = r.secondary_metrics['relaxation_time']
        assert relax > 1.0, f"Relaxation time {relax} too short (trivial?)"


class TestP25NegativeControls:
    """Negative controls: detector must NOT fire on these."""

    def test_diffusive_not_detected(self):
        """Diffusive dynamics (ICs diverge) → NOT detected."""
        model = DiffusiveDynamics(
            n_dims=10, n_ics=20, n_steps=200, seed=42,
        )
        h = model.simulate()
        det = P25EquifinalityDetector(n_permutations=99, seed=42)
        r = det.detect(h, model.get_metadata())

        assert not r.detected, f"P25 false positive on diffusive: {r.summary()}"

    def test_trivial_collapse_not_detected(self):
        """Trivial collapse (instant map to constant) → NOT detected.

        The trivial-collapse gate in the detector requires dynamics-driven
        convergence (relaxation_time > 1); instant collapse is rejected.
        """
        model = TrivialCollapse(
            n_dims=10, n_ics=20, n_steps=200, seed=42,
        )
        h = model.simulate()
        det = P25EquifinalityDetector(n_permutations=99, seed=42)
        r = det.detect(h, model.get_metadata())

        assert not r.detected, (
            f"P25 false positive on trivial collapse: {r.summary()}"
        )


class TestP25PerturbationRecovery:
    """Perturbation recovery: canalized system restores after displacement."""

    def test_restoration_accuracy(self):
        """Restoration accuracy ≥ 0.8 after mid-trajectory perturbation."""
        params = CanalizedLandscapeParams(
            n_dims=10, n_ics=20, n_steps=200, seed=42,
            perturbation_time=50, perturbation_magnitude=3.0,
        )
        model = CanalizedLandscape(params)
        h = model.simulate()
        det = P25EquifinalityDetector(n_permutations=199, seed=42)
        r = det.detect(h, model.get_metadata())

        assert r.detected, f"P25 not detected with perturbation: {r.summary()}"
        restore = r.secondary_metrics.get('restoration_accuracy')
        assert restore is not None, "No restoration_accuracy in secondaries"
        assert restore >= 0.8, f"Restoration accuracy {restore} < 0.8"


class TestP25Determinism:
    """Determinism: same seed → identical results."""

    def test_deterministic_results(self):
        """Two runs with identical seeds produce identical detection results."""
        params = CanalizedLandscapeParams(n_dims=10, n_ics=20, n_steps=200, seed=42)

        for _ in range(2):
            model = CanalizedLandscape(params)
            h = model.simulate()
            det = P25EquifinalityDetector(n_permutations=99, seed=42)
            r = det.detect(h, model.get_metadata())

        # Run a second time and compare
        model2 = CanalizedLandscape(params)
        h2 = model2.simulate()
        det2 = P25EquifinalityDetector(n_permutations=99, seed=42)
        r2 = det2.detect(h2, model2.get_metadata())

        assert r.tier == r2.tier
        assert r.confidence == r2.confidence
        assert r.null_p_value == r2.null_p_value


class TestP25ObservationBundle:
    """T1a observation bundle adapter."""

    def test_extract_observation_bundle(self):
        """Bundle extraction produces correct shapes."""
        params = CanalizedLandscapeParams(n_dims=5, n_ics=10, n_steps=50, seed=42)
        model = CanalizedLandscape(params)
        h = model.simulate()
        bundle = extract_observation_bundle(h)

        assert bundle['ics'].shape == (10, 5)
        assert bundle['finals'].shape == (10, 5)
        assert bundle['target'].shape == (5,)
        assert bundle['n_trials'] == 10
        assert bundle['n_dims'] == 5
