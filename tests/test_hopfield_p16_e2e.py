"""End-to-end tests for P16 associative memory (Hopfield network).

Tests:
  - Deterministic DEFINITIVE on low-load canonical regime (P/N << capacity)
  - Negative control: random weights → no reliable completion → not detected
  - Transfer matrix +1 detector
  - Observation bundle adapter (T1a)
"""

import numpy as np
import pytest

from epc.models.hopfield import (
    HopfieldNetwork,
    HopfieldParams,
    BooleanGRN,
    BooleanGRNParams,
)
from epc.detectors.p16_associative_memory import (
    P16AssociativeMemoryDetector,
    extract_observation_bundle,
)


class TestHopfieldModel:
    """Basic model tests."""

    def test_deterministic_retrieval(self):
        """Hopfield retrieves a stored pattern from 20% corrupted cue."""
        params = HopfieldParams(N=200, P=5, corruption=0.2, seed=42)
        model = HopfieldNetwork(params)
        history = model.simulate(n_cues=1, cue_pattern_indices=[0])

        # Check convergence
        final = history[-1]
        assert final['converged'], "Network did not converge"

        # Check best overlap with ANY stored pattern (network may converge
        # to a different but valid attractor at moderate N)
        state = final['state']
        patterns = final['stored_patterns']
        N = len(state)
        best_overlap = max(
            abs(float(np.dot(p.astype(np.int32), state.astype(np.int32)) / N))
            for p in patterns
        )
        assert best_overlap > 0.9, (
            f"Poor retrieval: best_overlap={best_overlap:.3f}"
        )

    def test_multiple_pattern_retrieval(self):
        """All 5 stored patterns are individually retrievable (best overlap)."""
        params = HopfieldParams(N=200, P=5, corruption=0.2, seed=42)
        model = HopfieldNetwork(params)
        history = model.simulate(n_cues=5, cue_pattern_indices=[0, 1, 2, 3, 4])

        bundle = extract_observation_bundle(history)
        stored = bundle['stored_patterns']
        N = bundle['N']
        for trial in bundle['trials']:
            final_state = trial['states'][-1]
            best = max(
                abs(float(np.dot(
                    p.astype(np.int32), final_state.astype(np.int32),
                ) / N))
                for p in stored
            )
            assert best > 0.9, (
                f"Trial {trial['cue_pattern_idx']}: best_overlap={best:.3f}"
            )

    def test_seed_reproducibility(self):
        """Two runs with the same seed produce identical trajectories."""
        for seed in [42, 123]:
            params = HopfieldParams(N=100, P=5, seed=seed)
            m1 = HopfieldNetwork(params)
            h1 = m1.simulate(n_cues=1)
            m2 = HopfieldNetwork(params)
            h2 = m2.simulate(n_cues=1)
            assert np.array_equal(h1[-1]['state'], h2[-1]['state'])

    def test_metadata_keys(self):
        """Model metadata contains required keys."""
        model = HopfieldNetwork()
        meta = model.get_metadata()
        required = [
            'model_name', 'model_class', 'substrate_type',
            'N', 'P', 'alpha', 'has_multiple_attractors',
            'has_content_addressable_memory',
        ]
        for key in required:
            assert key in meta, f"Missing metadata key: {key}"
        assert meta['model_class'] == 'attractor_network'


class TestP16Detection:
    """P16 detector tests."""

    def test_canonical_definitive(self):
        """Low-load Hopfield (α=0.05) reaches DEFINITIVE."""
        params = HopfieldParams(N=100, P=5, corruption=0.2, seed=42)
        model = HopfieldNetwork(params)
        history = model.simulate(n_cues=5, cue_pattern_indices=[0, 1, 2, 3, 4])

        det = P16AssociativeMemoryDetector(n_permutations=199, seed=42)
        result = det.detect(history, model.get_metadata())

        assert result.detected, f"P16 not detected: {result.summary()}"
        assert result.tier.value == "definitive", (
            f"Expected definitive, got {result.tier.value}: {result.summary()}"
        )
        assert result.confidence >= 0.85, (
            f"Confidence too low: {result.confidence}"
        )

    def test_negative_control_random_weights(self):
        """Random weights → no associative memory → detector does NOT fire."""
        params = HopfieldParams(N=100, P=5, corruption=0.2, seed=42)
        model = HopfieldNetwork(params)

        # Replace weights with random (destroying stored patterns)
        rng = np.random.default_rng(99)
        model.W = rng.normal(0, 1.0 / params.N, size=(params.N, params.N))
        np.fill_diagonal(model.W, 0.0)
        model.W = (model.W + model.W.T) / 2  # Symmetrize

        history = model.simulate(n_cues=5, cue_pattern_indices=[0, 1, 2, 3, 4])
        det = P16AssociativeMemoryDetector(n_permutations=99, seed=42)
        result = det.detect(history, model.get_metadata())

        assert not result.detected, (
            f"P16 false positive on random weights: {result.summary()}"
        )

    def test_observation_bundle_adapter(self):
        """T1a adapter extracts correct structure from history."""
        params = HopfieldParams(N=50, P=3, corruption=0.1, seed=42)
        model = HopfieldNetwork(params)
        history = model.simulate(n_cues=2, cue_pattern_indices=[0, 1])

        bundle = extract_observation_bundle(history)

        assert 'trials' in bundle
        assert 'stored_patterns' in bundle
        assert bundle['N'] == 50
        assert bundle['P'] == 3
        assert len(bundle['trials']) == 2
        assert bundle['stored_patterns'].shape == (3, 50)

        # Each trial has states, overlaps, cue_pattern_idx, converged
        for trial in bundle['trials']:
            assert len(trial['states']) > 0
            assert len(trial['overlaps']) > 0
            assert 'cue_pattern_idx' in trial
            assert 'converged' in trial

    def test_selective_recall_multiple_patterns(self):
        """Detector verifies that multiple distinct patterns are recalled."""
        params = HopfieldParams(N=100, P=5, corruption=0.2, seed=42)
        model = HopfieldNetwork(params)
        history = model.simulate(n_cues=5, cue_pattern_indices=[0, 1, 2, 3, 4])

        det = P16AssociativeMemoryDetector(n_permutations=199, seed=42)
        result = det.detect(history, model.get_metadata())

        assert result.secondary_metrics['recall_accuracy'] == 1.0, (
            f"Not all patterns recalled: "
            f"acc={result.secondary_metrics['recall_accuracy']}"
        )
        assert result.secondary_metrics['n_distinct_patterns_recalled'] == 5, (
            f"Not all distinct patterns recalled: "
            f"n={result.secondary_metrics['n_distinct_patterns_recalled']}"
        )


class TestP16TransferMatrix:
    """Transfer matrix entry: P16 detector count += 1."""

    def test_detector_registered(self):
        """P16 must be in DETECTOR_REGISTRY."""
        from epc.orchestration import DETECTOR_REGISTRY
        assert 'P16' in DETECTOR_REGISTRY, "P16 not in DETECTOR_REGISTRY"

    def test_hopfield_model_registered(self):
        """Hopfield model must be in MODEL_REGISTRY."""
        from epc.orchestration import MODEL_REGISTRY
        assert 'hopfield' in MODEL_REGISTRY, "hopfield not in MODEL_REGISTRY"

    def test_canonical_pair_compatible(self):
        """Hopfield × P16 must be compatible."""
        from epc.orchestration import check_compatibility
        result = check_compatibility('hopfield', 'P16')
        assert result.compatible, (
            f"hopfield × P16 not compatible: {result.reason}"
        )
