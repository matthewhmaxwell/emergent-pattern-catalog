"""End-to-end tests for P4 territoriality / exclusion boundaries.

Tests:
  - Deterministic: DEFINITIVE on canonical ScentMarkingModel
  - Negative control: PlainRandomWalkModel → NOT detected
  - Boundary persistence: stable partition over time window
  - Observation bundle adapter (T1a)
"""

from __future__ import annotations

import numpy as np
import pytest

from epc.models.territoriality import (
    ScentMarkingModel,
    ScentMarkingParams,
    PlainRandomWalkModel,
)
from epc.detectors.p4_territoriality import (
    P4TerritorialityDetector,
    extract_observation_bundle,
    _exclusivity_index,
    _pairwise_overlap,
    _boundary_persistence,
)


class TestP4Deterministic:
    """Deterministic seed-pinned P4 detection."""

    def test_scent_marking_definitive(self):
        """ScentMarkingModel → DEFINITIVE detection."""
        params = ScentMarkingParams(n_steps=20000, snapshot_interval=100, seed=42)
        model = ScentMarkingModel(params)
        history = model.simulate()
        metadata = model.get_metadata()

        det = P4TerritorialityDetector(n_permutations=199, seed=42)
        result = det.detect(history, metadata)

        assert result.detected
        assert result.tier.value == "definitive"
        assert result.confidence >= 0.75
        assert result.null_p_value <= 0.01

    def test_scent_marking_exclusivity_high(self):
        """Territorial model should have exclusivity > 0.85."""
        params = ScentMarkingParams(n_steps=20000, snapshot_interval=100, seed=42)
        model = ScentMarkingModel(params)
        history = model.simulate()

        final_occ = history[-1]['occupancy']
        excl = _exclusivity_index(final_occ)
        assert excl > 0.85, f"Exclusivity too low: {excl:.3f}"

    def test_scent_marking_low_overlap(self):
        """Territorial model should have low pairwise overlap."""
        params = ScentMarkingParams(n_steps=20000, snapshot_interval=100, seed=42)
        model = ScentMarkingModel(params)
        history = model.simulate()

        final_occ = history[-1]['occupancy']
        overlap = _pairwise_overlap(final_occ)
        assert overlap < 0.10, f"Overlap too high: {overlap:.3f}"


class TestP4NegativeControl:
    """PlainRandomWalkModel → NOT detected."""

    def test_random_walk_not_detected(self):
        """Plain random walk → no territorial detection."""
        model = PlainRandomWalkModel(
            n_agents=4, grid_size=48, n_steps=5000,
            snapshot_interval=50, seed=42,
        )
        history = model.simulate()
        metadata = model.get_metadata()

        det = P4TerritorialityDetector(n_permutations=199, seed=42)
        result = det.detect(history, metadata)

        assert not result.detected, (
            f"False positive on random walk: tier={result.tier}, "
            f"excl={result.primary_metric.get('exclusivity_index', '?')}"
        )


class TestP4BoundaryPersistence:
    """Territory boundaries should be stable over time."""

    def test_persistence_high(self):
        """Boundary persistence > 0.7 for territorial model."""
        params = ScentMarkingParams(n_steps=20000, snapshot_interval=100, seed=42)
        model = ScentMarkingModel(params)
        history = model.simulate()

        occ_list = [h['occupancy'] for h in history]
        n = len(occ_list)
        early = occ_list[n // 4]
        late = occ_list[-1]
        pers = _boundary_persistence(early, late)
        assert pers > 0.7, f"Boundary persistence too low: {pers:.3f}"


class TestP4ObservationBundle:
    """T1a observation bundle adapter tests."""

    def test_bundle_keys(self):
        """Bundle has required keys."""
        params = ScentMarkingParams(n_steps=500, snapshot_interval=50, seed=42)
        model = ScentMarkingModel(params)
        history = model.simulate()

        bundle = extract_observation_bundle(history)
        assert 'positions' in bundle
        assert 'scent_fields' in bundle
        assert 'occupancy' in bundle
        assert 'steps' in bundle
        assert 'n_agents' in bundle
        assert 'grid_size' in bundle

    def test_bundle_shapes(self):
        """Bundle arrays have correct shapes."""
        params = ScentMarkingParams(
            n_agents=3, grid_size=32, n_steps=500,
            snapshot_interval=50, seed=42,
        )
        model = ScentMarkingModel(params)
        history = model.simulate()

        bundle = extract_observation_bundle(history)
        assert bundle['n_agents'] == 3
        assert bundle['grid_size'] == 32
        assert bundle['positions'][0].shape == (3, 2)
        assert bundle['scent_fields'][0].shape == (3, 32, 32)
        assert bundle['occupancy'][0].shape == (3, 32, 32)


class TestP4Metadata:
    """Model metadata tests."""

    def test_scent_marking_metadata(self):
        """ScentMarkingModel metadata has required fields."""
        model = ScentMarkingModel()
        meta = model.get_metadata()
        assert meta['model_class'] == 'territorial_random_walk'
        assert meta['has_foreign_avoidance'] is True
        assert meta['has_scent_marking'] is True

    def test_random_walk_metadata(self):
        """PlainRandomWalkModel metadata indicates no avoidance."""
        model = PlainRandomWalkModel()
        meta = model.get_metadata()
        assert meta['has_foreign_avoidance'] is False


class TestP4TransferMatrix:
    """Transfer-matrix coverage for P4."""

    def test_p4_registered(self):
        """P4 is in DETECTOR_REGISTRY."""
        from epc.orchestration import DETECTOR_REGISTRY
        assert 'P4' in DETECTOR_REGISTRY

    def test_scent_marking_registered(self):
        """scent_marking_territory is in MODEL_REGISTRY."""
        from epc.orchestration import MODEL_REGISTRY
        assert 'scent_marking_territory' in MODEL_REGISTRY

    def test_pheromone_repulsion_registered(self):
        """pheromone_repulsion_territory is in MODEL_REGISTRY."""
        from epc.orchestration import MODEL_REGISTRY
        assert 'pheromone_repulsion_territory' in MODEL_REGISTRY
