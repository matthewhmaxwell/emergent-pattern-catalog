"""End-to-end tests for P29 trail / network formation.

Tests:
  - Deterministic detection on ant trail model (DEFINITIVE expected).
  - Negative control: no-reinforcement model must NOT be detected.
  - Observation bundle (T1a) adapter contract.
  - Metadata correctness.
  - Transfer-matrix registration check.
"""

import numpy as np
import pytest

from epc.models.trail_network import (
    AntTrailModel,
    AntTrailParams,
    NoReinforcementModel,
    NoReinforcementParams,
    PhysarumModel,
    PhysarumParams,
)
from epc.detectors.p29_trail_network import (
    P29TrailNetworkDetector,
    extract_observation_bundle,
)
from epc.orchestration import DETECTOR_REGISTRY, MODEL_REGISTRY


class TestP29Deterministic:
    """P29 detection on canonical ant trail model."""

    def test_ant_trail_detected(self):
        """Ant trail model should reach at least CONFIRMATION tier."""
        params = AntTrailParams(
            n_nodes=7, n_agents=50, grid_size=100,
            alpha=1.0, beta=2.0, deposition_rate=10.0,
            evaporation_rate=0.02, n_steps=500, snapshot_interval=50,
            node_layout="grid", seed=42,
        )
        model = AntTrailModel(params)
        history = model.simulate()
        detector = P29TrailNetworkDetector(n_permutations=199, seed=42)
        result = detector.detect(history, model.get_metadata())

        assert result.detected, f"P29 not detected on ant trail: {result.summary()}"
        assert result.tier.value in ("confirmation", "definitive"), (
            f"Expected confirmation+; got {result.tier.value}"
        )
        assert result.confidence >= 0.50

    def test_physarum_definitive(self):
        """Physarum model should reach DEFINITIVE tier."""
        params = PhysarumParams(
            n_nodes=5, grid_size=80, gamma=1.8, decay_rate=0.01,
            n_steps=2000, snapshot_interval=50,
            node_layout="grid", seed=42,
        )
        model = PhysarumModel(params)
        history = model.simulate()
        detector = P29TrailNetworkDetector(n_permutations=199, seed=42)
        result = detector.detect(history, model.get_metadata())

        assert result.detected, f"P29 not detected on Physarum: {result.summary()}"
        assert result.tier.value == "definitive", (
            f"Expected definitive; got {result.tier.value}"
        )
        assert result.confidence >= 0.75

    def test_ant_trail_length_ratio(self):
        """Ant trail network length should be within bounded factor of MST."""
        params = AntTrailParams(
            n_nodes=7, n_agents=50, grid_size=100,
            alpha=1.0, beta=2.0, deposition_rate=10.0,
            evaporation_rate=0.02, n_steps=500, snapshot_interval=50,
            node_layout="grid", seed=42,
        )
        model = AntTrailModel(params)
        history = model.simulate()
        detector = P29TrailNetworkDetector(n_permutations=99, seed=42)
        result = detector.detect(history, model.get_metadata())

        length_ratio = result.primary_metric['length_ratio']
        assert length_ratio < 2.5, f"Network length ratio too high: {length_ratio}"

    def test_ant_trail_connectivity(self):
        """Ant trail network should be well-connected."""
        params = AntTrailParams(
            n_nodes=7, n_agents=50, grid_size=100,
            alpha=1.0, beta=2.0, deposition_rate=10.0,
            evaporation_rate=0.02, n_steps=500, snapshot_interval=50,
            node_layout="grid", seed=42,
        )
        model = AntTrailModel(params)
        history = model.simulate()
        detector = P29TrailNetworkDetector(n_permutations=99, seed=42)
        result = detector.detect(history, model.get_metadata())

        assert result.primary_metric['connectivity'] >= 0.6


class TestP29NegativeControl:
    """No-reinforcement model must NOT be detected as P29."""

    def test_no_reinforcement_rejected(self):
        """No-reinforcement random walk should not produce efficient network."""
        params = NoReinforcementParams(
            n_nodes=7, n_agents=50, grid_size=100,
            n_steps=500, snapshot_interval=50,
            node_layout="grid", seed=42,
        )
        model = NoReinforcementModel(params)
        history = model.simulate()
        detector = P29TrailNetworkDetector(n_permutations=99, seed=42)
        result = detector.detect(history, model.get_metadata())

        assert not result.detected, (
            f"P29 false positive on no-reinforcement model: {result.summary()}"
        )


class TestP29ObservationBundle:
    """T1a adapter tests."""

    def test_bundle_keys(self):
        """Bundle should contain required keys."""
        params = AntTrailParams(
            n_nodes=4, n_agents=20, grid_size=50,
            n_steps=100, snapshot_interval=50, seed=0,
        )
        model = AntTrailModel(params)
        history = model.simulate()
        bundle = extract_observation_bundle(history)

        required = {'node_positions', 'edge_weights', 'pheromone_fields',
                     'steps', 'n_nodes', 'grid_size'}
        assert required.issubset(set(bundle.keys()))

    def test_bundle_shapes(self):
        """Bundle arrays should have correct shapes."""
        params = AntTrailParams(
            n_nodes=4, n_agents=20, grid_size=50,
            n_steps=100, snapshot_interval=50, seed=0,
        )
        model = AntTrailModel(params)
        history = model.simulate()
        bundle = extract_observation_bundle(history)

        T = len(bundle['steps'])
        N = bundle['n_nodes']
        G = bundle['grid_size']

        assert len(bundle['node_positions']) == T
        assert bundle['node_positions'][0].shape == (N, 2)
        assert len(bundle['edge_weights']) == T
        assert bundle['edge_weights'][0].shape == (N, N)
        assert len(bundle['pheromone_fields']) == T
        assert bundle['pheromone_fields'][0].shape == (G, G)


class TestP29Metadata:
    """Metadata correctness checks."""

    def test_ant_trail_metadata(self):
        model = AntTrailModel()
        meta = model.get_metadata()
        assert meta['has_pheromone_reinforcement'] is True
        assert meta['model_class'] == 'trail_network'
        assert meta['substrate'] == 'trail_network'

    def test_physarum_metadata(self):
        model = PhysarumModel()
        meta = model.get_metadata()
        assert meta['has_pheromone_reinforcement'] is True
        assert meta['model_class'] == 'trail_network'

    def test_no_reinforcement_metadata(self):
        model = NoReinforcementModel()
        meta = model.get_metadata()
        assert meta['has_pheromone_reinforcement'] is False
        assert meta['model_class'] == 'random_walk'


class TestP29TransferMatrix:
    """Registry and compatibility checks."""

    def test_p29_in_detector_registry(self):
        assert 'P29' in DETECTOR_REGISTRY

    def test_ant_trail_in_model_registry(self):
        assert 'ant_trail_network' in MODEL_REGISTRY

    def test_physarum_in_model_registry(self):
        assert 'physarum_network' in MODEL_REGISTRY

    def test_p29_compatible_with_ant_trail(self):
        from epc.orchestration import check_compatibility
        result = check_compatibility('ant_trail_network', 'P29')
        assert result.compatible, f"Not compatible: {result.reason}"

    def test_p29_compatible_with_physarum(self):
        from epc.orchestration import check_compatibility
        result = check_compatibility('physarum_network', 'P29')
        assert result.compatible, f"Not compatible: {result.reason}"

    def test_p29_incompatible_with_lattice_2d(self):
        from epc.orchestration import check_compatibility
        result = check_compatibility('greenberg_hastings', 'P29')
        assert not result.compatible
