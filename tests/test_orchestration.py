"""Orchestration layer tests — substrate-aware detector dispatch.

Tests the compatibility matrix: 6 substrate types, 14 models × 13 detectors,
43 compatible pairs (Sprint 13 added gray_scott × P3 on the new
lattice_2d_continuous substrate). Verifies substrate filtering, observable
guards, canonical positive/negative assignments, and that Sprint 5 models
(Nowak-May, HK) are correctly registered.

Architecture decision #25 (updated Sprint 13).
"""

import pytest

from epc.orchestration import (
    MODEL_REGISTRY,
    DETECTOR_REGISTRY,
    check_compatibility,
    get_compatible_pairs,
    get_compatibility_matrix,
)


# ---------------------------------------------------------------------------
# Substrate type counts
# ---------------------------------------------------------------------------

class TestSubstrateCounts:

    def test_lattice_1d_models_exist(self):
        models = [m for m in MODEL_REGISTRY.values()
                  if m.substrate_type == 'lattice_1d']
        assert len(models) >= 1

    def test_lattice_2d_models_exist(self):
        models = [m for m in MODEL_REGISTRY.values()
                  if m.substrate_type == 'lattice_2d']
        assert len(models) >= 4

    def test_continuous_2d_models_exist(self):
        models = [m for m in MODEL_REGISTRY.values()
                  if m.substrate_type == 'continuous_2d']
        assert len(models) >= 2

    def test_oscillator_models_exist(self):
        models = [m for m in MODEL_REGISTRY.values()
                  if m.substrate_type == 'oscillator']
        assert len(models) >= 1

    def test_opinion_space_models_exist(self):
        models = [m for m in MODEL_REGISTRY.values()
                  if m.substrate_type == 'opinion_space']
        assert len(models) >= 1

    def test_five_substrate_types_total(self):
        substrates = {m.substrate_type for m in MODEL_REGISTRY.values()}
        assert len(substrates) == 6, f"Expected 6 substrate types, got {substrates}"


class TestRegistryCounts:

    def test_model_count(self):
        assert len(MODEL_REGISTRY) == 16, \
            f"Expected 16 models, got {len(MODEL_REGISTRY)}: {list(MODEL_REGISTRY.keys())}"

    def test_detector_count(self):
        assert len(DETECTOR_REGISTRY) == 15, \
            f"Expected 15 detectors, got {len(DETECTOR_REGISTRY)}: {list(DETECTOR_REGISTRY.keys())}"


class TestCompatibility:

    def test_total_compatible_pairs(self):
        pairs = get_compatible_pairs()
        assert len(pairs) == 49, \
            f"Expected 49 compatible pairs, got {len(pairs)}: {pairs}"

    def test_total_cells(self):
        matrix = get_compatibility_matrix()
        total = sum(len(row) for row in matrix.values())
        assert total == 240, f"Expected 240 cells (16x15), got {total}"

    def test_mismatch_count(self):
        pairs = get_compatible_pairs()
        mismatches = 240 - len(pairs)
        assert mismatches == 191, f"Expected 191 mismatches, got {mismatches}"


class TestCanonicalPairs:

    @pytest.mark.parametrize("model_name,detector_id", [
        ('schelling', 'P1'),
        ('vicsek', 'P5'),
        ('dorsogna', 'P6'),
        ('kuramoto', 'P9'),
        ('greenberg_hastings', 'P13'),
        ('btw_sandpile', 'P14'),
        ('game_of_life', 'P15'),
        ('hegselmann_krause', 'P21'),
        ('nowak_may', 'P27'),
        ('zhang_sequential', 'P31'),
        ('sir_epidemic', 'P22'),
    ])
    def test_canonical_pair_compatible(self, model_name, detector_id):
        result = check_compatibility(model_name, detector_id)
        assert result.compatible, \
            f"{model_name} x {detector_id} should be compatible: {result.reason}"


class TestObservableGuards:

    def test_p27_requires_coop_fraction(self):
        p27 = DETECTOR_REGISTRY['P27']
        assert 'coop_fraction' in p27.required_observables

    def test_p27_not_compatible_with_gh(self):
        result = check_compatibility('greenberg_hastings', 'P27')
        assert not result.compatible, "P27 should NOT match GH"

    def test_p27_compatible_with_nowak_may(self):
        result = check_compatibility('nowak_may', 'P27')
        assert result.compatible, f"Nowak-May should match P27: {result.reason}"

    def test_p14_only_btw(self):
        p14_pairs = [(m, d) for m, d in get_compatible_pairs() if d == 'P14']
        assert len(p14_pairs) == 1 and p14_pairs[0][0] == 'btw_sandpile'

    def test_p21_only_hk(self):
        p21_pairs = [(m, d) for m, d in get_compatible_pairs() if d == 'P21']
        assert len(p21_pairs) == 1 and p21_pairs[0][0] == 'hegselmann_krause'


class TestSprint9Registrations:
    """Sprint 9: rps_spatial model + P12 detector are registered and
    compatible with the expected counterparts."""

    def test_rps_spatial_registered(self):
        assert 'rps_spatial' in MODEL_REGISTRY, \
            "rps_spatial model should be registered"
        m = MODEL_REGISTRY['rps_spatial']
        assert m.substrate_type == 'lattice_2d'
        assert 'grid' in m.observables
        assert 'P12' in m.primary_patterns

    def test_p12_registered(self):
        assert 'P12' in DETECTOR_REGISTRY, "P12 detector should be registered"
        d = DETECTOR_REGISTRY['P12']
        assert 'lattice_2d' in d.required_substrate
        assert 'grid' in d.required_observables
        assert d.observable_scope == 'state_history_only'

    def test_rps_p12_compatible(self):
        """The canonical Sprint 9 pair: rps_spatial × P12 must be compatible."""
        from epc.orchestration import check_compatibility
        r = check_compatibility('rps_spatial', 'P12')
        assert r.compatible, f"rps_spatial × P12 should be compatible: {r.reason}"

    def test_p12_compatible_with_all_lattice_2d_grid_models(self):
        """P12 is a lattice_2d detector requiring 'grid' observable. Every
        lattice_2d model exposing a 'grid' observable should pair.

        Note: btw_sandpile is lattice_2d but exposes only avalanche-level
        observables (no 'grid'), so it is correctly REJECTED by observable
        mismatch, not substrate mismatch. This is expected behavior.
        """
        p12_pairs = {m for m, d in get_compatible_pairs() if d == 'P12'}
        lattice_2d_grid_models = {
            name for name, reg in MODEL_REGISTRY.items()
            if reg.substrate_type == 'lattice_2d' and 'grid' in reg.observables
        }
        assert p12_pairs == lattice_2d_grid_models, (
            f"P12 pairs {p12_pairs} should match lattice_2d-with-grid models "
            f"{lattice_2d_grid_models}"
        )

    def test_p12_rejects_non_lattice_2d_substrates(self):
        """P12 should NOT be compatible with oscillator/continuous_2d/etc."""
        from epc.orchestration import check_compatibility
        for model in ['kuramoto', 'vicsek', 'dorsogna',
                      'hegselmann_krause', 'zhang_sequential']:
            r = check_compatibility(model, 'P12')
            assert not r.compatible, \
                f"{model} (non-lattice_2d) should NOT match P12"
