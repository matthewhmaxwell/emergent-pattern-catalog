"""Orchestration layer tests — substrate-aware detector dispatch.

Tests the compatibility matrix: 5 substrate types, 11 models × 10 detectors,
24 compatible pairs. Verifies substrate filtering, observable guards, canonical
positive/negative assignments, and that Sprint 5 models (Nowak-May, HK) are
correctly registered.

Architecture decision #25 (updated Sprint 5 verification).
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
        assert len(substrates) == 5, f"Expected 5 substrate types, got {substrates}"


class TestRegistryCounts:

    def test_model_count(self):
        assert len(MODEL_REGISTRY) == 11, \
            f"Expected 11 models, got {len(MODEL_REGISTRY)}: {list(MODEL_REGISTRY.keys())}"

    def test_detector_count(self):
        assert len(DETECTOR_REGISTRY) == 10, \
            f"Expected 10 detectors, got {len(DETECTOR_REGISTRY)}: {list(DETECTOR_REGISTRY.keys())}"


class TestCompatibility:

    def test_total_compatible_pairs(self):
        pairs = get_compatible_pairs()
        assert len(pairs) == 24, \
            f"Expected 24 compatible pairs, got {len(pairs)}: {pairs}"

    def test_total_cells(self):
        matrix = get_compatibility_matrix()
        total = sum(len(row) for row in matrix.values())
        assert total == 110, f"Expected 110 cells (11x10), got {total}"

    def test_mismatch_count(self):
        pairs = get_compatible_pairs()
        mismatches = 110 - len(pairs)
        assert mismatches == 86, f"Expected 86 mismatches, got {mismatches}"


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
