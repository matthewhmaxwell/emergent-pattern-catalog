"""Orchestration layer tests — substrate-aware detector dispatch.

Tests the compatibility matrix: 7 substrate types, 19 models × 18 detectors,
65 compatible pairs (Sprint 19 added the lotka_volterra model + P11 detector
that were shipped in Sprint 11 but never registered in the orchestration
layer — closes the oldest carry-forward). Verifies substrate filtering,
observable guards, canonical positive/negative assignments, and that
Sprint 5 models (Nowak-May, HK) are correctly registered.

Architecture decision #25 (updated Sprint 19).
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

    def test_scalar_wealth_models_exist(self):
        models = [m for m in MODEL_REGISTRY.values()
                  if m.substrate_type == 'scalar_wealth']
        assert len(models) >= 1

    def test_five_substrate_types_total(self):
        substrates = {m.substrate_type for m in MODEL_REGISTRY.values()}
        assert len(substrates) == 7, f"Expected 7 substrate types, got {substrates}"


class TestRegistryCounts:

    def test_model_count(self):
        assert len(MODEL_REGISTRY) == 19, \
            f"Expected 19 models, got {len(MODEL_REGISTRY)}: {list(MODEL_REGISTRY.keys())}"

    def test_detector_count(self):
        assert len(DETECTOR_REGISTRY) == 18, \
            f"Expected 18 detectors, got {len(DETECTOR_REGISTRY)}: {list(DETECTOR_REGISTRY.keys())}"


class TestCompatibility:

    def test_total_compatible_pairs(self):
        pairs = get_compatible_pairs()
        assert len(pairs) == 65, \
            f"Expected 65 compatible pairs, got {len(pairs)}: {pairs}"

    def test_total_cells(self):
        matrix = get_compatibility_matrix()
        total = sum(len(row) for row in matrix.values())
        assert total == 342, f"Expected 342 cells (19x18), got {total}"

    def test_mismatch_count(self):
        pairs = get_compatible_pairs()
        mismatches = 342 - len(pairs)
        assert mismatches == 277, f"Expected 277 mismatches, got {mismatches}"


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
        ('yard_sale', 'P28'),
        ('kuramoto_nonlocal', 'P10'),
        ('lotka_volterra', 'P11'),
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


class TestSprint17Registrations:
    """Sprint 17: yard_sale model + P28 detector are registered on the
    new scalar_wealth substrate and compatible with each other; P28
    rejects all 16 other substrates at substrate_mismatch."""

    def test_yard_sale_registered(self):
        assert 'yard_sale' in MODEL_REGISTRY, \
            "yard_sale model should be registered"
        m = MODEL_REGISTRY['yard_sale']
        assert m.substrate_type == 'scalar_wealth'
        assert 'wealth' in m.observables
        assert 'P28' in m.primary_patterns

    def test_yard_sale_metadata_rule_flags(self):
        """P28 mechanistic-null gate depends on four metadata flags.

        These MUST be declared in the yard_sale MODEL_REGISTRY entry so
        downstream detector code can rely on their presence (Decision 49).
        """
        m = MODEL_REGISTRY['yard_sale']
        for flag in ['has_conserved_resource', 'has_multiplicative_stake',
                     'has_saving_propensity', 'has_redistribution']:
            assert flag in m.metadata_keys, \
                f"yard_sale must declare metadata flag '{flag}' for P28"

    def test_p28_registered(self):
        assert 'P28' in DETECTOR_REGISTRY, "P28 detector should be registered"
        d = DETECTOR_REGISTRY['P28']
        assert 'scalar_wealth' in d.required_substrate
        assert 'wealth' in d.required_observables
        assert d.observable_scope == 'model_metadata_assisted'

    def test_yard_sale_p28_compatible(self):
        """The canonical Sprint 17 pair: yard_sale × P28 must be compatible."""
        r = check_compatibility('yard_sale', 'P28')
        assert r.compatible, f"yard_sale × P28 should be compatible: {r.reason}"

    def test_p28_rejects_all_non_wealth_substrates(self):
        """P28 should reject every substrate except scalar_wealth at
        substrate_mismatch. This exercises the full registry for the
        new substrate type."""
        for m in MODEL_REGISTRY:
            reg = MODEL_REGISTRY[m]
            r = check_compatibility(m, 'P28')
            if reg.substrate_type == 'scalar_wealth':
                assert r.compatible, f"{m} (scalar_wealth) should match P28"
            else:
                assert not r.compatible, \
                    f"{m} ({reg.substrate_type}) should NOT match P28"
                assert 'substrate_mismatch' in r.reason, \
                    f"{m} × P28 rejection reason should be substrate_mismatch"

    def test_p28_only_pairs_with_yard_sale(self):
        """P28 should have exactly one compatible model (yard_sale) at Sprint 17."""
        p28_pairs = [(m, d) for m, d in get_compatible_pairs() if d == 'P28']
        assert len(p28_pairs) == 1 and p28_pairs[0][0] == 'yard_sale'

    def test_yard_sale_only_pairs_with_p28(self):
        """yard_sale is on the new scalar_wealth substrate, so no other
        detector should match it at Sprint 17."""
        ys_pairs = [(m, d) for m, d in get_compatible_pairs() if m == 'yard_sale']
        assert len(ys_pairs) == 1 and ys_pairs[0][1] == 'P28'

    def test_scalar_wealth_substrate_exclusive(self):
        """scalar_wealth should have exactly one model (yard_sale)."""
        sw_models = [name for name, reg in MODEL_REGISTRY.items()
                     if reg.substrate_type == 'scalar_wealth']
        assert sw_models == ['yard_sale']


class TestSprint19Registrations:
    """Sprint 19: close the Sprint 11 registration gap.

    The Lotka-Volterra lattice model and P11 predator-prey oscillation
    detector shipped in Sprint 11 (commit 5f013c9) with direct-import
    tests in ``test_lv_p11_e2e.py`` but were never added to
    MODEL_REGISTRY / DETECTOR_REGISTRY. This block pins the Sprint 19
    registration so the gap cannot silently reopen.

    Scientific anchors (Sprint 11 findings, preserved for grep-ability):
    - LV occupies lattice_2d with the 'grid' observable (same substrate
      as RPS) and is substrate+observable-compatible with 6 detectors:
      P1, P11, P12, P13, P15, P22. Of these, P11 is canonical DEFINITIVE;
      P1 is confirmation (Sprint 11 §5); P12/P13/P15/P22 are
      substrate-compatible but correctly rejected at content level.
    - P11 is lattice_2d with 'grid' required; its bilateral-vs-cyclic
      discrimination vs P12 happens at the content-level
      n_unique_species_observed prerequisite (=2 for P11, ≥3 for P12),
      NOT at the substrate/observable guard. See Decisions 34–36.
    - P11 is substrate+observable-compatible with 7 models (all
      lattice_2d-with-grid); btw_sandpile is lattice_2d but lacks the
      'grid' observable and is correctly rejected at the observable
      guard (same pattern as P12).
    """

    def test_lotka_volterra_registered(self):
        assert 'lotka_volterra' in MODEL_REGISTRY, \
            "lotka_volterra model should be registered"
        m = MODEL_REGISTRY['lotka_volterra']
        assert m.substrate_type == 'lattice_2d'
        assert 'grid' in m.observables
        assert 'P11' in m.primary_patterns

    def test_lotka_volterra_species_metadata(self):
        """P11's prerequisite uses species identity; LV's n_species = 2
        metadata flag should be declared so downstream inspection /
        future species-hint auto-resolution can rely on it (Decision 34)."""
        m = MODEL_REGISTRY['lotka_volterra']
        for key in ['n_species', 'n_states', 'model_class', 'interaction_type']:
            assert key in m.metadata_keys, \
                f"lotka_volterra must declare metadata key '{key}' for P11"

    def test_p11_registered(self):
        assert 'P11' in DETECTOR_REGISTRY, "P11 detector should be registered"
        d = DETECTOR_REGISTRY['P11']
        assert 'lattice_2d' in d.required_substrate
        assert 'grid' in d.required_observables
        assert d.observable_scope == 'state_history_only'

    def test_lotka_volterra_p11_compatible(self):
        """The canonical Sprint 11 pair: lotka_volterra × P11 — now
        finally pinned in the orchestration layer."""
        r = check_compatibility('lotka_volterra', 'P11')
        assert r.compatible, \
            f"lotka_volterra × P11 should be compatible: {r.reason}"

    def test_p11_rejects_non_lattice_2d_substrates(self):
        """P11 should reject every non-lattice_2d substrate at
        substrate_mismatch. Parallel structure to TestSprint9Registrations'
        P12 test (P11 and P12 share substrate by design)."""
        for m in MODEL_REGISTRY:
            reg = MODEL_REGISTRY[m]
            r = check_compatibility(m, 'P11')
            if reg.substrate_type != 'lattice_2d':
                assert not r.compatible, \
                    f"{m} ({reg.substrate_type}) should NOT match P11"
                assert 'substrate_mismatch' in r.reason, \
                    f"{m} × P11 rejection reason should be substrate_mismatch"

    def test_p11_rejects_lattice_2d_without_grid(self):
        """btw_sandpile is lattice_2d but exposes only avalanche-level
        observables. P11 should reject it at the observable guard, not
        at substrate guard — exactly the pattern P12 exhibits."""
        r = check_compatibility('btw_sandpile', 'P11')
        assert not r.compatible
        assert 'missing_observable' in r.reason, \
            f"btw_sandpile × P11 should fail on observable, got: {r.reason}"

    def test_p11_compatible_with_all_lattice_2d_grid_models(self):
        """P11 is a lattice_2d detector requiring 'grid' observable. Every
        lattice_2d model exposing a 'grid' observable should pair at the
        orchestration layer. Discrimination vs positive-screening happens
        at content level in the detector, not at the registry gate.

        Sprint 11's negative-model tests in test_lv_p11_e2e.py exercise
        the content-level rejection for Schelling, Nowak-May, SIR, RPS.
        GH and GoL are registry-compatible but no content-level negative
        test exists yet — a carry-forward beyond this cleanup sprint.
        """
        p11_pairs = {m for m, d in get_compatible_pairs() if d == 'P11'}
        lattice_2d_grid_models = {
            name for name, reg in MODEL_REGISTRY.items()
            if reg.substrate_type == 'lattice_2d' and 'grid' in reg.observables
        }
        assert p11_pairs == lattice_2d_grid_models, (
            f"P11 pairs {p11_pairs} should equal lattice_2d-with-grid models "
            f"{lattice_2d_grid_models}"
        )

    def test_lotka_volterra_pairs_with_six_detectors(self):
        """LV is lattice_2d with 'grid': should pair with P1, P11, P12,
        P13, P15, P22 (six). P14 needs avalanche_sizes (rejected at
        observable); P27 needs coop_fraction (rejected at observable);
        all other detectors are on other substrates (rejected at
        substrate). This pins the full per-model column."""
        lv_dets = {d for m, d in get_compatible_pairs() if m == 'lotka_volterra'}
        assert lv_dets == {'P1', 'P11', 'P12', 'P13', 'P15', 'P22'}, \
            f"Expected LV to pair with 6 detectors, got {lv_dets}"

    def test_p11_and_p12_share_substrate_by_design(self):
        """P11 and P12 are the two predator-prey detectors: P11 for
        bilateral (2-species, LV), P12 for cyclic (3+ species, RPS).
        They MUST share substrate and observable requirements at the
        orchestration layer, because their bilateral-vs-cyclic
        discrimination is content-level (n_species prerequisite)."""
        p11 = DETECTOR_REGISTRY['P11']
        p12 = DETECTOR_REGISTRY['P12']
        assert p11.required_substrate == p12.required_substrate
        assert p11.required_observables == p12.required_observables
        # They should have identical compatible-model sets
        p11_pairs = {m for m, d in get_compatible_pairs() if d == 'P11'}
        p12_pairs = {m for m, d in get_compatible_pairs() if d == 'P12'}
        assert p11_pairs == p12_pairs, \
            f"P11 pairs {p11_pairs} should equal P12 pairs {p12_pairs}"
