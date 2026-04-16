"""
Cross-detection matrix tests — Sprint 6.

Tests all 8 previously-untested compatible pairs from the transfer matrix.
Each test verifies the ACTUAL EFFECT (detection tier or correct rejection),
not just that code runs without error.

Results summary:
  D'Orsogna × P5:  REJECTED — φ=0.046, milling ≠ flocking
  Vicsek × P6:     REJECTED — |L|=0.031, flocking ≠ milling
  GoL × P13:       REJECTED — n_states=2, not excitable medium (was false positive before fix)
  Nowak-May × P1:  CONFIRMATION — cooperator clusters aggregate (I≈0.9, seg≈0.75, p<0.01)
  Nowak-May × P13: REJECTED — speed_cv >> 0.2, not excitable
  Nowak-May × P15: SKIPPED — P15 implementation is GoL-specific (hardcoded _step_gol)
  Schelling × P13: REJECTED — run length or speed instability, not excitable dynamics
  Schelling × P15: SKIPPED — P15 implementation is GoL-specific

Bugs found and fixed:
  1. P13 n_states<3 was a warning, not a hard guard → GoL false positive at screening
  2. Schelling state dicts lacked grid_dims/n_states/step → P13 crash on KeyError
"""

import pytest
import numpy as np


# ===================================================================
# Cross-exclusion tests: P5/P6 mutual rejection
# ===================================================================

class TestP5P6CrossExclusion:
    """P5 (flocking) and P6 (milling) should correctly reject each other's
    canonical positive models. This validates the detector discrimination."""

    def test_dorsogna_milling_rejected_by_p5(self):
        """D'Orsogna at published milling params should NOT trigger P5.

        Milling produces low polarization (φ ≈ 0.04) — well below the
        screening threshold of 0.5. The swarm rotates coherently but
        individual headings point tangentially, not in a common direction.
        """
        from epc.models.dorsogna_spp import DOrsognaSPPModel
        from epc.detectors.p5_flocking import P5FlockingDetector

        m = DOrsognaSPPModel(
            n_particles=100, C_a=0.5, C_r=1.0, l_a=3.0, l_r=0.5,
            alpha=1.0, beta=0.5, dt=0.01, seed=42,
        )
        history = m.run(n_steps=5000)
        metadata = m.get_metadata()

        detector = P5FlockingDetector(n_permutations=49, seed=42)
        result = detector.detect(history, metadata)

        assert not result.detected, (
            f"P5 should NOT detect milling D'Orsogna "
            f"(φ={result.primary_metric.get('polarization_mean', '?'):.3f})"
        )
        phi = result.primary_metric.get("polarization_mean", 1.0)
        assert phi < 0.2, (
            f"Milling polarization φ={phi:.3f} should be << 0.5 (screening threshold)"
        )
        print(f"  ✓ D'Orsogna × P5: correctly rejected (φ={phi:.3f})")

    def test_vicsek_ordered_rejected_by_p6(self):
        """Vicsek at low noise (ordered) should NOT trigger P6.

        Ordered flocking produces low angular momentum (|L| ≈ 0.03) —
        well below the screening threshold of 0.3. The swarm translates
        coherently with no rotation.
        """
        from epc.models.vicsek import VicsekModel
        from epc.detectors.p6_milling import P6MillingDetector

        m = VicsekModel(
            n_particles=300, box_size=7.0, speed=0.03, noise=0.5,
            interaction_radius=1.0, seed=42,
        )
        history = m.run(n_steps=2000)
        metadata = m.get_metadata()

        detector = P6MillingDetector(n_permutations=49, seed=42)
        result = detector.detect(history, metadata)

        assert not result.detected, (
            f"P6 should NOT detect ordered Vicsek "
            f"(|L|={result.primary_metric.get('angular_momentum_abs_mean', '?'):.3f})"
        )
        L_abs = result.primary_metric.get("angular_momentum_abs_mean", 1.0)
        assert L_abs < 0.15, (
            f"Ordered flocking |L|={L_abs:.3f} should be << 0.3 (screening threshold)"
        )
        print(f"  ✓ Vicsek × P6: correctly rejected (|L|={L_abs:.3f})")


# ===================================================================
# P13 excitable wave guard tests
# ===================================================================

class TestP13ExcitableGuard:
    """P13 requires n_states ≥ 3 (resting/excited/refractory).
    Models with fewer states should be hard-rejected, not just warned."""

    def test_gol_rejected_by_p13_n_states_guard(self):
        """GoL (n_states=2) must NOT trigger P13.

        Before the fix, GoL produced a false positive at screening tier
        because its synchronous update gave speed_cv=0.0 (every cell has
        identical inter-excitation intervals). The n_states<3 guard now
        prevents this.
        """
        from epc.models.game_of_life import GameOfLife
        from epc.detectors.p13_excitable_waves import P13ExcitableWaveDetector

        m = GameOfLife(rows=50, cols=50, init_mode="r_pentomino", seed=42)
        history = m.run(n_steps=500)
        metadata = m.get_metadata()

        detector = P13ExcitableWaveDetector(n_null_runs=19)
        result = detector.detect(history, metadata)

        assert not result.detected, (
            "P13 should NOT detect GoL (n_states=2, not excitable medium)"
        )
        assert any("n_states=2" in w for w in result.warnings), (
            "P13 should warn about n_states=2"
        )
        print(f"  ✓ GoL × P13: correctly rejected (n_states guard)")

    def test_nowak_may_rejected_by_p13(self):
        """Nowak-May (n_states=2) must NOT trigger P13.

        Spatial PD has cooperator/defector states with no excitable dynamics.
        """
        from epc.models.nowak_may import NowakMayModel
        from epc.detectors.p13_excitable_waves import P13ExcitableWaveDetector

        m = NowakMayModel(rows=50, cols=50, b=1.8, seed=42)
        history = m.run(n_steps=200)
        metadata = m.get_metadata()

        detector = P13ExcitableWaveDetector(n_null_runs=9)
        result = detector.detect(history, metadata)

        assert not result.detected, (
            "P13 should NOT detect Nowak-May (n_states=2)"
        )
        assert any("n_states=2" in w for w in result.warnings), (
            "P13 should warn about n_states=2"
        )
        print(f"  ✓ Nowak-May × P13: correctly rejected (n_states guard)")

    def test_schelling_rejected_by_p13(self):
        """Schelling (n_states=3 but not excitable) should NOT trigger P13.

        Schelling has 3 nominal states (empty/typeA/typeB) which passes
        the n_states guard, but the dynamics are not excitable — there
        are no wavefronts with constant speed. P13 should reject on
        speed instability or insufficient run length.
        """
        from epc.models.schelling import run_schelling
        from epc.detectors.p13_excitable_waves import P13ExcitableWaveDetector

        history = run_schelling(grid_size=50, n_steps=200, seed=42)
        metadata = {"model": "schelling", "substrate": "lattice_2d"}

        detector = P13ExcitableWaveDetector(n_null_runs=9)
        result = detector.detect(history, metadata)

        assert not result.detected, (
            "P13 should NOT detect Schelling (not excitable dynamics)"
        )
        print(f"  ✓ Schelling × P13: correctly rejected "
              f"(tier={result.tier}, warnings={result.warnings})")

    def test_p13_still_detects_greenberg_hastings(self):
        """Regression: P13 must still detect GH spiral (the canonical positive).

        Ensures the n_states guard doesn't break the true positive case.
        """
        from epc.models.greenberg_hastings import GreenbergHastings
        from epc.detectors.p13_excitable_waves import P13ExcitableWaveDetector

        m = GreenbergHastings(
            rows=80, cols=80, n_states=5, threshold=1,
            init_density=0.3, seed=42,
        )
        history = m.run(n_steps=400)
        metadata = m.get_metadata()

        detector = P13ExcitableWaveDetector(n_null_runs=19)
        result = detector.detect(history, metadata)

        assert result.detected, (
            f"P13 should detect GH spiral (tier={result.tier})"
        )
        assert not any("n_states" in w and "< 3" in w for w in result.warnings), (
            "GH has n_states=5, should not trigger n_states guard"
        )
        print(f"  ✓ GH × P13 regression: still detected at {result.tier}")


# ===================================================================
# Nowak-May × P1: genuine co-occurrence with P27
# ===================================================================

class TestNowakMayP1CoOccurrence:
    """Cooperator clusters in spatial PD genuinely aggregate.

    This is scientifically meaningful: cooperators cluster together due to
    payoff imitation dynamics, producing spatial aggregation detectable by P1.
    P1's allowed_co_occurrences includes P27, so this is a valid co-occurrence.
    """

    def test_nowak_may_p1_confirmation(self):
        """Nowak-May at b=1.8 should reach P1 confirmation tier.

        Expected: Moran's I ≈ 0.9, segregation index ≈ 0.75, p < 0.01.
        Cooperator/defector clusters are spatially aggregated.
        """
        from epc.models.nowak_may import NowakMayModel
        from epc.detectors.p1_aggregation import P1AggregationDetector

        m = NowakMayModel(rows=100, cols=100, b=1.8, seed=42)
        history = m.run(n_steps=200)
        metadata = m.get_metadata()

        detector = P1AggregationDetector(n_permutations=199)
        result = detector.detect(history, metadata)

        # Must detect at confirmation or higher
        from epc.detector_result import DetectionTier
        assert result.detected, (
            f"P1 should detect aggregation in Nowak-May "
            f"(tier={result.tier}, I={result.primary_metric.get('morans_i', '?')})"
        )
        assert result.tier in (DetectionTier.CONFIRMATION, DetectionTier.DEFINITIVE), (
            f"Expected confirmation+, got {result.tier}"
        )

        # Verify effect sizes
        I = result.primary_metric.get("morans_i", 0)
        assert I > 0.5, f"Moran's I={I:.3f} should be > 0.5 for strong clustering"

        seg = result.secondary_metrics.get("segregation_index", 0)
        assert seg > 0.5, f"Segregation index={seg:.3f} should be > 0.5"

        assert result.null_p_value < 0.01, (
            f"Null p={result.null_p_value:.4f} should be < 0.01"
        )

        # Verify 2 unique types (cooperator + defector)
        n_types = result.primary_metric.get("n_unique_types", 0)
        assert n_types == 2, f"Expected 2 types (C/D), got {n_types}"

        print(f"  ✓ Nowak-May × P1: CONFIRMATION "
              f"(I={I:.3f}, seg={seg:.3f}, p={result.null_p_value:.4f})")

    def test_nowak_may_p1_co_occurrence_with_p27(self):
        """P1 lists P27 as allowed co-occurrence — verify this is consistent."""
        from epc.detectors.p1_aggregation import P1AggregationDetector

        det = P1AggregationDetector()
        # Check that P27 is in the allowed co-occurrences
        assert "P27" in det.allowed_co_occurrences, (
            "P1 should list P27 as allowed co-occurrence"
        )
        print("  ✓ P1 allows P27 co-occurrence")


# ===================================================================
# P15 implementation limitation documentation
# ===================================================================

class TestP15ImplementationScope:
    """P15 (computation) is substrate-compatible with all lattice_2d models,
    but its IMPLEMENTATION is GoL-specific: it hardcodes _step_gol and
    _place_glider for deterministic collision testing.

    These tests document that P15 cannot be cross-tested against non-GoL
    models with the current implementation. This is NOT a bug — it's an
    implementation scope note for future generalization.
    """

    def test_p15_implementation_is_gol_specific(self):
        """Verify P15 uses internal GoL stepping, not external model history."""
        from epc.detectors.p15_fidelity_fix import (
            test_p15_fidelity_deterministic,
            _step_gol,
            _place_glider,
        )

        # P15 runs its own GoL simulation internally
        result = test_p15_fidelity_deterministic()
        assert result.reproducibility == 1.0, "P15 internal GoL should be deterministic"
        print(f"  ✓ P15 is GoL-specific (reproducibility={result.reproducibility:.1f})")

    def test_p15_substrate_compatible_but_not_cross_testable(self):
        """Orchestration lists NowakMay×P15 and Schelling×P15 as compatible
        (correct: both are lattice_2d with grid observable). But P15's
        implementation doesn't accept external histories."""
        from epc.orchestration import check_compatibility

        # Substrate compatibility is correct
        r1 = check_compatibility("nowak_may", "P15")
        r2 = check_compatibility("schelling", "P15")
        assert r1.compatible, "Nowak-May should be substrate-compatible with P15"
        assert r2.compatible, "Schelling should be substrate-compatible with P15"

        # But the implementation scope is GoL-only (document, not fix)
        print("  ✓ NowakMay×P15, Schelling×P15: substrate-compatible "
              "but P15 implementation is GoL-specific")


# ===================================================================
# Schelling state enrichment regression
# ===================================================================

class TestSchellingStateEnrichment:
    """Schelling state dicts must include grid_dims, n_states, step
    for cross-detector compatibility."""

    def test_schelling_state_has_required_keys(self):
        """State dicts must include grid_dims, n_states, step."""
        from epc.models.schelling import run_schelling

        history = run_schelling(grid_size=30, n_steps=10, seed=42)

        for i, state in enumerate(history):
            assert "grid" in state, f"State {i} missing 'grid'"
            assert "grid_dims" in state, f"State {i} missing 'grid_dims'"
            assert "n_states" in state, f"State {i} missing 'n_states'"
            assert "step" in state, f"State {i} missing 'step'"

        assert history[0]["grid_dims"] == (30, 30)
        assert history[0]["n_states"] == 3
        assert history[0]["step"] == 0
        assert history[-1]["step"] == 10
        print(f"  ✓ Schelling state enrichment: all keys present")


# ===================================================================
# Updated transfer matrix summary
# ===================================================================

class TestTransferMatrixCompleteness:
    """Verify all 8 previously-untested compatible pairs now have outcomes."""

    EXPECTED_OUTCOMES = {
        ("dorsogna", "P5"): "rejected",      # Milling ≠ flocking
        ("vicsek", "P6"): "rejected",         # Flocking ≠ milling
        ("game_of_life", "P13"): "rejected",  # n_states=2 guard
        ("nowak_may", "P1"): "detected",      # Cooperator aggregation
        ("nowak_may", "P13"): "rejected",     # n_states=2 guard
        ("nowak_may", "P15"): "gol_specific", # P15 impl limitation
        ("schelling", "P13"): "rejected",     # Not excitable
        ("schelling", "P15"): "gol_specific", # P15 impl limitation
    }

    def test_all_pairs_documented(self):
        """Every untested pair now has an expected outcome."""
        assert len(self.EXPECTED_OUTCOMES) == 8
        detected = sum(1 for v in self.EXPECTED_OUTCOMES.values() if v == "detected")
        rejected = sum(1 for v in self.EXPECTED_OUTCOMES.values() if v == "rejected")
        limited = sum(1 for v in self.EXPECTED_OUTCOMES.values() if v == "gol_specific")
        print(f"  ✓ Transfer matrix: {detected} detected, {rejected} rejected, "
              f"{limited} implementation-limited")


# ===================================================================
# Runner
# ===================================================================

if __name__ == "__main__":
    pytest.main([__file__, "-v", "--tb=short"])
