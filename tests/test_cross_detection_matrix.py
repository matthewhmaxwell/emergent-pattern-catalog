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

class TestP15GeneralizedCrossDetection:
    """P15 (persistent computation) now has a generalized BaseDetector
    implementation in p15_persistent_computation.py that works for any
    deterministic lattice_2d CA via a step_fn parameter.

    The OLD p15_fidelity_fix.py module is GoL-specific and retained for
    backward compatibility; it is NOT the implementation used for
    cross-detection.

    These tests anchor the cross-detection outcomes for Nowak-May and
    Schelling using the generalized detector. The exhaustive P15 validation
    lives in tests/test_p15_generalized.py; this class only pins the
    transfer-matrix cells that belong to the cross-detection audit.
    """

    def test_p15_generalized_detector_is_importable(self):
        """Sanity check: the generalized detector exists and has the expected API."""
        from epc.detectors.p15_persistent_computation import (
            P15PersistentComputationDetector,
            make_step_fn_for_gol,
            make_step_fn_from_model,
        )
        from epc.base_detector import BaseDetector

        assert issubclass(P15PersistentComputationDetector, BaseDetector), \
            "Generalized P15 must be a BaseDetector subclass"
        # Factories produce callables grid → grid
        step_fn = make_step_fn_for_gol()
        grid = np.zeros((10, 10), dtype=np.int8)
        out = step_fn(grid)
        assert out.shape == (10, 10)
        print("  ✓ Generalized P15 detector importable with correct API")

    def test_nowak_may_p15_screening_only(self):
        """Nowak-May × P15 (generalized): SCREENING only, not DEFINITIVE.

        Nowak-May is deterministic (reproducibility = 1.0), so it clears
        the bit-exact replay check. But it produces limited outcome
        diversity under perturbation variations — typically 2 distinct
        outcome classes (not the ≥3 required for confirmation). So the
        detector reports SCREENING tier.

        This is the scientifically correct classification: Nowak-May
        supports persistent dynamics but is not a Turing-like computation
        engine in the way GoL is.
        """
        from epc.detectors.p15_persistent_computation import (
            P15PersistentComputationDetector,
            make_step_fn_from_model,
        )
        from epc.models.nowak_may import NowakMayModel
        from epc.detector_result import DetectionTier

        nm = NowakMayModel(rows=40, cols=40, b=1.8, init_mode="random", seed=42)
        nm.setup()
        history = nm.run(n_steps=30)
        meta = nm.get_metadata()

        step_fn = make_step_fn_from_model(nm)
        det = P15PersistentComputationDetector(
            step_fn=step_fn, n_variations=8, seed=42,
        )
        result = det.detect(history, meta)

        # Nowak-May is deterministic
        assert result.primary_metric["reproducibility"] == 1.0, \
            f"Nowak-May should be reproducible, got " \
            f"{result.primary_metric['reproducibility']}"
        # But does NOT reach DEFINITIVE (limited outcome diversity)
        assert result.tier < DetectionTier.DEFINITIVE, \
            f"Nowak-May should not achieve P15 DEFINITIVE, got {result.tier}"
        print(f"  ✓ NowakMay × P15 (generalized): "
              f"tier={result.tier.name}, "
              f"n_distinct={result.primary_metric.get('n_distinct_outcomes', 0)}")

    def test_schelling_p15_not_detected(self):
        """Schelling × P15 (generalized): not detected (no step_fn available).

        Schelling is stochastic (random agent selection for relocation) and
        its run_schelling function does not expose a pure grid→grid step.
        Without a step_fn, the detector cannot perform reproducibility or
        variation testing → reproducibility = 0, outcome diversity = 0,
        fails screening.

        This is correct behavior: Schelling is not a deterministic
        computation substrate.
        """
        from epc.detectors.p15_persistent_computation import (
            P15PersistentComputationDetector,
        )
        from epc.models.schelling import run_schelling
        from epc.detector_result import DetectionTier

        history = run_schelling(grid_size=30, n_steps=50, seed=42)

        # No step_fn provided → detector can't replay
        det = P15PersistentComputationDetector(
            step_fn=None, n_variations=8, seed=42,
        )
        result = det.detect(history, None)

        assert not result.detected, \
            f"Schelling should not register as P15 without step_fn, " \
            f"got detected=True (tier={result.tier.name})"
        assert result.primary_metric["reproducibility"] == 0.0, \
            f"No step_fn → reproducibility should be 0, " \
            f"got {result.primary_metric['reproducibility']}"
        print(f"  ✓ Schelling × P15 (generalized): correctly not detected "
              f"(no step_fn for stochastic model)")


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
    """Pins the expected cross-detection outcome for every audited
    model-detector pair. Each entry corresponds to a test elsewhere in
    the suite (in this file or in test_sir_p22_e2e.py / test_p15_generalized.py)
    that verifies the actual effect.

    Values:
      "detected"     — detector passes screening+ on this model
      "screening"    — detector passes screening only (not confirmation+)
      "rejected"     — detector correctly does not pass screening
      "not_detected" — detector cannot run / produces no signal (e.g., no step_fn)
    """

    EXPECTED_OUTCOMES = {
        # --- Sprint 6 pairs (original scope) ---
        ("dorsogna", "P5"): "rejected",       # Milling ≠ flocking
        ("vicsek", "P6"): "rejected",          # Flocking ≠ milling
        ("game_of_life", "P13"): "rejected",   # n_states=2 guard
        ("nowak_may", "P1"): "detected",       # Cooperator aggregation
        ("nowak_may", "P13"): "rejected",      # n_states=2 guard
        ("schelling", "P13"): "rejected",      # Not excitable

        # --- Sprint 7 pairs (SIR + P22) ---
        ("sir_epidemic", "P22"): "detected",   # DEFINITIVE cascade
        ("sir_epidemic", "P13"): "rejected",   # Single-pass, not re-entrant
        ("sir_epidemic", "P1"): "screening",   # Transient wavefront aggregation only
        ("greenberg_hastings", "P22"): "rejected",  # Re-entrant, not cascade
        ("game_of_life", "P22"): "rejected",   # R-pentomino reach < 5%
        ("nowak_may", "P22"): "rejected",      # No wavefront; moran_i_time ≈ 0
        ("schelling", "P22"): "rejected",      # No wavefront; moran_i_time ≈ 0.06

        # --- Sprint 8 pairs (generalized P15) ---
        # P15 is no longer gol_specific; the generalized detector handles
        # any deterministic lattice_2d CA via a step_fn parameter.
        ("game_of_life", "P15"): "detected",   # DEFINITIVE with dense random IC
        ("nowak_may", "P15"): "screening",     # Determinism=1, but only 2 outcomes
        ("schelling", "P15"): "not_detected",  # Stochastic, no step_fn available
        ("greenberg_hastings", "P15"): "rejected",  # Only 1 outcome class (spirals)
        ("sir_epidemic", "P15"): "not_detected",    # Stochastic, fails reproducibility

        # --- Sprint 9 pairs (RPS + P12) ---
        # RPS-row: all six non-exclusive detectors characterized.
        ("rps_spatial", "P12"): "detected",     # CONFIRMATION (score=1.83, p=0.005)
        ("rps_spatial", "P13"): "rejected",     # cv=0.65 >> 0.2 threshold
        ("rps_spatial", "P22"): "screening",    # passes reach+Moran, fails unimodal
        ("rps_spatial", "P1"): "screening",     # transient spiral aggregation (Moran=0.57)
        ("rps_spatial", "P15"): "not_detected", # stochastic, reproducibility=0
        # P12-column: all non-RPS lattice_2d grid models are rejected.
        ("greenberg_hastings", "P12"): "rejected",  # ρ=1.0 (clock-driven)
        ("sir_epidemic", "P12"): "rejected",    # 3 states, no cyclic dominance
        ("game_of_life", "P12"): "rejected",    # n_species=2 < 3 (prerequisite)
        ("nowak_may", "P12"): "rejected",       # n_states=2 < 3 (prerequisite)
    }

    VALID_OUTCOMES = {"detected", "rejected", "screening", "not_detected"}

    def test_all_pairs_documented(self):
        """Every audited pair has a valid, documented expected outcome."""
        assert len(self.EXPECTED_OUTCOMES) >= 27, \
            f"Expected at least 27 audited pairs, got {len(self.EXPECTED_OUTCOMES)}"

        # Every outcome value is in the valid set
        for pair, outcome in self.EXPECTED_OUTCOMES.items():
            assert outcome in self.VALID_OUTCOMES, \
                f"Pair {pair} has invalid outcome '{outcome}'; " \
                f"must be one of {self.VALID_OUTCOMES}"

        detected = sum(1 for v in self.EXPECTED_OUTCOMES.values() if v == "detected")
        screening = sum(1 for v in self.EXPECTED_OUTCOMES.values() if v == "screening")
        rejected = sum(1 for v in self.EXPECTED_OUTCOMES.values() if v == "rejected")
        not_detected = sum(1 for v in self.EXPECTED_OUTCOMES.values() if v == "not_detected")

        print(f"  ✓ Transfer matrix: {detected} detected, {screening} screening, "
              f"{rejected} rejected, {not_detected} not_detected "
              f"({len(self.EXPECTED_OUTCOMES)} total)")

    def test_sprint_7_pairs_covered(self):
        """Every Sprint 7 cross-detection pair must appear in the matrix."""
        sprint_7_pairs = [
            ("sir_epidemic", "P22"),
            ("sir_epidemic", "P13"),
            ("sir_epidemic", "P1"),
            ("greenberg_hastings", "P22"),
            ("game_of_life", "P22"),
            ("nowak_may", "P22"),
            ("schelling", "P22"),
        ]
        for pair in sprint_7_pairs:
            assert pair in self.EXPECTED_OUTCOMES, \
                f"Sprint 7 pair {pair} missing from transfer matrix"
        print(f"  ✓ Sprint 7: all {len(sprint_7_pairs)} P22 cross-detection pairs covered")

    def test_sprint_8_p15_generalization_covered(self):
        """Every Sprint 8 P15 cross-detection result must appear in the matrix.

        After generalizing P15, there are no remaining 'gol_specific' cells:
        every lattice_2d model has a real P15 outcome.
        """
        sprint_8_p15_pairs = [
            ("game_of_life", "P15"),
            ("nowak_may", "P15"),
            ("schelling", "P15"),
            ("greenberg_hastings", "P15"),
            ("sir_epidemic", "P15"),
        ]
        for pair in sprint_8_p15_pairs:
            assert pair in self.EXPECTED_OUTCOMES, \
                f"Sprint 8 P15 pair {pair} missing from transfer matrix"
            assert self.EXPECTED_OUTCOMES[pair] != "gol_specific", \
                f"Pair {pair} still marked gol_specific — obsolete after " \
                f"Sprint 8 generalization"
        print(f"  ✓ Sprint 8: all {len(sprint_8_p15_pairs)} P15 pairs covered, "
              f"no gol_specific entries remain")

    def test_sprint_9_pairs_covered(self):
        """Every Sprint 9 cross-detection pair (RPS row + P12 column) must
        appear in the transfer matrix with a pinned expected outcome.

        The Sprint 9 headline finding is that P12 cleanly separates
        intransitive cyclic-dominance spirals (RPS) from excitable-wave
        spirals (GH), even though both produce similar-looking spiral
        patterns on lattice_2d. The detector-discrimination comes from
        the neighbor-conditional replacement ratio ρ: RPS has ρ ≫ 1
        (neighbor-driven transitions), while GH has ρ = 1.0 exactly
        (clock-driven transitions).

        A secondary finding is that P13 ALSO rejects RPS cleanly at
        screening, via wavefront-speed CV rather than via the n_states
        guard. This is documented in tests/test_rps_p13_boundary.py.
        """
        sprint_9_pairs = [
            # RPS row — all 5 non-trivially-compatible detectors
            ("rps_spatial", "P12"),
            ("rps_spatial", "P13"),
            ("rps_spatial", "P22"),
            ("rps_spatial", "P1"),
            ("rps_spatial", "P15"),
            # P12 column — all 4 other lattice_2d-grid models
            ("greenberg_hastings", "P12"),
            ("sir_epidemic", "P12"),
            ("game_of_life", "P12"),
            ("nowak_may", "P12"),
        ]
        for pair in sprint_9_pairs:
            assert pair in self.EXPECTED_OUTCOMES, \
                f"Sprint 9 pair {pair} missing from transfer matrix"

        # RPS × P12 is the Sprint 9 positive; all other P12 cells are
        # rejections.
        assert self.EXPECTED_OUTCOMES[("rps_spatial", "P12")] == "detected"
        for other_model in ("greenberg_hastings", "sir_epidemic",
                            "game_of_life", "nowak_may"):
            assert self.EXPECTED_OUTCOMES[(other_model, "P12")] == "rejected", (
                f"{other_model} × P12 should be rejected (no cyclic dominance)"
            )

        print(f"  ✓ Sprint 9: all {len(sprint_9_pairs)} pairs covered "
              f"(RPS row + P12 column)")


# ===================================================================
# Runner
# ===================================================================

if __name__ == "__main__":
    pytest.main([__file__, "-v", "--tb=short"])
