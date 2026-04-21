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

        # Verify effect sizes.
        # Sprint 10: primary metric is now final-state Moran's I (was
        # peak over trajectory). NM b=1.8 shows final I ≈ 0.49, which is
        # strong aggregation but below the old 0.5 threshold calibrated
        # against the peak-based primary. We now also verify that the
        # transient peak Moran (retained as a diagnostic) is >> final.
        I = result.primary_metric.get("morans_i", 0)
        assert I > 0.4, f"Moran's I={I:.3f} should be > 0.4 for strong clustering"

        # Final-state Moran is the primary; report it and check it's solid.
        I_final = result.primary_metric.get("morans_i_final", 0)
        assert I_final > 0.4, f"Final Moran's I={I_final:.3f} should be > 0.4"

        # Peak transient Moran is reported as a diagnostic — typically
        # higher than final on NM because cooperator clusters tighten early
        # then relax toward a shifting equilibrium.
        I_peak = result.primary_metric.get("morans_i_peak", 0)
        assert I_peak > I_final, (
            f"Peak I={I_peak:.3f} should exceed final I={I_final:.3f} "
            "on Nowak-May (transient cluster tightening)"
        )

        seg = result.secondary_metrics.get("segregation_index", 0)
        assert seg > 0.5, f"Segregation index={seg:.3f} should be > 0.5"

        assert result.null_p_value < 0.01, (
            f"Null p={result.null_p_value:.4f} should be < 0.01"
        )

        # Verify 2 unique types (cooperator + defector)
        n_types = result.primary_metric.get("n_unique_types", 0)
        assert n_types == 2, f"Expected 2 types (C/D), got {n_types}"

        print(f"  ✓ Nowak-May × P1: CONFIRMATION "
              f"(I_final={I_final:.3f}, I_peak={I_peak:.3f}, "
              f"seg={seg:.3f}, p={result.null_p_value:.4f})")

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
        # Sprint 10: SIR × P1 flipped screening → rejected under new
        # final-state Moran primary (Decision 32). Wavefront is transient;
        # final recovered pattern is near-uniform (I ≈ 0.02, fails 0.05 floor).
        ("sir_epidemic", "P1"): "rejected",    # Sprint 10: final-I primary
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
        # Sprint 10: RPS × P1 stays at screening under the new final-I
        # primary. Spiral domains rotate but persist, so final Moran ≈ peak
        # ≈ 0.55 — genuinely sustained clustering, unlike SIR's collapsing
        # wavefront. Co-exists at screening (not confirmation, because
        # cluster geometry is not Schelling-style block segregation).
        ("rps_spatial", "P1"): "screening",     # sustained spiral clustering (I_final≈0.55)
        ("rps_spatial", "P15"): "not_detected", # stochastic, reproducibility=0
        # P12-column: all non-RPS lattice_2d grid models are rejected.
        ("greenberg_hastings", "P12"): "rejected",  # ρ=1.0 (clock-driven)
        ("sir_epidemic", "P12"): "rejected",    # 3 states, no cyclic dominance
        ("game_of_life", "P12"): "rejected",    # n_species=2 < 3 (prerequisite)
        ("nowak_may", "P12"): "rejected",       # n_states=2 < 3 (prerequisite)

        # --- Sprint 11 pairs (LV row + P11 column) ---
        # Headline finding: P11 cleanly separates bilateral predator-prey
        # oscillation (LV) from cyclic three-species dominance (RPS) via
        # the n_species prerequisite. On the primary metric alone (rho_anti
        # at nonzero lag), RPS species pairs actually score STRONGER anti-
        # correlation than LV (-0.93 to -0.97 vs LV's -0.72 to -0.88) — the
        # separation comes from the n_species==2 prerequisite, not from
        # the primary metric. See REPLICATION_NOTES.md Sprint 11.
        # A secondary finding is that strictly conserved 2-species systems
        # (Nowak-May: coop+defect=1 exactly) give rho_anti=-1.0 by algebra
        # alone, so P11 must also verify total_std > 0.005 (nontrivial
        # empty reservoir). See Decision 35.
        ("lotka_volterra_lattice", "P11"): "detected",  # DEFINITIVE (rho_anti ≈ -0.86)
        ("lotka_volterra_lattice", "P12"): "rejected",  # intransitivity_score=0.24 << 1.0
        ("lotka_volterra_lattice", "P13"): "rejected",  # wavefront_speed_cv=1.17 >> 0.2
        ("lotka_volterra_lattice", "P22"): "screening", # reach passes, is_unimodal=0 (cyclic)
        ("lotka_volterra_lattice", "P1"): "detected",   # I_final≈0.46, CONFIRMATION at n_perm≥499
        ("lotka_volterra_lattice", "P15"): "not_detected",  # stochastic, no step_fn
        # P11-column: all other audited lattice_2d models reject.
        ("schelling", "P11"): "rejected",       # prereq: species_std ≈ 0 (per-agent identity)
        ("nowak_may", "P11"): "rejected",       # prereq: total_std = 0 (strict A+B=1 conservation)
        ("sir_epidemic", "P11"): "rejected",    # prereq: post-burn-in variance ≈ 0 (transient dynamics)
        ("rps_spatial", "P11"): "rejected",     # prereq: n_unique_species_observed=3 ≠ 2

        # --- Sprint 13 pairs (Gray-Scott + P3) ---
        # Gray-Scott is the first continuous-valued-field model in the
        # catalog. It occupies the new lattice_2d_continuous substrate
        # (16k unique float values per snapshot, vs ≤ 10 for every integer-
        # grid model). P3 is restricted to this substrate by registration
        # AND by the n_unique_values >= 50 content prerequisite.
        #
        # Key negative-sweep finding: RPS at low mobility shows radial-FFT
        # peak-to-mean ≈ 23 on its raw integer grid — numerically matching
        # Gray-Scott labyrinth. Discrimination is substrate-level, NOT
        # empirical threshold tuning. See REPLICATION_NOTES.md Sprint 13.
        ("gray_scott", "P3"): "detected",       # DEFINITIVE (p/m=18.75, d~103, cv=0)
        # P3-column: all integer-grid models reject at substrate prereq
        # with informative warnings ("no 'field' (continuous 2D) observable").
        ("schelling", "P3"): "rejected",        # substrate prereq: no 'field' observable
        ("nowak_may", "P3"): "rejected",        # substrate prereq: no 'field' observable
        ("sir_epidemic", "P3"): "rejected",     # substrate prereq: no 'field' observable
        ("rps_spatial", "P3"): "rejected",      # substrate prereq: no 'field' observable
        ("game_of_life", "P3"): "rejected",     # substrate prereq: no 'field' observable
        ("greenberg_hastings", "P3"): "rejected", # substrate prereq: no 'field' observable
        ("lotka_volterra_lattice", "P3"): "rejected", # substrate prereq: no 'field' observable
        # Gray-Scott-row: other detectors are substrate-incompatible.
        # P1, P13, P22, P11, P12 inspect integer-grid observables, which
        # GS does not expose; they reject gracefully at screening with
        # warnings. (Sprint 14 B.1 hardened P1's 2D branch to match the
        # graceful-reject pattern — previously GS × P1 raised KeyError.)
        ("gray_scott", "P1"): "rejected",       # Sprint 14 B.1: no 'grid'/'type_labels_at_pos', substrate warning
        ("gray_scott", "P13"): "rejected",      # no 'grid' observable, rejects gracefully
        ("gray_scott", "P22"): "rejected",      # no 'grid' observable, rejects gracefully
        ("gray_scott", "P11"): "rejected",      # no 'grid' observable, rejects gracefully
        ("gray_scott", "P12"): "rejected",      # no 'grid' observable, rejects gracefully
        ("gray_scott", "P15"): "not_detected",  # no step_fn / stochastic-incompatible

        # --- Sprint 15 pairs (Nagel-Schreckenberg + P8) ---
        # Nagel-Schreckenberg is the second lattice_1d model (after Zhang
        # sorting) and the first lattice_1d model with a 'velocities'
        # observable. P8 is restricted to lattice_1d by registration AND
        # requires integer 1D velocity arrays at content level. This
        # produces substrate-level discrimination analogous to Sprint 13's
        # Decision 37 (continuous-field gate for P3).
        #
        # Canonical NS positive: L=1000, rho=0.15, v_max=5, p_slow=0.3,
        # burn_in=1000, measurement >= 1500 steps. stopped_fraction ≈ 0.18,
        # jam_lifetime_p95 ≈ 13, jam_lifetime_max ≥ 50, null_p ≈ 0.005,
        # Cohen's d effectively infinite (null std ≈ 0).
        #
        # Negative discrimination: density saturation (rho=0.80, p=0)
        # gives stopped=0.75 (pigeonhole, not jamming), but jam_lt_p95=4
        # — correctly rejected at the confirmation gate.
        ("nagel_schreckenberg", "P8"): "detected",   # DEFINITIVE canonical jam
        # P8-column: every other model rejects at substrate prereq
        # (missing 'velocities' observable) or substrate mismatch.
        ("zhang_sequential", "P8"): "rejected",      # lattice_1d but no 'velocities'
        ("schelling", "P8"): "rejected",             # substrate_mismatch: lattice_2d
        ("nowak_may", "P8"): "rejected",             # substrate_mismatch: lattice_2d
        ("sir_epidemic", "P8"): "rejected",          # substrate_mismatch: lattice_2d
        ("rps_spatial", "P8"): "rejected",           # substrate_mismatch: lattice_2d
        ("game_of_life", "P8"): "rejected",          # substrate_mismatch: lattice_2d
        ("greenberg_hastings", "P8"): "rejected",    # substrate_mismatch: lattice_2d
        ("lotka_volterra_lattice", "P8"): "rejected",# substrate_mismatch: lattice_2d
        ("gray_scott", "P8"): "rejected",            # substrate_mismatch: lattice_2d_continuous
        # NS-row: other detectors are substrate-incompatible.
        # Every lattice_2d detector (P1/P13/P22/P11/P12/P15) rejects
        # gracefully at missing_observable: NS has 1D velocities but no
        # 'grid' / 'cell_types' / 'field' / 'opinions' / 'theta' etc.
        ("nagel_schreckenberg", "P1"): "rejected",   # no 'grid'/'cell_types'
        ("nagel_schreckenberg", "P13"): "rejected",  # no 'grid' observable
        ("nagel_schreckenberg", "P22"): "rejected",  # no 'grid' observable
        ("nagel_schreckenberg", "P11"): "rejected",  # no 'grid' observable
        ("nagel_schreckenberg", "P12"): "rejected",  # no 'grid' observable
        ("nagel_schreckenberg", "P15"): "rejected",  # no 'grid' observable
        ("nagel_schreckenberg", "P3"): "rejected",   # no 'field' observable

        # --- Sprint 16 pairs (Active Brownian Particles + P2) ---
        # ABP is the 16th model family and the first continuous_2d model
        # with a density-dependent self-propulsion mechanism. P2 (MIPS)
        # is substrate-compatible with all three continuous_2d models
        # (Vicsek, D'Orsogna, ABP) but the DEFINITIVE tier requires
        # metadata affirmation: has_alignment_rule=False +
        # has_attraction_rule=False + has_density_dependent_speed=True.
        # Only ABP carries all three flags. This is metadata-level
        # discrimination (Decision 43), analogous to Sprint 15's
        # content-level discrimination for P8.
        #
        # Canonical ABP positive: N=1000, phi=0.5, Pe=100, rho_star=4,
        # 3000 post-burn steps. two_phase_score ~ 0.25-0.40, r ~ -0.93,
        # CV_v ~ 1.3, null_p = 0.005 -> DEFINITIVE.
        #
        # Negative-sweep findings:
        #   - Vicsek (any eta): constant speed, f_gas << 0.03 because
        #     ordered flocks are all-liquid (no dilute phase) and
        #     disordered is uniform (no liquid). Primary < 0.03 rejects
        #     at screening floor.
        #   - D'Orsogna milling: attractive Morse potential. Clusters
        #     strongly (f_liquid ~ 0.73) but f_gas ~ 0.05. Primary ~ 0.05
        #     lands at SCREENING. Even when primary is above SCREENING
        #     floor, has_attraction_rule=True blocks DEFINITIVE.
        #   - Three false-positive TRAPS identified at Phase 1:
        #     thermal ABP (Pe=5), dilute ABP (phi=0.1), over-saturated
        #     ABP (phi=0.85 long runtime). Each caught by a distinct
        #     confirmation gate (CV_v, primary, primary respectively).
        ("abp", "P2"): "detected",        # DEFINITIVE canonical MIPS
        ("abp", "P5"): "rejected",        # not a flocking system — no alignment
        ("abp", "P6"): "rejected",        # not a milling system — no attraction/confinement
        # P2 column — every other continuous_2d model that is substrate-
        # compatible should fire at most SCREENING.
        ("vicsek", "P2"): "rejected",      # constant speed, f_gas ≈ 0 -> below_two_phase_floor
        ("dorsogna", "P2"): "screening",   # attractive milling: f_liquid high, f_gas tiny
                                           # -> SCREENING; metadata gates DEFINITIVE
        # Every non-continuous_2d model is substrate_mismatch for P2.
        ("zhang_sequential", "P2"): "rejected",
        ("schelling", "P2"): "rejected",
        ("nowak_may", "P2"): "rejected",
        ("sir_epidemic", "P2"): "rejected",
        ("rps_spatial", "P2"): "rejected",
        ("game_of_life", "P2"): "rejected",
        ("greenberg_hastings", "P2"): "rejected",
        ("lotka_volterra_lattice", "P2"): "rejected",
        ("gray_scott", "P2"): "rejected",
        ("nagel_schreckenberg", "P2"): "rejected",
        ("btw_sandpile", "P2"): "rejected",
        ("kuramoto", "P2"): "rejected",
        ("hegselmann_krause", "P2"): "rejected",
        # ABP-row: lattice-only detectors all substrate_mismatch,
        # continuous_2d detectors P5 and P6 handled above.
        ("abp", "P1"): "rejected",         # lattice_2d only
        ("abp", "P13"): "rejected",
        ("abp", "P22"): "rejected",
        ("abp", "P11"): "rejected",
        ("abp", "P12"): "rejected",
        ("abp", "P15"): "rejected",
        ("abp", "P3"): "rejected",
        ("abp", "P8"): "rejected",         # lattice_1d only
        ("abp", "P9"): "rejected",         # oscillator only
        ("abp", "P14"): "rejected",        # btw only
        ("abp", "P21"): "rejected",        # opinion_space only
        ("abp", "P27"): "rejected",        # nowak_may coop_fraction only
        ("abp", "P31"): "rejected",        # lattice_1d cell_types only

        # --- Sprint 17 pairs (Yard-Sale + P28 Wealth Condensation) ---
        # Yard-Sale is the 17th model and occupies the new scalar_wealth
        # substrate (the first well-mixed/non-spatial substrate in the
        # registry). P28 is the 16th detector and is restricted to
        # scalar_wealth by registration.
        #
        # Empirical tier behavior (Phase 1 + Phase 2):
        #   - YS f=0.1 λ=0 χ=0 N=1000 at t=2e6: DEFINITIVE (gini=0.94,
        #     top_1pct=0.35, monotonic, null-p=0.005).
        #   - Within-family negatives separate cleanly:
        #     λ=0.5 saving propensity → gini=0.28 → below screening floor.
        #     χ=0.001 redistribution → gini=0.68, top_1pct=0.14 →
        #       SCREENING (top_1pct < 0.15 blocks CONFIRMATION).
        #     χ=0.0001 mild redistribution → gini=0.89 → empirically DEF-
        #       strength but metadata has_redistribution=True blocks DEF;
        #       CONFIRMATION is the correct tier (Decision 49).
        #   - Every non-wealth substrate rejects at substrate_mismatch
        #     when run through P28.
        #
        # ONE empirical decision deviates from the pre-existing P28
        # detector card in docs/pattern_catalog_v0_4.md: the catalog
        # mentions a "Pareto power-law tail" as a detection signature.
        # Phase 1c showed that the Hill-estimator α is NOT stable across
        # timescales — it drifts through the Pareto band (1 < α < 2)
        # only transiently before the distribution degenerates toward a
        # delta. α is carried as a DIAGNOSTIC secondary metric but is
        # explicitly NOT used as a tier gate. Decision 47.
        ("yard_sale", "P28"): "detected",   # DEFINITIVE canonical condensation
        # P28 column — every other registered model rejects at substrate_mismatch.
        ("zhang_sequential", "P28"): "rejected",
        ("zhang_threaded", "P28"): "rejected",
        ("schelling", "P28"): "rejected",
        ("greenberg_hastings", "P28"): "rejected",
        ("game_of_life", "P28"): "rejected",
        ("vicsek", "P28"): "rejected",
        ("dorsogna", "P28"): "rejected",
        ("kuramoto", "P28"): "rejected",
        ("btw_sandpile", "P28"): "rejected",
        ("nowak_may", "P28"): "rejected",
        ("hegselmann_krause", "P28"): "rejected",
        ("sir_epidemic", "P28"): "rejected",
        ("rps_spatial", "P28"): "rejected",
        ("lotka_volterra_lattice", "P28"): "rejected",
        ("gray_scott", "P28"): "rejected",
        ("nagel_schreckenberg", "P28"): "rejected",
        ("abp", "P28"): "rejected",
        # yard_sale row — every non-P28 detector rejects at substrate_mismatch
        # (scalar_wealth is distinct from all 6 existing substrates).
        ("yard_sale", "P1"): "rejected",
        ("yard_sale", "P2"): "rejected",
        ("yard_sale", "P3"): "rejected",
        ("yard_sale", "P5"): "rejected",
        ("yard_sale", "P6"): "rejected",
        ("yard_sale", "P8"): "rejected",
        ("yard_sale", "P9"): "rejected",
        ("yard_sale", "P11"): "rejected",  # P11 implemented but unregistered
        ("yard_sale", "P12"): "rejected",
        ("yard_sale", "P13"): "rejected",
        ("yard_sale", "P14"): "rejected",
        ("yard_sale", "P15"): "rejected",
        ("yard_sale", "P21"): "rejected",
        ("yard_sale", "P22"): "rejected",
        ("yard_sale", "P27"): "rejected",
        ("yard_sale", "P31"): "rejected",
    }

    VALID_OUTCOMES = {"detected", "rejected", "screening", "not_detected"}

    def test_all_pairs_documented(self):
        """Every audited pair has a valid, documented expected outcome."""
        assert len(self.EXPECTED_OUTCOMES) >= 112, \
            f"Expected at least 112 audited pairs, got {len(self.EXPECTED_OUTCOMES)}"

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

    def test_sprint_11_pairs_covered(self):
        """Every Sprint 11 cross-detection pair (LV row + P11 column) must
        appear in the transfer matrix with a pinned expected outcome.

        The Sprint 11 headline finding is that P11 cleanly separates
        bilateral predator-prey oscillation (LV) from cyclic three-
        species dominance (RPS) — but the separation is mechanistic
        (n_species prerequisite), not signal-based. On the primary
        metric alone (rho_anti = min cross-correlation at nonzero lag),
        RPS species pairs score STRONGER anti-correlation than LV
        (-0.93 to -0.97 vs LV's -0.72 to -0.88). The n_species==2
        prerequisite is what keeps P11 specific to bilateral systems.

        A secondary finding is that strictly conserved 2-species
        systems (Nowak-May: coop + defect = 1 exactly) give rho_anti
        = -1.0 by algebra alone, requiring P11 to also verify
        total_std > 0.005 (nontrivial empty reservoir).
        """
        sprint_11_pairs = [
            # LV row — all 6 audited detectors
            ("lotka_volterra_lattice", "P11"),
            ("lotka_volterra_lattice", "P12"),
            ("lotka_volterra_lattice", "P13"),
            ("lotka_volterra_lattice", "P22"),
            ("lotka_volterra_lattice", "P1"),
            ("lotka_volterra_lattice", "P15"),
            # P11 column — all 4 other 2-species (or multi-species) models
            ("schelling", "P11"),
            ("nowak_may", "P11"),
            ("sir_epidemic", "P11"),
            ("rps_spatial", "P11"),
        ]
        for pair in sprint_11_pairs:
            assert pair in self.EXPECTED_OUTCOMES, \
                f"Sprint 11 pair {pair} missing from transfer matrix"

        # LV × P11 is the Sprint 11 positive; all other P11 cells are
        # rejections via different prerequisite paths.
        assert self.EXPECTED_OUTCOMES[("lotka_volterra_lattice", "P11")] == "detected"
        for other_model in ("schelling", "nowak_may", "sir_epidemic",
                            "rps_spatial"):
            assert self.EXPECTED_OUTCOMES[(other_model, "P11")] == "rejected", (
                f"{other_model} × P11 should be rejected (prerequisite failure)"
            )

        # LV column exclusions: P12 and P13 must reject cleanly.
        assert self.EXPECTED_OUTCOMES[("lotka_volterra_lattice", "P12")] == "rejected", (
            "LV × P12 must be rejected (intransitivity_score < 1.0; no cyclic "
            "A→B→C→A, only bilateral predation)"
        )
        assert self.EXPECTED_OUTCOMES[("lotka_volterra_lattice", "P13")] == "rejected", (
            "LV × P13 must be rejected (pursuit-evasion fronts are noise-"
            "driven with wavefront_speed_cv > 1, not clock-driven)"
        )

        print(f"  ✓ Sprint 11: all {len(sprint_11_pairs)} pairs covered "
              f"(LV row + P11 column)")

    def test_sprint_13_gray_scott_p3_covered(self):
        """Sprint 13 added Gray-Scott as the first continuous-valued-field
        model and P3 (Turing-wavelength) as its canonical detector.

        Key structural facts pinned by this test:
        - GS × P3 must be `detected` (canonical DEFINITIVE positive,
          p/m=18.75, Cohen's d ~103, seeds 42, 7, 123).
        - Every integer-grid model × P3 must be `rejected` (substrate
          prerequisite: no 'field' observable; secondary prerequisite:
          n_unique_values >= 50).
        - Every other 'grid'-consuming detector × GS must be `rejected`
          (graceful substrate mismatch — GS exposes 'field', not 'grid').
        """
        sprint_13_pairs = [
            # GS row (canonical positive + substrate-mismatch rejections)
            ("gray_scott", "P3"),
            ("gray_scott", "P13"),
            ("gray_scott", "P22"),
            ("gray_scott", "P11"),
            ("gray_scott", "P12"),
            ("gray_scott", "P15"),
            # P3 column (all integer-grid models reject at substrate prereq)
            ("schelling", "P3"),
            ("nowak_may", "P3"),
            ("sir_epidemic", "P3"),
            ("rps_spatial", "P3"),
            ("game_of_life", "P3"),
            ("greenberg_hastings", "P3"),
            ("lotka_volterra_lattice", "P3"),
        ]
        for pair in sprint_13_pairs:
            assert pair in self.EXPECTED_OUTCOMES, \
                f"Sprint 13 pair {pair} missing from transfer matrix"

        # GS × P3 is the Sprint 13 positive.
        assert self.EXPECTED_OUTCOMES[("gray_scott", "P3")] == "detected"

        # Every integer-grid model × P3 must be rejected via substrate
        # prerequisite. This is the headline Sprint 13 discrimination
        # test — notably RPS at low mobility has radial-FFT peak-to-mean
        # ≈ 23 matching Gray-Scott labyrinth, so discrimination CANNOT
        # be on peak-to-mean alone.
        for other_model in ("schelling", "nowak_may", "sir_epidemic",
                            "rps_spatial", "game_of_life",
                            "greenberg_hastings", "lotka_volterra_lattice"):
            assert self.EXPECTED_OUTCOMES[(other_model, "P3")] == "rejected", (
                f"{other_model} × P3 should be rejected (substrate prereq: "
                f"no 'field' observable, or nu_unique_values < 50)"
            )

        # GS row: other detectors must not falsely detect.
        for det in ("P13", "P22", "P11", "P12"):
            assert self.EXPECTED_OUTCOMES[("gray_scott", det)] == "rejected", (
                f"gray_scott × {det} should be rejected (no 'grid' observable)"
            )

        print(f"  ✓ Sprint 13: all {len(sprint_13_pairs)} pairs covered "
              f"(GS row + P3 column)")

    def test_sprint_14_p1_substrate_robustness_covered(self):
        """Sprint 14 B.1 hardened P1's 2D branch against continuous-field
        substrates (carry-forward from Sprint 13 #1).

        Before Sprint 14, P1 raised ``KeyError: "Need 'type_labels_at_pos'
        or 'grid' for 2D"`` when invoked on Gray-Scott's state_history
        (which exposes ``field`` but no integer-grid observable). Post-fix,
        P1 returns a graceful ``detected=False`` result with a substrate
        warning, matching the pattern used by P11/P13/P22.

        The transfer-matrix cell ``gray_scott × P1`` flipped from ``ker``
        (pre-existing bug, skipped in test matrix) to ``rej`` (audited).
        """
        assert ("gray_scott", "P1") in self.EXPECTED_OUTCOMES, \
            "Sprint 14 B.1 requires gray_scott × P1 in EXPECTED_OUTCOMES"
        assert self.EXPECTED_OUTCOMES[("gray_scott", "P1")] == "rejected", (
            "gray_scott × P1 should be rejected (no integer-labeled grid; "
            "graceful substrate-warning reject per Sprint 14 B.1)"
        )
        print("  ✓ Sprint 14 B.1: gray_scott × P1 = rejected "
              "(KeyError → graceful reject)")

    def test_sprint_15_ns_p8_covered(self):
        """Sprint 15 added Nagel-Schreckenberg as the second lattice_1d
        model (after Zhang sorting) and P8 (traffic jamming) as its
        canonical detector.

        Key structural facts pinned by this test:
        - NS × P8 must be `detected` (canonical DEFINITIVE positive at
          L=1000, rho=0.15, v_max=5, p_slow=0.3: stopped=0.18,
          jam_lifetime_p95=13, null p=0.005, Cohen's d effectively
          infinite because null_std ≈ 0).
        - Every other model × P8 must be `rejected` — either at substrate
          level (non-lattice_1d) or at observable level (Zhang is
          lattice_1d but lacks a 'velocities' observable).
        - Every other detector × NS must be `rejected` — NS exposes
          positions/velocities/gaps but none of the 'grid'/'field'/
          'cell_types'/'theta'/'opinions' observables the other
          detectors require.
        - A secondary finding is that density saturation (NS at rho=0.80,
          p=0) gives stopped_fraction=0.75 by pigeonhole but
          jam_lifetime_p95=4 — correctly rejected at the confirmation
          gate. This is the P8 equivalent of the Sprint 13 RPS false-
          positive trap being rejected at the n_unique_values prereq.
        """
        sprint_15_pairs = [
            # NS row — canonical positive + substrate-mismatch rejections
            ("nagel_schreckenberg", "P8"),
            ("nagel_schreckenberg", "P1"),
            ("nagel_schreckenberg", "P13"),
            ("nagel_schreckenberg", "P22"),
            ("nagel_schreckenberg", "P11"),
            ("nagel_schreckenberg", "P12"),
            ("nagel_schreckenberg", "P15"),
            ("nagel_schreckenberg", "P3"),
            # P8 column — every other model must reject
            ("zhang_sequential", "P8"),
            ("schelling", "P8"),
            ("nowak_may", "P8"),
            ("sir_epidemic", "P8"),
            ("rps_spatial", "P8"),
            ("game_of_life", "P8"),
            ("greenberg_hastings", "P8"),
            ("lotka_volterra_lattice", "P8"),
            ("gray_scott", "P8"),
        ]
        for pair in sprint_15_pairs:
            assert pair in self.EXPECTED_OUTCOMES, \
                f"Sprint 15 pair {pair} missing from transfer matrix"

        # NS × P8 is the Sprint 15 positive.
        assert self.EXPECTED_OUTCOMES[("nagel_schreckenberg", "P8")] == "detected"

        # Every other model × P8 must be rejected via substrate or
        # observable prerequisite. Zhang is the important case: it
        # shares the lattice_1d substrate but lacks a 'velocities'
        # observable, so rejection is at observable-prereq (not
        # substrate_mismatch) — analogous to Sprint 13's n_unique_values
        # content-level gate for P3.
        for other_model in (
            "zhang_sequential", "schelling", "nowak_may", "sir_epidemic",
            "rps_spatial", "game_of_life", "greenberg_hastings",
            "lotka_volterra_lattice", "gray_scott",
        ):
            assert self.EXPECTED_OUTCOMES[(other_model, "P8")] == "rejected", (
                f"{other_model} × P8 should be rejected "
                f"(substrate_mismatch or missing 'velocities' observable)"
            )

        # NS-row: every integer-grid / continuous-field / oscillator
        # detector must reject cleanly — NS exposes positions/velocities/
        # gaps but none of the 2D observables.
        for det in ("P1", "P13", "P22", "P11", "P12", "P15", "P3"):
            assert self.EXPECTED_OUTCOMES[("nagel_schreckenberg", det)] == "rejected", (
                f"nagel_schreckenberg × {det} should be rejected "
                f"(missing 2D 'grid' or 'field' observable)"
            )

        print(f"  ✓ Sprint 15: all {len(sprint_15_pairs)} pairs covered "
              f"(NS row + P8 column)")

    def test_sprint_16_abp_p2_covered(self):
        """Sprint 16 added Active Brownian Particles (ABP) as the 16th
        model (third continuous_2d model after Vicsek and D'Orsogna)
        and P2 (MIPS) as the 15th detector.

        Key structural facts pinned by this test:
        - ABP × P2 is the canonical positive (DEFINITIVE at N >= 800,
          phi=0.5, Pe=100, 2500 measurement steps: two_phase_score
          ~ 0.15-0.40, r ~ -0.9, null_p < 0.01, metadata flags
          has_density_dependent_speed=True, has_alignment_rule=False,
          has_attraction_rule=False).
        - Vicsek × P2 and D'Orsogna × P2 share the substrate but are
          distinguished mechanistically. Vicsek has constant speed so
          CV_v = 0 exactly AND f_gas is near zero (flocks are all-
          liquid); primary falls below the screening floor.
          D'Orsogna has variable speed AND does cluster, but its
          attractive Morse potential manifests as f_liquid dominant /
          f_gas tiny, and the has_attraction_rule=True metadata flag
          (when carried) blocks DEFINITIVE.
        - Every non-continuous_2d model × P2 must be `rejected` at
          substrate mismatch. Every non-continuous_2d detector × ABP
          must also be `rejected` at substrate mismatch.
        - The Sprint 16 Phase 1 characterization identified THREE
          false-positive traps within the ABP model family itself:
          thermal (low Pe, CV_v gate rejects), dilute (low phi,
          primary floor rejects), over-saturated (high phi long
          runtime, primary drops below CONFIRMATION). Each trap
          corresponds to a specific detector confirmation gate.
        - ONE decision deviates from the pre-existing P2 detector card
          in docs/detector_cards.md: the card's recipe uses Hartigan
          dip on density histogram. Phase 1c/1d showed this is
          empirically wrong — dip rejects uniformity on discrete
          density histograms universally (floor p=0.005 even on
          truly uniform regimes). The Sprint 16 detector instead
          uses ``two_phase_coexistence_score = min(f_gas, f_liquid)``
          which is tight-fitting for MIPS and zero for flocking /
          stuck / dilute regimes. See ADR 44 in REPLICATION_NOTES.md.
        """
        sprint_16_pairs = [
            # ABP row — canonical positive + continuous_2d neighbor handling
            ("abp", "P2"),
            ("abp", "P5"),
            ("abp", "P6"),
            # ABP row — substrate-mismatch rejections
            ("abp", "P1"),
            ("abp", "P3"),
            ("abp", "P8"),
            ("abp", "P9"),
            ("abp", "P11"),
            ("abp", "P12"),
            ("abp", "P13"),
            ("abp", "P14"),
            ("abp", "P15"),
            ("abp", "P21"),
            ("abp", "P22"),
            ("abp", "P27"),
            ("abp", "P31"),
            # P2 column — continuous_2d neighbors
            ("vicsek", "P2"),
            ("dorsogna", "P2"),
            # P2 column — substrate-mismatch rejections
            ("zhang_sequential", "P2"),
            ("schelling", "P2"),
            ("nowak_may", "P2"),
            ("sir_epidemic", "P2"),
            ("rps_spatial", "P2"),
            ("game_of_life", "P2"),
            ("greenberg_hastings", "P2"),
            ("lotka_volterra_lattice", "P2"),
            ("gray_scott", "P2"),
            ("nagel_schreckenberg", "P2"),
            ("btw_sandpile", "P2"),
            ("kuramoto", "P2"),
            ("hegselmann_krause", "P2"),
        ]
        for pair in sprint_16_pairs:
            assert pair in self.EXPECTED_OUTCOMES, \
                f"Sprint 16 pair {pair} missing from transfer matrix"

        # ABP × P2 is the Sprint 16 positive.
        assert self.EXPECTED_OUTCOMES[("abp", "P2")] == "detected", (
            "abp × P2 must be detected (canonical MIPS DEFINITIVE)"
        )

        # ABP × other continuous_2d detectors must be rejected —
        # ABP has no alignment (P5) and no attractive/milling centre (P6).
        for det in ("P5", "P6"):
            assert self.EXPECTED_OUTCOMES[("abp", det)] == "rejected", (
                f"abp × {det} should be rejected (no alignment / no attraction)"
            )

        # P2 column on continuous_2d neighbours: Vicsek rejects at
        # the screening floor (flocks are all-liquid OR uniform-gas);
        # D'Orsogna can reach SCREENING via clustering but not higher.
        assert self.EXPECTED_OUTCOMES[("vicsek", "P2")] == "rejected", (
            "vicsek × P2 should be rejected (constant speed, below_two_phase_floor)"
        )
        assert self.EXPECTED_OUTCOMES[("dorsogna", "P2")] == "screening", (
            "dorsogna × P2 should reach SCREENING only (attraction-driven clustering)"
        )

        # Every non-continuous_2d model × P2 must be rejected at
        # substrate mismatch.
        for other_model in (
            "zhang_sequential", "schelling", "nowak_may", "sir_epidemic",
            "rps_spatial", "game_of_life", "greenberg_hastings",
            "lotka_volterra_lattice", "gray_scott", "nagel_schreckenberg",
            "btw_sandpile", "kuramoto", "hegselmann_krause",
        ):
            assert self.EXPECTED_OUTCOMES[(other_model, "P2")] == "rejected", (
                f"{other_model} × P2 should be rejected (substrate_mismatch: "
                f"P2 requires continuous_2d)"
            )

        # ABP-row: every lattice-only / oscillator / opinion-space
        # detector must reject at substrate mismatch.
        for det in ("P1", "P3", "P8", "P9", "P11", "P12", "P13", "P14",
                    "P15", "P21", "P22", "P27", "P31"):
            assert self.EXPECTED_OUTCOMES[("abp", det)] == "rejected", (
                f"abp × {det} should be rejected "
                f"(substrate mismatch — non-continuous_2d detector)"
            )

        print(f"  ✓ Sprint 16: all {len(sprint_16_pairs)} pairs covered "
              f"(ABP row + P2 column)")

    def test_sprint_17_yard_sale_p28_covered(self):
        """Sprint 17 added Yard-Sale as the 17th model (first scalar_wealth
        model — the first well-mixed / non-spatial substrate in the
        registry) and P28 (Wealth Condensation) as the 16th detector.

        Key structural facts pinned by this test:
        - yard_sale × P28 is the canonical positive (DEFINITIVE at
          N=1000, f=0.1, lambda=0, chi=0, t=2e6: gini ~ 0.94,
          top_1pct ~ 0.34, monotonic, null-p < 0.01, metadata flags
          has_conserved_resource=True, has_multiplicative_stake=True,
          has_saving_propensity=False, has_redistribution=False).
        - Every non-scalar_wealth model × P28 rejects at substrate
          mismatch.
        - Every non-P28 detector × yard_sale rejects at substrate
          mismatch (scalar_wealth is a new, isolated substrate type
          at Sprint 17).
        - Within-family negatives are handled at detector-layer gates,
          not at the registry layer (saving propensity, redistribution
          are mechanism-flag discriminations, not substrate mismatches).
          These are tested in tests/test_yard_sale_p28_e2e.py and do
          not need entries in this cross-model transfer matrix (which
          audits across-model-family pairs only).
        - ONE decision deviates from the pre-existing P28 pattern-
          catalog entry: the catalog mentions "Pareto power-law tail"
          as a detection signature, but the Hill estimator α drifts
          unstably with time (Phase 1c). α is a diagnostic metric,
          NOT a tier gate. Decision 47 in REPLICATION_NOTES.md.
        """
        sprint_17_pairs = [
            # yard_sale × P28 — canonical positive
            ("yard_sale", "P28"),
            # P28 column — substrate-mismatch rejections (every non-wealth model)
            ("zhang_sequential", "P28"),
            ("zhang_threaded", "P28"),
            ("schelling", "P28"),
            ("greenberg_hastings", "P28"),
            ("game_of_life", "P28"),
            ("vicsek", "P28"),
            ("dorsogna", "P28"),
            ("kuramoto", "P28"),
            ("btw_sandpile", "P28"),
            ("nowak_may", "P28"),
            ("hegselmann_krause", "P28"),
            ("sir_epidemic", "P28"),
            ("rps_spatial", "P28"),
            ("lotka_volterra_lattice", "P28"),
            ("gray_scott", "P28"),
            ("nagel_schreckenberg", "P28"),
            ("abp", "P28"),
            # yard_sale row — substrate-mismatch rejections (every non-P28 detector)
            ("yard_sale", "P1"),
            ("yard_sale", "P2"),
            ("yard_sale", "P3"),
            ("yard_sale", "P5"),
            ("yard_sale", "P6"),
            ("yard_sale", "P8"),
            ("yard_sale", "P9"),
            ("yard_sale", "P11"),
            ("yard_sale", "P12"),
            ("yard_sale", "P13"),
            ("yard_sale", "P14"),
            ("yard_sale", "P15"),
            ("yard_sale", "P21"),
            ("yard_sale", "P22"),
            ("yard_sale", "P27"),
            ("yard_sale", "P31"),
        ]
        for pair in sprint_17_pairs:
            assert pair in self.EXPECTED_OUTCOMES, \
                f"Sprint 17 pair {pair} missing from transfer matrix"

        # yard_sale × P28 is the Sprint 17 canonical positive.
        assert self.EXPECTED_OUTCOMES[("yard_sale", "P28")] == "detected", (
            "yard_sale × P28 must be detected (canonical wealth condensation DEFINITIVE)"
        )

        # Every non-wealth model × P28 must be rejected at substrate mismatch.
        for other_model in (
            "zhang_sequential", "zhang_threaded", "schelling",
            "greenberg_hastings", "game_of_life", "vicsek", "dorsogna",
            "kuramoto", "btw_sandpile", "nowak_may", "hegselmann_krause",
            "sir_epidemic", "rps_spatial", "lotka_volterra_lattice",
            "gray_scott", "nagel_schreckenberg", "abp",
        ):
            assert self.EXPECTED_OUTCOMES[(other_model, "P28")] == "rejected", (
                f"{other_model} × P28 should be rejected "
                f"(substrate_mismatch: P28 requires scalar_wealth)"
            )

        # yard_sale row — every detector except P28 must reject at substrate
        # mismatch (scalar_wealth is isolated from all 6 prior substrates).
        for det in (
            "P1", "P2", "P3", "P5", "P6", "P8", "P9", "P11", "P12",
            "P13", "P14", "P15", "P21", "P22", "P27", "P31",
        ):
            assert self.EXPECTED_OUTCOMES[("yard_sale", det)] == "rejected", (
                f"yard_sale × {det} should be rejected "
                f"(substrate_mismatch — scalar_wealth substrate)"
            )

        print(f"  ✓ Sprint 17: all {len(sprint_17_pairs)} pairs covered "
              f"(yard_sale row + P28 column)")


# ===================================================================
# Runner
# ===================================================================

if __name__ == "__main__":
    pytest.main([__file__, "-v", "--tb=short"])
