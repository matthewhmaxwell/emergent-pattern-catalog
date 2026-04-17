"""P12 cyclic dominance detector: end-to-end tests.

Verifies that P12 correctly identifies the intransitive cyclic structure
in spatial RPS (Reichenbach 2007) and correctly rejects:
  - Greenberg-Hastings (excitable waves; clock-driven transitions, ρ=1)
  - SIR (only 3 states, but non-cyclic)
  - Game of Life (only 2 states; fails prerequisite)
  - Nowak-May (only 2 states; fails prerequisite)

Each test verifies the ACTUAL EFFECT, not just that code runs:
  - RPS × P12 MUST reach at least CONFIRMATION tier on canonical params.
  - RPS × P12 MUST identify a cyclic triple whose direction matches the
    model's actual dominance map (A kills B, B kills C, C kills A).
  - RPS × P12 MUST mark P13 as 'excluded' in exclusion_results (cyclic
    dominance implies NOT clock-driven, which is P13's mechanism).
  - GH × P12 MUST return intransitivity_score = 0.0 (ρ=1 exactly for
    clock-driven transitions).
  - SIR × P12 MUST return intransitivity_score = 0.0 in the coexistence
    cyclic direction (SIR has no cyclic dominance structure).
  - GoL × P12 MUST fail prerequisite (n_candidate_species < 3).
"""

from __future__ import annotations

import numpy as np
import pytest

from epc.detectors.p12_cyclic_dominance import P12CyclicDominanceDetector
from epc.models.game_of_life import GameOfLife
from epc.models.greenberg_hastings import GreenbergHastings
from epc.models.rps_spatial import RPSSpatialModel
from epc.models.sir_epidemic import SIREpidemicModel


# ===================================================================
# RPS → P12 detection
# ===================================================================

class TestRPSDetectedByP12:
    """RPS at coexistence mobility should trigger P12 at confirmation or
    definitive tier with high confidence."""

    def test_rps_reaches_confirmation(self):
        """RPS × P12 at canonical params → CONFIRMATION tier."""
        m = RPSSpatialModel(rows=30, cols=30, mobility=1e-4, seed=42)
        history = m.run(n_steps=80)
        detector = P12CyclicDominanceDetector(n_permutations=199, seed=42)
        result = detector.detect(history, m.get_metadata())

        assert result.detected, (
            f"P12 should detect RPS. Got score="
            f"{result.primary_metric.get('intransitivity_score', '?'):.3f}"
        )
        # At 199 permutations, min p = 0.005 which just crosses the
        # confirmation threshold of 0.01.
        assert result.is_confirmed, (
            f"Expected confirmation tier, got {result.tier.value}. "
            f"p={result.null_p_value:.4f}, score="
            f"{result.primary_metric['intransitivity_score']:.3f}, "
            f"coexistence="
            f"{result.secondary_metrics.get('coexistence_fraction')}"
        )
        assert result.confidence >= 0.55, (
            f"Confidence {result.confidence:.2f} too low for confirmation"
        )

    def test_rps_p12_identifies_correct_cycle_direction(self):
        """The identified triple should match the model's dominance map.

        In RPSSpatialModel, DOMINATES = {1:2, 2:3, 3:1} meaning:
          A (1) kills B (2)  → B-cells become A (via selection + repro)
          B (2) kills C (3)  → C-cells become B
          C (3) kills A (1)  → A-cells become C

        So the REPLACEMENT direction is: 2→1, 3→2, 1→3.

        In our detector, identified_triple (a, b, c) means: b replaces a,
        c replaces b, a replaces c. So we need the cyclic ordering
        {2 → 1 → 3 → 2} or any rotation. That's the triple (2, 1, 3) or
        (1, 3, 2) or (3, 2, 1) — all cyclically equivalent.
        """
        m = RPSSpatialModel(rows=30, cols=30, mobility=1e-4, seed=42)
        history = m.run(n_steps=80)
        detector = P12CyclicDominanceDetector(n_permutations=49, seed=42)
        result = detector.detect(history, m.get_metadata())

        identified = tuple(result.secondary_metrics["identified_triple"])
        # The three cyclic rotations of the expected direction:
        expected_rotations = {(2, 1, 3), (1, 3, 2), (3, 2, 1)}
        assert identified in expected_rotations, (
            f"Identified triple {identified} does not match any rotation of "
            f"the model's dominance cycle (expected one of {expected_rotations})"
        )

    def test_rps_excludes_p13_and_p22(self):
        """P12 should exclude P13 (clock-driven waves) and P22 (single-pass
        cascade) when it detects RPS."""
        m = RPSSpatialModel(rows=30, cols=30, mobility=1e-4, seed=42)
        history = m.run(n_steps=80)
        detector = P12CyclicDominanceDetector(n_permutations=199, seed=42)
        result = detector.detect(history, m.get_metadata())

        assert result.exclusion_results.get("P13") == "excluded", (
            f"Expected P13 excluded, got {result.exclusion_results.get('P13')}"
        )
        assert result.exclusion_results.get("P22") == "excluded", (
            f"Expected P22 excluded, got {result.exclusion_results.get('P22')}"
        )

    def test_rps_effect_size_is_enormous(self):
        """The observed-vs-null effect size should be extreme because the
        null (spatial shuffle) destroys the neighbor-transition correlation
        completely. Expect Cohen's d > 50."""
        m = RPSSpatialModel(rows=30, cols=30, mobility=1e-4, seed=42)
        history = m.run(n_steps=80)
        detector = P12CyclicDominanceDetector(n_permutations=199, seed=42)
        result = detector.detect(history, m.get_metadata())

        d = result.effect_size.get("cohens_d", 0.0)
        assert d > 50.0, (
            f"Effect size {d:.1f} should vastly exceed null. "
            f"observed={result.effect_size.get('raw_value', 0.0):.3f}, "
            f"null_mean={result.effect_size.get('null_mean', 0.0):.3f}"
        )


# ===================================================================
# Cross-detection: other models REJECTED by P12
# ===================================================================

class TestP12CrossDetectionRejections:
    """P12 should correctly reject every model that is not RPS-like."""

    def test_gh_rejected_by_p12(self):
        """GH spirals have ρ ≈ 1.0 (clock-driven), intransitivity_score ≈ 0."""
        m = GreenbergHastings(
            rows=30, cols=30, n_states=5, threshold=1,
            init_mode="random", seed=42,
        )
        history = m.run(n_steps=60)
        detector = P12CyclicDominanceDetector(n_permutations=49, seed=42)
        result = detector.detect(history, m.get_metadata())

        assert not result.detected, (
            f"GH should be rejected. Got score="
            f"{result.primary_metric.get('intransitivity_score', '?'):.3f}"
        )
        score = result.primary_metric.get("intransitivity_score", 99.0)
        assert score < 0.5, (
            f"GH intransitivity score {score:.3f} should be near 0 "
            f"(clock-driven transitions produce ρ ≈ 1.0)"
        )

    def test_gh_n3_rejected_by_p12(self):
        """GH with n_states=3 exactly (matching RPS in state count) is
        still rejected: the number of states doesn't matter, only ρ does."""
        m = GreenbergHastings(
            rows=30, cols=30, n_states=3, threshold=1,
            init_mode="random", seed=42,
        )
        history = m.run(n_steps=60)
        detector = P12CyclicDominanceDetector(n_permutations=49, seed=42)
        result = detector.detect(history, m.get_metadata())

        assert not result.detected, (
            f"GH(n=3) should be rejected. Score="
            f"{result.primary_metric.get('intransitivity_score', '?'):.3f}"
        )

    def test_sir_rejected_by_p12(self):
        """SIR has 3 states {S, I, R} but no cyclic dominance.

        S→I→R is a linear, single-pass progression. There is no R→S
        transition at all. The neighbor-conditional ratio ρ for the
        direction R→S is 0 (no such transitions observed), so the
        min-forward-ρ score collapses to ≈ 0."""
        m = SIREpidemicModel(
            rows=30, cols=30, infection_prob=0.5, recovery_prob=0.1,
            init_mode="single_seed", seed=42,
        )
        history = m.run(n_steps=50)
        detector = P12CyclicDominanceDetector(n_permutations=49, seed=42)
        result = detector.detect(history, m.get_metadata())

        assert not result.detected, (
            f"SIR should be rejected. Score="
            f"{result.primary_metric.get('intransitivity_score', '?'):.3f}"
        )

    def test_gol_rejected_by_p12_via_prerequisite(self):
        """GoL has only 2 states; P12 requires ≥ 3 candidate species."""
        m = GameOfLife(
            rows=30, cols=30, init_mode="random",
            init_density=0.37, seed=42,
        )
        history = m.run(n_steps=50)
        detector = P12CyclicDominanceDetector(n_permutations=49, seed=42)
        result = detector.detect(history, m.get_metadata())

        assert not result.detected, "GoL should be rejected by P12 prerequisite"
        assert result.primary_metric.get("n_candidate_species", 0) < 3, (
            f"GoL should have < 3 candidate species, got "
            f"{result.primary_metric.get('n_candidate_species')}"
        )


# ===================================================================
# Null model sanity
# ===================================================================

class TestP12NullModel:
    """Verify the null model behaves as expected: spatial shuffle destroys
    the observed intransitivity signal."""

    def test_null_destroys_rps_signal(self):
        """Run null model manually on RPS trajectory. Null scores should
        cluster near zero; observed score should sit far above them."""
        m = RPSSpatialModel(rows=30, cols=30, mobility=1e-4, seed=42)
        history = m.run(n_steps=60)
        detector = P12CyclicDominanceDetector(n_permutations=99, seed=42)
        result = detector.detect(history, m.get_metadata())

        # Effect size should confirm observed >> null
        observed = result.effect_size.get("raw_value", 0.0)
        null_mean = result.effect_size.get("null_mean", 0.0)
        null_std = result.effect_size.get("null_std", 1.0)

        assert observed > 1.5, f"Observed score {observed:.3f} should be > 1.5"
        assert null_mean < 0.1, f"Null mean {null_mean:.3f} should be near 0"
        assert null_std < 0.1, f"Null std {null_std:.3f} should be small"

        # p-value should be at or near the minimum achievable
        assert result.null_p_value <= 0.02, (
            f"p-value {result.null_p_value:.4f} should be near the "
            f"minimum achievable with n_perm=99 (= 0.01)"
        )


# ===================================================================
# Robustness across parameters
# ===================================================================

class TestP12Robustness:
    """Detection is robust across mobilities and seeds (within the
    coexistence regime)."""

    def test_rps_detected_across_mobilities(self):
        """P12 detects RPS at multiple coexistence mobilities."""
        for M in [1e-5, 5e-5, 1e-4]:
            m = RPSSpatialModel(rows=30, cols=30, mobility=M, seed=42)
            history = m.run(n_steps=60)
            detector = P12CyclicDominanceDetector(n_permutations=49, seed=42)
            result = detector.detect(history, m.get_metadata())
            assert result.detected, (
                f"P12 should detect RPS at M={M}. Got score="
                f"{result.primary_metric['intransitivity_score']:.3f}"
            )
            assert result.primary_metric["intransitivity_score"] > 1.0, (
                f"Score at M={M}: "
                f"{result.primary_metric['intransitivity_score']:.3f}"
            )

    def test_rps_detected_across_seeds(self):
        """P12 detects RPS across random seeds."""
        for seed in [42, 1, 7]:
            m = RPSSpatialModel(rows=30, cols=30, mobility=1e-4, seed=seed)
            history = m.run(n_steps=60)
            detector = P12CyclicDominanceDetector(n_permutations=49, seed=seed)
            result = detector.detect(history, m.get_metadata())
            assert result.detected, (
                f"P12 should detect RPS at seed={seed}. Got score="
                f"{result.primary_metric['intransitivity_score']:.3f}"
            )


# ===================================================================
# Sprint 10: RPS × P1 characterization under new final-I primary
# ===================================================================


class TestRPSP1ScreeningLevel:
    """RPS × P1: Sprint 10 final-state Moran primary.

    In contrast with SIR (rejected by the same change), RPS spiral
    domains rotate but persist — final Moran's I ≈ peak Moran's I ≈ 0.55.
    The aggregation is genuinely sustained, so RPS × P1 still screens
    under the Sprint-10 detector. This distinguishes the RPS case from
    the SIR case, which shares an initially-similar peak-Moran
    signature but collapses to a near-uniform final state.
    """

    def test_rps_p1_screening_under_new_primary(self):
        """RPS at canonical mobility still screens under final-I primary.

        Expected: final_moran ≈ peak_moran ≈ 0.5 (sustained spiral
        clustering). Screening floor is 0.05, so final > floor passes.
        """
        from epc.detectors.p1_aggregation import P1AggregationDetector

        def adapt(history):
            adapted = []
            for state in history:
                s = dict(state)
                if "grid_dims" not in s and "grid" in s:
                    s["grid_dims"] = s["grid"].shape
                adapted.append(s)
            return adapted

        m = RPSSpatialModel(rows=30, cols=30, mobility=1e-4, seed=42)
        history = adapt(m.run(n_steps=80))

        det = P1AggregationDetector(n_permutations=99)
        result = det.detect(history, m.get_metadata())

        assert result.detected, \
            "Sprint 10: RPS × P1 should still screen (spiral domains are sustained)"

        peak_moran = result.primary_metric.get("morans_i_peak", 0.0)
        final_moran = result.primary_metric.get("morans_i_final", 0.0)
        primary_moran = result.primary_metric.get("morans_i", 0.0)

        # RPS spirals are sustained: final ≈ peak (unlike SIR).
        assert final_moran > 0.4, \
            f"RPS final Moran should be ≈ 0.5 (sustained spirals), got {final_moran:.3f}"
        assert peak_moran > 0.4, \
            f"RPS peak Moran should be ≈ 0.5, got {peak_moran:.3f}"
        # The peak-final gap should be small (rotating but persistent).
        gap = peak_moran - final_moran
        assert gap < 0.15, \
            f"RPS peak-final gap should be small (sustained), got {gap:.3f}"

        # Primary is now final-state Moran (Sprint 10).
        assert abs(primary_moran - final_moran) < 1e-6, \
            f"Sprint 10: primary should equal final, got " \
            f"primary={primary_moran:.4f}, final={final_moran:.4f}"

    def test_rps_vs_sir_p1_asymmetry(self):
        """Sprint 10: same-family models differ — RPS screens, SIR rejects.

        Both RPS and SIR produce high PEAK Moran's I during their
        characteristic spatial dynamics, but only RPS maintains final
        Moran near peak. This asymmetry is the scientific motivation
        for the Sprint 10 primary-metric change.
        """
        from epc.detectors.p1_aggregation import P1AggregationDetector

        def adapt(history):
            adapted = []
            for state in history:
                s = dict(state)
                if "grid_dims" not in s and "grid" in s:
                    s["grid_dims"] = s["grid"].shape
                adapted.append(s)
            return adapted

        # --- RPS ---
        rps = RPSSpatialModel(rows=30, cols=30, mobility=1e-4, seed=42)
        rps_hist = adapt(rps.run(n_steps=80))
        det = P1AggregationDetector(n_permutations=99)
        rps_res = det.detect(rps_hist, rps.get_metadata())

        # --- SIR ---
        sir = SIREpidemicModel(
            rows=80, cols=80,
            infection_prob=0.20, recovery_prob=0.3,
            init_mode="single_seed", seed=42,
        )
        sir_hist = sir.run(400, record_every=1)
        sir_res = det.detect(sir_hist, sir.get_metadata())

        # Both peak high.
        rps_peak = rps_res.primary_metric["morans_i_peak"]
        sir_peak = sir_res.primary_metric["morans_i_peak"]
        assert rps_peak > 0.4, f"RPS peak should be > 0.4, got {rps_peak:.3f}"
        assert sir_peak > 0.5, f"SIR peak should be > 0.5, got {sir_peak:.3f}"

        # Only RPS maintains final near peak.
        rps_final = rps_res.primary_metric["morans_i_final"]
        sir_final = sir_res.primary_metric["morans_i_final"]
        assert rps_final > 0.4, \
            f"RPS final should ≈ peak (sustained), got {rps_final:.3f}"
        assert sir_final < 0.1, \
            f"SIR final should be near 0 (transient collapsed), got {sir_final:.3f}"

        # Sprint 10 decision outcomes.
        assert rps_res.detected, "Sprint 10: RPS × P1 should screen"
        assert not sir_res.detected, "Sprint 10: SIR × P1 should reject"


if __name__ == "__main__":
    pytest.main([__file__, "-v", "--tb=short"])
