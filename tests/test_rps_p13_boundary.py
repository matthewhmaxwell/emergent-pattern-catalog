"""RPS × P13 boundary test — pins the Sprint 9 scientific finding.

The premise of this test file is subtle. The Sprint 9 prompt predicted
that P13 might give a FALSE POSITIVE on RPS: RPS has n_states ≥ 3
(passing P13's hard guard), generates persistent wavefronts at nonzero
speed, and re-excites cells (unlike SIR which dies out). On paper, RPS
should look like a P13 excitable medium.

The ACTUAL OBSERVED BEHAVIOR is different. P13 rejects RPS cleanly at
screening because the wavefront speed CV is far too high (0.59–0.68
across seeds and mobilities tested), comfortably above the 0.2
screening threshold.

Why? Excitable media (GH) have CLOCK-DRIVEN transitions: once excited,
a cell ticks through refractory → rest on a deterministic schedule. The
wavefront speed is set by this clock and is highly uniform (CV ≈ 0.05
for GH spirals). RPS transitions are entirely NEIGHBOR-DRIVEN and
STOCHASTIC: a cell changes state only when a specific-species neighbor
is selected by the Gillespie scheduler. This produces wavefronts that
LOOK spiral-like on visual inspection (same wavelength order) but have
much greater speed variability at the cell-level.

This finding matters because it validates that our P13 detector already
discriminates between excitable and cyclic-dominance spirals via a
mechanism (speed uniformity) that was not specifically engineered for
that discrimination. The test locks in this behavior as a regression
barrier: if a future change to the wavefront-speed estimator loosens
its selectivity, this test will catch the resulting RPS false-positive.

Tests below:
  1. RPS × P13 rejects at screening at coexistence mobility (M=1e-4)
  2. Rejection reason is speed CV > 0.2 (not died_out, not short run)
  3. Rejection is robust across mobilities (1e-5, 1e-4, 3e-4)
  4. Rejection is robust across seeds (42, 1, 7)
  5. Basic comparison vs GH: GH has CV that is dramatically lower
     (characterization test — shows the magnitude of the gap)
"""

from __future__ import annotations

import numpy as np
import pytest

from epc.models.rps_spatial import RPSSpatialModel
from epc.detectors.p13_excitable_waves import P13ExcitableWaveDetector


# ===================================================================
# Core rejection test
# ===================================================================

class TestRPSP13Rejection:
    """The central Sprint 9 finding: P13 rejects RPS cleanly at screening."""

    def test_rps_rejected_by_p13_at_coexistence(self):
        """P13 rejects RPS (M=1e-4 coexistence regime) at screening.

        Expected outcome: detected=False, tier=SCREENING, cv > 0.2.
        """
        m = RPSSpatialModel(rows=60, cols=60, mobility=1e-4, seed=42)
        history = m.run(n_steps=250)
        metadata = m.get_metadata()

        detector = P13ExcitableWaveDetector(n_null_runs=49)
        result = detector.detect(history, metadata)

        assert not result.detected, (
            f"P13 should REJECT RPS (not detect). "
            f"cv={result.primary_metric.get('wavefront_speed_cv', '?'):.3f}"
        )

        cv = result.primary_metric.get("wavefront_speed_cv", float("inf"))
        assert cv > 0.2, (
            f"Rejection must be on speed-CV basis. "
            f"Got cv={cv:.3f}, needed > 0.2 (screening threshold)"
        )

        # Make sure we're not rejecting because of a trivial reason
        # (short run, died out, or n_states<3).
        assert result.primary_metric.get("died_out", 1.0) < 0.5, (
            "Rejection should NOT be due to epidemic dying out; RPS is "
            "re-entrant and activity should persist"
        )
        assert result.primary_metric.get("longest_active_streak", 0) > 100, (
            "Rejection should NOT be due to short active streak"
        )
        # RPS has n_states=4 ≥ 3, so it passes the hard guard. Ensure no
        # n_states warning is present.
        n_states_warnings = [
            w for w in result.warnings if "n_states" in w
        ]
        assert not n_states_warnings, (
            f"RPS has n_states=4 ≥ 3, should NOT trigger n_states warning. "
            f"Got: {n_states_warnings}"
        )

    def test_rps_p13_rejection_robust_across_mobilities(self):
        """Rejection holds across the coexistence mobility range.

        Tests M ∈ {1e-5, 1e-4, 3e-4}. All should reject with cv > 0.2.
        """
        mobilities = [1e-5, 1e-4, 3e-4]
        for M in mobilities:
            m = RPSSpatialModel(rows=50, cols=50, mobility=M, seed=42)
            history = m.run(n_steps=200)
            detector = P13ExcitableWaveDetector(n_null_runs=19)
            result = detector.detect(history, m.get_metadata())

            cv = result.primary_metric.get("wavefront_speed_cv", float("inf"))
            assert not result.detected, (
                f"RPS should be rejected by P13 at M={M}. "
                f"Got detected={result.detected}, cv={cv:.3f}"
            )
            assert cv > 0.2, (
                f"At M={M}, expected cv > 0.2 (screening threshold), got {cv:.3f}"
            )

    def test_rps_p13_rejection_robust_across_seeds(self):
        """Rejection holds across random seeds — not a single-seed artifact."""
        seeds = [42, 1, 7]
        for seed in seeds:
            m = RPSSpatialModel(rows=50, cols=50, mobility=1e-4, seed=seed)
            history = m.run(n_steps=200)
            detector = P13ExcitableWaveDetector(n_null_runs=19)
            result = detector.detect(history, m.get_metadata())

            cv = result.primary_metric.get("wavefront_speed_cv", float("inf"))
            assert not result.detected, (
                f"RPS should be rejected at seed={seed}, got cv={cv:.3f}"
            )
            assert cv > 0.2, f"At seed={seed}, cv={cv:.3f} should be > 0.2"


# ===================================================================
# Mechanism characterization: RPS CV vs GH CV
# ===================================================================

class TestRPSvsGHCVDiscrimination:
    """Characterize the MAGNITUDE of the CV gap between RPS and GH.

    This is the test that documents WHY P13 rejects RPS: RPS wavefronts
    have far higher speed variability than GH spirals. GH is clock-driven
    (highly uniform); RPS is stochastic-neighbor-driven.
    """

    def test_rps_cv_substantially_higher_than_gh(self):
        """GH CV should be well below 0.2 screening threshold; RPS above.

        The literal P13 screening threshold is 0.2. GH typically has
        cv ≈ 0.03–0.1 for well-formed spirals. RPS has cv ≈ 0.5–0.7.
        We check: RPS cv > 3 × GH cv (substantial qualitative gap).
        """
        from epc.models.greenberg_hastings import GreenbergHastings

        # GH spiral-forming parameters (from existing GH literature)
        gh = GreenbergHastings(
            rows=50, cols=50, n_states=8, threshold=1,
            neighborhood="moore", init_mode="random", seed=42,
        )
        gh_history = gh.run(n_steps=200)
        detector = P13ExcitableWaveDetector(n_null_runs=19)
        gh_result = detector.detect(gh_history, gh.get_metadata())
        gh_cv = gh_result.primary_metric.get("wavefront_speed_cv", float("inf"))

        # RPS at coexistence mobility
        rps = RPSSpatialModel(rows=50, cols=50, mobility=1e-4, seed=42)
        rps_history = rps.run(n_steps=200)
        rps_result = detector.detect(rps_history, rps.get_metadata())
        rps_cv = rps_result.primary_metric.get("wavefront_speed_cv", float("inf"))

        # The qualitative finding: RPS cv is substantially larger than GH cv.
        # Empirically RPS cv ≈ 0.65, GH cv ≈ 0.05–0.15, so 3x is a generous
        # margin.
        assert rps_cv > 3.0 * gh_cv, (
            f"Expected RPS cv >> GH cv (at least 3x). "
            f"Got RPS cv={rps_cv:.3f}, GH cv={gh_cv:.3f}, ratio={rps_cv/max(gh_cv, 1e-9):.2f}"
        )

        # GH should be below screening threshold; RPS above
        assert rps_cv > 0.2, f"RPS cv={rps_cv:.3f} should exceed 0.2 threshold"
        # (We don't assert GH < 0.2 here — GH's detection is covered in
        # its own e2e tests. We only care about the RPS/GH gap.)


# ===================================================================
# n_states=4 guard passes for RPS
# ===================================================================

class TestRPSPassesNStatesGuard:
    """Regression test: RPS has n_states=4 (empty + 3 species), so it
    passes P13's hard n_states ≥ 3 guard. The rejection comes from
    speed CV, not from the guard."""

    def test_rps_metadata_has_n_states_ge_3(self):
        """Metadata and state snapshots both report n_states = 4."""
        m = RPSSpatialModel(rows=20, cols=20, mobility=1e-4, seed=42)
        history = m.run(n_steps=5)

        meta = m.get_metadata()
        assert meta["n_states"] == 4, f"Expected n_states=4, got {meta['n_states']}"

        for i, state in enumerate(history):
            assert state["n_states"] == 4, (
                f"State {i} has n_states={state['n_states']} != 4"
            )

    def test_rps_model_class_avoids_ca_substring(self):
        """P13's placeholder exclusion logic checks for 'ca' in model_class.

        If RPS reported model_class containing 'ca' or 'excitable', P13
        would mark P12 as 'excluded' — which would be scientifically wrong
        because RPS IS the canonical P12 positive. Guard that.
        """
        m = RPSSpatialModel(rows=20, cols=20, mobility=1e-4, seed=42)
        meta = m.get_metadata()
        mc = meta["model_class"]
        assert "ca" not in mc.lower(), (
            f"model_class '{mc}' contains 'ca' — will break P13 exclusion logic"
        )
        assert "excitable" not in mc.lower(), (
            f"model_class '{mc}' contains 'excitable' — will break P13 exclusion logic"
        )


if __name__ == "__main__":
    pytest.main([__file__, "-v", "--tb=short"])
