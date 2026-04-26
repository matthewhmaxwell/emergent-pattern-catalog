"""End-to-end tests for the voter model and P18 consensus detector (Sprint 20).

Test classes:
- TestVoterModel: model dynamics, init modes, neighborhoods.
- TestVoterMetrics: spatial metrics (Moran's I, wall density) on raw grids.
- TestP18Primary: detector primary metric extraction.
- TestP18Discrimination: detector tiers vs voter (target) and discriminators
    (GH random, GH broken_wave, GoL random, GoL r_pentomino).
- TestSchellingP18ContentLevel: Sprint 21 carry-forward #20.
    Schelling × P18 content-level negative test confirming pure-metric
    rejection without reliance on the metadata flag (Class 4 of §6.10).
- TestP18MultiSeed: detector behavior across seeds (run-to-run robustness).
- TestSprint20SlowReplication: finite-size slow tests at L ∈ {64, 128, 256}.
"""

from __future__ import annotations

import numpy as np
import pytest

from epc.detector_result import DetectionTier
from epc.detectors.p18_consensus import (
    P18ConsensusDetector,
    _moran_i_moore,
    _wall_density_moore,
)
from epc.models.greenberg_hastings import GreenbergHastings
from epc.models.game_of_life import GameOfLife
from epc.models.schelling import run_schelling
from epc.models.voter import VoterModel


# ============================================================================
# Voter model
# ============================================================================


class TestVoterModel:
    """Voter model basics: init, step, dynamics."""

    def test_setup_random(self):
        m = VoterModel(rows=16, cols=16, init_mode="random", seed=0)
        s = m.setup()
        assert s["grid"].shape == (16, 16)
        assert set(np.unique(s["grid"])) <= {0, 1}
        # Random init: roughly balanced
        frac1 = (s["grid"] == 1).mean()
        assert 0.30 < frac1 < 0.70

    def test_setup_half_and_half(self):
        m = VoterModel(rows=8, cols=8, init_mode="half_and_half", seed=0)
        s = m.setup()
        assert (s["grid"][:, :4] == 0).all()
        assert (s["grid"][:, 4:] == 1).all()

    def test_setup_biased(self):
        m = VoterModel(rows=32, cols=32, init_mode="biased",
                       init_fraction=0.8, seed=0)
        s = m.setup()
        frac1 = (s["grid"] == 1).mean()
        assert 0.70 < frac1 < 0.90

    def test_step_changes_grid(self):
        m = VoterModel(rows=16, cols=16, seed=0)
        m.setup()
        g0 = m.grid.copy()
        m.step()
        # At least some sites should have flipped (almost certainly)
        assert not np.array_equal(m.grid, g0)

    def test_state_dict_keys(self):
        m = VoterModel(rows=16, cols=16, seed=0)
        s = m.setup()
        for key in ("grid", "grid_dims", "step", "magnetization",
                    "abs_magnetization", "wall_density", "moran_i",
                    "consensus_reached", "n_states"):
            assert key in s, f"state dict missing key '{key}'"

    def test_consensus_terminates_run_until_consensus(self):
        # Small lattice should reach consensus quickly
        m = VoterModel(rows=8, cols=8, seed=0)
        h = m.run_until_consensus(max_steps=10000)
        # The final state must be consensus
        assert h[-1]["consensus_reached"]

    def test_metadata(self):
        m = VoterModel(rows=32, cols=32, neighborhood="moore", seed=0)
        md = m.get_metadata()
        assert md["model"] == "voter"
        assert md["substrate"] == "lattice_2d"
        assert md["update"] == "asynchronous_copy_neighbor"
        assert md["has_movement"] is False

    def test_von_neumann_neighborhood(self):
        m = VoterModel(rows=16, cols=16, neighborhood="von_neumann", seed=0)
        m.setup()
        m.step()
        assert m._n_neighbors == 4

    def test_moore_neighborhood(self):
        m = VoterModel(rows=16, cols=16, neighborhood="moore", seed=0)
        m.setup()
        m.step()
        assert m._n_neighbors == 8

    def test_invalid_neighborhood_rejected(self):
        with pytest.raises(ValueError):
            VoterModel(rows=8, cols=8, neighborhood="hexagonal")

    def test_invalid_init_mode_rejected(self):
        m = VoterModel(rows=8, cols=8, init_mode="lemur")
        with pytest.raises(ValueError):
            m.setup()

    def test_dynamics_reduces_wall_density(self):
        """Coarsening: wall density should decrease from random init."""
        m = VoterModel(rows=32, cols=32, seed=0)
        s0 = m.setup()
        for _ in range(50):
            s = m.step()
        assert s["wall_density"] < s0["wall_density"], \
            "voter dynamics should reduce wall density (coarsening)"


# ============================================================================
# Spatial metric helpers
# ============================================================================


class TestVoterMetrics:
    """Sanity-check the Moran's I and wall-density helpers."""

    def test_moran_constant_grid_returns_zero(self):
        g = np.ones((16, 16), dtype=np.int8)
        assert _moran_i_moore(g) == 0.0

    def test_moran_random_grid_near_zero(self):
        rng = np.random.default_rng(0)
        g = (rng.random((32, 32)) < 0.5).astype(np.int8)
        # Random binary grid: Moran's I should be small (< 0.1 typically)
        assert abs(_moran_i_moore(g)) < 0.1

    def test_moran_clustered_grid_positive(self):
        # Half-and-half spatial split: Moran's I should be very positive
        g = np.zeros((32, 32), dtype=np.int8)
        g[:, 16:] = 1
        assert _moran_i_moore(g) > 0.5

    def test_wall_density_constant_grid_zero(self):
        g = np.ones((16, 16), dtype=np.int8)
        assert _wall_density_moore(g) == 0.0

    def test_wall_density_random_near_half(self):
        rng = np.random.default_rng(0)
        g = (rng.random((64, 64)) < 0.5).astype(np.int8)
        # Random grid: half of neighbor pairs differ
        assert 0.45 < _wall_density_moore(g) < 0.55

    def test_wall_density_half_split_low(self):
        # Half-and-half spatial split has only the boundary as wall
        g = np.zeros((32, 32), dtype=np.int8)
        g[:, 16:] = 1
        # Only ~ 2 columns of wall (Moore = 8 neighbors ⇒ small fraction)
        assert _wall_density_moore(g) < 0.10


# ============================================================================
# P18 detector — primary metric extraction
# ============================================================================


class TestP18Primary:
    """Probe P18's primary metric extraction in isolation."""

    def test_extract_grids_returns_int8(self):
        det = P18ConsensusDetector()
        m = VoterModel(rows=16, cols=16, seed=0)
        h = m.run(n_steps=20)
        grids = det._extract_grids(h)
        assert len(grids) == len(h)
        for g in grids:
            assert g.dtype == np.int8

    def test_trajectory_metrics_keys(self):
        det = P18ConsensusDetector()
        m = VoterModel(rows=32, cols=32, seed=0)
        h = m.run(n_steps=100)
        metrics = det._trajectory_metrics(h, timescale=16.0)
        for k in ("moran_spearman_early", "wall_spearman_early",
                  "moran_final_qtr_mean", "wall_final_qtr_mean",
                  "minority_fraction_final"):
            assert k in metrics

    def test_voter_moran_grows_in_early_window(self):
        det = P18ConsensusDetector()
        m = VoterModel(rows=64, cols=64, seed=0)
        h = m.run(n_steps=200)
        metrics = det._trajectory_metrics(h, timescale=32.0)
        # Voter's early-window Moran Spearman must be strongly positive
        assert metrics["moran_spearman_early"] > 0.7

    def test_voter_wall_decays_in_early_window(self):
        det = P18ConsensusDetector()
        m = VoterModel(rows=64, cols=64, seed=0)
        h = m.run(n_steps=200)
        metrics = det._trajectory_metrics(h, timescale=32.0)
        # Voter's early-window wall Spearman must be strongly negative
        assert metrics["wall_spearman_early"] < -0.5

    def test_estimate_timescale_uses_grid_dims(self):
        det = P18ConsensusDetector()
        m = VoterModel(rows=64, cols=64, seed=0)
        h = m.run(n_steps=20)
        tau = det._estimate_timescale(h, None)
        # τ = max(rows, cols) / 2
        assert tau == 32.0


# ============================================================================
# P18 discrimination: voter target vs the four discriminators
# ============================================================================


class TestP18Discrimination:
    """The Sprint 20 §4.20 discriminator table, encoded as tests."""

    def test_voter_reaches_definitive(self):
        """Voter at L=64, 400 sweeps, seed=0: must reach DEFINITIVE."""
        det = P18ConsensusDetector(n_permutations=199, seed=0)
        m = VoterModel(rows=64, cols=64, seed=0)
        h = m.run(n_steps=400)
        result = det.detect(h, model_metadata=m.get_metadata())
        assert result.tier == DetectionTier.DEFINITIVE, \
            f"voter seed=0 should reach DEFINITIVE, got {result.tier.name}"
        assert result.confidence >= 0.85

    def test_voter_p_value_below_strict_gate(self):
        """Voter null p-value must be < 0.01 (the confirmation gate)."""
        det = P18ConsensusDetector(n_permutations=199, seed=0)
        m = VoterModel(rows=64, cols=64, seed=0)
        h = m.run(n_steps=400)
        r = det.detect(h, model_metadata=m.get_metadata())
        assert r.null_p_value < 0.01

    def test_gh_broken_wave_rejected_at_screening(self):
        """GH broken_wave init (P13 spiral): Moran is stationary, rejected."""
        det = P18ConsensusDetector(n_permutations=199, seed=0)
        gh = GreenbergHastings(rows=64, cols=64, n_states=8, threshold=2,
                               init_mode="broken_wave", seed=0)
        h = [gh.setup()] + [gh.step() for _ in range(400)]
        r = det.detect(h, model_metadata={"update": "excitable"})
        assert r.tier == DetectionTier.SCREENING
        assert not r.detected

    def test_gh_random_excluded_from_definitive(self):
        """GH random init: passes screening (Moran grows from rare excited
        cells clustering) but excluded from DEFINITIVE by wall_final < 0.05.
        Reaching CONFIRMATION is acceptable — the three-tier framework is
        designed for this kind of ambiguity at lower tiers."""
        det = P18ConsensusDetector(n_permutations=199, seed=0)
        gh = GreenbergHastings(rows=64, cols=64, n_states=8, threshold=2,
                               init_mode="random", seed=0)
        h = [gh.setup()] + [gh.step() for _ in range(400)]
        r = det.detect(h, model_metadata={"update": "excitable"})
        assert r.tier != DetectionTier.DEFINITIVE, \
            f"GH random should NOT reach DEFINITIVE, got {r.tier.name}"

    def test_gol_random_rejected_at_screening(self):
        """GoL random: Moran plateau ~0.27 < 0.30 screening gate."""
        det = P18ConsensusDetector(n_permutations=199, seed=0)
        gol = GameOfLife(rows=64, cols=64, init_mode="random", seed=0)
        h = [gol.setup()] + [gol.step() for _ in range(400)]
        r = det.detect(h, model_metadata={"update": "gol_b3s23"})
        assert r.tier == DetectionTier.SCREENING
        assert not r.detected

    def test_gol_rpent_rejected_at_screening(self):
        """GoL r_pentomino: high initial Moran, no early-time growth."""
        det = P18ConsensusDetector(n_permutations=199, seed=0)
        gol = GameOfLife(rows=64, cols=64, init_mode="r_pentomino", seed=0)
        h = [gol.setup()] + [gol.step() for _ in range(400)]
        r = det.detect(h, model_metadata={"update": "gol_b3s23"})
        assert r.tier == DetectionTier.SCREENING
        assert not r.detected

    def test_voter_exclusions_cleared(self):
        """At DEFINITIVE, voter must clear all three exclusions: P13, P15, P1."""
        det = P18ConsensusDetector(n_permutations=199, seed=0)
        m = VoterModel(rows=64, cols=64, seed=0)
        h = m.run(n_steps=400)
        r = det.detect(h, model_metadata=m.get_metadata())
        for excl in ("P13", "P15", "P1"):
            assert r.exclusion_results.get(excl) == "excluded", \
                f"voter should exclude {excl}, got {r.exclusion_results}"


# ============================================================================
# Sprint 21 #20: Schelling × P18 content-level negative test
# ============================================================================


class TestSchellingP18ContentLevel:
    """Sprint 21 carry-forward #20.

    Schelling × P18 content-level negative test. Confirms that the
    metric-level signature of Schelling segregation does not pass P18
    confirmation gates without relying on the metadata flag
    ``update = 'asynchronous_copy_neighbor'``.

    Sprint 20 §4.20 stated "Schelling P1: saturates to wall ~0.02,
    rejected here" but lacked a 5-seed empirical characterization. This
    class supplies that characterization and corrects the §4.20 record:
    Schelling's three-state grid {0, 1, 2} (with state 0 = empty) yields
    wall_final_qtr_mean ~0.36 — not 0.02 — because empty-cell boundaries
    are counted as walls by the Moore-neighborhood difference metric.
    The actual rejection mechanism is therefore not wall_final < 0.05
    but the joint failure of the screening and confirmation gates:
      - moran_final_qtr_mean ~0.27 fails the 0.30 screening floor
      - wall_final_qtr_mean ~0.36 fails the 0.30 confirmation ceiling
    on most seeds, with the worst case (seed 2) reaching SCREENING
    tier but stopping there because the wall plateau exceeds the
    confirmation threshold.

    The metadata flag in P18's _check_exclusions for P1 thus serves as
    defense-in-depth (Class 4 belt-and-suspenders) rather than the
    primary discriminator. This affirms §6.10's pure-metric within-
    substrate discrimination class for the Schelling × P18 pair.

    Test parameters: L = 64, n_steps = 300, seeds {0, 1, 2, 3, 4}.
    Schelling at L = 64 with default density 0.9 and threshold 0.375
    converges by ~step 150; n_steps = 300 leaves ample plateau.
    """

    SEEDS = [0, 1, 2, 3, 4]
    L = 64
    N_STEPS = 300

    @pytest.fixture(scope="class")
    def schelling_runs(self):
        """Cache Schelling runs across all 5 seeds once for the class.

        Each Schelling run at L=64, n_steps=300 takes ~12s; without
        caching the parametrized tests below would re-run every seed
        for every test (5 tests × 5 seeds × 12s = 300s). With this
        fixture, each seed runs once (~60s total) and the resulting
        history is shared.
        """
        out = {}
        for seed in self.SEEDS:
            h = run_schelling(grid_size=self.L, n_steps=self.N_STEPS, seed=seed)
            out[seed] = h
        return out

    @pytest.fixture(scope="class")
    def schelling_results(self, schelling_runs):
        """Cache (DetectorResult, trajectory_metrics) for each seed."""
        det = P18ConsensusDetector(n_permutations=199, seed=0)
        out = {}
        for seed, h in schelling_runs.items():
            r = det.detect(h, model_metadata=None)
            m = det._trajectory_metrics(h)
            out[seed] = (r, m)
        return out

    @pytest.fixture(scope="class")
    def schelling_results_with_metadata(self, schelling_runs):
        """Cache results when realistic Schelling metadata is supplied."""
        det = P18ConsensusDetector(n_permutations=199, seed=0)
        meta = {"threshold": 0.375, "density": 0.9}
        out = {}
        for seed, h in schelling_runs.items():
            r = det.detect(h, model_metadata=meta)
            out[seed] = r
        return out

    @pytest.mark.parametrize("seed", SEEDS)
    def test_schelling_does_not_reach_confirmation_without_metadata(
        self, seed, schelling_results
    ):
        """Pure-metric: Schelling must stay below CONFIRMATION across all seeds.

        This is the headline content-level negative test. Schelling is
        passed to P18 with model_metadata=None (no metadata path
        available for any exclusion). The detector must stop at SCREENING
        tier — either by failing one of the screening gates outright or
        by passing screening but failing one of the confirmation gates.
        """
        r, _ = schelling_results[seed]
        assert r.tier != DetectionTier.CONFIRMATION, (
            f"Schelling seed={seed} reached CONFIRMATION (false positive). "
            f"primary={r.primary_metric}, secondary={r.secondary_metrics}, "
            f"p={r.null_p_value:.4f}"
        )
        assert r.tier != DetectionTier.DEFINITIVE, (
            f"Schelling seed={seed} reached DEFINITIVE (false positive). "
            f"primary={r.primary_metric}"
        )

    @pytest.mark.parametrize("seed", SEEDS)
    def test_schelling_does_not_reach_confirmation_with_metadata(
        self, seed, schelling_results_with_metadata
    ):
        """Belt-and-suspenders: same assertion holds with realistic metadata.

        Schelling's registry metadata keys are ``threshold`` and
        ``density``. Neither contains 'copy', 'imitation', or 'voter',
        so the metadata-aided P1 exclusion in _check_exclusions returns
        "inconclusive" rather than "excluded". The pure-metric path is
        unchanged, so confirmation must still not be reached.
        """
        r = schelling_results_with_metadata[seed]
        assert r.tier != DetectionTier.CONFIRMATION, (
            f"Schelling (with metadata) seed={seed} reached CONFIRMATION. "
            f"primary={r.primary_metric}, secondary={r.secondary_metrics}"
        )
        assert r.tier != DetectionTier.DEFINITIVE

    @pytest.mark.parametrize("seed", SEEDS)
    def test_schelling_wall_plateau_corrects_section_4_20_record(
        self, seed, schelling_results
    ):
        """Empirical correction to §4.20 docstring claim of wall ~0.02.

        Schelling's three-state grid {0, 1, 2} yields wall_final_qtr_mean
        in [0.30, 0.45] across L=64 seeds, NOT ~0.02 as §4.20 states.
        The reason: ~10% of cells are empty (state 0) and every empty
        cell adjacent to an occupied cell counts as a wall under the
        Moore-neighborhood difference metric. The DEFINITIVE wall_final
        > 0.05 gate would therefore PASS for Schelling — confirming
        that the actual rejection mechanism is the screening floor on
        moran_final_qtr_mean and the confirmation ceiling on
        wall_final_qtr_mean (0.30), not the definitive wall floor.
        """
        _, m = schelling_results[seed]
        wall_final = m["wall_final_qtr_mean"]
        assert 0.30 <= wall_final <= 0.45, (
            f"Schelling seed={seed}: wall_final_qtr_mean={wall_final:.4f} "
            f"outside expected [0.30, 0.45]. §4.20 record needs update."
        )
        # Confirms the §4.20 claim of "wall ~0.02" is empirically wrong:
        assert wall_final > 0.05, (
            f"Schelling seed={seed}: §4.20 docstring claims wall < 0.05 but "
            f"empirical value is {wall_final:.4f}; the DEFINITIVE wall floor "
            f"would PASS for Schelling."
        )

    @pytest.mark.parametrize("seed", SEEDS)
    def test_schelling_moran_plateau_below_or_near_screening_floor(
        self, seed, schelling_results
    ):
        """Schelling final Moran is in [0.20, 0.32] — at or below the 0.30 floor.

        This is one of the two metric-level discriminators (the other
        is wall_final_qtr_mean exceeding the confirmation ceiling). On
        most seeds, moran_final_qtr_mean < 0.30 fails screening
        outright; on seeds where it edges above (e.g. seed 2 ~0.301),
        the wall plateau still blocks confirmation.
        """
        _, m = schelling_results[seed]
        moran_final = m["moran_final_qtr_mean"]
        assert 0.20 <= moran_final <= 0.32, (
            f"Schelling seed={seed}: moran_final_qtr_mean={moran_final:.4f} "
            f"outside expected [0.20, 0.32]"
        )

    def test_schelling_at_least_four_of_five_seeds_fail_screening(
        self, schelling_results
    ):
        """Aggregate robustness across 5 seeds: ≥4 must fail SCREENING.

        On the Sprint 20 characterization run at L=64, n_steps=300:
          - 4 of 5 seeds (0, 1, 3, 4) fail screening outright with
            detected=False (moran_final < 0.30 floor or moran_growth
            below 0.20).
          - 1 of 5 seeds (seed 2) passes screening with detected=True
            but stops at SCREENING tier because wall_final_qtr_mean
            ~0.36 > 0.30 confirmation ceiling.
        Either way, no seed reaches CONFIRMATION. This aggregate test
        catches drift in the marginal seed boundary.
        """
        screening_failures = 0
        for seed in self.SEEDS:
            r, _ = schelling_results[seed]
            if r.tier == DetectionTier.SCREENING and not r.detected:
                screening_failures += 1
        assert screening_failures >= 4, (
            f"Expected ≥4 of 5 seeds to fail screening outright; "
            f"got {screening_failures}. The marginal seed-2-style edge "
            f"case may have shifted; review the screening-floor calibration."
        )


# ============================================================================
# P18 multi-seed robustness
# ============================================================================


class TestP18MultiSeed:
    """Run-to-run robustness: voter should reach DEFINITIVE across most seeds."""

    @pytest.mark.parametrize("seed", [0, 1, 2, 3, 4])
    def test_voter_seed_reaches_definitive(self, seed):
        det = P18ConsensusDetector(n_permutations=199, seed=0)
        m = VoterModel(rows=64, cols=64, seed=seed)
        h = m.run(n_steps=400)
        r = det.detect(h, model_metadata=m.get_metadata())
        assert r.tier == DetectionTier.DEFINITIVE, \
            f"voter seed={seed} should reach DEFINITIVE, got {r.tier.name}"

    @pytest.mark.parametrize("seed", [0, 1, 2, 3, 4])
    def test_gh_random_does_not_reach_definitive(self, seed):
        det = P18ConsensusDetector(n_permutations=199, seed=0)
        gh = GreenbergHastings(rows=64, cols=64, n_states=8, threshold=2,
                               init_mode="random", seed=seed)
        h = [gh.setup()] + [gh.step() for _ in range(400)]
        r = det.detect(h, model_metadata={"update": "excitable"})
        assert r.tier != DetectionTier.DEFINITIVE

    @pytest.mark.parametrize("seed", [0, 1, 2, 3, 4])
    def test_gol_random_does_not_reach_screening_pass(self, seed):
        det = P18ConsensusDetector(n_permutations=199, seed=0)
        gol = GameOfLife(rows=64, cols=64, init_mode="random", seed=seed)
        h = [gol.setup()] + [gol.step() for _ in range(400)]
        r = det.detect(h, model_metadata={"update": "gol_b3s23"})
        # GoL random fails screening (Moran plateau too low).
        assert r.tier == DetectionTier.SCREENING and not r.detected


# ============================================================================
# Slow tests: finite-size robustness
# ============================================================================


@pytest.mark.slow
class TestSprint20SlowReplication:
    """Finite-size slow tests at L ∈ {64, 128, 256} × 3 seeds.

    Pinned per the Sprint 19 finite-size-as-routine convention. P18
    must reach DEFINITIVE at all three L for the seeds tested. If a
    detector improvement broadens the basin further this test will
    still pass; if a future change degrades performance below the
    L=256 ceiling, this test fails (a regression to fix).
    """

    @pytest.mark.parametrize("L,seed", [
        (64, 0), (64, 1), (64, 42),
        (128, 0), (128, 1), (128, 42),
        (256, 0), (256, 42),
    ])
    def test_voter_definitive_across_L(self, L, seed):
        det = P18ConsensusDetector(n_permutations=199, seed=0)
        # Use a fixed sweep budget that scales modestly with L. At L=64
        # this gives 4× plateau time; at L=256 with 300 sweeps the run
        # is just past plateau onset, which is sufficient for the detector
        # since the early-window primary metric (t ≤ τ = 128) is what
        # matters and the confirmation gates don't need long plateau time.
        n_sweeps = {64: 400, 128: 400, 256: 300}[L]
        m = VoterModel(rows=L, cols=L, seed=seed)
        h = m.run(n_steps=n_sweeps)
        r = det.detect(h, model_metadata=m.get_metadata())
        assert r.tier == DetectionTier.DEFINITIVE, (
            f"voter L={L} seed={seed} should reach DEFINITIVE, "
            f"got {r.tier.name}, p={r.null_p_value:.4f}, "
            f"primary={r.primary_metric}"
        )
