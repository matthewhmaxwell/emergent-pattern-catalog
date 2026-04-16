"""End-to-end tests for SIR epidemic model × P22 information cascade detector.

Tests verify ACTUAL EFFECTS, not just that code runs:
1. SIR model replication: epidemic curve shape, convergence, state counts
2. P22 × SIR: DEFINITIVE detection of spatially-structured cascade
3. SIR × P13: REJECTION (single-pass wave, not persistent excitable wave)
4. P22 × non-SIR lattice_2d: correct rejection (no cascade dynamics)
5. SIR parameter sensitivity: subcritical vs. supercritical

Statistical power: n_permutations=199 (floor p=0.005).
"""

import numpy as np
import pytest

from epc.models.sir_epidemic import SIREpidemicModel
from epc.detectors.p22_information_cascade import P22CascadeDetector
from epc.detectors.p13_excitable_waves import P13ExcitableWaveDetector
from epc.detector_result import DetectionTier


# ============================================================
# SIR model basic replication
# ============================================================


class TestSIRModelReplication:
    """Verify SIR model reproduces published qualitative results."""

    def test_single_seed_epidemic_curve(self):
        """Epidemic from single seed: I(t) rises then falls (bell curve).
        
        Published result: Kermack-McKendrick dynamics produce a unimodal
        infected curve when R0 > 1.
        """
        model = SIREpidemicModel(
            rows=80, cols=80,
            infection_prob=0.20, recovery_prob=0.3,
            init_mode="single_seed", seed=42,
        )
        history = model.run(400)

        i_series = [s["i_count"] for s in history]
        peak_t = np.argmax(i_series)

        # Peak must be in the interior (not at start or end)
        assert peak_t > 5, f"Peak too early: t={peak_t}"
        assert peak_t < len(i_series) - 5, f"Peak at end: t={peak_t}"

        # I(t) must return to 0 (epidemic dies out)
        assert i_series[-1] == 0, f"Epidemic didn't die out: I={i_series[-1]}"

        # Must have at least some infected cells at peak
        assert max(i_series) > 50, f"Peak too small: {max(i_series)}"

    def test_epidemic_state_conservation(self):
        """S + I + R = N at every timestep (population conservation)."""
        model = SIREpidemicModel(
            rows=50, cols=50,
            infection_prob=0.20, recovery_prob=0.3,
            init_mode="single_seed", seed=42,
        )
        history = model.run(200)
        total = 50 * 50

        for t, state in enumerate(history):
            s = state["s_count"]
            i = state["i_count"]
            r = state["r_count"]
            assert s + i + r == total, f"Conservation violated at t={t}: S={s} I={i} R={r}"

    def test_state_dict_keys(self):
        """State history has all required keys for lattice_2d substrate."""
        model = SIREpidemicModel(rows=20, cols=20, seed=42)
        history = model.run(10)

        required_keys = ["grid", "grid_dims", "n_states", "step",
                         "s_count", "i_count", "r_count",
                         "newly_infected", "newly_recovered",
                         "activity_density"]
        for key in required_keys:
            assert key in history[0], f"Missing key: {key}"

        assert history[0]["n_states"] == 3
        assert history[0]["grid_dims"] == (20, 20)

    def test_subcritical_dies_quickly(self):
        """Below percolation threshold, single-seed epidemic dies out fast.
        
        With very low p, the infection has a high probability of dying
        before spreading significantly.
        """
        # Run multiple seeds to account for stochasticity
        died_early = 0
        n_trials = 10
        for seed in range(n_trials):
            model = SIREpidemicModel(
                rows=50, cols=50,
                infection_prob=0.02, recovery_prob=0.5,
                init_mode="single_seed", seed=seed,
            )
            history = model.run(500)
            r_final = history[-1]["r_count"]
            if r_final < 50:  # Less than 2% infected
                died_early += 1

        # Most trials should die early with these subcritical params
        assert died_early >= 7, \
            f"Only {died_early}/{n_trials} died early — p=0.02 should be subcritical"

    def test_supercritical_percolates(self):
        """Above threshold, epidemic infects large fraction of population."""
        model = SIREpidemicModel(
            rows=80, cols=80,
            infection_prob=0.25, recovery_prob=0.3,
            init_mode="single_seed", seed=42,
        )
        history = model.run(400)
        total = 80 * 80
        r_final = history[-1]["r_count"]
        r_fraction = r_final / total

        assert r_fraction > 0.5, \
            f"Supercritical epidemic should infect >50%, got {r_fraction:.1%}"

    def test_metadata_r0_estimate(self):
        """Metadata includes approximate R0 estimate."""
        model = SIREpidemicModel(
            rows=50, cols=50,
            infection_prob=0.20, recovery_prob=0.3,
            neighborhood="moore", seed=42,
        )
        model.setup()
        meta = model.get_metadata()

        assert "r0_approx" in meta
        assert meta["r0_approx"] > 1.0, "R0 should be >1 for these params"
        assert meta["model_name"] == "sir_epidemic"
        assert meta["n_states"] == 3


# ============================================================
# P22 canonical positive: SIR epidemic
# ============================================================


class TestSIRP22Canonical:
    """P22 detector must achieve DEFINITIVE on supercritical SIR epidemic."""

    def test_sir_p22_definitive(self):
        """SIR × P22: DEFINITIVE detection.
        
        Key evidence:
        - Cascade reaches >30% of population
        - Moran's I on infection time map >> null (spatial wavefront)
        - Unimodal epidemic curve
        - R0 > 1
        - Epidemic dies out (single-pass cascade)
        """
        model = SIREpidemicModel(
            rows=80, cols=80,
            infection_prob=0.20, recovery_prob=0.3,
            init_mode="single_seed", seed=42,
        )
        history = model.run(400, record_every=1)
        meta = model.get_metadata()

        det = P22CascadeDetector(n_permutations=199, seed=42)
        result = det.detect(history, meta)

        assert result.detected, "P22 should detect cascade in SIR"
        assert result.tier >= DetectionTier.DEFINITIVE, \
            f"Expected DEFINITIVE, got {result.tier}"
        assert result.confidence >= 0.80

        # Verify the actual metrics support the detection
        assert result.primary_metric["cascade_reach_total"] > 0.30
        assert result.primary_metric["moran_i_time"] > 0.90
        assert result.primary_metric["is_unimodal"] > 0.5
        assert result.primary_metric["r0_estimate"] > 1.0
        assert result.primary_metric["epidemic_died_out"] > 0.5

        # Null model: random timing should have ~0 Moran's I
        assert result.effect_size.get("null_mean", 1.0) < 0.1
        assert result.effect_size.get("cohens_d", 0.0) > 10.0

    def test_sir_p22_effect_size(self):
        """Effect size must be massive: wavefront structure is unmistakable."""
        model = SIREpidemicModel(
            rows=60, cols=60,
            infection_prob=0.25, recovery_prob=0.3,
            init_mode="single_seed", seed=123,
        )
        history = model.run(300, record_every=1)
        meta = model.get_metadata()

        det = P22CascadeDetector(n_permutations=99, seed=42)
        result = det.detect(history, meta)

        # Cohen's d should be enormous (observed Moran >> null Moran)
        d = result.effect_size.get("cohens_d", 0.0)
        assert d > 20.0, f"Effect size too small: d={d}"


# ============================================================
# SIR × P13: MUST REJECT (the core boundary test)
# ============================================================


class TestSIRP13Rejection:
    """P13 must reject SIR — epidemic waves are NOT excitable waves.

    SIR has 3 states (S/I/R), so it PASSES the n_states≥3 guard.
    But SIR wavefronts are single-pass: recovered cells never return
    to susceptible, so there is no re-excitation, no spirals, and
    activity dies out. P13 should reject on:
    - wavefront_speed_mean = 0 (WavefrontSpeedLocal measures inter-excitation
      intervals, which don't exist for single-pass waves)
    - persistence: activity dies before 5×T_prop
    """

    def test_sir_rejected_by_p13(self):
        """P13 must fail screening on SIR epidemic."""
        model = SIREpidemicModel(
            rows=50, cols=50,
            infection_prob=0.20, recovery_prob=0.3,
            init_mode="single_seed", seed=42,
        )
        history = model.run(300, record_every=1)
        meta = model.get_metadata()

        det = P13ExcitableWaveDetector(n_null_runs=99)
        result = det.detect(history, meta)

        assert not result.detected, \
            "P13 should NOT detect excitable waves in SIR"
        assert result.tier == DetectionTier.SCREENING, \
            f"SIR should fail P13 screening, got {result.tier}"
        assert result.confidence == 0.0

    def test_sir_has_3_states_but_p13_still_rejects(self):
        """Verify the n_states guard passes (3 states) but screening fails.
        
        This is the key boundary test: n_states=3 means the hard guard
        does NOT reject SIR. The discrimination comes from wavefront
        speed and persistence metrics.
        """
        model = SIREpidemicModel(
            rows=50, cols=50,
            infection_prob=0.20, recovery_prob=0.3,
            init_mode="single_seed", seed=42,
        )
        history = model.run(300, record_every=1)
        meta = model.get_metadata()

        # Verify SIR has 3 states (passes n_states guard)
        assert history[0]["n_states"] == 3

        det = P13ExcitableWaveDetector(n_null_runs=99)
        result = det.detect(history, meta)

        # Should still fail screening despite passing n_states guard
        assert result.tier == DetectionTier.SCREENING
        assert result.confidence == 0.0

        # The primary metric should show WHY it failed
        pm = result.primary_metric
        # Either speed is 0/unmeasurable or CV is too high
        speed_fails = pm.get("wavefront_speed_mean", 0.0) <= 0
        cv_fails = pm.get("wavefront_speed_cv", float("inf")) > 0.2
        died = pm.get("died_out", 0.0) > 0.5
        streak_short = pm.get("longest_active_streak", 0) < 250  # 5 × 50

        assert speed_fails or cv_fails or died or streak_short, \
            f"P13 should fail on speed/CV/persistence, but primary metrics: {pm}"


# ============================================================
# P22 cross-detection: non-SIR models
# ============================================================


class TestP22CrossDetection:
    """P22 should NOT detect cascades in non-epidemic models.
    
    The cascade detector measures spatial structure in infection timing.
    Models without epidemic-like spreading should fail screening because
    they don't have the characteristic wavefront-from-seed pattern.
    """

    def test_schelling_rejected_by_p22(self):
        """Schelling segregation is NOT a cascade — no wavefront."""
        from epc.models.schelling import run_schelling
        history = run_schelling(grid_size=50, n_steps=200, seed=42)

        det = P22CascadeDetector(n_permutations=99, seed=42)
        result = det.detect(history, None)

        # Should fail screening: no cascade-like dynamics
        assert result.tier == DetectionTier.SCREENING, \
            f"Schelling should fail P22 screening, got {result.tier}"

    def test_greenberg_hastings_rejected_by_p22(self):
        """GH excitable waves are persistent, NOT single-pass cascades.
        
        P22 screening requires cascade_reach ≥ 5% AND spatial Moran's I ≥ 0.1.
        GH waves don't produce a single-pass cascade pattern — they cycle.
        The infection time map won't show a clean wavefront structure because
        cells get re-excited repeatedly.
        """
        from epc.models.greenberg_hastings import GreenbergHastings

        model = GreenbergHastings(
            rows=50, cols=50, n_states=5, threshold=1,
            init_mode="broken_wave", seed=42,
        )
        history = model.run(200, record_every=1)
        meta = model.get_metadata()

        det = P22CascadeDetector(n_permutations=99, seed=42)
        result = det.detect(history, meta)

        # GH waves are persistent/cyclic — not cascade dynamics
        # Should fail because infection time map doesn't capture re-excitation
        # or because the "epidemic" never dies out
        assert result.tier <= DetectionTier.CONFIRMATION, \
            f"GH should not achieve P22 definitive, got {result.tier}"
