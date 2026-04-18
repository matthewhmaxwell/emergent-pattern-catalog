"""P3 Turing-wavelength detector: end-to-end tests.

Verifies that P3 correctly identifies stationary periodic spatial
patterns from reaction-diffusion instability (Gray-Scott), and correctly
rejects every existing integer-grid model — even when those models
produce high radial-FFT peak-to-mean ratios (notably RPS at low mobility
where p/m ~ 23 matches Gray-Scott labyrinth numerically).

Design principle (Sprint 13): discrimination is by SUBSTRATE TYPE
(continuous-valued field vs integer-labelled grid), not by empirical
threshold tuning. The n_unique_values >= 50 prerequisite cleanly
separates Gray-Scott (~16k unique float values) from every integer-grid
model (<= 10 unique values).

Each test verifies the ACTUAL EFFECT, not just that code runs:
  - Gray-Scott labyrinth (F=0.037, k=0.060) MUST reach DEFINITIVE across
    seeds 42, 7, 123.
  - Gray-Scott spots (F=0.030, k=0.062) MUST reach at least CONFIRMATION.
  - Gray-Scott uniform decay regime (F=0.10, k=0.10) MUST reject at
    prereq (field_std = 0).
  - All integer-grid models (RPS, GH, Schelling, Nowak-May, GoL, SIR,
    LV) MUST reject at the substrate prerequisite with informative
    warnings — they don't expose a 'field' observable.
  - Adversarial cast: an RPS grid manually relabelled as 'field' MUST
    reject at n_unique_values prerequisite (nu = 4 << 50).

Timing: the slowest test is the N=128 T=4000 canonical positive at
~3-5 seconds per seed. Three seeds = ~15s. Negative sweep trivial.
"""

from __future__ import annotations

import numpy as np
import pytest

from epc.detector_result import DetectionTier
from epc.detectors.p3_turing_wavelength import P3TuringWavelengthDetector
from epc.models.gray_scott import GrayScott
from epc.models.greenberg_hastings import GreenbergHastings
from epc.models.lotka_volterra_lattice import LotkaVolterraLattice
from epc.models.nowak_may import NowakMayModel
from epc.models.rps_spatial import RPSSpatialModel
from epc.models.schelling import run_schelling
from epc.models.sir_epidemic import SIREpidemicModel
from epc.models.game_of_life import GameOfLife


# ===================================================================
# Canonical positive: Gray-Scott labyrinth
# ===================================================================


class TestGrayScottLabyrinthP3Definitive:
    """Canonical P3 positive: Gray-Scott labyrinth (F=0.037, k=0.060),
    N=128, T>=4000. Must reach DEFINITIVE across multiple seeds.
    """

    @pytest.mark.parametrize("seed", [42, 7, 123])
    def test_labyrinth_reaches_definitive(self, seed):
        m = GrayScott(
            rows=128, cols=128,
            feed_rate=0.037, kill_rate=0.060,
            seed=seed,
        )
        history = m.run(n_steps=4000)

        det = P3TuringWavelengthDetector(n_permutations=199, seed=42)
        result = det.detect(history, m.get_metadata())

        assert result.detected, (
            f"GS labyrinth seed={seed} should be detected; got tier={result.tier.name}, "
            f"primary={result.primary_metric}"
        )
        assert result.tier == DetectionTier.DEFINITIVE, (
            f"GS labyrinth seed={seed} should reach DEFINITIVE; got {result.tier.name}. "
            f"p/m={result.primary_metric.get('peak_to_mean')}, "
            f"d={result.effect_size.get('cohens_d')}, "
            f"p={result.null_p_value}, "
            f"cv={result.secondary_metrics.get('peak_k_cv')}"
        )

    def test_labyrinth_wavelength_physical(self):
        """Wavelength should be ~12 pixels (Sprint 13 characterization
        finding; invariant across grid sizes N=64..192)."""
        m = GrayScott(rows=128, cols=128, feed_rate=0.037, kill_rate=0.060, seed=42)
        history = m.run(4000)
        det = P3TuringWavelengthDetector(n_permutations=199, seed=42)
        result = det.detect(history, m.get_metadata())

        lam = result.primary_metric["wavelength_pixels"]
        assert 10.0 <= lam <= 16.0, (
            f"Labyrinth wavelength should be ~12 px; got {lam}"
        )

    def test_labyrinth_peak_k_stable(self):
        """peak_k must be stable across late-trajectory snapshots."""
        m = GrayScott(rows=128, cols=128, feed_rate=0.037, kill_rate=0.060, seed=42)
        history = m.run(4000)
        det = P3TuringWavelengthDetector(n_permutations=199, seed=42)
        result = det.detect(history, m.get_metadata())

        cv = result.secondary_metrics.get("peak_k_cv", 999)
        assert cv < 0.05, (
            f"Labyrinth peak_k should be highly stable; cv={cv}"
        )

    def test_labyrinth_strong_effect_size(self):
        """Cohen's d vs shuffle null must be enormous (characterization
        showed d ~100). Anything < 30 suggests the null is wrong or the
        signal isn't Turing."""
        m = GrayScott(rows=128, cols=128, feed_rate=0.037, kill_rate=0.060, seed=42)
        history = m.run(4000)
        det = P3TuringWavelengthDetector(n_permutations=199, seed=42)
        result = det.detect(history, m.get_metadata())

        d = result.effect_size.get("cohens_d", 0)
        assert d > 30.0, (
            f"Labyrinth d vs shuffle null should be >> 30; got {d}"
        )

    def test_labyrinth_exclusions_cleared(self):
        """P1 and P13 should both be marked 'excluded' (continuous-field
        substrate excludes their canonical positives by construction)."""
        m = GrayScott(rows=128, cols=128, feed_rate=0.037, kill_rate=0.060, seed=42)
        history = m.run(4000)
        det = P3TuringWavelengthDetector(n_permutations=199, seed=42)
        result = det.detect(history, m.get_metadata())

        assert result.exclusion_results.get("P1") == "excluded"
        assert result.exclusion_results.get("P13") == "excluded"


# ===================================================================
# Canonical positive: Gray-Scott spots
# ===================================================================


class TestGrayScottSpotsP3Confirmation:
    """Gray-Scott spots regime (F=0.030, k=0.062) reaches at least
    CONFIRMATION — p/m is lower than labyrinth (~13 vs ~19)."""

    def test_spots_reaches_confirmation(self):
        m = GrayScott(rows=128, cols=128, feed_rate=0.030, kill_rate=0.062, seed=42)
        history = m.run(4000)
        det = P3TuringWavelengthDetector(n_permutations=199, seed=42)
        result = det.detect(history, m.get_metadata())

        assert result.detected
        assert result.tier >= DetectionTier.CONFIRMATION, (
            f"GS spots should reach >= CONFIRMATION; got {result.tier.name}. "
            f"p/m={result.primary_metric.get('peak_to_mean')}, "
            f"d={result.effect_size.get('cohens_d')}"
        )


# ===================================================================
# Canonical negative: non-Turing Gray-Scott regimes
# ===================================================================


class TestNonTuringGrayScottRegimesRejected:
    """Gray-Scott parameters OUTSIDE the Turing window must be rejected.
    These test the prereq_gate, not the substrate distinction."""

    def test_uniform_decay_rejected_by_field_std_prereq(self):
        """(F=0.10, k=0.10): v decays to 0 uniformly. field_std = 0 should
        fail the prereq."""
        m = GrayScott(rows=64, cols=64, feed_rate=0.10, kill_rate=0.10, seed=42)
        history = m.run(2000)
        det = P3TuringWavelengthDetector(n_permutations=99, seed=42)
        result = det.detect(history, m.get_metadata())

        assert not result.detected, (
            f"Uniform-decay GS should not be detected; got tier={result.tier.name}"
        )
        assert result.primary_metric["field_std"] < 0.01, (
            f"Uniform decay should have near-zero std; got "
            f"{result.primary_metric['field_std']}"
        )


# ===================================================================
# Substrate-prerequisite rejections: integer-grid models
# ===================================================================


class TestIntegerGridModelsRejectedBySubstrate:
    """All integer-grid models must reject at the substrate prerequisite.

    They don't expose a 'field' observable — their state carries a 'grid'
    observable instead. P3 must NOT consume 'grid' (would be a silent
    substrate violation). A substrate warning must be emitted.
    """

    def _check_rejected_with_substrate_warning(self, history, metadata):
        det = P3TuringWavelengthDetector(n_permutations=99, seed=42)
        result = det.detect(history, metadata)
        assert not result.detected, (
            f"Integer-grid model should reject; got tier={result.tier.name}"
        )
        assert any(
            "continuous" in w.lower() or "field" in w.lower()
            for w in result.warnings
        ), (
            f"Expected substrate warning; got warnings={result.warnings}"
        )

    def test_rps_rejected(self):
        """RPS at low mobility has radial-FFT p/m ~23 on the raw grid —
        this is the critical false-positive risk P3 must handle."""
        m = RPSSpatialModel(rows=100, cols=100, mobility=1e-4, seed=42)
        history = m.run(1200)
        self._check_rejected_with_substrate_warning(history, m.get_metadata())

    def test_gh_rejected(self):
        m = GreenbergHastings(rows=100, cols=100, n_states=8, threshold=2,
                              init_mode="random", seed=42)
        history = m.run(200)
        self._check_rejected_with_substrate_warning(history, m.get_metadata())

    def test_lv_rejected(self):
        m = LotkaVolterraLattice(
            rows=80, cols=80,
            predation_rate=4.0,
            prey_reproduction_rate=1.0,
            predator_death_rate=1.0,
            seed=42,
        )
        history = m.run(500)
        self._check_rejected_with_substrate_warning(history, m.get_metadata())

    def test_schelling_rejected(self):
        history = run_schelling(grid_size=60, density=0.9, threshold=0.375,
                                 n_steps=50, seed=42)
        self._check_rejected_with_substrate_warning(
            history, {"model_name": "schelling"}
        )

    def test_nowak_may_rejected(self):
        m = NowakMayModel(rows=80, cols=80, b=1.8, seed=42)
        history = m.run(100)
        self._check_rejected_with_substrate_warning(history, m.get_metadata())

    def test_gol_rejected(self):
        m = GameOfLife(rows=80, cols=80, init_mode="random", init_density=0.37,
                       seed=42)
        history = m.run(150)
        self._check_rejected_with_substrate_warning(history, m.get_metadata())

    def test_sir_rejected(self):
        m = SIREpidemicModel(
            rows=80, cols=80,
            infection_prob=0.3, recovery_prob=0.1,
            init_mode="single_seed", seed=42,
        )
        history = m.run(80)
        self._check_rejected_with_substrate_warning(history, m.get_metadata())


# ===================================================================
# Adversarial rejection: n_unique_values prerequisite
# ===================================================================


class TestAdversarialDiscreteFieldRejected:
    """If an integer-grid model's grid is manually re-labelled as 'field',
    P3 must still reject via the n_unique_values prerequisite.

    This exercises the core discrimination mechanism: substrate detection
    is content-based (n_unique_values), not just key-based. RPS raw grid
    has p/m ≈ 23 — higher than Gray-Scott labyrinth — so ANY detector
    gating only on peak-to-mean would false-positive.
    """

    def test_rps_with_field_alias_rejected(self):
        """Adversarial: inject RPS grid as 'field' key."""
        m = RPSSpatialModel(rows=100, cols=100, mobility=1e-4, seed=42)
        history = m.run(1200)

        adversarial = []
        for snap in history:
            s = dict(snap)
            s["field"] = snap["grid"].astype(float)
            adversarial.append(s)

        det = P3TuringWavelengthDetector(n_permutations=99, seed=42)
        result = det.detect(adversarial, m.get_metadata())

        assert not result.detected, (
            f"Adversarial RPS-as-field should reject; got tier={result.tier.name}. "
            f"n_unique={result.primary_metric.get('n_unique_values')}, "
            f"p/m={result.primary_metric.get('peak_to_mean')}"
        )
        # nu should be small (4 for 4-state RPS)
        assert result.primary_metric["n_unique_values"] < 50, (
            f"RPS should have nu < 50; got "
            f"{result.primary_metric['n_unique_values']}"
        )

    def test_binary_grid_rejected(self):
        """Adversarial: binary (0/1) field masquerading as continuous."""
        rng = np.random.default_rng(42)
        grid = (rng.random((64, 64)) < 0.5).astype(float)
        adversarial = [
            {"field": grid, "grid_dims": grid.shape, "step": 0},
            {"field": grid, "grid_dims": grid.shape, "step": 1},
        ]
        det = P3TuringWavelengthDetector(n_permutations=99, seed=42)
        result = det.detect(adversarial, {"model_name": "fake"})
        assert not result.detected, (
            f"Binary field should reject on nu prereq; got tier={result.tier.name}"
        )


# ===================================================================
# Scaling: peak_k scales with N, wavelength is invariant
# ===================================================================


class TestTuringWavelengthScalesWithGrid:
    """The Turing wavelength is a physical property of the system, not
    a grid artifact. peak_k must scale linearly with N while
    wavelength in pixels remains invariant.

    Verified in Sprint 13 characterization: peak_k = 5 at N=64,
    peak_k = 7 at N=96, peak_k = 10 at N=128, peak_k = 16 at N=192.
    All give lambda ≈ 12 px.
    """

    @pytest.mark.parametrize("N,expected_k,expected_tier", [
        (64, 5, DetectionTier.DEFINITIVE),
        (96, 7, DetectionTier.DEFINITIVE),
        (128, 10, DetectionTier.DEFINITIVE),
    ])
    def test_peak_k_scaling(self, N, expected_k, expected_tier):
        # At smaller N the pattern selects faster. T=3000 suffices for N<=96,
        # T=4000 for N=128.
        T = 3000 if N <= 96 else 4000
        m = GrayScott(rows=N, cols=N, feed_rate=0.037, kill_rate=0.060, seed=42)
        history = m.run(T)
        det = P3TuringWavelengthDetector(n_permutations=199, seed=42)
        result = det.detect(history, m.get_metadata())

        assert result.detected
        # Allow ±2 on peak_k due to radial binning discretization
        peak_k = result.primary_metric["peak_k"]
        assert abs(peak_k - expected_k) <= 2, (
            f"N={N} peak_k should be ~{expected_k}; got {peak_k}"
        )
        # Wavelength in pixels must be ~12 px regardless of N
        lam = result.primary_metric["wavelength_pixels"]
        assert 10.0 <= lam <= 16.0, (
            f"N={N} wavelength should be ~12 px; got {lam}"
        )
        assert result.tier >= expected_tier, (
            f"N={N} should reach >= {expected_tier.name}; got {result.tier.name}"
        )


# ===================================================================
# Sanity: empty / too-short histories handled gracefully
# ===================================================================


class TestEdgeCases:
    """Degenerate inputs must not crash."""

    def test_empty_history(self):
        det = P3TuringWavelengthDetector(seed=42)
        result = det.detect([], None)
        assert not result.detected
        assert any("empty" in w.lower() for w in result.warnings)

    def test_no_field_observable(self):
        det = P3TuringWavelengthDetector(seed=42)
        fake_history = [
            {"positions": np.zeros((10, 2)), "step": 0},
            {"positions": np.zeros((10, 2)), "step": 1},
        ]
        result = det.detect(fake_history, None)
        assert not result.detected
        assert any("field" in w.lower() for w in result.warnings)

    def test_constant_field_rejected(self):
        """Uniform field (std=0) must reject at prereq."""
        flat = np.full((64, 64), 0.5)
        history = [
            {"field": flat, "grid_dims": flat.shape, "step": i}
            for i in range(200)
        ]
        det = P3TuringWavelengthDetector(seed=42)
        result = det.detect(history, {"model_name": "flat"})
        assert not result.detected
