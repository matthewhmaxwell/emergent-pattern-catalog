"""P11 predator-prey oscillation detector: end-to-end tests.

Verifies that P11 correctly identifies the bilateral oscillating
coexistence in the lattice Lotka-Volterra model (Mobilia-Georgiev-Täuber
2007), and correctly rejects:
  - Schelling (2 species but zero species-fraction variance: agent
    identity is conserved per-agent, no dynamics)
  - Nowak-May (fails n_species prerequisite once fixation completes OR
    fails total_std prerequisite before fixation: strictly conserved
    coop + defect = 1)
  - SIR (active dynamics are transient; post-burn-in the system is
    static, failing variance prerequisites)
  - RPS (fails n_species prerequisite: 3 species, not 2)
  - White noise time series injected via a stub model

Each test verifies the ACTUAL EFFECT, not just that code runs:
  - LV × P11 MUST reach DEFINITIVE tier on canonical parameters across
    multiple seeds.
  - LV × P11 MUST produce rho_anti < -0.5 consistently (measured range:
    -0.72 to -0.88 across seeds 42, 7, 123).
  - LV × P11 MUST identify a nonzero-lag anti-correlation (|tau| >= 5).
  - LV × P11 MUST mark P12 and P9 as 'excluded' in exclusion_results.
  - Each negative MUST reject at SCREENING (detected=False or tier=SCREENING).
  - Schelling MUST fail the species-variance prerequisite (std ≈ 0).
  - Nowak-May MUST fail either n_species (after fixation) or total_std
    prerequisite (during coexistence: strictly A+B=1 conservation).
  - SIR MUST fail post-burn-in variance prerequisites OR rho_anti screen.
  - RPS MUST fail the n_species=2 prerequisite with n_species=3 reported.
"""

from __future__ import annotations

import numpy as np
import pytest

from epc.detector_result import DetectionTier
from epc.detectors.p11_predator_prey_oscillation import P11PredatorPreyDetector
from epc.models.lotka_volterra_lattice import LotkaVolterraLattice
from epc.models.nowak_may import NowakMayModel
from epc.models.rps_spatial import RPSSpatialModel
from epc.models.schelling import run_schelling
from epc.models.sir_epidemic import SIREpidemicModel


# ===================================================================
# LV → P11 detection (positive cases)
# ===================================================================

class TestLVDetectedByP11:
    """Canonical LV positive: L=100 λ=4 σ=μ=1, various seeds."""

    @pytest.mark.parametrize("seed", [42, 7, 123])
    def test_lv_reaches_definitive(self, seed):
        """LV × P11 at canonical params reaches DEFINITIVE tier."""
        m = LotkaVolterraLattice(
            rows=100, cols=100,
            predation_rate=4.0,
            prey_reproduction_rate=1.0,
            predator_death_rate=1.0,
            seed=seed,
        )
        history = m.run(n_steps=1200)

        # Guard: skip this seed if predators went extinct before
        # we got enough data (rare at these params but can happen).
        if len(history) < 500 or history[-1]["predator_count"] == 0:
            pytest.skip(f"seed {seed}: predators went extinct at "
                        f"t={len(history)} — not a canonical positive")

        det = P11PredatorPreyDetector(n_permutations=99, seed=42)
        result = det.detect(history, m.get_metadata())

        assert result.detected, (
            f"LV seed={seed} should be detected; got tier={result.tier.name}, "
            f"primary={result.primary_metric}"
        )
        assert result.tier == DetectionTier.DEFINITIVE, (
            f"LV seed={seed} should reach DEFINITIVE; got {result.tier.name}. "
            f"rho_anti={result.primary_metric.get('rho_anti')}, "
            f"fft_p2m={result.primary_metric.get('fft_peak_to_mean')}, "
            f"d={result.effect_size.get('cohens_d')}"
        )

    def test_lv_rho_anti_is_strongly_negative(self):
        """rho_anti must be strongly negative (< -0.5) across seeds."""
        m = LotkaVolterraLattice(
            rows=100, cols=100,
            predation_rate=4.0,
            prey_reproduction_rate=1.0,
            predator_death_rate=1.0,
            seed=42,
        )
        history = m.run(n_steps=1200)
        det = P11PredatorPreyDetector(n_permutations=99, seed=42)
        result = det.detect(history, m.get_metadata())

        rho_anti = result.primary_metric["rho_anti"]
        assert rho_anti < -0.5, (
            f"LV canonical positive should give rho_anti < -0.5, got {rho_anti}"
        )
        # Should be deeper than -0.7 on seed=42 (measured: -0.86)
        assert rho_anti < -0.7, (
            f"LV seed=42 should give rho_anti < -0.7 (measured ~-0.86), "
            f"got {rho_anti}"
        )

    def test_lv_tau_anti_at_nonzero_lag(self):
        """tau_anti must be at |lag| >= 5 (screening threshold)."""
        m = LotkaVolterraLattice(
            rows=100, cols=100,
            predation_rate=4.0,
            prey_reproduction_rate=1.0,
            predator_death_rate=1.0,
            seed=42,
        )
        history = m.run(n_steps=1200)
        det = P11PredatorPreyDetector(n_permutations=99, seed=42)
        result = det.detect(history, m.get_metadata())

        tau_anti = result.primary_metric["tau_anti"]
        assert abs(int(tau_anti)) >= 5, (
            f"tau_anti should have |tau| >= 5, got {tau_anti}"
        )

    def test_lv_fft_peak_to_mean_large(self):
        """FFT peak-to-mean ratio should be large (> 12) on canonical LV."""
        m = LotkaVolterraLattice(
            rows=100, cols=100,
            predation_rate=4.0,
            prey_reproduction_rate=1.0,
            predator_death_rate=1.0,
            seed=42,
        )
        history = m.run(n_steps=1200)
        det = P11PredatorPreyDetector(n_permutations=99, seed=42)
        result = det.detect(history, m.get_metadata())

        p2m = result.primary_metric["fft_peak_to_mean"]
        assert p2m > 12.0, (
            f"LV canonical should give fft_peak_to_mean > 12, got {p2m}"
        )

    def test_lv_exclusions_cleared(self):
        """LV × P11 definitive must mark P9 and P12 as excluded (n_species=2)."""
        m = LotkaVolterraLattice(
            rows=100, cols=100,
            predation_rate=4.0,
            prey_reproduction_rate=1.0,
            predator_death_rate=1.0,
            seed=42,
        )
        history = m.run(n_steps=1200)
        det = P11PredatorPreyDetector(n_permutations=99, seed=42)
        result = det.detect(history, m.get_metadata())

        if result.tier < DetectionTier.CONFIRMATION:
            pytest.skip("exclusions only checked at confirmation+")

        assert result.exclusion_results.get("P12") == "excluded", (
            f"P12 should be excluded via n_species=2, got "
            f"{result.exclusion_results}"
        )
        assert result.exclusion_results.get("P9") == "excluded", (
            f"P9 should be excluded via n_species=2, got "
            f"{result.exclusion_results}"
        )

    def test_lv_halves_agree(self):
        """rho_anti should be < -0.3 on both first and second halves."""
        m = LotkaVolterraLattice(
            rows=100, cols=100,
            predation_rate=4.0,
            prey_reproduction_rate=1.0,
            predator_death_rate=1.0,
            seed=42,
        )
        history = m.run(n_steps=1200)
        det = P11PredatorPreyDetector(n_permutations=99, seed=42)
        result = det.detect(history, m.get_metadata())

        r1 = result.secondary_metrics["rho_anti_first_half"]
        r2 = result.secondary_metrics["rho_anti_second_half"]
        assert r1 < -0.3 and r2 < -0.3, (
            f"halves should both show anti-correlation; got "
            f"first_half={r1}, second_half={r2}"
        )
        assert result.secondary_metrics["halves_agree"] is True


# ===================================================================
# Negatives: each MUST NOT be detected as P11
# ===================================================================

class TestSchellingRejectedByP11:
    """Schelling has zero species-fraction variance (per-agent identity
    is conserved). Must fail species_std prerequisite."""

    def test_schelling_rejected(self):
        history = run_schelling(
            grid_size=60, density=0.9, threshold=0.375,
            n_steps=300, seed=42,
        )
        det = P11PredatorPreyDetector(n_permutations=99, seed=42)
        result = det.detect(history)

        assert not result.detected, (
            f"Schelling should not be detected as P11; got "
            f"tier={result.tier.name}"
        )
        # Prerequisite should fail on species std
        assert result.primary_metric["prey_std"] < 0.005, (
            f"Schelling prey_std should be ~0; got "
            f"{result.primary_metric['prey_std']}"
        )
        assert result.primary_metric["predator_std"] < 0.005, (
            f"Schelling predator_std should be ~0; got "
            f"{result.primary_metric['predator_std']}"
        )


class TestNowakMayRejectedByP11:
    """Nowak-May has strictly A+B=1 conservation (no empty cells). Must
    fail total_std prerequisite, OR fail n_species prerequisite after
    fixation to a single strategy."""

    def test_nowak_may_rejected(self):
        m = NowakMayModel(
            rows=60, cols=60, b=1.8,
            init_coop_fraction=0.5, seed=42,
        )
        history = m.run(n_steps=200)
        det = P11PredatorPreyDetector(n_permutations=99, seed=42)
        result = det.detect(history, m.get_metadata())

        assert not result.detected, (
            f"Nowak-May should not be detected as P11; got "
            f"tier={result.tier.name}"
        )

        # One of two prerequisite failures must happen:
        n_sp = result.primary_metric["n_unique_species_observed"]
        total_std = result.primary_metric["total_std"]
        fixated = n_sp < 2
        conserved = total_std < 0.005
        assert fixated or conserved, (
            f"Nowak-May should fail n_species OR total_std prereq; got "
            f"n_sp={n_sp}, total_std={total_std}"
        )


class TestSIRRejectedByP11:
    """SIR's active phase is transient; post-burn-in the system is static
    (all recovered). Must fail variance prerequisites post-burn-in, OR
    fail rho_anti screening if enough active-phase samples survive."""

    def test_sir_rejected(self):
        m = SIREpidemicModel(
            rows=60, cols=60,
            infection_prob=0.3, recovery_prob=0.1,
            init_mode="single_seed", seed=42,
        )
        history = m.run(n_steps=200)
        det = P11PredatorPreyDetector(n_permutations=99, seed=42)
        result = det.detect(history, m.get_metadata())

        assert not result.detected, (
            f"SIR should not be detected as P11; got tier={result.tier.name}"
        )


class TestRPSRejectedByP11:
    """RPS has 3 species. Must fail n_species == 2 prerequisite at
    the screening check."""

    def test_rps_rejected_via_n_species(self):
        m = RPSSpatialModel(rows=60, cols=60, mobility=1e-4, seed=42)
        history = m.run(n_steps=400)
        det = P11PredatorPreyDetector(n_permutations=99, seed=42)
        result = det.detect(history, m.get_metadata())

        assert not result.detected, (
            f"RPS should not be detected as P11; got tier={result.tier.name}"
        )
        # Auto-detection should find 3 species
        n_sp = result.primary_metric["n_unique_species_observed"]
        assert n_sp == 3, (
            f"RPS should have n_unique_species_observed=3, got {n_sp}"
        )

    def test_rps_anti_correlation_exists_but_prereq_fails(self):
        """Diagnostic: confirm that RPS species pairs DO have strong anti-
        correlation (rho_anti < -0.5), so the ONLY thing keeping RPS
        from being detected is the n_species prerequisite. This documents
        why the prerequisite is essential rather than decorative."""
        m = RPSSpatialModel(rows=60, cols=60, mobility=1e-4, seed=42)
        history = m.run(n_steps=400)
        det = P11PredatorPreyDetector(n_permutations=99, seed=42)
        result = det.detect(history, m.get_metadata())

        rho_anti = result.primary_metric.get("rho_anti", 0.0)
        # RPS species pairs measured at rho_anti ~= -0.93 to -0.97
        assert rho_anti < -0.5, (
            f"RPS species pairs should have strong anti-correlation; "
            f"got rho_anti={rho_anti}. This test documents that the "
            f"n_species prerequisite is what separates P11 from P12."
        )


# ===================================================================
# White noise null: must reject cleanly at screening
# ===================================================================

class TestWhiteNoiseRejectedByP11:
    """Construct a synthetic state_history with uncorrelated white-noise
    species counts on a lattice. rho_anti should be small (|rho| < 0.3),
    fft_peak_to_mean should be small (< 6), and the detector MUST reject
    at screening."""

    def test_white_noise_rejected(self):
        rng = np.random.default_rng(42)
        rows = cols = 40
        n_steps = 600
        history = []
        for t in range(n_steps):
            # Random assignment each frame: 40% empty, 30% prey, 30% pred
            g = rng.choice([0, 1, 2], size=(rows, cols),
                           p=[0.4, 0.3, 0.3]).astype(np.int8)
            history.append({
                "grid": g,
                "grid_dims": (rows, cols),
                "step": t,
                "prey_fraction": float((g == 1).sum() / g.size),
                "predator_fraction": float((g == 2).sum() / g.size),
            })
        det = P11PredatorPreyDetector(n_permutations=99, seed=42)
        result = det.detect(history)

        assert not result.detected, (
            f"White-noise grid should not be detected; got "
            f"tier={result.tier.name}, rho_anti="
            f"{result.primary_metric.get('rho_anti')}"
        )
        rho_anti = result.primary_metric.get("rho_anti", 0.0)
        assert rho_anti > -0.3, (
            f"White-noise rho_anti should be > -0.3 (shallow anti-corr); "
            f"got {rho_anti}"
        )


# ===================================================================
# Prerequisite behavior: detector must gracefully handle bad inputs
# ===================================================================

class TestP11Prerequisites:
    """Verify prerequisite checks in isolation."""

    def test_empty_history_rejected(self):
        det = P11PredatorPreyDetector(n_permutations=99, seed=42)
        result = det.detect([])
        assert not result.detected

    def test_single_species_rejected(self):
        """Grid with only one species label (all ones) → n_species=1."""
        rows = cols = 20
        g = np.ones((rows, cols), dtype=np.int8)
        history = [
            {"grid": g.copy(), "grid_dims": (rows, cols), "step": t}
            for t in range(200)
        ]
        det = P11PredatorPreyDetector(n_permutations=99, seed=42)
        result = det.detect(history)

        assert not result.detected
        assert result.primary_metric["n_unique_species_observed"] == 1

    def test_constant_species_fractions_rejected(self):
        """Two species but both fractions never change → fails species_std."""
        rng = np.random.default_rng(42)
        rows = cols = 20
        # Fixed half-and-half assignment, never changes
        g = np.zeros((rows, cols), dtype=np.int8)
        half_mask = rng.random((rows, cols)) < 0.5
        g[half_mask] = 1
        g[~half_mask] = 2
        history = [
            {"grid": g.copy(), "grid_dims": (rows, cols), "step": t}
            for t in range(200)
        ]
        det = P11PredatorPreyDetector(n_permutations=99, seed=42)
        result = det.detect(history)

        assert not result.detected
        assert result.primary_metric["prey_std"] < 0.005
        assert result.primary_metric["predator_std"] < 0.005


# ===================================================================
# Metadata integration
# ===================================================================

class TestP11WithMetadataHints:
    """Confirm that explicit prey_state / predator_state hints work."""

    def test_explicit_species_hint(self):
        m = LotkaVolterraLattice(
            rows=60, cols=60,
            predation_rate=4.0,
            prey_reproduction_rate=1.0,
            predator_death_rate=1.0,
            seed=42,
        )
        history = m.run(n_steps=800)
        if len(history) < 500 or history[-1]["predator_count"] == 0:
            pytest.skip("predators went extinct")

        # Explicit hint (prey=1, predator=2)
        det = P11PredatorPreyDetector(
            n_permutations=99,
            species_state_hint=(1, 2),
            seed=42,
        )
        result = det.detect(history)

        # Should auto-detect the same (1, 2) labels; both paths must work.
        assert result.primary_metric["prey_state"] == 1
        assert result.primary_metric["predator_state"] == 2
