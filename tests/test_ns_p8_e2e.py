"""End-to-end tests for Nagel-Schreckenberg + P8 traffic-jamming detector.

Coverage:
- Canonical positive (rho=0.15, p=0.3) → DEFINITIVE
- Confirmation-tier regime (rho=0.12, p=0.3) → CONFIRMATION
- Free-flow negative (rho=0.05, p=0.3) → screening_rejection
- Deterministic NS (p=0) at moderate density → screening_rejection
- Pigeonhole density saturation (rho=0.80, p=0) → SCREENING only
  (stopped_fraction is high but jam_lifetime_p95 ≤ 4 rejects at
  confirmation — the designed-for discriminator)
- Fundamental-diagram replication (flow peak ~0.46 near rho=0.10-0.12)
- Broad negative sweep: every other model rejects gracefully
- (SLOW) Replication-quality canonical run at L=1000, 2000 measurement

Timing budget:
- Fast-half target: each test < 15s, total ~60s.
- Slow-half: single replication-quality test, ~3 min.
"""

from __future__ import annotations

import numpy as np
import pytest

from epc.detector_result import DetectionTier
from epc.detectors.p8_traffic_jamming import P8TrafficJammingDetector
from epc.models.nagel_schreckenberg import NagelSchreckenberg


# -----------------------------------------------------------------------------
# Fast-half canonical regimes
# -----------------------------------------------------------------------------


def _run_p8(rho: float, p_slow: float, L: int, n_steps: int, burn_in: int,
            n_permutations: int = 199, seed: int = 42,
            init_mode: str = "uniform") -> tuple:
    """Helper: build model, run, apply P8, return (result, model_metadata)."""
    m = NagelSchreckenberg(
        L=L, density=rho, p_slow=p_slow, v_max=5,
        init_mode=init_mode, seed=seed,
    )
    hist = m.run(n_steps)
    det = P8TrafficJammingDetector(
        n_permutations=n_permutations, burn_in=burn_in, seed=42,
    )
    return det.detect(hist, model_metadata=m.get_metadata()), m.get_metadata()


class TestCanonicalRegimes:
    """P8 tier assignments across canonical NS regimes.

    Fast-half sizing: L=500 with 400 burn-in + 1000 measurement. At L=500
    the statistics are noisier than the L=1000 Phase 1 characterization
    but the DEFINITIVE/CONFIRMATION/SCREENING boundaries are preserved.
    """

    def test_canonical_positive_is_definitive(self):
        """rho=0.15, p=0.3 at L=500 → DEFINITIVE.

        The canonical NS jamming regime — primary signature of the Nagel-
        Schreckenberg 1992 paper. Stopped-fraction ≈ 0.18, jam_lifetime
        heavy-tailed (p95 > 5, max > 20), null p ≤ 0.01.
        """
        result, _ = _run_p8(rho=0.15, p_slow=0.3, L=500, n_steps=1400,
                            burn_in=400, n_permutations=199)
        assert result.detected, "canonical jam must be detected"
        assert result.tier == DetectionTier.DEFINITIVE, (
            f"canonical jam should be DEFINITIVE, got {result.tier.name} "
            f"(stopped={result.primary_metric['stopped_fraction']:.3f}, "
            f"lt_p95={result.secondary_metrics.get('jam_lifetime_p95', 0):.1f}, "
            f"lt_max={result.secondary_metrics.get('jam_lifetime_max', 0)}, "
            f"null_p={result.null_p_value:.4f})"
        )
        # Specific quantitative anchors — pin these so a threshold or
        # algorithm regression forces explicit review.
        sf = result.primary_metric["stopped_fraction"]
        assert 0.15 < sf < 0.25, f"stopped_fraction={sf:.3f} outside canonical range [0.15, 0.25]"
        lt_p95 = result.secondary_metrics["jam_lifetime_p95"]
        assert lt_p95 > 5.0, f"jam_lifetime_p95={lt_p95} must exceed confirmation threshold 5"
        lt_max = result.secondary_metrics["jam_lifetime_max"]
        assert lt_max > 20, f"jam_lifetime_max={lt_max} must exceed definitive threshold 20"
        assert result.null_p_value < 0.01, f"null_p={result.null_p_value:.4f} not < 0.01"
        print(f"  ✓ canonical (rho=0.15, p=0.3): {result.tier.name}, "
              f"stopped={sf:.3f}, lt_p95={lt_p95:.1f}, lt_max={lt_max}, "
              f"null_p={result.null_p_value:.4f}")

    def test_deep_jam_is_definitive(self):
        """rho=0.30, p=0.3 → DEFINITIVE. Stopped fraction ≈ 0.43."""
        result, _ = _run_p8(rho=0.30, p_slow=0.3, L=500, n_steps=1400,
                            burn_in=400, n_permutations=199)
        assert result.detected
        assert result.tier == DetectionTier.DEFINITIVE, (
            f"deep jam should be DEFINITIVE, got {result.tier.name}"
        )
        sf = result.primary_metric["stopped_fraction"]
        assert 0.38 < sf < 0.50, f"stopped_fraction={sf:.3f} outside expected [0.38, 0.50]"
        print(f"  ✓ deep jam (rho=0.30): {result.tier.name}, stopped={sf:.3f}")

    def test_near_transition_is_confirmation(self):
        """rho=0.12, p=0.3 → CONFIRMATION at L=1000.

        Stopped_fraction is below the DEFINITIVE threshold (0.15) but
        above the SCREENING threshold (0.05), and jam_lifetime_p95 > 5
        still passes the confirmation gate. This is the Sprint 15
        'canonical CONFIRMATION' example (analogous to GS spots at
        Sprint 14.6).

        NB: This test REQUIRES L=1000. At L=500 the onset regime is
        contaminated by finite-size effects (only ~60 cars at rho=0.12)
        and stopped_fraction fluctuates below the 0.05 screening floor
        for some seeds. The CONFIRMATION window at L=500 is effectively
        empty: the transition is so sharp that rho jumps from SCREENING
        (stopped ≤ 0.05) to DEFINITIVE (stopped ≥ 0.15) without a
        CONFIRMATION regime in between.
        """
        result, _ = _run_p8(rho=0.12, p_slow=0.3, L=1000, n_steps=1400,
                            burn_in=400, n_permutations=199)
        assert result.detected, "near-transition regime must be detected"
        assert result.tier == DetectionTier.CONFIRMATION, (
            f"rho=0.12 at L=1000 should be CONFIRMATION (not DEFINITIVE), "
            f"got {result.tier.name}. "
            f"stopped={result.primary_metric['stopped_fraction']:.3f}, "
            f"lt_p95={result.secondary_metrics.get('jam_lifetime_p95', 0):.1f}, "
            f"lt_max={result.secondary_metrics.get('jam_lifetime_max', 0)}"
        )
        sf = result.primary_metric["stopped_fraction"]
        # Must be BELOW definitive threshold (0.15) and ABOVE screening (0.05)
        assert 0.05 < sf < 0.15, (
            f"rho=0.12 stopped={sf:.3f} outside CONFIRMATION window (0.05, 0.15)"
        )
        lt_p95 = result.secondary_metrics["jam_lifetime_p95"]
        assert lt_p95 > 5.0, f"lt_p95={lt_p95} below confirmation gate"
        print(f"  ✓ near-transition (rho=0.12, L=1000): {result.tier.name}, "
              f"stopped={sf:.3f}, lt_p95={lt_p95:.1f}")

    def test_free_flow_rejected_at_screening(self):
        """rho=0.05, p=0.3 → SCREENING (no detection).

        In deep free-flow, no car ever stops after burn-in, so
        stopped_fraction ≈ 0 fails the screening floor. This is the
        primary P8 negative-at-screening case.
        """
        result, _ = _run_p8(rho=0.05, p_slow=0.3, L=500, n_steps=1400,
                            burn_in=400, n_permutations=99)
        assert not result.detected, "free flow must not be detected"
        assert result.tier == DetectionTier.SCREENING
        sf = result.primary_metric["stopped_fraction"]
        assert sf < 0.01, f"free-flow stopped_fraction={sf:.4f} should be ~0"
        reason = result.primary_metric.get("screening_rejection_reason", "?")
        assert reason == "below_stopped_floor", (
            f"expected below_stopped_floor, got {reason!r}"
        )
        print(f"  ✓ free flow (rho=0.05): SCREENING, stopped={sf:.4f}, "
              f"reason={reason}")

    def test_deterministic_at_low_density_rejected(self):
        """p=0 at rho=0.15 → SCREENING.

        Deterministic NS at sub-saturation density has no randomization,
        so no spontaneous jams form — cars reach a stable moving pattern.
        Stopped fraction ≈ 0.
        """
        result, _ = _run_p8(rho=0.15, p_slow=0.0, L=500, n_steps=1400,
                            burn_in=400, n_permutations=99)
        assert not result.detected
        assert result.tier == DetectionTier.SCREENING
        sf = result.primary_metric["stopped_fraction"]
        assert sf < 0.01
        print(f"  ✓ deterministic (rho=0.15, p=0): SCREENING, stopped={sf:.4f}")

    def test_density_saturation_rejected_at_confirmation(self):
        """rho=0.80, p=0 → SCREENING (primary passes, confirmation fails).

        This is the DESIGNED-FOR discriminator of the P8 detector:
        pigeonhole density-saturation at rho > 1/(v_max+1) produces a
        high stopped_fraction (~0.75) because cars physically cannot all
        move, BUT jam_lifetime_p95 ≈ 4 — short uniform stops, no heavy
        tail — so confirmation rejects. True NS jamming (Sprint 15
        canonical) produces lt_p95 = 13. The ~3x separation is the
        discriminator.

        This mirrors the Sprint 13 RPS-vs-GrayScott false-positive trap
        where raw integer-grid FFT gave matching peak-to-mean and the
        discrimination had to happen at substrate prerequisite level.
        Here, discrimination happens at the confirmation secondary
        prerequisite level.
        """
        result, _ = _run_p8(rho=0.80, p_slow=0.0, L=500, n_steps=1400,
                            burn_in=400, n_permutations=99)
        # Screening should pass (stopped_fraction is high):
        assert result.detected, (
            "density-saturation regime passes screening (stopped_fraction high)"
        )
        # But confirmation must fail:
        assert result.tier == DetectionTier.SCREENING, (
            f"density saturation should STOP AT SCREENING, got {result.tier.name}. "
            f"This is the key P8 discriminator — pigeonhole high stopped_fraction "
            f"must not pass the jam_lifetime_p95 gate."
        )
        sf = result.primary_metric["stopped_fraction"]
        assert sf > 0.5, f"density saturation stopped={sf:.3f} should be > 0.5"
        lt_p95 = result.secondary_metrics.get("jam_lifetime_p95", 0.0)
        assert lt_p95 <= 5.0, (
            f"density saturation lt_p95={lt_p95} must not exceed 5 "
            f"(the confirmation gate). If this fails, the P8 discriminator "
            f"is no longer separating pigeonhole from true jamming."
        )
        print(f"  ✓ density saturation (rho=0.80, p=0): correctly stops at "
              f"SCREENING, stopped={sf:.3f}, lt_p95={lt_p95:.1f} (≤5 "
              f"confirmation gate)")


class TestFundamentalDiagramReplication:
    """Verify the NS model replicates the published fundamental diagram.

    Published anchors (Nagel-Schreckenberg 1992; Bette et al. 2017;
    Wikipedia):
      - Free-flow mean velocity: <v> = v_max - p in the dilute limit.
      - Peak flow at rho ≈ 0.10-0.15 with flow ~ 0.46-0.47 at p=0.3, v_max=5.
      - Stopped fraction ~ 0 below rho=0.10, > 0.15 above rho=0.15.
    """

    def test_dilute_limit_mean_velocity(self):
        """At rho=0.02, p=0.3, v_max=5: mean velocity → v_max - p = 4.7 (analytic)."""
        m = NagelSchreckenberg(L=500, density=0.02, p_slow=0.3, v_max=5, seed=42)
        hist = m.run(1500)
        post_burn = hist[400:]
        mean_v = np.mean([h["mean_velocity"] for h in post_burn])
        assert 4.60 < mean_v < 4.75, (
            f"dilute-limit <v>={mean_v:.3f} should be ~4.7 (v_max - p)"
        )
        print(f"  ✓ dilute limit <v>={mean_v:.3f} (expected 4.70 = v_max - p)")

    def test_flow_peak_near_expected_density(self):
        """Sweep rho at L=1000 and check peak flow is at rho ∈ [0.08, 0.15]
        with value ~ 0.46 (published Nagel-Schreckenberg 1992 value).

        NB: uses L=1000 because the fundamental-diagram peak magnitude
        has significant finite-size correction — at L=500 the peak is
        ~0.50 rather than the asymptotic ~0.46. Total runtime ~1.5s.
        """
        densities = [0.05, 0.08, 0.10, 0.12, 0.15, 0.20, 0.30]
        flows = []
        for rho in densities:
            m = NagelSchreckenberg(L=1000, density=rho, p_slow=0.3, v_max=5, seed=42)
            hist = m.run(1200)
            post_burn = hist[400:]
            flow = np.mean([h["flow"] for h in post_burn])
            flows.append(flow)
        peak_idx = int(np.argmax(flows))
        peak_rho = densities[peak_idx]
        peak_flow = flows[peak_idx]
        # Peak must be in [0.08, 0.15]
        assert 0.08 <= peak_rho <= 0.15, (
            f"peak flow at rho={peak_rho}, should be in [0.08, 0.15]. "
            f"Full sweep: {dict(zip(densities, [f'{f:.3f}' for f in flows]))}"
        )
        # Peak flow magnitude ≈ 0.46 at L=1000 (published NS1992 / Wikipedia)
        assert 0.42 < peak_flow < 0.50, (
            f"peak flow={peak_flow:.3f} should be ~0.46 (published ~0.46 at L=1000)"
        )
        print(f"  ✓ flow peak: rho={peak_rho}, flow={peak_flow:.3f} "
              f"(expected rho ∈ [0.08, 0.15], flow ~ 0.46)")


class TestSubstrateRobustness:
    """P8 must reject every substrate that is not NS-like."""

    def test_p8_rejects_without_velocities(self):
        """Synthetic state history without 'velocities' → substrate_mismatch."""
        fake_hist = [
            {"grid": np.zeros((20, 20), dtype=np.int32), "step": t}
            for t in range(200)
        ]
        det = P8TrafficJammingDetector(n_permutations=19, burn_in=0, seed=42)
        result = det.detect(fake_hist, model_metadata=None)
        assert not result.detected
        assert result.tier == DetectionTier.SCREENING
        assert result.primary_metric["screening_rejection_reason"] == "substrate_mismatch"

    def test_p8_rejects_2d_velocities(self):
        """Vicsek/D'Orsogna-style 2D float velocities → substrate_mismatch."""
        fake_hist = [
            {"velocities": np.random.random((50, 2)), "step": t}
            for t in range(200)
        ]
        det = P8TrafficJammingDetector(n_permutations=19, burn_in=0, seed=42)
        result = det.detect(fake_hist, model_metadata=None)
        assert not result.detected
        assert result.primary_metric["screening_rejection_reason"] == "substrate_mismatch"

    def test_p8_rejects_float_velocities(self):
        """1D float velocities → non_integer_velocities rejection reason."""
        fake_hist = [
            {"velocities": np.random.random(50).astype(np.float32), "step": t}
            for t in range(200)
        ]
        det = P8TrafficJammingDetector(n_permutations=19, burn_in=0, seed=42)
        result = det.detect(fake_hist, model_metadata=None)
        assert not result.detected
        assert result.primary_metric["screening_rejection_reason"] == "non_integer_velocities"

    def test_p8_rejects_out_of_range_velocities(self):
        """1D integer velocities with absurd range → velocity_range_out_of_bounds.

        Guards against mis-labeled observables where some OTHER integer-
        valued 1D array accidentally gets called 'velocities'.
        """
        rng = np.random.default_rng(42)
        fake_hist = [
            {"velocities": rng.integers(0, 10000, 50, dtype=np.int64), "step": t}
            for t in range(200)
        ]
        det = P8TrafficJammingDetector(n_permutations=19, burn_in=0, seed=42)
        result = det.detect(fake_hist, model_metadata=None)
        assert not result.detected
        assert result.primary_metric["screening_rejection_reason"] == "velocity_range_out_of_bounds"

    def test_p8_rejects_too_few_cars(self):
        """n_cars < 20 → too_few_cars."""
        rng = np.random.default_rng(42)
        fake_hist = [
            {"velocities": rng.integers(0, 6, 5, dtype=np.int8), "step": t}
            for t in range(200)
        ]
        det = P8TrafficJammingDetector(n_permutations=19, burn_in=0, seed=42)
        result = det.detect(fake_hist, model_metadata=None)
        assert not result.detected
        assert result.primary_metric["screening_rejection_reason"] == "too_few_cars"

    def test_p8_rejects_short_run(self):
        """Post-burn-in run < 100 steps → run_too_short."""
        rng = np.random.default_rng(42)
        fake_hist = [
            {"velocities": rng.integers(0, 6, 50, dtype=np.int8), "step": t}
            for t in range(50)
        ]
        det = P8TrafficJammingDetector(n_permutations=19, burn_in=0, seed=42)
        result = det.detect(fake_hist, model_metadata=None)
        assert not result.detected
        assert result.primary_metric["screening_rejection_reason"] == "run_too_short"


class TestModelInvariants:
    """Verify the NS model itself preserves its invariants."""

    def test_car_count_conserved(self):
        """N_cars must be exactly conserved across all steps."""
        m = NagelSchreckenberg(L=200, density=0.3, p_slow=0.3, v_max=5, seed=42)
        hist = m.run(300)
        n_initial = hist[0]["n_cars"]
        for t, snap in enumerate(hist):
            assert snap["n_cars"] == n_initial, (
                f"car count broken at step {t}: {snap['n_cars']} != {n_initial}"
            )
            assert len(snap["velocities"]) == n_initial
            assert len(snap["positions"]) == n_initial

    def test_no_collisions(self):
        """All positions must be distinct at every step (no two cars
        in the same cell)."""
        m = NagelSchreckenberg(L=200, density=0.3, p_slow=0.3, v_max=5, seed=42)
        hist = m.run(300)
        for t, snap in enumerate(hist):
            positions = snap["positions"]
            assert len(np.unique(positions)) == len(positions), (
                f"duplicate positions at step {t}: {positions}"
            )

    def test_velocity_range(self):
        """Velocities must always lie in [0, v_max]."""
        m = NagelSchreckenberg(L=200, density=0.3, p_slow=0.3, v_max=5, seed=42)
        hist = m.run(300)
        for t, snap in enumerate(hist):
            v = snap["velocities"]
            assert v.min() >= 0, f"negative velocity at step {t}"
            assert v.max() <= 5, f"velocity > v_max at step {t}: {v.max()}"

    def test_gap_consistency(self):
        """Gaps must equal (next_position - position) mod L - 1."""
        m = NagelSchreckenberg(L=200, density=0.3, p_slow=0.3, v_max=5, seed=42)
        hist = m.run(200)
        for t, snap in enumerate(hist):
            positions = snap["positions"]
            gaps = snap["gaps"]
            if len(positions) <= 1:
                continue
            expected = (np.roll(positions, -1) - positions) % 200 - 1
            assert np.array_equal(gaps, expected), (
                f"gap inconsistency at step {t}"
            )

    def test_reproducibility_via_seed(self):
        """Same seed → identical state history."""
        m1 = NagelSchreckenberg(L=100, density=0.2, p_slow=0.3, v_max=5, seed=1234)
        m2 = NagelSchreckenberg(L=100, density=0.2, p_slow=0.3, v_max=5, seed=1234)
        h1 = m1.run(200)
        h2 = m2.run(200)
        for t, (s1, s2) in enumerate(zip(h1, h2)):
            assert np.array_equal(s1["velocities"], s2["velocities"]), (
                f"divergence at step {t}"
            )


# -----------------------------------------------------------------------------
# Slow-half: replication-quality canonical run
# -----------------------------------------------------------------------------


@pytest.mark.slow
class TestReplicationQuality:
    """Replication-quality canonical run matching the Phase 1 characterization.

    L=1000, v_max=5, p_slow=0.3, rho=0.15; 1000 burn-in + 2000 measurement
    (~7 seconds per seed with 199 null permutations; total < 3 minutes
    with 3 seeds). Pinned quantitative anchors come directly from Phase 1
    characterization at Sprint 15 HEAD.
    """

    @pytest.mark.parametrize("seed", [42, 123, 2024])
    def test_canonical_l1000_is_definitive(self, seed):
        """L=1000 canonical: stopped ≈ 0.18, lt_p95 ≈ 13, tier DEFINITIVE."""
        result, _ = _run_p8(rho=0.15, p_slow=0.3, L=1000, n_steps=3000,
                            burn_in=1000, n_permutations=199, seed=seed)
        assert result.detected
        assert result.tier == DetectionTier.DEFINITIVE
        sf = result.primary_metric["stopped_fraction"]
        # Tight Phase 1 anchors: seed-to-seed sigma ~ 0.003 at L=1000
        assert 0.16 < sf < 0.21, (
            f"seed={seed}: stopped={sf:.3f} outside Phase 1 range [0.16, 0.21]"
        )
        lt_p95 = result.secondary_metrics["jam_lifetime_p95"]
        assert lt_p95 >= 10, f"seed={seed}: lt_p95={lt_p95} below Phase 1 anchor ~13"
        print(f"  ✓ L=1000 seed={seed}: {result.tier.name}, stopped={sf:.3f}, "
              f"lt_p95={lt_p95:.1f}")


# -----------------------------------------------------------------------------
# Slow-half: finite-size scaling (Sprint 19 — closes carry-forward #15)
# -----------------------------------------------------------------------------


@pytest.mark.slow
class TestFiniteSizeScaling:
    """Finite-size scaling test pinning the minimum L for reliable
    CONFIRMATION/DEFINITIVE tier across seeds.

    Sprint 15 REPLICATION_NOTES documented that at L = 500 near the
    critical density ρ = 0.12, stopped_fraction fluctuates near the
    screening threshold (0.05), so some seeds land at SCREENING rather
    than CONFIRMATION. The canonical CONFIRMATION demonstration requires
    L = 1000. This test pins that finding by running the canonical
    ρ = 0.15, p_slow = 0.3 regime across L ∈ {250, 500, 1000} with three
    seeds at each L, and asserting the expected finite-size behavior:

    - At L ≥ 500, all seeds reach DEFINITIVE (the canonical regime is
      sufficiently above the transition that finite-size variance is
      bounded).
    - At L = 250, we relax to at least CONFIRMATION across all seeds;
      SCREENING-only outcomes at L = 250 would indicate an emerging
      finite-size sensitivity.
    - stopped_fraction spread across seeds contracts monotonically with
      L (σ(sf) decreasing as L increases), confirming the standard
      1/√L finite-size scaling for this lattice_1d order parameter.

    Runtime: ~4–8 minutes total for 9 runs (3 L-values × 3 seeds) at
    n_permutations = 199 (the canonical P8 setting — n_perm = 99 hits
    the null-p floor of 0.01 exactly and fails the strict `< 0.01`
    confirmation gate, per the Sprint 8b null-p-floor finding).
    """

    @pytest.mark.parametrize("L", [250, 500, 1000])
    @pytest.mark.parametrize("seed", [42, 123, 2024])
    def test_definitive_across_L_and_seed(self, L, seed):
        """Canonical ρ=0.15 regime reaches at least CONFIRMATION at every
        tested (L, seed); DEFINITIVE at L ≥ 500."""
        # Scale burn_in + measurement proportionally so total cars × steps
        # gives comparable per-trajectory statistics across L.
        n_steps = 1500 if L <= 500 else 2000
        burn_in = 500 if L <= 500 else 800
        result, _ = _run_p8(rho=0.15, p_slow=0.3, L=L, n_steps=n_steps,
                            burn_in=burn_in, n_permutations=199, seed=seed)
        assert result.detected, (
            f"L={L} seed={seed}: expected at least SCREENING, got {result.tier.name}"
        )
        # Every tested (L, seed) should clear CONFIRMATION at the canonical
        # ρ=0.15 regime; L=250 is the lower bound we're pinning as safe.
        assert result.tier.value >= DetectionTier.CONFIRMATION.value, (
            f"L={L} seed={seed}: finite-size regression — "
            f"canonical ρ=0.15 should reach CONFIRMATION, got {result.tier.name}"
        )
        if L >= 500:
            assert result.tier == DetectionTier.DEFINITIVE, (
                f"L={L} seed={seed}: L≥500 should give DEFINITIVE, got "
                f"{result.tier.name}; stopped_fraction="
                f"{result.primary_metric['stopped_fraction']:.3f}"
            )

    def test_stopped_fraction_spread_contracts_with_L(self):
        """σ(stopped_fraction) across seeds should decrease with L
        (standard 1/√L finite-size scaling for ring order parameters).

        This is the structural test that motivated the carry-forward:
        Sprint 15 observed the near-transition regime had growing σ at
        small L; we confirm the canonical regime shows the expected
        shrinking σ as L grows.
        """
        sigmas = {}
        for L in [250, 500, 1000]:
            n_steps = 1500 if L <= 500 else 2000
            burn_in = 500 if L <= 500 else 800
            sf_values = []
            for seed in [42, 123, 2024]:
                result, _ = _run_p8(rho=0.15, p_slow=0.3, L=L, n_steps=n_steps,
                                    burn_in=burn_in, n_permutations=199,
                                    seed=seed)
                sf_values.append(result.primary_metric["stopped_fraction"])
            sigmas[L] = float(np.std(sf_values, ddof=1))
            print(f"  L={L}: stopped_fraction = {sf_values}, σ={sigmas[L]:.4f}")

        # Monotone shrinkage is the strong claim, but 3 seeds give noisy
        # σ estimates; assert the weaker but robust claim that σ at L=1000
        # is no larger than σ at L=250.
        assert sigmas[1000] <= sigmas[250] + 0.01, (
            f"σ(stopped_fraction) at L=1000 ({sigmas[1000]:.4f}) exceeds "
            f"σ at L=250 ({sigmas[250]:.4f}) + 0.01 — finite-size "
            f"scaling regression"
        )
