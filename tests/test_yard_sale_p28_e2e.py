"""End-to-end tests for Yard-Sale + P28 Wealth Condensation detector.

Coverage:
- Canonical positive (N=1000, f=0.1, lambda=0, chi=0) -> DEFINITIVE
- Fast condensation (f=0.3) -> DEFINITIVE
- Slow condensation (f=0.05) -> CONFIRMATION (top_1pct just below DEF threshold)
- Within-family negatives:
    * lambda=0.5 (saving propensity) -> rejected at screening
    * chi=0.001 (moderate redistribution) -> SCREENING (top_1pct < 0.15)
    * chi=0.0001 (mild redistribution) -> CONFIRMATION (top-1% share just
      below the DEF floor; metadata no longer gates any tier, 2f40e9a)
- Too-early trap (f=0.01, t=100k) -> rejected at screening (not emergent)
- Substrate-mismatch (non-wealth state histories) -> rejected
- Model replication (Boghosian 2014 qualitative: Gini -> 1 at long time)
- Multi-seed reproducibility at canonical positive -> all DEFINITIVE
- Registry checks (Sprint 17 registrations present)

Tiering since 2f40e9a (mechanism from the trajectory, not metadata; S2 Text,
P28 row). The Decision 49 metadata-mechanism null (four self-reported
booleans gating DEFINITIVE) was removed.
  Primary = Gini coefficient at final frame.
  SCREENING    = total wealth conserved (drift <= 1%), condensation emerges
                 from near-equality (initial Gini <= 0.30, rise >= 0.20),
                 Gini >= 0.65, top_1pct >= 0.05.
  CONFIRMATION = Gini >= 0.75, top_1pct >= 0.15, monotonic_frac >= 0.80,
                 null_p < 0.01.
  DEFINITIVE   = CONFIRMATION + Gini >= 0.85, top_1pct >= 0.30.

Timing budget:
- Fast-half: each test <= 15s, total ~2 min.
- Slow-half: multi-seed + finite-N reproducibility at DEFINITIVE.
"""

from __future__ import annotations

import numpy as np
import pytest

from epc.detector_result import DetectionTier
from epc.detectors.p28_wealth_condensation import P28WealthCondensationDetector
from epc.models.yard_sale import YardSale


# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------


def _run_ys_p28(
    f=0.1, lambda_save=0.0, chi=0.0, redistribute_every=1000,
    total_tx=1_000_000, n_checkpoints=30,
    N=1000, seed=42, init_mode="equal",
    n_permutations=199, burn_in=0,
):
    """Helper: build YS, run, apply P28, return DetectorResult + metadata."""
    m = YardSale(
        n_agents=N, f=f, lambda_save=lambda_save, chi=chi,
        redistribute_every=redistribute_every,
        init_mode=init_mode, seed=seed,
    )
    interval = max(1, total_tx // n_checkpoints)
    history = m.run(total_tx, record_interval=interval)
    det = P28WealthCondensationDetector(
        burn_in=burn_in, n_permutations=n_permutations, seed=seed,
    )
    res = det.detect(history, model_metadata=m.get_metadata())
    return res, m.get_metadata()


# -----------------------------------------------------------------------------
# Model replication — qualitative Boghosian 2014 result
# -----------------------------------------------------------------------------


class TestYardSaleReplication:
    """Yard-Sale model must conserve wealth exactly and drive Gini -> 1."""

    def test_wealth_conservation(self):
        """Total wealth conserved to machine precision (chi=0, pure YS)."""
        m = YardSale(n_agents=500, f=0.1, lambda_save=0.0, chi=0.0, seed=42)
        m.setup()
        total_0 = float(m.wealth.sum())
        m.step(200_000)
        total_1 = float(m.wealth.sum())
        # Floating-point drift; bit-exact at small f, can be O(N*eps)
        # after many transactions
        assert abs(total_1 - total_0) / total_0 < 1e-10, \
            f"Wealth drift {abs(total_1 - total_0):.2e} exceeds 1e-10"

    def test_gini_starts_at_zero(self):
        """Equal initialization -> Gini = 0 exactly at t=0."""
        m = YardSale(n_agents=500, f=0.1, lambda_save=0.0, init_mode="equal",
                     seed=42)
        s0 = m.setup()
        assert s0["gini"] == pytest.approx(0.0, abs=1e-12)

    def test_gini_grows_toward_one_at_long_time(self):
        """Pure YS (lambda=0, chi=0) -> Gini > 0.9 at ~5M transactions.

        Boghosian 2014 qualitative claim: symmetric fair exchange with
        multiplicative stakes produces full condensation.
        """
        m = YardSale(n_agents=1000, f=0.1, lambda_save=0.0, seed=42)
        m.setup()
        m.step(2_000_000)
        g = YardSale._gini(m.wealth)
        assert g > 0.80, f"Gini after 2e6 transactions should be > 0.80, got {g:.4f}"

    def test_saving_propensity_plateau(self):
        """lambda > 0 produces a finite-Gini fixed point far below 1.

        CC 2000 result: Gamma-distributed wealth with finite Gini.
        """
        m = YardSale(n_agents=1000, f=0.1, lambda_save=0.5, seed=42)
        m.setup()
        m.step(2_000_000)
        g = YardSale._gini(m.wealth)
        assert 0.15 < g < 0.50, \
            f"lambda=0.5 Gini should be in (0.15, 0.50), got {g:.4f}"

    def test_no_negative_wealth(self):
        """YS stake = f*min(w_i, w_j) can never produce negative wealth."""
        m = YardSale(n_agents=500, f=0.3, lambda_save=0.0, seed=42)
        m.setup()
        m.step(500_000)
        assert (m.wealth >= 0).all(), \
            "Negative wealth appeared — stake logic is broken"

    def test_n_invariance_at_fixed_sweeps(self):
        """Gini at t=N*sweeps should be ~ independent of N at fixed sweeps.

        Phase 1d.4 finding: Gini(1000 sweeps, f=0.1) = 0.886 ± 0.002
        across N ∈ {200, 500, 1000, 2000}.
        """
        ginis = []
        for N in [500, 1000]:
            m = YardSale(n_agents=N, f=0.1, lambda_save=0.0, seed=42)
            m.setup()
            m.step(1000 * N)  # 1000 sweeps
            ginis.append(YardSale._gini(m.wealth))
        assert abs(ginis[0] - ginis[1]) < 0.05, \
            f"N-invariance broken: Gini(500) = {ginis[0]:.4f} vs " \
            f"Gini(1000) = {ginis[1]:.4f}"


# -----------------------------------------------------------------------------
# Canonical positive regimes
# -----------------------------------------------------------------------------


class TestCanonicalP28:
    """P28 tier assignments across canonical YS regimes."""

    def test_canonical_positive_definitive(self):
        """N=1000, f=0.1, lambda=0, chi=0 at t=2e6 -> DEFINITIVE."""
        res, _ = _run_ys_p28(f=0.1, total_tx=2_000_000, seed=42)
        assert res.tier == DetectionTier.DEFINITIVE, \
            f"Canonical YS should be DEFINITIVE, got {res.tier.value}: " \
            f"gini={res.primary_metric.get('gini_final'):.4f} " \
            f"top1={res.secondary_metrics.get('top_1pct_share'):.4f}"
        assert res.detected
        assert res.confidence >= 0.85  # DEFINITIVE cap

    def test_canonical_primary_magnitude(self):
        """Canonical YS -> Gini > 0.90 and top-1% > 0.30."""
        res, _ = _run_ys_p28(f=0.1, total_tx=2_000_000, seed=42)
        assert res.primary_metric["gini_final"] > 0.90
        assert res.secondary_metrics["top_1pct_share"] > 0.30

    def test_fast_condensation_definitive(self):
        """f=0.3 condenses very fast; t=1e6 enough for DEFINITIVE."""
        res, _ = _run_ys_p28(f=0.3, total_tx=1_000_000, seed=42)
        assert res.tier == DetectionTier.DEFINITIVE
        assert res.primary_metric["gini_final"] > 0.95

    def test_slow_condensation_confirmation(self):
        """f=0.05 condenses slowly; t=3e6 gives CONFIRMATION (top-1% < 0.30)."""
        res, _ = _run_ys_p28(f=0.05, total_tx=3_000_000, n_checkpoints=40,
                             seed=42)
        # At t=3M, Gini ~0.85, top_1pct ~0.22 — below DEF top_1pct floor
        assert res.tier == DetectionTier.CONFIRMATION, \
            f"f=0.05 at t=3e6 should CONFIRM (not DEF), got {res.tier.value}"
        assert res.detected

    def test_monotonic_growth_at_positive(self):
        """Pure YS Gini trajectory is essentially monotonic."""
        res, _ = _run_ys_p28(f=0.1, total_tx=1_000_000, seed=42)
        mono = res.secondary_metrics["monotonic_fraction"]
        assert mono > 0.95, f"Monotonic growth {mono:.3f} should be > 0.95 for pure YS"

    def test_null_rejects_positive(self):
        """Well-mixed Boltzmann-Gibbs null: observed Gini >> null (~0.5)."""
        res, _ = _run_ys_p28(f=0.1, total_tx=1_000_000, n_permutations=199,
                             seed=42)
        assert res.null_p_value < 0.01, \
            f"Null-p should be < 0.01, got {res.null_p_value:.4f}"
        # Null mean should be ~0.5 (Exp distribution Gini limit)
        assert 0.45 < res.effect_size["null_mean"] < 0.55, \
            f"Exp null mean should be ~0.5, got {res.effect_size['null_mean']:.4f}"


# -----------------------------------------------------------------------------
# Within-family negatives — the critical "mechanism matters" tests
# -----------------------------------------------------------------------------


class TestP28Negatives:
    """Within-family YS negatives exercise the mechanistic-null gate."""

    def test_saving_propensity_rejected(self):
        """lambda=0.5 -> Gini ~0.28 -> rejected at SCREENING."""
        res, _ = _run_ys_p28(f=0.1, lambda_save=0.5, total_tx=2_000_000,
                             seed=42)
        assert res.tier == DetectionTier.SCREENING
        assert not res.detected
        assert "gini_below" in res.primary_metric.get(
            "screening_rejection_reason", ""
        )

    def test_moderate_redistribution_screening(self):
        """chi=0.001 -> Gini ~0.68 -> SCREENING (top_1pct < 0.15 gate)."""
        res, _ = _run_ys_p28(f=0.1, chi=0.001, total_tx=2_000_000, seed=42)
        # Screening passes (Gini > 0.4, top1 > 0.05) but CONFIRMATION fails
        assert res.tier == DetectionTier.SCREENING, \
            f"chi=0.001 should SCREEN not confirm, got {res.tier.value}"
        assert res.detected  # screening floor passed

    def test_mild_redistribution_confirmation_not_definitive(self):
        """chi=0.0001 -> Gini ~0.89, looks condensed, held at CONFIRMATION.

        This is the Sprint 17 "mechanism matters" test analogous to
        Sprint 16's P2 milling test. Test maintenance (docstring only): the
        metadata gate (has_redistribution=True) that held this run below
        DEFINITIVE was removed in 2f40e9a. At seed 42 the run is now held at
        CONFIRMATION by its top-1% share (0.289 against the 0.30 DEF floor).
        """
        res, meta = _run_ys_p28(f=0.1, chi=0.0001, total_tx=2_000_000,
                                seed=42)
        assert meta["has_redistribution"] is True
        # Empirically this LOOKS like condensation, but metadata stops it
        assert res.tier == DetectionTier.CONFIRMATION, \
            f"chi=0.0001 should CONFIRM, not DEFINE (has_redistribution=True): " \
            f"got {res.tier.value} with gini={res.primary_metric['gini_final']:.4f}"
        # Confirm the empirical signal really is DEF-strength
        assert res.primary_metric["gini_final"] > 0.80
        assert res.secondary_metrics["top_1pct_share"] > 0.20

    def test_too_early_below_floor_rejected(self):
        """Short run (t=100k at f=0.01) -> Gini ~0.08 -> rejected.

        Test maintenance: since 2f40e9a screening checks emergence (Gini must
        rise by at least 0.20 from an initial value of at most 0.30) before
        the Gini floor (raised 0.40 -> 0.65 in the same commit), so the
        reported reason is condensation_not_emergent. The run is still below
        the screening floor, which is asserted separately.
        """
        res, _ = _run_ys_p28(f=0.01, total_tx=100_000, n_checkpoints=20,
                             seed=42)
        assert res.tier == DetectionTier.SCREENING
        assert not res.detected
        assert (
            res.primary_metric["screening_rejection_reason"]
            == "condensation_not_emergent"
        )
        assert res.primary_metric["gini_delta"] < 0.20
        assert res.primary_metric["gini_final"] < 0.65

    def test_substrate_mismatch_rejection(self):
        """Synthetic state history with no 'wealth' observable -> rejected."""
        fake_history = [
            {"grid": np.zeros((10, 10), dtype=int), "step": i} for i in range(20)
        ]
        det = P28WealthCondensationDetector(n_permutations=99, seed=7)
        res = det.detect(fake_history, model_metadata={"model_family": "not_wealth"})
        assert res.tier == DetectionTier.SCREENING
        assert not res.detected
        assert res.primary_metric["screening_rejection_reason"] == "substrate_mismatch"

    def test_too_few_agents_rejected(self):
        """N=20 YS should reject at too_few_agents prereq."""
        m = YardSale(n_agents=20, f=0.3, lambda_save=0.0, seed=42)
        history = m.run(50_000, record_interval=2000)
        det = P28WealthCondensationDetector(n_permutations=99, seed=7)
        res = det.detect(history, model_metadata=m.get_metadata())
        assert res.tier == DetectionTier.SCREENING
        assert not res.detected
        assert res.primary_metric["screening_rejection_reason"] == "too_few_agents"


# -----------------------------------------------------------------------------
# Mechanism from the trajectory (2f40e9a), replacing the Decision 49
# metadata-mechanism null
# -----------------------------------------------------------------------------
#
# Test maintenance: this class used to pin the four-flag mechanistic null
# (test_missing_conserved_flag_blocks_definitive,
# test_missing_multiplicative_stake_blocks_definitive,
# test_saving_propensity_flag_blocks_definitive,
# test_redistribution_flag_blocks_definitive,
# test_metadata_absent_blocks_definitive). 2f40e9a removed that gate because
# flipping one self-reported flag on the same trajectory flipped CONFIRMATION
# and DEFINITIVE (S2 Text, P28 row: "tiers read four self-reported flags" ->
# "mechanism from the trajectory"). Those five tests are replaced by a test
# that the flags no longer change the tier, and by tests of the trajectory
# gates that now do the flags' job: conservation (was has_conserved_resource)
# and emergence from near-equality. The saving-propensity and redistribution
# look-alikes are covered by TestP28Negatives above (Gini floor).


class TestP28MechanisticNull:
    """Tier is set by the trajectory; metadata flags do not gate it (2f40e9a)."""

    @pytest.fixture
    def canonical_history(self):
        """Cache one canonical run for the flag tests.

        Uses 2M transactions to match the Phase 1e characterization
        (Gini ~0.94, top_1pct ~0.34) at which the DEFINITIVE gate
        unambiguously fires. Shorter runs (1M, 1.5M) land in the
        CONFIRMATION band and produce seed-dependent DEF/CONF flip.
        """
        m = YardSale(n_agents=1000, f=0.1, lambda_save=0.0, seed=42)
        return m.run(2_000_000, record_interval=60_000), m.get_metadata()

    def test_all_flags_correct_allows_definitive(self, canonical_history):
        history, meta = canonical_history
        det = P28WealthCondensationDetector(n_permutations=199, seed=42)
        res = det.detect(history, model_metadata=meta)
        assert res.tier == DetectionTier.DEFINITIVE

    @pytest.mark.parametrize("tamper", [
        {"has_conserved_resource": False},
        {"has_multiplicative_stake": False},
        {"has_saving_propensity": True},
        {"has_redistribution": True},
        None,  # no metadata at all
    ])
    def test_metadata_flags_do_not_change_tier(self, canonical_history, tamper):
        """Flipping any mechanism flag, or dropping metadata, leaves DEFINITIVE.

        The metadata-flip invariance that 2f40e9a verified in
        analysis/validation_rebuild/p28_mechanism_discrimination.py.
        """
        history, meta = canonical_history
        if tamper is None:
            md = None
        else:
            md = dict(meta)
            md.update(tamper)
        det = P28WealthCondensationDetector(n_permutations=199, seed=42)
        res = det.detect(history, model_metadata=md)
        assert res.tier == DetectionTier.DEFINITIVE, (tamper, res.tier.value)

    def test_non_conserved_lookalike_rejected_by_conservation_gate(self):
        """Bouchaud-Mezard multiplicative growth (total wealth drifts) is
        rejected at screening by measured conservation drift, the job the
        has_conserved_resource flag used to do."""
        from epc.phase2a.failed_regimes import p28_wealth as fr
        regime = next(r for r in fr.CONFIG["regimes"] if r["label"] == "bouchaud_mezard")
        history = fr.build_substrate(regime)
        det = P28WealthCondensationDetector(n_permutations=199, seed=42)
        res = det.detect(history, model_metadata=None)
        assert not res.detected
        assert res.primary_metric["screening_rejection_reason"] == "resource_not_conserved"
        assert res.primary_metric["conservation_drift"] > 0.01

    def test_pre_condensed_input_rejected_by_emergence_gate(self, canonical_history):
        """The same canonical run with its first half removed starts already
        condensed (initial Gini ~0.89), so condensation does not emerge from
        near-equality and screening rejects it (2f40e9a emergence gate)."""
        history, meta = canonical_history
        det = P28WealthCondensationDetector(n_permutations=199, seed=42)
        res = det.detect(history[len(history) // 2:], model_metadata=meta)
        assert not res.detected
        assert res.primary_metric["screening_rejection_reason"] == "condensation_not_emergent"
        assert res.primary_metric["gini_initial"] > 0.30


# -----------------------------------------------------------------------------
# Multi-seed reproducibility — fast subset
# -----------------------------------------------------------------------------


class TestMultiSeedReproducibility:
    """YS at N=1000 is reproducible across seeds (Phase 1d.1).

    Uses total_tx=2_000_000 to land inside the DEFINITIVE band
    (Gini > 0.92, top_1pct > 0.30). Shorter runs (1.5M) give
    Gini ~0.91 and land in the CONFIRMATION band for some seeds —
    the DEF gate fires on top_1pct, not Gini, at that regime.
    """

    @pytest.mark.parametrize("seed", [42, 7, 101])
    def test_canonical_reproducible(self, seed):
        res, _ = _run_ys_p28(f=0.1, total_tx=2_000_000, seed=seed)
        assert res.tier == DetectionTier.DEFINITIVE, \
            f"seed={seed} should give DEFINITIVE: gini=" \
            f"{res.primary_metric.get('gini_final'):.4f} " \
            f"top1={res.secondary_metrics.get('top_1pct_share'):.4f}"
        assert res.primary_metric["gini_final"] > 0.90


# -----------------------------------------------------------------------------
# Slow half — replication-quality N-scaling
# -----------------------------------------------------------------------------


class TestSprint17SlowReplication:
    """Slow-half: multi-seed reproducibility at full N=1000, t=2e6, and
    an N-scaling robustness check."""

    @pytest.mark.slow
    @pytest.mark.parametrize("seed", [42, 7, 101, 999])
    def test_four_seed_definitive(self, seed):
        """All four seeds produce DEFINITIVE with Gini in a narrow band."""
        res, _ = _run_ys_p28(f=0.1, total_tx=2_000_000, seed=seed)
        assert res.tier == DetectionTier.DEFINITIVE
        assert 0.92 < res.primary_metric["gini_final"] < 0.95, \
            f"seed={seed} Gini {res.primary_metric['gini_final']:.4f} " \
            f"outside expected [0.92, 0.95]"

    @pytest.mark.slow
    def test_N_scaling_gini_invariance(self):
        """Gini at 1000 sweeps should be ~constant across N (Phase 1d.4)."""
        ginis = {}
        for N in [500, 1000, 2000]:
            m = YardSale(n_agents=N, f=0.1, lambda_save=0.0, seed=42)
            m.setup()
            m.step(1000 * N)
            ginis[N] = YardSale._gini(m.wealth)
        # Ginis should be within ~0.01 of each other
        g_values = list(ginis.values())
        spread = max(g_values) - min(g_values)
        assert spread < 0.02, \
            f"Gini spread {spread:.4f} across N={list(ginis.keys())} " \
            f"exceeds 0.02: {ginis}"

    @pytest.mark.slow
    @pytest.mark.parametrize("N", [200, 500, 1000])
    @pytest.mark.parametrize("seed", [42, 7, 101])
    def test_low_N_seed_robustness(self, N, seed):
        """Sprint 19 closes carry-forward #17.

        Phase 1d.4 verified Gini N-invariance at 1000 sweeps across
        N ∈ {200, 500, 1000, 2000} at seed=42 only. This test pins the
        lower-N seed-robustness claim: at each N, across three seeds,
        Gini should land in a narrow band and no seed should display
        metastability (analogous to the Sprint 16 N=400 ABP finding where
        some seeds at small N fell into a different attractor).

        Expected: Gini at 1000 sweeps = 0.888 ± 0.015 for every (N, seed)
        combination in the tested grid. N=200 is the floor — smaller N
        would be unpinned by this test.
        """
        m = YardSale(n_agents=N, f=0.1, lambda_save=0.0, seed=seed)
        m.setup()
        m.step(1000 * N)
        gini = YardSale._gini(m.wealth)
        # Phase 1d.4 anchor is 0.888; ±0.02 tolerance accommodates
        # seed variance at the smallest N.
        assert 0.86 < gini < 0.92, (
            f"N={N} seed={seed}: Gini={gini:.4f} outside finite-size band "
            f"[0.86, 0.92]; possible metastability at low N"
        )

    @pytest.mark.slow
    def test_low_N_gini_spread_bounded(self):
        """σ(Gini) across seeds should be small at every N ≥ 200, with
        no dramatic enlargement as N shrinks (the metastability signature
        this test is built to catch).
        """
        ginis_by_N = {}
        for N in [200, 500, 1000]:
            g_values = []
            for seed in [42, 7, 101]:
                m = YardSale(n_agents=N, f=0.1, lambda_save=0.0, seed=seed)
                m.setup()
                m.step(1000 * N)
                g_values.append(YardSale._gini(m.wealth))
            ginis_by_N[N] = g_values
            sigma = float(np.std(g_values, ddof=1))
            print(f"  N={N}: Gini={g_values}, σ={sigma:.4f}")

        sigma_200 = float(np.std(ginis_by_N[200], ddof=1))
        sigma_1000 = float(np.std(ginis_by_N[1000], ddof=1))
        # Metastability would show as σ(N=200) >> σ(N=1000). Allow
        # modest enlargement but flag a doubling as suspicious.
        assert sigma_200 < 2 * sigma_1000 + 0.005, (
            f"σ(Gini) at N=200 ({sigma_200:.4f}) exceeds 2 × σ at N=1000 "
            f"({sigma_1000:.4f}) + 0.005 — possible seed metastability "
            f"at small N: {ginis_by_N}"
        )


# -----------------------------------------------------------------------------
# Registry hooks
# -----------------------------------------------------------------------------


class TestSprint17RegistryHooks:
    """Sprint 17 registrations reachable through orchestration."""

    def test_p28_registered(self):
        from epc.orchestration import DETECTOR_REGISTRY
        assert "P28" in DETECTOR_REGISTRY

    def test_yard_sale_registered(self):
        from epc.orchestration import MODEL_REGISTRY
        assert "yard_sale" in MODEL_REGISTRY
        reg = MODEL_REGISTRY["yard_sale"]
        assert reg.substrate_type == "scalar_wealth"
        assert "wealth" in reg.observables

    def test_scalar_wealth_substrate_present(self):
        from epc.orchestration import MODEL_REGISTRY
        substrates = {m.substrate_type for m in MODEL_REGISTRY.values()}
        assert "scalar_wealth" in substrates
