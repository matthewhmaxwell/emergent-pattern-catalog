"""End-to-end tests for non-local Kuramoto + P10 chimera detector.

Coverage:
- Canonical positive (N=128, A=0.995, β=0.05, asymmetric IC, seed=0) -> DEFINITIVE
- Seed robustness across the chimera basin (β=0.05 multi-seed) -> DEFINITIVE
- Sync basin (β=0.10 seed=0) -> screening (no_coexistence)
- Ordinary Kuramoto (K=0.3, 1.0, 2.0, 4.0) -> screening rejection
- Random phases -> screening rejection (no_coexistence)
- Short run -> screening rejection (run_too_short)
- Small N -> screening rejection (too_few_oscillators)
- Missing theta observable -> screening rejection (substrate_mismatch)
- Metadata: DEFINITIVE requires has_nonlocal_coupling=True
           AND has_frequency_heterogeneity != True.
- Registry: both model and detector registered correctly.

Timing budget:
- Fast-half: each test <= ~8s, total ~45s.
- Slow-half: 3-4 seeds at larger N or longer T for replication-quality checks.
"""

from __future__ import annotations

import numpy as np
import pytest

from epc.detector_result import DetectionTier
from epc.detectors.p10_chimera import P10ChimeraDetector
from epc.models.kuramoto import KuramotoModel, KuramotoParams
from epc.models.kuramoto_nonlocal import KuramotoNonlocal


# -----------------------------------------------------------------------------
# Helpers
# -----------------------------------------------------------------------------


def _run_chimera_p10(seed=0, beta=0.05, N=128, n_frames=50,
                     n_permutations=199, A=0.995):
    """Build a non-local Kuramoto, run, apply P10, return DetectorResult."""
    m = KuramotoNonlocal(N=N, A=A, beta=beta, seed=seed)
    hist = m.run(n_frames=n_frames)
    det = P10ChimeraDetector(n_permutations=n_permutations, seed=42)
    return det.detect(hist, model_metadata=m.get_metadata())


def _run_ordinary_kuramoto_p10(K=1.0, seed=42, N=128, n_steps=4000,
                                equilibration=2000, record_every=40,
                                n_permutations=199):
    """Build an ordinary all-to-all Kuramoto, run, apply P10."""
    params = KuramotoParams(N=N, K=K, gamma=0.5, dt=0.05, seed=seed,
                             freq_dist='lorentzian')
    m = KuramotoModel(params)
    hist = m.run(n_steps=n_steps, record_every=record_every,
                  equilibration=equilibration)
    det = P10ChimeraDetector(n_permutations=n_permutations, seed=42)
    return det.detect(hist, model_metadata=m.get_metadata())


# -----------------------------------------------------------------------------
# Canonical positive
# -----------------------------------------------------------------------------


class TestCanonicalChimera:
    """P10 tier assignments on the canonical Abrams-Strogatz chimera regime.

    Fast-half sizing: N=128, 50 frames, n_permutations=199.

    Phase 1j calibration at these parameters:
      pos_vel_ac[4] ≈ 0.84–0.86  across seeds {0, 1, 42, 200, 500}
      Chimera null (position-shuffle, n=199) gives floor p = 0.005.
      Cohen's d on (observed - null mean) / null std ~ 9.
    """

    def test_canonical_positive_is_definitive(self):
        """seed=0, β=0.05, N=128, T=50 -> DEFINITIVE."""
        result = _run_chimera_p10(seed=0, beta=0.05, N=128, n_frames=50,
                                   n_permutations=199)
        assert result.detected, (
            f"canonical chimera must be detected; got "
            f"pos_vel_ac={result.primary_metric.get('pos_vel_ac'):.3f}"
        )
        assert result.tier == DetectionTier.DEFINITIVE, (
            f"canonical chimera must be DEFINITIVE, got {result.tier}; "
            f"primary={result.primary_metric}"
        )
        # Primary metric band
        pva = result.primary_metric["pos_vel_ac"]
        assert 0.75 < pva < 1.0, (
            f"pos_vel_ac[4] should be in chimera band, got {pva:.3f}"
        )
        # Coexistence must have been established
        assert result.primary_metric["coexistence"] is True
        assert result.primary_metric["n_persistent_coh"] >= 1
        assert result.primary_metric["n_persistent_incoh"] >= 1
        # Null passed at < 0.01
        assert result.null_p_value < 0.01, (
            f"expected null_p < 0.01, got {result.null_p_value}"
        )
        # Effect size is very large
        assert result.effect_size["cohens_d"] > 2.0, (
            f"expected strong effect, got d={result.effect_size['cohens_d']:.2f}"
        )
        # Confidence at DEFINITIVE floor
        assert result.confidence >= 0.75

    @pytest.mark.parametrize("seed", [1, 42, 200])
    def test_multi_seed_chimera_is_definitive(self, seed):
        """β=0.05 with seeds {1, 42, 200} -> DEFINITIVE.

        At β=0.05 the chimera basin is wide (5/5 tested seeds in Phase 1j
        produced chimeras). seed=500 is also a positive but slow.
        """
        result = _run_chimera_p10(seed=seed, beta=0.05, N=128, n_frames=50,
                                   n_permutations=199)
        assert result.tier == DetectionTier.DEFINITIVE, (
            f"seed={seed} should give DEFINITIVE, got {result.tier}; "
            f"pos_vel_ac={result.primary_metric.get('pos_vel_ac'):.3f}"
        )

    def test_canonical_secondary_metrics_in_chimera_band(self):
        """r_global in (0.3, 0.95), velocity_neighbor_corr > 0, persistence_corr > 0."""
        result = _run_chimera_p10(seed=0, beta=0.05, n_permutations=199)
        assert 0.3 < result.secondary_metrics["r_global_mean"] < 0.95, (
            f"r_global should be partial-sync, got "
            f"{result.secondary_metrics['r_global_mean']:.3f}"
        )
        assert result.secondary_metrics["velocity_neighbor_corr"] > 0.0


# -----------------------------------------------------------------------------
# Content-level rejections (same substrate, different regime)
# -----------------------------------------------------------------------------


class TestNegativeOnSameSubstrate:
    """Ordinary all-to-all Kuramoto (same substrate as kuramoto_nonlocal).

    These runs all have ``theta`` observable; P10 must reject them at
    content level via coexistence and/or pos_vel_ac gates — NOT at
    substrate_mismatch. See Phase 1j: Kuramoto K=1.0 max pos_vel_ac
    across 6 seeds was 0.448, safely below the 0.55 screening floor.
    """

    @pytest.mark.parametrize("K", [0.3, 1.0, 2.0, 4.0])
    def test_ordinary_kuramoto_rejected_at_screening(self, K):
        result = _run_ordinary_kuramoto_p10(K=K, seed=42, n_permutations=99)
        assert not result.detected, (
            f"ordinary Kuramoto K={K} must NOT be detected; "
            f"got pos_vel_ac={result.primary_metric.get('pos_vel_ac'):.3f}"
        )
        assert result.tier == DetectionTier.SCREENING, (
            f"expected SCREENING tier (rejection), got {result.tier}"
        )
        srr = result.primary_metric.get("screening_rejection_reason", "none")
        assert srr in ("no_coexistence", "pos_vel_ac_below_floor"), (
            f"unexpected rejection reason: {srr}"
        )
        pva = result.primary_metric.get("pos_vel_ac", 0.0)
        assert pva < 0.55, (
            f"ordinary Kuramoto K={K} primary must be below screening floor, "
            f"got {pva:.3f}"
        )

    def test_nonlocal_sync_basin_rejected_at_coexistence(self):
        """β=0.10 seed=0 lands in sync basin -> no_coexistence rejection.

        Note: at β=0.10 the chimera basin is narrow; different seeds may
        still produce chimeras. seed=0 reliably escapes to full sync
        (Phase 1c); seed=42 at β=0.10 actually finds the chimera basin.
        Only seed=0 is safe for this negative test.
        """
        result = _run_chimera_p10(seed=0, beta=0.10, N=128, n_frames=50,
                                   n_permutations=99)
        assert not result.detected
        assert result.tier == DetectionTier.SCREENING
        assert result.primary_metric.get("screening_rejection_reason") \
            == "no_coexistence"
        # All windows should be coherent (full sync)
        assert result.primary_metric["n_persistent_coh"] >= 12
        assert result.primary_metric["n_persistent_incoh"] == 0


# -----------------------------------------------------------------------------
# Prerequisite rejections
# -----------------------------------------------------------------------------


class TestPrerequisiteRejections:
    """Screening gates: too small, too short, missing observables."""

    def test_random_phases_rejected(self):
        """Random-phase snapshots with no ring structure -> screening."""
        rng = np.random.default_rng(0)
        hist = [{"theta": rng.uniform(0, 2 * np.pi, 128)} for _ in range(60)]
        det = P10ChimeraDetector(n_permutations=99, seed=42)
        result = det.detect(hist, model_metadata=None)
        assert not result.detected
        assert result.tier == DetectionTier.SCREENING
        # Either no_coexistence (usual) or pos_vel_ac_below_floor
        srr = result.primary_metric.get("screening_rejection_reason")
        assert srr in ("no_coexistence", "pos_vel_ac_below_floor")

    def test_missing_theta_observable_substrate_mismatch(self):
        """State history without `theta` -> substrate_mismatch."""
        hist = [{"grid": np.zeros((10, 10))} for _ in range(60)]
        det = P10ChimeraDetector(n_permutations=49, seed=42)
        result = det.detect(hist, model_metadata={"model_family": "foo"})
        assert not result.detected
        assert result.primary_metric.get("screening_rejection_reason") \
            == "substrate_mismatch"

    def test_short_run_rejected(self):
        """Post-burn < 30 frames -> run_too_short."""
        m = KuramotoNonlocal(N=128, seed=0)
        hist = m.run(n_frames=20)  # post-burn=14 < 30
        det = P10ChimeraDetector(n_permutations=49, seed=42)
        result = det.detect(hist, model_metadata=m.get_metadata())
        assert not result.detected
        assert result.primary_metric.get("screening_rejection_reason") \
            == "run_too_short"

    def test_small_N_rejected(self):
        """N < 32 -> too_few_oscillators."""
        m = KuramotoNonlocal(N=16, seed=0)
        hist = m.run(n_frames=50)
        det = P10ChimeraDetector(n_permutations=49, seed=42)
        result = det.detect(hist, model_metadata=m.get_metadata())
        assert not result.detected
        assert result.primary_metric.get("screening_rejection_reason") \
            == "too_few_oscillators"


# -----------------------------------------------------------------------------
# Metadata-gate tests (DEFINITIVE depends on metadata)
# -----------------------------------------------------------------------------


class TestDefinitiveMetadataGate:
    """DEFINITIVE requires has_nonlocal_coupling=True in metadata.

    The mechanistic-null gate follows the Sprint 16 P2 / Sprint 17 P28
    pattern — content-level signal alone can only reach CONFIRMATION
    without metadata confirmation of the underlying mechanism.
    """

    def test_definitive_requires_nonlocal_coupling_metadata(self):
        """Same chimera data but with has_nonlocal_coupling flag removed
        -> caps at CONFIRMATION."""
        m = KuramotoNonlocal(seed=0)
        hist = m.run(n_frames=50)
        meta = m.get_metadata()
        # Strip the mechanism flag
        meta = {k: v for k, v in meta.items()
                if k not in ("has_nonlocal_coupling",
                             "has_frequency_heterogeneity")}
        det = P10ChimeraDetector(n_permutations=199, seed=42)
        result = det.detect(hist, model_metadata=meta)
        assert result.detected
        assert result.tier == DetectionTier.CONFIRMATION, (
            f"without metadata gate, expected CONFIRMATION, got {result.tier}"
        )

    def test_definitive_blocked_by_frequency_heterogeneity_flag(self):
        """Chimera data but metadata claims has_frequency_heterogeneity=True
        -> caps at CONFIRMATION (this claim contradicts the ring model's
        identical ω)."""
        m = KuramotoNonlocal(seed=0)
        hist = m.run(n_frames=50)
        meta = dict(m.get_metadata())
        meta["has_frequency_heterogeneity"] = True  # override
        det = P10ChimeraDetector(n_permutations=199, seed=42)
        result = det.detect(hist, model_metadata=meta)
        assert result.detected
        assert result.tier == DetectionTier.CONFIRMATION

    def test_no_metadata_caps_at_confirmation(self):
        """model_metadata=None -> caps at CONFIRMATION."""
        m = KuramotoNonlocal(seed=0)
        hist = m.run(n_frames=50)
        det = P10ChimeraDetector(n_permutations=199, seed=42)
        result = det.detect(hist, model_metadata=None)
        assert result.detected
        assert result.tier == DetectionTier.CONFIRMATION


# -----------------------------------------------------------------------------
# Permutation-count floor
# -----------------------------------------------------------------------------


class TestPermutationFloor:
    """n_permutations < 199 does not reach CONFIRMATION on canonical data.

    With n_permutations=99, floor p = 1 / 100 = 0.01, which FAILS the
    `< 0.01` confirmation gate. This is the same trap as Sprint 15 P8
    and Sprint 16 P2.
    """

    def test_n_perm_99_caps_at_screening(self):
        """n_permutations=99 -> floor p = 0.01, fails confirmation gate."""
        result = _run_chimera_p10(seed=0, beta=0.05, n_permutations=99)
        assert result.detected, "should still clear screening"
        assert result.null_p_value == pytest.approx(0.01, abs=1e-6), (
            f"floor p for n_perm=99 should be 0.01, got {result.null_p_value}"
        )
        # With floor p = 0.01, the `null_p < 0.01` confirmation gate fails,
        # so tier caps at SCREENING.
        assert result.tier == DetectionTier.SCREENING


# -----------------------------------------------------------------------------
# Registry
# -----------------------------------------------------------------------------


class TestRegistry:
    """Both kuramoto_nonlocal and P10 must be registered and compatible."""

    def test_kuramoto_nonlocal_registered(self):
        from epc.orchestration import MODEL_REGISTRY
        assert "kuramoto_nonlocal" in MODEL_REGISTRY
        reg = MODEL_REGISTRY["kuramoto_nonlocal"]
        assert reg.substrate_type == "oscillator"
        assert "theta" in reg.observables
        assert "P10" in reg.primary_patterns

    def test_p10_detector_registered(self):
        from epc.orchestration import DETECTOR_REGISTRY
        assert "P10" in DETECTOR_REGISTRY
        reg = DETECTOR_REGISTRY["P10"]
        assert "oscillator" in reg.required_substrate
        assert "theta" in reg.required_observables

    def test_p10_only_matches_oscillator_substrate(self):
        from epc.orchestration import MODEL_REGISTRY, check_compatibility
        compatible_models = [
            m for m in MODEL_REGISTRY
            if check_compatibility(m, "P10").compatible
        ]
        # Only kuramoto and kuramoto_nonlocal
        assert set(compatible_models) == {"kuramoto", "kuramoto_nonlocal"}

    def test_kuramoto_nonlocal_matches_p9_and_p10_only(self):
        from epc.orchestration import DETECTOR_REGISTRY, check_compatibility
        compatible_dets = [
            d for d in DETECTOR_REGISTRY
            if check_compatibility("kuramoto_nonlocal", d).compatible
        ]
        assert set(compatible_dets) == {"P9", "P10"}


# -----------------------------------------------------------------------------
# SLOW replication-quality tests
# -----------------------------------------------------------------------------


@pytest.mark.slow
class TestSlowReplication:
    """Larger-N and longer-T replication checks. Run with -m slow."""

    def test_N_256_chimera_is_definitive(self):
        """N=256 at β=0.05 still gives DEFINITIVE.

        Phase 1e N-scaling: chimera structure holds at N=256 with
        pos_vel_ac[4] ~ 0.91 and the basin is still wide at β=0.05.
        """
        result = _run_chimera_p10(seed=0, beta=0.05, N=256, n_frames=50,
                                   n_permutations=199)
        assert result.tier == DetectionTier.DEFINITIVE, (
            f"N=256 should give DEFINITIVE, got {result.tier}; "
            f"pos_vel_ac={result.primary_metric.get('pos_vel_ac'):.3f}"
        )

    def test_long_run_T100_chimera_is_definitive(self):
        """T=100 (long-run stability) -> DEFINITIVE."""
        result = _run_chimera_p10(seed=0, beta=0.05, N=128, n_frames=100,
                                   n_permutations=199)
        assert result.tier == DetectionTier.DEFINITIVE

    def test_paper_beta_018_chimera_basin(self):
        """β=0.18 (paper value) chimera basin. Seed=42 lands in the
        chimera basin; Phase 1j confirmed DEFINITIVE at N=128."""
        result = _run_chimera_p10(seed=42, beta=0.18, N=128, n_frames=50,
                                   n_permutations=199)
        assert result.tier == DetectionTier.DEFINITIVE, (
            f"seed=42 β=0.18 should reach DEFINITIVE (paper parameters); "
            f"got {result.tier}, "
            f"pos_vel_ac={result.primary_metric.get('pos_vel_ac'):.3f}"
        )


# -----------------------------------------------------------------------------
# Slow-half: finite-size robustness (Sprint 19 — closes carry-forward #22)
# -----------------------------------------------------------------------------


@pytest.mark.slow
class TestFiniteSizeRobustness:
    """Finite-size robustness test pinning the lower-N bound of the
    chimera basin at β=0.05 across multiple seeds.

    Sprint 18 Phase 1e characterized pos_vel_ac[4] at seed=0 across
    N ∈ {64, 128, 256}; Phase 1j confirmed all five seeds
    {0, 1, 42, 200, 500} reach the chimera basin at β=0.05 for N=128.
    This test combines the two: at each N in {128, 256}, run all
    five chimera-basin seeds and require DEFINITIVE tier at every
    (N, seed) pair.

    Sprint 19 finding (added when this test was written): at N=64 the
    chimera basin is NOT robust across the Phase-1j seed set. A clean
    run across the five anchored seeds at N=64, β=0.05 gives
    pos_vel_ac[4] ∈ [−0.008, 0.59] — five SCREENING outcomes with no
    DEFINITIVE. The Sprint 18 Phase 1e claim "pos_vel_ac[4] = 0.863 at
    N=64" holds for a specific seed-IC configuration that does not
    generalize. N=64 is therefore below the seed-robust chimera basin
    floor; N=128 is the true lower bound pinned by this test. The
    discrepancy is documented as a Sprint 19 scientific finding for
    the Sprint 20 transfer prompt to carry forward.

    Runtime: ~4–7 minutes total for 10 runs (2 N-values × 5 seeds) at
    n_permutations=199.
    """

    @pytest.mark.parametrize("N", [128, 256])
    @pytest.mark.parametrize("seed", [0, 1, 42, 200, 500])
    def test_chimera_basin_across_N_and_seed(self, N, seed):
        """Canonical β=0.05 regime reaches DEFINITIVE at every tested
        (N, seed) combination for N ≥ 128. Failure of any case would
        indicate that the chimera basin narrows further than expected
        at N=128 for some specific seeds."""
        result = _run_chimera_p10(seed=seed, beta=0.05, N=N, n_frames=50,
                                   n_permutations=199)
        assert result.tier == DetectionTier.DEFINITIVE, (
            f"N={N} seed={seed} β=0.05: expected DEFINITIVE, got "
            f"{result.tier}; pos_vel_ac="
            f"{result.primary_metric.get('pos_vel_ac'):.3f}. "
            f"Chimera basin regression at N ≥ 128."
        )

    def test_pos_vel_ac_N_invariance(self):
        """pos_vel_ac[lag=4] averaged across the five chimera-basin
        seeds should be tight across N ∈ {128, 256} (within ±0.10) —
        the structural claim that the chimera fixed point is a
        continuum attractor whose dimensionless order parameter is
        genuinely N-invariant above the finite-size basin-shrinkage
        regime.
        """
        means = {}
        for N in [128, 256]:
            vals = []
            for seed in [0, 1, 42, 200, 500]:
                result = _run_chimera_p10(seed=seed, beta=0.05, N=N,
                                           n_frames=50, n_permutations=199)
                vals.append(float(result.primary_metric["pos_vel_ac"]))
            means[N] = float(np.mean(vals))
            print(f"  N={N}: pos_vel_ac[4] per-seed = "
                  f"{[f'{v:.3f}' for v in vals]}, mean={means[N]:.3f}")

        m_values = list(means.values())
        spread = max(m_values) - min(m_values)
        assert spread < 0.10, (
            f"pos_vel_ac[4] spread {spread:.3f} across N={list(means.keys())} "
            f"exceeds 0.10 — finite-size variation larger than expected "
            f"for a continuum attractor: {means}"
        )

    @pytest.mark.parametrize("seed", [0, 1, 42, 200, 500])
    def test_N64_below_robust_basin(self, seed):
        """N=64 is below the seed-robust chimera basin floor.

        Pinned as a scientific finding discovered during Sprint 19
        cleanup: at N=64, β=0.05, none of the five Phase-1j anchored
        seeds reach DEFINITIVE (pos_vel_ac[4] ≤ 0.59 for all five).
        Future detector improvements that extend the robust basin
        below N=128 should make THIS test fail (flipping from
        SCREENING to DEFINITIVE); that failure would be a scientific
        finding, not a regression.
        """
        result = _run_chimera_p10(seed=seed, beta=0.05, N=64, n_frames=50,
                                   n_permutations=199)
        # The empirical upper bound observed across the 5 seeds at N=64
        # is pos_vel_ac ≈ 0.59. Any value exceeding 0.70 would mean the
        # finite-size basin has improved materially — worth investigating.
        pv = float(result.primary_metric["pos_vel_ac"])
        assert pv < 0.70, (
            f"N=64 seed={seed}: pos_vel_ac={pv:.3f} exceeds 0.70 floor. "
            f"If this is a detector improvement, celebrate and update "
            f"the seed-robust basin lower bound from N=128 to N=64."
        )
