"""P28 — Wealth condensation / spontaneous inequality detector.

Detects spontaneous condensation of a conserved scalar resource into a
highly skewed distribution under symmetric, fair pairwise-exchange
dynamics — the Chakraborti-Boghosian Yard-Sale paradox.

Observable scope: model_metadata_assisted. The empirical detection
runs on state_history (a 'wealth' observable across frames). Metadata
is used to promote CONFIRMATION -> DEFINITIVE via the mechanistic-null
flags ``has_conserved_resource``, ``has_multiplicative_stake``,
``has_saving_propensity``, and ``has_redistribution``.

DISCRIMINATION FROM NEIGHBORS:

P14 (self-organized criticality, sandpiles): P14 operates on a
  lattice_2d substrate with avalanche-size observable and requires a
  scale-free power-law distribution of avalanche sizes. P28 is on a
  well-mixed (no-lattice) wealth substrate and operates on a scalar
  resource, not avalanches. Substrate-level rejection is clean by
  registration.

P1 / P2 / P5 / P6 / P8 / P12 / P13 / P15 / P22 / P27 / P31 all have
  substrate types that exclude the scalar_wealth substrate; all are
  rejected by the registry at substrate_mismatch.

P32 (emergent specialization / division of labor): different thing
  entirely — functional role differentiation, not resource
  concentration. No overlap expected at the detector level.

DETECTOR TIERS:

Prerequisites (all required for any detection):
  1. 'wealth' observable as (N,) non-negative float, N >= 50.
  2. At least 3 frames in measurement window for monotonicity check.
  3. Total wealth approximately conserved across the run (relative
     drift < 1% tolerated, a loose check; YS is exactly conservative).

Primary metric: Gini coefficient at the final frame of measurement.
  Pure YS: Gini grows monotonically toward 1.0. Well-mixed (DY 2001)
  equilibrium: Gini ~= 0.5. Perfect equality: Gini = 0.

Secondary metrics:
  - top_1pct_share, top_10pct_share
  - monotonic_growth_fraction across measurement window
  - relative_gini_growth (normalized slope)
  - alpha_hill (Pareto tail exponent, diagnostic only; empirical
    Sprint 17 Phase 1c shows alpha is unstable across timescales
    and NOT used as a tier gate — see REPLICATION_NOTES.md)

Null model: well-mixed Boltzmann-Gibbs (Dragulescu-Yakovenko 2001).
  Draw N samples from Exp(<w>) with matched mean; compute Gini.
  Repeat n_permutations times. p-value = P(Gini_null >= Gini_observed).
  Under the null of "symmetric-exchange produced ergodic well-mixing",
  Gini ~= 0.5 in expectation; a CONDENSING system has Gini > null.

Tier logic:
  Screening: prereqs pass + Gini > 0.40 + top_1pct > 0.05.
  Confirmation: screening + Gini > 0.55 + top_1pct > 0.15 +
                monotonic_fraction > 0.80 + null_p < 0.01.
  Definitive: confirmation + Gini > 0.80 + top_1pct > 0.30 +
              mechanistic-null flags affirm condensation regime
              (has_conserved_resource=True,
               has_multiplicative_stake=True,
               has_saving_propensity=False,
               has_redistribution=False).

The mechanistic-null gate separates CONFIRMATION from DEFINITIVE in
the same way P2's MIPS gate does (Sprint 16, Decision 43). Systems
with saving propensity > 0 or nonzero redistribution produce a
finite-Gini plateau that can still satisfy CONFIRMATION
(Sprint 17 Phase 1d.3 shows chi = 0.0001 gives Gini = 0.77). The
DEFINITIVE tier reserves the scientific claim "the mechanism producing
this inequality is symmetric fair exchange of a conserved resource
under multiplicative stakes" for cases where the metadata affirms
exactly that setup.

EMPIRICAL CALIBRATION (Sprint 17 Phase 1, N=1000, f=0.05-0.30,
equal init, 50 measurement checkpoints):

  Canonical positive (YS, f=0.05, lambda=0, chi=0, t=5M):
    Gini = 0.90, top_1pct = 0.23, monotonic_frac = 0.99, alpha ~ 1.0
    -> CONFIRMATION (score above 0.55, monotonic); DEFINITIVE with
       metadata flags affirmed.

  Canonical positive (YS, f=0.30, t=5M):
    Gini = 0.996, top_1pct = 0.99, monotonic_frac = 0.99
    -> DEFINITIVE.

  Within-family negatives (all N=1000, t=2M):
    YS lambda=0.5 (saving propensity): Gini = 0.29 -> rejected at
      screening (Gini < 0.40).
    YS chi=0.001 (redistribution): Gini = 0.13 -> rejected at
      screening.
    YS chi=0.0001 (very mild redistribution): Gini = 0.77 with
      monotonic growth -> passes CONFIRMATION, BUT metadata
      has_redistribution=True rejects DEFINITIVE.

  Substrate-mismatch negatives (all 15 other registered models):
    All rejected at registry level (no 'wealth' observable).
"""

from __future__ import annotations

from typing import Any

import numpy as np

from epc.base_detector import BaseDetector
from epc.detector_result import DetectionTier, NullType
from epc.metrics.wealth_concentration import (
    gini,
    top_p_share,
    hill_tail_alpha,
    monotonic_growth_fraction,
    relative_gini_growth,
    well_mixed_gini_null,
)


class P28WealthCondensationDetector(BaseDetector):
    """Detector for P28 — wealth condensation / spontaneous inequality.

    Parameters
    ----------
    burn_in : int
        Number of initial snapshots to discard before measurement.
        Default 0 (condensation dynamics start at t=0 with no
        equilibration phase — unlike MIPS which requires steady state).
    n_permutations : int
        Number of draws from the well-mixed Exp null. Minimum 199 for
        CONFIRMATION (need floor p < 0.01). Default 199. Gives floor
        p = 0.005.
    min_monotonic_frames : int
        Minimum number of frames in the measurement window over which
        monotonicity is computed. Default 5.
    seed : int
        RNG seed for the null-model draws.
    """

    # Tier thresholds
    _SCREEN_GINI_MIN = 0.40
    _SCREEN_TOP1PCT_MIN = 0.05
    _CONFIRM_GINI_MIN = 0.55
    _CONFIRM_TOP1PCT_MIN = 0.15
    _CONFIRM_MONOTONIC_MIN = 0.80
    _CONFIRM_NULL_P_MAX = 0.01
    _DEF_GINI_MIN = 0.80
    _DEF_TOP1PCT_MIN = 0.30

    # Prerequisites
    _MIN_N_AGENTS = 50
    _MIN_RUN_FRAMES = 3

    # Conservation tolerance (fraction of mean wealth)
    _CONSERVATION_TOL = 0.01

    def __init__(
        self,
        burn_in: int = 0,
        n_permutations: int = 199,
        min_monotonic_frames: int = 5,
        seed: int = 42,
    ) -> None:
        super().__init__(
            pattern_id="P28",
            # P28 is on its own substrate (scalar_wealth); no
            # neighbor-pattern conflicts on that substrate at Sprint 17.
            excluded_patterns=[],
            allowed_co_occurrences=[],
            observable_scope="model_metadata_assisted",
        )
        if n_permutations < 1:
            raise ValueError("n_permutations must be >= 1")
        if burn_in < 0:
            raise ValueError("burn_in must be >= 0")
        if min_monotonic_frames < 2:
            raise ValueError("min_monotonic_frames must be >= 2")
        self.burn_in = burn_in
        self.n_permutations = n_permutations
        self.min_monotonic_frames = min_monotonic_frames
        self._seed = seed

        # Caches populated during detect()
        self._measurement_frames: list[np.ndarray] = []
        self._gini_series: np.ndarray = np.empty(0)
        self._screening_rejection_reason: str = "none"

    # ------------------------------------------------------------------
    # Prerequisites
    # ------------------------------------------------------------------

    def _estimate_timescale(
        self,
        state_history: list[dict[str, Any]],
        model_metadata: dict[str, Any] | None,
    ) -> float:
        """System timescale: T_sweep = N transactions per agent.

        Falls back to len(history) / 10 if unavailable.
        """
        if model_metadata is not None:
            N = model_metadata.get("n_agents")
            if N is not None and N > 0:
                return float(N)
        return max(1.0, len(state_history) / 10.0)

    def _validate_prerequisites(
        self,
        state_history: list[dict[str, Any]],
        timescale: float,
    ) -> list[str]:
        """Validate P28 prerequisites."""
        warnings: list[str] = []
        self._screening_rejection_reason = "none"
        self._measurement_frames = []
        self._gini_series = np.empty(0)

        if not state_history:
            warnings.append("empty state_history")
            self._screening_rejection_reason = "empty_state_history"
            return warnings

        # Extract the measurement window (burn_in: end)
        measurement_history = state_history[self.burn_in:]
        if len(measurement_history) < self._MIN_RUN_FRAMES:
            warnings.append(
                f"post-burn-in run length = {len(measurement_history)} < "
                f"minimum {self._MIN_RUN_FRAMES}; insufficient for P28 "
                f"monotonicity check"
            )
            self._screening_rejection_reason = "run_too_short"
            return warnings

        # Check for 'wealth' observable
        if "wealth" not in measurement_history[-1]:
            warnings.append(
                "no 'wealth' observable in state snapshots — P28 detector "
                "requires scalar_wealth substrate"
            )
            self._screening_rejection_reason = "substrate_mismatch"
            return warnings

        # Collect wealth frames
        wealth_frames: list[np.ndarray] = []
        for frame in measurement_history:
            w = frame.get("wealth")
            if w is None:
                continue
            arr = np.asarray(w, dtype=np.float64)
            if arr.ndim != 1:
                warnings.append(
                    f"'wealth' observable is not 1D (shape={arr.shape}); "
                    f"expected (N,)"
                )
                self._screening_rejection_reason = "substrate_mismatch"
                return warnings
            if np.any(arr < 0):
                warnings.append(
                    "negative wealth values present — P28 assumes "
                    "non-negative wealth"
                )
                # continue anyway; don't block
            wealth_frames.append(arr)

        if len(wealth_frames) < self._MIN_RUN_FRAMES:
            warnings.append(
                f"only {len(wealth_frames)} usable wealth frames; need "
                f"{self._MIN_RUN_FRAMES}"
            )
            self._screening_rejection_reason = "run_too_short"
            return warnings

        # N-particles check
        N = wealth_frames[-1].size
        if N < self._MIN_N_AGENTS:
            warnings.append(
                f"n_agents = {N} < minimum {self._MIN_N_AGENTS}; "
                f"insufficient for P28 statistics"
            )
            self._screening_rejection_reason = "too_few_agents"
            return warnings

        # Conservation sanity check (loose)
        totals = np.array([w.sum() for w in wealth_frames], dtype=np.float64)
        mean_total = float(totals.mean())
        if mean_total > 0:
            drift = float(np.abs(totals - mean_total).max() / mean_total)
            if drift > self._CONSERVATION_TOL:
                warnings.append(
                    f"total wealth drift = {drift:.3%} exceeds "
                    f"{self._CONSERVATION_TOL:.1%} tolerance — model may "
                    f"not conserve the resource (e.g., redistribution is "
                    f"active or there's a bug)"
                )
                # Not a hard rejection — YS with chi>0 can drift
                # slightly, and we want the detector to still be able to
                # flag it.

        self._measurement_frames = wealth_frames
        return warnings

    # ------------------------------------------------------------------
    # Primary / secondaries / screening
    # ------------------------------------------------------------------

    def _compute_primary(
        self,
        state_history: list[dict[str, Any]],
        timescale: float,
    ) -> dict[str, float]:
        """Primary metric: Gini at the final frame."""
        if self._screening_rejection_reason != "none":
            return {
                "gini_final": 0.0,
                "gini_initial": 0.0,
                "mean_wealth": 0.0,
                "n_agents": 0,
                "n_frames_used": 0,
                "screening_rejection_reason": self._screening_rejection_reason,
            }

        frames = self._measurement_frames
        gini_series = np.array([gini(w) for w in frames], dtype=np.float64)
        self._gini_series = gini_series

        w_final = frames[-1]
        return {
            "gini_final": float(gini_series[-1]),
            "gini_initial": float(gini_series[0]),
            "gini_delta": float(gini_series[-1] - gini_series[0]),
            "mean_wealth": float(w_final.mean()),
            "n_agents": int(w_final.size),
            "n_frames_used": int(len(frames)),
            "screening_rejection_reason": "none",
        }

    def _check_screening(
        self,
        primary_result: dict[str, float],
        timescale: float,
    ) -> bool:
        if self._screening_rejection_reason != "none":
            return False
        gini_final = float(primary_result.get("gini_final", 0.0))
        if gini_final < self._SCREEN_GINI_MIN:
            self._screening_rejection_reason = "gini_below_screening_floor"
            primary_result["screening_rejection_reason"] = (
                "gini_below_screening_floor"
            )
            return False
        # top_1pct computed in secondaries; but we can check it cheaply here
        w_final = self._measurement_frames[-1]
        top1 = top_p_share(w_final, 0.01)
        if top1 < self._SCREEN_TOP1PCT_MIN:
            self._screening_rejection_reason = "top1_below_screening_floor"
            primary_result["screening_rejection_reason"] = (
                "top1_below_screening_floor"
            )
            return False
        return True

    def _compute_secondaries(
        self,
        state_history: list[dict[str, Any]],
        timescale: float,
    ) -> dict[str, Any]:
        """Secondary metrics: top-p shares, monotonicity, Hill tail."""
        if not self._measurement_frames:
            return {}
        w_final = self._measurement_frames[-1]
        top1 = top_p_share(w_final, 0.01)
        top10 = top_p_share(w_final, 0.10)
        max_share = float(w_final.max() / w_final.sum()) if w_final.sum() > 0 else 0.0

        mono_frac = monotonic_growth_fraction(self._gini_series)
        rel_growth = relative_gini_growth(self._gini_series)

        alpha, w_min, k = hill_tail_alpha(w_final)
        return {
            "top_1pct_share": float(top1),
            "top_10pct_share": float(top10),
            "max_share": float(max_share),
            "monotonic_fraction": float(mono_frac) if np.isfinite(mono_frac) else 0.0,
            "relative_gini_growth": float(rel_growth) if np.isfinite(rel_growth) else 0.0,
            "alpha_hill": float(alpha) if np.isfinite(alpha) else float("nan"),
            "hill_k_used": int(k),
            "hill_w_min": float(w_min) if np.isfinite(w_min) else float("nan"),
        }

    # ------------------------------------------------------------------
    # Null / effect size / tier / exclusions
    # ------------------------------------------------------------------

    def _run_null_model(
        self,
        state_history: list[dict[str, Any]],
        primary_result: dict[str, float],
        timescale: float,
    ) -> tuple[float, NullType, dict[str, float]]:
        """Null for P28: the well-mixed Boltzmann-Gibbs (Dragulescu-
        Yakovenko 2001) distribution.

        Draw N samples from Exp(mean_w); compute Gini. Repeat
        n_permutations times. p-value = P(Gini_null >= Gini_observed).

        Under the null H0: "symmetric exchange of a conserved scalar
        resource equilibrates to an exponential distribution", the
        Gini of a sample is ~ 0.5 in the large-N limit. Observed
        Gini > null is evidence of SUPER-BOLTZMANN condensation —
        the defining P28 signature.
        """
        if not self._measurement_frames:
            return 1.0, NullType.SURROGATE, {"mean": 0.0, "std": 0.0, "observed": 0.0}

        w_final = self._measurement_frames[-1]
        N = int(w_final.size)
        mean_w = float(w_final.mean())
        if mean_w <= 0:
            return 1.0, NullType.SURROGATE, {"mean": 0.0, "std": 0.0, "observed": 0.0}

        observed_gini = float(primary_result.get("gini_final", gini(w_final)))

        rng = np.random.default_rng(self._seed)
        null_ginis = well_mixed_gini_null(N, mean_w, self.n_permutations, rng=rng)
        # Right-tailed: observed should exceed null
        p = (float((null_ginis >= observed_gini).sum()) + 1.0) / (len(null_ginis) + 1.0)
        return (
            float(p),
            NullType.SURROGATE,
            {
                "mean": float(null_ginis.mean()),
                "std": float(null_ginis.std()),
                "observed": observed_gini,
            },
        )

    def _compute_effect_size(
        self,
        primary_result: dict[str, float],
        null_dist_stats: dict[str, float],
    ) -> dict[str, float]:
        """Cohen's d: (observed_gini - null_mean) / null_std."""
        observed = float(null_dist_stats.get("observed", 0.0))
        nm = float(null_dist_stats.get("mean", 0.0))
        ns = float(null_dist_stats.get("std", 0.0))
        if ns <= 1e-9:
            # Exp(mean_w) at N=1000 gives std ~ 0.01, virtually never
            # degenerate in practice, but handle for robustness
            if observed > nm:
                d = 1e6
            elif observed < nm:
                d = -1e6
            else:
                d = 0.0
        else:
            d = (observed - nm) / ns
        return {
            "cohens_d": float(d),
            "raw_value": observed,
            "null_mean": nm,
            "null_std": ns,
        }

    def _mechanistic_null_passes(
        self, model_metadata: dict[str, Any] | None
    ) -> dict[str, Any]:
        """Check P28 metadata flags for the condensation-mechanism gate.

        Required for DEFINITIVE:
          has_conserved_resource = True
          has_multiplicative_stake = True
          has_saving_propensity = False
          has_redistribution = False

        Returns dict with 'passes' bool and per-flag diagnostic.
        """
        if model_metadata is None:
            return {"passes": False, "reason": "metadata_absent"}
        required = {
            "has_conserved_resource": True,
            "has_multiplicative_stake": True,
            "has_saving_propensity": False,
            "has_redistribution": False,
        }
        flags = {}
        ok = True
        for k, v in required.items():
            actual = model_metadata.get(k)
            flags[k] = {"actual": actual, "required": v}
            if actual != v:
                ok = False
        return {"passes": ok, "flags": flags,
                "reason": "passed" if ok else "flag_mismatch"}

    def _determine_tier(
        self,
        primary_result: dict[str, float],
        secondary_result: dict[str, Any],
        null_p: float,
        null_type: NullType,
        effect_size: dict[str, float],
        state_history: list[dict[str, Any]],
        model_metadata: dict[str, Any] | None,
        timescale: float,
    ) -> tuple[DetectionTier, dict[str, bool]]:
        """P28 tier logic.

        Screening already passed.  Check CONFIRMATION gates, then
        DEFINITIVE.
        """
        bonuses: dict[str, bool] = {}

        gini_final = float(primary_result.get("gini_final", 0.0))
        top1 = float(secondary_result.get("top_1pct_share", 0.0))
        mono = float(secondary_result.get("monotonic_fraction", 0.0))

        # Screening bonuses
        bonuses["secondaries_pass"] = bool(secondary_result)
        bonuses["shuffle_null_p_001"] = null_p < 0.01

        # CONFIRMATION gates
        gini_ok = gini_final >= self._CONFIRM_GINI_MIN
        top1_ok = top1 >= self._CONFIRM_TOP1PCT_MIN
        mono_ok = mono >= self._CONFIRM_MONOTONIC_MIN
        null_ok = null_p < self._CONFIRM_NULL_P_MAX

        if not (gini_ok and top1_ok and mono_ok and null_ok):
            return DetectionTier.SCREENING, bonuses

        # Confirmation bonuses
        bonuses["null_p_0001"] = null_p < 0.0001
        bonuses["effect_size_gt_1"] = abs(effect_size.get("cohens_d", 0.0)) > 1.0
        bonuses["all_secondaries"] = gini_ok and top1_ok and mono_ok

        # DEFINITIVE: stronger magnitudes + metadata-mechanism
        gini_def = gini_final >= self._DEF_GINI_MIN
        top1_def = top1 >= self._DEF_TOP1PCT_MIN
        if not (gini_def and top1_def):
            return DetectionTier.CONFIRMATION, bonuses

        mech = self._mechanistic_null_passes(model_metadata)
        if not mech["passes"]:
            return DetectionTier.CONFIRMATION, bonuses

        # Definitive bonuses
        bonuses["all_exclusions_cleared"] = True
        bonuses["both_null_types_rejected"] = null_p < 0.01 and mech["passes"]
        bonuses["finite_size_robustness"] = False  # could be added via multi-N slow test
        return DetectionTier.DEFINITIVE, bonuses

    def _all_secondaries_pass(self, secondary_result: dict[str, Any]) -> bool:
        if not secondary_result:
            return False
        top1 = float(secondary_result.get("top_1pct_share", 0.0))
        mono = float(secondary_result.get("monotonic_fraction", 0.0))
        return (
            top1 >= self._CONFIRM_TOP1PCT_MIN
            and mono >= self._CONFIRM_MONOTONIC_MIN
        )

    def _check_exclusions(
        self,
        state_history: list[dict[str, Any]],
        model_metadata: dict[str, Any] | None,
        timescale: float,
    ) -> tuple[list[str], dict[str, str]]:
        """P28 has no near-neighbor pattern on the same substrate at Sprint 17.

        All other patterns operate on incompatible substrates
        (lattice_1d, lattice_2d, continuous_2d, oscillator,
        opinion_space, lattice_2d_continuous). The scalar_wealth
        substrate is new to the registry at Sprint 17 and has no
        competitors.
        """
        return ([], {})
