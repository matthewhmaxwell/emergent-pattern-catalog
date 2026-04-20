"""P8 — traffic jamming detector.

Detects emergent stop-and-go jamming in a 1D traffic-flow cellular
automaton (Nagel-Schreckenberg family).

Observable scope: state_history_only. The detector reads the 'velocities'
observable (integer-valued 1D array, values in 0..v_max) from each
snapshot; it is substrate-restricted to ``lattice_1d`` via the orchestration
registry AND to integer velocity arrays via a content-level prerequisite.

DISCRIMINATION FROM NEIGHBORS:

P14 (SOC / self-organized criticality): NS jams have a long-but-not-
  power-law lifetime distribution (the tail truncates at finite ring size),
  while BTW sandpile avalanche distributions show clean power laws over
  decades. P14 operates on 2D substrate (avalanche sizes as metadata)
  and does not overlap substrate-wise. Substrate-level rejection is
  clean.

P4 (territoriality): NS has no territorial mechanism: no "home range"
  or resource-claiming per car. Cars move uniformly forward. P4 is not
  currently implemented in EPC, but even if introduced, a velocities-
  only primary cannot be confounded.

Density saturation (NOT a pattern): at very high density (rho >
  1/(v_max+1)) even a noise-free NS has high stopped-fraction simply
  from pigeonhole: there isn't room for every car to move. This is a
  TRIVIAL consequence of density, not emergent jamming. The confirmation
  gate `jam_lifetime_p95 > 5` cleanly rejects this case: deterministic
  pigeonhole produces only short uniform stops (p95 <= 4 at rho=0.80,
  p=0), while true NS jamming produces heavy-tailed lifetimes.

MECHANISM:

NS jamming arises from the interaction of the randomization rule (cars
occasionally over-slow) with the slowing-down rule (cars must brake for
the car ahead). When a car slows unnecessarily, the car behind has to
brake, propagating a shockwave upstream. At low density the wave
dissipates; above a critical density rho_c ~ 0.1 - 0.15 (at v_max = 5,
p = 0.3), shockwaves persist and spontaneous jams appear. See Bette
et al. (2017) for the full jamming-rate x lifetime x size decomposition.

DETECTOR TIERS:

Prerequisites (all required for any detection):
  1. 'velocities' observable present as a 1D integer array in every
     snapshot (rejects Zhang sorting, Schelling, GoL, SIR, RPS, LV,
     GH, NowakMay, GS, BTW — all of which either have no velocities
     or 2D velocity vectors).
  2. Integer-valued velocities in [0, 64] (rejects Vicsek/D'Orsogna's
     continuous 2D velocity arrays should substrate gating ever fail).
  3. n_cars >= 20 and run length >= 100 post-burn-in steps.

Primary metric: stopped_fraction = <1[v_i(t) = 0]>_{t > burn_in, i} —
  time-and-car-averaged fraction of v=0 states. (Bette et al. 2017 P(v=0)
  order parameter.)

Secondary metrics:
  - jam_lifetime_p95, jam_lifetime_max: tail statistics of per-car
    consecutive-stopped runs. Heavy-tailed in true jamming.
  - gap_cv, zero_gap_fraction: spatial inhomogeneity of gaps. Bimodal
    in jammed regimes (some gaps = 0 inside jams, large gaps between).
  - n_jam_events: total number of stopped runs.

Screening:    prereqs pass + stopped_fraction > 0.05.
Confirmation: screening + jam_lifetime_p95 > 5 + null p < 0.01 on
              jam_lifetime_p95 (temporal-shuffle null). The jam_lifetime
              gate distinguishes emergent NS jamming from deterministic
              density saturation.
Definitive:   confirmation + stopped_fraction > 0.15 + jam_lifetime_max
              > 20 + exclusions cleared.

Null model: per-car independent temporal shuffle of v(t). Preserves the
  marginal v-distribution per car (so stopped_fraction is unchanged),
  destroys the persistence that makes real jam lifetimes heavy-tailed.
  Under the null, p95 collapses to the geometric-distribution p95 at
  the observed stopped marginal (typically 2-4 steps).

EMPIRICAL CALIBRATION (Sprint 15, L=1000, v_max=5, 3-seed mean,
1000 burn-in + 2000 measurement):

  Canonical positives:
    rho=0.15, p=0.3: stopped=0.18, lt_p95=13, lt_max~60, p<0.01 -> DEFINITIVE
    rho=0.30, p=0.3: stopped=0.43, lt_p95=13, lt_max~85, p<0.01 -> DEFINITIVE

  Near-transition (CONFIRMATION / DEFINITIVE boundary):
    rho=0.12, p=0.3: stopped=0.08, lt_p95=12, lt_max~50, p<0.01 -> CONFIRMATION

  Negatives:
    rho=0.05, p=0.3:  stopped=0.00 -> rejected at screening
    rho=0.15, p=0.0:  stopped=0.00 -> rejected at screening
    rho=0.80, p=0.0:  stopped=0.75, lt_p95=4  -> rejected at confirmation
                                                  (density saturation)

References:
  Nagel, K. & Schreckenberg, M. (1992). J. Phys. I France 2, 2221.
  Bette, H. M., Habel, L., Emig, T. & Schreckenberg, M. (2017).
    Phys. Rev. E 95, 012311.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from epc.base_detector import BaseDetector
from epc.detector_result import DetectionTier, DetectorResult, NullType
from epc.metrics.traffic_jamming import (
    collect_gap_history,
    collect_velocity_history,
    gap_cv as compute_gap_cv,
    jam_lifetime_stats,
    stopped_fraction,
    temporal_shuffle_null_jam_p95,
    zero_gap_fraction,
)


class P8TrafficJammingDetector(BaseDetector):
    """Detector for P8 — emergent traffic jamming.

    Parameters
    ----------
    n_permutations : int
        Number of temporal-shuffle null permutations. Default 199 (floor
        p = 1/200 = 0.005). Confirmation requires p < 0.01 which needs
        >= 199.
    burn_in : int
        Number of initial snapshots to discard before measurement.
        Default 1000 (matches the canonical measurement protocol at
        L = 1000, v_max = 5).
    null_max_cars : int
        For the temporal-shuffle null, we subsample to at most this many
        cars to keep the permutation cost bounded. Default 80. The null
        is a local (per-car) shuffle so the p95 statistic is extensively
        insensitive to car count above ~ 50.
    seed : int
        RNG seed for the null permutations.
    """

    # Prerequisite thresholds
    _MIN_CARS = 20
    _MIN_RUN_LENGTH_POST_BURN = 100
    _MAX_V = 64  # reject continuous-valued or absurd velocity observables

    # Screening thresholds
    _SCREEN_STOPPED_MIN = 0.05

    # Confirmation thresholds (in addition to screening)
    _CONFIRM_JAM_LT_P95_MIN = 5.0
    _CONFIRM_NULL_P_MAX = 0.01

    # Definitive thresholds (in addition to confirmation)
    _DEF_STOPPED_MIN = 0.15
    _DEF_JAM_LT_MAX_MIN = 20

    def __init__(
        self,
        n_permutations: int = 199,
        burn_in: int = 1000,
        null_max_cars: int = 80,
        seed: int = 42,
    ) -> None:
        super().__init__(
            pattern_id="P8",
            # No P* neighbors currently share substrate with P8.
            # Leaving the list empty rather than listing P14 (SOC) avoids
            # creating a spurious exclusion check — P14 is on 2D substrate
            # and cannot co-occur by registration.
            excluded_patterns=[],
            allowed_co_occurrences=[],
            observable_scope="state_history_only",
        )
        if n_permutations < 1:
            raise ValueError("n_permutations must be >= 1")
        if burn_in < 0:
            raise ValueError("burn_in must be >= 0")
        if null_max_cars < 5:
            raise ValueError("null_max_cars must be >= 5")
        self.n_permutations = n_permutations
        self.burn_in = burn_in
        self.null_max_cars = null_max_cars
        self._seed = seed

        # Caches populated during detect()
        self._velocity_history: np.ndarray | None = None
        self._gap_history: np.ndarray | None = None
        self._screening_rejection_reason: str = "none"

    # ------------------------------------------------------------------
    # Timescale / prerequisites
    # ------------------------------------------------------------------

    def _estimate_timescale(
        self,
        state_history: list[dict[str, Any]],
        model_metadata: dict[str, Any] | None,
    ) -> float:
        """System timescale: L / v_max (one ring traversal).

        Falls back to run_length / 10 if metadata is absent, matching
        the base-class default.
        """
        if model_metadata is not None:
            L = model_metadata.get("L")
            v_max = model_metadata.get("v_max")
            if L is not None and v_max is not None and v_max > 0:
                return float(L) / float(v_max)
        if state_history:
            last = state_history[-1]
            L = last.get("L")
            v_max = last.get("v_max")
            if L is not None and v_max is not None and v_max > 0:
                return float(L) / float(v_max)
        return max(1.0, len(state_history) / 10.0)

    def _validate_prerequisites(
        self,
        state_history: list[dict[str, Any]],
        timescale: float,
    ) -> list[str]:
        """Populate detector caches while checking substrate prerequisites."""
        warnings: list[str] = []
        self._screening_rejection_reason = "none"

        if not state_history:
            warnings.append("empty state_history")
            self._screening_rejection_reason = "empty_state_history"
            return warnings

        last = state_history[-1]
        if "velocities" not in last:
            warnings.append(
                "no 'velocities' observable in state_history — P8 requires "
                "an integer 1D velocity array per car per snapshot"
            )
            self._screening_rejection_reason = "substrate_mismatch"
            return warnings

        # Collect post-burn-in velocity history
        burn = min(self.burn_in, max(0, len(state_history) - 1))
        vh = collect_velocity_history(state_history, burn_in=burn)
        if vh is None:
            warnings.append(
                "velocities observable is non-1D or shape-inconsistent across "
                "snapshots — P8 requires per-car integer 1D velocity arrays"
            )
            self._screening_rejection_reason = "substrate_mismatch"
            return warnings

        T, N = vh.shape
        if N < self._MIN_CARS:
            warnings.append(
                f"n_cars = {N} < minimum {self._MIN_CARS}; insufficient for "
                f"stable P8 statistics"
            )
            self._screening_rejection_reason = "too_few_cars"
            return warnings
        if T < self._MIN_RUN_LENGTH_POST_BURN:
            warnings.append(
                f"post-burn-in run length = {T} < minimum "
                f"{self._MIN_RUN_LENGTH_POST_BURN}; insufficient for stable "
                f"P8 statistics"
            )
            self._screening_rejection_reason = "run_too_short"
            return warnings

        # Integer-typed velocities with reasonable range
        if not np.issubdtype(vh.dtype, np.integer):
            warnings.append(
                f"velocities dtype is {vh.dtype} (non-integer) — P8 requires "
                f"integer velocities in [0, v_max]"
            )
            self._screening_rejection_reason = "non_integer_velocities"
            return warnings
        if vh.min() < 0 or vh.max() > self._MAX_V:
            warnings.append(
                f"velocity range [{vh.min()}, {vh.max()}] outside "
                f"[0, {self._MAX_V}] — P8 is calibrated for traffic-CA-style "
                f"small integer velocities"
            )
            self._screening_rejection_reason = "velocity_range_out_of_bounds"
            return warnings

        self._velocity_history = vh
        # Gap history is optional (secondary metrics); fall back to None.
        self._gap_history = collect_gap_history(state_history, burn_in=burn)

        # Run-length sanity vs timescale
        recommended_min = int(5 * timescale)
        if T < recommended_min:
            warnings.append(
                f"post-burn-in run length {T} < recommended {recommended_min} "
                f"(5 * timescale {timescale:.0f}) — primary metric may be "
                f"under-converged"
            )

        return warnings

    # ------------------------------------------------------------------
    # Primary / secondaries / screening
    # ------------------------------------------------------------------

    def _compute_primary(
        self,
        state_history: list[dict[str, Any]],
        timescale: float,
    ) -> dict[str, float]:
        """Primary metric: stopped_fraction. Also reports key diagnostic
        fields for downstream tier logic.
        """
        if self._velocity_history is None:
            return {
                "stopped_fraction": 0.0,
                "mean_velocity": 0.0,
                "n_cars": 0,
                "n_post_burn_steps": 0,
                "screening_rejection_reason": self._screening_rejection_reason,
            }
        vh = self._velocity_history
        T, N = vh.shape
        return {
            "stopped_fraction": stopped_fraction(vh),
            "mean_velocity": float(vh.mean()),
            "n_cars": int(N),
            "n_post_burn_steps": int(T),
            "screening_rejection_reason": self._screening_rejection_reason,
        }

    def _check_screening(
        self,
        primary_result: dict[str, float],
        timescale: float,
    ) -> bool:
        """Screening gate: substrate prereqs PASSED and stopped_fraction > 0.05.

        Note: when screening fails specifically on the stopped-fraction
        floor, we mutate ``primary_result['screening_rejection_reason']``
        in-place so that the final DetectorResult surfaces the specific
        failure mode (below_stopped_floor) rather than leaving the
        prereq-level 'none'.
        """
        if self._screening_rejection_reason != "none":
            # Substrate-level prereq already failed; leave reason alone.
            return False
        sf = primary_result.get("stopped_fraction", 0.0)
        if sf <= self._SCREEN_STOPPED_MIN:
            self._screening_rejection_reason = "below_stopped_floor"
            # Propagate to the primary_result dict that DetectorResult
            # will carry forward (mutated in-place; safe because
            # primary_result is a local dict on this turn only).
            primary_result["screening_rejection_reason"] = "below_stopped_floor"
            return False
        return True

    def _compute_secondaries(
        self,
        state_history: list[dict[str, Any]],
        timescale: float,
    ) -> dict[str, Any]:
        """Secondaries: jam-lifetime distribution statistics + gap stats."""
        if self._velocity_history is None:
            return {}
        vh = self._velocity_history
        jl = jam_lifetime_stats(vh)
        out: dict[str, Any] = {
            "jam_lifetime_mean": float(jl["mean"]),
            "jam_lifetime_median": float(jl["median"]),
            "jam_lifetime_p95": float(jl["p95"]),
            "jam_lifetime_max": int(jl["max"]),
            "n_jam_events": int(jl["n_jam_events"]),
        }
        if self._gap_history is not None:
            out["gap_cv"] = float(compute_gap_cv(self._gap_history))
            out["zero_gap_fraction"] = float(zero_gap_fraction(self._gap_history))
        return out

    # ------------------------------------------------------------------
    # Null / effect size / tier / exclusions
    # ------------------------------------------------------------------

    def _run_null_model(
        self,
        state_history: list[dict[str, Any]],
        primary_result: dict[str, float],
        timescale: float,
    ) -> tuple[float, NullType, dict[str, float]]:
        """Null: per-car independent temporal shuffle on jam_lifetime_p95.

        This is the detector's confirmation-grade null. An observed
        jam_lt_p95 >> null-distribution p95 is the signature of genuine
        NS jamming (persistence) as opposed to independent per-car
        fluctuations. Right-tailed test.
        """
        if self._velocity_history is None:
            return 1.0, NullType.SHUFFLE, {"mean": 0.0, "std": 0.0}

        vh = self._velocity_history
        T, N = vh.shape

        # Subsample cars for the null to keep cost bounded. Per-car shuffle
        # is extensively local: p95 converges quickly in N.
        if N > self.null_max_cars:
            rng = np.random.default_rng(self._seed)
            car_idx = rng.choice(N, size=self.null_max_cars, replace=False)
            vh_sub = vh[:, car_idx]
        else:
            vh_sub = vh

        observed_p95 = float(jam_lifetime_stats(vh_sub)["p95"])
        nulls = temporal_shuffle_null_jam_p95(
            vh_sub,
            n_permutations=self.n_permutations,
            seed=self._seed,
        )
        if observed_p95 <= 0.0:
            return 1.0, NullType.SHUFFLE, {
                "mean": float(nulls.mean()),
                "std": float(nulls.std()),
                "observed": 0.0,
            }
        p = (float((nulls >= observed_p95).sum()) + 1.0) / (len(nulls) + 1.0)
        return (
            float(p),
            NullType.SHUFFLE,
            {
                "mean": float(nulls.mean()),
                "std": float(nulls.std()),
                "observed": observed_p95,
            },
        )

    def _compute_effect_size(
        self,
        primary_result: dict[str, float],
        null_dist_stats: dict[str, float],
    ) -> dict[str, float]:
        """Cohen's d on jam_lifetime_p95 vs temporal-shuffle null.

        Observed jam_lt_p95 is compared against the null distribution's
        mean and std. Positive d means observed is HIGHER than null —
        the direction indicating real jamming persistence.

        Note: null std can be exactly zero at canonical parameters (all
        shuffles give p95 = 2 because the shuffled marginal is identical
        across perms). We fall back to a small epsilon so d is reported
        as a large finite number rather than NaN / inf — this matches
        the Gray-Scott+P3 precedent.
        """
        observed = float(null_dist_stats.get("observed", 0.0))
        nm = float(null_dist_stats.get("mean", 0.0))
        ns = float(null_dist_stats.get("std", 0.0))
        if ns <= 1e-9:
            # Effectively infinite separation; report as large finite number
            # for numerical sanity.
            if observed > nm:
                d = 1e6
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

    def _check_confirmation(
        self,
        primary_result: dict[str, float],
        secondary_result: dict[str, Any],
        null_p: float,
        timescale: float,
    ) -> bool:
        """Confirmation: screening + jam_lifetime_p95 > 5 + null p < 0.01."""
        if not secondary_result:
            return False
        p95 = float(secondary_result.get("jam_lifetime_p95", 0.0))
        if p95 <= self._CONFIRM_JAM_LT_P95_MIN:
            return False
        if null_p >= self._CONFIRM_NULL_P_MAX:
            return False
        return True

    def _check_definitive(
        self,
        primary_result: dict[str, float],
        secondary_result: dict[str, Any],
        null_p: float,
        null_type: NullType,
        state_history: list[dict[str, Any]],
        model_metadata: dict[str, Any] | None,
        timescale: float,
    ) -> bool:
        """Definitive: confirmation + stopped_fraction > 0.15 +
        jam_lifetime_max > 20.

        No mechanistic-intervention null is required for P8: the
        temporal-shuffle null already delivers large effect sizes
        (Cohen's d > 1000 at canonical parameters due to null std ~ 0).
        This matches the Gray-Scott+P3 precedent (Decision 37/38).
        """
        sf = float(primary_result.get("stopped_fraction", 0.0))
        if sf <= self._DEF_STOPPED_MIN:
            return False
        lt_max = int(secondary_result.get("jam_lifetime_max", 0))
        if lt_max <= self._DEF_JAM_LT_MAX_MIN:
            return False
        return True

    def _all_secondaries_pass(self, secondary_result: dict[str, Any]) -> bool:
        """All P8 secondaries pass = jam_lifetime_p95 > 5."""
        if not secondary_result:
            return False
        return float(secondary_result.get("jam_lifetime_p95", 0.0)) > self._CONFIRM_JAM_LT_P95_MIN

    def _check_exclusions(
        self,
        state_history: list[dict[str, Any]],
        model_metadata: dict[str, Any] | None,
        timescale: float,
    ) -> tuple[list[str], dict[str, str]]:
        """No currently-implemented neighbor patterns share substrate with P8.

        P8's substrate prerequisite (1D integer velocities) cleanly excludes
        all other EPC patterns by construction. Future pattern additions
        (e.g., a territorial P4 on 1D occupancy, or a wave P13 variant on
        1D excitable media) would trigger an update here.
        """
        return [], {}
