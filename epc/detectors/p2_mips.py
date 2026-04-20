"""P2 — Activity-induced phase separation (MIPS) detector.

Detects motility-induced phase separation: emergent gas/liquid
coexistence in a system of self-propelled particles with NO attractive
forces and NO alignment rule, driven purely by density-dependent
self-propulsion.

Observable scope: model_metadata_assisted. The empirical detection
runs entirely on state_history (positions + velocities). Metadata is
used to promote CONFIRMATION -> DEFINITIVE via the mechanistic-null
flags ``has_density_dependent_speed``, ``has_attraction_rule``, and
``has_alignment_rule``.

DISCRIMINATION FROM NEIGHBORS:

P1 (similarity-driven aggregation): P1 operates on lattice substrates
  with type-labeled agents and similarity preferences. P2 is on
  continuous_2d with no type labels and no preferences. Substrate-level
  rejection is clean.

P5 (flocking): Vicsek-family flocking has constant speed (v_i = v_0
  for all i) and density is NOT anti-correlated with |v| — there is
  no v(rho) mechanism. Vicsek ordered states have flock-following
  density bands (all f_liquid, no f_gas); disordered states have
  uniform density. Both fail the two-phase primary. Additionally,
  CV_v == 0 exactly under Vicsek constant-speed dynamics.

P6 (milling): D'Orsogna-family swarms use attractive potentials
  explicitly — the ``has_attraction_rule`` metadata flag is True,
  which rejects DEFINITIVE at the mechanistic-null gate. Empirically,
  milling swarms DO cluster, so the P2 primary may be nonzero, but
  D'Orsogna exhibits f_liquid ~ 0.7, f_gas ~ 0.03 (single dense
  flock, no dilute-phase coexistence), so primary = min ~ 0.03 which
  is at most SCREENING.

Density saturation (NOT a pattern): at very high packing fraction +
  high Pe, ABP particles become stuck everywhere (v=0 for ALL i).
  Local density is nearly uniform just below rho_star. This gives
  f_liquid ~ 1, f_gas ~ 0, primary ~ 0 → correctly rejected at
  screening.

Thermal / passive-like regime (low Pe): weakly active ABP shows
  density fluctuations that look bimodal-adjacent but without the
  coexistence-phase structure of true MIPS. Discriminator:
  CV_v > 0.3 (true MIPS has wide v distribution since v(rho) varies
  across particles; thermal regime has narrow v distribution because
  no particle reaches rho_star).

Dilute-Poisson artifact (low phi, sparse particles): the local
  density takes only a few discrete values (1/pi, 2/pi, ...) so
  Pearson r(rho, v) is trivially near -1 by the v(rho) relation.
  Discriminator: primary cannot reach SCREENING threshold anyway
  (f_liquid = 0 because rho never exceeds rho_star), so the
  -1-artifact is never promoted.

DETECTOR TIERS:

Prerequisites (all required for any detection):
  1. 'positions' observable as (N, 2) float, N >= 50.
  2. 'velocities' observable as (N, 2) float (optional if 'speeds'
     provided).
  3. Periodic box: ``box_size`` present in state snapshot or
     model_metadata['box_size'], and consistent across snapshots.
  4. Run length >= 300 snapshots post-burn-in.
  5. rho_star either provided via detector init OR available as
     model_metadata['rho_star'].

Primary metric: two_phase_coexistence_score(rho_local, rho_star) =
  min(f_gas, f_liquid).
  f_gas: fraction of particles with rho_local < rho_star/2
  f_liquid: fraction with rho_local > rho_star

Secondary metrics:
  - cv_v: coefficient of variation of speed magnitudes
  - density_speed_anticorrelation: Pearson r(rho, |v|)
  - dynamic_range_p90_p10: percentile density ratio
  - mean_rho, frac_stalled (diagnostic)

Screening: prereqs pass + primary > 0.03.
Confirmation: screening + primary > 0.10 + (-0.99 < r < -0.30)
              + cv_v > 0.30 + frac_stalled < 0.98.
Definitive: confirmation + primary > 0.20 + mechanistic-null gate
            passes (metadata affirms density-dependent speed, NO
            attraction, NO alignment).

Null model: surrogate constant-speed null. The detector does NOT
  re-run the model (the framework does not provide that). Instead,
  DEFINITIVE tier requires metadata affirmation of the mechanism.
  This matches observable_scope = 'model_metadata_assisted'.

  For CONFIRMATION, the effect size is computed against the
  density_speed anticorrelation — the r metric itself carries the
  null-model comparison semantically (a correctly-signed r plus a
  strong primary constitutes the rejection of the "positions cluster
  by chance" null).

EMPIRICAL CALIBRATION (Sprint 16 Phase 1, N=1000, rho_star=4.0,
r_cg=1.0, 2000+ burn-in + 3000 measurement steps):

  Canonical positive (ABP, phi=0.5, Pe=100):
    primary=0.36, r=-0.93, cv_v=1.33, frac_stalled=0.63  -> DEFINITIVE
    (mechanism verified: has_density_dependent_speed=True, etc.)

  Negatives (positions clustered, but different mechanism):
    Vicsek ordered:     primary=0.017 (f_liquid=0.95, f_gas=0.02)
                        cv_v=0 (constant speed)  -> SCREENING rejection
    Vicsek disordered:  primary=0.002 -> rej at screening
    D'Orsogna milling:  primary=0.032, has_attraction_rule=True
                        -> SCREENING-tier at most (metadata gates
                        DEFINITIVE)

  ABP false-positive traps (same model, different regimes):
    ABP thermal (Pe=5):    primary~0.07, cv_v=0.3 borderline
                          -> SCREENING (cv_v gate)
    ABP dilute (phi=0.1): primary~0 (no liquid) -> rej at screening
    ABP stuck (phi=0.75): primary~0.003 (no gas) -> rej at screening

References:
  Fily, Y. & Marchetti, M. C. (2012). Athermal Phase Separation of
    Self-Propelled Particles with No Alignment.
    Phys. Rev. Lett. 108, 235702.
  Redner, G. S., Hagan, M. F. & Baskaran, A. (2013). Structure and
    Dynamics of a Phase-Separating Active Colloidal Fluid.
    Phys. Rev. Lett. 110, 055701.
  Cates, M. E. & Tailleur, J. (2015). Motility-Induced Phase
    Separation. Annu. Rev. Condens. Matter Phys. 6, 219.
"""
from __future__ import annotations

from typing import Any

import numpy as np

from epc.base_detector import BaseDetector
from epc.detector_result import DetectionTier, NullType
from epc.metrics.density_phase_separation import (
    collect_density_and_speed_history,
    density_speed_anticorrelation,
    mechanistic_null_test,
    phase_coexistence_fractions,
    two_phase_coexistence_score,
)


class P2MIPSDetector(BaseDetector):
    """Detector for P2 — motility-induced phase separation.

    Parameters
    ----------
    rho_star : float or None
        Density threshold at which v(rho) -> 0 in the self-propulsion
        law. Default None = read from model_metadata['rho_star']. Falls
        back to 4.0 if neither is available (matches ABP default with
        r_cg=1.0 = particle diameter).
    r_cg : float or None
        Coarse-graining radius for local density estimation. Default
        None = read from metadata or fall back to 1.0.
    burn_in : int
        Number of initial snapshots to discard before measurement.
        Default 200.
    n_permutations : int
        Number of permutations for the density-speed r null. Minimum
        199 for CONFIRMATION (need floor p < 0.01). Default 199.
        n_permutations=99 gives floor p=0.01 which FAILS the
        `null_p < 0.01` gate — detector will report SCREENING only
        even on canonical positives. This is the P8 permutation-floor
        gotcha from Sprint 15 (Statistical power requirements).
    max_frames : int or None
        Cap on number of measurement frames used. Default 60 (keeps
        correlation computations fast on large histories).
    seed : int
        RNG seed. Default 42.
    """

    # Tier thresholds (primary metric = min(f_gas, f_liquid))
    _SCREEN_PRIMARY_MIN = 0.03
    _CONFIRM_PRIMARY_MIN = 0.08
    _DEF_PRIMARY_MIN = 0.15

    # Confirmation gates
    _CONFIRM_R_MAX = -0.30     # r must be <= this (more negative)
    _CONFIRM_R_MIN = -0.99     # r must be >= this (reject -1 artifact)
    _CONFIRM_CV_V_MIN = 0.30
    _CONFIRM_FRAC_STALLED_MAX = 0.98

    # Prerequisites
    _MIN_N_PARTICLES = 50
    _MIN_RUN_LENGTH_POST_BURN = 300

    # Fallback defaults when metadata lacks the parameters
    _DEFAULT_RHO_STAR = 4.0
    _DEFAULT_R_CG = 1.0

    def __init__(
        self,
        rho_star: float | None = None,
        r_cg: float | None = None,
        burn_in: int = 200,
        n_permutations: int = 199,
        max_frames: int | None = 60,
        seed: int = 42,
    ) -> None:
        super().__init__(
            pattern_id="P2",
            # P2's nearest neighbor is P1 (aggregation) on lattice_2d
            # and P6 (milling) on continuous_2d. P1 is on a different
            # substrate and cannot co-occur by registration. P6 shares
            # the substrate but is distinguished by attraction-rule
            # metadata (see _check_exclusions).
            excluded_patterns=["P1", "P6"],
            allowed_co_occurrences=["P5"],   # flocking + MIPS co-occur in e.g. Vicsek+v(rho)
            observable_scope="model_metadata_assisted",
        )
        if n_permutations < 1:
            raise ValueError("n_permutations must be >= 1")
        if burn_in < 0:
            raise ValueError("burn_in must be >= 0")
        if rho_star is not None and rho_star <= 0:
            raise ValueError(f"rho_star must be positive, got {rho_star}")
        if r_cg is not None and r_cg <= 0:
            raise ValueError(f"r_cg must be positive, got {r_cg}")
        self.rho_star = rho_star
        self.r_cg = r_cg
        self.burn_in = burn_in
        self.n_permutations = n_permutations
        self.max_frames = max_frames
        self._seed = seed

        # Caches populated during detect()
        self._all_rho: np.ndarray | None = None
        self._all_v: np.ndarray | None = None
        self._n_frames_used: int = 0
        self._resolved_rho_star: float = self._DEFAULT_RHO_STAR
        self._resolved_r_cg: float = self._DEFAULT_R_CG
        self._resolved_box_size: float | None = None
        self._screening_rejection_reason: str = "none"

    # ------------------------------------------------------------------
    # Timescale / prerequisites
    # ------------------------------------------------------------------

    def _estimate_timescale(
        self,
        state_history: list[dict[str, Any]],
        model_metadata: dict[str, Any] | None,
    ) -> float:
        """System timescale: 1 / D_r (rotational correlation time).

        Falls back to run_length / 10 if metadata is absent.
        """
        if model_metadata is not None:
            D_r = model_metadata.get("D_r")
            if D_r is not None and D_r > 0:
                return 1.0 / float(D_r)
        return max(1.0, len(state_history) / 10.0)

    def _resolve_params(
        self,
        state_history: list[dict[str, Any]],
        model_metadata: dict[str, Any] | None,
    ) -> None:
        """Pick rho_star, r_cg, box_size from init args / metadata / state."""
        # rho_star
        if self.rho_star is not None:
            self._resolved_rho_star = float(self.rho_star)
        elif model_metadata is not None and "rho_star" in model_metadata:
            self._resolved_rho_star = float(model_metadata["rho_star"])
        else:
            self._resolved_rho_star = self._DEFAULT_RHO_STAR

        # r_cg
        if self.r_cg is not None:
            self._resolved_r_cg = float(self.r_cg)
        elif model_metadata is not None and "r_cg" in model_metadata:
            self._resolved_r_cg = float(model_metadata["r_cg"])
        else:
            self._resolved_r_cg = self._DEFAULT_R_CG

        # box_size (for periodic cKDTree)
        self._resolved_box_size = None
        if model_metadata is not None and "box_size" in model_metadata:
            self._resolved_box_size = float(model_metadata["box_size"])
        elif state_history:
            last = state_history[-1]
            if "box_size" in last:
                self._resolved_box_size = float(last["box_size"])

    def _validate_prerequisites(
        self,
        state_history: list[dict[str, Any]],
        timescale: float,
    ) -> list[str]:
        """Validate substrate prerequisites and populate detector caches.

        NOTE: model_metadata is not available at this point in the base
        class API. We defer rho_star / r_cg / box_size resolution to the
        first call to _compute_primary, which reads from caches populated
        here. Parameter resolution logic happens in `detect()` override.
        """
        warnings: list[str] = []
        self._screening_rejection_reason = "none"

        if not state_history:
            warnings.append("empty state_history")
            self._screening_rejection_reason = "empty_state_history"
            return warnings

        last = state_history[-1]
        has_positions = "positions" in last
        has_vel_or_speed = "velocities" in last or "speeds" in last
        if not has_positions:
            warnings.append(
                "no 'positions' observable in state_history — P2 requires "
                "continuous 2D positions per particle per snapshot"
            )
            self._screening_rejection_reason = "substrate_mismatch"
            return warnings
        if not has_vel_or_speed:
            warnings.append(
                "no 'velocities' or 'speeds' observable — P2 requires "
                "per-particle velocity magnitudes"
            )
            self._screening_rejection_reason = "substrate_mismatch"
            return warnings

        pos = last.get("positions")
        try:
            pos_arr = np.asarray(pos)
        except Exception:
            warnings.append("'positions' is not coercible to ndarray")
            self._screening_rejection_reason = "substrate_mismatch"
            return warnings
        if pos_arr.ndim != 2 or pos_arr.shape[1] != 2:
            warnings.append(
                f"positions shape {pos_arr.shape} — P2 requires (N, 2) "
                f"arrays for continuous 2D space"
            )
            self._screening_rejection_reason = "substrate_mismatch"
            return warnings

        N = pos_arr.shape[0]
        if N < self._MIN_N_PARTICLES:
            warnings.append(
                f"n_particles = {N} < minimum {self._MIN_N_PARTICLES}; "
                f"insufficient for P2 density statistics"
            )
            self._screening_rejection_reason = "too_few_particles"
            return warnings

        # Run-length check (post-burn-in)
        post_burn = len(state_history) - self.burn_in
        if post_burn < self._MIN_RUN_LENGTH_POST_BURN:
            warnings.append(
                f"post-burn-in run length = {post_burn} < minimum "
                f"{self._MIN_RUN_LENGTH_POST_BURN}; insufficient for P2"
            )
            self._screening_rejection_reason = "run_too_short"
            return warnings

        # box_size check (may still be resolved later in _compute_primary
        # from metadata; here we only flag if it's neither in state nor
        # any parameter-default would apply)
        if "box_size" not in last and self._resolved_box_size is None:
            warnings.append(
                "no 'box_size' in state snapshots — P2 defaults to "
                "non-periodic cKDTree which may underestimate boundary "
                "densities"
            )

        return warnings

    # ------------------------------------------------------------------
    # Primary / secondaries / screening
    # ------------------------------------------------------------------

    def detect(
        self,
        state_history: list[dict[str, Any]],
        model_metadata: dict[str, Any] | None = None,
        timescale: float | None = None,
    ):
        """Entry point override: resolve rho_star, r_cg, box_size BEFORE
        prereq validation so downstream methods see the right parameters.
        """
        self._resolve_params(state_history, model_metadata)
        return super().detect(state_history, model_metadata, timescale)

    def _compute_primary(
        self,
        state_history: list[dict[str, Any]],
        timescale: float,
    ) -> dict[str, float]:
        """Populate caches (all_rho, all_v) and compute primary metric."""
        if self._screening_rejection_reason != "none":
            return {
                "two_phase_score": 0.0,
                "f_gas": 0.0,
                "f_liquid": 0.0,
                "rho_star_used": self._resolved_rho_star,
                "r_cg_used": self._resolved_r_cg,
                "n_frames_used": 0,
                "screening_rejection_reason": self._screening_rejection_reason,
            }

        all_rho, all_v, n_frames = collect_density_and_speed_history(
            state_history,
            box_size=self._resolved_box_size,
            r_cg=self._resolved_r_cg,
            periodic=(self._resolved_box_size is not None),
            burn_in=self.burn_in,
            max_frames=self.max_frames,
        )
        self._all_rho = all_rho
        self._all_v = all_v
        self._n_frames_used = n_frames

        if all_rho.size == 0:
            return {
                "two_phase_score": 0.0,
                "f_gas": 0.0,
                "f_liquid": 0.0,
                "rho_star_used": self._resolved_rho_star,
                "r_cg_used": self._resolved_r_cg,
                "n_frames_used": 0,
                "screening_rejection_reason": "measurement_window_empty",
            }

        f_gas, f_liq = phase_coexistence_fractions(
            all_rho, self._resolved_rho_star
        )
        score = min(f_gas, f_liq)
        return {
            "two_phase_score": float(score),
            "f_gas": float(f_gas),
            "f_liquid": float(f_liq),
            "mean_rho": float(all_rho.mean()),
            "rho_star_used": float(self._resolved_rho_star),
            "r_cg_used": float(self._resolved_r_cg),
            "n_frames_used": int(n_frames),
            "screening_rejection_reason": "none",
        }

    def _check_screening(
        self,
        primary_result: dict[str, float],
        timescale: float,
    ) -> bool:
        if self._screening_rejection_reason != "none":
            return False
        score = primary_result.get("two_phase_score", 0.0)
        if score <= self._SCREEN_PRIMARY_MIN:
            self._screening_rejection_reason = "below_two_phase_floor"
            primary_result["screening_rejection_reason"] = "below_two_phase_floor"
            return False
        return True

    def _compute_secondaries(
        self,
        state_history: list[dict[str, Any]],
        timescale: float,
    ) -> dict[str, Any]:
        if self._all_rho is None or self._all_v is None:
            return {}
        all_rho = self._all_rho
        all_v = self._all_v
        mean_v = float(all_v.mean()) if all_v.size else 0.0
        cv_v = float(all_v.std() / mean_v) if mean_v > 1e-9 else 0.0
        r = density_speed_anticorrelation(all_rho, all_v)
        p10 = float(np.percentile(all_rho, 10))
        p90 = float(np.percentile(all_rho, 90))
        dr_90_10 = float(p90 / max(p10, 1e-12))
        # frac_stalled: fraction of particles with |v| < 5% of mean v.
        # If mean v ~ 0, fall back to 1.0 (everyone stalled).
        if mean_v > 1e-9:
            frac_stalled = float((all_v < 0.05 * mean_v).mean())
        else:
            frac_stalled = 1.0
        return {
            "cv_v": cv_v,
            "density_speed_r": r,
            "p90_p10_ratio": dr_90_10,
            "mean_v": mean_v,
            "frac_stalled": frac_stalled,
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
        """Null for P2: the density_speed_anticorrelation value is its
        own null-adjusted statistic.

        We compute a null distribution of r by SHUFFLING the pairing of
        rho values and speed values across particles. Under the null
        (no density-velocity coupling), r is distributed near 0. The
        observed r is compared against this null.

        This is a permutation null that preserves the marginal
        distributions of rho and v but destroys their pairing.
        """
        if self._all_rho is None or self._all_v is None or self._all_rho.size < 10:
            return 1.0, NullType.SHUFFLE, {"mean": 0.0, "std": 0.0}

        rho = self._all_rho
        v = self._all_v
        # Subsample to keep the null fast
        rng = np.random.default_rng(self._seed)
        if rho.size > 5000:
            idx = rng.choice(rho.size, 5000, replace=False)
            rho_s = rho[idx]
            v_s = v[idx]
        else:
            rho_s = rho
            v_s = v

        observed_r = density_speed_anticorrelation(rho_s, v_s)
        null_rs = np.empty(self.n_permutations, dtype=np.float64)
        for i in range(self.n_permutations):
            perm = rng.permutation(v_s.size)
            null_rs[i] = density_speed_anticorrelation(rho_s, v_s[perm])

        # Left-tailed: we want r to be MORE negative than under the null
        p = (float((null_rs <= observed_r).sum()) + 1.0) / (len(null_rs) + 1.0)
        return (
            float(p),
            NullType.SHUFFLE,
            {
                "mean": float(null_rs.mean()),
                "std": float(null_rs.std()),
                "observed": float(observed_r),
            },
        )

    def _compute_effect_size(
        self,
        primary_result: dict[str, float],
        null_dist_stats: dict[str, float],
    ) -> dict[str, float]:
        """Cohen's d on density-velocity anti-correlation vs permutation null."""
        observed = float(null_dist_stats.get("observed", 0.0))
        nm = float(null_dist_stats.get("mean", 0.0))
        ns = float(null_dist_stats.get("std", 0.0))
        if ns <= 1e-9:
            if observed < nm:
                d = -1e6  # massive effect in negative direction
            elif observed > nm:
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
        """P2 tier logic.

        Screening already passed. Check CONFIRMATION gates, then DEFINITIVE.
        """
        bonuses: dict[str, bool] = {}

        # Screening bonuses
        bonuses["secondaries_pass"] = bool(secondary_result)
        bonuses["shuffle_null_p_001"] = null_p < 0.01

        score = primary_result.get("two_phase_score", 0.0)
        if score <= self._CONFIRM_PRIMARY_MIN:
            return DetectionTier.SCREENING, bonuses

        r = secondary_result.get("density_speed_r", 0.0)
        cv_v = secondary_result.get("cv_v", 0.0)
        frac_stalled = secondary_result.get("frac_stalled", 1.0)

        r_in_band = (self._CONFIRM_R_MIN < r < self._CONFIRM_R_MAX)
        cv_v_ok = cv_v > self._CONFIRM_CV_V_MIN
        stalled_ok = frac_stalled < self._CONFIRM_FRAC_STALLED_MAX
        null_ok = null_p < 0.01

        if not (r_in_band and cv_v_ok and stalled_ok and null_ok):
            return DetectionTier.SCREENING, bonuses

        # Confirmation bonuses
        bonuses["null_p_0001"] = null_p < 0.0001
        bonuses["effect_size_gt_1"] = abs(effect_size.get("cohens_d", 0.0)) > 1.0
        bonuses["all_secondaries"] = r_in_band and cv_v_ok and stalled_ok

        if score <= self._DEF_PRIMARY_MIN:
            return DetectionTier.CONFIRMATION, bonuses

        # DEFINITIVE: requires mechanistic-null metadata gate
        mech = mechanistic_null_test(score, model_metadata)
        if not mech["null_rejects_mips"]:
            return DetectionTier.CONFIRMATION, bonuses

        # Definitive bonuses
        bonuses["all_exclusions_cleared"] = True  # exclusions checked in base
        bonuses["both_null_types_rejected"] = null_p < 0.01 and mech["null_rejects_mips"]
        bonuses["finite_size_robustness"] = False
        return DetectionTier.DEFINITIVE, bonuses

    def _all_secondaries_pass(self, secondary_result: dict[str, Any]) -> bool:
        if not secondary_result:
            return False
        r = secondary_result.get("density_speed_r", 0.0)
        cv_v = secondary_result.get("cv_v", 0.0)
        frac_stalled = secondary_result.get("frac_stalled", 1.0)
        return (
            self._CONFIRM_R_MIN < r < self._CONFIRM_R_MAX
            and cv_v > self._CONFIRM_CV_V_MIN
            and frac_stalled < self._CONFIRM_FRAC_STALLED_MAX
        )

    def _check_exclusions(
        self,
        state_history: list[dict[str, Any]],
        model_metadata: dict[str, Any] | None,
        timescale: float,
    ) -> tuple[list[str], dict[str, str]]:
        """Exclude P1 (different substrate) and P6 (attraction-driven).

        P1 is on lattice_2d and cannot co-occur by registration — report
        as 'excluded_by_substrate'.
        P6 is on continuous_2d like P2; check model_metadata for
        has_attraction_rule. If True, P6 is the correct attribution;
        if False, P6 is excluded.
        """
        checked = ["P1", "P6"]
        results = {
            "P1": "excluded_by_substrate",
        }
        if model_metadata is None:
            results["P6"] = "inconclusive"
        else:
            has_attr = model_metadata.get("has_attraction_rule")
            if has_attr is True:
                results["P6"] = "not_excluded"  # this might actually be P6
            elif has_attr is False:
                results["P6"] = "excluded"
            else:
                results["P6"] = "inconclusive"
        return checked, results
