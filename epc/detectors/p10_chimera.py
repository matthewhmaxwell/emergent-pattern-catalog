"""P10 — Chimera-state detector.

Detects coexistence of coherent and incoherent oscillator subpopulations
on a non-locally coupled oscillator ring (Abrams-Strogatz 2004 /
Kuramoto-Battogtokh 2002).

Observable scope: ``model_metadata_assisted``.

DISCRIMINATION FROM NEIGHBORS
------------------------------

P9 (temporal synchronization): P9 already has a P10 exclusion gate
  (``p10_excluded = local_r_cv < 0.2``) that prevents P9 from reaching
  DEFINITIVE on chimera runs. P10 is the symmetric detector that
  asserts coexistence affirmatively.

Ordinary all-to-all Kuramoto (no chimera): every superficial chimera
  signature — "gap between max and min local r", "persistent coexisting
  windows", "spatial window ordering" — gives false positives on
  ordinary Kuramoto at K ~ K_c. Partial sync on a mean-field model
  produces the same window-level appearance after the model's internal
  ω-sort: center of the ω distribution entrains (coherent), tails drift
  (incoherent). See Phase 1 tables in REPLICATION_NOTES Sprint 18.

  The clean discriminator is ``pos_vel_ac[lag=4]`` — spatial
  autocorrelation of time-averaged per-oscillator phase velocity on the
  ring. A chimera has neighbors drifting together (common coupling
  field); ordinary Kuramoto has each oscillator drifting at its own ω,
  so neighbors on the index-sorted-by-ω ring have uncorrelated long-run
  velocities. Phase 1j N=128 separation:

    Chimera (coex-passing, multi-seed β=0.18): 0.929 ± 0.007
    Kuramoto K=1.0, 6 seeds:                   0.312 ± 0.130
    Separation gap:  chim min 0.919 - Kur max 0.448 = +0.47

  DEFINITIVE is additionally gated by the ``has_nonlocal_coupling``
  metadata flag — only models that explicitly implement non-local
  coupling can ever reach DEFINITIVE, even if the observable-level gates
  all pass (they do not, but the metadata gate provides a second line
  of defense).

DETECTOR TIERS
--------------

Prerequisites (all required for any detection):
  1. ``theta`` observable (N,) phases in state snapshots.
  2. N >= 32 for windowing.
  3. ``positions`` observable present OR ring-position assumed from
     index. Ring-position from index is the default — the windowing
     uses index contiguity regardless.
  4. Post-burn run length >= 30 frames (needed for persistence
     statistics).

Screening prerequisites (hard gates, any failure -> not detected):
  A. Coexistence: n_persistent_coh >= 1 AND n_persistent_incoh >= 1.
     Rejects full sync (n_persistent_incoh = 0) and full incoherence
     (n_persistent_coh = 0).
  B. ``pos_vel_ac[lag=4] > 0.55`` (screening floor).

Primary metric: ``pos_vel_ac[lag=4]``.

Secondary metrics:
  - r_global_mean, r_global_std (Kuramoto order parameter time-average)
  - gap_timeavg (max - min of time-averaged local r)
  - persistence_corr (diagnostic only)
  - velocity_spatial_neighbor_corr (complementary time-series measure)
  - n_persistent_coh, n_persistent_incoh

Confirmation gate: primary >= 0.75 + null p < 0.01 (via position-shuffle
  surrogate, n_permutations >= 199).

Definitive gate: CONFIRMATION + ``has_nonlocal_coupling = True`` in
  metadata. Absence of the flag or True-valued
  ``has_frequency_heterogeneity`` (i.e., ordinary Kuramoto) blocks
  DEFINITIVE and caps at CONFIRMATION.

NULL TYPE: SURROGATE (position-shuffle).

EMPIRICAL CALIBRATION (Sprint 18 Phase 1j, N=128, A=0.995, β=0.05,
dt=0.025, record_dt=1.0, 50 frames, burn 30%):

  Canonical positive (seed=0, asymmetric_gaussian IC):
    primary (pos_vel_ac[4]) = 0.844
    null_p = 0.005 (floor at n_perm=199)
    cohens_d > 9
    coex = True, n_coh = 2, n_incoh = 4
    has_nonlocal_coupling = True
    -> DEFINITIVE

  Same model, β = 0.10 (chimera basin escapes to sync):
    coex = False -> screening rejection (full sync)

  Ordinary Kuramoto (K = 1.0, heterogeneous ω, hardest negative):
    primary <= 0.45 -> screening rejection

References
----------
Abrams, D. M. & Strogatz, S. H. (2004). Chimera states for coupled
    oscillators. Phys. Rev. Lett. 93, 174102.
Kuramoto, Y. & Battogtokh, D. (2002). Coexistence of coherence and
    incoherence in non-locally coupled phase oscillators. Nonlinear
    Phenomena in Complex Systems 5, 380–385.
Panaggio, M. J. & Abrams, D. M. (2015). Chimera states: coexistence of
    coherence and incoherence in networks of coupled oscillators.
    Nonlinearity 28, R67.
"""

from __future__ import annotations

from typing import Any

import numpy as np

from epc.base_detector import BaseDetector
from epc.detector_result import DetectionTier, NullType
from epc.metrics.chimera_coexistence import (
    coexistence_stats,
    global_order_parameter,
    local_r_time_matrix,
    persistence_corr,
    phase_velocities,
    pos_vel_ac,
    shuffle_null_pos_vel_ac,
    velocity_spatial_neighbor_corr,
)


class P10ChimeraDetector(BaseDetector):
    """Detector for P10 — chimera states.

    Parameters
    ----------
    n_windows : int
        Number of contiguous ring-arc windows used for local r. Default 16.
    coh_thresh : float
        Local r threshold above which a window is "coherent". Default 0.85.
    incoh_thresh : float
        Local r threshold below which a window is "incoherent". Default 0.70.
    persistence_fraction : float
        Fraction of post-burn frames a window must spend above/below the
        thresholds to count as persistent. Default 0.90.
    lag : int
        Ring-index lag for pos_vel_ac primary. Default 4.
    burn_fraction : float
        Fraction of initial frames discarded before measurement. Default 0.30.
    record_dt : float or None
        Wall-time between recorded frames. Default None = read from
        model_metadata['record_dt']. Falls back to 1.0 if neither is
        available.
    n_permutations : int
        Position-shuffle surrogate null permutations. Minimum 199 for
        CONFIRMATION (need floor p < 0.01 at n_permutations = 199 which
        gives floor p = (0+1)/(199+1) = 0.005). Default 199.
    seed : int
        RNG seed for the surrogate null. Default 42.
    """

    # Tier thresholds (primary = pos_vel_ac[lag=4])
    _SCREEN_PRIMARY_MIN = 0.55
    _CONFIRM_PRIMARY_MIN = 0.75

    # Prerequisites
    _MIN_N_OSCILLATORS = 32
    _MIN_RUN_LENGTH_POST_BURN = 30

    # Defaults when metadata is absent
    _DEFAULT_RECORD_DT = 1.0

    def __init__(
        self,
        n_windows: int = 16,
        coh_thresh: float = 0.85,
        incoh_thresh: float = 0.70,
        persistence_fraction: float = 0.90,
        lag: int = 4,
        burn_fraction: float = 0.30,
        record_dt: float | None = None,
        n_permutations: int = 199,
        seed: int = 42,
    ) -> None:
        super().__init__(
            pattern_id="P10",
            # P9 (global sync) is the nearest neighbor. P10 shares the
            # oscillator substrate with P9, so exclusion is done at
            # content level (coexistence gate and pos_vel_ac), not by
            # substrate mismatch.
            excluded_patterns=["P9"],
            allowed_co_occurrences=[],
            observable_scope="model_metadata_assisted",
        )
        if n_windows < 4:
            raise ValueError(f"n_windows must be >= 4, got {n_windows}")
        if not (0.0 < coh_thresh <= 1.0):
            raise ValueError(f"coh_thresh must be in (0, 1], got {coh_thresh}")
        if not (0.0 <= incoh_thresh < 1.0):
            raise ValueError(
                f"incoh_thresh must be in [0, 1), got {incoh_thresh}"
            )
        if incoh_thresh >= coh_thresh:
            raise ValueError(
                "incoh_thresh must be < coh_thresh, got "
                f"{incoh_thresh} vs {coh_thresh}"
            )
        if not (0.0 < persistence_fraction <= 1.0):
            raise ValueError(
                f"persistence_fraction must be in (0, 1], "
                f"got {persistence_fraction}"
            )
        if lag < 1:
            raise ValueError(f"lag must be >= 1, got {lag}")
        if not (0.0 <= burn_fraction < 1.0):
            raise ValueError(
                f"burn_fraction must be in [0, 1), got {burn_fraction}"
            )
        if record_dt is not None and record_dt <= 0:
            raise ValueError(f"record_dt must be positive, got {record_dt}")
        if n_permutations < 1:
            raise ValueError(
                f"n_permutations must be >= 1, got {n_permutations}"
            )

        self.n_windows = n_windows
        self.coh_thresh = coh_thresh
        self.incoh_thresh = incoh_thresh
        self.persistence_fraction = persistence_fraction
        self.lag = lag
        self.burn_fraction = burn_fraction
        self.record_dt = record_dt
        self.n_permutations = n_permutations
        self._seed = seed

        # Caches populated during detect()
        self._theta_post: np.ndarray | None = None
        self._lr_matrix: np.ndarray | None = None
        self._vel: np.ndarray | None = None
        self._resolved_record_dt: float = self._DEFAULT_RECORD_DT
        self._screening_rejection_reason: str = "none"

    # ------------------------------------------------------------------
    # Entry-point override — resolve parameters before base class runs
    # ------------------------------------------------------------------

    def detect(
        self,
        state_history: list[dict[str, Any]],
        model_metadata: dict[str, Any] | None = None,
        timescale: float | None = None,
    ):
        """Resolve record_dt before prereq validation runs."""
        if self.record_dt is not None:
            self._resolved_record_dt = float(self.record_dt)
        elif model_metadata is not None and "record_dt" in model_metadata:
            self._resolved_record_dt = float(model_metadata["record_dt"])
        else:
            self._resolved_record_dt = self._DEFAULT_RECORD_DT
        # Reset caches between calls
        self._theta_post = None
        self._lr_matrix = None
        self._vel = None
        self._screening_rejection_reason = "none"
        return super().detect(state_history, model_metadata, timescale)

    # ------------------------------------------------------------------
    # Timescale / prerequisites
    # ------------------------------------------------------------------

    def _estimate_timescale(
        self,
        state_history: list[dict[str, Any]],
        model_metadata: dict[str, Any] | None,
    ) -> float:
        """Use record_dt as the unit timescale; 1 frame = 1 τ."""
        return 1.0

    def _validate_prerequisites(
        self,
        state_history: list[dict[str, Any]],
        timescale: float,
    ) -> list[str]:
        warnings: list[str] = []
        self._screening_rejection_reason = "none"

        if not state_history:
            warnings.append("empty state_history")
            self._screening_rejection_reason = "empty_state_history"
            return warnings

        last = state_history[-1]
        if "theta" not in last:
            warnings.append(
                "no 'theta' observable in state_history — P10 requires "
                "per-oscillator phases"
            )
            self._screening_rejection_reason = "substrate_mismatch"
            return warnings

        theta0 = np.asarray(last["theta"])
        if theta0.ndim != 1:
            warnings.append(
                f"theta shape {theta0.shape} — P10 requires 1-D phase arrays"
            )
            self._screening_rejection_reason = "substrate_mismatch"
            return warnings

        N = int(theta0.size)
        if N < self._MIN_N_OSCILLATORS:
            warnings.append(
                f"N = {N} < minimum {self._MIN_N_OSCILLATORS} for P10"
            )
            self._screening_rejection_reason = "too_few_oscillators"
            return warnings

        T = len(state_history)
        burn = int(T * self.burn_fraction)
        post_burn = T - burn
        if post_burn < self._MIN_RUN_LENGTH_POST_BURN:
            warnings.append(
                f"post-burn run length = {post_burn} < minimum "
                f"{self._MIN_RUN_LENGTH_POST_BURN} for P10"
            )
            self._screening_rejection_reason = "run_too_short"
            return warnings

        # Populate theta_post cache for downstream computations
        try:
            theta_post = np.stack([
                np.asarray(h["theta"], dtype=np.float64)
                for h in state_history[burn:]
            ])
        except Exception as exc:  # heterogeneous shapes
            warnings.append(
                f"could not stack theta arrays into (T, N) matrix: {exc}"
            )
            self._screening_rejection_reason = "substrate_mismatch"
            return warnings

        if theta_post.ndim != 2 or theta_post.shape[1] != N:
            warnings.append(
                f"theta history shape {theta_post.shape} is not (T, N)"
            )
            self._screening_rejection_reason = "substrate_mismatch"
            return warnings

        self._theta_post = theta_post
        return warnings

    # ------------------------------------------------------------------
    # Primary
    # ------------------------------------------------------------------

    def _compute_primary(
        self,
        state_history: list[dict[str, Any]],
        timescale: float,
    ) -> dict[str, float]:
        if (
            self._screening_rejection_reason != "none"
            or self._theta_post is None
        ):
            return {
                "pos_vel_ac": 0.0,
                "lag": self.lag,
                "n_windows": self.n_windows,
                "n_persistent_coh": 0,
                "n_persistent_incoh": 0,
                "coexistence": False,
                "screening_rejection_reason": self._screening_rejection_reason,
            }

        theta_post = self._theta_post

        # Coexistence gate
        lr = local_r_time_matrix(theta_post, n_windows=self.n_windows)
        self._lr_matrix = lr
        cx = coexistence_stats(
            lr,
            coh_thresh=self.coh_thresh,
            incoh_thresh=self.incoh_thresh,
            persistence_fraction=self.persistence_fraction,
        )

        # Velocities
        try:
            vel = phase_velocities(theta_post, record_dt=self._resolved_record_dt)
        except ValueError:
            return {
                "pos_vel_ac": 0.0,
                "lag": self.lag,
                "n_windows": self.n_windows,
                **{k: cx[k] for k in ("n_persistent_coh",
                                       "n_persistent_incoh",
                                       "coexistence",
                                       "gap_timeavg",
                                       "local_r_min",
                                       "local_r_max")},
                "screening_rejection_reason": "run_too_short",
            }
        self._vel = vel
        ac = pos_vel_ac(vel, lag=self.lag)

        return {
            "pos_vel_ac": float(ac),
            "lag": int(self.lag),
            "n_windows": int(self.n_windows),
            "n_persistent_coh": int(cx["n_persistent_coh"]),
            "n_persistent_incoh": int(cx["n_persistent_incoh"]),
            "per_frame_coexistence_fraction": float(
                cx["per_frame_coexistence_fraction"]
            ),
            "coexistence": bool(cx["coexistence"]),
            "gap_timeavg": float(cx["gap_timeavg"]),
            "local_r_min": float(cx["local_r_min"]),
            "local_r_max": float(cx["local_r_max"]),
            "above_frac_max": float(cx["above_frac_max"]),
            "below_frac_max": float(cx["below_frac_max"]),
            "screening_rejection_reason": "none",
        }

    def _check_screening(
        self,
        primary_result: dict[str, float],
        timescale: float,
    ) -> bool:
        if self._screening_rejection_reason != "none":
            return False
        if not primary_result.get("coexistence", False):
            self._screening_rejection_reason = "no_coexistence"
            primary_result["screening_rejection_reason"] = "no_coexistence"
            return False
        ac = float(primary_result.get("pos_vel_ac", 0.0))
        if ac <= self._SCREEN_PRIMARY_MIN:
            self._screening_rejection_reason = "pos_vel_ac_below_floor"
            primary_result["screening_rejection_reason"] = "pos_vel_ac_below_floor"
            return False
        return True

    # ------------------------------------------------------------------
    # Secondaries
    # ------------------------------------------------------------------

    def _compute_secondaries(
        self,
        state_history: list[dict[str, Any]],
        timescale: float,
    ) -> dict[str, Any]:
        if self._theta_post is None or self._vel is None or self._lr_matrix is None:
            return {}
        r_mean, r_std = global_order_parameter(self._theta_post)
        return {
            "r_global_mean": r_mean,
            "r_global_std": r_std,
            "persistence_corr": persistence_corr(self._lr_matrix),
            "velocity_neighbor_corr": velocity_spatial_neighbor_corr(self._vel),
        }

    # ------------------------------------------------------------------
    # Null
    # ------------------------------------------------------------------

    def _run_null_model(
        self,
        state_history: list[dict[str, Any]],
        primary_result: dict[str, float],
        timescale: float,
    ) -> tuple[float, NullType, dict[str, float]]:
        """Position-shuffle surrogate null.

        Permute oscillator indices randomly; recompute velocities and
        pos_vel_ac. Preserves per-frame phase distributions and
        per-oscillator trajectories; destroys ring adjacency — which is
        the signal the primary depends on.

        Right-tailed p: fraction of null values >= observed. Floor at
        1 / (n_permutations + 1).
        """
        if self._theta_post is None:
            return 1.0, NullType.SURROGATE, {"mean": 0.0, "std": 0.0}

        observed = float(primary_result.get("pos_vel_ac", 0.0))
        rng = np.random.default_rng(self._seed)
        nulls = shuffle_null_pos_vel_ac(
            self._theta_post,
            record_dt=self._resolved_record_dt,
            n_permutations=self.n_permutations,
            lag=self.lag,
            rng=rng,
        )
        # Right-tailed with +1 smoothing (floor p = 1 / (n+1))
        p = (float((nulls >= observed).sum()) + 1.0) / (
            float(len(nulls)) + 1.0
        )
        return (
            float(p),
            NullType.SURROGATE,
            {
                "mean": float(nulls.mean()),
                "std": float(nulls.std()),
                "observed": float(observed),
            },
        )

    def _compute_effect_size(
        self,
        primary_result: dict[str, float],
        null_dist_stats: dict[str, float],
    ) -> dict[str, float]:
        observed = float(null_dist_stats.get("observed",
                                              primary_result.get("pos_vel_ac", 0.0)))
        nm = float(null_dist_stats.get("mean", 0.0))
        ns = float(null_dist_stats.get("std", 0.0))
        if ns <= 1e-9:
            d = 1e6 if observed > nm else (-1e6 if observed < nm else 0.0)
        else:
            d = (observed - nm) / ns
        return {
            "cohens_d": float(d),
            "raw_value": observed,
            "null_mean": nm,
            "null_std": ns,
        }

    # ------------------------------------------------------------------
    # Tier logic
    # ------------------------------------------------------------------

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
        bonuses: dict[str, bool] = {}

        # Screening bonuses
        bonuses["secondaries_pass"] = bool(secondary_result)
        bonuses["shuffle_null_p_001"] = null_p < 0.01

        ac = float(primary_result.get("pos_vel_ac", 0.0))
        if ac < self._CONFIRM_PRIMARY_MIN or null_p >= 0.01:
            return DetectionTier.SCREENING, bonuses

        # Confirmation bonuses
        bonuses = {
            "null_p_0001": null_p < 0.001,
            "effect_size_gt_1": (
                abs(effect_size.get("cohens_d", 0.0)) > 1.0
            ),
            "all_secondaries": self._all_secondaries_pass(secondary_result),
        }

        # DEFINITIVE gate: non-local coupling metadata flag.
        if model_metadata is None:
            return DetectionTier.CONFIRMATION, bonuses
        has_nl = model_metadata.get("has_nonlocal_coupling")
        has_het = model_metadata.get("has_frequency_heterogeneity")
        # DEFINITIVE requires BOTH:
        #   has_nonlocal_coupling True  (presence of chimera-enabling mechanism)
        #   has_frequency_heterogeneity False or absent (identical oscillators)
        if has_nl is not True:
            return DetectionTier.CONFIRMATION, bonuses
        if has_het is True:
            return DetectionTier.CONFIRMATION, bonuses

        bonuses["all_exclusions_cleared"] = True
        bonuses["both_null_types_rejected"] = null_p < 0.01
        # "finite_size_robustness" — conservatively False unless a larger-N
        # companion run has been verified by the caller. P10 does not
        # rerun the model.
        bonuses["finite_size_robustness"] = False
        return DetectionTier.DEFINITIVE, bonuses

    def _all_secondaries_pass(self, secondary_result: dict[str, Any]) -> bool:
        if not secondary_result:
            return False
        # Require r_global mean in the partial-sync band (not full sync,
        # not full incoherence) — a sanity check that the observed
        # dynamics are chimeric. The coexistence gate already enforces
        # this at the window level; this is a redundant aggregate check.
        r_mean = float(secondary_result.get("r_global_mean", 0.0))
        return 0.3 < r_mean < 0.95

    # ------------------------------------------------------------------
    # Exclusions
    # ------------------------------------------------------------------

    def _check_exclusions(
        self,
        state_history: list[dict[str, Any]],
        model_metadata: dict[str, Any] | None,
        timescale: float,
    ) -> tuple[list[str], dict[str, str]]:
        """Exclude P9 (global sync).

        P9 and P10 share the oscillator substrate, so exclusion happens
        at content level. The coexistence gate and pos_vel_ac primary
        BOTH have to pass; P9's canonical positive (K > K_c ordinary
        Kuramoto) fails both (no coexistence AND low pos_vel_ac). So if
        we reach CONFIRMATION+, P9 is already content-excluded.
        """
        return ["P9"], {"P9": "excluded"}
