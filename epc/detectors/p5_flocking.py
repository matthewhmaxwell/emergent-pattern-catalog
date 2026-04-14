"""
P5 — Translational alignment (flocking) detector.

Implements the full P5 detector card specification:
    - Primary metric: Polarization order parameter φ
    - Secondary metrics: Group speed ratio R, angular momentum |L|,
      heading autocorrelation, heading distribution
    - Null model: Heading-shuffle permutation test
    - Exclusions: P6 (milling), P7 (lanes)
    - Three detection tiers: screening, confirmation, definitive

Detector card reference:
    Screening: φ_mean > 0.5 over ≥ 5 × T_cross
    Confirmation: φ_mean > 0.7 over ≥ 10 × T_cross AND R > 0.5
                  AND above shuffle null (p < 0.01)
    Definitive: Confirmation + P6/P7/P8 exclusions
                + zero-coupling null shows φ ≈ 1/√N

Power requirements:
    - Shuffle null: ≥ 199 permutations (floor p = 0.005)
    - History must cover ≥ 5 × T_cross for screening,
      ≥ 10 × T_cross for confirmation
"""

from __future__ import annotations

import numpy as np
from dataclasses import dataclass, field
from typing import Any, Optional

from epc.metrics.collective_motion import (
    PolarizationMetric,
    GroupSpeedRatioMetric,
    AngularMomentumMetric,
    HeadingDistributionMetric,
)


# ===================================================================
# DetectorResult (matches project-wide schema)
# ===================================================================

@dataclass
class DetectorResult:
    """Standardized detector output."""

    pattern_id: str
    detected: bool
    tier: str  # "none" | "screening" | "confirmation" | "definitive"
    confidence: float
    primary_metric: dict
    secondary_metrics: dict
    effect_size: dict
    null_p_value: float
    null_type: str
    exclusions_checked: list
    exclusion_results: dict
    co_occurrence_candidates: list
    metadata_available: bool
    warnings: list = field(default_factory=list)
    notes: str = ""


# ===================================================================
# Null model: heading-shuffle permutation test
# ===================================================================

def heading_shuffle_null(
    history: list[dict[str, Any]],
    n_permutations: int = 199,
    measurement_window: Optional[tuple[int, int]] = None,
    rng: Optional[np.random.Generator] = None,
) -> dict[str, Any]:
    """Heading-shuffle null model for polarization.

    Replaces all headings with independent uniform random directions
    while keeping positions and speeds. This tests whether the
    observed alignment could arise from random heading assignment.

    Under this null, φ ≈ 1/√N (central limit theorem for unit vectors).

    Args:
        history: State history with 'velocities' key.
        n_permutations: Number of Monte Carlo samples. ≥199 for
            confirmation-tier p < 0.005.
        measurement_window: (start, end) indices. If None, use full history.
        rng: Random generator for reproducibility.

    Returns:
        Dict with:
            null_phi_values: (n_permutations,) mean φ per sample
            observed_phi: mean φ of original data
            p_value: fraction of null values ≥ observed
            effect_size_d: Cohen's d (observed - null_mean) / null_std
    """
    if rng is None:
        rng = np.random.default_rng(42)

    # Select measurement window
    if measurement_window is not None:
        start, end = measurement_window
        history = history[start:end]

    # Observed polarization time series
    observed_phis = np.array([
        PolarizationMetric.compute_instant(s) for s in history
    ])
    observed_phi = np.mean(observed_phis)

    # Get N from first state
    N = len(history[0]["velocities"])
    T = len(history)

    # Generate null distribution
    # For each permutation: at each timestep, assign independent random
    # headings to all particles, compute φ, then average over time.
    null_phi_values = np.empty(n_permutations)

    for perm in range(n_permutations):
        # Draw random headings for each (timestep, particle)
        random_headings = rng.uniform(-np.pi, np.pi, size=(T, N))
        # φ at each timestep: |(1/N) Σ (cos θ, sin θ)|
        mean_cos = np.mean(np.cos(random_headings), axis=1)  # (T,)
        mean_sin = np.mean(np.sin(random_headings), axis=1)  # (T,)
        phi_t = np.sqrt(mean_cos**2 + mean_sin**2)           # (T,)
        null_phi_values[perm] = np.mean(phi_t)

    # p-value: fraction of null values ≥ observed (one-sided, upper)
    p_value = (np.sum(null_phi_values >= observed_phi) + 1) / (n_permutations + 1)

    # Effect size (Cohen's d)
    null_mean = np.mean(null_phi_values)
    null_std = np.std(null_phi_values)
    if null_std > 1e-12:
        cohens_d = (observed_phi - null_mean) / null_std
    else:
        cohens_d = float("inf") if observed_phi > null_mean else 0.0

    return {
        "null_phi_values": null_phi_values,
        "null_mean": float(null_mean),
        "null_std": float(null_std),
        "observed_phi": float(observed_phi),
        "p_value": float(p_value),
        "cohens_d": float(cohens_d),
        "n_permutations": n_permutations,
    }


# ===================================================================
# P5 Flocking Detector
# ===================================================================

class P5FlockingDetector:
    """Detector for P5 — Translational alignment (flocking).

    Follows the three-tier detection protocol:
        1. Screening: φ > 0.5 sustained over ≥ 5 T_cross
        2. Confirmation: φ > 0.7 over ≥ 10 T_cross, R > 0.5,
           shuffle null p < 0.01
        3. Definitive: Confirmation + P6/P7 exclusions +
           zero-coupling null

    Args:
        n_permutations: Permutation count for shuffle null.
            Default 199 → floor p = 0.005.
        T_cross: Domain crossing time in timesteps.
            If None, estimated from metadata or history.
        box_size: Domain size for periodic-aware angular momentum.
            If None, skips periodic correction.
        seed: Random seed for permutation test.
    """

    # --- Detector card thresholds ---
    SCREENING_PHI = 0.5
    SCREENING_T_CROSS_MIN = 5  # multiples of T_cross

    CONFIRMATION_PHI = 0.7
    CONFIRMATION_T_CROSS_MIN = 10
    CONFIRMATION_R_MIN = 0.5
    CONFIRMATION_P_MAX = 0.01

    P6_EXCLUSION_L_THRESHOLD = 0.5
    P7_EXCLUSION_ANTIPARALLEL_THRESHOLD = 0.2

    def __init__(
        self,
        n_permutations: int = 199,
        T_cross: Optional[float] = None,
        box_size: Optional[float] = None,
        seed: int = 42,
    ):
        self.n_permutations = n_permutations
        self.T_cross = T_cross
        self.box_size = box_size
        self.rng = np.random.default_rng(seed)

    def detect(
        self,
        history: list[dict[str, Any]],
        metadata: Optional[dict[str, Any]] = None,
    ) -> DetectorResult:
        """Run full P5 detection pipeline.

        Args:
            history: State history — list of dicts with 'positions'
                and 'velocities' keys.
            metadata: Optional model metadata dict. Used for T_cross
                estimation and zero-coupling null context.

        Returns:
            DetectorResult with tier, confidence, and all metrics.
        """
        warnings = []
        notes_parts = []

        # --- Determine T_cross ---
        T_cross = self._resolve_T_cross(history, metadata, warnings)

        # --- Determine measurement window ---
        T = len(history)
        T_cross_steps = int(T_cross)

        # Check minimum run length
        min_steps_screening = self.SCREENING_T_CROSS_MIN * T_cross_steps
        min_steps_confirmation = self.CONFIRMATION_T_CROSS_MIN * T_cross_steps

        if T < min_steps_screening:
            warnings.append(
                f"Run length ({T}) < {self.SCREENING_T_CROSS_MIN} T_cross "
                f"({min_steps_screening}). Cannot reach screening tier."
            )

        # Use second half of trajectory as measurement window
        # (first half as equilibration)
        meas_start = max(T // 2, 0)
        meas_end = T
        meas_length = meas_end - meas_start
        meas_history = history[meas_start:meas_end]

        n_T_cross = meas_length / T_cross if T_cross > 0 else 0
        notes_parts.append(
            f"Measurement window: steps {meas_start}–{meas_end} "
            f"({n_T_cross:.1f} T_cross)"
        )

        # --- Compute primary metric: polarization ---
        phi_result = PolarizationMetric.compute(meas_history)
        phi_mean = phi_result["phi_mean"]

        # --- Compute secondary metrics ---
        R_result = GroupSpeedRatioMetric.compute(meas_history)
        R_mean = R_result["R_mean"]

        L_result = AngularMomentumMetric.compute(meas_history, self.box_size)
        L_abs_mean = L_result["L_abs_mean"]

        heading_dist = HeadingDistributionMetric.compute(meas_history)

        # --- Null model (shuffle) ---
        null_result = heading_shuffle_null(
            meas_history,
            n_permutations=self.n_permutations,
            rng=self.rng,
        )
        null_p = null_result["p_value"]
        cohens_d = null_result["cohens_d"]

        # --- Exclusion checks ---
        exclusion_results = {}

        # P6: |L| > 0.5 → milling
        p6_excluded = L_abs_mean > self.P6_EXCLUSION_L_THRESHOLD
        exclusion_results["P6"] = "TRIGGERED" if p6_excluded else "excluded"

        # P7: bimodal antiparallel heading → lanes
        p7_excluded = heading_dist["is_bimodal_antiparallel"]
        exclusion_results["P7"] = "TRIGGERED" if p7_excluded else "excluded"

        # P8: backward density waves → jamming (not implemented yet)
        exclusion_results["P8"] = "not_checked"
        if exclusion_results["P8"] == "not_checked":
            warnings.append("P8 exclusion (jamming) not yet implemented")

        any_exclusion_triggered = p6_excluded or p7_excluded
        all_implemented_exclusions_clear = (
            not p6_excluded and not p7_excluded
        )

        # --- Tier determination ---
        tier = "none"
        detected = False
        confidence = 0.0

        # Screening: φ_mean > 0.5, ≥ 5 T_cross
        screening_pass = (
            phi_mean > self.SCREENING_PHI
            and n_T_cross >= self.SCREENING_T_CROSS_MIN
        )

        # Confirmation: φ_mean > 0.7, R > 0.5, ≥ 10 T_cross, null p < 0.01
        confirmation_pass = (
            phi_mean > self.CONFIRMATION_PHI
            and R_mean > self.CONFIRMATION_R_MIN
            and n_T_cross >= self.CONFIRMATION_T_CROSS_MIN
            and null_p < self.CONFIRMATION_P_MAX
        )

        # Definitive: confirmation + exclusions clear
        definitive_pass = (
            confirmation_pass
            and all_implemented_exclusions_clear
        )

        if definitive_pass:
            tier = "definitive"
            detected = True
            confidence = 0.75  # base
            if all_implemented_exclusions_clear:
                confidence += 0.10
            if null_p < 0.001:
                confidence += 0.10
            # Cap at 1.0
            confidence = min(confidence, 1.0)
            # Note: zero-coupling null not run in this context
            # (would require model access). Documented as limitation.
            notes_parts.append(
                "Zero-coupling null requires model access — "
                "not run in state-history-only mode."
            )

        elif confirmation_pass:
            tier = "confirmation"
            detected = True
            confidence = 0.55  # base
            if null_p < 0.001:
                confidence += 0.15
            if cohens_d > 1.0:
                confidence += 0.10
            if R_mean > self.CONFIRMATION_R_MIN:
                confidence += 0.05
            confidence = min(confidence, 0.85)

            if any_exclusion_triggered:
                notes_parts.append(
                    f"Exclusion(s) triggered: {exclusion_results}. "
                    f"Confirmation tier, not definitive."
                )

        elif screening_pass:
            tier = "screening"
            detected = True
            confidence = 0.35  # base
            if R_mean > self.CONFIRMATION_R_MIN:
                confidence += 0.15  # secondary pass
            if null_p < 0.01:
                confidence += 0.10
            confidence = min(confidence, 0.60)

        else:
            tier = "none"
            detected = False
            confidence = 0.0

        # --- Assemble result ---
        N = len(history[0]["velocities"])
        return DetectorResult(
            pattern_id="P5",
            detected=detected,
            tier=tier,
            confidence=round(confidence, 3),
            primary_metric={"polarization_mean": phi_mean},
            secondary_metrics={
                "group_speed_ratio": R_mean,
                "angular_momentum_abs": L_abs_mean,
                "heading_resultant_length": heading_dist["resultant_length"],
                "heading_circular_variance": heading_dist["circular_variance"],
                "antiparallel_fraction": heading_dist["antiparallel_fraction"],
                "nematic_order": heading_dist["nematic_order"],
                "measurement_T_cross": n_T_cross,
            },
            effect_size={
                "cohens_d": cohens_d,
                "phi_observed": phi_mean,
                "phi_null_mean": null_result["null_mean"],
                "phi_null_std": null_result["null_std"],
            },
            null_p_value=null_p,
            null_type="shuffle",
            exclusions_checked=["P6", "P7", "P8"],
            exclusion_results=exclusion_results,
            co_occurrence_candidates=["P17", "P19", "P9"],
            metadata_available=metadata is not None,
            warnings=warnings,
            notes=" | ".join(notes_parts),
        )

    def _resolve_T_cross(
        self,
        history: list[dict[str, Any]],
        metadata: Optional[dict[str, Any]],
        warnings: list[str],
    ) -> float:
        """Determine T_cross from metadata, explicit setting, or estimation.

        Priority: explicit > metadata > estimation from history.
        """
        if self.T_cross is not None:
            return self.T_cross

        if metadata is not None:
            if "box_size" in metadata and "speed" in metadata:
                v = metadata["speed"]
                L = metadata["box_size"]
                if v > 0:
                    return L / v

        # Estimate from history: T_cross = L_est / v_mean
        # where L_est = domain size estimated from position range
        if len(history) > 0:
            all_pos = np.array([s["positions"] for s in history[:100]])
            pos_range = np.max(all_pos, axis=(0, 1)) - np.min(all_pos, axis=(0, 1))
            L_est = np.max(pos_range) * 1.1  # rough estimate

            all_vel = np.array([s["velocities"] for s in history[:100]])
            v_mean = np.mean(np.linalg.norm(all_vel, axis=2))

            if v_mean > 1e-12:
                T_est = L_est / v_mean
                warnings.append(
                    f"T_cross estimated from history: {T_est:.1f} steps. "
                    f"Provide metadata or set T_cross explicitly for accuracy."
                )
                return T_est

        warnings.append("Could not determine T_cross. Using T=100 as fallback.")
        return 100.0
