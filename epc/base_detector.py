"""Base detector interface.

Every pattern detector subclasses BaseDetector and implements:
- detect(): run full detection pipeline, return DetectorResult
- _compute_primary(): compute primary metric from state history
- _compute_secondaries(): compute secondary metrics
- _run_null_model(): generate null distribution and p-value
- _check_exclusions(): test nearest-neighbor exclusions

The base class handles tier logic, confidence scoring, and result assembly.
Subclasses focus on pattern-specific computation.
"""

from __future__ import annotations

from abc import ABC, abstractmethod
from typing import Any

import numpy as np

from epc.detector_result import (
    CONFIDENCE_BONUSES,
    DetectionTier,
    DetectorResult,
    NullType,
    compute_confidence,
)


class BaseDetector(ABC):
    """Abstract base for all pattern detectors.

    Parameters
    ----------
    pattern_id : str
        Pattern identifier, e.g. 'P1'.
    excluded_patterns : list[str]
        Nearest-neighbor patterns that must be checked for exclusion
        before declaring definitive detection.
    allowed_co_occurrences : list[str]
        Patterns that can legitimately co-occur with this one.
    observable_scope : str
        One of 'state_history_only', 'model_metadata_assisted',
        'model_metadata_required'.
    """

    def __init__(
        self,
        pattern_id: str,
        excluded_patterns: list[str] | None = None,
        allowed_co_occurrences: list[str] | None = None,
        observable_scope: str = "state_history_only",
    ) -> None:
        self.pattern_id = pattern_id
        self.excluded_patterns = excluded_patterns or []
        self.allowed_co_occurrences = allowed_co_occurrences or []
        self.observable_scope = observable_scope

    def detect(
        self,
        state_history: list[dict[str, Any]],
        model_metadata: dict[str, Any] | None = None,
        timescale: float | None = None,
    ) -> DetectorResult:
        """Run full detection pipeline.

        Parameters
        ----------
        state_history : list[dict[str, Any]]
            Sequence of state snapshots from a model run.
            Keys are standardized observable names (positions, velocities,
            types, grid, etc.).
        model_metadata : dict, optional
            Model rules, parameters, interaction structure. Required for
            some detectors (see observable_scope).
        timescale : float, optional
            System-intrinsic timescale τ in timesteps. If not provided,
            the detector will attempt to estimate it.

        Returns
        -------
        DetectorResult
            Structured detection assessment.
        """
        warnings: list[str] = []

        # Check metadata availability
        metadata_available = model_metadata is not None
        if self.observable_scope == "model_metadata_required" and not metadata_available:
            warnings.append("metadata missing for structural check")

        # Estimate timescale if not provided
        if timescale is None:
            timescale = self._estimate_timescale(state_history, model_metadata)

        # Validate prerequisites
        prereq_warnings = self._validate_prerequisites(state_history, timescale)
        warnings.extend(prereq_warnings)

        # Primary metric
        primary_result = self._compute_primary(state_history, timescale)

        # Screening check
        screening_pass = self._check_screening(primary_result, timescale)
        if not screening_pass:
            return DetectorResult(
                pattern_id=self.pattern_id,
                detected=False,
                tier=DetectionTier.SCREENING,
                confidence=0.0,
                primary_metric=primary_result,
                secondary_metrics={},
                effect_size={},
                null_p_value=1.0,
                null_type=NullType.SHUFFLE,
                metadata_available=metadata_available,
                warnings=warnings,
            )

        # Secondary metrics
        secondary_result = self._compute_secondaries(state_history, timescale)

        # Null model
        null_p, null_type, null_dist_stats = self._run_null_model(
            state_history, primary_result, timescale
        )

        # Effect size (Cohen's d equivalent)
        effect_size = self._compute_effect_size(primary_result, null_dist_stats)

        # Determine tier
        tier, bonuses = self._determine_tier(
            primary_result=primary_result,
            secondary_result=secondary_result,
            null_p=null_p,
            null_type=null_type,
            effect_size=effect_size,
            state_history=state_history,
            model_metadata=model_metadata,
            timescale=timescale,
        )

        # Confidence
        confidence = compute_confidence(tier, bonuses)

        # Exclusions (only computed if confirmation+ tier)
        exclusions_checked: list[str] = []
        exclusion_results: dict[str, str] = {}
        if tier >= DetectionTier.CONFIRMATION:
            exclusions_checked, exclusion_results = self._check_exclusions(
                state_history, model_metadata, timescale
            )

        return DetectorResult(
            pattern_id=self.pattern_id,
            detected=True,
            tier=tier,
            confidence=confidence,
            primary_metric=primary_result,
            secondary_metrics=secondary_result,
            effect_size=effect_size,
            null_p_value=null_p,
            null_type=null_type,
            exclusions_checked=exclusions_checked,
            exclusion_results=exclusion_results,
            co_occurrence_candidates=self.allowed_co_occurrences,
            metadata_available=metadata_available,
            warnings=warnings,
        )

    # --- Abstract methods for subclasses ---

    @abstractmethod
    def _compute_primary(
        self,
        state_history: list[dict[str, Any]],
        timescale: float,
    ) -> dict[str, float]:
        """Compute the primary detection metric.

        Returns dict, e.g. {'morans_i': 0.84, 'z_score': 3.2}.
        """
        ...

    @abstractmethod
    def _check_screening(
        self,
        primary_result: dict[str, float],
        timescale: float,
    ) -> bool:
        """Check if primary metric passes screening threshold."""
        ...

    @abstractmethod
    def _compute_secondaries(
        self,
        state_history: list[dict[str, Any]],
        timescale: float,
    ) -> dict[str, Any]:
        """Compute secondary metrics."""
        ...

    @abstractmethod
    def _run_null_model(
        self,
        state_history: list[dict[str, Any]],
        primary_result: dict[str, float],
        timescale: float,
    ) -> tuple[float, NullType, dict[str, float]]:
        """Run null model and return (p_value, null_type, null_distribution_stats).

        null_distribution_stats should include at minimum:
        {'mean': float, 'std': float} for effect size computation.
        """
        ...

    @abstractmethod
    def _check_exclusions(
        self,
        state_history: list[dict[str, Any]],
        model_metadata: dict[str, Any] | None,
        timescale: float,
    ) -> tuple[list[str], dict[str, str]]:
        """Test nearest-neighbor exclusions.

        Returns (patterns_checked, results_dict).
        results_dict values: 'excluded', 'not_excluded', 'inconclusive'.
        """
        ...

    # --- Overridable methods with defaults ---

    def _estimate_timescale(
        self,
        state_history: list[dict[str, Any]],
        model_metadata: dict[str, Any] | None,
    ) -> float:
        """Estimate system-intrinsic timescale from state history.

        Default: total run length / 10 (conservative). Subclasses should
        override with pattern-appropriate estimation.
        """
        return max(1.0, len(state_history) / 10.0)

    def _validate_prerequisites(
        self,
        state_history: list[dict[str, Any]],
        timescale: float,
    ) -> list[str]:
        """Check data sufficiency and return warning strings.

        Default checks: minimum history length relative to timescale.
        """
        warnings = []
        run_length = len(state_history)
        if run_length < 5 * timescale:
            warnings.append(f"run length ({run_length}) < 5τ ({5 * timescale:.0f})")
        return warnings

    def _compute_effect_size(
        self,
        primary_result: dict[str, float],
        null_dist_stats: dict[str, float],
    ) -> dict[str, float]:
        """Compute effect size (Cohen's d equivalent).

        Primary metric value minus null mean, divided by null SD.
        """
        if not null_dist_stats or null_dist_stats.get("std", 0) == 0:
            return {}

        # Use the first numeric value in primary_result
        primary_values = [v for v in primary_result.values() if isinstance(v, (int, float))]
        if not primary_values:
            return {}

        primary_val = primary_values[0]
        null_mean = null_dist_stats.get("mean", 0.0)
        null_std = null_dist_stats.get("std", 1.0)

        cohens_d = (primary_val - null_mean) / null_std if null_std > 0 else 0.0

        return {
            "cohens_d": cohens_d,
            "raw_value": primary_val,
            "null_mean": null_mean,
            "null_std": null_std,
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
        """Determine highest achieved tier and which bonuses are earned.

        Default implementation uses the standard tier criteria from detector_cards.md.
        Subclasses override _check_confirmation() and _check_definitive() for
        pattern-specific tier logic.
        """
        # Start at screening (we already passed screening check)
        bonuses: dict[str, bool] = {}

        # Check screening bonuses
        has_secondaries = bool(secondary_result)
        bonuses["secondaries_pass"] = has_secondaries
        bonuses["shuffle_null_p_001"] = null_p < 0.01

        # Try confirmation
        if self._check_confirmation(
            primary_result, secondary_result, null_p, timescale
        ):
            tier = DetectionTier.CONFIRMATION
            bonuses = {
                "null_p_0001": null_p < 0.001,
                "effect_size_gt_1": effect_size.get("cohens_d", 0) > 1.0,
                "all_secondaries": self._all_secondaries_pass(secondary_result),
            }

            # Try definitive
            if self._check_definitive(
                primary_result, secondary_result, null_p, null_type,
                state_history, model_metadata, timescale
            ):
                tier = DetectionTier.DEFINITIVE
                bonuses = {
                    "all_exclusions_cleared": True,  # Updated after exclusion check
                    "both_null_types_rejected": (
                        null_type == NullType.MECHANISTIC_INTERVENTION and null_p < 0.01
                    ),
                    "finite_size_robustness": self._check_finite_size_robustness(
                        state_history, timescale
                    ),
                }
        else:
            tier = DetectionTier.SCREENING

        return tier, bonuses

    def _check_confirmation(
        self,
        primary_result: dict[str, float],
        secondary_result: dict[str, Any],
        null_p: float,
        timescale: float,
    ) -> bool:
        """Check confirmation tier criteria. Override per pattern."""
        # Default: primary passes + at least one secondary + null p < 0.01
        return bool(secondary_result) and null_p < 0.01

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
        """Check definitive tier criteria. Override per pattern."""
        # Default: confirmation + mechanistic null available
        return null_p < 0.01 and null_type == NullType.MECHANISTIC_INTERVENTION

    def _all_secondaries_pass(self, secondary_result: dict[str, Any]) -> bool:
        """Check whether all secondary metrics pass. Override per pattern."""
        return bool(secondary_result)

    def _check_finite_size_robustness(
        self,
        state_history: list[dict[str, Any]],
        timescale: float,
    ) -> bool:
        """Check finite-size and robustness conditions. Override per pattern."""
        return len(state_history) >= 10 * timescale
