"""Detector result schema and confidence scoring.

Every detector returns a DetectorResult — a structured assessment, not just
pass/fail. Confidence is tier-capped: screening can never exceed 0.60,
confirmation 0.85, definitive 1.00. A strong screening signal must pass
confirmation criteria to unlock the higher confidence range.

See docs/detector_cards.md § Detector Output Schema for the specification.
"""

from __future__ import annotations

from dataclasses import dataclass, field
from enum import Enum
from typing import Any


class DetectionTier(Enum):
    """Detection evidence tiers, in ascending order of rigor.

    - SCREENING: Exploratory flagging. Primary metric above relaxed threshold.
    - CONFIRMATION: Standard detection. Primary + secondary + null model p < 0.01.
    - DEFINITIVE: Publication-grade. All confirmation + exclusions + both null types.
    """
    SCREENING = "screening"
    CONFIRMATION = "confirmation"
    DEFINITIVE = "definitive"

    @property
    def confidence_cap(self) -> float:
        return _TIER_CAPS[self]

    def __lt__(self, other: DetectionTier) -> bool:
        return _TIER_ORDER[self] < _TIER_ORDER[other]

    def __le__(self, other: DetectionTier) -> bool:
        return _TIER_ORDER[self] <= _TIER_ORDER[other]

    def __gt__(self, other: DetectionTier) -> bool:
        return _TIER_ORDER[self] > _TIER_ORDER[other]

    def __ge__(self, other: DetectionTier) -> bool:
        return _TIER_ORDER[self] >= _TIER_ORDER[other]


_TIER_CAPS = {
    DetectionTier.SCREENING: 0.60,
    DetectionTier.CONFIRMATION: 0.85,
    DetectionTier.DEFINITIVE: 1.00,
}

_TIER_ORDER = {
    DetectionTier.SCREENING: 0,
    DetectionTier.CONFIRMATION: 1,
    DetectionTier.DEFINITIVE: 2,
}

# Base confidence scores and bonus increments per tier.
# From detector_cards.md § Confidence score table.
CONFIDENCE_BASES = {
    DetectionTier.SCREENING: 0.35,
    DetectionTier.CONFIRMATION: 0.55,
    DetectionTier.DEFINITIVE: 0.75,
}

CONFIDENCE_BONUSES = {
    DetectionTier.SCREENING: {
        "secondaries_pass": 0.15,
        "shuffle_null_p_001": 0.10,
    },
    DetectionTier.CONFIRMATION: {
        "null_p_0001": 0.15,
        "effect_size_gt_1": 0.10,
        "all_secondaries": 0.05,
    },
    DetectionTier.DEFINITIVE: {
        "all_exclusions_cleared": 0.10,
        "both_null_types_rejected": 0.10,
        "finite_size_robustness": 0.05,
    },
}


def compute_confidence(
    tier: DetectionTier,
    bonuses_earned: dict[str, bool],
) -> float:
    """Compute tier-capped confidence score.

    Parameters
    ----------
    tier : DetectionTier
        The highest tier whose criteria are met.
    bonuses_earned : dict[str, bool]
        Keys matching CONFIDENCE_BONUSES[tier]; True if the bonus condition is met.

    Returns
    -------
    float
        Confidence in [0.0, tier.confidence_cap].
    """
    base = CONFIDENCE_BASES[tier]
    bonus_defs = CONFIDENCE_BONUSES[tier]
    total_bonus = sum(
        bonus_defs[k] for k, earned in bonuses_earned.items()
        if earned and k in bonus_defs
    )
    return min(base + total_bonus, tier.confidence_cap)


class NullType(Enum):
    """Null model evidence strength categories.

    - SHUFFLE: label/timing permutation of observed data.
    - SURROGATE: synthetic data matching marginal statistics.
    - MECHANISTIC_INTERVENTION: requires modifying model rules.
    """
    SHUFFLE = "shuffle"
    SURROGATE = "surrogate"
    MECHANISTIC_INTERVENTION = "mechanistic_intervention"


@dataclass
class DetectorResult:
    """Structured result from any pattern detector.

    This is the universal output schema. Every field is populated by the
    detector; downstream analysis (orchestration, transfer matrix, reports)
    consumes this without knowing which detector produced it.
    """

    # Identity
    pattern_id: str
    """Pattern identifier, e.g. 'P1', 'P31'."""

    detected: bool
    """Whether the pattern passes at least screening tier."""

    tier: DetectionTier
    """Highest tier achieved."""

    confidence: float
    """Tier-capped confidence score in [0.0, 1.0]."""

    # Measurements
    primary_metric: dict[str, float]
    """Primary detection metric, e.g. {'morans_i': 0.84}."""

    secondary_metrics: dict[str, Any]
    """Supporting measurements, pattern-specific."""

    effect_size: dict[str, float]
    """Cross-model comparison. Recommended: cohens_d, raw metric vs null mean."""

    # Null model assessment
    null_p_value: float
    """P-value from primary null model."""

    null_type: NullType
    """Which null model produced the p-value."""

    # Exclusion and co-occurrence
    exclusions_checked: list[str] = field(default_factory=list)
    """Neighbor patterns tested for exclusion, e.g. ['P2', 'P3', 'P30']."""

    exclusion_results: dict[str, str] = field(default_factory=dict)
    """Exclusion outcomes, e.g. {'P2': 'excluded', 'P3': 'excluded'}."""

    co_occurrence_candidates: list[str] = field(default_factory=list)
    """Compatible patterns also detected in the same run."""

    # Context
    metadata_available: bool = False
    """Whether model metadata (rules, parameters) was available."""

    warnings: list[str] = field(default_factory=list)
    """Flags: low N, short run, boundary artifacts, etc."""

    notes: str = ""
    """Free text for edge cases."""

    def __post_init__(self) -> None:
        if self.confidence > self.tier.confidence_cap:
            self.confidence = self.tier.confidence_cap
        if self.confidence < 0.0:
            self.confidence = 0.0

    @property
    def is_confirmed(self) -> bool:
        return self.detected and self.tier >= DetectionTier.CONFIRMATION

    @property
    def is_definitive(self) -> bool:
        return self.detected and self.tier >= DetectionTier.DEFINITIVE

    def summary(self) -> str:
        """One-line human-readable summary."""
        status = "DETECTED" if self.detected else "not detected"
        return (
            f"{self.pattern_id}: {status} at {self.tier.value} "
            f"(confidence={self.confidence:.2f}, "
            f"p={self.null_p_value:.4f})"
        )
