"""
P19 — Emergent leadership / minority guidance detector.

Implements the P19 detector card specification:
    - Primary metric: Group directional accuracy (alignment with preferred dir)
    - Confirmation metric: Influence asymmetry via directional pull —
      informed agents' mean heading is consistently biased toward the
      preferred direction relative to naive agents
    - Definitive metric: Guidance efficacy disproportionate to informed share
    - Null model: label-shuffle permutation test
    - Three detection tiers: screening, confirmation, definitive

Detection logic:
    The defining P19 signature is that a small informed minority steers the
    whole group without explicit signaling (Couzin et al. 2005, Nature).

    - Screening (≤0.60): Group mean heading aligns with the preferred
      direction above a threshold (accuracy > 0.3).
    - Confirmation (≤0.85): Informed agents' mean heading is closer to the
      preferred direction than naive agents' mean heading, significantly
      more often than expected by label-shuffle permutation null. This
      directly measures the directional influence asymmetry: the informed
      minority "leads" the group by being consistently biased toward the
      preferred direction.
    - Definitive (≤1.00): Guidance efficacy is disproportionate to informed
      share (accuracy / informed_fraction >> 1) AND the informed fraction
      is a genuine minority (ρ ≤ 0.25).

Canonical parameters (Couzin 2005):
    N=200, ρ=0.1, ω=0.3, preferred_direction=0.0, noise=0.1

Power requirements:
    - Label-shuffle null: ≥ 99 permutations
    - History length: ≥ 100 steps

References:
    Couzin, I. D., Krause, J., Franks, N. R., & Levin, S. A. (2005).
    Effective leadership and decision-making in animal groups on the move.
    Nature, 433(7025), 513–516.
"""

from __future__ import annotations

import numpy as np
from dataclasses import dataclass, field
from typing import Any, Optional


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


def _compute_influence_asymmetry(
    history: list[dict[str, Any]],
    informed_mask: np.ndarray,
    preferred_direction: float,
    window: tuple[int, int],
    n_permutations: int = 199,  # 99 floors p at 0.01 so the p<0.01 gates would be unreachable
    seed: int = 42,
) -> dict[str, Any]:
    """Compute influence asymmetry: how much closer are informed agents'
    headings to the preferred direction vs naive agents' headings.

    At each timestep, compute:
      pull = cos(θ_informed_mean - θ_pref) - cos(θ_naive_mean - θ_pref)

    If pull > 0 consistently, informed agents are leading (biased toward
    preferred direction). Test significance via label-shuffle permutation:
    randomly reassign which agents are "informed" and check if the pull
    is as strong.

    Args:
        history: State history with 'headings' key.
        informed_mask: Boolean mask for informed agents.
        preferred_direction: Target direction (radians).
        window: (start, end) indices.
        n_permutations: Permutations for significance test.
        seed: RNG seed.

    Returns:
        Dict with pull_mean, pull_fraction, p_value.
    """
    naive_mask = ~informed_mask
    n_informed = int(np.sum(informed_mask))
    n_total = len(informed_mask)

    if n_informed == 0 or n_informed == n_total:
        return {
            "pull_mean": 0.0,
            "pull_fraction": 0.0,
            "p_value": 1.0,
        }

    start, end = window
    T = end - start

    # Compute pull at each timestep
    pulls = np.zeros(T)
    for i, t in enumerate(range(start, end)):
        h = history[t]["headings"]
        # Circular mean for informed
        inf_h = h[informed_mask]
        inf_mean = np.arctan2(np.mean(np.sin(inf_h)), np.mean(np.cos(inf_h)))
        # Circular mean for naive
        nai_h = h[naive_mask]
        nai_mean = np.arctan2(np.mean(np.sin(nai_h)), np.mean(np.cos(nai_h)))
        # Pull: informed closer to preferred than naive
        pulls[i] = (
            np.cos(inf_mean - preferred_direction)
            - np.cos(nai_mean - preferred_direction)
        )

    observed_pull_mean = float(np.mean(pulls))
    observed_pull_fraction = float(np.mean(pulls > 0))

    # Permutation test: shuffle informed/naive labels
    rng = np.random.default_rng(seed)
    null_pull_means = np.zeros(n_permutations)

    for p in range(n_permutations):
        perm = rng.permutation(n_total)
        fake_informed = np.zeros(n_total, dtype=bool)
        fake_informed[perm[:n_informed]] = True
        fake_naive = ~fake_informed

        null_pulls = np.zeros(T)
        for i, t in enumerate(range(start, end)):
            h = history[t]["headings"]
            fi_h = h[fake_informed]
            fi_mean = np.arctan2(np.mean(np.sin(fi_h)), np.mean(np.cos(fi_h)))
            fn_h = h[fake_naive]
            fn_mean = np.arctan2(np.mean(np.sin(fn_h)), np.mean(np.cos(fn_h)))
            null_pulls[i] = (
                np.cos(fi_mean - preferred_direction)
                - np.cos(fn_mean - preferred_direction)
            )
        null_pull_means[p] = np.mean(null_pulls)

    p_value = float(np.mean(null_pull_means >= observed_pull_mean))
    if p_value == 0:
        p_value = 1.0 / (n_permutations + 1)

    return {
        "pull_mean": observed_pull_mean,
        "pull_fraction": observed_pull_fraction,
        "p_value": p_value,
    }


def detect_p19(
    history: list[dict[str, Any]],
    metadata: dict[str, Any],
    n_permutations: int = 199,  # 99 floors p at 0.01 so the p<0.01 gates would be unreachable
    seed: int = 42,
) -> DetectorResult:
    """Run the P19 emergent leadership / minority guidance detector.

    Args:
        history: State history from an informed-minority model.
            Must contain 'headings' and 'informed_mask' per state.
        metadata: Model metadata. Must contain 'preferred_direction',
            'informed_fraction'.
        n_permutations: Number of permutations for label-shuffle null.
        seed: RNG seed.

    Returns:
        DetectorResult with tier, confidence, and metrics.
    """
    warnings: list[str] = []

    # --- Prerequisites ---
    if len(history) < 10:
        return DetectorResult(
            pattern_id="P19",
            detected=False,
            tier="none",
            confidence=0.0,
            primary_metric={"group_accuracy": 0.0},
            secondary_metrics={},
            effect_size={},
            null_p_value=1.0,
            null_type="none",
            exclusions_checked=[],
            exclusion_results={},
            co_occurrence_candidates=[],
            metadata_available=bool(metadata),
            warnings=["History too short (< 10 steps)"],
        )

    # Check for informed_mask in history
    if "informed_mask" not in history[0]:
        return DetectorResult(
            pattern_id="P19",
            detected=False,
            tier="none",
            confidence=0.0,
            primary_metric={"group_accuracy": 0.0},
            secondary_metrics={},
            effect_size={},
            null_p_value=1.0,
            null_type="none",
            exclusions_checked=[],
            exclusion_results={},
            co_occurrence_candidates=[],
            metadata_available=bool(metadata),
            warnings=["No informed_mask in history"],
        )

    informed_mask = history[0]["informed_mask"]
    n_informed = int(np.sum(informed_mask))
    n_total = len(informed_mask)
    informed_fraction = n_informed / n_total if n_total > 0 else 0.0

    preferred_direction = metadata.get("preferred_direction", 0.0)

    # --- Measurement windows ---
    T = len(history)
    # Accuracy: last 50% (steady-state)
    warmup = T // 2
    acc_window = (warmup, T)
    # Influence: use convergence phase (first 20% to 80%) for pull measurement
    pull_start = max(1, T // 5)
    pull_end = min(T, int(T * 0.8))
    pull_window = (pull_start, pull_end)

    # === PRIMARY METRIC: Group directional accuracy ===
    accuracies = []
    for t in range(acc_window[0], acc_window[1]):
        h = history[t]["headings"]
        mean_vx = np.mean(np.cos(h))
        mean_vy = np.mean(np.sin(h))
        group_heading = np.arctan2(mean_vy, mean_vx)
        acc = float(np.cos(group_heading - preferred_direction))
        accuracies.append(acc)

    mean_accuracy = float(np.mean(accuracies))
    std_accuracy = float(np.std(accuracies))

    # === Influence asymmetry (label-shuffle null): the CAUSAL discriminator.
    #     Informed agents must STEER the group, not merely co-move with a supplied
    #     preferred direction. A leaderless flock with a random informed-mask gives
    #     pull ~ 0 and is not significant. Computed BEFORE screening so it gates it.
    influence = _compute_influence_asymmetry(
        history, informed_mask, preferred_direction, pull_window,
        n_permutations=n_permutations, seed=seed,
    )
    pull_mean = influence["pull_mean"]
    pull_p_value = influence["p_value"]
    pull_significant = (pull_mean > 0.0) and (pull_p_value < 0.01)  # 0.05 let chance-significant leaderless flocks through

    # === SCREENING: group steered toward preferred AND the informed minority
    #     demonstrably causes it (significant pull) ===
    screening_threshold = 0.3
    passes_screening = (mean_accuracy > screening_threshold) and pull_significant

    if not passes_screening:
        return DetectorResult(
            pattern_id="P19",
            detected=False,
            tier="none",
            confidence=0.0,
            primary_metric={
                "group_accuracy_mean": mean_accuracy,
                "group_accuracy_std": std_accuracy,
                "screening_threshold": screening_threshold,
            },
            secondary_metrics={"informed_fraction": informed_fraction,
                               "pull_mean": pull_mean, "pull_p_value": pull_p_value},
            effect_size={},
            null_p_value=pull_p_value,
            null_type="label_shuffle",
            exclusions_checked=[],
            exclusion_results={},
            co_occurrence_candidates=[],
            metadata_available=bool(metadata),
            warnings=warnings,
            notes=("Screening failed: no significant informed-minority influence "
                   "(accuracy=%.2f, pull_mean=%.3f, p=%.3f)" % (mean_accuracy, pull_mean, pull_p_value)),
        )

    # === CONFIRMATION: robust influence (stronger significance) ===
    passes_confirmation = (pull_p_value < 0.01) and (influence["pull_fraction"] > 0.5)

    if not passes_confirmation:
        return DetectorResult(
            pattern_id="P19",
            detected=True,
            tier="screening",
            confidence=0.50,
            primary_metric={
                "group_accuracy_mean": mean_accuracy,
                "group_accuracy_std": std_accuracy,
            },
            secondary_metrics={
                "informed_fraction": informed_fraction,
                "pull_mean": pull_mean,
                "pull_fraction": influence["pull_fraction"],
                "pull_p_value": pull_p_value,
            },
            effect_size={"accuracy_effect": mean_accuracy},
            null_p_value=pull_p_value,
            null_type="label_shuffle",
            exclusions_checked=[],
            exclusion_results={},
            co_occurrence_candidates=["P5"],
            metadata_available=bool(metadata),
            warnings=warnings,
            notes="Screening pass, confirmation fail: influence asymmetry not significant",
        )

    # === DEFINITIVE: disproportionate guidance efficacy ===
    is_minority = informed_fraction <= 0.25
    guidance_efficacy = (
        mean_accuracy / informed_fraction
        if informed_fraction > 0
        else 0.0
    )
    passes_definitive = (
        is_minority
        and guidance_efficacy > 2.0
        and pull_p_value < 0.01
    )

    if passes_definitive:
        tier = "definitive"
        confidence = 0.90
    else:
        tier = "confirmation"
        confidence = 0.70

    return DetectorResult(
        pattern_id="P19",
        detected=True,
        tier=tier,
        confidence=confidence,
        primary_metric={
            "group_accuracy_mean": mean_accuracy,
            "group_accuracy_std": std_accuracy,
        },
        secondary_metrics={
            "informed_fraction": informed_fraction,
            "n_informed": n_informed,
            "n_total": n_total,
            "pull_mean": pull_mean,
            "pull_fraction": influence["pull_fraction"],
            "pull_p_value": pull_p_value,
            "guidance_efficacy": guidance_efficacy,
            "is_minority": is_minority,
        },
        effect_size={
            "accuracy_effect": mean_accuracy,
            "pull_mean": pull_mean,
            "guidance_efficacy": guidance_efficacy,
        },
        null_p_value=pull_p_value,
        null_type="label_shuffle",
        exclusions_checked=["P5"],
        exclusion_results={
            "P5": "P19 distinguished from P5 by informed_mask + directional bias",
        },
        co_occurrence_candidates=["P5"],
        metadata_available=bool(metadata),
        warnings=warnings,
        notes=f"P19 {tier}: accuracy={mean_accuracy:.3f}, "
              f"pull={pull_mean:.4f} (p={pull_p_value:.4f}), "
              f"efficacy={guidance_efficacy:.1f}",
    )
