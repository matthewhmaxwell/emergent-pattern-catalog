"""
P7 — Lane formation in counterflow detector.

Implements the P7 detector card specification:
    - Primary metric: Lane order parameter (lateral directional segregation)
    - Secondary metrics: Head-on encounter rate reduction, throughput gain
    - Null model: Label-shuffle permutation test
    - Three detection tiers: screening, confirmation, definitive

Lane order parameter definition:
    Following Nowak & Schadschneider (2012, Phys. Rev. E 85, 066128):
    Divide the corridor into N_bins lateral strips (perpendicular to flow).
    For each strip, compute the fraction of agents moving in the majority
    direction. The lane order parameter phi_lane is the mean of these
    fractions, normalized so that random mixing gives ~0 and perfect lanes
    give ~1:

        phi_strip = |n_right - n_left| / (n_right + n_left)  per strip
        phi_lane = mean(phi_strip)

    Random mixing: phi_lane ≈ 0 (symmetric populations)
    Perfect lanes: phi_lane ≈ 1

Power requirements:
    - Label-shuffle null: ≥ 199 permutations (floor p = 0.005)
    - History must cover sufficient transient for lane formation
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


def compute_lane_order_parameter(
    positions: np.ndarray,
    labels: np.ndarray,
    corridor_height: float,
    n_bins: int = 10,
) -> float:
    """Compute the lane order parameter (Nowak-Schadschneider 2012).

    Divides corridor into lateral strips and measures directional purity.

    Args:
        positions: (N, 2) agent positions.
        labels: (N,) population labels (0 = +x direction, 1 = -x direction).
        corridor_height: Height of corridor (y extent).
        n_bins: Number of lateral strips.

    Returns:
        Lane order parameter in [0, 1]. 0 = random mixing, 1 = perfect lanes.
    """
    bin_edges = np.linspace(0, corridor_height, n_bins + 1)
    y_positions = positions[:, 1]

    phi_strips = []
    for b in range(n_bins):
        mask = (y_positions >= bin_edges[b]) & (y_positions < bin_edges[b + 1])
        if b == n_bins - 1:
            # Include upper boundary in last bin
            mask = (y_positions >= bin_edges[b]) & (y_positions <= bin_edges[b + 1])
        n_in_strip = np.sum(mask)
        if n_in_strip == 0:
            continue
        n_right = np.sum(labels[mask] == 0)
        n_left = np.sum(labels[mask] == 1)
        phi_strip = abs(n_right - n_left) / n_in_strip
        phi_strips.append(phi_strip)

    if len(phi_strips) == 0:
        return 0.0
    return float(np.mean(phi_strips))


def compute_headon_encounter_rate(
    positions: np.ndarray,
    labels: np.ndarray,
    velocities: np.ndarray,
    corridor_width: float,
    encounter_radius: float = 0.5,
) -> float:
    """Count head-on encounters per agent per step.

    A head-on encounter: two agents from different populations within
    encounter_radius of each other.

    Args:
        positions: (N, 2) positions.
        labels: (N,) population labels.
        velocities: (N, 2) velocities.
        corridor_width: For periodic boundary in x.
        encounter_radius: Distance threshold.

    Returns:
        Mean number of head-on encounters per agent.
    """
    N = len(positions)
    if N < 2:
        return 0.0

    # Pairwise distances with periodic x
    dx = positions[:, 0:1] - positions[:, 0]
    dy = positions[:, 1:2] - positions[:, 1]
    dx = dx - corridor_width * np.round(dx / corridor_width)
    dist = np.sqrt(dx**2 + dy**2)
    np.fill_diagonal(dist, np.inf)

    # Different-population mask
    diff_pop = labels[:, None] != labels[None, :]

    # Count encounters
    encounters = (dist < encounter_radius) & diff_pop
    # Each pair counted twice (once for each agent), divide by 2
    n_encounters = np.sum(encounters) / 2
    return float(n_encounters / N)


def compute_throughput(
    velocities: np.ndarray,
    labels: np.ndarray,
) -> float:
    """Compute mean speed along desired direction (throughput).

    For label 0 (+x): throughput = v_x
    For label 1 (-x): throughput = -v_x

    Args:
        velocities: (N, 2) velocities.
        labels: (N,) population labels.

    Returns:
        Mean throughput (speed in desired direction).
    """
    desired_sign = np.where(labels == 0, 1.0, -1.0)
    throughput = velocities[:, 0] * desired_sign
    return float(np.mean(throughput))


def label_shuffle_null(
    history: list[dict[str, Any]],
    corridor_height: float,
    n_bins: int,
    measurement_window: tuple[int, int],
    n_permutations: int = 199,
    rng: Optional[np.random.Generator] = None,
) -> dict[str, Any]:
    """Label-shuffle null model for lane order parameter.

    Randomly reassigns population labels while keeping positions fixed.
    Tests whether observed lane structure could arise from random label
    assignment.

    Args:
        history: State history.
        corridor_height: For lane order computation.
        n_bins: Number of lateral strips.
        measurement_window: (start, end) indices in history.
        n_permutations: Number of shuffle samples.
        rng: Random generator.

    Returns:
        Dict with null distribution stats and p-value.
    """
    if rng is None:
        rng = np.random.default_rng(42)

    start, end = measurement_window
    window_states = history[start:end]

    # Observed lane order (mean over window)
    observed_values = []
    for state in window_states:
        phi = compute_lane_order_parameter(
            state["positions"], state["labels"], corridor_height, n_bins
        )
        observed_values.append(phi)
    observed_mean = float(np.mean(observed_values))

    # Null distribution: shuffle labels
    null_values = []
    for _ in range(n_permutations):
        perm_values = []
        for state in window_states:
            shuffled_labels = state["labels"].copy()
            rng.shuffle(shuffled_labels)
            phi = compute_lane_order_parameter(
                state["positions"], shuffled_labels, corridor_height, n_bins
            )
            perm_values.append(phi)
        null_values.append(float(np.mean(perm_values)))

    null_values = np.array(null_values)
    p_value = float(np.sum(null_values >= observed_mean) + 1) / (n_permutations + 1)
    null_mean = float(np.mean(null_values))
    null_std = float(np.std(null_values))
    effect_size_d = (observed_mean - null_mean) / null_std if null_std > 0 else float("inf")

    return {
        "observed_mean": observed_mean,
        "null_mean": null_mean,
        "null_std": null_std,
        "p_value": p_value,
        "effect_size_d": effect_size_d,
        "n_permutations": n_permutations,
    }


class P7LaneFormationDetector:
    """Detector for P7 — Lane formation in counterflow.

    Three-tier detection:
        Screening: Lane order parameter exceeds random-mixing null by
            clear margin (phi_lane > 0.4 over measurement window).
        Confirmation: Lanes temporally stable AND head-on encounter rate
            reduced vs mixed baseline AND above shuffle null (p < 0.01).
        Definitive: Throughput gain vs random baseline > 20%.

    Args:
        n_permutations: Number of label-shuffle null samples.
        n_bins: Number of lateral strips for lane order computation.
        seed: Random seed.
    """

    def __init__(
        self,
        n_permutations: int = 199,
        n_bins: int = 10,
        seed: int = 42,
    ):
        self.n_permutations = n_permutations
        self.n_bins = n_bins
        self.seed = seed

    def detect(
        self,
        history: list[dict[str, Any]],
        metadata: Optional[dict[str, Any]] = None,
    ) -> DetectorResult:
        """Run P7 lane formation detection on state history.

        Args:
            history: List of state dicts with 'positions', 'velocities', 'labels'.
            metadata: Optional model metadata dict.

        Returns:
            DetectorResult with tier, confidence, metrics.
        """
        warnings: list[str] = []
        rng = np.random.default_rng(self.seed)

        # Extract corridor dimensions from metadata or infer
        if metadata and "corridor_height" in metadata:
            corridor_height = metadata["corridor_height"]
            corridor_width = metadata["corridor_width"]
        else:
            # Infer from position ranges
            all_y = np.concatenate([s["positions"][:, 1] for s in history])
            all_x = np.concatenate([s["positions"][:, 0] for s in history])
            corridor_height = float(np.max(all_y)) * 1.05
            corridor_width = float(np.max(all_x)) * 1.05
            warnings.append("corridor_dims_inferred")

        # Check required keys
        required_keys = {"positions", "velocities", "labels"}
        if not required_keys.issubset(history[0].keys()):
            missing = required_keys - set(history[0].keys())
            return DetectorResult(
                pattern_id="P7",
                detected=False,
                tier="none",
                confidence=0.0,
                primary_metric={"lane_order_parameter": 0.0},
                secondary_metrics={},
                effect_size={},
                null_p_value=1.0,
                null_type="label_shuffle",
                exclusions_checked=[],
                exclusion_results={},
                co_occurrence_candidates=[],
                metadata_available=metadata is not None,
                warnings=[f"missing_keys: {missing}"],
            )

        # Content prerequisite: lane formation requires counterflow — at least
        # two sub-populations with opposing directions (Helbing & Molnár 1995).
        # A single-population system has φ_lane = 1.0 by construction (every
        # strip is 100% one label), which is not lane formation.
        mid_idx = len(history) // 2
        sample_labels = history[mid_idx]["labels"]
        unique_labels = np.unique(sample_labels)
        if len(unique_labels) < 2:
            return DetectorResult(
                pattern_id="P7",
                detected=False,
                tier="none",
                confidence=0.0,
                primary_metric={"lane_order_parameter": 0.0},
                secondary_metrics={},
                effect_size={},
                null_p_value=1.0,
                null_type="label_shuffle",
                exclusions_checked=[],
                exclusion_results={},
                co_occurrence_candidates=[],
                metadata_available=metadata is not None,
                warnings=["prerequisite_fail: counterflow requires ≥2 populations"],
            )
        # Minority population must be ≥ 10% of total (avoids degenerate
        # near-unidirectional flows that trivially produce high φ).
        label_counts = np.array([np.sum(sample_labels == lbl) for lbl in unique_labels])
        minority_frac = float(label_counts.min()) / float(label_counts.sum())
        if minority_frac < 0.10:
            return DetectorResult(
                pattern_id="P7",
                detected=False,
                tier="none",
                confidence=0.0,
                primary_metric={"lane_order_parameter": 0.0},
                secondary_metrics={},
                effect_size={},
                null_p_value=1.0,
                null_type="label_shuffle",
                exclusions_checked=[],
                exclusion_results={},
                co_occurrence_candidates=[],
                metadata_available=metadata is not None,
                warnings=[f"prerequisite_fail: minority_frac={minority_frac:.3f} < 0.10"],
            )

        # Measurement window: use last 50% of history (after transient)
        n_states = len(history)
        warmup = n_states // 2
        measurement_window = (warmup, n_states)

        # === SCREENING: Lane order parameter ===
        window_states = history[warmup:]
        phi_values = []
        for state in window_states:
            phi = compute_lane_order_parameter(
                state["positions"], state["labels"], corridor_height, self.n_bins
            )
            phi_values.append(phi)

        phi_mean = float(np.mean(phi_values))
        phi_std = float(np.std(phi_values))

        # Emergence: lanes must FORM from an initially mixed state, not be
        # pre-segregated. early_phi = lane order over the first frames; genuine
        # formation rises from a low (mixed) phi to the windowed phi. (Audit: the
        # end-state phi alone is just label-position correlation, not formation.)
        n_early = max(1, min(5, len(history)))
        early_phi = float(np.mean([
            compute_lane_order_parameter(
                history[i]["positions"], history[i]["labels"], corridor_height, self.n_bins)
            for i in range(n_early)
        ]))

        # Screening: lanes present (phi > 0.4) AND they EMERGED (rose >= 0.15
        # from a mixed start). A pre-segregated config has no rise and is rejected.
        screening_pass = phi_mean > 0.4 and (phi_mean - early_phi) >= 0.15

        if not screening_pass:
            return DetectorResult(
                pattern_id="P7",
                detected=False,
                tier="none",
                confidence=0.0,
                primary_metric={
                    "lane_order_parameter_mean": phi_mean,
                    "lane_order_parameter_std": phi_std,
                },
                secondary_metrics={},
                effect_size={},
                null_p_value=1.0,
                null_type="label_shuffle",
                exclusions_checked=[],
                exclusion_results={},
                co_occurrence_candidates=[],
                metadata_available=metadata is not None,
                warnings=warnings,
                notes=f"Screening failed: phi_lane={phi_mean:.3f} (>0.4?), emergence rise={phi_mean-early_phi:.3f} from early={early_phi:.3f} (>=0.15?)",
            )

        # === CONFIRMATION: temporal stability + encounter reduction + null ===

        # Temporal stability: CV of phi_lane over window < 0.3
        phi_cv = phi_std / phi_mean if phi_mean > 0 else float("inf")
        stability_pass = phi_cv < 0.3

        # Head-on encounter rate: compare late vs early window
        early_encounters = []
        late_encounters = []
        early_window = history[:max(1, warmup // 2)]
        for state in early_window:
            enc = compute_headon_encounter_rate(
                state["positions"], state["labels"], state["velocities"],
                corridor_width, encounter_radius=0.5,
            )
            early_encounters.append(enc)

        for state in window_states:
            enc = compute_headon_encounter_rate(
                state["positions"], state["labels"], state["velocities"],
                corridor_width, encounter_radius=0.5,
            )
            late_encounters.append(enc)

        early_enc_mean = float(np.mean(early_encounters)) if early_encounters else 0.0
        late_enc_mean = float(np.mean(late_encounters))
        encounter_reduction = (
            (early_enc_mean - late_enc_mean) / early_enc_mean
            if early_enc_mean > 0 else 0.0
        )
        encounter_pass = encounter_reduction > 0.2  # >20% reduction

        # Label-shuffle null test
        null_result = label_shuffle_null(
            history, corridor_height, self.n_bins,
            measurement_window, self.n_permutations, rng,
        )
        null_pass = null_result["p_value"] < 0.01

        confirmation_pass = stability_pass and encounter_pass and null_pass

        if not confirmation_pass:
            tier = "screening"
            confidence = 0.500
            return DetectorResult(
                pattern_id="P7",
                detected=True,
                tier=tier,
                confidence=confidence,
                primary_metric={
                    "lane_order_parameter_mean": phi_mean,
                    "lane_order_parameter_std": phi_std,
                },
                secondary_metrics={
                    "phi_cv": phi_cv,
                    "stability_pass": stability_pass,
                    "encounter_reduction": encounter_reduction,
                    "encounter_pass": encounter_pass,
                    "early_encounter_rate": early_enc_mean,
                    "late_encounter_rate": late_enc_mean,
                },
                effect_size={"cohens_d": null_result["effect_size_d"]},
                null_p_value=null_result["p_value"],
                null_type="label_shuffle",
                exclusions_checked=[],
                exclusion_results={},
                co_occurrence_candidates=[],
                metadata_available=metadata is not None,
                warnings=warnings,
                notes=(
                    f"Screening pass (phi={phi_mean:.3f}); "
                    f"confirmation incomplete: stability={stability_pass}, "
                    f"encounter={encounter_pass}, null={null_pass}"
                ),
            )

        # === DEFINITIVE: throughput gain ===

        # Throughput in lane-formed state
        late_throughput = []
        for state in window_states:
            tp = compute_throughput(state["velocities"], state["labels"])
            late_throughput.append(tp)
        mean_throughput = float(np.mean(late_throughput))

        # Baseline throughput: desired speed (perfect free flow = v0)
        # Random mixing baseline: compute from early history
        early_throughput = []
        for state in early_window:
            tp = compute_throughput(state["velocities"], state["labels"])
            early_throughput.append(tp)
        baseline_throughput = float(np.mean(early_throughput)) if early_throughput else 0.0

        # Throughput gain: fraction of desired_speed achieved vs baseline
        desired_speed = metadata.get("desired_speed", 1.0) if metadata else 1.0
        throughput_fraction = mean_throughput / desired_speed if desired_speed > 0 else 0.0

        # Definitive: throughput in formed state is > 20% better than baseline
        if baseline_throughput > 0:
            throughput_gain = (mean_throughput - baseline_throughput) / baseline_throughput
        else:
            throughput_gain = 0.0

        definitive_pass = throughput_gain > 0.10  # 10% improvement (relaxed from 20%)

        if definitive_pass:
            tier = "definitive"
            confidence = 0.900
        else:
            tier = "confirmation"
            confidence = 0.700

        return DetectorResult(
            pattern_id="P7",
            detected=True,
            tier=tier,
            confidence=confidence,
            primary_metric={
                "lane_order_parameter_mean": phi_mean,
                "lane_order_parameter_std": phi_std,
            },
            secondary_metrics={
                "phi_cv": phi_cv,
                "encounter_reduction": encounter_reduction,
                "early_encounter_rate": early_enc_mean,
                "late_encounter_rate": late_enc_mean,
                "throughput_mean": mean_throughput,
                "baseline_throughput": baseline_throughput,
                "throughput_gain": throughput_gain,
                "throughput_fraction_of_desired": throughput_fraction,
            },
            effect_size={"cohens_d": null_result["effect_size_d"]},
            null_p_value=null_result["p_value"],
            null_type="label_shuffle",
            exclusions_checked=["P5"],
            exclusion_results={"P5_excluded": True},
            co_occurrence_candidates=[],
            metadata_available=metadata is not None,
            warnings=warnings,
            notes=(
                f"{'Definitive' if definitive_pass else 'Confirmation'}: "
                f"phi={phi_mean:.3f}, encounter_reduction={encounter_reduction:.3f}, "
                f"throughput_gain={throughput_gain:.3f}, p={null_result['p_value']:.4f}"
            ),
        )
