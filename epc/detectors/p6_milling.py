"""
P6 — Milling / vortex formation detector.

Implements the full P6 detector card specification:
    - Primary metric: Angular momentum |L|
    - Secondary metrics: Group speed ratio R, ring density profile,
      rotation coherence
    - Null model: Heading-shuffle (random headings)
    - Exclusions: P5 (R > 0.5 → flocking)
    - Three detection tiers: screening, confirmation, definitive

Detector card reference:
    Screening: |L| > 0.3 over ≥ 5 × T_cross
    Confirmation: |L| > 0.5 AND R < 0.3 over ≥ 10 × T_cross
    Definitive: Confirmation + ring density confirmed + P5 exclusion
"""

from __future__ import annotations

import numpy as np
from dataclasses import dataclass, field
from typing import Any, Optional

from epc.metrics.collective_motion import (
    PolarizationMetric,
    GroupSpeedRatioMetric,
    AngularMomentumMetric,
)


@dataclass
class DetectorResult:
    """Standardized detector output."""
    pattern_id: str
    detected: bool
    tier: str
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


def ring_density_profile(
    state: dict[str, Any],
    n_bins: int = 10,
) -> dict[str, Any]:
    """Compute radial density profile from center of mass.

    Returns hollowness = ρ(0) / ρ(r*) where r* is the peak radius.
    Hollowness < 0.5 indicates a ring (hollow core).

    Args:
        state: State dict with 'positions' key.
        n_bins: Number of radial bins.

    Returns:
        Dict with radial profile and hollowness.
    """
    positions = state["positions"]
    com = np.mean(positions, axis=0)
    radii = np.linalg.norm(positions - com, axis=1)

    max_r = np.max(radii) * 1.1
    bin_edges = np.linspace(0, max_r, n_bins + 1)
    bin_centers = 0.5 * (bin_edges[:-1] + bin_edges[1:])
    bin_widths = np.diff(bin_edges)

    counts, _ = np.histogram(radii, bins=bin_edges)
    # Normalize by annular area: 2πr·dr
    areas = 2 * np.pi * bin_centers * bin_widths
    areas = np.maximum(areas, 1e-12)
    density = counts / areas

    # Hollowness: density at center / density at peak
    peak_idx = np.argmax(density)
    peak_density = density[peak_idx]
    core_density = density[0] if len(density) > 0 else 0

    if peak_density > 0:
        hollowness = core_density / peak_density
    else:
        hollowness = 1.0  # no clear structure

    return {
        "bin_centers": bin_centers,
        "density": density,
        "peak_radius": float(bin_centers[peak_idx]),
        "hollowness": float(hollowness),
        "mean_radius": float(np.mean(radii)),
        "std_radius": float(np.std(radii)),
    }


def milling_null(
    history: list[dict[str, Any]],
    n_permutations: int = 199,
    box_size: Optional[float] = None,
    rng: Optional[np.random.Generator] = None,
) -> dict[str, Any]:
    """Random-heading null model for angular momentum.

    Replaces all headings with uniform random directions.
    Under this null, |L| ≈ 0 (no rotational coherence).

    Args:
        history: State history.
        n_permutations: Monte Carlo samples.
        box_size: For periodic-aware COM.
        rng: Random generator.

    Returns:
        Dict with null distribution and p-value.
    """
    if rng is None:
        rng = np.random.default_rng(42)

    # Observed |L| time series
    observed_L = np.array([
        abs(AngularMomentumMetric.compute_instant(s, box_size))
        for s in history
    ])
    observed_L_mean = float(np.mean(observed_L))

    N = len(history[0]["velocities"])
    T = len(history)

    null_L_values = np.empty(n_permutations)

    for perm in range(n_permutations):
        perm_L = []
        for state in history:
            # Randomize heading directions, keep positions and speeds
            pos = state["positions"]
            vel = state["velocities"]
            speeds = np.linalg.norm(vel, axis=1)
            random_headings = rng.uniform(-np.pi, np.pi, size=N)
            random_vel = np.column_stack([
                speeds * np.cos(random_headings),
                speeds * np.sin(random_headings),
            ])
            fake_state = {"positions": pos, "velocities": random_vel}
            perm_L.append(abs(
                AngularMomentumMetric.compute_instant(fake_state, box_size)
            ))
        null_L_values[perm] = np.mean(perm_L)

    p_value = (np.sum(null_L_values >= observed_L_mean) + 1) / (n_permutations + 1)
    null_mean = float(np.mean(null_L_values))
    null_std = float(np.std(null_L_values))

    if null_std > 1e-12:
        cohens_d = (observed_L_mean - null_mean) / null_std
    else:
        cohens_d = float("inf") if observed_L_mean > null_mean else 0.0

    return {
        "null_L_values": null_L_values,
        "null_mean": null_mean,
        "null_std": null_std,
        "observed_L_abs_mean": observed_L_mean,
        "p_value": float(p_value),
        "cohens_d": float(cohens_d),
        "n_permutations": n_permutations,
    }


class P6MillingDetector:
    """Detector for P6 — Milling / vortex formation.

    Detection tiers:
        Screening: |L| > 0.3 over ≥ 5 T_cross
        Confirmation: |L| > 0.5 AND R < 0.3 over ≥ 10 T_cross
        Definitive: Confirmation + ring density (hollowness < 0.5) + P5 exclusion

    Args:
        n_permutations: Permutation count for null test.
        T_cross: Domain crossing time in timesteps.
        box_size: For periodic-aware COM (None for open space).
        seed: Random seed.
    """

    SCREENING_L = 0.3
    SCREENING_T_CROSS_MIN = 5

    CONFIRMATION_L = 0.5
    CONFIRMATION_R_MAX = 0.3
    CONFIRMATION_T_CROSS_MIN = 10

    DEFINITIVE_HOLLOWNESS_MAX = 0.5
    P5_EXCLUSION_R_THRESHOLD = 0.5

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
        """Run full P6 detection pipeline."""
        warnings = []
        notes_parts = []

        # --- T_cross ---
        T_cross = self.T_cross if self.T_cross else 100.0
        if metadata and "v_eq" in metadata and "init_radius" in metadata:
            T_cross = 2 * metadata["init_radius"] / metadata["v_eq"]

        T = len(history)
        T_cross_steps = max(T_cross, 1.0)

        # Measurement window: second half of trajectory
        meas_start = T // 2
        meas_end = T
        meas_length = meas_end - meas_start
        meas_history = history[meas_start:meas_end]

        # Compute actual time covered by measurement window
        # Use step numbers from state dicts if available
        dt_per_state = 1.0  # default: each history entry = 1 unit
        if metadata and "dt" in metadata and len(meas_history) >= 2:
            step_first = meas_history[0].get("step", 0)
            step_last = meas_history[-1].get("step", len(meas_history) - 1)
            step_span = max(step_last - step_first, 1)
            dt_per_state = step_span * metadata["dt"] / max(len(meas_history) - 1, 1)
            total_time = dt_per_state * (len(meas_history) - 1)
        elif metadata and "dt" in metadata:
            total_time = meas_length * metadata["dt"]
        else:
            total_time = float(meas_length)

        n_T_cross = total_time / T_cross if T_cross > 0 else 0

        notes_parts.append(
            f"Measurement window: steps {meas_start}–{meas_end} "
            f"({n_T_cross:.1f} T_cross)"
        )

        # --- Primary metric: |L| ---
        L_result = AngularMomentumMetric.compute(meas_history, self.box_size)
        L_abs_mean = L_result["L_abs_mean"]

        # --- Secondary metrics ---
        R_result = GroupSpeedRatioMetric.compute(meas_history)
        R_mean = R_result["R_mean"]

        phi_result = PolarizationMetric.compute(meas_history)
        phi_mean = phi_result["phi_mean"]

        # Ring density from last state
        ring = ring_density_profile(meas_history[-1])
        hollowness = ring["hollowness"]

        # --- Null model ---
        # Subsample history for null (speed optimization)
        subsample = max(1, len(meas_history) // 50)
        null_history = meas_history[::subsample]
        null_result = milling_null(
            null_history,
            n_permutations=self.n_permutations,
            box_size=self.box_size,
            rng=self.rng,
        )
        null_p = null_result["p_value"]
        cohens_d = null_result["cohens_d"]

        # --- Exclusion checks ---
        exclusion_results = {}

        # P5: R > 0.5 → flocking
        p5_excluded = R_mean > self.P5_EXCLUSION_R_THRESHOLD
        exclusion_results["P5"] = "TRIGGERED" if p5_excluded else "excluded"

        # --- Tier determination ---
        screening_pass = (
            L_abs_mean > self.SCREENING_L
            and n_T_cross >= self.SCREENING_T_CROSS_MIN
        )

        confirmation_pass = (
            L_abs_mean > self.CONFIRMATION_L
            and R_mean < self.CONFIRMATION_R_MAX
            and n_T_cross >= self.CONFIRMATION_T_CROSS_MIN
        )

        definitive_pass = (
            confirmation_pass
            and hollowness < self.DEFINITIVE_HOLLOWNESS_MAX
            and not p5_excluded
        )

        tier = "none"
        detected = False
        confidence = 0.0

        if definitive_pass:
            tier = "definitive"
            detected = True
            confidence = 0.75
            if not p5_excluded:
                confidence += 0.10
            if null_p < 0.001:
                confidence += 0.10
            confidence = min(confidence, 1.0)

        elif confirmation_pass:
            tier = "confirmation"
            detected = True
            confidence = 0.55
            if null_p < 0.001:
                confidence += 0.15
            if cohens_d > 1.0:
                confidence += 0.10
            confidence = min(confidence, 0.85)

        elif screening_pass:
            tier = "screening"
            detected = True
            confidence = 0.35
            if null_p < 0.01:
                confidence += 0.10
            confidence = min(confidence, 0.60)

        return DetectorResult(
            pattern_id="P6",
            detected=detected,
            tier=tier,
            confidence=round(confidence, 3),
            primary_metric={"angular_momentum_abs_mean": L_abs_mean},
            secondary_metrics={
                "group_speed_ratio": R_mean,
                "polarization": phi_mean,
                "hollowness": hollowness,
                "peak_radius": ring["peak_radius"],
                "mean_radius": ring["mean_radius"],
                "measurement_T_cross": n_T_cross,
            },
            effect_size={
                "cohens_d": cohens_d,
                "L_observed": L_abs_mean,
                "L_null_mean": null_result["null_mean"],
            },
            null_p_value=null_p,
            null_type="shuffle",
            exclusions_checked=["P5"],
            exclusion_results=exclusion_results,
            co_occurrence_candidates=["P9"],
            metadata_available=metadata is not None,
            warnings=warnings,
            notes=" | ".join(notes_parts),
        )
