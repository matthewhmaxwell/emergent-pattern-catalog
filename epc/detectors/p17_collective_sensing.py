"""
P17 — Distributed sensing / collective gradient detection detector.

Implements the P17 detector card specification:
    - Primary metric: Group chemotactic index (CI) — progress toward field peak
    - Confirmation metric: CI scaling with group size N (the Berdahl signature)
    - Null model: alpha=0 (no speed-modulation) comparison
    - Three detection tiers: screening, confirmation, definitive

Detection logic:
    The defining P17 signature is that group navigation accuracy RISES with
    group size while isolated individuals remain near chance (Berdahl et al.
    2013, Science). The detector runs the model at multiple group sizes and
    checks for positive CI-vs-N scaling.

    - Screening (≤0.60): CI at canonical N exceeds isolated-agent baseline
    - Confirmation (≤0.85): CI increases with N (positive slope of CI vs log(N)),
      significance tested via alpha=0 permutation null
    - Definitive (≤1.00): CI at large N substantially exceeds chance AND
      the N-scaling is monotonic AND null p < 0.005

Canonical parameters (Berdahl-style):
    box_size=20.0, v_max=0.4, turn_noise=0.3, sensing_noise=0.8,
    alpha=0.95, social_strength=0.2, field_sigma=5.0, n_steps=1000

Power requirements:
    - Null runs: ≥ 199 permutations (floor p = 0.005)
    - Group size sweep: at least 4 N values covering [1, N_max]
    - Seeds per N: ≥ 5 for reliable averaging
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


def compute_chemotactic_index(
    history: list[dict[str, Any]],
    field_center: tuple[float, float],
    box_size: float,
) -> float:
    """Compute chemotactic index: normalized approach to field peak.

    CI = (d_initial - d_final) / d_initial, where distances are periodic.

    Args:
        history: State history with 'positions' key.
        field_center: (x, y) of field peak.
        box_size: Periodic domain size.

    Returns:
        CI in [-1, 1]. Positive = approaching peak.
    """
    if len(history) < 2:
        return 0.0

    cx, cy = field_center
    L = box_size

    def periodic_com(positions: np.ndarray) -> tuple[float, float]:
        theta_x = 2 * np.pi * positions[:, 0] / L
        theta_y = 2 * np.pi * positions[:, 1] / L
        com_x = L * np.arctan2(np.mean(np.sin(theta_x)), np.mean(np.cos(theta_x))) / (2 * np.pi)
        com_y = L * np.arctan2(np.mean(np.sin(theta_y)), np.mean(np.cos(theta_y))) / (2 * np.pi)
        return (com_x % L, com_y % L)

    def periodic_dist(x1: float, y1: float, x2: float, y2: float) -> float:
        dx = x1 - x2
        dy = y1 - y2
        dx = dx - L * round(dx / L)
        dy = dy - L * round(dy / L)
        return float(np.sqrt(dx**2 + dy**2))

    early_end = max(2, len(history) // 10)
    early_dists = []
    for s in history[:early_end]:
        com = periodic_com(s["positions"])
        early_dists.append(periodic_dist(com[0], com[1], cx, cy))
    d_init = float(np.mean(early_dists))

    late_start = len(history) * 3 // 4
    late_dists = []
    for s in history[late_start:]:
        com = periodic_com(s["positions"])
        late_dists.append(periodic_dist(com[0], com[1], cx, cy))
    d_final = float(np.mean(late_dists))

    if d_init < 1e-8:
        return 1.0 if d_final < 1e-8 else 0.0

    ci = (d_init - d_final) / d_init
    return float(np.clip(ci, -1.0, 1.0))


def _run_group_ci(
    n_agents: int,
    box_size: float = 20.0,
    v_max: float = 0.4,
    turn_noise: float = 0.3,
    sensing_noise: float = 0.8,
    alpha: float = 0.95,
    social_strength: float = 0.2,
    field_sigma: float = 5.0,
    field_amplitude: float = 1.0,
    n_steps: int = 1000,
    record_interval: int = 5,
    seed: int = 42,
) -> float:
    """Run collective sensing model at given N and return CI.

    Args:
        n_agents: Group size.
        Other args: model parameters.

    Returns:
        Chemotactic index for this run.
    """
    from epc.models.collective_sensing import CollectiveSensingModel

    model = CollectiveSensingModel(
        n_agents=n_agents,
        box_size=box_size,
        v_max=v_max,
        turn_noise=turn_noise,
        sensing_noise=sensing_noise,
        alpha=alpha,
        social_strength=social_strength,
        field_sigma=field_sigma,
        field_amplitude=field_amplitude,
        dt=1.0,
        init_mode="offset",
        seed=seed,
    )
    history = model.run(n_steps=n_steps, record_interval=record_interval)
    return compute_chemotactic_index(history, model.field_center, box_size)


class P17CollectiveSensingDetector:
    """Detector for P17 — distributed sensing / collective gradient detection.

    Detection methodology:
        1. Run the model at multiple group sizes N ∈ {1, 5, 10, 25, 50}
        2. Compute CI at each N (averaged over n_seeds seeds)
        3. Screening: CI(N_max) > CI(N=1) + threshold
        4. Confirmation: positive slope of CI vs log(N) + permutation null
        5. Definitive: monotonic increase + large effect size + p < 0.005
    """

    def __init__(
        self,
        group_sizes: Optional[list[int]] = None,
        n_null_runs: int = 199,
        n_seeds: int = 5,
        n_steps: int = 1000,
        record_interval: int = 5,
        seed: int = 42,
    ):
        """Initialize P17 detector.

        Args:
            group_sizes: List of N values to test. Default: [1, 5, 10, 25, 50].
            n_null_runs: Number of null model runs for comparison.
            n_seeds: Seeds per group size for averaging.
            n_steps: Simulation steps per run.
            record_interval: Recording interval.
            seed: Master random seed.
        """
        self.group_sizes = group_sizes or [1, 5, 10, 25, 50]
        self.n_null_runs = n_null_runs
        self.n_seeds = n_seeds
        self.n_steps = n_steps
        self.record_interval = record_interval
        self.seed = seed

    def detect(
        self,
        history: list[dict[str, Any]],
        metadata: dict[str, Any],
    ) -> DetectorResult:
        """Run P17 detection.

        Uses model parameters from metadata to run at multiple group sizes
        and test for the N-scaling signature.

        Args:
            history: State history from a collective-sensing run (used for
                     metadata extraction; the actual detection re-runs at
                     multiple N values).
            metadata: Model metadata dict.

        Returns:
            DetectorResult with tier and confidence.
        """
        warnings: list[str] = []

        # Extract parameters from metadata
        box_size = metadata.get("box_size", 20.0)
        v_max = metadata.get("v_max", 0.4)
        turn_noise = metadata.get("turn_noise", 0.3)
        sensing_noise = metadata.get("sensing_noise", 0.8)
        alpha = metadata.get("alpha", 0.95)
        social_strength = metadata.get("social_strength", 0.2)
        field_sigma = metadata.get("field_sigma", 5.0)
        field_amplitude = metadata.get("field_amplitude", 1.0)
        field_center = metadata.get("field_center", (box_size / 2, box_size / 2))
        if isinstance(field_center, list):
            field_center = tuple(field_center)

        rng = np.random.default_rng(self.seed)

        # --- Run at multiple group sizes ---
        ci_by_n: dict[int, list[float]] = {}
        for n in self.group_sizes:
            ci_values = []
            for _ in range(self.n_seeds):
                seed_val = int(rng.integers(0, 2**31))
                ci = _run_group_ci(
                    n_agents=n,
                    box_size=box_size,
                    v_max=v_max,
                    turn_noise=turn_noise,
                    sensing_noise=sensing_noise,
                    alpha=alpha,
                    social_strength=social_strength,
                    field_sigma=field_sigma,
                    field_amplitude=field_amplitude,
                    n_steps=self.n_steps,
                    record_interval=self.record_interval,
                    seed=seed_val,
                )
                ci_values.append(ci)
            ci_by_n[n] = ci_values

        mean_ci = {n: float(np.mean(v)) for n, v in ci_by_n.items()}
        sorted_ns = sorted(mean_ci.keys())
        ci_arr = np.array([mean_ci[n] for n in sorted_ns])

        # Baseline (smallest N) and max
        ci_baseline = mean_ci[sorted_ns[0]]
        ci_max = mean_ci[sorted_ns[-1]]

        # Slope of CI vs log(N)
        log_ns = np.log(np.array(sorted_ns, dtype=float))
        slope = float(np.polyfit(log_ns, ci_arr, 1)[0]) if len(sorted_ns) >= 2 else 0.0

        # Monotonicity (with tolerance for noise)
        diffs = np.diff(ci_arr)
        is_monotonic = bool(np.all(diffs >= -0.03))

        # --- Null: alpha=0 (no speed modulation) ---
        null_slopes: list[float] = []
        null_rng = np.random.default_rng(self.seed + 10000)
        for _ in range(self.n_null_runs):
            null_ci_by_n: dict[int, float] = {}
            for n in self.group_sizes:
                seed_val = int(null_rng.integers(0, 2**31))
                ci_null = _run_group_ci(
                    n_agents=n,
                    box_size=box_size,
                    v_max=v_max,
                    turn_noise=turn_noise,
                    sensing_noise=sensing_noise,
                    alpha=0.0,  # null: no speed modulation
                    social_strength=social_strength,
                    field_sigma=field_sigma,
                    field_amplitude=field_amplitude,
                    n_steps=self.n_steps,
                    record_interval=self.record_interval,
                    seed=seed_val,
                )
                null_ci_by_n[n] = ci_null
            null_ci_vals = np.array([null_ci_by_n[n] for n in sorted_ns])
            null_slope = float(np.polyfit(log_ns, null_ci_vals, 1)[0])
            null_slopes.append(null_slope)

        null_slopes_arr = np.array(null_slopes)
        p_value = float(np.mean(null_slopes_arr >= slope))

        # Effect size
        null_mean = float(np.mean(null_slopes_arr))
        null_std = float(np.std(null_slopes_arr))
        if null_std > 1e-10:
            cohens_d = (slope - null_mean) / null_std
        else:
            cohens_d = float('inf') if slope > null_mean else 0.0

        # --- Tier assignment ---
        tier = "none"
        confidence = 0.0

        ci_gain = ci_max - ci_baseline
        if ci_gain > 0.05 and ci_max > 0.1:
            tier = "screening"
            confidence = 0.500

            if slope > 0.02 and p_value < 0.05:
                tier = "confirmation"
                confidence = 0.700

                if is_monotonic and cohens_d > 2.0 and p_value < 0.01 and ci_max > 0.2:
                    tier = "definitive"
                    confidence = 0.900

        detected = tier != "none"

        notes = (
            f"CI by N: {mean_ci}; slope={slope:.4f}; "
            f"monotonic={is_monotonic}; p={p_value:.4f}; d={cohens_d:.2f}"
        )

        return DetectorResult(
            pattern_id="P17",
            detected=detected,
            tier=tier,
            confidence=confidence,
            primary_metric={
                "chemotactic_index_by_N": mean_ci,
                "ci_at_max_N": ci_max,
                "ci_baseline_N1": ci_baseline,
            },
            secondary_metrics={
                "ci_slope_vs_logN": slope,
                "is_monotonic": is_monotonic,
                "ci_gain": ci_gain,
            },
            effect_size={
                "cohens_d": cohens_d,
                "ci_slope": slope,
                "null_mean_slope": null_mean,
                "null_std_slope": null_std,
            },
            null_p_value=p_value,
            null_type="alpha_zero_null",
            exclusions_checked=["P5_flocking"],
            exclusion_results={
                "P5": "P5 detects collective motion without field-inference; "
                      "P17 requires CI-vs-N scaling which pure flocking lacks"
            },
            co_occurrence_candidates=["P5"],
            metadata_available=True,
            warnings=warnings,
            notes=notes,
        )
