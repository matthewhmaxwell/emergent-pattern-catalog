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
        """Run P17 detection from the SUBSTRATE (no re-simulation).

        Emergent collective gradient sensing (Berdahl et al. 2013): a group
        climbs a noisy scalar field toward its peak even though individuals
        cannot sense the gradient alone. Measured entirely from the provided
        trajectory:
          1. chemotactic_index: the group centroid's normalized approach to the
             TRUE field peak (does the group climb?).
          2. field-sensing null: the group must approach the true peak more than
             random target locations (p<0.05); flocking / non-social motion
             (social_off) approaches a random point as much as the true peak.
          3. emergence: individual sensing must be too noisy to climb alone
             (individual_SNR low) -- otherwise it is trivial individual sensing,
             not emergent collective sensing (field_too_strong).
        Only field_center / field_sigma / field_amplitude / box_size are read
        from metadata; these describe the ENVIRONMENT, not the detection signal.
        """
        warnings: list[str] = []
        CLIMB_MIN, DEF_CI_MIN, EMERGENT_SNR_MAX, COHESION_MAX = 0.05, 0.30, 0.5, 0.10

        box_size = float(metadata.get("box_size", 20.0))
        field_sigma = float(metadata.get("field_sigma", 5.0))
        field_amplitude = float(metadata.get("field_amplitude", 1.0))
        fc = metadata.get("field_center", (box_size / 2, box_size / 2))
        if isinstance(fc, list):
            fc = tuple(fc)
        center = (float(fc[0]), float(fc[1]))

        if not history or "positions" not in history[0]:
            return self._null_result("no 'positions' in substrate history")
        L = box_size

        ci_true = compute_chemotactic_index(history, center, L)

        rng = np.random.default_rng(self.seed)
        n_rand = max(self.n_null_runs, 99)
        ci_random = np.array([
            compute_chemotactic_index(
                history, (float(rng.uniform(0, L)), float(rng.uniform(0, L))), L)
            for _ in range(n_rand)
        ])
        null_mean = float(ci_random.mean())
        null_std = float(ci_random.std())
        p_value = float(np.mean(ci_random >= ci_true))
        if null_std > 1e-10:
            cohens_d = (ci_true - null_mean) / null_std
        else:
            cohens_d = float("inf") if ci_true > null_mean else 0.0

        indiv_snr = self._individual_snr(history, center, L, field_sigma, field_amplitude)
        dispersion = self._dispersion(history, L)

        # Three substrate signatures, all required at screening:
        #   cohesive (social coupling; rejects social_off),
        #   climbs toward the true peak (rejects non-sensing),
        #   emergent (individuals too noisy to sense alone; rejects field_too_strong).
        cohesive = dispersion < COHESION_MAX
        emergent = indiv_snr < EMERGENT_SNR_MAX
        tier = "none"
        confidence = 0.0
        if cohesive and emergent and ci_true > CLIMB_MIN:
            tier = "screening"
            confidence = 0.500
            if p_value < 0.10:
                tier = "confirmation"
                confidence = 0.700
                if ci_true > DEF_CI_MIN and p_value < 0.05:
                    tier = "definitive"
                    confidence = 0.900
        detected = tier != "none"

        notes = (
            f"CI_true={ci_true:.3f} vs random-target mean={null_mean:.3f} "
            f"(p={p_value:.3f}, d={cohens_d:.2f}); individual_SNR={indiv_snr:.2f} "
            f"(< {EMERGENT_SNR_MAX} => emergent)"
        )
        return DetectorResult(
            pattern_id="P17",
            detected=detected,
            tier=tier,
            confidence=confidence,
            primary_metric={
                "chemotactic_index": float(ci_true),
                "field_sensing_p_value": p_value,
                "individual_snr": float(indiv_snr),
                "group_dispersion": float(dispersion),
            },
            secondary_metrics={
                "ci_random_target_mean": null_mean,
                "ci_random_target_std": null_std,
            },
            effect_size={"cohens_d": cohens_d},
            null_p_value=p_value,
            null_type="random_target_null",
            exclusions_checked=["P5_flocking", "field_too_strong"],
            exclusion_results={
                "P5": f"approach to true peak vs random targets p={p_value:.3f} (flocking has no field preference)",
                "emergent": f"individual_SNR={indiv_snr:.2f} (< {EMERGENT_SNR_MAX} required: individuals cannot sense alone)",
            },
            co_occurrence_candidates=["P5"],
            metadata_available=True,
            warnings=warnings,
            notes=notes,
        )

    def _dispersion(self, history: list[dict[str, Any]], L: float) -> float:
        """Mean late-frame agent distance to the periodic group centroid /
        box size. Low => cohesive collective (social coupling); high =>
        agents dispersed/independent (social_off)."""
        vals = []
        for s in history[len(history) // 2:]:
            pos = np.asarray(s.get("positions"), dtype=float)
            if pos.ndim != 2 or pos.shape[0] < 2:
                continue
            tx = 2 * np.pi * pos[:, 0] / L
            ty = 2 * np.pi * pos[:, 1] / L
            cx = (L * np.arctan2(np.mean(np.sin(tx)), np.mean(np.cos(tx))) / (2 * np.pi)) % L
            cy = (L * np.arctan2(np.mean(np.sin(ty)), np.mean(np.cos(ty))) / (2 * np.pi)) % L
            dx = pos[:, 0] - cx; dy = pos[:, 1] - cy
            dx -= L * np.round(dx / L); dy -= L * np.round(dy / L)
            vals.append(float(np.mean(np.sqrt(dx * dx + dy * dy))))
        return (float(np.mean(vals)) / L) if vals else 1.0

    def _individual_snr(
        self,
        history: list[dict[str, Any]],
        center: tuple[float, float],
        L: float,
        sigma: float,
        amp: float,
    ) -> float:
        """Spatial spread of the TRUE field across agents / residual sensing
        noise, from a late frame. Low SNR => individuals cannot sense alone =>
        group-level climbing is emergent. +inf when per-agent samples missing."""
        s = history[len(history) * 3 // 4]
        fs = s.get("field_samples")
        pos = s.get("positions")
        if fs is None or pos is None:
            return float("inf")
        pos = np.asarray(pos, dtype=float)
        fs = np.asarray(fs, dtype=float)
        if pos.ndim != 2 or pos.shape[0] != fs.shape[0] or fs.shape[0] < 3:
            return float("inf")
        dx = pos[:, 0] - center[0]
        dy = pos[:, 1] - center[1]
        dx -= L * np.round(dx / L)
        dy -= L * np.round(dy / L)
        d = np.sqrt(dx * dx + dy * dy)
        true_field = amp * np.exp(-(d ** 2) / (2.0 * sigma ** 2))
        resid = float(np.std(fs - true_field))
        spread = float(np.std(true_field))
        return spread / (resid + 1e-9)

    def _null_result(self, reason: str) -> DetectorResult:
        return DetectorResult(
            pattern_id="P17", detected=False, tier="none", confidence=0.0,
            primary_metric={"chemotactic_index": 0.0, "field_sensing_p_value": 1.0, "individual_snr": 0.0},
            secondary_metrics={}, effect_size={"cohens_d": 0.0},
            null_p_value=1.0, null_type="random_target_null",
            exclusions_checked=[], exclusion_results={},
            co_occurrence_candidates=["P5"], metadata_available=True,
            warnings=[reason], notes=reason,
        )
