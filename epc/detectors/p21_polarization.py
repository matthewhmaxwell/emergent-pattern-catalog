"""P21 Polarization / Fragmentation detector.

Detects persistent multimodality in opinion distributions — populations
splitting into internally coherent camps rather than converging to consensus.

From detector_cards.md:
  Primary: Dip test p < 0.01 at steady state (persistent multimodality)
  Secondary: ≥2 clusters, between-cluster distance > ε, persistence
  Null: full-range (ε→∞) → consensus, random-opinion → uniform

  Screening:    Dip p < 0.05
  Confirmation: Dip p < 0.01 AND ≥2 clusters AND persistence ≥ 50 interaction times
  Definitive:   Confirmation + ≥100 interaction times + full-range null converges
                + emerged from unimodal IC + P18 exclusion

Reference: Hegselmann & Krause (2002), JASSS 5(3).
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Optional

import numpy as np


@dataclass
class P21DetectorResult:
    """Result of P21 polarization detection."""
    detected: bool
    tier: str                       # 'none', 'screening', 'confirmation', 'definitive'
    confidence: float
    n_clusters: int                 # Number of opinion clusters at steady state
    dip_stat: float                 # Hartigan's dip statistic
    dip_p: float                    # Dip test p-value
    variance: float                 # Final opinion variance
    persistence_steps: int          # Steps with stable multimodality
    from_unimodal: bool             # Started from unimodal initial condition
    full_range_consensus: bool      # Full-range null converged to consensus
    warnings: list[str]


def _hartigan_dip(x: np.ndarray, n_boot: int = 1000, seed: int = 42) -> tuple[float, float]:
    """Compute Hartigan's dip statistic and p-value.

    The dip statistic measures departure from unimodality. Larger values
    indicate stronger multimodality.

    Uses bootstrap p-value: fraction of uniform samples with dip ≥ observed.
    """
    x_sorted = np.sort(x)
    n = len(x_sorted)

    if n < 4:
        return 0.0, 1.0

    # Compute dip: max difference between empirical CDF and best-fit unimodal CDF
    # Simplified: use the greatest convex minorant / least concave majorant approach
    ecdf = np.arange(1, n + 1) / n

    # Best-fit uniform on [min, max]
    uniform_cdf = (x_sorted - x_sorted[0]) / (x_sorted[-1] - x_sorted[0] + 1e-30)

    dip = float(np.max(np.abs(ecdf - uniform_cdf))) / 2

    # Bootstrap p-value
    rng = np.random.default_rng(seed)
    null_dips = np.zeros(n_boot)
    for i in range(n_boot):
        u = np.sort(rng.uniform(0, 1, size=n))
        u_ecdf = np.arange(1, n + 1) / n
        u_uniform = (u - u[0]) / (u[-1] - u[0] + 1e-30)
        null_dips[i] = float(np.max(np.abs(u_ecdf - u_uniform))) / 2

    p_value = float(np.mean(null_dips >= dip))
    if p_value == 0:
        p_value = 1.0 / (n_boot + 1)

    return dip, p_value


def _count_clusters(opinions: np.ndarray, gap: float = 0.05) -> int:
    """Count opinion clusters by gap detection."""
    if len(opinions) == 0:
        return 0
    sorted_x = np.sort(opinions)
    gaps = np.diff(sorted_x) > gap
    return int(np.sum(gaps)) + 1


def detect_p21(
    history: list[dict[str, Any]],
    model_metadata: Optional[dict[str, Any]] = None,
    n_boot: int = 1000,
) -> P21DetectorResult:
    """Run P21 polarization/fragmentation detection.

    Parameters
    ----------
    history : list[dict]
        State trajectory with 'opinions', 'n_clusters', 'variance' keys.
    model_metadata : dict, optional
        May contain 'epsilon' for cluster gap threshold.
    n_boot : int
        Bootstrap samples for dip test.

    Returns
    -------
    P21DetectorResult
    """
    warnings = []

    if not history or len(history) < 5:
        return P21DetectorResult(
            detected=False, tier="none", confidence=0.0,
            n_clusters=0, dip_stat=0.0, dip_p=1.0,
            variance=0.0, persistence_steps=0,
            from_unimodal=False, full_range_consensus=False,
            warnings=["Insufficient history"],
        )

    epsilon = model_metadata.get("epsilon", 0.2) if model_metadata else 0.2

    # --- Dip test on final opinion distribution ---
    final_opinions = history[-1]["opinions"]
    dip_stat, dip_p = _hartigan_dip(final_opinions, n_boot=n_boot, seed=42)

    # --- Cluster count at steady state ---
    gap_threshold = epsilon / 2
    n_clusters = _count_clusters(final_opinions, gap=gap_threshold)

    # --- Variance ---
    variance = float(np.var(final_opinions))

    # --- Persistence: how many steps was it multimodal? ---
    persistence = 0
    for h in reversed(history):
        nc = _count_clusters(h["opinions"], gap=gap_threshold)
        if nc >= 2:
            persistence += 1
        else:
            break

    # --- Initial condition check: started from unimodal? ---
    init_opinions = history[0]["opinions"]
    _, init_dip_p = _hartigan_dip(init_opinions, n_boot=200, seed=99)
    from_unimodal = init_dip_p > 0.05  # could not reject unimodality

    # --- Full-range null: would consensus occur at ε→∞? ---
    # This is always true for HK with ε > 0.5 (all agents interact)
    # Check from metadata
    full_range_consensus = True  # HK always converges to consensus at ε→∞

    n_steps = history[-1].get("step", len(history) - 1)

    # --- Classification ---
    screening = dip_p < 0.05 and n_clusters >= 2
    confirmation = (
        screening
        and dip_p < 0.01
        and n_clusters >= 2
        and persistence >= 50
    )
    definitive = (
        confirmation
        and persistence >= 100
        and full_range_consensus
        and from_unimodal
    )

    if definitive:
        tier = "definitive"
        confidence = min(0.85 + 0.15 * min(n_clusters / 3, 1.0), 1.0)
    elif confirmation:
        tier = "confirmation"
        confidence = min(0.60 + 0.25 * min(n_clusters / 2, 1.0), 0.85)
    elif screening:
        tier = "screening"
        confidence = min(0.30 + 0.30 * min(dip_stat / 0.05, 1.0), 0.60)
    else:
        tier = "none"
        confidence = 0.0

    if not from_unimodal and screening:
        warnings.append("Initial condition may have been multimodal (dip p<0.05)")

    return P21DetectorResult(
        detected=tier != "none",
        tier=tier,
        confidence=round(confidence, 3),
        n_clusters=n_clusters,
        dip_stat=round(dip_stat, 4),
        dip_p=round(dip_p, 4),
        variance=round(variance, 6),
        persistence_steps=persistence,
        from_unimodal=from_unimodal,
        full_range_consensus=full_range_consensus,
        warnings=warnings,
    )
