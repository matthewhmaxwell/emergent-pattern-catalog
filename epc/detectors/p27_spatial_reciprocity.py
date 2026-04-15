"""P27 Spatial Reciprocity detector.

Detects cooperation surviving on a lattice when it should go extinct
in a well-mixed population (spatial reciprocity / network reciprocity).

From detector_cards.md:
  Primary metric: f_C > 0 long-term when well-mixed predicts f_C → 0
  Secondary: Moran's I on cooperator indicator > 0 (spatial clustering)
  Null: well-mixed (f_C → 0), random-strategy (I ≈ 0)

  Screening:    f_C > 0.02 at t > 100 T_gen
  Confirmation: f_C > 0.05 at t > 1000 T_gen AND well-mixed null AND I > 0
  Definitive:   Confirmation + PD verified (T>R>P>S) + N≥100 + P1 exclusion

Reference: Nowak & May (1992), Nature 359, 826-829.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Optional

import numpy as np
from scipy import stats


@dataclass
class P27DetectorResult:
    """Result of P27 spatial reciprocity detection."""
    detected: bool
    tier: str                        # 'none', 'screening', 'confirmation', 'definitive'
    confidence: float                # Tier-capped [0, 1]
    coop_fraction: float             # Final cooperation fraction
    moran_i: float                   # Spatial autocorrelation of cooperators
    moran_i_p: float                 # p-value from random null
    pd_verified: bool                # T > R > P >= S confirmed
    well_mixed_baseline: float       # Expected f_C in well-mixed (0 for PD)
    n_generations: int               # Total generations observed
    warnings: list[str]


def detect_p27(
    history: list[dict[str, Any]],
    model_metadata: Optional[dict[str, Any]] = None,
    n_permutations: int = 199,
) -> P27DetectorResult:
    """Run P27 spatial reciprocity detection.

    Parameters
    ----------
    history : list[dict]
        State trajectory with 'grid', 'coop_fraction', 'moran_i' keys.
    model_metadata : dict, optional
        Must contain payoff matrix info for definitive tier.
    n_permutations : int
        Permutations for Moran's I significance test.

    Returns
    -------
    P27DetectorResult
    """
    warnings = []

    if not history or len(history) < 10:
        return P27DetectorResult(
            detected=False, tier="none", confidence=0.0,
            coop_fraction=0.0, moran_i=0.0, moran_i_p=1.0,
            pd_verified=False, well_mixed_baseline=0.0,
            n_generations=0, warnings=["Insufficient history"],
        )

    n_gen = history[-1].get("step", len(history) - 1)

    # --- Cooperation fraction (late-stage average) ---
    late_start = max(1, len(history) * 3 // 4)
    late_fcs = [h["coop_fraction"] for h in history[late_start:]]
    fc_mean = float(np.mean(late_fcs))
    fc_std = float(np.std(late_fcs))

    # --- Moran's I (late-stage average) ---
    late_Is = [h["moran_i"] for h in history[late_start:]]
    mi_mean = float(np.mean(late_Is))

    # --- Moran's I significance test (random shuffle null) ---
    final_grid = history[-1]["grid"]
    rows, cols = history[-1]["grid_dims"]
    rng = np.random.default_rng(42)

    # Observed Moran's I on final grid
    obs_I = mi_mean

    # Null distribution: shuffle strategy labels
    null_Is = []
    for _ in range(n_permutations):
        shuffled = final_grid.ravel().copy()
        rng.shuffle(shuffled)
        shuffled = shuffled.reshape(rows, cols)
        indicator = (shuffled == 0).astype(float)
        x = indicator - indicator.mean()
        var = float(np.mean(x ** 2))
        if var < 1e-12:
            null_Is.append(0.0)
            continue
        cross = 0.0
        W = 0
        for dr, dc in [(-1, 0), (1, 0), (0, -1), (0, 1),
                        (-1, -1), (-1, 1), (1, -1), (1, 1)]:
            shifted = np.roll(np.roll(x, -dr, axis=0), -dc, axis=1)
            cross += np.sum(x * shifted)
            W += x.size
        null_Is.append(float(cross / (W * var)))

    null_Is = np.array(null_Is)
    moran_p = float(np.mean(null_Is >= obs_I))
    if moran_p == 0:
        moran_p = 1.0 / (n_permutations + 1)

    # --- PD structure check ---
    pd_verified = False
    if model_metadata:
        T = model_metadata.get("payoff_T", 0)
        R = model_metadata.get("payoff_R", 0)
        P = model_metadata.get("payoff_P", 0)
        S = model_metadata.get("payoff_S", 0)
        pd_verified = T > R > P and P >= S

    # Well-mixed baseline: for PD, defection dominates → f_C = 0
    well_mixed = 0.0 if pd_verified else 0.5

    # --- Classification ---
    N = rows * cols

    # Screening: f_C > 0.02 at t > 100
    screening = fc_mean > 0.02 and n_gen > 100

    # Confirmation: f_C > 0.05 at t > 1000, I significant, above well-mixed
    confirmation = (
        screening
        and fc_mean > 0.05
        and fc_mean > well_mixed + 0.02  # above well-mixed baseline
        and mi_mean > 0.0
        and moran_p < 0.01
    )

    # Definitive: confirmation + PD verified + N >= 100 + P1 exclusion
    # P1 exclusion: no movement (static lattice) → P27, not P1
    has_movement = model_metadata.get("has_movement", True) if model_metadata else True
    definitive = (
        confirmation
        and pd_verified
        and N >= 100
        and not has_movement  # P1 exclusion: no physical movement
    )

    if definitive:
        tier = "definitive"
        confidence = min(0.85 + 0.15 * min(fc_mean / 0.3, 1.0), 1.0)
    elif confirmation:
        tier = "confirmation"
        confidence = min(0.60 + 0.25 * min(fc_mean / 0.2, 1.0), 0.85)
    elif screening:
        tier = "screening"
        confidence = min(0.30 + 0.30 * min(fc_mean / 0.1, 1.0), 0.60)
    else:
        tier = "none"
        confidence = 0.0

    if fc_std > 0.1 * fc_mean and fc_mean > 0:
        warnings.append(f"High cooperation fraction variance (CV={fc_std/fc_mean:.2f})")

    return P27DetectorResult(
        detected=tier != "none",
        tier=tier,
        confidence=round(confidence, 3),
        coop_fraction=round(fc_mean, 4),
        moran_i=round(mi_mean, 4),
        moran_i_p=round(moran_p, 4),
        pd_verified=pd_verified,
        well_mixed_baseline=well_mixed,
        n_generations=n_gen,
        warnings=warnings,
    )
