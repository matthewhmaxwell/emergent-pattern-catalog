"""P33 — Active nematic order with ±1/2 topological defects (detector).

Discriminating metric: COHERENT half-integer defect density
  D* = half_def_density   IF (S_loc ≥ S_MIN AND φ ≤ PHI_MAX AND |L| ≤ L_MAX)
     = 0                   otherwise
high only for an active nematic, ~0 for every lookalike (polar flock fails φ;
milling fails |L|; isotropic fails S_loc; ordered nematic fails defect-presence).
The ±1/2 winding is the nematic fingerprint no polar/integer-vortex system has.

Validated to the Phase-2a bar (analysis/discovery/p33_active_nematic.py): positives
5/5 definitive, TNR 20/20 = 1.000, continuous d(D*) = 4.81; real P5 Vicsek rejected;
finite-size definitive at G = 64/96/128.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

import numpy as np

from epc.metrics.nematic_order import (
    angular_momentum, director_field, half_integer_defect_density,
    local_nematic_order, polar_order,
)

# --- Detector-card thresholds (calibrated against the real positive/negative
#     distributions: pos S_loc 0.43-0.67, disordered-Vicsek 0.35, isotropic 0.15) ---
S_MIN = 0.40        # local nematic order gate
PHI_MAX = 0.35      # polar-order ceiling (apolar; flock φ≈0.99)
L_MAX = 0.25        # angular-momentum ceiling (not milling; milling |L|≈0.69)
D_SCREEN = 0.004    # coherent half-defect density, screening
D_CONF = 0.010      # confirmation (positives D*≈0.05-0.14)


@dataclass
class DetectorResult:
    pattern_id: str
    detected: bool
    tier: str                      # "none"|"screening"|"confirmation"|"definitive"
    confidence: float
    primary_metric: dict
    secondary_metrics: dict
    null_p_value: float
    null_type: str = "shuffle(D*)-diagnostic"
    exclusions_checked: list = field(default_factory=list)
    exclusion_results: dict = field(default_factory=dict)
    notes: str = ""


def _none(pid: str, note: str) -> DetectorResult:
    return DetectorResult(pid, False, "none", 0.0,
                          {"coherent_half_defect_density": 0.0}, {}, 1.0, notes=note)


class P33ActiveNematicDetector:
    """Detector for P33 — active nematic order with ±1/2 topological defects.

    Tiers are deterministic + faithful (the 20-case Phase-2a negative panel is the
    discrimination rigor; the shuffle null on D* is reported as a diagnostic — it is
    fragile on the most-ordered seeds whose local order survives a permutation)."""

    def __init__(self, n_permutations: int = 199, seed: int = 42) -> None:
        self.n_permutations = n_permutations
        self.seed = seed

    def detect(self, history: List[Dict[str, Any]],
               metadata: Optional[Dict[str, Any]] = None) -> DetectorResult:
        rng = np.random.default_rng(self.seed)
        if not history or "velocities" not in history[0]:
            return _none("P33", "substrate mismatch — P33 needs an orientation field "
                                "(velocities/headings or theta_field).")
        T = len(history)
        meas = history[T // 2:]

        def fstats(f):
            theta = director_field(f)
            if theta is None:
                return None
            sloc = local_nematic_order(theta)
            hdef, idef = half_integer_defect_density(theta)
            ang = f["headings"] if "headings" in f else np.arctan2(
                f["velocities"][:, 1], f["velocities"][:, 0])
            phi = polar_order(ang)
            L = (angular_momentum(f["positions"], np.asarray(f.get("headings", ang)),
                                  float(f.get("box_size", theta.shape[0])))
                 if "positions" in f else 0.0)
            return sloc, phi, hdef, idef, L, theta

        stats = [s for s in (fstats(f) for f in meas) if s is not None]
        if not stats:
            return _none("P33", "no director field could be formed.")
        S_loc = float(np.mean([s[0] for s in stats]))
        phi = float(np.mean([s[1] for s in stats]))
        hdef = float(np.mean([s[2] for s in stats]))
        idef = float(np.mean([s[3] for s in stats]))
        L = float(np.mean([s[4] for s in stats]))
        min_hdef = float(min(s[2] for s in stats))
        init = fstats(history[0])
        S_loc_init = init[0] if init is not None else 0.0

        gates = (S_loc >= S_MIN) and (phi <= PHI_MAX) and (L <= L_MAX)
        Dstar = hdef if gates else 0.0

        # diagnostic null on D*: shuffle a representative (median-defect) frame
        hdefs = [s[2] for s in stats]
        theta_obs = stats[int(np.argsort(hdefs)[len(hdefs) // 2])][5]
        flat = theta_obs.ravel()
        null_D = np.empty(self.n_permutations)
        for k in range(self.n_permutations):
            th = rng.permutation(flat).reshape(theta_obs.shape)
            sl = local_nematic_order(th)
            hd, _ = half_integer_defect_density(th)
            null_D[k] = hd if (sl >= S_MIN and phi <= PHI_MAX and L <= L_MAX) else 0.0
        p_val = float((np.sum(null_D >= Dstar) + 1) / (self.n_permutations + 1))

        emerged = S_loc > S_loc_init + 0.15
        sustained = min_hdef > D_SCREEN

        tier = "none"; detected = False; conf = 0.0
        if Dstar > D_SCREEN:
            tier, detected, conf = "screening", True, 0.4
        if detected and Dstar > D_CONF and emerged and sustained:
            tier, conf = "confirmation", 0.6
        if tier == "confirmation" and phi <= PHI_MAX and L <= L_MAX and idef < hdef:
            tier = "definitive"; conf = 0.85 + (0.1 if p_val < 0.01 else 0.0)

        return DetectorResult(
            pattern_id="P33", detected=detected, tier=tier, confidence=round(conf, 3),
            primary_metric={"coherent_half_defect_density": Dstar},
            secondary_metrics={"S_loc": S_loc, "polar_phi": phi,
                               "half_def_density": hdef, "integer_def_density": idef,
                               "angular_momentum": L, "S_loc_init": S_loc_init,
                               "min_hdef": min_hdef, "emerged": emerged,
                               "sustained": sustained, "gates_pass": gates},
            null_p_value=p_val,
            exclusions_checked=["P5", "P6"],
            exclusion_results={"P5_polar": phi <= PHI_MAX, "P6_milling": L <= L_MAX},
            notes=f"D*={Dstar:.4f} via coherent ±1/2 defect gas; gates "
                  f"S_loc={S_loc:.2f}/φ={phi:.2f}/|L|={L:.2f}.")
