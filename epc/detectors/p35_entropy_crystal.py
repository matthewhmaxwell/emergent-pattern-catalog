"""P35 — Entropy-driven crystallization (detector).

Discriminating metric: HEXATIC-ORDER EMERGENCE (psi6 gain)
  D* = psi6_gain (= psi6_late - psi6_early)   IF psi6_late >= PSI6_MIN
     = 0                                       otherwise
positive ONLY for crystallization. Rejects aggregation/clumping (RCP, no gain), a
fluid/random pattern (psi6 below the hexatic floor), and a static pre-formed crystal
(high psi6 but ZERO gain). The GAIN — hexatic order EMERGING from a disordered start
— is the clean discriminator (mirrors P34 modularity-gain / P3 wavelength-growth).

Validated (analysis/discovery/p35_entropy_crystal.py): positives 5/5 DEFINITIVE,
TNR 30/30 = 1.000, continuous d = 8.31.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

import numpy as np

from epc.metrics.positional_order import local_psi6

PSI6_MIN = 0.42       # must reach hexatic order (fluid/clump background ~0.30-0.40)
GAIN_SCREEN = 0.04    # hexatic order emerged from a disordered start
GAIN_CONF = 0.07      # confirmation (positives gain ~0.08-0.13)


@dataclass
class DetectorResult:
    pattern_id: str
    detected: bool
    tier: str
    confidence: float
    primary_metric: dict
    secondary_metrics: dict
    null_p_value: float = 1.0
    null_type: str = "panel-discrimination"
    exclusions_checked: list = field(default_factory=list)
    exclusion_results: dict = field(default_factory=dict)
    notes: str = ""


def _none(note: str) -> DetectorResult:
    return DetectorResult("P35", False, "none", 0.0, {"psi6_gain": 0.0}, {}, notes=note)


class P35EntropyCrystalDetector:
    """Detector for P35 — emergent hexatic (bond-orientational) order from packing."""

    def __init__(self, seed: int = 42) -> None:
        self.seed = seed

    def detect(self, history: List[Dict[str, Any]],
               metadata: Optional[Dict[str, Any]] = None) -> DetectorResult:
        ps = [np.asarray(f["positions"]) for f in history if "positions" in f]
        if len(ps) < 6:
            return _none("substrate mismatch — P35 needs a 'positions' time series.")
        L = (metadata or {}).get("box_size")
        if not L:
            allp = np.vstack(ps)
            L = float(allp.max() - allp.min()) or None
        early = ps[1:4]
        late = ps[-5:]                  # equilibrated final state (annealing ends late)
        psi6_early = float(np.mean([local_psi6(p, L) for p in early]))
        psi6_late = float(np.mean([local_psi6(p, L) for p in late]))
        gain = psi6_late - psi6_early
        sustained = min(local_psi6(p, L) for p in late) > PSI6_MIN - 0.06

        gates = psi6_late >= PSI6_MIN
        Dstar = gain if gates else 0.0
        tier = "none"; detected = False; conf = 0.0
        if Dstar > GAIN_SCREEN:
            tier, detected, conf = "screening", True, 0.4
        if detected and Dstar > GAIN_CONF and sustained:
            tier, conf = "confirmation", 0.6
        if tier == "confirmation" and psi6_late >= 0.45:
            tier, conf = "definitive", 0.9

        return DetectorResult(
            pattern_id="P35", detected=detected, tier=tier, confidence=round(conf, 3),
            primary_metric={"psi6_gain": Dstar},
            secondary_metrics={"psi6_early": psi6_early, "psi6_late": psi6_late,
                               "psi6_gain": gain, "gates_pass": gates,
                               "sustained": sustained},
            exclusions_checked=["P1", "P2"],
            exclusion_results={"hexatic_not_clump": gates, "emergent_not_static": Dstar > 0},
            notes=f"psi6 {psi6_early:.2f}->{psi6_late:.2f} (gain {gain:+.2f}).")
