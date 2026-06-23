"""P34 — Adaptive-network fragmentation (detector).

Discriminating metric: COEVOLUTION-DRIVEN MODULARITY GAIN
  D* = Q_gain (= Q_late - Q_early)   IF (Q_late >= Q_MIN AND Q_early <= Q_EARLY_MAX)
     = 0                              otherwise
high only when community structure EMERGES from a low-modularity connected start.
Rejects random rewiring (no community), a static modular graph (no gain), a
fixed-topology voter (no topology change), static ER (no structure), a
pre-fragmented graph (no emergence). Confirmed by giant-component fragmentation
and opinion-topology assortativity (state-driven, not random).

Validated (analysis/discovery/p34_adaptive_network.py): positives 5/5 DEFINITIVE,
TNR 25/25 = 1.000, continuous d(D*) = 16.74.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

import numpy as np

from epc.metrics.network_emergence import (
    binarize, giant_fraction, modularity, state_assortativity,
)

Q_MIN = 0.30          # final modularity must reach a community state
Q_EARLY_MAX = 0.22    # started from a low-modularity (connected) graph
GAIN_SCREEN = 0.08    # coevolution-driven modularity gain, screening
GAIN_CONF = 0.15      # confirmation


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
    return DetectorResult("P34", False, "none", 0.0, {"modularity_gain": 0.0}, {}, notes=note)


class P34AdaptiveNetworkDetector:
    """Detector for P34 — emergent community fragmentation in a coevolving network."""

    def __init__(self, seed: int = 42) -> None:
        self.seed = seed

    def detect(self, history: List[Dict[str, Any]],
               metadata: Optional[Dict[str, Any]] = None) -> DetectorResult:
        frames = [(f, binarize(f)) for f in history]
        frames = [(f, A) for f, A in frames if A is not None and A.sum() > 0]
        if len(frames) < 4:
            return _none("substrate mismatch — P34 needs an 'adjacency' time series.")
        T = len(frames)
        early = frames[:3]                  # genuine initial state (rewiring is fast)
        late = frames[(3 * T) // 4:]

        Q_early = float(np.mean([modularity(A) for _, A in early]))
        Q_late = float(np.mean([modularity(A) for _, A in late]))
        Q_gain = Q_late - Q_early
        giant_early = float(np.mean([giant_fraction(A) for _, A in early]))
        giant_late = float(np.mean([giant_fraction(A) for _, A in late]))
        la = [(f, A) for f, A in late if "opinions" in f]
        assort = (float(np.mean([state_assortativity(A, np.asarray(f["opinions"]))
                                 for f, A in la])) if la else None)

        gates = (Q_late >= Q_MIN) and (Q_early <= Q_EARLY_MAX)
        Dstar = Q_gain if gates else 0.0
        fragmented = giant_late < giant_early - 0.10
        assort_ok = (assort is None) or (assort >= 0.85)

        tier = "none"; detected = False; conf = 0.0
        if Dstar > GAIN_SCREEN:
            tier, detected, conf = "screening", True, 0.4
        if detected and Dstar > GAIN_CONF and fragmented and assort_ok:
            tier, conf = "confirmation", 0.6
        if tier == "confirmation" and Q_late >= 0.35 and giant_late < 0.8:
            tier, conf = "definitive", 0.9

        return DetectorResult(
            pattern_id="P34", detected=detected, tier=tier, confidence=round(conf, 3),
            primary_metric={"modularity_gain": Dstar},
            secondary_metrics={"Q_early": Q_early, "Q_late": Q_late, "Q_gain": Q_gain,
                               "giant_early": giant_early, "giant_late": giant_late,
                               "assortativity": assort, "fragmented": fragmented,
                               "gates_pass": gates},
            exclusions_checked=["P18", "P21", "P29"],
            exclusion_results={"fixed_topology": fragmented, "emergent_not_static": gates},
            notes=f"Q {Q_early:.2f}->{Q_late:.2f} (gain {Q_gain:.2f}); giant "
                  f"{giant_early:.2f}->{giant_late:.2f}; assort={assort}.")
