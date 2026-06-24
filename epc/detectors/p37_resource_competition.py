"""P37 — Field-mediated resource competition: sustained competitive oscillation (detector).

Discriminating metric: COEXISTENCE OSCILLATION
  D* = osc_cv (median per-species late-window temporal CV)  IF n_persist >= 3
     = 0                                                      otherwise
high only when several species COEXIST while their abundances keep oscillating
(never reach a fixed point) — the Huisman-Weissing signature of competition
mediated through shared depletable resources. Rejects competitive exclusion (few
survivors), STABLE multi-species coexistence (same coexistence, but settles to a
fixed point -> osc_cv ~ 0, the key control), and single-resource exclusion.
Requires a 'resources' observable (the field), so direct-interaction oscillators
(P11 Lotka-Volterra) and phase oscillators (P9/P10) are substrate mismatches.

Validated (analysis/discovery/p37_resource_competition.py): TNR 1.0, d large.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

from epc.metrics.competition_dynamics import coexistence_oscillation

OSC_SCREEN = 0.05     # sustained-oscillation screening
OSC_CONF = 0.08      # confirmation
OSC_DEF = 0.10       # definitive
MIN_PERSIST = 3      # multi-species coexistence (the exclusion-violation essence)


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
    return DetectorResult("P37", False, "none", 0.0, {"coexistence_oscillation": 0.0}, {}, notes=note)


class P37ResourceCompetitionDetector:
    """Detector for P37 — sustained oscillatory coexistence under resource competition."""

    def __init__(self, seed: int = 42) -> None:
        self.seed = seed

    def detect(self, history: List[Dict[str, Any]],
               metadata: Optional[Dict[str, Any]] = None) -> DetectorResult:
        m = coexistence_oscillation(history)
        if m is None:
            return _none("substrate mismatch — P37 needs 'abundances'+'resources' time series.")
        osc = float(m["osc_cv"]); npr = int(m["n_persist"])
        nres = int(m["n_resources"]); nsp = int(m["n_species"])
        gates = npr >= MIN_PERSIST
        Dstar = osc if gates else 0.0

        tier, detected, conf = "none", False, 0.0
        if gates and osc > OSC_SCREEN:
            tier, detected, conf = "screening", True, 0.4
        if detected and osc > OSC_CONF and npr >= 4:
            tier, conf = "confirmation", 0.6
        if tier == "confirmation" and osc > OSC_DEF:
            tier, conf = "definitive", 0.9

        return DetectorResult(
            pattern_id="P37", detected=detected, tier=tier, confidence=round(conf, 3),
            primary_metric={"coexistence_oscillation": round(Dstar, 4)},
            secondary_metrics={"osc_cv": round(osc, 4), "n_persist": npr,
                               "n_resources": nres, "n_species": nsp,
                               "multi_species": npr >= MIN_PERSIST},
            exclusions_checked=["P11", "P9"],
            exclusion_results={"multi_species_coexistence": npr >= MIN_PERSIST,
                               "sustained_oscillation": osc > OSC_SCREEN},
            notes=f"osc_cv={osc:.3f} persist={npr}/{nsp} resources={nres}.")
