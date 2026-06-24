"""P36 — Heterogeneous-agent kinetic exchange: Pareto power-law wealth tail (detector).

Discriminating metric: PARETO-OVER-EXPONENTIAL TAIL ADVANTAGE
  D* = (ks_exp - ks_pareto)   IF [ alpha_hill in Pareto range AND gini high AND
        not single-agent condensed AND Pareto fit good ]
     = 0                       otherwise
high only when the steady-state wealth tail is a clean power law — the
Chatterjee-Chakrabarti-Manna (2004) signature of agent heterogeneity. Rejects the
homogeneous-lambda Gamma, the lambda=0 Boltzmann exponential, Yard-Sale
condensation (P28), equal wealth, and a static lognormal: all are
exponential-class or condensed, not Pareto. Gini alone is NOT specific (P28 and
the Gamma share high/moderate Gini), so the discriminator is the tail SHAPE.

Validated (analysis/discovery/p36_pareto_exchange.py): TNR 1.0, continuous d(D*) large.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Any, Dict, List, Optional

import numpy as np

from epc.metrics.wealth_concentration import (
    gini as _gini, hill_tail_alpha, pareto_ks_distance, exp_ks_distance,
)

A_LO, A_HI = 0.6, 1.8     # Pareto tail index (CCM 2004 ~ 1); look-alikes alpha >= 2.1
G_MIN = 0.45              # real inequality (rejects Gamma gini~0.27 and equal)
MS_MAX = 0.20             # not single-agent condensed (rejects full Yard-Sale)
KS_GATE = 0.08            # Pareto fit must be good (rejects condensation ks~0.11)
ADV_SCREEN = 0.12        # power-law-over-exponential advantage, screening
ADV_CONF = 0.18         # confirmation
ADV_DEF = 0.22          # definitive


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
    return DetectorResult("P36", False, "none", 0.0, {"pareto_advantage": 0.0}, {}, notes=note)


class P36ParetoExchangeDetector:
    """Detector for P36 — Pareto power-law wealth tail from heterogeneous exchange."""

    def __init__(self, seed: int = 42) -> None:
        self.seed = seed

    def detect(self, history: List[Dict[str, Any]],
               metadata: Optional[Dict[str, Any]] = None) -> DetectorResult:
        ws = [np.asarray(f["wealth"], dtype=np.float64) for f in history
              if isinstance(f, dict) and "wealth" in f]
        ws = [w for w in ws if w.size >= 100 and w.sum() > 0]
        if len(ws) < 2:
            return _none("substrate mismatch — P36 needs a 'wealth' (conserved scalar) time series.")
        T = len(ws)
        late = ws[(3 * T) // 4:]
        alphas, ksp, kse, gn, ms = [], [], [], [], []
        for w in late:
            a, wmin, _ = hill_tail_alpha(w)
            gn.append(_gini(w)); ms.append(float(w.max() / w.sum()))
            if np.isfinite(a) and np.isfinite(wmin) and wmin > 0:
                alphas.append(a)
                ksp.append(pareto_ks_distance(w, a, wmin))
                kse.append(exp_ks_distance(w, wmin))
        if not alphas:
            return _none("degenerate tail — no finite Pareto fit.")
        alpha = float(np.median(alphas))
        ks_par = float(np.nanmedian(ksp))
        ks_exp = float(np.nanmedian(kse))
        gini = float(np.mean(gn))
        max_share = float(np.mean(ms))
        adv = ks_exp - ks_par

        gates = ((A_LO <= alpha <= A_HI) and (gini > G_MIN)
                 and (max_share < MS_MAX) and (ks_par < KS_GATE))
        Dstar = adv if gates else 0.0

        tier, detected, conf = "none", False, 0.0
        if gates and adv > ADV_SCREEN:
            tier, detected, conf = "screening", True, 0.4
        if detected and adv > ADV_CONF and ks_par < 0.06:
            tier, conf = "confirmation", 0.6
        if tier == "confirmation" and adv > ADV_DEF and 0.8 <= alpha <= 1.5:
            tier, conf = "definitive", 0.9

        return DetectorResult(
            pattern_id="P36", detected=detected, tier=tier, confidence=round(conf, 3),
            primary_metric={"pareto_advantage": round(Dstar, 4)},
            secondary_metrics={"alpha_hill": round(alpha, 3), "ks_pareto": round(ks_par, 4),
                               "ks_exp": round(ks_exp, 4), "gini": round(gini, 3),
                               "max_share": round(max_share, 4), "gates_pass": gates},
            exclusions_checked=["P28", "P14"],
            exclusion_results={"not_condensed": max_share < MS_MAX, "power_law_over_exp": adv > 0},
            notes=f"alpha={alpha:.2f} ks_par={ks_par:.3f} ks_exp={ks_exp:.3f} "
                  f"adv={adv:+.3f} gini={gini:.2f} max_share={max_share:.3f}.")
