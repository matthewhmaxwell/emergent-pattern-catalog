"""T2a — Calibration layer for the EPC instrument.

Turns a detector's raw canonical metric value into a CALIBRATED, cross-detector
score using the faithful reference distributions from the validation rebuild
(positive seeds + the negative panel; the per-pattern scalar named in
`epc/phase2a/continuous_metrics.py`). This is what lets the battery be pointed at
an arbitrary system and report a comparable per-pattern profile, instead of an
uncalibrated raw tier.

Every pattern's scalar is first ORIENTED (multiplied by its direction) so that
"higher = more pattern-like" holds uniformly. For an observation's oriented
value ``ox`` against reference positives P and negatives N:

  pos_percentile        fraction of P with value <= ox  -> where it sits among
                        genuine positives (0.5 = typical positive, 1.0 = beyond
                        all positives, 0.0 = below all positives).
  neg_tail_p            (1 + #{n in N : n >= ox}) / (1 + |N|)  -> conformal tail
                        probability: how likely a NEGATIVE is at least this
                        pattern-like. Small = unlikely under the negatives.
  calibrated_confidence P(ox more pattern-like than a random negative), the
                        AUC-style separation score in [0,1] (1 = clearly
                        positive-side, ~0.5 = ambiguous, 0 = negative-side).
                        This is the single comparable instrument score.

`neg_tail_p` is the primitive T2b ("none-of-the-above") thresholds on; the
calibrated_confidence is the headline per-pattern score for T2d profiles.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Any, Dict, List, Optional

import numpy as np


@dataclass
class Calibrator:
    pattern_id: str
    metric: str
    direction: int          # +1 higher-is-positive, -1 lower-is-positive
    pos: np.ndarray         # ORIENTED positive reference values
    neg: np.ndarray         # ORIENTED negative reference values

    @property
    def saturated(self) -> bool:
        return float(np.ptp(self.pos)) < 1e-12 and float(np.ptp(self.neg)) < 1e-12

    def calibrate(self, value: Optional[float]) -> Dict[str, Any]:
        if value is None or not np.isfinite(value):
            return {"pattern_id": self.pattern_id, "metric": self.metric,
                    "value": None, "pos_percentile": None, "neg_tail_p": None,
                    "calibrated_confidence": None, "verdict": "uncomputable"}
        ox = self.direction * float(value)
        pos_pct = float(np.mean(self.pos <= ox)) if self.pos.size else float("nan")
        n = self.neg.size
        neg_tail_p = (1.0 + float(np.sum(self.neg >= ox))) / (1.0 + n) if n else float("nan")
        # AUC-style separation: fraction of negatives strictly below, ties half.
        if n:
            conf = float(np.mean(self.neg < ox) + 0.5 * np.mean(self.neg == ox))
        else:
            conf = float("nan")
        # PROVISIONAL verdict from the robust AUC-style separation (the
        # conformal neg_tail_p floors at ~1/(|N|+1), too coarse for hard p<=0.05
        # gates at the current reference sizes). T2c sets final operating points
        # from held-out precision/recall; these tiers are placeholders.
        if not np.isfinite(conf):
            verdict = "uncomputable"
        elif conf >= 0.98:
            verdict = "strong"
        elif conf >= 0.90:
            verdict = "present"
        elif conf >= 0.70:
            verdict = "weak"
        else:
            verdict = "absent"
        return {"pattern_id": self.pattern_id, "metric": self.metric,
                "value": float(value), "oriented_value": ox,
                "pos_percentile": pos_pct, "neg_tail_p": neg_tail_p,
                "calibrated_confidence": conf, "saturated": self.saturated,
                "verdict": verdict}


def build_from_panel(panel_summary: Dict[str, Any]) -> Optional[Calibrator]:
    """Build a Calibrator from a panel JSON's stored reference value arrays."""
    s = panel_summary.get("summary", panel_summary)
    pos = s.get("continuous_pos_values")
    neg = s.get("continuous_neg_values")
    if not pos or not neg:
        return None
    direction = -1 if s.get("continuous_metric_direction") == "lower" else 1
    return Calibrator(
        pattern_id=panel_summary.get("pattern_id", s.get("pattern_id", "?")),
        metric=s.get("continuous_metric", "?"),
        direction=direction,
        pos=direction * np.asarray(pos, dtype=float),
        neg=direction * np.asarray(neg, dtype=float),
    )
