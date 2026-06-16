"""Battery profiler — run ALL detectors on one observation -> calibrated profile.

The instrument keystone (Milestone C). For a single observation it runs every
detector, extracts each pattern's canonical scalar, calibrates it against that
pattern's reference distribution, and returns the ranked per-pattern profile plus
the generic-emergence score and the three-way verdict (MATCH /
EMERGENT-UNCLASSIFIED / NO-EMERGENCE). Detectors whose substrate the observation
doesn't match reject at prerequisite (or error) -> metric absent -> calibrated
'absent', which is correct: only the patterns actually present light up.

This is what T2c (held-out systems) and T2d (external demo) call.
"""
from __future__ import annotations

from typing import Any, Dict, List, Optional

import analysis.run_phase2a_panel as R
from epc.phase2a.battery import Battery
from epc.phase2a.continuous_metrics import CONTINUOUS_METRIC, extract_continuous
from epc.phase2a.emergence import generic_emergence, three_way_verdict
from epc.phase2a.panel import _detected, _verdict

_TIER_RANK = {'definitive': 3, 'confirmation': 2, 'screening': 1, 'none': 0}

PATTERNS = [f"p{i}" for i in range(1, 31)] + ["p32"]


def build_detector_fns() -> Dict[str, Any]:
    fns: Dict[str, Any] = {}
    for pid in PATTERNS:
        mk = getattr(R, f"make_{pid}_detector_fn", None)
        if mk is not None:
            fns[pid.upper()] = mk()
    return fns


def profile_observation(history: List[Dict[str, Any]],
                        metadata: Optional[Dict[str, Any]] = None,
                        battery: Optional[Battery] = None,
                        detector_fns: Optional[Dict[str, Any]] = None,
                        seed: int = 0) -> Dict[str, Any]:
    battery = battery if battery is not None else Battery.load()
    detector_fns = detector_fns if detector_fns is not None else build_detector_fns()
    rows: List[Dict[str, Any]] = []
    for pid_up, fn in detector_fns.items():
        entry = CONTINUOUS_METRIC.get(pid_up)
        cal = battery.calibrators.get(pid_up)
        raw = None
        detected = False
        tier = "none"
        try:
            res = fn(history, metadata)
            raw = extract_continuous(res, entry)
            detected = _detected(res)
            tier = _verdict(res)
        except Exception:
            raw, detected, tier = None, False, "none"
        row = cal.calibrate(raw) if cal is not None else {
            "pattern_id": pid_up, "calibrated_confidence": None, "verdict": "no-calibrator"}
        row["detected"] = bool(detected)
        row["detector_tier"] = tier if detected else "none"
        row["tier_rank"] = _TIER_RANK.get(tier, 0) if detected else 0
        row["fired_detected"] = bool(detected)
        rows.append(row)
    # Rank by the detector's actual tier (its multi-criterion gate + exclusions
    # disambiguate within-family overlap that a single calibrated scalar cannot),
    # then by calibrated confidence as the comparable magnitude.
    rows.sort(key=lambda r: (1 if r.get("fired_detected") else 0,
                             r.get("tier_rank", 0),
                             r.get("calibrated_confidence") or -1.0), reverse=True)
    em = generic_emergence(history, seed=seed)
    top = rows[0] if rows else None
    # MATCH iff the top detector ACTUALLY FIRED (det=True) via its full gate; the
    # tier is the strength. If nothing fired, the verdict is emergence-driven
    # (EMERGENT-UNCLASSIFIED vs NO-EMERGENCE).
    if top is not None and top.get("fired_detected"):
        verdict = {"verdict": "MATCH", "pattern_id": top["pattern_id"],
                   "detector_tier": top["detector_tier"],
                   "calibrated_confidence": top.get("calibrated_confidence"),
                   "emergence_score": em["score"]}
    else:
        verdict = three_way_verdict(None, em["score"])
    return {"profile": rows, "emergence": em, "verdict": verdict,
            "top": {"pattern_id": top["pattern_id"],
                    "detector_tier": top.get("detector_tier"),
                    "calibrated_confidence": top.get("calibrated_confidence")} if top else None}
