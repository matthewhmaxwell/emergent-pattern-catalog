"""Battery — load all per-pattern calibrators and produce a calibrated
emergent-pattern profile (Milestone C, T2a deliverable; feeds T2d).

The instrument view: given a system's observation, each detector yields its
canonical continuous scalar, the matching Calibrator maps it to a CALIBRATED,
cross-detector-comparable confidence, and the battery returns the ranked
profile. Detectors whose substrate the observation does not match reject at
prerequisite -> calibrate as 'absent', which is the correct battery behavior:
only the patterns actually present light up.
"""
from __future__ import annotations

import glob
import json
import os
from typing import Any, Dict, List, Optional

from epc.phase2a.calibration import Calibrator, build_from_panel


class Battery:
    def __init__(self, calibrators: Dict[str, Calibrator]) -> None:
        self.calibrators = calibrators

    @classmethod
    def load(cls, outputs_dir: str = "analysis/outputs") -> "Battery":
        """Load one Calibrator per pattern from the Phase-2a panel JSONs."""
        cals: Dict[str, Calibrator] = {}
        missing: List[str] = []
        for f in sorted(glob.glob(os.path.join(outputs_dir, "*_phase2a_panel.json"))):
            try:
                summ = json.load(open(f))
            except Exception:
                continue
            cal = build_from_panel(summ)
            if cal is not None:
                cals[cal.pattern_id] = cal
            else:
                missing.append(os.path.basename(f))
        if missing:
            # reference arrays not present (panel predates T2a) — surfaced, not hidden
            cals.setdefault("_missing", None)  # type: ignore
            cals.pop("_missing", None)
            cls._missing = missing  # type: ignore
        return cls(cals)

    def profile_from_values(self, values: Dict[str, Optional[float]]) -> List[Dict[str, Any]]:
        """Calibrated profile from a dict of {pattern_id: raw canonical metric}.

        Patterns absent from ``values`` (or None) calibrate as uncomputable and
        sort last. Result is ranked by calibrated_confidence descending.
        """
        out: List[Dict[str, Any]] = []
        for pid, cal in self.calibrators.items():
            out.append(cal.calibrate(values.get(pid)))
        out.sort(key=lambda r: (r.get("calibrated_confidence") is not None,
                                r.get("calibrated_confidence") or -1.0),
                 reverse=True)
        return out

    def top(self, values: Dict[str, Optional[float]], k: int = 5) -> List[Dict[str, Any]]:
        return self.profile_from_values(values)[:k]

    def __len__(self) -> int:
        return len(self.calibrators)
