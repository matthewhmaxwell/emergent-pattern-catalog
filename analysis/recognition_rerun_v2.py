"""Recognition re-run after the P33-P37 merge: self-recognition for all 36 panel positives (P1-P37 except P31),
five seeds each, through the full set of 37 detectors (P37 at its post-panel-fix detector), plus the same cross-model
runs as recognition_rerun. Writes analysis/outputs/recognition_rerun_v2.json; recognition_rerun.json is kept.

  python -m analysis.recognition_rerun_v2
"""
import json
import analysis.recognition_rerun as RR

RR.PATTERNS = [f"p{i}" for i in range(1, 38) if i != 31]
_dump = json.dump


def _redirect(obj, fh, **kw):
    fh.close()
    with open("analysis/outputs/recognition_rerun_v2.json", "w") as f:
        _dump(obj, f, **kw)


if __name__ == "__main__":
    RR.json.dump = _redirect
    RR.main()
