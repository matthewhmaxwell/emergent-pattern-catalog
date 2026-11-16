"""Recognition re-run after the P33-P37 merge: self-recognition for all 36 panel positives (P1-P37 except P31),
five seeds each, through the 36 single-run detectors in the battery (P31 has none; P37 at its post-panel-fix
detector), plus the same cross-model runs as recognition_rerun. Writes analysis/outputs/recognition_rerun_v2.json and
leaves recognition_rerun.json untouched.

  python -m analysis.recognition_rerun_v2

The first version of this script (ac4d6ca) redirected json.dump, but recognition_rerun.main() had already opened
recognition_rerun.json for writing, so that run emptied the v1 record. The v2 results are unaffected; the v1 record was
restored from the copy kept with the paper. This version redirects the open call instead.
"""
import builtins
import analysis.recognition_rerun as RR

RR.PATTERNS = [f"p{i}" for i in range(1, 38) if i != 31]
_SRC = "analysis/outputs/recognition_rerun.json"
_DST = "analysis/outputs/recognition_rerun_v2.json"


def _open(path, *args, **kwargs):
    return builtins.open(_DST if path == _SRC else path, *args, **kwargs)


if __name__ == "__main__":
    RR.open = _open
    RR.main()
