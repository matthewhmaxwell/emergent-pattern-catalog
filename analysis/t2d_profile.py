"""T2d profiler — point the catalog battery at LLM-swarm trace(s), blind.

Loads one or more trace JSONs written by analysis/t2d_llm_swarm.py and runs
analysis.battery_profile.profile_observation on each, in both the native
(screening) mode and the recommended OOD mode (match_min_tier='confirmation').
Reports the verdict, the ranked top-5 profile, and the generic-emergence score —
without telling the battery which condition produced the trace.
"""
from __future__ import annotations

import json
import sys

import numpy as np

from analysis.battery_profile import build_detector_fns, profile_observation
from epc.phase2a.battery import Battery


def load(path):
    d = json.load(open(path))
    hist = []
    for f in d["history"]:
        nf = {"positions": np.asarray(f["positions"], dtype=float)}
        if "velocities" in f:
            nf["velocities"] = np.asarray(f["velocities"], dtype=float)
        if "types" in f:
            nf["types"] = np.asarray(f["types"], dtype=int)
        hist.append(nf)
    return d, hist


def main():
    bat = Battery.load()
    fns = build_detector_fns()
    for path in sys.argv[1:]:
        d, hist = load(path)
        cond = d.get("condition", "?")
        print(f"\n=== {cond}  ({len(hist)} frames, n={d.get('params',{}).get('n')}) ===", flush=True)
        for mt in ("screening", "confirmation"):
            out = profile_observation(hist, None, battery=bat, detector_fns=fns,
                                      match_min_tier=mt)
            v = out["verdict"]
            em = out["emergence"]
            line = (f"  [{mt:<12}] {v.get('verdict'):<22} "
                    f"pat={v.get('pattern_id')} demoted={v.get('demoted_match')} "
                    f"em={em['score']} kind={em['kind']}")
            print(line, flush=True)
            if mt == "confirmation":
                top = [(r["pattern_id"], r.get("verdict"),
                        round(r.get("calibrated_confidence") or -1, 2),
                        r.get("detected"), r.get("detector_tier"))
                       for r in out["profile"][:5]]
                print("     top5:", top, flush=True)


if __name__ == "__main__":
    main()
