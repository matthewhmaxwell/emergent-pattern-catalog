"""Exit 0 only when every task behind a group B gating rule has a result file; print what is still missing.

Gating: P5 N=400 old_velocity (90 runs); P6 T6.1/T6.2 (40 ring runs at the dt chosen by the dt check);
P10 Stage A N=256 (10), Stage B for any paper-IC chimera, and the pre-registered fallback (attempts and their
Stage B ramps) until it has nothing left to run. Reported-only tasks may still be queued.

Usage (worktree root): sudo -u matthewhmaxwell env PYTHONPATH=. /home/matthewhmaxwell/epc-venv/bin/python \
    analysis/reproductions/repro_v2/gating_done.py
"""
from __future__ import annotations

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))

import p5_vicsek1995 as p5  # noqa: E402
import p6_dorsogna2006 as p6  # noqa: E402
import p10_abrams_strogatz2004 as p10  # noqa: E402


def missing(rawdir: Path, tasks: list[dict]) -> list[str]:
    return [t["id"] for t in tasks if not (rawdir / f"{t['id']}.json").exists()]


def main() -> int:
    out = {}
    out["P5 N=400 old_velocity"] = missing(p5.RAWDIR, [t for t in p5.all_tasks(False)
                                                        if t["N"] == 400 and t["update"] == "old_velocity"])
    dtc = p6.dt_check_outcome()
    out["P6 T6.1/T6.2"] = (missing(p6.RAWDIR, p6.ring_tasks_main(dtc["dt_for_T61"])) if dtc
                           else ["dt check not finished"])
    a = missing(p10.RAWDIR, p10.tasks_stage_a(256, 0.18))
    if a:
        out["P10 T10.1/T10.2"] = a
    else:
        tb, _ = p10.tasks_stage_b(256)
        tf, _ = p10.tasks_fallback()
        out["P10 T10.1/T10.2"] = missing(p10.RAWDIR, tb) + [t["id"] for t in tf]
        need = p10.fallback_need()
        att = [k for k in range(p10.FALLBACK_MAX_ATTEMPTS) if (p10.RAWDIR / f"P10_F_N256_k{k}.json").exists()]
        print(f"P10 fallback: needed {need}, attempts finished {len(att)} (seeds {[100 + k for k in att]}), "
              f"chimera realizations in use {len(p10.fallback_used_ks())}"
              + (" -- ATTEMPT CAP REACHED BEFORE 10 CHIMERAS" if need and len(att) >= p10.FALLBACK_MAX_ATTEMPTS
                 and len(p10.fallback_used_ks()) < need else ""))
    errs = sorted(str(f.name) for f in p5.RAWDIR.parent.glob("*/*.error.json"))
    done = all(not v for v in out.values())
    for k, v in out.items():
        print(f"{k}: {'done' if not v else f'{len(v)} missing, e.g. ' + ', '.join(v[:3])}")
    if errs:
        print("ERROR files:", errs)
    print("GATING_DONE" if done else "GATING_NOT_DONE")
    return 0 if done else 1


if __name__ == "__main__":
    sys.exit(main())
