"""Ring-2 Stage 1 — particle-life sourcing pilot (a LENS TEST, not a hunt).

Samples N_SAMPLES random parameter sets (the K x K force matrix), runs each particle-
life system, pushes it through the instrument (named lenses via generic_emergence + the
model-free novelty tripwire), and records one row per system. Then answers the three
questions that tell us whether the lens stack is trustworthy before we invest further:
  Q1 separation  — do different parameter sets light different lenses / spread the descriptor?
  Q2 tripwire    — does anything come up COMPLEX but UNCLASSIFIED (a novelty lead)?
  Q3 null floor  — do the dead/trivial sets stay quiet (no false trips)?

    python -m analysis.ring2.stage1_particle_life_pilot
"""
from __future__ import annotations

import json
import signal
from collections import Counter

import numpy as np

from epc.models.particle_life import particle_life
from epc.phase2a.novelty_tripwire import novelty_tripwire

N_SAMPLES = 300
OUT = "analysis/outputs/ring2_stage1_particle_life.jsonl"
PER_SAMPLE_TIMEOUT = 30   # robustness: no single pathological parameter set can stall the run


class _Timeout(Exception):
    pass


def _alarm(s, f):
    raise _Timeout()


signal.signal(signal.SIGALRM, _alarm)


def _raw_descriptor(history):
    """Cheap interpretable structural axes from the late window."""
    late = history[len(history) // 2:]
    L = 10.0
    pol, spd, clus = [], [], []
    for f in late:
        v = np.asarray(f["velocities"], float); p = np.asarray(f["positions"], float)
        sp = np.linalg.norm(v, axis=1)
        spd.append(float(sp.mean()))
        if sp.mean() > 1e-6:
            u = v / (np.linalg.norm(v, axis=1, keepdims=True) + 1e-12)
            pol.append(float(np.linalg.norm(u.mean(0))))     # polar order (flocking)
        h = np.histogram2d(p[:, 0], p[:, 1], bins=10, range=[[0, L], [0, L]])[0]
        clus.append(float(h.std() / (h.mean() + 1e-9)))       # density CV (clumping)
    return {"polar": round(float(np.mean(pol)) if pol else 0.0, 3),
            "speed": round(float(np.mean(spd)), 3),
            "clustering": round(float(np.mean(clus)), 3)}


def run():
    rows = []
    for i in range(N_SAMPLES):
        rng = np.random.default_rng(i)
        K = int(rng.integers(3, 7))
        F = rng.uniform(-1.0, 1.0, size=(K, K))
        try:
            signal.alarm(PER_SAMPLE_TIMEOUT)
            hist, meta = particle_life(i, K=K, F=F, N=200, steps=250)
            tw = novelty_tripwire(hist)
            desc = _raw_descriptor(hist)
            signal.alarm(0)
            rows.append({"i": i, "K": K, "F": F.round(3).tolist(),
                         "tripped": tw["tripped"], "is_complex": tw["is_complex"],
                         "classified": tw["classified"], "C": tw["C"], "psi": tw["psi"],
                         "em": tw["em_score"], "em_kind": tw["em_kind"], **desc})
        except _Timeout:
            signal.alarm(0)
            rows.append({"i": i, "K": K, "error": f"timeout>{PER_SAMPLE_TIMEOUT}s"})
        except Exception as e:
            signal.alarm(0)
            rows.append({"i": i, "K": K, "error": repr(e)[:80]})
        if (i + 1) % 50 == 0:
            print(f"  ...{i + 1}/{N_SAMPLES}", flush=True)

    import os
    os.makedirs("analysis/outputs", exist_ok=True)
    with open(OUT, "w") as fh:
        for r in rows:
            fh.write(json.dumps(r) + "\n")

    ok = [r for r in rows if "error" not in r]
    errs = [r for r in rows if "error" in r]
    classified = [r for r in ok if r["classified"]]
    tripped = [r for r in ok if r["tripped"]]
    dead = [r for r in ok if not r["classified"] and not r["is_complex"]]
    kinds = Counter(r["em_kind"] for r in classified)

    print(f"\n=== Stage-1 pilot: {len(ok)} systems ({len(errs)} errors) -> {OUT} ===")
    print(f"\nQ1 SEPARATION — lenses lit across the sample (em>=0.50):")
    for k, c in kinds.most_common():
        print(f"    {str(k):<34} {c}")
    print(f"    classified total: {len(classified)}/{len(ok)}")
    em = np.array([r["em"] for r in ok]); C = np.array([r["C"] or 0 for r in ok])
    print(f"    em_score range {em.min():.2f}-{em.max():.2f} (mean {em.mean():.2f});  "
          f"C range {C.min():.3f}-{C.max():.3f} (mean {C.mean():.3f})")
    pol = np.array([r["polar"] for r in ok]); cl = np.array([r["clustering"] for r in ok])
    print(f"    polar range {pol.min():.2f}-{pol.max():.2f};  clustering range {cl.min():.2f}-{cl.max():.2f}")

    print(f"\nQ2 TRIPWIRE — complex-but-unclassified (novelty leads): {len(tripped)}")
    for r in tripped[:12]:
        print(f"    i={r['i']} K={r['K']} C={r['C']} psi={r['psi']} em={r['em']} "
              f"polar={r['polar']} clustering={r['clustering']}")

    print(f"\nQ3 NULL FLOOR — dead/trivial (unclassified + not complex): {len(dead)}/{len(ok)} "
          f"({100*len(dead)//max(len(ok),1)}%)")
    null_complex = [r for r in dead if r["is_complex"]]
    print(f"    of the dead, falsely flagged complex: {len(null_complex)} (want 0)")
    print("\nNOTE: Stage 1 is a lens test. Tripwire hits are LEADS to inspect, not claims.")
    return 0


if __name__ == "__main__":
    raise SystemExit(run())
