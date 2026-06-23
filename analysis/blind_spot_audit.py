"""Blind-spot / detection-recall audit.

Runs every probe through three detection channels:
  CURRENT  — the instrument's generic_emergence indicator (4 morphology channels)
  PSI      — Psi_CE synergy / causal emergence (best over a generic macro bank)
  MPR-C    — Martin-Plastino-Rosso statistical complexity (structure ⟂ entropy)

Calibrates the PSI and MPR-C decision thresholds against the NULL probes, then
reports per-probe hits and a per-family recall table: what the current
instrument catches vs misses, and how much PSI + MPR-C recover.
"""
from __future__ import annotations

import json
import numpy as np

from analysis.blind_spot_probes import PROBES
from analysis.emergence_channels import mpr_complexity, psi_ce_best
from epc.phase2a.emergence import generic_emergence

CUR_THR = 0.50          # generic_emergence verdict threshold (EMERGENT-UNCLASSIFIED)


def _candidates(p):
    micro = np.asarray(p["micro"], float)
    macro = np.asarray(p["macro"], float)
    c = {"provided": macro, "mean": micro.mean(1), "std": micro.std(1)}
    u = np.unique(np.round(micro, 6))
    if set(u.tolist()).issubset({0.0, 1.0}):            # binary micro → add collective features
        c["parity"] = micro.astype(int).sum(1) % 2
        c["majority"] = (micro.mean(1) > 0.5).astype(int)
    return c


def run():
    rows = []
    for fn in PROBES:
        try:
            p = fn(0)
            em = generic_emergence(p["history"])
            psi, feat = psi_ce_best(p["micro"], _candidates(p))
            mpr = mpr_complexity(p["macro"])
            rows.append({"name": fn.__name__, "family": p["family"], "truth": p["truth"],
                         "cur": round(float(em["score"]), 3),
                         "psi": (round(float(psi), 3) if psi == psi else None), "psi_feat": feat,
                         "C": round(float(mpr["C"]), 3), "H": round(float(mpr["H"]), 3)})
        except Exception as e:
            rows.append({"name": fn.__name__, "family": "?", "truth": "?", "error": repr(e)[:140]})

    nulls = [r for r in rows if r.get("truth") == "null"]
    psi_thr = max(0.05, max([(r.get("psi") or -9) for r in nulls], default=0) + 0.02)
    C_thr = max([r.get("C", 0) for r in nulls], default=0) + 0.03
    for r in rows:
        if "error" in r:
            continue
        r["cur_hit"] = r["cur"] >= CUR_THR
        r["psi_hit"] = (r.get("psi") or -9) > psi_thr
        r["C_hit"] = r["C"] > C_thr
        r["any_new"] = r["psi_hit"] or r["C_hit"]

    print(f"thresholds: PSI>{psi_thr:.3f}  MPR-C>{C_thr:.3f}\n")
    hdr = f"{'probe':<22}{'family':<22}{'truth':<9}{'CUR':>6}{'hit':>4} {'PSI':>7}{'hit':>4} {'MPR-C':>7}{'hit':>4}"
    print(hdr); print("-" * len(hdr))
    for r in rows:
        if "error" in r:
            print(f"{r['name']:<22}ERROR: {r['error']}"); continue
        print(f"{r['name']:<22}{r['family']:<22}{r['truth']:<9}"
              f"{r['cur']:>6}{'Y' if r['cur_hit'] else '.':>4} "
              f"{(r['psi'] if r['psi'] is not None else float('nan')):>7.3f}{'Y' if r['psi_hit'] else '.':>4} "
              f"{r['C']:>7.3f}{'Y' if r['C_hit'] else '.':>4}")

    em_rows = [r for r in rows if r.get("truth") == "emergent"]
    cur_recall = np.mean([r["cur_hit"] for r in em_rows]) if em_rows else 0
    new_recall = np.mean([r["cur_hit"] or r["any_new"] for r in em_rows]) if em_rows else 0
    null_fp_cur = np.mean([r["cur_hit"] for r in nulls]) if nulls else 0
    null_fp_new = np.mean([(r["cur_hit"] or r["any_new"]) for r in nulls]) if nulls else 0
    missed_cur = [r["family"] for r in em_rows if not r["cur_hit"]]
    recovered = [r["family"] for r in em_rows if not r["cur_hit"] and r["any_new"]]
    still_missed = [r["family"] for r in em_rows if not r["cur_hit"] and not r["any_new"]]

    print("\n=== RECALL on emergent probes ===")
    print(f"  CURRENT instrument:           {cur_recall:.2f}  ({sum(r['cur_hit'] for r in em_rows)}/{len(em_rows)})")
    print(f"  CURRENT + PSI + MPR-C:        {new_recall:.2f}  ({sum(r['cur_hit'] or r['any_new'] for r in em_rows)}/{len(em_rows)})")
    print(f"  null false-positives  current:{null_fp_cur:.2f}  with new channels:{null_fp_new:.2f}")
    print(f"\n  families MISSED by current : {sorted(set(missed_cur))}")
    print(f"  RECOVERED by PSI/MPR-C     : {sorted(set(recovered))}")
    print(f"  still missed (build needed): {sorted(set(still_missed))}")
    json.dump(rows, open("analysis/outputs/blind_spot_audit.json", "w"), indent=1)
    print("\nwrote analysis/outputs/blind_spot_audit.json")


if __name__ == "__main__":
    run()
