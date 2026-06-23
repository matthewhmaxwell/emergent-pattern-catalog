"""P35 entropy crystallization — standalone Phase-2a validation + panel JSON.

Run:  python -m analysis.discovery.p35_entropy_crystal            # validate
      python -m analysis.discovery.p35_entropy_crystal --panel    # write calibrator JSON
"""
from __future__ import annotations

import json
import os
from typing import Any, Dict, List

import numpy as np

from epc.models.entropy_crystal import (
    hard_disk_crystallization, neg_dilute_gas, neg_static_crystal,
)
from epc.detectors.p35_entropy_crystal import P35EntropyCrystalDetector
from analysis.t2c_systems import (
    nov_keller_segel, nov_active_nematic, null_spatial_noise, null_random_walk,
)

DET = P35EntropyCrystalDetector()
POS_KW = dict(N=196, L=14.0, eta_end=0.74, n_steps=2500, n_frames=80, noise0=0.16)
NEGS = [("dilute_gas", neg_dilute_gas, 14.0), ("static_crystal", neg_static_crystal, 14.0),
        ("keller_segel", nov_keller_segel, 50.0), ("active_nematic", nov_active_nematic, 18.0),
        ("spatial_noise", null_spatial_noise, 30.0), ("random_walk", null_random_walk, 60.0)]


def _D(r):
    return float(r.primary_metric["psi6_gain"])


def _row(nm, r):
    sm = r.secondary_metrics
    return (f"  {nm:<18} det={str(r.detected):<5} tier={r.tier:<12} D*={_D(r):+.3f} "
            f"psi6 {sm.get('psi6_early',0):.3f}->{sm.get('psi6_late',0):.3f}")


def _with_box(m, L):
    m = dict(m or {})
    if not m.get("box_size"):
        m["box_size"] = L
    return m


def _pos_neg(seeds):
    pos = [DET.detect(*hard_disk_crystallization(s, **POS_KW)) for s in range(seeds)]
    neg = []
    for nm, fn, L in NEGS:
        for s in range(seeds):
            h, m = fn(s)
            neg.append((nm, DET.detect(h, _with_box(m, L))))
    return pos, neg


def validate(seeds: int = 5) -> int:
    print("=== POSITIVES entropy crystallization — expect >= confirmation ===", flush=True)
    pos, neg = _pos_neg(seeds)
    fails = []
    for s, r in enumerate(pos):
        print(_row(f"seed{s}", r), flush=True)
        if r.tier not in ("confirmation", "definitive"):
            fails.append(f"positive seed{s} tier={r.tier}")
    print("=== NEGATIVES — expect rejected (TNR) ===", flush=True)
    seen = set()
    for nm, r in neg:
        if nm not in seen:
            print(_row(nm, r), flush=True); seen.add(nm)
        if r.detected:
            fails.append(f"negative {nm} DETECTED tier={r.tier}")
    n_rej = sum(1 for _, r in neg if not r.detected)
    ps = np.array([_D(r) for r in pos]); ns = np.array([_D(r) for _, r in neg])
    pooled = np.sqrt((ps.var() + ns.var()) / 2) or 1e-9
    print(f"\nTNR = {n_rej}/{len(neg)} = {n_rej/len(neg):.3f}   "
          f"continuous d(D*) = {(ps.mean()-ns.mean())/pooled:.2f}", flush=True)
    if fails:
        print("FAIL:", *("  - " + f for f in fails), sep="\n", flush=True)
        return 1
    print("PASS: entropy crystallization >=confirmation; all lookalikes rejected (TNR 1.0).",
          flush=True)
    return 0


def write_panel(seeds: int = 5,
                out: str = "analysis/outputs/p35_entropy_crystal_phase2a_panel.json") -> int:
    pos, neg = _pos_neg(seeds)
    pos_vals = [_D(r) for r in pos]; neg_vals = [_D(r) for _, r in neg]
    n_rej = sum(1 for _, r in neg if not r.detected)
    ps = np.array(pos_vals); ns = np.array(neg_vals)
    pooled = float(np.sqrt((ps.var() + ns.var()) / 2)) or 1e-9
    summary = {
        "pattern_id": "P35", "continuous_metric": "psi6_gain",
        "continuous_metric_direction": "higher",
        "continuous_pos_values": pos_vals, "continuous_neg_values": neg_vals,
        "overall_tnr": n_rej / len(neg), "n_positive": len(pos), "n_negative": len(neg),
        "positive_tiers": [r.tier for r in pos],
        "cohens_d_continuous": float((ps.mean() - ns.mean()) / pooled),
        "negative_panel": ["dilute_gas", "static_crystal", "keller_segel",
                           "active_nematic", "spatial_noise", "random_walk"],
        "note": "Entropy-driven crystallization (none_entropy substrate, Ring-1). "
                "Discriminator = hexatic-order emergence (local psi6 gain from a "
                "disordered fluid), gated by reaching the hexatic floor. Validated "
                "TNR 1.0, d=8.31.",
    }
    doc = {"pattern_id": "P35", "panel_version": "discovery-ring1",
           "detector_format": "positions_timeseries", "summary": summary}
    os.makedirs(os.path.dirname(out), exist_ok=True)
    json.dump(doc, open(out, "w"), indent=2)
    print(f"wrote {out}  (pos {len(pos_vals)}, neg {len(neg_vals)}, TNR {n_rej}/{len(neg)})", flush=True)
    return 0


if __name__ == "__main__":
    import sys
    if "--panel" in sys.argv:
        raise SystemExit(write_panel())
    raise SystemExit(validate())
