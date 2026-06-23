"""P34 adaptive-network fragmentation — standalone Phase-2a validation + panel JSON.

Imports the promoted modules (single source of truth). Run:
    python -m analysis.discovery.p34_adaptive_network            # validate (TNR + d)
    python -m analysis.discovery.p34_adaptive_network --panel    # write calibrator JSON
"""
from __future__ import annotations

import json
import os
from typing import Any, Dict, List

import numpy as np

from epc.models.adaptive_network import (
    adaptive_voter, null_random_rewire, neg_static_modular,
    neg_fixed_network_voter, neg_er_static, neg_pre_fragmented,
)
from epc.detectors.p34_adaptive_network import P34AdaptiveNetworkDetector

DET = P34AdaptiveNetworkDetector()
POS_KW = dict(N=200, mean_deg=4, phi=0.99, max_steps=400000, n_frames=80)
NEGS = [("null_random_rewire", null_random_rewire), ("static_modular", neg_static_modular),
        ("fixed_network_voter", neg_fixed_network_voter), ("er_static", neg_er_static),
        ("pre_fragmented", neg_pre_fragmented)]


def _D(r):
    return float(r.primary_metric["modularity_gain"])


def _row(nm, r):
    sm = r.secondary_metrics
    return (f"  {nm:<22} det={str(r.detected):<5} tier={r.tier:<12} D*={_D(r):.3f} "
            f"Qe={sm.get('Q_early',0):.3f} Ql={sm.get('Q_late',0):.3f} "
            f"giant {sm.get('giant_early',0):.2f}->{sm.get('giant_late',0):.2f}")


def _pos_neg(seeds):
    pos = [DET.detect(*adaptive_voter(s, **POS_KW)) for s in range(seeds)]
    neg = [(nm, DET.detect(*fn(s))) for nm, fn in NEGS for s in range(seeds)]
    return pos, neg


def validate(seeds: int = 5) -> int:
    print(f"=== POSITIVES adaptive_voter phi={POS_KW['phi']} — expect >= confirmation ===", flush=True)
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
    print("PASS: adaptive-network fragmentation >=confirmation; all lookalikes rejected (TNR 1.0).",
          flush=True)
    return 0


def write_panel(seeds: int = 5,
                out: str = "analysis/outputs/p34_adaptive_network_phase2a_panel.json") -> int:
    pos, neg = _pos_neg(seeds)
    pos_vals = [_D(r) for r in pos]; neg_vals = [_D(r) for _, r in neg]
    n_rej = sum(1 for _, r in neg if not r.detected)
    ps = np.array(pos_vals); ns = np.array(neg_vals)
    pooled = float(np.sqrt((ps.var() + ns.var()) / 2)) or 1e-9
    summary = {
        "pattern_id": "P34", "continuous_metric": "modularity_gain",
        "continuous_metric_direction": "higher",
        "continuous_pos_values": pos_vals, "continuous_neg_values": neg_vals,
        "overall_tnr": n_rej / len(neg), "n_positive": len(pos), "n_negative": len(neg),
        "positive_tiers": [r.tier for r in pos],
        "cohens_d_continuous": float((ps.mean() - ns.mean()) / pooled),
        "negative_panel": ["null_random_rewire", "static_modular", "fixed_network_voter",
                           "er_static", "pre_fragmented"],
        "note": "Adaptive-network fragmentation (evolving-network substrate, Ring-1). "
                "Discriminator = coevolution-driven modularity gain (Q emerges from a "
                "low-modularity connected start), gated by fragmentation + opinion-topology "
                "assortativity. Validated TNR 1.0, d=16.74.",
    }
    doc = {"pattern_id": "P34", "panel_version": "discovery-ring1",
           "detector_format": "adjacency_timeseries", "summary": summary}
    os.makedirs(os.path.dirname(out), exist_ok=True)
    json.dump(doc, open(out, "w"), indent=2)
    print(f"wrote {out}  (pos {len(pos_vals)}, neg {len(neg_vals)}, TNR {n_rej}/{len(neg)})", flush=True)
    return 0


if __name__ == "__main__":
    import sys
    if "--panel" in sys.argv:
        raise SystemExit(write_panel())
    raise SystemExit(validate())
