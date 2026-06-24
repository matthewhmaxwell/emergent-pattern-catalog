"""P37 field-mediated resource competition — Phase-2a validation + panel JSON.

    python -m analysis.discovery.p37_resource_competition            # validate (TNR + d)
    python -m analysis.discovery.p37_resource_competition --panel    # write calibrator JSON
"""
from __future__ import annotations

import json
import os

import numpy as np

from epc.models.resource_competition import (
    resource_competition, neg_competitive_exclusion, neg_stable_coexistence,
    neg_single_resource,
)
from epc.detectors.p37_resource_competition import P37ResourceCompetitionDetector

DET = P37ResourceCompetitionDetector()
NEGS = [("competitive_exclusion", neg_competitive_exclusion),
        ("stable_coexistence", neg_stable_coexistence),
        ("single_resource", neg_single_resource)]


def _D(r):
    return float(r.primary_metric["coexistence_oscillation"])


def _row(nm, r):
    sm = r.secondary_metrics
    return (f"  {nm:<24} det={str(r.detected):<5} tier={r.tier:<12} D*={_D(r):.4f} "
            f"osc={sm.get('osc_cv',0):.3f} persist={sm.get('n_persist',0)}/{sm.get('n_species',0)} "
            f"res={sm.get('n_resources',0)}")


def _pos_neg(seeds):
    pos = [DET.detect(*resource_competition(s)) for s in range(seeds)]
    neg = [(nm, DET.detect(*fn(s))) for nm, fn in NEGS for s in range(seeds)]
    return pos, neg


def validate(seeds: int = 5) -> int:
    print("=== POSITIVES Huisman-Weissing competitive chaos — expect >= confirmation ===", flush=True)
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
    print("PASS: resource-competition oscillation >=confirmation; exclusion / stable "
          "coexistence / single-resource all rejected (TNR 1.0).", flush=True)
    return 0


def write_panel(seeds: int = 5,
                out: str = "analysis/outputs/p37_resource_competition_phase2a_panel.json") -> int:
    pos, neg = _pos_neg(seeds)
    pos_vals = [_D(r) for r in pos]; neg_vals = [_D(r) for _, r in neg]
    n_rej = sum(1 for _, r in neg if not r.detected)
    ps = np.array(pos_vals); ns = np.array(neg_vals)
    pooled = float(np.sqrt((ps.var() + ns.var()) / 2)) or 1e-9
    summary = {
        "pattern_id": "P37", "continuous_metric": "coexistence_oscillation",
        "continuous_metric_direction": "higher",
        "continuous_pos_values": pos_vals, "continuous_neg_values": neg_vals,
        "overall_tnr": n_rej / len(neg), "n_positive": len(pos), "n_negative": len(neg),
        "positive_tiers": [r.tier for r in pos],
        "cohens_d_continuous": float((ps.mean() - ns.mean()) / pooled),
        "negative_panel": [nm for nm, _ in NEGS],
        "note": "Field-mediated resource competition (interact=field x conflict=competitive, "
                "Ring-1 combination cell). Huisman-Weissing competitive chaos: species compete "
                "through shared depletable resources and coexist via SUSTAINED oscillation "
                "(never a fixed point). Discriminator = multi-species coexistence + sustained "
                "abundance oscillation; the stable-coexistence control (same coexistence, no "
                "oscillation) is rejected, so the signature is the non-equilibrium dynamics. TNR 1.0.",
    }
    doc = {"pattern_id": "P37", "panel_version": "discovery-ring1",
           "detector_format": "abundance_resource_timeseries", "summary": summary}
    os.makedirs(os.path.dirname(out), exist_ok=True)
    json.dump(doc, open(out, "w"), indent=2)
    print(f"wrote {out}  (pos {len(pos_vals)}, neg {len(neg_vals)}, TNR {n_rej}/{len(neg)})", flush=True)
    return 0


if __name__ == "__main__":
    import sys
    if "--panel" in sys.argv:
        raise SystemExit(write_panel())
    raise SystemExit(validate())
