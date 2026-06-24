"""P36 heterogeneous kinetic exchange (Pareto tail) — Phase-2a validation + panel JSON.

    python -m analysis.discovery.p36_pareto_exchange            # validate (TNR + d)
    python -m analysis.discovery.p36_pareto_exchange --panel    # write calibrator JSON
"""
from __future__ import annotations

import json
import os

import numpy as np

from epc.models.kinetic_exchange import (
    kinetic_exchange, neg_uniform_saving, neg_boltzmann, neg_condensation,
    neg_equal, neg_lognormal,
)
from epc.detectors.p36_pareto_exchange import P36ParetoExchangeDetector

DET = P36ParetoExchangeDetector()
NEGS = [("uniform_saving_gamma", neg_uniform_saving), ("boltzmann_exp", neg_boltzmann),
        ("yard_sale_condensation", neg_condensation), ("equal_wealth", neg_equal),
        ("lognormal_static", neg_lognormal)]


def _D(r):
    return float(r.primary_metric["pareto_advantage"])


def _row(nm, r):
    sm = r.secondary_metrics
    return (f"  {nm:<24} det={str(r.detected):<5} tier={r.tier:<12} D*={_D(r):+.3f} "
            f"alpha={sm.get('alpha_hill',0):.2f} ks_par={sm.get('ks_pareto',0):.3f} "
            f"gini={sm.get('gini',0):.2f} ms={sm.get('max_share',0):.3f}")


def _pos_neg(seeds):
    pos = [DET.detect(*kinetic_exchange(s)) for s in range(seeds)]
    neg = [(nm, DET.detect(*fn(s))) for nm, fn in NEGS for s in range(seeds)]
    return pos, neg


def validate(seeds: int = 5) -> int:
    print("=== POSITIVES kinetic_exchange (heterogeneous lambda) — expect >= confirmation ===", flush=True)
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
    print("PASS: Pareto-tail exchange >=confirmation; all wealth look-alikes rejected (TNR 1.0).",
          flush=True)
    return 0


def write_panel(seeds: int = 5,
                out: str = "analysis/outputs/p36_pareto_exchange_phase2a_panel.json") -> int:
    pos, neg = _pos_neg(seeds)
    pos_vals = [_D(r) for r in pos]; neg_vals = [_D(r) for _, r in neg]
    n_rej = sum(1 for _, r in neg if not r.detected)
    ps = np.array(pos_vals); ns = np.array(neg_vals)
    pooled = float(np.sqrt((ps.var() + ns.var()) / 2)) or 1e-9
    summary = {
        "pattern_id": "P36", "continuous_metric": "pareto_advantage",
        "continuous_metric_direction": "higher",
        "continuous_pos_values": pos_vals, "continuous_neg_values": neg_vals,
        "overall_tnr": n_rej / len(neg), "n_positive": len(pos), "n_negative": len(neg),
        "positive_tiers": [r.tier for r in pos],
        "cohens_d_continuous": float((ps.mean() - ns.mean()) / pooled),
        "negative_panel": [nm for nm, _ in NEGS],
        "note": "Heterogeneous-agent kinetic exchange (resource_exchange x heterogeneous, "
                "Ring-1 combination cell). Discriminator = Pareto-over-exponential tail "
                "advantage (ks_exp - ks_pareto), gated by Pareto tail index ~1, high Gini, "
                "not condensed. The power-law tail is the CCM-2004 signature of agent "
                "heterogeneity; Gini alone is shared with P28 condensation. TNR 1.0.",
    }
    doc = {"pattern_id": "P36", "panel_version": "discovery-ring1",
           "detector_format": "wealth_timeseries", "summary": summary}
    os.makedirs(os.path.dirname(out), exist_ok=True)
    json.dump(doc, open(out, "w"), indent=2)
    print(f"wrote {out}  (pos {len(pos_vals)}, neg {len(neg_vals)}, TNR {n_rej}/{len(neg)})", flush=True)
    return 0


if __name__ == "__main__":
    import sys
    if "--panel" in sys.argv:
        raise SystemExit(write_panel())
    raise SystemExit(validate())
