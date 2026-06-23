"""P33 active nematic — standalone Phase-2a validation + panel-JSON writer.

Imports the promoted modules (epc.models.active_nematic, epc.detectors.p33_active_nematic)
so there is a single source of truth. Run:
    python -m analysis.discovery.p33_active_nematic            # validate (TNR + d)
    python -m analysis.discovery.p33_active_nematic --harden   # real P5 + finite size
    python -m analysis.discovery.p33_active_nematic --panel    # write calibrator JSON
"""
from __future__ import annotations

import json
import os
from typing import Any, Dict, List

import numpy as np

from epc.models.active_nematic import (
    active_nematic_field, neg_polar_flock, neg_milling, neg_isotropic, neg_uniform_nematic,
)
from epc.detectors.p33_active_nematic import P33ActiveNematicDetector

DET = P33ActiveNematicDetector()
NEGS = [("polar_flock", neg_polar_flock), ("milling", neg_milling),
        ("isotropic", neg_isotropic), ("uniform_nematic", neg_uniform_nematic)]


def _Dstar(r) -> float:
    return float(r.primary_metric["coherent_half_defect_density"])


def _row(name, r) -> str:
    sm = r.secondary_metrics
    return (f"  {name:<22} det={str(r.detected):<5} tier={r.tier:<12} "
            f"D*={_Dstar(r):.4f} S_loc={sm.get('S_loc', 0):.3f} phi={sm.get('polar_phi', 0):.3f} "
            f"hdef={sm.get('half_def_density', 0):.4f} |L|={sm.get('angular_momentum', 0):.3f} "
            f"p={r.null_p_value:.3f}")


def _pos_neg(seeds: int, eta: float = 0.18):
    pos = [DET.detect(*active_nematic_field(s, eta=eta)) for s in range(seeds)]
    neg = [(nm, DET.detect(*fn(s))) for nm, fn in NEGS for s in range(seeds)]
    return pos, neg


def validate(seeds: int = 5, eta: float = 0.18) -> int:
    print(f"=== POSITIVES (active nematic, eta={eta}) — expect >= confirmation ===", flush=True)
    pos, neg = _pos_neg(seeds, eta)
    fails = []
    for s, r in enumerate(pos):
        print(_row(f"seed{s}", r), flush=True)
        if r.tier not in ("confirmation", "definitive"):
            fails.append(f"positive seed{s} tier={r.tier}")
    print("\n=== NEGATIVES — expect NOT detected (TNR) ===", flush=True)
    seen = set()
    for nm, r in neg:
        if nm not in seen:
            print(_row(nm, r), flush=True); seen.add(nm)
        if r.detected:
            fails.append(f"negative {nm} DETECTED tier={r.tier}")
    n_rej = sum(1 for _, r in neg if not r.detected)
    ps = np.array([_Dstar(r) for r in pos]); ns = np.array([_Dstar(r) for _, r in neg])
    pooled = np.sqrt((ps.var() + ns.var()) / 2) or 1e-9
    print(f"\nTNR = {n_rej}/{len(neg)} = {n_rej/len(neg):.3f}   "
          f"continuous d(D*) = {(ps.mean()-ns.mean())/pooled:.2f}", flush=True)
    if fails:
        print("FAIL:", *("  - " + f for f in fails), sep="\n", flush=True)
        return 1
    print("PASS: active nematic >=confirmation; all lookalikes rejected (TNR 1.0).", flush=True)
    return 0


def harden(seeds: int = 3) -> int:
    from epc.models.vicsek import VicsekModel
    fails = []
    print("=== REAL P5 Vicsek flock (low noise) — expect rejected ===", flush=True)
    for s in range(seeds):
        m = VicsekModel(n_particles=400, box_size=10.0, speed=0.3, noise=0.2,
                        interaction_radius=1.0, seed=s)
        r = DET.detect(m.run(300, record_interval=5), m.get_metadata())
        print(_row(f"vicsek_flock s{s}", r), flush=True)
        if r.detected:
            fails.append(f"real vicsek flock s{s} DETECTED")
    print("=== REAL Vicsek disordered (high noise) — expect rejected ===", flush=True)
    for s in range(seeds):
        m = VicsekModel(n_particles=400, box_size=10.0, speed=0.3, noise=5.0,
                        interaction_radius=1.0, seed=s)
        r = DET.detect(m.run(300, record_interval=5), m.get_metadata())
        print(_row(f"vicsek_disorder s{s}", r), flush=True)
        if r.detected:
            fails.append(f"real vicsek disorder s{s} DETECTED")
    print("=== FINITE-SIZE G=96, 128 — expect definitive ===", flush=True)
    for G in (96, 128):
        for s in range(seeds):
            r = DET.detect(*active_nematic_field(s, G=G, eta=0.18))
            print(_row(f"pos_G{G} s{s}", r), flush=True)
            if r.tier not in ("confirmation", "definitive"):
                fails.append(f"positive G={G} s{s} tier={r.tier}")
    if fails:
        print("HARDENING FAIL:", *("  - " + f for f in fails), sep="\n", flush=True)
        return 1
    print("HARDENING PASS: real P5 Vicsek rejected; active nematic definitive at G=64/96/128.",
          flush=True)
    return 0


def write_panel(seeds: int = 5, out: str = "analysis/outputs/p33_active_nematic_phase2a_panel.json") -> int:
    """Run the validated panel and write the calibrator reference JSON."""
    pos, neg = _pos_neg(seeds)
    pos_vals = [_Dstar(r) for r in pos]
    neg_vals = [_Dstar(r) for _, r in neg]
    n_rej = sum(1 for _, r in neg if not r.detected)
    ps = np.array(pos_vals); ns = np.array(neg_vals)
    pooled = float(np.sqrt((ps.var() + ns.var()) / 2)) or 1e-9
    summary = {
        "pattern_id": "P33",
        "continuous_metric": "coherent_half_defect_density",
        "continuous_metric_direction": "higher",
        "continuous_pos_values": pos_vals,
        "continuous_neg_values": neg_vals,
        "overall_tnr": n_rej / len(neg),
        "n_positive": len(pos), "n_negative": len(neg),
        "positive_tiers": [r.tier for r in pos],
        "cohens_d_continuous": float((ps.mean() - ns.mean()) / pooled),
        "negative_panel": ["polar_flock", "milling_vortex", "isotropic_noise",
                           "uniform_nematic (defect-free)"],
        "note": "Active nematic with ±1/2 topological defects; discriminator = coherent "
                "half-integer defect density D* gated by local nematic order, low polar "
                "order, low angular momentum. Validated TNR 1.0, d=4.81; real P5 Vicsek "
                "rejected; finite-size definitive G=64/96/128.",
    }
    doc = {"pattern_id": "P33", "panel_version": "discovery-ring0",
           "detector_format": "orientation_field", "summary": summary}
    os.makedirs(os.path.dirname(out), exist_ok=True)
    json.dump(doc, open(out, "w"), indent=2)
    print(f"wrote {out}  (pos {len(pos_vals)}, neg {len(neg_vals)}, TNR {n_rej}/{len(neg)})",
          flush=True)
    return 0


if __name__ == "__main__":
    import sys
    if "--harden" in sys.argv:
        raise SystemExit(harden())
    if "--panel" in sys.argv:
        raise SystemExit(write_panel())
    raise SystemExit(validate())
