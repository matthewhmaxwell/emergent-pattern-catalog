"""Out-of-TNR diagnostics for the P33-P37 full panels (v12).

1. SKIPPED Class A substrates. Where a detector declares a primary-metric
   invariance (epc/phase2a/detector_invariance.py), run_panel skips the
   corresponding positive-derived Class A substrate (spec v1.2). The flags were
   set from the metric definitions before the panels ran; this script still
   scores every skipped substrate (5 generator seeds against canonical positive
   seed 0, exactly as the harness would build it) so the write-up can say what
   the detector does on them. These results are NOT part of any TNR.
2. P35 box handling. run_p35 scores each substrate in its own periodic box
   (frames' box_size). This re-scores the P35 continuous_2d catalog mates with the
   harness default instead (canonical metadata, L=14 for every substrate) to show
   whether that choice changes any verdict.
3. Every negative that FIRED in a v12 panel JSON is rebuilt exactly as the
   harness built it (same generator / seed / adapter) and re-scored, recording
   the detector's secondary metrics so the write-up can say which gate let it
   through. Deterministic: the verdict must match the panel JSON (checked).

Run:  PYTHONPATH=. python analysis/validation_rebuild/p33_p37_panel_diagnostics.py
Writes analysis/outputs/p33_p37_panel_diagnostics.json
"""
from __future__ import annotations

import json
import os
from typing import Any, Dict

import analysis.run_phase2a_panel as R
from epc.phase2a import catalog as C
from epc.phase2a import structured as ST
from epc.phase2a import synthetic as S
from epc.phase2a.detector_invariance import get_flags

OUT = "analysis/outputs/p33_p37_panel_diagnostics.json"

SPEC = {
    "P33": (R.build_p33_positives, R.make_p33_detector_fn, "orientation_field"),
    "P34": (R.build_p34_positives, R.make_p34_detector_fn, "adjacency_timeseries"),
    "P35": (R.build_p35_positives, lambda: R._p35_substrate_box(R.make_p35_detector_fn()), "positions"),
    "P36": (R.build_p36_positives, R.make_p36_detector_fn, "scalar_wealth"),
    "P37": (R.build_p37_positives, R.make_p37_detector_fn, "abundance_resource_timeseries"),
}


def _num(v):
    try:
        return float(v)
    except (TypeError, ValueError):
        return v


def _row(r) -> Dict[str, Any]:
    return {"detected": bool(r.detected), "tier": str(r.tier), "confidence": float(r.confidence),
            "primary_metric": {k: _num(v) for k, v in r.primary_metric.items()},
            "secondary_metrics": {k: _num(v) for k, v in (r.secondary_metrics or {}).items()}}


FAILED = {"P33": R.p33_failed, "P34": R.p34_failed, "P35": R.p35_failed,
          "P36": R.p36_failed, "P37": R.p37_failed}


def fired_negatives(out: Dict[str, Any]) -> None:
    out["fired_negatives"] = {}
    for pid, (bp, mk, fmt) in SPEC.items():
        path = R.P33_P37_OUT[pid]
        if not os.path.exists(path):
            continue
        J = json.load(open(path))
        pos, md = bp(n_seeds=1)
        det = mk()
        rows = []
        for r in J["synthetic"]:
            if r.get("detected"):
                sid = r["substrate"]
                g = S.SYNTHETIC_GENERATORS[sid]
                h = (g(fmt, r["seed"], positive=pos[0]) if sid in ("permutation_shuffled", "time_shuffled")
                     else g(fmt, r["seed"], n_steps=200))
                rows.append(dict(cls="A", substrate=sid, seed=r["seed"], panel_tier=r["verdict"], **_row(det(h, md))))
        for r in J["catalog"]:
            if r.get("detected"):
                if r["kind"] == "catalog_mate":
                    h = C.load_catalog_substrate_for_format(r["substrate"], fmt, target_steps=200,
                                                            target_shape=(32, 32), target_n=300)
                else:
                    h = ST.SUPPLEMENT_BUILDERS[r["substrate"]](r["seed"])
                rows.append(dict(cls="B", substrate=r["substrate"], kind=r["kind"], seed=r.get("seed"),
                                 panel_tier=r["verdict"], **_row(det(h, md))))
        for r in J["failed_regime"]:
            if r.get("detected"):
                reg = next(g for g in FAILED[pid].CONFIG["regimes"]
                           if g["label"] == r["substrate"] and g["seed"] == r["seed"])
                rows.append(dict(cls="C", substrate=r["substrate"], seed=r["seed"], panel_tier=r["verdict"],
                                 **_row(det(FAILED[pid].build_substrate(reg), md))))
        for x in rows:
            x["reproduces_panel_tier"] = (x["tier"] == x["panel_tier"])
            print(pid, x["cls"], x["substrate"], x.get("seed"), x["tier"], x["secondary_metrics"], flush=True)
        out["fired_negatives"][pid] = rows


def main() -> int:
    out: Dict[str, Any] = {"skipped_class_a": {}, "p35_canonical_box_mates": {}}
    for pid, (bp, mk, fmt) in SPEC.items():
        f = get_flags(pid)
        skipped = ([("permutation_shuffled", "primary_metric_permutation_invariant")]
                   if f.permutation_invariant else []) + \
                  ([("time_shuffled", "primary_metric_time_shuffle_invariant")]
                   if f.time_shuffle_invariant else [])
        if not skipped:
            continue
        pos, md = bp(n_seeds=1)
        det = mk()
        for sid, why in skipped:
            rows = []
            for s in range(5):
                h = S.SYNTHETIC_GENERATORS[sid](fmt, s, positive=pos[0])
                rows.append(dict(seed=s, **_row(det(h, md))))
            out["skipped_class_a"][f"{pid}:{sid}"] = {"skip_reason": why, "runs": rows,
                                                     "n_detected": sum(r["detected"] for r in rows)}
            print(pid, sid, [(r["detected"], r["tier"]) for r in rows], flush=True)

    # P35 catalog mates scored with the harness-default canonical box (L=14).
    pos, md = R.build_p35_positives(n_seeds=1)
    raw = R.make_p35_detector_fn()
    own = R._p35_substrate_box(raw)
    for m in C.class_b_for_pattern("P35")["catalog_mates"]:
        h = C.load_catalog_substrate_for_format(m, "positions", target_steps=200)
        out["p35_canonical_box_mates"][m] = {
            "own_box": float(h[0]["box_size"]),
            "scored_with_own_box": _row(own(h, md)),
            "scored_with_canonical_L14": _row(raw(h, md)),
        }
        print("P35 mate", m, out["p35_canonical_box_mates"][m], flush=True)

    fired_negatives(out)

    os.makedirs(os.path.dirname(OUT), exist_ok=True)
    json.dump(out, open(OUT, "w"), indent=2, default=str)
    print(f"wrote {OUT}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
