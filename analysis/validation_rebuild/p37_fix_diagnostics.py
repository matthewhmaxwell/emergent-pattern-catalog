"""Per-negative detail for the P37 post-panel-fix panel (v2 detector).

Rebuilds every negative of analysis/outputs/p37_resource_competition_phase2a_panel_v12_fixed.json
exactly as run_panel built it (same generator, seed, adapter), re-scores it with
the v2 detector, and records the secondary metrics (oscillation p/z, coupling p,
persistence, amplitudes) plus the v1 score, so the write-up can say which test
rejected or passed each one. Also scores the Class A permutation_shuffled
substrate the harness skips (5 seeds), out of TNR.

Run:  PYTHONPATH=. python analysis/validation_rebuild/p37_fix_diagnostics.py
Writes analysis/outputs/p37_fix_diagnostics.json
"""
from __future__ import annotations

import json

import analysis.run_phase2a_panel as R
from epc.phase2a import catalog as C
from epc.phase2a import structured as ST
from epc.phase2a import synthetic as S

FMT = "abundance_resource_timeseries"
OUT = "analysis/outputs/p37_fix_diagnostics.json"


def _row(r):
    sm = dict(r.secondary_metrics)
    return {"detected": bool(r.detected), "tier": r.tier,
            "sustained_oscillation": r.primary_metric["sustained_oscillation"],
            "v1_coexistence_oscillation": r.primary_metric["coexistence_oscillation"],
            **{k: sm.get(k) for k in ("n_persist", "n_material", "n_sustained", "median_amp",
                                      "p_osc", "z_osc", "p_coupling", "coupling_corr", "n_window",
                                      "species_p", "species_amp")}}


def main() -> int:
    J = json.load(open(R.P33_P37_OUT["P37"]))
    pos, md = R.build_p37_positives(1)
    det = R.make_p37_detector_fn()
    rows = []
    for r in J["synthetic"]:
        sid = r["substrate"]
        if r.get("verdict") == "SKIPPED-degenerate-by-construction":
            for s in range(5):
                h = S.SYNTHETIC_GENERATORS[sid](FMT, s, positive=pos[0])
                rows.append(dict(cls="A-skipped(out of TNR)", substrate=sid, seed=s, **_row(det(h, md))))
            continue
        g = S.SYNTHETIC_GENERATORS[sid]
        h = (g(FMT, r["seed"], positive=pos[0]) if sid in ("permutation_shuffled", "time_shuffled")
             else g(FMT, r["seed"], n_steps=200))
        rows.append(dict(cls="A", substrate=sid, seed=r["seed"], panel_tier=r["verdict"], **_row(det(h, md))))
    for r in J["catalog"]:
        h = (C.load_catalog_substrate_for_format(r["substrate"], FMT, target_steps=200)
             if r["kind"] == "catalog_mate" else ST.SUPPLEMENT_BUILDERS[r["substrate"]](r["seed"]))
        rows.append(dict(cls="B", substrate=r["substrate"], seed=r.get("seed"), panel_tier=r["verdict"],
                         **_row(det(h, md))))
    for r in J["failed_regime"]:
        reg = next(g for g in R.p37_failed.CONFIG["regimes"]
                   if g["label"] == r["substrate"] and g["seed"] == r["seed"])
        rows.append(dict(cls="C", substrate=r["substrate"], seed=r["seed"], panel_tier=r["verdict"],
                         **_row(det(R.p37_failed.build_substrate(reg), md))))
    for x in rows:
        print(f"{x['cls']:22s} {x['substrate']:32s} {str(x.get('seed')):>11s} {x['tier']:12s} "
              f"persist={x['n_persist']} material={x['n_material']} p_osc={x['p_osc']} z={x['z_osc']} "
              f"p_cpl={x['p_coupling']} v1={x['v1_coexistence_oscillation']}", flush=True)
    json.dump({"panel": R.P33_P37_OUT["P37"], "rows": rows}, open(OUT, "w"), indent=2, default=str)
    print("wrote", OUT)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
