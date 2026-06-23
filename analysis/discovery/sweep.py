"""Discovery sweep (Ring 0) — run out-of-catalog emergent systems through the
comprehensive instrument and surface EMERGENT-UNCLASSIFIED candidates, each tagged
with the firing channel and its closest catalog pattern, for vetting toward P33+.

OOD operating point: match_min_tier="confirmation" — the T2c-recommended gate for
UNKNOWN systems. A screening-only firing is demoted to the emergence-driven verdict
(surfaced as `demoted_to`), so a MATCH here is a high-confidence over-claim worth
investigating, not screening noise.

Per-system readout:
  verdict        MATCH / EMERGENT-UNCLASSIFIED / NO-EMERGENCE
  emergence      generic-emergence score (0-1)
  channel        which emergence channel fired (em["kind"]): morphology kind,
                 orientation, synergy(psi_ce:…), temporal(spectral-peak),
                 heavy-tail(power-law), network(modularity/scale-free)
  closest        the top-ranked catalog detector + its tier (the nearest neighbour)
  demoted_to     a screening-only firing that the confirmation gate demoted

Usage:
    python -m analysis.discovery.sweep --seeds 3
    python -m analysis.discovery.sweep --seeds 1            # quick smoke test
"""
from __future__ import annotations

import argparse
import json
import os
import time
from collections import Counter
from typing import Any, Dict, List, Optional

from analysis.battery_profile import build_detector_fns, profile_observation
from epc.phase2a.battery import Battery


def run_sweep(systems: List[Dict[str, Any]], seeds: int = 3,
              match_min_tier: str = "confirmation") -> List[Dict[str, Any]]:
    battery = Battery.load()
    fns = build_detector_fns()
    rows: List[Dict[str, Any]] = []
    print(f"battery={len(battery)} detectors={len(fns)} systems={len(systems)} "
          f"seeds={seeds} gate={match_min_tier}", flush=True)
    for s in systems:
        reps = seeds if s.get("stochastic", True) else 1
        for rep in range(reps):
            t0 = time.time()
            try:
                hist, meta = s["build"](rep)
                out = profile_observation(hist, meta, battery=battery,
                                          detector_fns=fns, match_min_tier=match_min_tier)
                v = out["verdict"]; em = out["emergence"]; top = out.get("top") or {}
                demoted = v.get("demoted_match") or {}
                row = {"name": s["name"], "family": s.get("family"), "ref": s.get("ref"),
                       "seed": rep, "verdict": v.get("verdict"),
                       "matched": v.get("pattern_id"),
                       "emergence": em.get("score"), "channel": em.get("kind"),
                       "order_gain": em.get("order_gain"),
                       "closest": top.get("pattern_id"),
                       "closest_tier": top.get("detector_tier"),
                       "demoted_to": demoted.get("pattern_id"),
                       "n_frames": len(hist), "secs": round(time.time() - t0, 1)}
            except Exception as ex:
                row = {"name": s["name"], "family": s.get("family"), "seed": rep,
                       "verdict": "ERROR", "error": f"{type(ex).__name__}: {ex}",
                       "secs": round(time.time() - t0, 1)}
            rows.append(row)
            print(f"  {row['name']:<26} s{rep} -> {str(row.get('verdict')):<22} "
                  f"em={row.get('emergence')} ch={str(row.get('channel')):<26} "
                  f"closest={row.get('closest')}/{row.get('closest_tier')} "
                  f"demoted={row.get('demoted_to')} ({row['secs']}s)", flush=True)
    return rows


def _agg(rows: List[Dict[str, Any]]) -> List[str]:
    names: List[str] = []
    for r in rows:
        if r["name"] not in names:
            names.append(r["name"])
    L = ["| system | family | verdict(s) | mean em | channel(s) | closest | ref |",
         "|---|---|---|---|---|---|---|"]
    for nm in names:
        rs = [r for r in rows if r["name"] == nm]
        vc = Counter(r["verdict"] for r in rs)
        vstr = ", ".join(f"{k}×{v}" if v > 1 else k for k, v in vc.most_common())
        ems = [r["emergence"] for r in rs if isinstance(r.get("emergence"), (int, float))]
        emm = round(sum(ems) / len(ems), 3) if ems else "-"
        chs = Counter(r.get("channel") for r in rs if r.get("channel"))
        chstr = ", ".join(k for k, _ in chs.most_common())
        cl = Counter(f"{r.get('closest')}/{r.get('closest_tier')}" for r in rs
                     if r.get("closest"))
        clstr = ", ".join(k for k, _ in cl.most_common(2))
        ref = rs[0].get("ref") or "—"
        fam = rs[0].get("family") or "—"
        L.append(f"| {nm} | {fam} | {vstr} | {emm} | {chstr} | {clstr} | {ref} |")
    return L


def summarize(rows: List[Dict[str, Any]]) -> Dict[str, Any]:
    cand = [r for r in rows if r["family"] != "null"]
    nulls = [r for r in rows if r["family"] == "null"]
    by_verdict = Counter(r["verdict"] for r in cand)
    null_verdict = Counter(r["verdict"] for r in nulls)
    false_match = [r for r in cand if r["verdict"] == "MATCH"]
    abstained = [r for r in cand if r["verdict"] == "EMERGENT-UNCLASSIFIED"]
    missed = [r for r in cand if r["verdict"] == "NO-EMERGENCE"]
    null_leak = [r for r in nulls if r["verdict"] != "NO-EMERGENCE"]
    return {"candidate_verdicts": dict(by_verdict),
            "null_verdicts": dict(null_verdict),
            "n_abstained": len(abstained), "n_false_match": len(false_match),
            "n_missed_emergence": len(missed), "n_null_leak": len(null_leak),
            "false_match": [(r["name"], r["matched"], r["seed"]) for r in false_match],
            "missed": [(r["name"], r["seed"]) for r in missed],
            "null_leak": [(r["name"], r["verdict"], r["seed"]) for r in null_leak]}


def write_report(rows: List[Dict[str, Any]], summary: Dict[str, Any],
                 seeds: int, gate: str, md_path: str) -> None:
    L: List[str] = []
    L.append("# Discovery Sweep — Ring 0 (out-of-catalog emergent systems)\n")
    L.append(f"**Seeds/system:** {seeds}. **OOD gate:** match_min_tier=`{gate}`. "
             "Candidates are canonical emergent phenomena absent from the 32-pattern "
             "catalog; the target read is EMERGENT-UNCLASSIFIED (emergence detected, "
             "no false MATCH) with a named firing channel. Null controls must read "
             "NO-EMERGENCE.\n")
    L.append("## Outcome summary\n")
    L.append(f"- Candidate verdicts: `{summary['candidate_verdicts']}`")
    L.append(f"- EMERGENT-UNCLASSIFIED (correct abstention): **{summary['n_abstained']}**")
    L.append(f"- false MATCH (over-claim — instrument finding): **{summary['n_false_match']}** "
             f"{summary['false_match']}")
    L.append(f"- NO-EMERGENCE (missed — coverage gap): **{summary['n_missed_emergence']}** "
             f"{summary['missed']}")
    L.append(f"- Null controls: `{summary['null_verdicts']}` — leaks: {summary['null_leak']}\n")
    L.append("## Candidates\n")
    L.extend(_agg([r for r in rows if r["family"] != "null"]))
    L.append("\n## Null controls\n")
    L.extend(_agg([r for r in rows if r["family"] == "null"]))
    errs = [r for r in rows if r["verdict"] == "ERROR"]
    if errs:
        L.append("\n## Errors\n")
        for r in errs:
            L.append(f"- `{r['name']}` s{r['seed']}: {r.get('error')}")
    open(md_path, "w").write("\n".join(L) + "\n")
    print(f"\nwrote {md_path}", flush=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", type=int, default=3)
    ap.add_argument("--gate", type=str, default="confirmation")
    ap.add_argument("--json", type=str, default="analysis/outputs/ring0_discovery.json")
    ap.add_argument("--md", type=str,
                    default="docs/validation_rebuild/ring0_discovery_sweep.md")
    args = ap.parse_args()
    from analysis.discovery.ring0_systems import ALL
    t0 = time.time()
    rows = run_sweep(ALL, seeds=args.seeds, match_min_tier=args.gate)
    summary = summarize(rows)
    os.makedirs("analysis/outputs", exist_ok=True)
    json.dump({"rows": rows, "summary": summary, "seeds": args.seeds, "gate": args.gate},
              open(args.json, "w"), indent=2, default=str)
    write_report(rows, summary, args.seeds, args.gate, args.md)
    print(f"\n=== Ring-0 sweep: {summary['n_abstained']} abstained, "
          f"{summary['n_false_match']} false-MATCH, {summary['n_missed_emergence']} missed, "
          f"{summary['n_null_leak']} null-leak  (total {round(time.time() - t0, 1)}s)",
          flush=True)


if __name__ == "__main__":
    main()
