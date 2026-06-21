"""T2c — OOD validation suite for the EPC instrument.

Points the full detector battery (`analysis/battery_profile.profile_observation`)
at held-out systems it was never built on, and reports the instrument's
discrimination as precision/recall-style rates:

  recognition recall   independent re-implementations of catalog phenomena that
                       fire the CORRECT pattern (true positive).
  novelty abstention   out-of-catalog emergent systems flagged EMERGENT-UNCLASSIFIED
                       (correct "none-of-the-above").
  null specificity     non-emergent systems flagged NO-EMERGENCE.
  false-MATCH rate     novelty/null systems the catalog WRONGLY claims as a pattern
                       (over-claiming — the cardinal instrument sin).
  false-novel rate     null systems wrongly flagged emergent.

PROBE systems expose known emergence-indicator coverage gaps and are reported
separately, NOT folded into the headline rates.

Usage:
    python -m analysis.t2c_ood_suite --seeds 3
    python -m analysis.t2c_ood_suite --seeds 1 --only null,novelty   # quick subset
"""
from __future__ import annotations

import argparse
import json
import time
from collections import Counter
from typing import Any, Dict, List, Optional

from analysis.battery_profile import build_detector_fns, profile_observation
from analysis.t2c_systems import SYSTEMS
from epc.phase2a.battery import Battery


def _classify(cls: str, verdict: str, matched: Optional[str],
              top_pid: Optional[str], expect: Optional[str]) -> str:
    """Map a single run's verdict to a per-class outcome label."""
    if cls == "recognition":
        if verdict == "MATCH" and matched == expect:
            return "recognized"
        if verdict == "MATCH":
            return "wrong_match"
        if top_pid == expect:
            return "top1_only"          # ranked correct #1 but below firing threshold
        return "missed"
    if cls == "novelty":
        if verdict == "EMERGENT-UNCLASSIFIED":
            return "abstained"
        if verdict == "MATCH":
            return "false_match"
        return "missed_emergence"
    if cls == "null":
        if verdict == "NO-EMERGENCE":
            return "correct"
        if verdict == "EMERGENT-UNCLASSIFIED":
            return "false_novel"
        return "false_match"
    return verdict  # probe: pass through raw verdict


def run_suite(seeds: int = 3, only: Optional[List[str]] = None) -> Dict[str, Any]:
    battery = Battery.load()
    fns = build_detector_fns()
    systems = [s for s in SYSTEMS if (only is None or s["cls"] in only)]
    rows: List[Dict[str, Any]] = []
    print(f"battery={len(battery)} detectors={len(fns)} systems={len(systems)} seeds={seeds}",
          flush=True)
    for s in systems:
        reps = seeds if s.get("stochastic", True) else 1
        for rep in range(reps):
            t0 = time.time()
            try:
                hist, meta = s["build"](rep)
                out = profile_observation(hist, meta, battery=battery, detector_fns=fns)
                v = out["verdict"]
                verdict = v.get("verdict")
                matched = v.get("pattern_id")
                top = out.get("top") or {}
                top_pid = top.get("pattern_id")
                em = out["emergence"]["score"]
                outcome = _classify(s["cls"], verdict, matched, top_pid, s.get("expect"))
                row = {"name": s["name"], "cls": s["cls"], "seed": rep,
                       "expect": s.get("expect"), "verdict": verdict,
                       "matched": matched, "top_pid": top_pid,
                       "top_conf": top.get("calibrated_confidence"),
                       "top_tier": top.get("detector_tier"),
                       "emergence": em, "outcome": outcome,
                       "n_frames": len(hist), "secs": round(time.time() - t0, 1)}
            except Exception as ex:
                row = {"name": s["name"], "cls": s["cls"], "seed": rep,
                       "expect": s.get("expect"), "verdict": "ERROR",
                       "outcome": "error", "error": f"{type(ex).__name__}: {ex}",
                       "secs": round(time.time() - t0, 1)}
            rows.append(row)
            print(f"  {row['name']:<28} s{rep} -> {row.get('verdict'):<22} "
                  f"{row.get('outcome'):<16} em={row.get('emergence')} "
                  f"({row['secs']}s)", flush=True)
    return {"rows": rows, "rates": _rates(rows), "seeds": seeds}


def _rate(num: int, den: int) -> Optional[float]:
    return round(num / den, 3) if den else None


_STRONG = {"confirmation", "definitive"}


def _strict_outcome(r: Dict[str, Any], emergence_threshold: float = 0.5) -> str:
    """Re-score a run under a tier-gated MATCH rule: a MATCH only counts if the
    firing detector reached >= confirmation tier. A screening-only detection is
    demoted to the emergence-driven verdict (EMERGENT-UNCLASSIFIED / NO-EMERGENCE).
    This isolates how much of the raw false-MATCH rate is screening-tier noise vs
    a genuine high-confidence over-claim, without mutating the instrument keystone."""
    verdict = r.get("verdict")
    if verdict == "ERROR":
        return "error"
    em = r.get("emergence") or 0.0
    strong = verdict == "MATCH" and r.get("top_tier") in _STRONG
    if strong:
        sv = "MATCH"
    elif em >= emergence_threshold:
        sv = "EMERGENT-UNCLASSIFIED"
    else:
        sv = "NO-EMERGENCE"
    return _classify(r["cls"], sv, r.get("matched") if strong else None,
                     r.get("top_pid"), r.get("expect"))


def _rates_from(rows: List[Dict[str, Any]], getter) -> Dict[str, Any]:
    def oc(cls: str) -> Counter:
        return Counter(getter(r) for r in rows if r["cls"] == cls)
    rec, nov, nul = oc("recognition"), oc("novelty"), oc("null")
    n_rec = sum(rec.values()); n_nov = sum(nov.values()); n_nul = sum(nul.values())
    nn_false_match = nov["false_match"] + nul["false_match"]
    return {
        "recognition_recall": _rate(rec["recognized"], n_rec),
        "recognition_top1": _rate(rec["recognized"] + rec["top1_only"], n_rec),
        "novelty_abstention": _rate(nov["abstained"], n_nov),
        "null_specificity": _rate(nul["correct"], n_nul),
        "false_match_rate": _rate(nn_false_match, n_nov + n_nul),
        "false_novel_rate": _rate(nul["false_novel"], n_nul),
        "counts": {"recognition": dict(rec), "novelty": dict(nov), "null": dict(nul)},
        "n": {"recognition": n_rec, "novelty": n_nov, "null": n_nul},
    }


def _rates(rows: List[Dict[str, Any]]) -> Dict[str, Any]:
    return {
        "raw": _rates_from(rows, lambda r: r["outcome"]),
        "strict_confirmation": _rates_from(rows, _strict_outcome),
    }


def _agg_table(rows: List[Dict[str, Any]], cls: str) -> List[str]:
    names = [s["name"] for s in SYSTEMS if s["cls"] == cls]
    lines = ["| system | expect | verdicts (across seeds) | outcome | mean em |",
             "|---|---|---|---|---|"]
    for nm in names:
        rs = [r for r in rows if r["name"] == nm]
        if not rs:
            continue
        vc = Counter(r["verdict"] for r in rs)
        vstr = ", ".join(f"{k}×{v}" if v > 1 else k for k, v in vc.most_common())
        oc = Counter(r["outcome"] for r in rs).most_common(1)[0][0]
        ems = [r["emergence"] for r in rs if isinstance(r.get("emergence"), (int, float))]
        emm = round(sum(ems) / len(ems), 3) if ems else "-"
        exp = rs[0].get("expect") or "—"
        # surface which pattern a false MATCH landed on
        mm = Counter(r["matched"] for r in rs if r.get("matched"))
        if mm:
            vstr += f"  [→ {', '.join(k for k, _ in mm.most_common())}]"
        lines.append(f"| {nm} | {exp} | {vstr} | {oc} | {emm} |")
    return lines


def _rate_rows(rt: Dict[str, Any]) -> List[str]:
    return [
        "| metric | value | reading |",
        "|---|---|---|",
        f"| recognition recall | {rt['recognition_recall']} | alt-model fires the correct pattern (MATCH) |",
        f"| recognition top-1 | {rt['recognition_top1']} | correct pattern ranked #1 (incl. below-threshold) |",
        f"| novelty abstention | {rt['novelty_abstention']} | out-of-catalog emergence → EMERGENT-UNCLASSIFIED |",
        f"| null specificity | {rt['null_specificity']} | non-emergent → NO-EMERGENCE |",
        f"| **false-MATCH rate** | **{rt['false_match_rate']}** | novelty/null wrongly claimed as a pattern (lower=better) |",
        f"| false-novel rate | {rt['false_novel_rate']} | null wrongly flagged emergent (lower=better) |",
    ]


def write_report(result: Dict[str, Any], md_path: str) -> None:
    raw = result["rates"]["raw"]
    strict = result["rates"]["strict_confirmation"]
    L: List[str] = []
    L.append("# T2c — OOD Validation Suite (instrument generalization)\n")
    L.append(f"**Seeds per system:** {result['seeds']}. Battery pointed at held-out "
             "systems via `analysis/battery_profile.profile_observation`. Recognition "
             "systems are independent re-implementations the detectors were never tuned "
             "on; novelty systems are emergent but out-of-catalog; null systems are "
             "non-emergent. Probes (reported last) expose known emergence-indicator "
             "coverage gaps and are excluded from the headline rates.\n")
    L.append("## Headline rates — raw (MATCH at any tier, current keystone)\n")
    L.extend(_rate_rows(raw))
    L.append(f"\nCounts: {json.dumps(raw['counts'])}\n")
    L.append("## Headline rates — strict (MATCH requires ≥ confirmation tier)\n")
    L.append("Screening-tier detections are demoted to the emergence-driven verdict. "
             "This is a pure instrument-layer gate (no detector is changed) and shows "
             "how much raw false-MATCH is low-rigor screening noise.\n")
    L.extend(_rate_rows(strict))
    L.append(f"\nCounts: {json.dumps(strict['counts'])}\n")
    for cls, title in [("recognition", "Recognition arm (held-out re-implementations)"),
                       ("novelty", "Novelty arm (out-of-catalog emergent)"),
                       ("null", "Null arm (non-emergent)"),
                       ("probe", "Probe arm (coverage-limit demonstrations — not scored)")]:
        L.append(f"\n## {title}\n")
        L.extend(_agg_table(result["rows"], cls))
    errs = [r for r in result["rows"] if r["outcome"] == "error"]
    if errs:
        L.append("\n## Errors\n")
        for r in errs:
            L.append(f"- `{r['name']}` s{r['seed']}: {r.get('error')}")
    open(md_path, "w").write("\n".join(L) + "\n")
    print(f"\nwrote {md_path}", flush=True)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--seeds", type=int, default=3)
    ap.add_argument("--only", type=str, default=None,
                    help="comma-separated classes: null,novelty,recognition,probe")
    ap.add_argument("--json", type=str, default="analysis/outputs/t2c_ood_results.json")
    ap.add_argument("--md", type=str, default="docs/validation_rebuild/t2c_ood_validation.md")
    args = ap.parse_args()
    only = args.only.split(",") if args.only else None
    t0 = time.time()
    result = run_suite(seeds=args.seeds, only=only)
    json.dump(result, open(args.json, "w"), indent=2, default=str)
    write_report(result, args.md)
    raw = result["rates"]["raw"]; strict = result["rates"]["strict_confirmation"]
    print(f"\nRAW    : rec={raw['recognition_recall']} nov={raw['novelty_abstention']} "
          f"null-spec={raw['null_specificity']} false-MATCH={raw['false_match_rate']}", flush=True)
    print(f"STRICT : rec={strict['recognition_recall']} nov={strict['novelty_abstention']} "
          f"null-spec={strict['null_specificity']} false-MATCH={strict['false_match_rate']} "
          f"(total {round(time.time() - t0, 1)}s)", flush=True)


if __name__ == "__main__":
    main()
