"""Sprint 24 Phase 1 — synthetic candidate-fix dry-run grader.

Loads the baseline JSON from `docs/sprint24/baseline_voter_schelling.json`
and applies six candidate-fix rules to the saved metrics. Predicts the
post-fix tier for each run WITHOUT modifying any detector code. This
allows grading multiple candidate fixes against one shared baseline in
a single pass — much cheaper than the Sprint 21 approach of
characterizing each candidate against a freshly modified detector.

Run from repository root:
  PYTHONPATH=. python scripts/sprint24_grade_candidates.py

Output: docs/sprint24/candidate_grades.json + a printed summary table.

Six candidates evaluated:

  C1: tighten CONFIRMATION_WALL_FINAL_MAX 0.30 → 0.25.
      Schelling thr ∈ {0.43, 0.5} has wall_final ≈ 0.275, would
      screen-fail. Voter wall_final at L=128 max = 0.238 (margin
      +0.012 — thin).

  C2: P1-aware definitive downgrade. Modify _check_definitive to
      reject when exclusion_results['P1'] != 'excluded'. Schelling
      lacks an 'update' metadata key so P1=inconclusive, downgrades
      to CONFIRMATION. Voter has update='asynchronous_copy_neighbor'
      so P1=excluded, unaffected.

  C3: Schelling registry update token. Add metadata 'update'='move'.
      Standalone, this is a NO-OP because the Sprint 23 detector's
      _check_definitive does not consult exclusion results — the
      bonus dict's all_exclusions_cleared flag is hardcoded to True
      regardless. C3 is meaningful only as plumbing for C2.

  C4 = C1 + C2.

  C5: raise DEFINITIVE_MORAN_FINAL_MIN from 0.30 to 0.45. Cleanest
      empirical separation: voter moran_final ∈ [0.499, 0.663] vs
      Schelling thr ∈ {0.43, 0.5} ∈ [0.375, 0.410] — a 0.089-wide
      gap with ~0.05 margin both sides of 0.45.

  C6 = C5 + C2.  *** RECOMMENDED, IMPLEMENTED IN SPRINT 24 ***
"""
from __future__ import annotations

import json
from pathlib import Path
from typing import Any

REPO_ROOT = Path(__file__).resolve().parent.parent
BASELINE_PATH = REPO_ROOT / "docs" / "sprint24" / "baseline_voter_schelling.json"
OUTPUT_PATH = REPO_ROOT / "docs" / "sprint24" / "candidate_grades.json"

# Definitive metric gates that the candidates do NOT modify (used to
# infer the metric-only post-candidate outcome for boundary cases).
DEFINITIVE_MORAN_FINAL_MAX = 0.75
DEFINITIVE_WALL_FINAL_MIN = 0.05
DEFINITIVE_MINORITY_FINAL_MIN = 0.05

# Candidate constants
C1_WALL_FINAL_MAX = 0.25  # tightened from 0.30
C5_MORAN_FINAL_MIN = 0.45  # raised from 0.30


def apply_candidate(run: dict, candidate: str) -> str:
    """Compute the tier outcome under a candidate fix.

    Operates on a frozen detector state: takes the run's already-
    computed primary_metric, secondary_metrics, exclusion_results,
    null_p, etc., and re-derives the tier under the modified gate.

    Tier ladder:
      BELOW_SCREENING (screening primary fails)
      SCREENING        (passes screening, fails confirmation)
      CONFIRMATION     (passes confirmation, fails definitive)
      DEFINITIVE       (passes definitive)
    """
    baseline = run["tier"]
    pm = run.get("primary_metric", {})
    sm = run.get("secondary_metrics", {})
    excl = run.get("exclusion_results", {})

    if baseline == "BELOW_SCREENING":
        return "BELOW_SCREENING"

    moran_final = pm.get("moran_final_qtr_mean", 0.0)
    wall_final = sm.get("wall_final_qtr_mean", 0.0)
    p1_excl = excl.get("P1") or "inconclusive"

    def passes_confirmation(wall_max: float) -> bool:
        return wall_final < wall_max

    if candidate == "BASELINE":
        return baseline

    if candidate == "C1":
        if baseline in ("CONFIRMATION", "DEFINITIVE"):
            if not passes_confirmation(C1_WALL_FINAL_MAX):
                return "SCREENING"
        return baseline

    if candidate == "C2":
        if baseline == "DEFINITIVE" and p1_excl != "excluded":
            return "CONFIRMATION"
        return baseline

    if candidate == "C3":
        # Standalone: no effect (definitive doesn't consult exclusions
        # in the Sprint 23 detector).
        return baseline

    if candidate == "C2_AND_C3":
        if baseline == "DEFINITIVE" and p1_excl != "excluded":
            return "CONFIRMATION"
        return baseline

    if candidate == "C4_C1_AND_C2":
        new_tier = baseline
        if baseline in ("CONFIRMATION", "DEFINITIVE"):
            if not passes_confirmation(C1_WALL_FINAL_MAX):
                new_tier = "SCREENING"
        if new_tier == "DEFINITIVE" and p1_excl != "excluded":
            new_tier = "CONFIRMATION"
        return new_tier

    if candidate == "C5":
        if baseline == "DEFINITIVE":
            if moran_final < C5_MORAN_FINAL_MIN:
                return "CONFIRMATION"
        return baseline

    if candidate == "C6_C5_AND_C2":
        new_tier = baseline
        if new_tier == "DEFINITIVE" and moran_final < C5_MORAN_FINAL_MIN:
            new_tier = "CONFIRMATION"
        if new_tier == "DEFINITIVE" and p1_excl != "excluded":
            new_tier = "CONFIRMATION"
        return new_tier

    raise ValueError(f"Unknown candidate: {candidate}")


def grade_runs(runs: list[dict], candidate: str, label: str) -> dict:
    counts: dict[str, int] = {}
    detail = []
    for r in runs:
        new = apply_candidate(r, candidate)
        counts[new] = counts.get(new, 0) + 1
        detail.append({
            "model": r.get("model"),
            "L": r.get("L"),
            "seed": r.get("seed"),
            "threshold": r.get("threshold"),
            "metadata_path": r.get("metadata_path"),
            "baseline_tier": r["tier"],
            "candidate_tier": new,
        })
    return {"label": label, "candidate": candidate, "counts": counts,
            "n": len(runs), "detail": detail}


def fmt_grade(grade: dict) -> str:
    counts = grade["counts"]
    parts = [f"{k}={counts.get(k, 0)}" for k in
             ["DEFINITIVE", "CONFIRMATION", "SCREENING", "BELOW_SCREENING"]]
    return "  ".join(parts) + f"  (n={grade['n']})"


def main() -> None:
    if not BASELINE_PATH.exists():
        raise SystemExit(
            f"Baseline JSON not found at {BASELINE_PATH}. "
            "Run scripts/sprint24_baseline.py first."
        )
    with open(BASELINE_PATH) as f:
        data = json.load(f)
    voter = data["voter_runs"]
    sch = data["schelling_runs"]
    sch_with_meta = [r for r in sch if r.get("metadata_path") == "with_metadata"]

    print("=" * 78)
    print("SPRINT 24 PHASE 1 — synthetic candidate-fix dry-run grading")
    print("=" * 78)
    print(f"\nBaseline: {BASELINE_PATH}")
    print(f"Voter runs: {len(voter)} (canonical positive)")
    print(f"Schelling runs (with_metadata): {len(sch_with_meta)}")

    candidates = ["BASELINE", "C1", "C2", "C3", "C2_AND_C3",
                  "C4_C1_AND_C2", "C5", "C6_C5_AND_C2"]

    print("\n" + "-" * 78)
    print("VOTER (must hold DEFINITIVE — any regression is a FAIL)")
    print("-" * 78)
    voter_grades = {}
    for cand in candidates:
        g = grade_runs(voter, cand, "voter")
        voter_grades[cand] = g
        print(f"  {cand:<14} {fmt_grade(g)}")

    print("\n" + "-" * 78)
    print("SCHELLING (with metadata, by threshold) — must NOT reach DEFINITIVE")
    print("-" * 78)
    sch_grades_by_thr: dict[float, dict] = {}
    for thr in [0.30, 0.375, 0.43, 0.5]:
        runs_thr = [r for r in sch_with_meta if r.get("threshold") == thr]
        print(f"\n  threshold = {thr}  (n={len(runs_thr)} seeds)")
        for cand in candidates:
            g = grade_runs(runs_thr, cand, f"sch_thr={thr}")
            sch_grades_by_thr.setdefault(thr, {})[cand] = g
            marker = ""
            if g["counts"].get("DEFINITIVE", 0) > 0:
                marker = " [FALSE POSITIVE]"
            print(f"    {cand:<14} {fmt_grade(g)}{marker}")

    print("\n" + "-" * 78)
    print("VOTER wall_final by L (C1 risk assessment, candidate ceiling = 0.25)")
    print("-" * 78)
    by_L: dict[int, list[float]] = {}
    for r in voter:
        L = r.get("L")
        w = r.get("secondary_metrics", {}).get("wall_final_qtr_mean", float("nan"))
        by_L.setdefault(L, []).append(w)
    for L in sorted(by_L):
        ws = by_L[L]
        print(f"  L={L:<4}: min={min(ws):.3f}  max={max(ws):.3f}  "
              f"mean={sum(ws)/len(ws):.3f}  margin_to_0.25={0.25 - max(ws):+.3f}")

    print("\n" + "-" * 78)
    print("VOTER vs SCHELLING moran_final (C5 separation, candidate floor = 0.45)")
    print("-" * 78)
    voter_morans = [r["primary_metric"]["moran_final_qtr_mean"] for r in voter]
    sch_target_morans = [
        r["primary_metric"]["moran_final_qtr_mean"]
        for r in sch_with_meta
        if r.get("threshold") in (0.43, 0.5)
    ]
    print(f"  voter (n={len(voter_morans)}): "
          f"min={min(voter_morans):.3f} max={max(voter_morans):.3f}")
    print(f"  schelling thr in {{0.43, 0.5}} (n={len(sch_target_morans)}): "
          f"min={min(sch_target_morans):.3f} max={max(sch_target_morans):.3f}")
    gap = min(voter_morans) - max(sch_target_morans)
    midpoint = (min(voter_morans) + max(sch_target_morans)) / 2
    print(f"  separation gap: {gap:.3f}  midpoint: {midpoint:.3f}  "
          f"voter_margin_to_0.45: {min(voter_morans) - 0.45:+.3f}  "
          f"schelling_margin_to_0.45: {0.45 - max(sch_target_morans):+.3f}")

    grades_out = {
        "schema_version": 1,
        "generated_by": "scripts/sprint24_grade_candidates.py",
        "baseline_source": str(BASELINE_PATH.relative_to(REPO_ROOT)),
        "candidate_constants": {
            "C1_WALL_FINAL_MAX": C1_WALL_FINAL_MAX,
            "C5_MORAN_FINAL_MIN": C5_MORAN_FINAL_MIN,
        },
        "voter": voter_grades,
        "schelling_by_threshold": {str(k): v for k, v in
                                   sch_grades_by_thr.items()},
        "voter_wall_final_by_L": {str(L): {"values": by_L[L],
                                           "min": min(by_L[L]),
                                           "max": max(by_L[L]),
                                           "mean": sum(by_L[L]) / len(by_L[L])}
                                  for L in by_L},
    }
    with open(OUTPUT_PATH, "w") as f:
        json.dump(grades_out, f, indent=2, default=str)
    print(f"\nSaved {OUTPUT_PATH}")


if __name__ == "__main__":
    main()
