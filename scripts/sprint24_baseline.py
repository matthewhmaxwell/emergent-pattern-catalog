"""Sprint 24 Phase 1 — baseline characterization driver.

Reproducibility script for #20b recovery. Runs voter (canonical positive)
and Schelling (negative across parameter regimes) against the current
P18ConsensusDetector and saves the full per-run primary/secondary/
exclusion record to a JSON archive used by:

  - the Phase 1 candidate-fix dry-run grader
    (`scripts/sprint24_grade_candidates.py`)
  - the Phase 1 analysis document
    (`docs/sprint24/phase1_baseline.md`)
  - the Phase 2 verification step (re-run against the C6-modified
    detector, compare against Phase 1 dry-run predictions)

Workflow (look-before-touching):
  1. Voter L ∈ {64, 128, 256} × seeds {0, 1, 42, 200, 500} = 15 runs.
     n_steps follows the Sprint 20 TestSprint20SlowReplication
     convention {64: 400, 128: 400, 256: 300}.
  2. Schelling thresholds {0.30, 0.375, 0.43, 0.5} × seeds
     {0, 1, 2, 3, 4} = 20 runs at L = 64, density = 0.9, n_steps = 300.
     Each Schelling run graded twice (with realistic metadata
     `{'threshold': τ, 'density': 0.9}` and with metadata=None) for
     diagnostic comparison; the two paths produce identical tier
     outcomes.
  3. For each run: tier, primary metrics, secondary metrics, null p,
     P13/P15/P1 exclusion outcomes, effect size, warnings.
  4. Save merged baseline to `docs/sprint24/baseline_voter_schelling.json`.

Run from repository root:
  PYTHONPATH=. python scripts/sprint24_baseline.py

Wall-clock budget: ~6 minutes. Voter L=256 × 5 seeds dominates (~75s);
Schelling at 5 seeds × 4 thresholds adds ~75s; detector-side
permutation null adds ~30s.

This script is read-only against the detector and registry — it
performs no edits. Re-running it after any detector or registry
change re-derives the baseline against the new code state. The
Phase 1 candidate-fix dry-run grader operates on the saved JSON,
not on a live detector, so it remains valid for the original baseline
even after the detector is modified.
"""
from __future__ import annotations

import json
import time
from pathlib import Path
from typing import Any

import numpy as np

from epc.detector_result import DetectionTier
from epc.detectors.p18_consensus import P18ConsensusDetector
from epc.models.schelling import run_schelling
from epc.models.voter import VoterModel


# --- Configuration ----------------------------------------------------------

VOTER_LS = [64, 128, 256]
VOTER_SEEDS = [0, 1, 42, 200, 500]
VOTER_NSTEPS = {64: 400, 128: 400, 256: 300}

SCHELLING_THRESHOLDS = [0.30, 0.375, 0.43, 0.5]
SCHELLING_SEEDS = [0, 1, 2, 3, 4]
SCHELLING_L = 64
SCHELLING_NSTEPS = 300
SCHELLING_DENSITY = 0.9

# Output: docs/sprint24/baseline_voter_schelling.json (resolved relative
# to the repo root, not this script's location, so the script can be
# invoked from either)
REPO_ROOT = Path(__file__).resolve().parent.parent
OUTPUT_DIR = REPO_ROOT / "docs" / "sprint24"


# --- Helpers ----------------------------------------------------------------

def tier_label(tier: DetectionTier, detected: bool) -> str:
    """Display label combining tier and detected flag.

    The detector returns tier=SCREENING with detected=False when the
    screening primary fails. We label that case "BELOW_SCREENING" so
    the baseline JSON unambiguously distinguishes "passed screening,
    didn't reach confirmation" from "failed screening".
    """
    if not detected:
        return "BELOW_SCREENING"
    return {
        DetectionTier.SCREENING: "SCREENING",
        DetectionTier.CONFIRMATION: "CONFIRMATION",
        DetectionTier.DEFINITIVE: "DEFINITIVE",
    }[tier]


def run_voter(L: int, seed: int, n_steps: int) -> tuple[list[dict], dict]:
    """Run voter; return (history, metadata)."""
    m = VoterModel(rows=L, cols=L, seed=seed, init_mode="random",
                   neighborhood="moore", boundary="periodic")
    history = m.run(n_steps=n_steps)
    return history, m.get_metadata()


def schelling_realistic_metadata(threshold: float) -> dict:
    """Realistic Schelling metadata as supplied at runtime.

    Schelling's registered metadata_keys = ['threshold', 'density'] —
    no 'update' key. This matches what the runtime passes to the detector.
    """
    return {"threshold": threshold, "density": SCHELLING_DENSITY}


def detect_one(history: list[dict], metadata: dict | None,
               det: P18ConsensusDetector) -> dict[str, Any]:
    """Run detector + collect all relevant fields into a flat dict."""
    r = det.detect(history, model_metadata=metadata)

    # Pull trajectory metrics directly so the JSON record is self-
    # contained (some derived metrics don't appear in primary or
    # secondary on every code path).
    tm = det._trajectory_metrics(history)

    return {
        "tier": tier_label(r.tier, r.detected),
        "detected": bool(r.detected),
        "tier_value": str(r.tier.value) if hasattr(r.tier, "value") else str(r.tier),
        "confidence": float(r.confidence) if r.confidence is not None else None,
        "null_p": float(r.null_p_value) if r.null_p_value is not None else None,
        "null_type": str(r.null_type) if r.null_type is not None else None,
        "primary_metric": {k: float(v) if isinstance(v, (int, float, np.floating)) else v
                           for k, v in (r.primary_metric or {}).items()},
        "secondary_metrics": {k: float(v) if isinstance(v, (int, float, np.floating)) else v
                              for k, v in (r.secondary_metrics or {}).items()},
        "trajectory_metrics": {k: float(v) if isinstance(v, (int, float, np.floating)) else v
                               for k, v in (tm or {}).items()
                               if not isinstance(v, (list, np.ndarray, dict))},
        "exclusions_checked": list(r.exclusions_checked or []),
        "exclusion_results": dict(r.exclusion_results or {}),
        "effect_size": dict(r.effect_size or {}),
        "metadata_available": bool(r.metadata_available),
        "warnings": list(r.warnings or []),
    }


# --- Main -------------------------------------------------------------------

def main() -> None:
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    out: dict[str, Any] = {
        "schema_version": 1,
        "generated_by": "scripts/sprint24_baseline.py",
        "detector_constants": {
            "SCREENING_MORAN_SPEARMAN_MIN":
                P18ConsensusDetector.SCREENING_MORAN_SPEARMAN_MIN,
            "SCREENING_MORAN_FINAL_MIN":
                P18ConsensusDetector.SCREENING_MORAN_FINAL_MIN,
            "SCREENING_MORAN_GROWTH_MIN":
                P18ConsensusDetector.SCREENING_MORAN_GROWTH_MIN,
            "CONFIRMATION_WALL_SPEARMAN_MAX":
                P18ConsensusDetector.CONFIRMATION_WALL_SPEARMAN_MAX,
            "CONFIRMATION_WALL_FINAL_MAX":
                P18ConsensusDetector.CONFIRMATION_WALL_FINAL_MAX,
            "CONFIRMATION_WALL_DECAY_MIN":
                P18ConsensusDetector.CONFIRMATION_WALL_DECAY_MIN,
            "DEFINITIVE_MORAN_FINAL_MIN":
                P18ConsensusDetector.DEFINITIVE_MORAN_FINAL_MIN,
            "DEFINITIVE_MORAN_FINAL_MAX":
                P18ConsensusDetector.DEFINITIVE_MORAN_FINAL_MAX,
            "DEFINITIVE_WALL_FINAL_MIN":
                P18ConsensusDetector.DEFINITIVE_WALL_FINAL_MIN,
            "DEFINITIVE_MINORITY_FINAL_MIN":
                P18ConsensusDetector.DEFINITIVE_MINORITY_FINAL_MIN,
        },
        "voter_runs": [],
        "schelling_runs": [],
    }

    det = P18ConsensusDetector(n_permutations=199, seed=0)

    print("=== VOTER baseline (canonical positive) ===", flush=True)
    for L in VOTER_LS:
        n_steps = VOTER_NSTEPS[L]
        for seed in VOTER_SEEDS:
            t0 = time.time()
            history, meta = run_voter(L, seed, n_steps)
            t_run = time.time() - t0
            t1 = time.time()
            res = detect_one(history, meta, det)
            t_det = time.time() - t1
            res["L"] = L
            res["seed"] = seed
            res["n_steps"] = n_steps
            res["model"] = "voter"
            res["model_metadata_used"] = meta
            res["timing_s"] = {"run": round(t_run, 2),
                               "detect": round(t_det, 2)}
            out["voter_runs"].append(res)
            print(f"  voter L={L} seed={seed:3d}: tier={res['tier']:<14} "
                  f"p={res['null_p']:.4f}  "
                  f"moran_final={res['primary_metric'].get('moran_final_qtr_mean', float('nan')):.3f} "
                  f"wall_final={res['secondary_metrics'].get('wall_final_qtr_mean', float('nan')):.3f} "
                  f"P1_excl={res['exclusion_results'].get('P1', '-')}  "
                  f"({t_run:.0f}+{t_det:.0f}s)", flush=True)

    print("\n=== SCHELLING baseline (negative across parameter regimes) ===",
          flush=True)
    for threshold in SCHELLING_THRESHOLDS:
        meta = schelling_realistic_metadata(threshold)
        for seed in SCHELLING_SEEDS:
            t0 = time.time()
            history = run_schelling(grid_size=SCHELLING_L,
                                    density=SCHELLING_DENSITY,
                                    threshold=threshold,
                                    n_steps=SCHELLING_NSTEPS,
                                    seed=seed)
            t_run = time.time() - t0
            t1 = time.time()
            res_meta = detect_one(history, meta, det)
            res_none = detect_one(history, None, det)
            t_det = time.time() - t1
            for tag, res in [("with_metadata", res_meta),
                             ("none_metadata", res_none)]:
                res["L"] = SCHELLING_L
                res["seed"] = seed
                res["n_steps"] = SCHELLING_NSTEPS
                res["density"] = SCHELLING_DENSITY
                res["threshold"] = threshold
                res["model"] = "schelling"
                res["metadata_path"] = tag
                res["model_metadata_used"] = meta if tag == "with_metadata" else None
                res["timing_s"] = {"run": round(t_run, 2),
                                   "detect": round(t_det / 2, 2)}
                out["schelling_runs"].append(res)
            res = res_meta
            print(f"  schelling thr={threshold:.3f} seed={seed}: "
                  f"tier={res['tier']:<14} p={res['null_p']:.4f}  "
                  f"moran_final={res['primary_metric'].get('moran_final_qtr_mean', float('nan')):.3f} "
                  f"wall_final={res['secondary_metrics'].get('wall_final_qtr_mean', float('nan')):.3f} "
                  f"P1_excl={res['exclusion_results'].get('P1', '-')}  "
                  f"({t_run:.0f}+{t_det:.0f}s)", flush=True)

    out_path = OUTPUT_DIR / "baseline_voter_schelling.json"
    with open(out_path, "w") as f:
        json.dump(out, f, indent=2, default=str)
    print(f"\nSaved {out_path} ({out_path.stat().st_size} bytes)")
    print(f"Voter runs: {len(out['voter_runs'])}  "
          f"Schelling runs (incl. both metadata paths): "
          f"{len(out['schelling_runs'])}")


if __name__ == "__main__":
    main()
