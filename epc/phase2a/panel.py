"""Phase-2a panel runner.

The harness runs a detector against:
- 5 seeds of its canonical positive (positive distribution),
- 10 synthetic substrates (Class A),
- 10 catalog-derived non-positives (Class B; 1 self-replacement → P12 fallback),
- 10 pattern-specific failed regimes (Class C),

scores each, computes per-class TNR, overall TNR, and Cohen's d between
positive and pooled-negative score distributions, and writes the result to
a JSON file matching the schema in ``docs/phase2a_panel_spec.md``.

A detector is a callable ``(history, metadata) -> result`` where ``result``
exposes ``.detected`` (bool) and ``.confidence`` (float).
"""

from __future__ import annotations

import json
import os
import subprocess
import time
from dataclasses import dataclass
from typing import Any, Callable, Dict, List, Optional

import numpy as np

from epc.phase2a import PANEL_VERSION
from epc.phase2a import catalog as catalog_mod
from epc.phase2a import synthetic as synth_mod


# --- Score / verdict helpers ------------------------------------------------

def _score(result: Any) -> float:
    """Confidence is the universal score across detector result types."""
    return float(getattr(result, "confidence", 0.0))


def _detected(result: Any) -> bool:
    return bool(getattr(result, "detected", False))


def _verdict(result: Any) -> str:
    """Tier name as a stable string."""
    tier = getattr(result, "tier", None)
    if tier is None:
        return "none"
    name = getattr(tier, "name", None)
    if name is not None:
        return name.lower()
    return str(tier).lower()


# --- Stats ------------------------------------------------------------------

def compute_tnr(verdicts: List[bool]) -> float:
    """TNR = fraction of negatives correctly rejected (detected==False).

    ``verdicts`` is a list of ``detected`` booleans on negative substrates.
    """
    if not verdicts:
        return float("nan")
    n_correct = sum(1 for v in verdicts if v is False)
    return n_correct / len(verdicts)


def cohens_d(positive_scores: List[float], negative_scores: List[float]) -> float:
    """Cohen's d between two independent score distributions.

    d = (mean_pos − mean_neg) / pooled_std
    pooled_std = sqrt(((n_pos−1)·var_pos + (n_neg−1)·var_neg) / (n_pos + n_neg − 2))
    """
    pos = np.asarray(positive_scores, dtype=float)
    neg = np.asarray(negative_scores, dtype=float)
    if pos.size < 2 or neg.size < 2:
        return float("nan")
    var_pos = float(np.var(pos, ddof=1))
    var_neg = float(np.var(neg, ddof=1))
    pooled_var = ((pos.size - 1) * var_pos + (neg.size - 1) * var_neg) / (pos.size + neg.size - 2)
    if pooled_var <= 0.0:
        return float("inf") if pos.mean() > neg.mean() else float("-inf")
    return (float(pos.mean()) - float(neg.mean())) / float(np.sqrt(pooled_var))


def overall_verdict(
    overall_tnr: float,
    per_class_tnr: Dict[str, float],
    overall_tnr_threshold: float = 0.95,
    weak_class_threshold: float = 0.90,
) -> str:
    """Apply the spec's PASS criteria.

    PASS               : overall_tnr ≥ 0.95 and every reported class ≥ 0.90.
    PASS-with-weakness : overall_tnr ≥ 0.95 but some class < 0.90.
    PARTIAL            : overall_tnr < 0.95.
    """
    if overall_tnr >= overall_tnr_threshold:
        weak = [c for c, t in per_class_tnr.items() if t < weak_class_threshold]
        if weak:
            return "PASS-with-weakness"
        return "PASS"
    return "PARTIAL"


# --- Runner -----------------------------------------------------------------

@dataclass
class SubstrateResult:
    substrate: str
    score: float
    verdict: str
    detected: bool


def _run_one(
    detector_fn: Callable[[List[Dict[str, Any]], Optional[Dict[str, Any]]], Any],
    history: List[Dict[str, Any]],
    metadata: Optional[Dict[str, Any]] = None,
) -> Any:
    return detector_fn(history, metadata)


def _git_head() -> str:
    try:
        return subprocess.check_output(
            ["git", "rev-parse", "HEAD"], stderr=subprocess.DEVNULL
        ).decode().strip()
    except Exception:
        return "unknown"


def run_panel(
    detector_fn: Callable[[List[Dict[str, Any]], Optional[Dict[str, Any]]], Any],
    *,
    pattern_id: str,
    detector_format: str,
    canonical_positive_runs: List[Dict[str, Any]],
    canonical_metadata: Optional[Dict[str, Any]] = None,
    failed_regime_module: Any,
    output_path: Optional[str] = None,
    target_steps: int = 200,
    target_shape: tuple = (32, 32),
    target_n: int = 300,
    seed: int = 42,
    verbose: bool = False,
) -> Dict[str, Any]:
    """Run the Phase-2a panel against ``detector_fn`` and return summary dict.

    ``canonical_positive_runs`` is a list of pre-built canonical-positive
    histories (one per seed). Recommended size: 5.
    """
    t0 = time.time()
    rng = np.random.default_rng(seed)

    # 1) Canonical positive distribution.
    positive_results: List[Dict[str, Any]] = []
    positive_scores: List[float] = []
    canonical_first = canonical_positive_runs[0]
    for i, history in enumerate(canonical_positive_runs):
        result = _run_one(detector_fn, history, canonical_metadata)
        positive_scores.append(_score(result))
        positive_results.append({
            "seed_index": i,
            "score": _score(result),
            "verdict": _verdict(result),
            "detected": _detected(result),
        })
        if verbose:
            print(f"  [pos {i}] verdict={_verdict(result)} score={_score(result):.3f}")

    # 2) Class A — synthetic.
    synthetic_results: List[Dict[str, Any]] = []
    class_a_kwargs = {"n_steps": target_steps}
    if detector_format == "grid":
        class_a_kwargs["shape"] = target_shape
    elif detector_format == "phases":
        class_a_kwargs["n"] = target_n

    for i, (sub_id, gen) in enumerate(synth_mod.SYNTHETIC_GENERATORS.items()):
        s = int(rng.integers(0, 2**31 - 1))
        if sub_id in ("permutation_shuffled", "time_shuffled"):
            history = gen(detector_format, s, positive=canonical_first)
        else:
            history = gen(detector_format, s, **class_a_kwargs)
        result = _run_one(detector_fn, history, canonical_metadata)
        synthetic_results.append({
            "substrate": sub_id,
            "seed": s,
            "score": _score(result),
            "verdict": _verdict(result),
            "detected": _detected(result),
        })
        if verbose:
            print(f"  [syn {i:2d} {sub_id:24s}] verdict={_verdict(result)} score={_score(result):.3f}")

    # 3) Class B — catalog-derived (with self-replacement → fallback).
    catalog_results: List[Dict[str, Any]] = []
    catalog_ids = catalog_mod.catalog_ids_for_pattern(pattern_id)
    for i, sub_id in enumerate(catalog_ids):
        history = catalog_mod.load_catalog_substrate_for_format(
            sub_id, detector_format,
            target_steps=target_steps, target_shape=target_shape, target_n=target_n,
        )
        result = _run_one(detector_fn, history, canonical_metadata)
        catalog_results.append({
            "substrate": sub_id,
            "score": _score(result),
            "verdict": _verdict(result),
            "detected": _detected(result),
        })
        if verbose:
            print(f"  [cat {i:2d} {sub_id:24s}] verdict={_verdict(result)} score={_score(result):.3f}")

    # 4) Class C — failed regimes.
    failed_results: List[Dict[str, Any]] = []
    config = failed_regime_module.CONFIG
    for i, regime in enumerate(config["regimes"]):
        history = failed_regime_module.build_substrate(regime)
        result = _run_one(detector_fn, history, canonical_metadata)
        failed_results.append({
            "substrate": regime["label"],
            "params": regime["params"],
            "seed": regime["seed"],
            "score": _score(result),
            "verdict": _verdict(result),
            "detected": _detected(result),
        })
        if verbose:
            print(f"  [fai {i:2d} {regime['label']:24s}] verdict={_verdict(result)} score={_score(result):.3f}")

    # 5) Aggregates.
    syn_detected = [r["detected"] for r in synthetic_results]
    cat_detected = [r["detected"] for r in catalog_results]
    fai_detected = [r["detected"] for r in failed_results]
    all_neg_detected = syn_detected + cat_detected + fai_detected
    all_neg_scores = [r["score"] for r in synthetic_results + catalog_results + failed_results]

    syn_tnr = compute_tnr(syn_detected)
    cat_tnr = compute_tnr(cat_detected)
    fai_tnr = compute_tnr(fai_detected)
    overall_tnr = compute_tnr(all_neg_detected)
    d = cohens_d(positive_scores, all_neg_scores)

    per_class = {
        "synthetic_tnr": syn_tnr,
        "catalog_tnr": cat_tnr,
        "failed_regime_tnr": fai_tnr,
    }
    verdict = overall_verdict(overall_tnr, per_class)

    summary = {
        "pattern_id": pattern_id,
        "head_commit": _git_head(),
        "panel_version": PANEL_VERSION,
        "detector_format": detector_format,
        "elapsed_seconds": round(time.time() - t0, 2),
        "canonical_positive": {
            "n_seeds": len(positive_scores),
            "scores": positive_scores,
            "mean_score": float(np.mean(positive_scores)) if positive_scores else float("nan"),
            "results": positive_results,
        },
        "synthetic": synthetic_results,
        "catalog": catalog_results,
        "failed_regime": failed_results,
        "summary": {
            "n_negatives": len(all_neg_detected),
            "overall_tnr": overall_tnr,
            "synthetic_tnr": syn_tnr,
            "catalog_tnr": cat_tnr,
            "failed_regime_tnr": fai_tnr,
            "cohens_d_positive_vs_panel": d,
            "verdict": verdict,
        },
    }

    if output_path is not None:
        os.makedirs(os.path.dirname(output_path), exist_ok=True)
        with open(output_path, "w") as f:
            json.dump(summary, f, indent=2, default=_json_safe)

    return summary


def _json_safe(o: Any) -> Any:
    if isinstance(o, (np.integer,)):
        return int(o)
    if isinstance(o, (np.floating,)):
        return float(o)
    if isinstance(o, np.ndarray):
        return o.tolist()
    raise TypeError(f"Object of type {type(o).__name__} is not JSON serializable")
