"""Phase-2a panel runner (v1.1).

The harness runs a detector against:
- N seeds of its canonical positive (positive distribution; recommended 5),
- 10 synthetic substrates (Class A, fixed),
- substrate-typed catalog-derived non-positives + Class B' synthetic
  supplements (Class B; variable size driven by
  :func:`epc.phase2a.catalog.class_b_for_pattern`),
- pattern-specific failed regimes if Class C is populated for the pattern,
  OR skipped entirely if Class C is declared N/A in the failed-regime config.

It scores each substrate, computes per-class TNR, overall TNR, and Cohen's d
between the positive and pooled-negative score distributions, and writes the
result to a JSON file matching the schema in
``docs/phase2a_panel_spec.md`` (v1.1).

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
from epc.phase2a import structured as structured_mod
from epc.phase2a import synthetic as synth_mod
from epc.phase2a.detector_invariance import get_flags as _get_invariance_flags


# Class A substrate ids whose semantics are tied to a primary-metric invariance.
# v1.2 spec §"Change 1": when the detector under test has the corresponding
# invariance flag set to True, the substrate is degenerate-by-construction
# (it preserves the canonical positive's aggregate distribution) and the
# harness skips it.
PERM_SHUFFLED_SUBSTRATE = "permutation_shuffled"
TIME_SHUFFLED_SUBSTRATE = "time_shuffled"
SKIP_VERDICT = "SKIPPED-degenerate-by-construction"


# Class size below which per-class TNR is reported as advisory only and does
# not gate the PASS verdict (v1.1 spec §"PASS criterion").
ADVISORY_CLASS_SIZE = 5


def _adapt_supplement_history_to_opinions(
    history: List[Dict[str, Any]],
    target_n: int = 100,
) -> List[Dict[str, Any]]:
    """Adapt a supplement's grid-format history to opinions format.

    Network supplements (random_graph_evolution, network_random_walks) return
    grid-format histories.  P21's detector requires ``opinions`` key.  This
    adapter flattens each grid frame and applies a rank transform to [0, 1] so
    the resulting distribution is uniform (unimodal) — P21 should reject.
    """
    import numpy as _np

    def _rank_uniform(arr: _np.ndarray) -> _np.ndarray:
        n = len(arr)
        if n <= 1:
            return arr.astype(float)
        order = _np.argsort(arr)
        result = _np.empty(n, dtype=float)
        result[order] = _np.arange(n, dtype=float) / (n - 1)
        return result

    adapted = []
    for snap in history:
        if "opinions" in snap:
            adapted.append(snap)
            continue
        if "grid" in snap:
            flat = _np.asarray(snap["grid"]).flatten().astype(float)
        else:
            # Fallback: uniform null.
            flat = _np.linspace(0.0, 1.0, max(target_n, 4))
        if len(flat) >= target_n:
            flat = flat[:target_n]
        else:
            reps = target_n // len(flat) + 1
            flat = _np.tile(flat, reps)[:target_n]
        opinions = _rank_uniform(flat)
        adapted.append({"opinions": opinions, "step": snap.get("step", 0)})
    return adapted


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
    """TNR = fraction of negatives correctly rejected (detected==False)."""
    if not verdicts:
        return float("nan")
    n_correct = sum(1 for v in verdicts if v is False)
    return n_correct / len(verdicts)


def cohens_d(positive_scores: List[float], negative_scores: List[float]) -> float:
    """Cohen's d between two independent score distributions.

    d = (mean_pos − mean_neg) / pooled_std with sample-variance pooling.

    When both score arrays are effectively constant (pooled variance below
    the degeneracy threshold ``1e-12``) the comparison is undefined in the
    standard formula. We disambiguate from the difference of means:
    +inf / -inf for non-zero mean separation, 0 when means are also equal.
    The degeneracy threshold absorbs ULP-level numerical noise that
    accumulates in ``np.var`` on long constant arrays — without it,
    constant-vs-constant inputs spuriously produce non-trivial finite d.
    """
    pos = np.asarray(positive_scores, dtype=float)
    neg = np.asarray(negative_scores, dtype=float)
    if pos.size < 2 or neg.size < 2:
        return float("nan")
    var_pos = float(np.var(pos, ddof=1))
    var_neg = float(np.var(neg, ddof=1))
    pooled_var = ((pos.size - 1) * var_pos + (neg.size - 1) * var_neg) / (pos.size + neg.size - 2)
    mean_diff = float(pos.mean()) - float(neg.mean())
    if pooled_var <= 1e-12:
        if abs(mean_diff) <= 1e-9:
            return 0.0
        return float("inf") if mean_diff > 0 else float("-inf")
    return mean_diff / float(np.sqrt(pooled_var))


def overall_verdict(
    overall_tnr: float,
    per_class_tnr_gating: Dict[str, float],
    cohens_d_value: Optional[float] = None,
    *,
    overall_tnr_threshold: float = 0.95,
    weak_class_threshold: float = 0.90,
    cohens_d_pass_threshold: float = 1.0,
    cohens_d_partial_threshold: float = 0.5,
) -> str:
    """v1.1 verdict labels (spec §"PASS criterion").

    - **PASS** — overall_tnr ≥ 0.95 AND every gating class ≥ 0.90 AND d ≥ 1.0.
    - **PASS-with-weakness** — overall_tnr ≥ 0.95 AND d ≥ 1.0 but a gating class < 0.90.
    - **PARTIAL** — overall_tnr < 0.95 (or d < 1.0 with TNR ≥ 0.95) AND d ≥ 0.5.
    - **FAIL** — overall_tnr < 0.95 AND d < 0.5.

    ``per_class_tnr_gating`` only contains classes large enough to gate
    (≥ ADVISORY_CLASS_SIZE substrates); advisory-only classes are excluded
    by the caller.
    """
    # None → nan so it fails every threshold; ±inf are kept as-is so a
    # truly degenerate-perfect or degenerate-broken detector still gates correctly.
    if cohens_d_value is None:
        cohens_d_value = float("nan")

    if overall_tnr >= overall_tnr_threshold:
        if cohens_d_value >= cohens_d_pass_threshold:
            weak = [c for c, t in per_class_tnr_gating.items() if t < weak_class_threshold]
            return "PASS-with-weakness" if weak else "PASS"
        # TNR passes but effect size is below the PASS gate — treat as PARTIAL
        # if there is at least signal (d ≥ partial threshold), else FAIL.
        if cohens_d_value >= cohens_d_partial_threshold:
            return "PARTIAL"
        return "FAIL"

    # overall_tnr below the PASS gate — verdict split by Cohen's d.
    if cohens_d_value >= cohens_d_partial_threshold:
        return "PARTIAL"
    return "FAIL"


def _per_class(
    detected_list: List[bool],
    *,
    advisory_size: int = ADVISORY_CLASS_SIZE,
) -> Dict[str, Any]:
    """Compute per-class TNR + advisory status."""
    n = len(detected_list)
    if n == 0:
        return {"n": 0, "tnr": None, "advisory": True}
    return {
        "n": n,
        "tnr": compute_tnr(detected_list),
        "advisory": n < advisory_size,
    }


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
    canonical_positive_runs: List[List[Dict[str, Any]]],
    canonical_metadata: Optional[Dict[str, Any]] = None,
    failed_regime_module: Any,
    output_path: Optional[str] = None,
    target_steps: int = 200,
    target_shape: tuple = (32, 32),
    target_n: int = 300,
    seed: int = 42,
    verbose: bool = False,
) -> Dict[str, Any]:
    """Run the v1.1 Phase-2a panel against ``detector_fn`` and return summary."""
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

    # 2) Class A — synthetic (with v1.2 invariance-driven skipping).
    invariance = _get_invariance_flags(pattern_id)
    synthetic_results: List[Dict[str, Any]] = []
    class_a_kwargs: Dict[str, Any] = {"n_steps": target_steps}
    if detector_format == "grid":
        class_a_kwargs["shape"] = target_shape
    elif detector_format == "phases":
        class_a_kwargs["n"] = target_n
    elif detector_format == "sequence":
        # n = road/array length; target_n is reused as the sequence dimension.
        class_a_kwargs["n"] = target_n
    elif detector_format == "avalanches":
        # Avalanche-format generators only take n_avalanches; n_steps unused.
        class_a_kwargs = {}
    elif detector_format == "opinions":
        # Opinions generators use their own OPINIONS_DEFAULT_N; n_steps still applies.
        pass  # n_steps already in class_a_kwargs
    elif detector_format == "scalar_timeseries":
        # Scalar-timeseries generators use their own defaults; n_steps applies.
        pass  # n_steps already in class_a_kwargs
    elif detector_format == "choice_timeseries":
        # Choice-timeseries generators use their own defaults; n_steps applies.
        pass  # n_steps already in class_a_kwargs
    elif detector_format == "state_vector":
        # State-vector generators use their own defaults (N, P, trials).
        pass  # n_steps already in class_a_kwargs; generators ignore it
    elif detector_format == "noise_sweep":
        # Noise-sweep generators use their own defaults; n_steps unused.
        class_a_kwargs = {}
    elif detector_format == "density_sweep":
        # Density-sweep generators use their own defaults (n_levels, density range).
        pass  # n_steps already in class_a_kwargs
    elif detector_format == "canalization_bundle":
        # Canalization-bundle generators use their own defaults (n_ics, n_dims).
        pass  # n_steps already in class_a_kwargs

    class_a_total = len(synth_mod.SYNTHETIC_GENERATORS)
    for i, (sub_id, gen) in enumerate(synth_mod.SYNTHETIC_GENERATORS.items()):
        # v1.2: skip degenerate-by-construction substrates for invariant detectors.
        skip_reason: Optional[str] = None
        if sub_id == PERM_SHUFFLED_SUBSTRATE and invariance.permutation_invariant:
            skip_reason = "primary_metric_permutation_invariant"
        elif sub_id == TIME_SHUFFLED_SUBSTRATE and invariance.time_shuffle_invariant:
            skip_reason = "primary_metric_time_shuffle_invariant"

        if skip_reason is not None:
            synthetic_results.append({
                "substrate": sub_id,
                "verdict": SKIP_VERDICT,
                "skip_reason": skip_reason,
            })
            if verbose:
                print(f"  [syn {i:2d} {sub_id:24s}] SKIPPED ({skip_reason})")
            continue

        s = int(rng.integers(0, 2**31 - 1))
        if sub_id in (PERM_SHUFFLED_SUBSTRATE, TIME_SHUFFLED_SUBSTRATE):
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

    class_a_evaluated = sum(1 for r in synthetic_results if r["verdict"] != SKIP_VERDICT)

    # 3) Class B (v1.1) — substrate-typed catalog mates + Class B' supplements.
    class_b = catalog_mod.class_b_for_pattern(pattern_id)
    catalog_mates: List[str] = class_b["catalog_mates"]
    synthetic_supplements: List[str] = class_b["synthetic_supplements"]
    class_b_substrate_type: Optional[str] = class_b["substrate_type"]

    catalog_results: List[Dict[str, Any]] = []
    for sub_id in catalog_mates:
        history = catalog_mod.load_catalog_substrate_for_format(
            sub_id, detector_format,
            target_steps=target_steps, target_shape=target_shape, target_n=target_n,
        )
        result = _run_one(detector_fn, history, canonical_metadata)
        catalog_results.append({
            "substrate": sub_id,
            "kind": "catalog_mate",
            "score": _score(result),
            "verdict": _verdict(result),
            "detected": _detected(result),
        })
        if verbose:
            print(f"  [cat {sub_id:30s}] verdict={_verdict(result)} score={_score(result):.3f}")
    for supp_id in synthetic_supplements:
        builder = structured_mod.SUPPLEMENT_BUILDERS[supp_id]
        s = int(rng.integers(0, 2**31 - 1))
        history = builder(s)
        # Supplement builders return their substrate type's native format (e.g.
        # "network" supplements return grid-format histories).  For opinions-
        # format detectors adapt via rank transform → uniform → unimodal → TN.
        if detector_format == "opinions" and history and "opinions" not in history[0]:
            history = _adapt_supplement_history_to_opinions(history, target_n=target_n)
        result = _run_one(detector_fn, history, canonical_metadata)
        catalog_results.append({
            "substrate": supp_id,
            "kind": "synthetic_supplement",
            "seed": s,
            "score": _score(result),
            "verdict": _verdict(result),
            "detected": _detected(result),
        })
        if verbose:
            print(f"  [sup {supp_id:30s}] verdict={_verdict(result)} score={_score(result):.3f}")

    # 4) Class C — failed regimes (with N/A escape hatch).
    fr_config = getattr(failed_regime_module, "CONFIG", {})
    class_c_status = "populated"
    class_c_n_a_reason: Optional[str] = None
    failed_results: List[Dict[str, Any]] = []

    if fr_config.get("status") == "N/A":
        class_c_status = "N/A"
        class_c_n_a_reason = fr_config.get("n_a_reason", "")
        if verbose:
            print(f"  [class C N/A] {class_c_n_a_reason}")
    else:
        for regime in fr_config.get("regimes", []):
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
                print(f"  [fai {regime['label']:30s}] verdict={_verdict(result)} score={_score(result):.3f}")

    # 5) Aggregates. Skipped Class A substrates are excluded from TNR + d.
    syn_detected = [r["detected"] for r in synthetic_results if r["verdict"] != SKIP_VERDICT]
    cat_detected = [r["detected"] for r in catalog_results]
    fai_detected = [r["detected"] for r in failed_results]
    syn_scores = [r["score"] for r in synthetic_results if r["verdict"] != SKIP_VERDICT]
    cat_scores = [r["score"] for r in catalog_results]
    fai_scores = [r["score"] for r in failed_results]

    if class_c_status == "N/A":
        all_neg_detected = syn_detected + cat_detected
        all_neg_scores = syn_scores + cat_scores
    else:
        all_neg_detected = syn_detected + cat_detected + fai_detected
        all_neg_scores = syn_scores + cat_scores + fai_scores

    syn_class = _per_class(syn_detected)
    cat_class = _per_class(cat_detected)
    fai_class = _per_class(fai_detected) if class_c_status == "populated" else {"n": 0, "tnr": None, "advisory": True}

    overall_tnr = compute_tnr(all_neg_detected)
    d = cohens_d(positive_scores, all_neg_scores)

    # Build the gating class dict (only ≥ ADVISORY_CLASS_SIZE classes count).
    per_class_gating: Dict[str, float] = {}
    if syn_class["tnr"] is not None and not syn_class["advisory"]:
        per_class_gating["synthetic_tnr"] = syn_class["tnr"]
    if cat_class["tnr"] is not None and not cat_class["advisory"]:
        per_class_gating["catalog_tnr"] = cat_class["tnr"]
    if fai_class["tnr"] is not None and not fai_class["advisory"]:
        per_class_gating["failed_regime_tnr"] = fai_class["tnr"]

    verdict = overall_verdict(overall_tnr, per_class_gating, d)

    summary = {
        "pattern_id": pattern_id,
        "head_commit": _git_head(),
        "panel_version": PANEL_VERSION,
        "detector_format": detector_format,
        "elapsed_seconds": round(time.time() - t0, 2),
        "detector_invariance": {
            "permutation_invariant": invariance.permutation_invariant,
            "time_shuffle_invariant": invariance.time_shuffle_invariant,
            "primary_metric": invariance.primary_metric,
        },
        "canonical_positive": {
            "n_seeds": len(positive_scores),
            "scores": positive_scores,
            "mean_score": float(np.mean(positive_scores)) if positive_scores else float("nan"),
            "results": positive_results,
        },
        "synthetic": synthetic_results,
        "catalog": catalog_results,
        "failed_regime": failed_results,
        "class_b_composition": {
            "substrate_type": class_b_substrate_type,
            "catalog_mates": catalog_mates,
            "synthetic_supplements": synthetic_supplements,
        },
        "class_c_status": class_c_status,
        "class_c_n_a_reason": class_c_n_a_reason,
        "summary": {
            "n_negatives": len(all_neg_detected),
            "overall_tnr": overall_tnr,
            "class_a_size_total": class_a_total,
            "class_a_size_evaluated": class_a_evaluated,
            "synthetic": syn_class,
            "catalog": cat_class,
            "failed_regime": fai_class,
            "synthetic_tnr": syn_class["tnr"],
            "catalog_tnr": cat_class["tnr"],
            "catalog_tnr_advisory": cat_class["advisory"],
            "failed_regime_tnr": fai_class["tnr"],
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
