"""Run the Phase-2a panel against P18 and/or P9 and write JSON outputs.

Usage::

    PYTHONPATH=. python3.12 analysis/run_phase2a_panel.py p18
    PYTHONPATH=. python3.12 analysis/run_phase2a_panel.py p9
    PYTHONPATH=. python3.12 analysis/run_phase2a_panel.py both

Outputs:
    analysis/outputs/p18_phase2a_panel.json
    analysis/outputs/p9_phase2a_panel.json
"""

from __future__ import annotations

import sys
from typing import Any, Dict, List, Optional

import numpy as np

from epc.phase2a.panel import run_panel
from epc.phase2a.failed_regimes import p18_voter as p18_failed
from epc.phase2a.failed_regimes import p9_kuramoto as p9_failed


# --- Canonical positives -----------------------------------------------------

def build_p18_positives(n_seeds: int = 5) -> tuple[List[List[Dict[str, Any]]], Dict[str, Any]]:
    from epc.models.voter import VoterModel
    runs: List[List[Dict[str, Any]]] = []
    metadata: Dict[str, Any] = {}
    for seed in range(n_seeds):
        m = VoterModel(rows=64, cols=64, seed=seed)
        runs.append(m.run(n_steps=400))
        if seed == 0:
            metadata = m.get_metadata()
    return runs, metadata


def build_p9_positives(n_seeds: int = 5) -> tuple[List[List[Dict[str, Any]]], Dict[str, Any]]:
    from epc.models.kuramoto import KuramotoModel, KuramotoParams
    runs: List[List[Dict[str, Any]]] = []
    metadata: Dict[str, Any] = {}
    for seed in range(n_seeds):
        # n_steps=6000 → 600 records → n_T_osc ≈ 12 ≥ 10 (P9 screening prerequisite).
        m = KuramotoModel(KuramotoParams(N=300, K=8.0, gamma=0.5, dt=0.05, seed=seed))
        runs.append(m.run(n_steps=6000, record_every=10))
        if seed == 0:
            metadata = m.get_metadata()
    return runs, metadata


# --- Detector adapters (history, metadata) → result --------------------------

def make_p18_detector_fn(n_permutations: int = 99, seed: int = 42):
    from epc.detectors.p18_consensus import P18ConsensusDetector
    detector = P18ConsensusDetector(n_permutations=n_permutations, seed=seed)
    def fn(history, metadata=None):
        return detector.detect(history, model_metadata=metadata)
    return fn


def make_p9_detector_fn(n_null_runs: int = 99, seed: int = 42):
    from epc.detectors.p9_synchronization import detect_p9
    def fn(history, metadata=None):
        return detect_p9(history, n_null_runs=n_null_runs, seed=seed, model_metadata=metadata)
    # The P9 result tier is a string ("definitive" / "screening" / etc.) rather
    # than an enum; the panel runner already handles both via `_verdict()`.
    return fn


# --- Entry points ------------------------------------------------------------

def run_p18(out_path: str = "analysis/outputs/p18_phase2a_panel.json", verbose: bool = True) -> Dict[str, Any]:
    print(f"--- Running P18 panel → {out_path}")
    positives, metadata = build_p18_positives(n_seeds=5)
    detector_fn = make_p18_detector_fn(n_permutations=99, seed=42)
    return run_panel(
        detector_fn,
        pattern_id="P18",
        detector_format="grid",
        canonical_positive_runs=positives,
        canonical_metadata=metadata,
        failed_regime_module=p18_failed,
        output_path=out_path,
        target_steps=200,
        target_shape=(32, 32),
        verbose=verbose,
    )


def run_p9(out_path: str = "analysis/outputs/p9_phase2a_panel.json", verbose: bool = True) -> Dict[str, Any]:
    print(f"--- Running P9 panel → {out_path}")
    positives, metadata = build_p9_positives(n_seeds=5)
    detector_fn = make_p9_detector_fn(n_null_runs=99, seed=42)
    return run_panel(
        detector_fn,
        pattern_id="P9",
        detector_format="phases",
        canonical_positive_runs=positives,
        canonical_metadata=metadata,
        failed_regime_module=p9_failed,
        output_path=out_path,
        target_steps=600,  # cadence=10 → n_T_osc ≈ 12 (≥10 screening prereq)
        target_n=300,
        verbose=verbose,
    )


def main(argv: Optional[List[str]] = None) -> int:
    argv = argv if argv is not None else sys.argv[1:]
    which = argv[0] if argv else "both"

    summaries: Dict[str, Dict[str, Any]] = {}
    if which in ("p18", "both"):
        summaries["P18"] = run_p18()
    if which in ("p9", "both"):
        summaries["P9"] = run_p9()

    def _fmt(x):
        return "  N/A " if x is None else f"{x:>5.3f}"

    print()
    print("=" * 80)
    print(f"{'pattern':<8} {'overall':>7} {'syn':>6} {'cat':>6} {'fai':>6} {'d':>6} {'verdict':<22}")
    print("=" * 80)
    for pid, summary in summaries.items():
        s = summary["summary"]
        print(
            f"{pid:<8} {_fmt(s['overall_tnr']):>7} {_fmt(s['synthetic_tnr']):>6} "
            f"{_fmt(s['catalog_tnr']):>6} {_fmt(s['failed_regime_tnr']):>6} "
            f"{_fmt(s['cohens_d_positive_vs_panel']):>6} {s['verdict']:<22}"
        )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
