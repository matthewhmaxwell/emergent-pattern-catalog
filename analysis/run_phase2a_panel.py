"""Run the Phase-2a panel against P18, P9, P22, P27, etc. and write JSON outputs.

Usage::

    PYTHONPATH=. python3.12 analysis/run_phase2a_panel.py p18
    PYTHONPATH=. python3.12 analysis/run_phase2a_panel.py p9
    PYTHONPATH=. python3.12 analysis/run_phase2a_panel.py both
    PYTHONPATH=. python3.12 analysis/run_phase2a_panel.py p22
    PYTHONPATH=. python3.12 analysis/run_phase2a_panel.py p27
    PYTHONPATH=. python3.12 analysis/run_phase2a_panel.py p1

Outputs:
    analysis/outputs/p18_phase2a_panel.json
    analysis/outputs/p9_phase2a_panel.json
    analysis/outputs/p22_phase2a_panel.json
    analysis/outputs/p27_phase2a_panel.json
    analysis/outputs/p1_phase2a_panel.json
"""

from __future__ import annotations

import sys
from typing import Any, Dict, List, Optional

import numpy as np

from epc.phase2a.panel import run_panel
from epc.phase2a.failed_regimes import p18_voter as p18_failed
from epc.phase2a.failed_regimes import p9_kuramoto as p9_failed
from epc.phase2a.failed_regimes import p15_gol as p15_failed
from epc.phase2a.failed_regimes import p14_btw as p14_failed
from epc.phase2a.failed_regimes import p27_nowak_may as p27_failed
from epc.phase2a.failed_regimes import p22_sir as p22_failed
from epc.phase2a.failed_regimes import p1_schelling as p1_failed


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


def build_p15_positives(n_seeds: int = 5) -> tuple[List[List[Dict[str, Any]]], Dict[str, Any]]:
    """P15 canonical positive: dense random GoL (density=0.37, L=40, n_steps=300).

    Per ``tests/test_p15_generalized.py::test_gol_dense_definitive`` this is the
    canonical positive that reaches DEFINITIVE under P15's multi-variation
    reproducibility test. R-pentomino is too sparse for P15's structural-
    diversity screening.
    """
    from epc.models.game_of_life import GameOfLife
    runs: List[List[Dict[str, Any]]] = []
    metadata: Dict[str, Any] = {}
    for seed in range(n_seeds):
        m = GameOfLife(
            rows=40, cols=40, init_mode="random", init_density=0.37,
            boundary="periodic", seed=42 + seed,
        )
        runs.append(m.run(n_steps=300, record_every=1))
        if seed == 0:
            metadata = m.get_metadata()
    return runs, metadata


def build_p14_positives(n_seeds: int = 5) -> tuple[List[List[Dict[str, Any]]], Dict[str, Any]]:
    """P14 canonical positive: BTW sandpile at L=32, n_drive=10000.

    Each seed produces an avalanche-size array; we wrap as a single-element
    "history" so the panel runner's iteration model still works.
    """
    from epc.models.btw_sandpile import run_sandpile, BTWSandpileParams
    runs: List[List[Dict[str, Any]]] = []
    for seed in range(n_seeds):
        params = BTWSandpileParams(L=32, n_drive=10_000, n_burn=1_000, seed=seed)
        result = run_sandpile(params)
        runs.append([{
            "avalanche_sizes": result.avalanche_sizes,
            "avalanche_durations": result.avalanche_durations,
            "activity": result.activity,
            "energy": result.energy_history,
            "step": 0,
        }])
    metadata = {"is_self_tuned": True, "model": "btw_sandpile"}
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


def make_p15_detector_fn(seed: int = 42, n_variations: int = 8):
    """P15 detector with the canonical GoL step_fn so positives reach DEFINITIVE
    via the multi-variation reproducibility test (Sprint 8 generalization)."""
    from epc.detectors.p15_persistent_computation import (
        P15PersistentComputationDetector, make_step_fn_for_gol,
    )
    step_fn = make_step_fn_for_gol()
    detector = P15PersistentComputationDetector(
        step_fn=step_fn, n_variations=n_variations, seed=seed,
    )
    def fn(history, metadata=None):
        return detector.detect(history, model_metadata=metadata)
    return fn


def make_p14_detector_fn():
    """P14 detector wrapper: unwrap the single-element avalanches "history"
    and pass arrays directly to detect_p14."""
    from epc.detectors.p14_soc import detect_p14
    def fn(history, metadata=None):
        h0 = history[0] if isinstance(history, list) else history
        return detect_p14(
            avalanche_sizes=h0["avalanche_sizes"],
            avalanche_durations=h0.get("avalanche_durations"),
            activity=h0.get("activity"),
            energy=h0.get("energy"),
            null_sizes=h0.get("null_sizes"),
            is_self_tuned=metadata.get("is_self_tuned", True) if metadata else True,
        )
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


def run_p15(out_path: str = "analysis/outputs/p15_phase2a_panel.json", verbose: bool = True) -> Dict[str, Any]:
    print(f"--- Running P15 panel → {out_path}")
    positives, metadata = build_p15_positives(n_seeds=5)
    detector_fn = make_p15_detector_fn(seed=42)
    return run_panel(
        detector_fn,
        pattern_id="P15",
        detector_format="grid",
        canonical_positive_runs=positives,
        canonical_metadata=metadata,
        failed_regime_module=p15_failed,
        output_path=out_path,
        target_steps=200,
        target_shape=(32, 32),
        verbose=verbose,
    )


def run_p14(out_path: str = "analysis/outputs/p14_phase2a_panel.json", verbose: bool = True) -> Dict[str, Any]:
    print(f"--- Running P14 panel → {out_path}")
    positives, metadata = build_p14_positives(n_seeds=5)
    detector_fn = make_p14_detector_fn()
    return run_panel(
        detector_fn,
        pattern_id="P14",
        detector_format="avalanches",
        canonical_positive_runs=positives,
        canonical_metadata=metadata,
        failed_regime_module=p14_failed,
        output_path=out_path,
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


def build_p27_positives(n_seeds: int = 5) -> tuple[List[List[Dict[str, Any]]], Dict[str, Any]]:
    """P27 canonical positive: Nowak-May at b=1.8 (chaotic coexistence regime).

    L=50 grid, n_steps=200 so n_gen > 100 satisfies P27 screening prerequisite.
    """
    from epc.models.nowak_may import NowakMayModel
    runs: List[List[Dict[str, Any]]] = []
    metadata: Dict[str, Any] = {}
    for seed in range(n_seeds):
        m = NowakMayModel(rows=50, cols=50, b=1.8, init_coop_fraction=0.5, seed=seed)
        runs.append(m.run(n_steps=200))
        if seed == 0:
            metadata = m.get_metadata()
    return runs, metadata


def build_p22_positives(n_seeds: int = 5) -> tuple[List[List[Dict[str, Any]]], Dict[str, Any]]:
    """P22 canonical positive: SIR epidemic above percolation threshold.

    infection_prob=0.4 is well above threshold; single_seed init produces
    the characteristic circular wavefront and epidemic curve.
    """
    from epc.models.sir_epidemic import SIREpidemicModel
    runs: List[List[Dict[str, Any]]] = []
    metadata: Dict[str, Any] = {}
    for seed in range(n_seeds):
        m = SIREpidemicModel(
            rows=64, cols=64,
            infection_prob=0.4, recovery_prob=0.1,
            init_mode="single_seed", seed=seed,
        )
        runs.append(m.run(n_steps=200))
        if seed == 0:
            metadata = m.get_metadata()
    return runs, metadata


def _augment_history_p27(history: List[Dict[str, Any]]) -> List[Dict[str, Any]]:
    """Pass-through: Sprint 40 removed the coop_fraction augmentation.

    Previously this function computed coop_fraction = (grid == 0).mean()
    for any grid history that lacked the key, so that detect_p27 would not
    crash on generic substrates. Sprint 40 added a prerequisite guard
    directly in detect_p27 that short-circuits before accessing coop_fraction
    when the key is absent. Augmenting every grid history was the root cause
    of P27's spurious screening-tier fires on non-Nowak-May substrates
    (Sprint 39 panel finding: TNR_A=0.111). With the guard in place,
    non-Nowak-May histories are caught early and returned as not-applicable;
    augmentation is neither necessary nor correct here.
    """
    return history


def make_p27_detector_fn(n_permutations: int = 99, seed: int = 42):
    from epc.detectors.p27_spatial_reciprocity import detect_p27
    def fn(history, metadata=None):
        aug = _augment_history_p27(history)
        return detect_p27(aug, model_metadata=metadata, n_permutations=n_permutations)
    return fn


def make_p22_detector_fn(n_permutations: int = 99, seed: int = 42):
    from epc.detectors.p22_information_cascade import P22CascadeDetector
    detector = P22CascadeDetector(n_permutations=n_permutations, seed=seed)
    def fn(history, metadata=None):
        return detector.detect(history, model_metadata=metadata)
    return fn


def run_p27(out_path: str = "analysis/outputs/p27_phase2a_panel.json", verbose: bool = True) -> Dict[str, Any]:
    print(f"--- Running P27 panel → {out_path}")
    positives, metadata = build_p27_positives(n_seeds=5)
    detector_fn = make_p27_detector_fn(n_permutations=99, seed=42)
    return run_panel(
        detector_fn,
        pattern_id="P27",
        detector_format="grid",
        canonical_positive_runs=positives,
        canonical_metadata=metadata,
        failed_regime_module=p27_failed,
        output_path=out_path,
        target_steps=200,
        target_shape=(32, 32),
        verbose=verbose,
    )


def run_p22(out_path: str = "analysis/outputs/p22_phase2a_panel.json", verbose: bool = True) -> Dict[str, Any]:
    print(f"--- Running P22 panel → {out_path}")
    positives, metadata = build_p22_positives(n_seeds=5)
    detector_fn = make_p22_detector_fn(n_permutations=99, seed=42)
    return run_panel(
        detector_fn,
        pattern_id="P22",
        detector_format="grid",
        canonical_positive_runs=positives,
        canonical_metadata=metadata,
        failed_regime_module=p22_failed,
        output_path=out_path,
        target_steps=200,
        target_shape=(32, 32),
        verbose=verbose,
    )


def build_p1_positives(n_seeds: int = 5) -> tuple[List[List[Dict[str, Any]]], Dict[str, Any]]:
    """P1 canonical positive: Schelling segregation at threshold=0.375, 64×64.

    threshold=0.375 is the canonical segregating regime per Schelling (1971).
    n_steps=200 ensures the system reaches a stable segregated configuration.
    Returns a dummy metadata dict (run_schelling is a function with no get_metadata).
    """
    from epc.models.schelling import run_schelling
    runs: List[List[Dict[str, Any]]] = []
    for seed in range(n_seeds):
        runs.append(run_schelling(
            grid_size=64, density=0.9, threshold=0.375,
            n_steps=200, seed=seed,
        ))
    metadata: Dict[str, Any] = {
        "model": "schelling_segregation",
        "model_class": "schelling",
        "substrate_type": "lattice_2d",
    }
    return runs, metadata


def make_p1_detector_fn(n_permutations: int = 999):
    from epc.detectors.p1_aggregation import P1AggregationDetector
    detector = P1AggregationDetector(n_permutations=n_permutations)
    def fn(history, metadata=None):
        return detector.detect(history, model_metadata=metadata)
    return fn


def run_p1(out_path: str = "analysis/outputs/p1_phase2a_panel.json", verbose: bool = True) -> Dict[str, Any]:
    print(f"--- Running P1 panel → {out_path}")
    positives, metadata = build_p1_positives(n_seeds=5)
    detector_fn = make_p1_detector_fn(n_permutations=999)
    return run_panel(
        detector_fn,
        pattern_id="P1",
        detector_format="grid",
        canonical_positive_runs=positives,
        canonical_metadata=metadata,
        failed_regime_module=p1_failed,
        output_path=out_path,
        target_steps=200,
        target_shape=(32, 32),
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
    if which in ("p15",):
        summaries["P15"] = run_p15()
    if which in ("p14",):
        summaries["P14"] = run_p14()
    if which in ("p27",):
        summaries["P27"] = run_p27()
    if which in ("p22",):
        summaries["P22"] = run_p22()
    if which in ("p1",):
        summaries["P1"] = run_p1()

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
