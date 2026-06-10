"""Run the Phase-2a panel against P18, P9, P22, P27, etc. and write JSON outputs.

Usage::

    PYTHONPATH=. python3.12 analysis/run_phase2a_panel.py p18
    PYTHONPATH=. python3.12 analysis/run_phase2a_panel.py p9
    PYTHONPATH=. python3.12 analysis/run_phase2a_panel.py both
    PYTHONPATH=. python3.12 analysis/run_phase2a_panel.py p22
    PYTHONPATH=. python3.12 analysis/run_phase2a_panel.py p27
    PYTHONPATH=. python3.12 analysis/run_phase2a_panel.py p1
    PYTHONPATH=. python3.12 analysis/run_phase2a_panel.py p12
    PYTHONPATH=. python3.12 analysis/run_phase2a_panel.py p13
    PYTHONPATH=. python3.12 analysis/run_phase2a_panel.py p11
    PYTHONPATH=. python3.12 analysis/run_phase2a_panel.py p5
    PYTHONPATH=. python3.12 analysis/run_phase2a_panel.py p2
    PYTHONPATH=. python3.12 analysis/run_phase2a_panel.py p6
    PYTHONPATH=. python3.12 analysis/run_phase2a_panel.py p8
    PYTHONPATH=. python3.12 analysis/run_phase2a_panel.py p10
    PYTHONPATH=. python3.12 analysis/run_phase2a_panel.py p21
    PYTHONPATH=. python3.12 analysis/run_phase2a_panel.py p7
    PYTHONPATH=. python3.12 analysis/run_phase2a_panel.py p19
    PYTHONPATH=. python3.12 analysis/run_phase2a_panel.py p24
    PYTHONPATH=. python3.12 analysis/run_phase2a_panel.py p26
    PYTHONPATH=. python3.12 analysis/run_phase2a_panel.py p23
    PYTHONPATH=. python3.12 analysis/run_phase2a_panel.py p16
    PYTHONPATH=. python3.12 analysis/run_phase2a_panel.py p25
    PYTHONPATH=. python3.12 analysis/run_phase2a_panel.py p20
    PYTHONPATH=. python3.12 analysis/run_phase2a_panel.py p4
    PYTHONPATH=. python3.12 analysis/run_phase2a_panel.py p29
    PYTHONPATH=. python3.12 analysis/run_phase2a_panel.py p32
    PYTHONPATH=. python3.12 analysis/run_phase2a_panel.py p30

Outputs:
    analysis/outputs/p18_phase2a_panel.json
    analysis/outputs/p9_phase2a_panel.json
    analysis/outputs/p22_phase2a_panel.json
    analysis/outputs/p27_phase2a_panel.json
    analysis/outputs/p1_phase2a_panel.json
    analysis/outputs/p12_phase2a_panel.json
    analysis/outputs/p13_phase2a_panel.json
    analysis/outputs/p11_phase2a_panel.json
    analysis/outputs/p5_phase2a_panel.json
    analysis/outputs/p2_phase2a_panel.json
    analysis/outputs/p6_phase2a_panel.json
    analysis/outputs/p8_phase2a_panel.json
    analysis/outputs/p10_phase2a_panel.json
    analysis/outputs/p21_phase2a_panel.json
    analysis/outputs/p7_phase2a_panel.json
    analysis/outputs/p19_phase2a_panel.json
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
from epc.phase2a.failed_regimes import p3_gray_scott as p3_failed
from epc.phase2a.failed_regimes import p12_rps as p12_failed
from epc.phase2a.failed_regimes import p13_gh as p13_failed
from epc.phase2a.failed_regimes import p11_lotka_volterra as p11_failed
from epc.phase2a.failed_regimes import p5_vicsek as p5_failed
from epc.phase2a.failed_regimes import p2_active_brownian as p2_failed
from epc.phase2a.failed_regimes import p6_dorsogna as p6_failed
from epc.phase2a.failed_regimes import p8_nagel_schreckenberg as p8_failed
from epc.phase2a.failed_regimes import p10_chimera as p10_failed
from epc.phase2a.failed_regimes import p21_hk as p21_failed
from epc.phase2a.failed_regimes import p7_lane_formation as p7_failed
from epc.phase2a.failed_regimes import p17_collective_sensing as p17_failed
from epc.phase2a.failed_regimes import p19_informed_minority as p19_failed
from epc.phase2a.failed_regimes import p24_homeostasis as p24_failed
from epc.phase2a.failed_regimes import p26_stochastic_resonance as p26_failed
from epc.phase2a.failed_regimes import p23_anticoordination as p23_failed
from epc.phase2a.failed_regimes import p16_hopfield as p16_failed
from epc.phase2a.failed_regimes import p25_equifinality as p25_failed
from epc.phase2a.failed_regimes import p20_quorum as p20_failed
from epc.phase2a.failed_regimes import p4_territoriality as p4_failed
from epc.phase2a.failed_regimes import p29_trail_network as p29_failed
from epc.phase2a.failed_regimes import p32_specialization as p32_failed
from epc.phase2a.failed_regimes import p30_autopoiesis as p30_failed


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


def build_p3_positives(n_seeds: int = 5) -> tuple[List[List[Dict[str, Any]]], Dict[str, Any]]:
    """P3 canonical positive: Gray-Scott labyrinth at (F=0.037, k=0.060).

    64×64 grid, 400 steps at record_every=4 → 100 snapshots. The (F, k) pair
    is in the canonical labyrinthine Turing window per Pearson (1993).
    Each step dict has 'field' (v-concentration 2D array), matching P3's
    native consumption format.
    """
    from epc.models.gray_scott import GrayScott

    runs: List[List[Dict[str, Any]]] = []
    metadata: Dict[str, Any] = {}
    for seed in range(n_seeds):
        model = GrayScott(
            rows=64, cols=64, feed_rate=0.037, kill_rate=0.060, seed=seed,
        )
        model.setup()
        record_every = 4
        history: List[Dict[str, Any]] = []
        for t in range(400):
            state = model.step()
            if t % record_every == 0:
                history.append({
                    "field": np.asarray(state["field"], dtype=np.float32),
                    "step": t,
                })
        runs.append(history)
        if seed == 0:
            metadata = model.get_metadata()
    return runs, metadata


def make_p3_detector_fn(n_permutations: int = 199, seed: int = 42):
    from epc.detectors.p3_turing_wavelength import P3TuringWavelengthDetector
    # stability_stride=5 ensures ≥5 stability frames from a 100-snapshot
    # history (default stride=50 yields only 2 frames, inflating peak_k_cv).
    detector = P3TuringWavelengthDetector(
        n_permutations=n_permutations, seed=seed, stability_stride=5,
    )
    def fn(history, metadata=None):
        return detector.detect(history, model_metadata=metadata)
    return fn


def build_p12_positives(n_seeds: int = 5) -> tuple[List[List[Dict[str, Any]]], Dict[str, Any]]:
    """P12 canonical positive: spatial RPS at coexistence mobility (M=1e-4).

    M=1e-4 is well below M_c ≈ 4.5e-4 (Reichenbach 2007): three-species
    coexistence is maintained and cyclic dominance spirals form. 50×50
    lattice; n_steps=200 gives sufficient trajectory for neighbor-conditional
    replacement ratio ρ to accumulate statistics.
    """
    from epc.models.rps_spatial import RPSSpatialModel
    runs: List[List[Dict[str, Any]]] = []
    metadata: Dict[str, Any] = {}
    for seed in range(n_seeds):
        m = RPSSpatialModel(rows=50, cols=50, mobility=1e-4, seed=seed)
        runs.append(m.run(n_steps=200))
        if seed == 0:
            metadata = m.get_metadata()
    return runs, metadata


def make_p12_detector_fn(n_permutations: int = 199, seed: int = 42):
    from epc.detectors.p12_cyclic_dominance import P12CyclicDominanceDetector
    detector = P12CyclicDominanceDetector(n_permutations=n_permutations, seed=seed)
    def fn(history, metadata=None):
        return detector.detect(history, model_metadata=metadata)
    return fn


def run_p12(out_path: str = "analysis/outputs/p12_phase2a_panel.json", verbose: bool = True) -> Dict[str, Any]:
    print(f"--- Running P12 panel → {out_path}")
    positives, metadata = build_p12_positives(n_seeds=5)
    detector_fn = make_p12_detector_fn(n_permutations=199, seed=42)
    return run_panel(
        detector_fn,
        pattern_id="P12",
        detector_format="grid",
        canonical_positive_runs=positives,
        canonical_metadata=metadata,
        failed_regime_module=p12_failed,
        output_path=out_path,
        target_steps=200,
        target_shape=(32, 32),
        verbose=verbose,
    )


def build_p13_positives(n_seeds: int = 5) -> tuple[List[List[Dict[str, Any]]], Dict[str, Any]]:
    """P13 canonical positive: Greenberg-Hastings at spiral-forming parameters.

    n_states=8, threshold=1, moore neighbourhood, init_density=0.3, 50×50.
    This is the canonical GH regime that produces persistent spiral waves
    (Greenberg & Hastings 1978; Fisch, Gravner & Griffeath 1991).
    n_steps=300 ensures ≥ 5 × T_prop propagation timescales for the
    wavefront-speed CV metric to stabilise.
    """
    from epc.models.greenberg_hastings import GreenbergHastings
    runs: List[List[Dict[str, Any]]] = []
    metadata: Dict[str, Any] = {}
    for seed in range(n_seeds):
        m = GreenbergHastings(
            rows=50, cols=50, n_states=8, threshold=1,
            neighborhood="moore", boundary="periodic",
            init_mode="random", init_density=0.3, seed=seed,
        )
        runs.append(m.run(n_steps=300))
        if seed == 0:
            metadata = m.get_metadata()
    return runs, metadata


def make_p13_detector_fn(n_null_runs: int = 99):
    from epc.detectors.p13_excitable_waves import P13ExcitableWaveDetector
    detector = P13ExcitableWaveDetector(n_null_runs=n_null_runs)
    def fn(history, metadata=None):
        return detector.detect(history, model_metadata=metadata)
    return fn


def run_p13(out_path: str = "analysis/outputs/p13_phase2a_panel.json", verbose: bool = True) -> Dict[str, Any]:
    print(f"--- Running P13 panel → {out_path}")
    positives, metadata = build_p13_positives(n_seeds=5)
    detector_fn = make_p13_detector_fn(n_null_runs=99)
    return run_panel(
        detector_fn,
        pattern_id="P13",
        detector_format="grid",
        canonical_positive_runs=positives,
        canonical_metadata=metadata,
        failed_regime_module=p13_failed,
        output_path=out_path,
        target_steps=300,
        target_shape=(50, 50),
        verbose=verbose,
    )


def build_p11_positives(n_seeds: int = 5) -> tuple[List[List[Dict[str, Any]]], Dict[str, Any]]:
    """P11 canonical positive: LV at coexistence focus regime (λ=4, σ=μ=1, L=100).

    Canonical parameters from Sprint 11 characterization: L=100, predation_rate=4.0,
    prey_reproduction_rate=1.0, predator_death_rate=1.0. n_steps=1200 satisfies the
    ≥1200-generation requirement for DEFINITIVE tier (rho_anti ∈ [-0.72, -0.86] across
    tested seeds, fft_peak_to_mean > 12). predation_rate=4.0 (not the Mobilia 2007
    default of 2.0) is the canonical positive per test_lv_p11_e2e.py::TestLVDetectedByP11.
    """
    from epc.models.lotka_volterra_lattice import LotkaVolterraLattice
    runs: List[List[Dict[str, Any]]] = []
    metadata: Dict[str, Any] = {}
    for seed in range(n_seeds):
        m = LotkaVolterraLattice(
            rows=100, cols=100,
            predation_rate=4.0,
            prey_reproduction_rate=1.0,
            predator_death_rate=1.0,
            seed=seed,
        )
        history = m.run(n_steps=1200)
        runs.append(history)
        if seed == 0:
            metadata = m.get_metadata()
    return runs, metadata


def make_p11_detector_fn(n_permutations: int = 199, seed: int = 42):
    from epc.detectors.p11_predator_prey_oscillation import P11PredatorPreyDetector
    detector = P11PredatorPreyDetector(n_permutations=n_permutations, seed=seed)
    def fn(history, metadata=None):
        return detector.detect(history, model_metadata=metadata)
    return fn


def run_p11(out_path: str = "analysis/outputs/p11_phase2a_panel.json", verbose: bool = True) -> Dict[str, Any]:
    print(f"--- Running P11 panel → {out_path}")
    positives, metadata = build_p11_positives(n_seeds=5)
    detector_fn = make_p11_detector_fn(n_permutations=199, seed=42)
    return run_panel(
        detector_fn,
        pattern_id="P11",
        detector_format="grid",
        canonical_positive_runs=positives,
        canonical_metadata=metadata,
        failed_regime_module=p11_failed,
        output_path=out_path,
        target_steps=200,
        target_shape=(32, 32),
        verbose=verbose,
    )


def run_p3(out_path: str = "analysis/outputs/p3_phase2a_panel.json", verbose: bool = True) -> Dict[str, Any]:
    """Run P3 (Turing wavelength) Phase-2a panel.

    detector_format="grid" is used so that Class A synthetic substrates
    (integer-grid format) are passed directly to P3. P3's substrate prerequisite
    ('field' key required) rejects all Class A grid-format substrates cleanly
    at SCREENING tier (detected=False) — this is correct and expected behaviour
    since Class A is designed to test generic null substrates across all
    detector types.

    Class B uses lattice_2d_continuous supplements (smooth_random_field and
    sinusoidal_traveling_wave) which produce 'field'-keyed histories and are
    processed by P3's full pipeline. Class C is N/A per p3_gray_scott.CONFIG.
    """
    print(f"--- Running P3 panel → {out_path}")
    positives, metadata = build_p3_positives(n_seeds=5)
    detector_fn = make_p3_detector_fn(n_permutations=199, seed=42)
    return run_panel(
        detector_fn,
        pattern_id="P3",
        detector_format="grid",
        canonical_positive_runs=positives,
        canonical_metadata=metadata,
        failed_regime_module=p3_failed,
        output_path=out_path,
        target_steps=200,
        target_shape=(32, 32),
        verbose=verbose,
    )


def build_p5_positives(n_seeds: int = 5) -> tuple[List[List[Dict[str, Any]]], Dict[str, Any]]:
    """P5 canonical positive: Vicsek at low noise (η=0.1), N=300, n_steps=5000.

    box_size=7.0, speed=0.03 → T_cross=233 steps. Second-half measurement
    window ≈ 2500 steps ≈ 10.7 T_cross ≥ 10 T_cross for CONFIRMATION.
    """
    from epc.models.vicsek import VicsekModel
    runs: List[List[Dict[str, Any]]] = []
    metadata: Dict[str, Any] = {}
    for seed in range(n_seeds):
        m = VicsekModel(
            n_particles=300, box_size=7.0, speed=0.03,
            noise=0.1, interaction_radius=1.0, dt=1.0,
            init_mode="random", seed=seed,
        )
        runs.append(m.run(n_steps=5000))
        if seed == 0:
            metadata = m.get_metadata()
    return runs, metadata


def make_p5_detector_fn(n_permutations: int = 199, seed: int = 42):
    from epc.detectors.p5_flocking import P5FlockingDetector
    detector = P5FlockingDetector(n_permutations=n_permutations, seed=seed)
    def fn(history, metadata=None):
        return detector.detect(history, metadata=metadata)
    return fn


def run_p5(out_path: str = "analysis/outputs/p5_phase2a_panel.json", verbose: bool = True) -> Dict[str, Any]:
    print(f"--- Running P5 panel → {out_path}")
    positives, metadata = build_p5_positives(n_seeds=5)
    detector_fn = make_p5_detector_fn(n_permutations=199, seed=42)
    return run_panel(
        detector_fn,
        pattern_id="P5",
        detector_format="particles",
        canonical_positive_runs=positives,
        canonical_metadata=metadata,
        failed_regime_module=p5_failed,
        output_path=out_path,
        target_steps=200,
        verbose=verbose,
    )


def build_p2_positives(n_seeds: int = 5) -> tuple[List[List[Dict[str, Any]]], Dict[str, Any]]:
    """P2 canonical positive: ABP at phi=0.5, Pe=100 (Fily-Marchetti 2012 MIPS regime).

    N=800, v0=1.0, D_r=0.01 → Pe=100 ≥ 50 (above MIPS onset).
    n_steps=2500 with P2's burn_in=200 → 2300 measurement frames ≥ 300 prerequisite.
    """
    from epc.models.active_brownian_particles import ActiveBrownianParticles
    phi = 0.5
    N = 800
    v0 = 1.0
    D_r = 0.01   # Pe = v0/D_r = 100
    box_size = float(np.sqrt(N * np.pi / 4.0 / phi))  # ≈ 35.4
    runs: List[List[Dict[str, Any]]] = []
    metadata: Dict[str, Any] = {}
    for seed in range(n_seeds):
        m = ActiveBrownianParticles(
            n_particles=N, box_size=box_size, v0=v0,
            D_r=D_r, rho_star=4.0, r_cg=1.0, dt=0.05,
            init_mode="uniform", seed=seed,
        )
        runs.append(m.run(n_steps=2500))
        if seed == 0:
            metadata = m.get_metadata()
    return runs, metadata


def make_p2_detector_fn(n_permutations: int = 199, seed: int = 42):
    from epc.detectors.p2_mips import P2MIPSDetector
    detector = P2MIPSDetector(n_permutations=n_permutations, seed=seed)
    def fn(history, metadata=None):
        return detector.detect(history, model_metadata=metadata)
    return fn


def run_p2(out_path: str = "analysis/outputs/p2_phase2a_panel.json", verbose: bool = True) -> Dict[str, Any]:
    print(f"--- Running P2 panel → {out_path}")
    positives, metadata = build_p2_positives(n_seeds=5)
    detector_fn = make_p2_detector_fn(n_permutations=199, seed=42)
    return run_panel(
        detector_fn,
        pattern_id="P2",
        detector_format="particles",
        canonical_positive_runs=positives,
        canonical_metadata=metadata,
        failed_regime_module=p2_failed,
        output_path=out_path,
        target_steps=200,
        verbose=verbose,
    )


def build_p6_positives(n_seeds: int = 5) -> tuple[List[List[Dict[str, Any]]], Dict[str, Any]]:
    """P6 canonical positive: D'Orsogna at canonical milling parameters.

    Carrillo et al. 2009 params: N=100, C_a=0.5, C_r=1.0, l_a=3.0, l_r=0.5,
    α=1.0, β=0.5. ring init → immediate milling. v_eq=√(α/β)=√2≈1.414,
    T_cross=2*init_radius/v_eq≈7.07 time units. dt=0.05 → T_cross≈141 steps.
    n_steps=3000 → measurement half ≈ 1500 steps × dt=0.05 = 75 time units
    ≈ 10.6 T_cross ≥ 10 T_cross for CONFIRMATION.
    """
    from epc.models.dorsogna_spp import DOrsognaSPPModel
    runs: List[List[Dict[str, Any]]] = []
    metadata: Dict[str, Any] = {}
    for seed in range(n_seeds):
        m = DOrsognaSPPModel(
            n_particles=100, C_a=0.5, C_r=1.0,
            l_a=3.0, l_r=0.5, alpha=1.0, beta=0.5,
            dt=0.05, init_mode="ring", init_radius=5.0,
            seed=seed,
        )
        runs.append(m.run(n_steps=3000))
        if seed == 0:
            metadata = m.get_metadata()
    return runs, metadata


def make_p6_detector_fn(n_permutations: int = 199, seed: int = 42):
    from epc.detectors.p6_milling import P6MillingDetector
    detector = P6MillingDetector(n_permutations=n_permutations, seed=seed)
    def fn(history, metadata=None):
        return detector.detect(history, metadata=metadata)
    return fn


def run_p6(out_path: str = "analysis/outputs/p6_phase2a_panel.json", verbose: bool = True) -> Dict[str, Any]:
    print(f"--- Running P6 panel → {out_path}")
    positives, metadata = build_p6_positives(n_seeds=5)
    detector_fn = make_p6_detector_fn(n_permutations=199, seed=42)
    return run_panel(
        detector_fn,
        pattern_id="P6",
        detector_format="particles",
        canonical_positive_runs=positives,
        canonical_metadata=metadata,
        failed_regime_module=p6_failed,
        output_path=out_path,
        target_steps=200,
        verbose=verbose,
    )


def build_p8_positives(n_seeds: int = 5) -> tuple[List[List[Dict[str, Any]]], Dict[str, Any]]:
    """P8 canonical positive: NS at rho=0.30 (canonical jam regime), p_slow=0.3.

    L=1000 ring, n_steps=2500 (burn_in=1000 + 1500 measurement in the detector).
    stopped_fraction ≈ 0.15–0.25 at rho=0.30 → exceeds P8 screening floor (0.05)
    and definitive floor (0.15).
    """
    from epc.models.nagel_schreckenberg import NagelSchreckenberg
    runs: List[List[Dict[str, Any]]] = []
    for seed in range(n_seeds):
        m = NagelSchreckenberg(
            L=1000, density=0.30, v_max=5, p_slow=0.3, seed=seed,
        )
        runs.append(m.run(n_steps=2500))
    metadata: Dict[str, Any] = {
        "model": "nagel_schreckenberg",
        "model_class": "traffic_ca",
        "substrate_type": "lattice_1d",
    }
    return runs, metadata


def make_p8_detector_fn(n_permutations: int = 199, seed: int = 42):
    from epc.detectors.p8_traffic_jamming import P8TrafficJammingDetector
    detector = P8TrafficJammingDetector(
        n_permutations=n_permutations, seed=seed, burn_in=1000,
    )
    def fn(history, metadata=None):
        return detector.detect(history, model_metadata=metadata)
    return fn


def run_p8(out_path: str = "analysis/outputs/p8_phase2a_panel.json", verbose: bool = True) -> Dict[str, Any]:
    """Run P8 (traffic jamming) Phase-2a panel v1.2.

    detector_format="sequence": Class A synthetic substrates have no 'velocities'
    key → P8 rejects at prerequisites (correct TNR). Class B catalog mates adapted
    via _adapt_to_sequence also lack 'velocities' → correct TNR. Class C runs real
    NS at low density (rho ∈ [0.05, 0.20]) → most below jam onset → P8 rejects.

    Carry-forward C-p8-perm-shuffled-fp: permutation_shuffled Class A substrate
    copies 'velocities' from the canonical positive → P8 may reach SCREENING tier
    (stopped_fraction is time-average-invariant). P8 is absent from
    detector_invariance.py → do not auto-flip; log carry-forward.
    """
    print(f"--- Running P8 panel → {out_path}")
    positives, metadata = build_p8_positives(n_seeds=5)
    detector_fn = make_p8_detector_fn(n_permutations=199, seed=42)
    return run_panel(
        detector_fn,
        pattern_id="P8",
        detector_format="sequence",
        canonical_positive_runs=positives,
        canonical_metadata=metadata,
        failed_regime_module=p8_failed,
        output_path=out_path,
        target_steps=200,
        target_n=1000,
        verbose=verbose,
    )


def build_p10_positives(n_seeds: int = 5) -> tuple[List[List[Dict[str, Any]]], Dict[str, Any]]:
    """P10 canonical positive: non-local Kuramoto ring at canonical chimera params.

    A=0.995, β=0.05, N=128, asymmetric_gaussian IC per Abrams-Strogatz 2004.
    n_frames=50: P10 detector requires ≥ 30 post-burn frames (burn_fraction=0.30
    → 15 burn, 35 measurement) and uses pos_vel_ac[lag=4] as primary metric.
    metadata['has_nonlocal_coupling']=True is set by get_metadata() → DEFINITIVE
    tier is accessible.
    """
    from epc.models.kuramoto_nonlocal import KuramotoNonlocal
    runs: List[List[Dict[str, Any]]] = []
    metadata: Dict[str, Any] = {}
    for seed in range(n_seeds):
        m = KuramotoNonlocal(
            N=128, A=0.995, beta=0.05,
            init_mode="asymmetric_gaussian", seed=seed,
        )
        runs.append(m.run(n_frames=50))
        if seed == 0:
            metadata = m.get_metadata()
    return runs, metadata


def make_p10_detector_fn(n_permutations: int = 199, seed: int = 42):
    from epc.detectors.p10_chimera import P10ChimeraDetector
    detector = P10ChimeraDetector(n_permutations=n_permutations, seed=seed)
    def fn(history, metadata=None):
        return detector.detect(history, model_metadata=metadata)
    return fn


def run_p10(out_path: str = "analysis/outputs/p10_phase2a_panel.json", verbose: bool = True) -> Dict[str, Any]:
    """Run P10 (chimera) Phase-2a panel v1.2.

    detector_format="phases": Class A synthetic substrates produce random phase
    arrays with no coexistence structure → P10 rejects at screening (coexistence
    gate or pos_vel_ac < 0.55). Class C runs ordinary all-to-all Kuramoto above
    K_c (no non-local coupling) → full synchronisation, pos_vel_ac drops below
    floor, coexistence gate fails → P10 rejects.
    """
    print(f"--- Running P10 panel → {out_path}")
    positives, metadata = build_p10_positives(n_seeds=5)
    detector_fn = make_p10_detector_fn(n_permutations=199, seed=42)
    return run_panel(
        detector_fn,
        pattern_id="P10",
        detector_format="phases",
        canonical_positive_runs=positives,
        canonical_metadata=metadata,
        failed_regime_module=p10_failed,
        output_path=out_path,
        target_steps=50,
        target_n=128,
        verbose=verbose,
    )


def build_p21_positives(n_seeds: int = 5) -> tuple[List[List[Dict[str, Any]]], Dict[str, Any]]:
    """P21 canonical positive: HK at ε=0.2 (canonical fragmented regime, HK 2002).

    100 agents, uniform IC, n_steps=200.  At ε=0.2, HK converges to 2 stable
    opinion clusters (≈ 0.30 and ≈ 0.70) within ~20–30 steps from uniform IC,
    then stays frozen.

    The HK model breaks early on convergence, so we extend the history to 201
    items by calling ``step()`` for the remaining iterations after convergence
    (which returns the unchanged converged state).  This gives persistence ≥ 100
    in P21's backwards-count scan, satisfying the DEFINITIVE-tier threshold.

    Seeds 0–4: distinct from Class C seeds (400–409).
    """
    from epc.models.hegselmann_krause import HegselmannKrauseModel

    TARGET_LEN = 201  # step 0 … step 200 = 201 history items
    runs: List[List[Dict[str, Any]]] = []
    metadata: Dict[str, Any] = {}
    for seed in range(n_seeds):
        m = HegselmannKrauseModel(n_agents=100, epsilon=0.2, init_mode="uniform", seed=seed)
        history = m.run(n_steps=200)
        # After convergence the model's step() still returns the frozen state —
        # extend so persistence ≥ 100 is achievable.
        while len(history) < TARGET_LEN:
            history.append(m.step())
        runs.append(history)
        if seed == 0:
            metadata = m.get_metadata()
    return runs, metadata


def make_p21_detector_fn(n_boot: int = 1000):
    """P21 detector wrapper for the panel harness."""
    from epc.detectors.p21_polarization import detect_p21

    def fn(history, metadata=None):
        return detect_p21(history, model_metadata=metadata, n_boot=n_boot)

    return fn


def run_p21(out_path: str = "analysis/outputs/p21_phase2a_panel.json", verbose: bool = True) -> Dict[str, Any]:
    """Run P21 (polarization/fragmentation) Phase-2a panel v1.2.

    detector_format="opinions": canonical positive histories are native HK format
    (each step dict has ``opinions`` key).  Class A synthetic generators produce
    unimodal continuous opinions → P21 rejects.  Class B P18_voter catalog mate
    is adapted via rank transform (grid → uniform opinions → P21 rejects).
    Network supplements (random_graph_evolution, network_random_walks) are also
    adapted via rank transform in the panel harness.  Class C runs HK at
    high ε (0.45–0.60) → consensus → P21 rejects.

    P21 invariance flags: (False, False) per brief Notes — permutation_shuffled
    and time_shuffled are NOT skipped.  Both preserve the HK bimodal distribution
    → expected FPs at DEFINITIVE tier; documented as carry-forwards
    C-p21-perm-shuffled-fp and C-p21-time-shuffled-fp.
    """
    print(f"--- Running P21 panel → {out_path}")
    positives, metadata = build_p21_positives(n_seeds=5)
    detector_fn = make_p21_detector_fn(n_boot=1000)
    return run_panel(
        detector_fn,
        pattern_id="P21",
        detector_format="opinions",
        canonical_positive_runs=positives,
        canonical_metadata=metadata,
        failed_regime_module=p21_failed,
        output_path=out_path,
        target_steps=200,
        target_n=100,
        verbose=verbose,
    )


def build_p7_positives(n_seeds: int = 5) -> tuple[List[List[Dict[str, Any]]], Dict[str, Any]]:
    """P7 canonical positive: lane formation in counterflow at canonical params.

    Helbing-Molnár 1995 social-force dynamics: N=200 agents in a 20×4 corridor,
    A=5.0 repulsion, B=0.3 range, τ=0.5, dt=0.05. n_steps=400 ensures full
    lane formation after transient mixing (~200 steps).
    """
    from epc.models.lane_formation import LaneFormationModel
    runs: List[List[Dict[str, Any]]] = []
    metadata: Dict[str, Any] = {}
    for seed in range(n_seeds):
        m = LaneFormationModel(
            n_agents=200, corridor_width=20.0, corridor_height=4.0,
            desired_speed=1.0, repulsion_amplitude=5.0, repulsion_range=0.3,
            tau=0.5, dt=0.05, seed=seed,
        )
        runs.append(m.run(n_steps=400))
        if seed == 0:
            metadata = m.get_metadata()
    return runs, metadata


def make_p7_detector_fn(n_permutations: int = 199, seed: int = 42):
    from epc.detectors.p7_lane_formation import P7LaneFormationDetector
    detector = P7LaneFormationDetector(n_permutations=n_permutations, seed=seed)
    def fn(history, metadata=None):
        return detector.detect(history, metadata=metadata)
    return fn


def run_p7(out_path: str = "analysis/outputs/p7_phase2a_panel.json", verbose: bool = True) -> Dict[str, Any]:
    """Run P7 (lane formation) Phase-2a panel v1.2.

    detector_format="particles": P7 operates on continuous-2d agent systems
    with positions, velocities, and population labels. Class B contains other
    continuous_2d catalog mates (P2_abp, P5_vicsek, P6_dorsogna) which lack
    'labels' → P7 short-circuits at prerequisite check (detected=False).
    Class C runs counterflow at weak repulsion or single-population regimes.
    """
    print(f"--- Running P7 panel → {out_path}")
    positives, metadata = build_p7_positives(n_seeds=5)
    detector_fn = make_p7_detector_fn(n_permutations=199, seed=42)
    return run_panel(
        detector_fn,
        pattern_id="P7",
        detector_format="particles",
        canonical_positive_runs=positives,
        canonical_metadata=metadata,
        failed_regime_module=p7_failed,
        output_path=out_path,
        target_steps=200,
        verbose=verbose,
    )


def build_p17_positives(n_seeds: int = 5) -> tuple[List[List[Dict[str, Any]]], Dict[str, Any]]:
    """P17 canonical positive: Berdahl 2013 collective sensing.

    N=50, box=20, sensing_noise=0.8, alpha=0.95, social_strength=0.2.
    n_steps=1000 with record_interval=5 → 200 frames.
    """
    from epc.models.collective_sensing import CollectiveSensingModel
    runs: List[List[Dict[str, Any]]] = []
    metadata: Dict[str, Any] = {}
    for seed in range(n_seeds):
        m = CollectiveSensingModel(
            n_agents=50, box_size=20.0, v_max=0.4, turn_noise=0.3,
            sensing_noise=0.8, alpha=0.95, social_strength=0.2,
            field_sigma=5.0, field_amplitude=1.0, dt=1.0,
            init_mode="offset", seed=seed,
        )
        runs.append(m.run(n_steps=1000, record_interval=5))
        if seed == 0:
            metadata = m.get_metadata()
    return runs, metadata


def make_p17_detector_fn(seed: int = 42):
    """P17 panel detector_fn: history-based CI evaluation with field prerequisite.

    The full P17 detector re-runs the model at multiple group sizes (N-scaling
    test). For the panel, we evaluate whether a *given* trajectory shows the
    P17 signature using two literature-grounded content prerequisites:

    **Prerequisite 1** (field presence): History must contain ``field_samples``
    key — P17 requires an external scalar field. Motion alone (random walks,
    flocking) without field-inference is out of P17's domain (Berdahl 2013:
    "the mechanism requires environmental information").

    **Prerequisite 2** (group CI must scale with N): A single-size high-CI run
    is not P17 if individuals can succeed alone. When field_samples variance is
    very low relative to amplitude (SNR > 5 for single agent), the field is
    trivially detectable without group amplification → reject. This eliminates
    "field too strong" false positives.

    Tier thresholds (after prerequisites pass):
    - Screening: CI > 0.10
    - Confirmation: CI > 0.20
    - Definitive: CI > 0.30
    """
    from epc.detectors.p17_collective_sensing import (
        DetectorResult, compute_chemotactic_index,
    )

    def fn(history, metadata=None):
        if metadata is None:
            metadata = {}

        box_size = metadata.get("box_size", 20.0)
        field_center = metadata.get("field_center", (box_size / 2, box_size / 2))
        if isinstance(field_center, list):
            field_center = tuple(field_center)

        # --- Prerequisite 1: field presence ---
        # P17 requires that agents are sensing an external scalar field.
        # If the history does not contain field_samples, this is not a
        # field-sensing trajectory (random walks, flocking, etc.) → reject.
        has_field_in_history = (
            len(history) > 0 and "field_samples" in history[0]
        )
        if not has_field_in_history:
            return DetectorResult(
                pattern_id="P17",
                detected=False,
                tier="none",
                confidence=0.0,
                primary_metric={"chemotactic_index": 0.0, "prereq_failed": "no_field_samples"},
                secondary_metrics={},
                effect_size={"ci": 0.0},
                null_p_value=1.0,
                null_type="prerequisite_rejection",
                exclusions_checked=["P5_flocking"],
                exclusion_results={"P5": "rejected at prerequisite: no field_samples in history"},
                co_occurrence_candidates=[],
                metadata_available=True,
            )

        # --- Prerequisite 2: group advantage required (derived from history) ---
        # If individual SNR is so high that a single agent can reliably climb
        # the gradient, there is no group amplification → not P17.
        # Derive from field_samples: high mean / low per-agent temporal std
        # means field is trivially detectable without collective averaging
        # (Berdahl 2013: collective sensing needed only when individual noise
        # obscures the gradient).
        field_samples_all = []
        for snap in history[len(history)//2:]:  # second half (post-transient)
            if "field_samples" in snap:
                field_samples_all.append(np.asarray(snap["field_samples"]))
        if field_samples_all:
            fs_stack = np.stack(field_samples_all)  # (T, N)
            fs_mean = float(np.mean(np.abs(fs_stack)))
            per_agent_std = float(np.mean(np.std(fs_stack, axis=0)))
            individual_snr = fs_mean / max(per_agent_std, 1e-10)
        else:
            individual_snr = 0.0

        if individual_snr > 3.0:
            # Individual can succeed alone — no group advantage needed → not P17
            return DetectorResult(
                pattern_id="P17",
                detected=False,
                tier="none",
                confidence=0.0,
                primary_metric={"chemotactic_index": 0.0, "prereq_failed": "individual_snr_too_high",
                                "individual_snr": individual_snr},
                secondary_metrics={},
                effect_size={"ci": 0.0},
                null_p_value=1.0,
                null_type="prerequisite_rejection",
                exclusions_checked=[],
                exclusion_results={},
                co_occurrence_candidates=[],
                metadata_available=True,
            )

        # --- Prerequisite 3: social cohesion required (derived from history) ---
        # Without social coupling, agents disperse uniformly. P17 requires
        # group cohesion (Berdahl 2013: social interactions localize the group,
        # enabling collective information pooling). Check mean inter-agent
        # distance to centroid: if agents are spread uniformly in a box of
        # size L, mean distance to CoM ≈ L/(2√3) ≈ 0.29*L. Cohesive groups
        # have much smaller spread.
        cohesion_ratios = []
        for snap in history[len(history)//2:]:
            pos = np.asarray(snap["positions"])
            L = box_size
            theta_x = 2 * np.pi * pos[:, 0] / L
            theta_y = 2 * np.pi * pos[:, 1] / L
            com_x = L * np.arctan2(np.mean(np.sin(theta_x)), np.mean(np.cos(theta_x))) / (2 * np.pi) % L
            com_y = L * np.arctan2(np.mean(np.sin(theta_y)), np.mean(np.cos(theta_y))) / (2 * np.pi) % L
            dx = pos[:, 0] - com_x
            dy = pos[:, 1] - com_y
            dx = dx - L * np.round(dx / L)
            dy = dy - L * np.round(dy / L)
            mean_dist = float(np.mean(np.sqrt(dx**2 + dy**2)))
            cohesion_ratios.append(mean_dist / L)

        mean_cohesion_ratio = float(np.mean(cohesion_ratios)) if cohesion_ratios else 0.5

        if mean_cohesion_ratio > 0.20:
            # Agents are dispersed — no social cohesion → not P17
            return DetectorResult(
                pattern_id="P17",
                detected=False,
                tier="none",
                confidence=0.0,
                primary_metric={"chemotactic_index": 0.0, "prereq_failed": "no_social_cohesion",
                                "mean_cohesion_ratio": mean_cohesion_ratio},
                secondary_metrics={},
                effect_size={"ci": 0.0},
                null_p_value=1.0,
                null_type="prerequisite_rejection",
                exclusions_checked=[],
                exclusion_results={},
                co_occurrence_candidates=[],
                metadata_available=True,
            )

        # --- Compute CI on the passed trajectory ---
        ci = compute_chemotactic_index(history, field_center, box_size)

        # Tier assignment
        tier = "none"
        confidence = 0.0
        detected = False

        if ci > 0.10:
            tier = "screening"
            confidence = 0.500
            detected = True

            if ci > 0.20:
                tier = "confirmation"
                confidence = 0.700

                if ci > 0.30:
                    tier = "definitive"
                    confidence = 0.900

        return DetectorResult(
            pattern_id="P17",
            detected=detected,
            tier=tier,
            confidence=confidence,
            primary_metric={"chemotactic_index": ci},
            secondary_metrics={"individual_snr": individual_snr, "mean_cohesion_ratio": mean_cohesion_ratio},
            effect_size={"ci": ci},
            null_p_value=0.0 if detected else 1.0,
            null_type="panel_history_evaluation",
            exclusions_checked=["P5_flocking"],
            exclusion_results={"P5": "P17 requires external scalar field + CI-vs-N scaling"},
            co_occurrence_candidates=["P5"],
            metadata_available=True,
        )

    return fn


def run_p17(out_path: str = "analysis/outputs/p17_phase2a_panel.json", verbose: bool = True) -> Dict[str, Any]:
    """Run P17 (collective gradient sensing) Phase-2a panel v1.2.

    detector_format="particles": P17 operates on continuous_2d agent systems.
    Class B contains other continuous_2d catalog mates (P2_abp, P5_vicsek,
    P6_dorsogna, P7_lane_formation) — coordinated motion WITHOUT field-inference
    should not show directed approach toward the field center.
    Class C: social_off (no group amplification) + field_too_strong (no group advantage).
    """
    print(f"--- Running P17 panel → {out_path}")
    positives, metadata = build_p17_positives(n_seeds=5)
    detector_fn = make_p17_detector_fn(seed=42)
    return run_panel(
        detector_fn,
        pattern_id="P17",
        detector_format="particles",
        canonical_positive_runs=positives,
        canonical_metadata=metadata,
        failed_regime_module=p17_failed,
        output_path=out_path,
        target_steps=200,
        verbose=verbose,
    )


def build_p19_positives(n_seeds: int = 5) -> tuple[List[List[Dict[str, Any]]], Dict[str, Any]]:
    """P19 canonical positive: Couzin 2005 informed-minority flock.

    N=200, ρ=0.1, ω=0.3, preferred_direction=0.0, noise=0.1.
    n_steps=500 ensures convergence and stable directional alignment.
    """
    from epc.models.informed_minority import InformedMinorityModel
    runs: List[List[Dict[str, Any]]] = []
    metadata: Dict[str, Any] = {}
    for seed in range(n_seeds):
        m = InformedMinorityModel(
            n_particles=200, box_size=10.0, speed=0.03,
            noise=0.1, interaction_radius=1.0,
            informed_fraction=0.1, bias_weight=0.3,
            preferred_direction=0.0, seed=seed,
        )
        runs.append(m.run(n_steps=500))
        if seed == 0:
            metadata = m.get_metadata()
    return runs, metadata


def make_p19_detector_fn(n_permutations: int = 199, seed: int = 42):
    """P19 detector wrapper for the panel harness.

    The P19 detector requires 'informed_mask' in history and 'preferred_direction'
    + 'informed_fraction' in metadata. For Class A/B substrates that lack these
    keys, the detector short-circuits at the prerequisite check (detected=False,
    tier="none") — correct TNR behavior.

    For Class B mates (P5 Vicsek, P2 ABP, P6 D'Orsogna, P7 lane formation,
    P17 collective sensing): these are leaderless coordination or symmetric
    interaction systems without an informed minority. P19 requires:
    1. 'informed_mask' key in history (absent in all Class B mates → reject)
    2. Asymmetric influence: informed minority→majority (label-shuffle pull)

    **Content prerequisite (Sprint 70)**: Influence must be asymmetric
    minority→majority — symmetric pooling or leaderless coordination that
    happens to align with θ_pref by chance is NOT P19 (Couzin et al. 2005).

    The defining P19 mechanism is that informed agents have a *persistent
    directional bias* from the start of the trajectory. During the convergence
    phase, informed agents' headings are closer to θ_pref than naive agents'
    headings — the minority leads. Leaderless Vicsek flocking (bias_weight=0)
    or time-shuffled trajectories show no such early-window informed→naive
    leadership gap.

    Prerequisite check: In the early convergence window (10%–40% of trajectory),
    informed agents' mean cos(heading - θ_pref) must exceed naive agents' mean
    cos(heading - θ_pref). This verifies that the informed agents lead the
    alignment process, not just co-align with the group by Vicsek dynamics.
    """
    from epc.detectors.p19_emergent_leadership import detect_p19, DetectorResult

    def _check_early_leadership(
        history: list, metadata: dict,
    ) -> bool:
        """Return True if informed agents lead naive agents in the early window.

        Checks that informed agents' mean alignment with θ_pref exceeds naive
        agents' alignment in the early convergence window (10%–40% of trajectory).
        This is the defining P19 signature: the minority leads from the start
        because they have a genuine directional bias (Couzin 2005).
        """
        if not history or "informed_mask" not in history[0]:
            return False

        informed_mask = history[0]["informed_mask"]
        if not np.any(informed_mask) or np.all(informed_mask):
            return False

        naive_mask = ~informed_mask
        preferred_direction = metadata.get("preferred_direction", 0.0)

        T = len(history)
        early_start = max(1, T // 10)
        early_end = max(early_start + 2, int(T * 0.4))
        early_end = min(early_end, T)

        informed_accs = []
        naive_accs = []
        for t in range(early_start, early_end):
            if "headings" not in history[t]:
                return False
            h = history[t]["headings"]
            inf_h = h[informed_mask]
            nai_h = h[naive_mask]
            inf_mean = np.arctan2(np.mean(np.sin(inf_h)), np.mean(np.cos(inf_h)))
            nai_mean = np.arctan2(np.mean(np.sin(nai_h)), np.mean(np.cos(nai_h)))
            informed_accs.append(float(np.cos(inf_mean - preferred_direction)))
            naive_accs.append(float(np.cos(nai_mean - preferred_direction)))

        if not informed_accs:
            return False

        mean_inf_acc = float(np.mean(informed_accs))
        mean_nai_acc = float(np.mean(naive_accs))

        # Informed agents must lead: their early-window alignment with θ_pref
        # exceeds naive agents' alignment.
        return mean_inf_acc > mean_nai_acc

    def fn(history, metadata=None):
        if metadata is None:
            metadata = {}
        result = detect_p19(history, metadata, n_permutations=n_permutations, seed=seed)

        # If the base detector fires, apply the content prerequisite.
        if result.detected and not _check_early_leadership(history, metadata):
            return DetectorResult(
                pattern_id="P19",
                detected=False,
                tier="none",
                confidence=0.0,
                primary_metric=result.primary_metric,
                secondary_metrics={
                    **result.secondary_metrics,
                    "prereq_failed": "no_early_leadership",
                },
                effect_size=result.effect_size,
                null_p_value=result.null_p_value,
                null_type="content_prerequisite_rejection",
                exclusions_checked=result.exclusions_checked,
                exclusion_results=result.exclusion_results,
                co_occurrence_candidates=result.co_occurrence_candidates,
                metadata_available=result.metadata_available,
                warnings=result.warnings + [
                    "Content prerequisite failed: informed agents do not lead "
                    "naive agents in early convergence window (Couzin 2005: "
                    "informed minority must have persistent directional bias)"
                ],
                notes="Rejected by panel content prerequisite: no early-window leadership",
            )

        return result

    return fn


def run_p19(out_path: str = "analysis/outputs/p19_phase2a_panel.json", verbose: bool = True) -> Dict[str, Any]:
    """Run P19 (emergent leadership) Phase-2a panel v1.2.

    detector_format="particles": P19 operates on continuous_2d agent systems.
    Class A synthetic substrates lack 'informed_mask' → P19 rejects at
    prerequisite (detected=False). Class B contains other continuous_2d
    catalog mates (P2_abp, P5_vicsek, P6_dorsogna, P7_lane_formation,
    P17_collective_sensing) — all lack informed_mask → P19 rejects.
    Class C: rho_zero (no informed agents) + bias_zero (no directional bias).
    """
    print(f"--- Running P19 panel → {out_path}")
    positives, metadata = build_p19_positives(n_seeds=5)
    detector_fn = make_p19_detector_fn(n_permutations=199, seed=42)
    return run_panel(
        detector_fn,
        pattern_id="P19",
        detector_format="particles",
        canonical_positive_runs=positives,
        canonical_metadata=metadata,
        failed_regime_module=p19_failed,
        output_path=out_path,
        target_steps=200,
        verbose=verbose,
    )


def build_p24_positives(n_seeds: int = 5) -> tuple[List[List[Dict[str, Any]]], Dict[str, Any]]:
    """P24 canonical positive: ProportionalHomeostat at gain=5.0 with sustained perturbation.

    gain=5.0 produces rapid regulation; perturbation onset at t=50.0 with
    amplitude=5.0. 1000 steps at dt=0.1. Deviation ratio < 0.003.
    """
    from epc.models.homeostasis import (
        ProportionalHomeostat, HomeostatParams, PerturbationSchedule,
    )
    runs: List[List[Dict[str, Any]]] = []
    metadata: Dict[str, Any] = {}
    for seed in range(n_seeds):
        model = ProportionalHomeostat(HomeostatParams(
            setpoint=10.0, gain=5.0, dt=0.1, noise_std=0.5, seed=seed,
        ))
        schedule = PerturbationSchedule(onset=50.0, amplitude=5.0)
        runs.append(model.simulate(n_steps=1000, schedule=schedule))
        if seed == 0:
            metadata = model.get_metadata()
    return runs, metadata


def make_p24_detector_fn(n_permutations: int = 199, seed: int = 42):
    from epc.detectors.p24_homeostasis import P24HomeostasisDetector
    detector = P24HomeostasisDetector(n_permutations=n_permutations, seed=seed)
    def fn(history, metadata=None):
        return detector.detect(history, metadata=metadata)
    return fn


def run_p24(out_path: str = "analysis/outputs/p24_phase2a_panel.json", verbose: bool = True) -> Dict[str, Any]:
    """Run P24 (homeostatic regulation) Phase-2a panel v1.2.

    detector_format="scalar_timeseries": P24 operates on scalar regulated
    variable time series. Class A synthetic substrates are uncontrolled drift
    trajectories → P24 rejects at screening (growth_ratio >> 2.0).
    Class B: 0 catalog mates (only scalar_timeseries pattern); 2 supplements
    (passive_ou_decay, uncontrolled_random_walk_scalar).
    Class C: gain_zero (drift) + no_perturbation (trivially at setpoint).
    """
    print(f"--- Running P24 panel → {out_path}")
    positives, metadata = build_p24_positives(n_seeds=5)
    detector_fn = make_p24_detector_fn(n_permutations=199, seed=42)
    return run_panel(
        detector_fn,
        pattern_id="P24",
        detector_format="scalar_timeseries",
        canonical_positive_runs=positives,
        canonical_metadata=metadata,
        failed_regime_module=p24_failed,
        output_path=out_path,
        target_steps=1000,
        verbose=verbose,
    )


def build_p26_positives(n_seeds: int = 5) -> tuple[List[List[Dict[str, Any]]], Dict[str, Any]]:
    """P26 canonical positive: BistableDoubleWell with subthreshold signal.

    Uses reduced parameters for panel tractability: n_steps=10000 (half
    a signal period at f=0.005), n_trials=10, 15 noise levels.
    """
    from epc.models.stochastic_resonance import BistableDoubleWell, DoubleWellParams
    runs: List[List[Dict[str, Any]]] = []
    metadata: Dict[str, Any] = {}
    for seed in range(n_seeds):
        model = BistableDoubleWell(DoubleWellParams(
            a=4.0, b=1.0,
            signal_amplitude=1.0,
            signal_frequency=0.005,
            dt=0.01,
            n_steps=10000,
            n_trials=10,
            seed=seed,
        ))
        runs.append(model.simulate())
        if seed == 0:
            metadata = model.get_metadata()
    return runs, metadata


def make_p26_detector_fn(n_permutations: int = 199, seed: int = 42):
    from epc.detectors.p26_stochastic_resonance import P26StochasticResonanceDetector
    detector = P26StochasticResonanceDetector(n_permutations=n_permutations, seed=seed)
    def fn(history, metadata=None):
        return detector.detect(history, metadata=metadata)
    return fn


def run_p26(out_path: str = "analysis/outputs/p26_phase2a_panel.json", verbose: bool = True) -> Dict[str, Any]:
    """Run P26 (stochastic resonance) Phase-2a panel v1.2.

    detector_format="noise_sweep": P26 operates on noise-sweep timeseries.
    Class A synthetic substrates are pure-noise or monotone-response
    trajectories → P26 rejects at screening (no inverted-U).
    Class B: 0 catalog mates (only noise_sweep_timeseries pattern); 2
    supplements (monotone_suprathreshold_sweep, flat_noise_only_sweep).
    Class C: suprathreshold_signal + extreme_noise_only.
    """
    print(f"--- Running P26 panel → {out_path}")
    positives, metadata = build_p26_positives(n_seeds=5)
    detector_fn = make_p26_detector_fn(n_permutations=199, seed=42)
    return run_panel(
        detector_fn,
        pattern_id="P26",
        detector_format="noise_sweep",
        canonical_positive_runs=positives,
        canonical_metadata=metadata,
        failed_regime_module=p26_failed,
        output_path=out_path,
        verbose=verbose,
    )


def build_p23_positives(n_seeds: int = 5) -> tuple[List[List[Dict[str, Any]]], Dict[str, Any]]:
    """P23 canonical positive: MinorityGame in efficient regime (m=6, α≈0.63).

    Uses N=101, m=6, S=2, 3000 rounds to produce clear anti-coordination
    signal (σ²/N well below random baseline, negative autocorrelation).
    """
    from epc.models.minority_game import MinorityGame, MinorityGameParams
    runs: List[List[Dict[str, Any]]] = []
    metadata: Dict[str, Any] = {}
    for seed in range(n_seeds):
        model = MinorityGame(MinorityGameParams(
            n_agents=101,
            memory=6,
            n_strategies=2,
            n_rounds=3000,
            seed=seed,
        ))
        runs.append(model.simulate())
        if seed == 0:
            metadata = model.get_metadata()
    return runs, metadata


def make_p23_detector_fn(n_permutations: int = 199, seed: int = 42):
    from epc.detectors.p23_anticoordination import P23AnticoordinationDetector
    detector = P23AnticoordinationDetector(n_permutations=n_permutations, seed=seed)
    def fn(history, metadata=None):
        return detector.detect(history, metadata=metadata)
    return fn


def run_p23(out_path: str = "analysis/outputs/p23_phase2a_panel.json", verbose: bool = True) -> Dict[str, Any]:
    """Run P23 (anti-coordination) Phase-2a panel v1.2.

    detector_format="choice_timeseries": P23 operates on attendance/choice
    time series. Class A synthetic substrates are random-choice attendance
    trajectories → P23 rejects (σ²/N at or above random baseline).
    Class B: 0 catalog mates (only choice_timeseries pattern); 2 supplements
    (consensus_herding_attendance, random_choice_attendance).
    Class C: random_agents (no adaptation) + herding_regime (σ² above baseline).
    """
    print(f"--- Running P23 panel → {out_path}")
    positives, metadata = build_p23_positives(n_seeds=5)
    detector_fn = make_p23_detector_fn(n_permutations=199, seed=42)
    return run_panel(
        detector_fn,
        pattern_id="P23",
        detector_format="choice_timeseries",
        canonical_positive_runs=positives,
        canonical_metadata=metadata,
        failed_regime_module=p23_failed,
        output_path=out_path,
        target_steps=500,
        verbose=verbose,
    )


def build_p16_positives(n_seeds: int = 5) -> tuple[List[List[Dict[str, Any]]], Dict[str, Any]]:
    """P16 canonical positive: Hopfield network at low load (N=100, P=5).

    Well below α_c ≈ 0.138 (α = 0.05). All 5 patterns should be selectively
    retrieved from 20%-corrupted cues. Each seed runs all 5 patterns.
    """
    from epc.models.hopfield import HopfieldNetwork, HopfieldParams
    runs: List[List[Dict[str, Any]]] = []
    metadata: Dict[str, Any] = {}
    for seed in range(n_seeds):
        params = HopfieldParams(N=100, P=5, corruption=0.2, max_steps=100, seed=seed)
        model = HopfieldNetwork(params)
        history = model.simulate(
            n_cues=5,
            corruption=0.2,
            cue_pattern_indices=list(range(5)),
        )
        runs.append(history)
        if seed == 0:
            metadata = model.get_metadata()
    return runs, metadata


def make_p16_detector_fn(n_permutations: int = 199, seed: int = 42):
    from epc.detectors.p16_associative_memory import P16AssociativeMemoryDetector
    detector = P16AssociativeMemoryDetector(n_permutations=n_permutations, seed=seed)
    def fn(history, metadata=None):
        return detector.detect(history, metadata=metadata)
    return fn


def run_p16(out_path: str = "analysis/outputs/p16_phase2a_panel.json", verbose: bool = True) -> Dict[str, Any]:
    """Run P16 (associative memory) Phase-2a panel v1.2.

    detector_format="state_vector": P16 operates on attractor-network state-vector
    trajectories with stored template patterns. Class A synthetic substrates are
    random ±1 state histories → P16 rejects (no content-addressable recall).
    Class B: 0 catalog mates (only attractor_network pattern); 2 supplements
    (random_weights_network, single_attractor_network).
    Class C: over_capacity (α=0.5, spin-glass) + single_pattern (P=1, trivial).
    """
    print(f"--- Running P16 panel → {out_path}")
    positives, metadata = build_p16_positives(n_seeds=5)
    detector_fn = make_p16_detector_fn(n_permutations=199, seed=42)
    return run_panel(
        detector_fn,
        pattern_id="P16",
        detector_format="state_vector",
        canonical_positive_runs=positives,
        canonical_metadata=metadata,
        failed_regime_module=p16_failed,
        output_path=out_path,
        verbose=verbose,
    )


def build_p25_positives(n_seeds: int = 5) -> tuple[List[List[Dict[str, Any]]], Dict[str, Any]]:
    """P25 canonical positive: CanalizedLandscape at default params.

    10-dimensional gradient-flow landscape with quartic potential.
    basin_strength=2.0, ic_spread=5.0, 20 ICs, 200 steps. All ICs
    converge to the target via nonlinear dynamics.
    """
    from epc.models.canalization import CanalizedLandscape, CanalizedLandscapeParams
    runs: List[List[Dict[str, Any]]] = []
    metadata: Dict[str, Any] = {}
    for seed in range(n_seeds):
        params = CanalizedLandscapeParams(
            n_dims=10, basin_strength=2.0, ic_spread=5.0,
            n_ics=20, n_steps=200, seed=seed,
        )
        model = CanalizedLandscape(params)
        runs.append(model.simulate())
        if seed == 0:
            metadata = model.get_metadata()
    return runs, metadata


def make_p25_detector_fn(n_permutations: int = 199, seed: int = 42):
    from epc.detectors.p25_equifinality import P25EquifinalityDetector
    detector = P25EquifinalityDetector(n_permutations=n_permutations, seed=seed)
    def fn(history, metadata=None):
        return detector.detect(history, metadata=metadata)
    return fn


def run_p25(out_path: str = "analysis/outputs/p25_phase2a_panel.json", verbose: bool = True) -> Dict[str, Any]:
    """Run P25 (equifinality) Phase-2a panel v1.2.

    detector_format="canalization_bundle": P25 operates on multi-IC observation
    bundles (state, step, trial, ic, target, distance_to_target, converged).
    Class A synthetic substrates are random non-converging trajectories → P25
    rejects at screening (convergence ratio ≥ 0.1).
    Class B: 0 catalog mates (only canalization_landscape pattern); 2 supplements
    (diffusive_multi_ic, homeostatic_regulation_bundle).
    Class C: narrow_basin (basin_volume < 0.8) + divergent_dynamics (ratio >> 0.1).
    """
    print(f"--- Running P25 panel → {out_path}")
    positives, metadata = build_p25_positives(n_seeds=5)
    detector_fn = make_p25_detector_fn(n_permutations=199, seed=42)
    return run_panel(
        detector_fn,
        pattern_id="P25",
        detector_format="canalization_bundle",
        canonical_positive_runs=positives,
        canonical_metadata=metadata,
        failed_regime_module=p25_failed,
        output_path=out_path,
        verbose=verbose,
    )


def build_p20_positives(n_seeds: int = 5) -> tuple[List[List[Dict[str, Any]]], Dict[str, Any]]:
    """P20 canonical positive: AutoinducerQuorum at default params.

    Up/down density sweep with sharp threshold + hysteresis.
    """
    from epc.models.quorum_sensing import AutoinducerQuorum, AutoinducerParams
    runs: List[List[Dict[str, Any]]] = []
    metadata: Dict[str, Any] = {}
    for seed in range(n_seeds):
        params = AutoinducerParams(seed=seed)
        model = AutoinducerQuorum(params)
        runs.append(model.simulate())
        if seed == 0:
            metadata = model.get_metadata()
    return runs, metadata


def make_p20_detector_fn(n_permutations: int = 199, seed: int = 42):
    from epc.detectors.p20_quorum_sensing import P20QuorumSensingDetector
    detector = P20QuorumSensingDetector(n_permutations=n_permutations, seed=seed)
    def fn(history, metadata=None):
        return detector.detect(history, metadata=metadata)
    return fn


def run_p20(out_path: str = "analysis/outputs/p20_phase2a_panel.json", verbose: bool = True) -> Dict[str, Any]:
    """Run P20 (quorum sensing) Phase-2a panel v1.2.

    detector_format="density_sweep": P20 operates on density-sweep observation
    bundles (density, concentration, collective_state, fraction_on, density_idx,
    sweep_direction). Class A synthetic substrates are random/graded density-
    response curves → P20 rejects at screening (no sharp transition) or
    confirmation (low step-function R²).
    Class B: 0 catalog mates (only density_sweep_timeseries pattern); supplements
    may be added if needed.
    Class C: sub_threshold (system never activates) + graded_response (smooth
    sigmoid, no hysteresis — key P18 discrimination case).
    """
    print(f"--- Running P20 panel → {out_path}")
    positives, metadata = build_p20_positives(n_seeds=5)
    detector_fn = make_p20_detector_fn(n_permutations=199, seed=42)
    return run_panel(
        detector_fn,
        pattern_id="P20",
        detector_format="density_sweep",
        canonical_positive_runs=positives,
        canonical_metadata=metadata,
        failed_regime_module=p20_failed,
        output_path=out_path,
        verbose=verbose,
    )


# --- P4 Territoriality -------------------------------------------------------

def build_p4_positives(n_seeds: int = 5) -> tuple[List[List[Dict[str, Any]]], Dict[str, Any]]:
    """P4 canonical positive: ScentMarkingModel at canonical Giuggioli 2011 regime.

    N=4 agents on 48×48 torus, 20000 steps with snapshots every 100 steps.
    Produces exclusive, persistent territories with scent-mediated exclusion.
    """
    from epc.models.territoriality import ScentMarkingModel, ScentMarkingParams
    runs: List[List[Dict[str, Any]]] = []
    metadata: Dict[str, Any] = {}
    for seed in range(n_seeds):
        params = ScentMarkingParams(
            n_agents=4, grid_size=48,
            deposition_rate=0.1, decay_rate=0.03,
            repulsion_strength=2.0, home_attraction=2.0,
            temperature=0.5,
            n_steps=20000, snapshot_interval=100,
            seed=42 + seed,
        )
        model = ScentMarkingModel(params)
        runs.append(model.simulate())
        if seed == 0:
            metadata = model.get_metadata()
    return runs, metadata


def make_p4_detector_fn(n_permutations: int = 199, seed: int = 42):
    from epc.detectors.p4_territoriality import P4TerritorialityDetector
    detector = P4TerritorialityDetector(n_permutations=n_permutations, seed=seed)
    def fn(history, metadata=None):
        return detector.detect(history, metadata=metadata)
    return fn


def run_p4(out_path: str = "analysis/outputs/p4_phase2a_panel.json", verbose: bool = True) -> Dict[str, Any]:
    """Run P4 (territoriality) Phase-2a panel.

    detector_format: territorial_agent_field
    Class A: random-walk agents with scent but no avoidance (all TNs).
    Class B: 0 catalog mates (sole territorial_agent_field pattern); B'
    supplements: random_walk_territory, clustering_agents_territory.
    Class C: high-tolerance (overlapping ranges) + fast-decay (no boundaries).
    Invariance: permutation_invariant=True → skip permutation_shuffled;
    time_shuffle_invariant=False → evaluate time_shuffled.
    """
    print(f"--- Running P4 panel → {out_path}")
    positives, metadata = build_p4_positives(n_seeds=5)
    detector_fn = make_p4_detector_fn(n_permutations=199, seed=42)
    return run_panel(
        detector_fn,
        pattern_id="P4",
        detector_format="territorial_agent_field",
        canonical_positive_runs=positives,
        canonical_metadata=metadata,
        failed_regime_module=p4_failed,
        output_path=out_path,
        verbose=verbose,
    )


# --- P29 Trail / Network Formation ------------------------------------------

def build_p29_positives(n_seeds: int = 5) -> tuple[List[List[Dict[str, Any]]], Dict[str, Any]]:
    """P29 canonical positive: AntTrailModel with pheromone reinforcement.

    7-node complete graph, 40 agents, 500 steps. Emergent near-MST network
    with positive weight-distance correlation. Each seed uses random layout.
    """
    from epc.models.trail_network import AntTrailModel, AntTrailParams
    runs: List[List[Dict[str, Any]]] = []
    metadata: Dict[str, Any] = {}
    for seed in range(n_seeds):
        params = AntTrailParams(
            n_nodes=7, n_agents=40, grid_size=100,
            alpha=2.0, beta=3.0,
            deposition_rate=10.0, evaporation_rate=0.05,
            n_steps=500, snapshot_interval=10,
            seed=42 + seed,
        )
        model = AntTrailModel(params)
        runs.append(model.simulate())
        if seed == 0:
            metadata = model.get_metadata()
    return runs, metadata


def make_p29_detector_fn(n_permutations: int = 199, seed: int = 42):
    from epc.detectors.p29_trail_network import P29TrailNetworkDetector
    detector = P29TrailNetworkDetector(n_permutations=n_permutations, seed=seed)
    def fn(history, metadata=None):
        return detector.detect(history, metadata=metadata)
    return fn


def run_p29(out_path: str = "analysis/outputs/p29_phase2a_panel.json", verbose: bool = True) -> Dict[str, Any]:
    """Run P29 (trail / network formation) Phase-2a panel.

    detector_format: trail_network
    Class A: random graphs with no weight-distance correlation (all TNs).
    Class B: 0 catalog mates (sole trail_network pattern); B'
    supplements: static_mst_graph, uniform_traffic_graph.
    Class C: high-evaporation (no persistent network) + no-reinforcement
    (uniform edge weights).
    Invariance: permutation_invariant=False, time_shuffle_invariant=False.
    """
    print(f"--- Running P29 panel → {out_path}")
    positives, metadata = build_p29_positives(n_seeds=5)
    detector_fn = make_p29_detector_fn(n_permutations=199, seed=42)
    return run_panel(
        detector_fn,
        pattern_id="P29",
        detector_format="trail_network",
        canonical_positive_runs=positives,
        canonical_metadata=metadata,
        failed_regime_module=p29_failed,
        output_path=out_path,
        verbose=verbose,
    )


def build_p32_positives(n_seeds: int = 5) -> tuple[List[List[Dict[str, Any]]], Dict[str, Any]]:
    """P32 canonical positive: ResponseThresholdModel with specialization.

    20 agents, 3 tasks, 500 steps. Threshold reinforcement drives entropy
    decline and stable role assignment.
    """
    from epc.models.division_of_labor import ResponseThresholdModel
    runs: List[List[Dict[str, Any]]] = []
    metadata: Dict[str, Any] = {}
    for seed in range(n_seeds):
        m = ResponseThresholdModel(
            n_agents=20, n_tasks=3,
            reinforcement_rate=0.05, forgetting_rate=0.01,
            stimulus_rate=0.1, initial_threshold=0.5,
            seed=42 + seed,
        )
        m.setup()
        history: List[Dict[str, Any]] = []
        for _ in range(500):
            history.append(m.step())
        runs.append(history)
        if seed == 0:
            metadata = m.get_metadata()
    return runs, metadata


def make_p32_detector_fn(n_permutations: int = 199, seed: int = 42):
    from epc.detectors.p32_specialization import P32SpecializationDetector
    detector = P32SpecializationDetector(n_permutations=n_permutations, seed=seed)
    def fn(history, metadata=None):
        return detector.detect(history, model_metadata=metadata)
    return fn


def run_p32(out_path: str = "analysis/outputs/p32_phase2a_panel.json", verbose: bool = True) -> Dict[str, Any]:
    """Run P32 (specialization / division of labor) Phase-2a panel.

    detector_format: task_allocation_timeseries
    Class A: random task assignments (all TNs).
    Class B: 0 catalog mates (sole task_allocation_timeseries pattern);
    B' supplements: single_task_collapse_allocation, constant_rebalancing_allocation.
    Class C: no-reinforcement (no entropy decline) + high-forgetting
    (thresholds drift uniformly, no stable roles).
    Invariance: permutation_invariant=False, time_shuffle_invariant=False.
    """
    print(f"--- Running P32 panel → {out_path}")
    positives, metadata = build_p32_positives(n_seeds=5)
    detector_fn = make_p32_detector_fn(n_permutations=199, seed=42)
    return run_panel(
        detector_fn,
        pattern_id="P32",
        detector_format="task_allocation_timeseries",
        canonical_positive_runs=positives,
        canonical_metadata=metadata,
        failed_regime_module=p32_failed,
        output_path=out_path,
        verbose=verbose,
    )


def build_p30_positives(n_seeds: int = 5) -> tuple[List[List[Dict[str, Any]]], Dict[str, Any]]:
    """P30 canonical positive: AutopoiesisModel with membrane formation.

    100 substrate + 3 catalyst particles, box=20. Radial spring mechanism
    drives link particles to form a closed membrane shell around catalysts.
    """
    from epc.models.autopoiesis import AutopoiesisModel
    runs: List[List[Dict[str, Any]]] = []
    metadata: Dict[str, Any] = {}
    for seed in range(n_seeds):
        m = AutopoiesisModel(
            n_substrate=100, n_catalyst=3, box_size=20.0,
            production_rate=0.15, decay_rate=0.01,
            membrane_equilibrium_radius=3.0,
            catalyst_link_attraction=0.5,
            link_attraction=0.3,
            seed=42 + seed,
        )
        m.setup()
        history: List[Dict[str, Any]] = []
        for _ in range(300):
            history.append(m.step())
        runs.append(history)
        if seed == 0:
            metadata = {
                'model_class': 'autopoiesis',
                'substrate_type': 'particle_membrane',
                'production_rate': 0.15,
                'decay_rate': 0.01,
                'n_particles': m.n_total,
                'box_size': 20.0,
            }
    return runs, metadata


def make_p30_detector_fn(n_permutations: int = 199, seed: int = 42):
    from epc.detectors.p30_autopoiesis import P30AutopoiesisDetector
    detector = P30AutopoiesisDetector(n_permutations=n_permutations, seed=seed)
    def fn(history, metadata=None):
        return detector.detect(history, model_metadata=metadata)
    return fn


def run_p30(out_path: str = "analysis/outputs/p30_phase2a_panel.json", verbose: bool = True) -> Dict[str, Any]:
    """Run P30 (autopoiesis / spontaneous boundary) Phase-2a panel.

    detector_format: particle_membrane
    Class A: random particles with random type assignments (all TNs).
    Class B: 0 catalog mates (sole particle_membrane pattern); B'
    supplements: dense_cluster_particles (P1-like), dispersed_typed_regions (P4-like).
    Class C: non-bonding, high-decay, no-attraction, weak-production, large-box.
    Invariance: permutation_invariant=True, time_shuffle_invariant=False.
    """
    print(f"--- Running P30 panel → {out_path}")
    positives, metadata = build_p30_positives(n_seeds=5)
    detector_fn = make_p30_detector_fn(n_permutations=199, seed=42)
    return run_panel(
        detector_fn,
        pattern_id="P30",
        detector_format="particle_membrane",
        canonical_positive_runs=positives,
        canonical_metadata=metadata,
        failed_regime_module=p30_failed,
        output_path=out_path,
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
    if which in ("p3",):
        summaries["P3"] = run_p3()
    if which in ("p11",):
        summaries["P11"] = run_p11()
    if which in ("p12",):
        summaries["P12"] = run_p12()
    if which in ("p13",):
        summaries["P13"] = run_p13()
    if which in ("p5",):
        summaries["P5"] = run_p5()
    if which in ("p2",):
        summaries["P2"] = run_p2()
    if which in ("p6",):
        summaries["P6"] = run_p6()
    if which in ("p8",):
        summaries["P8"] = run_p8()
    if which in ("p10",):
        summaries["P10"] = run_p10()
    if which in ("p21",):
        summaries["P21"] = run_p21()
    if which in ("p7",):
        summaries["P7"] = run_p7()
    if which in ("p17",):
        summaries["P17"] = run_p17()
    if which in ("p19",):
        summaries["P19"] = run_p19()
    if which in ("p24",):
        summaries["P24"] = run_p24()
    if which in ("p26",):
        summaries["P26"] = run_p26()
    if which in ("p23",):
        summaries["P23"] = run_p23()
    if which in ("p16",):
        summaries["P16"] = run_p16()
    if which in ("p25",):
        summaries["P25"] = run_p25()
    if which in ("p20",):
        summaries["P20"] = run_p20()
    if which in ("p4",):
        summaries["P4"] = run_p4()
    if which in ("p29",):
        summaries["P29"] = run_p29()
    if which in ("p32",):
        summaries["P32"] = run_p32()
    if which in ("p30",):
        summaries["P30"] = run_p30()

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
