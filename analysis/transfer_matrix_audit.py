"""Full transfer-matrix tier audit (the real 40 x 32 audit, no carry-forward).

For every substrate-and-observable-compatible (model, detector) cell, the model
is run with FIVE seeds and the detector is run on each; the cell's reported tier
is the MEDIAN over the five seeds (robust to a single lucky/unlucky seed in
either direction — the seed-fragility of P2/P27/P29 documented in §4A.4 would
otherwise distort a single-seed read). The achieved tier is one of
definitive / confirmation / screening / none.

Canonical builders are the ALREADY-VALIDATED ones: the Phase-2a panel positives
(analysis.run_phase2a_panel.build_pXX_positives, five seeds) for the 31 pattern
primaries, the T2c cross-model builders (analysis.t2c_systems.recog_*) for the
cross-model alternatives, and canonical constructors for the remaining variants.

Special cases recorded honestly, not hidden:
  - P31 has no single-run panel detector (validated by the separate multi-run
    non-redundancy test) -> 'n/a_nonredundancy'.
  - A detector that raises on input it cannot parse (e.g. P1 on Zhang's
    non-integer 'cell_types') -> 'input_rejected' (the substrate-content
    prerequisite is not met; the honest verdict is a rejection).

Emits analysis/outputs/transfer_matrix_audit.json (per-cell median tier + the
five per-seed tiers) and prints a summary.
"""
from __future__ import annotations

import json
import time
from collections import Counter
from typing import Any, Dict, List, Optional, Tuple

import analysis.run_phase2a_panel as R
from analysis import t2c_systems as T
from analysis.battery_profile import build_detector_fns
from epc.orchestration import DETECTOR_REGISTRY, MODEL_REGISTRY, get_compatible_pairs
from epc.phase2a.panel import _detected, _verdict

N_SEEDS = 5
TIER_RANK = {"input_rejected": 0, "none": 1, "screening": 2, "confirmation": 3, "definitive": 4}
RANK_TIER = {v: k for k, v in TIER_RANK.items()}
History = Any
SeedRuns = List[Tuple[History, Optional[dict]]]


def _md(m) -> Optional[Dict[str, Any]]:
    try:
        md = m.get_metadata()
        return md if isinstance(md, dict) else None
    except Exception:
        return None


def _bp(fn_name: str):
    """Five Phase-2a panel positives (one per seed)."""
    def f() -> SeedRuns:
        runs, meta = getattr(R, fn_name)(n_seeds=N_SEEDS)
        m = meta if isinstance(meta, dict) else None
        return [(r, m) for r in runs]
    return f


def _recog(fn_name: str):
    def f() -> SeedRuns:
        out: SeedRuns = []
        for s in range(N_SEEDS):
            h, meta = getattr(T, fn_name)(s)
            out.append((h, meta if isinstance(meta, dict) else None))
        return out
    return f


def _zhang() -> SeedRuns:
    from epc.models.cell_view_sorting import CellViewSorting
    out: SeedRuns = []
    for s in range(N_SEEDS):
        m = CellViewSorting(n=60, algorithm="insertion", seed=s)
        out.append((m.run_to_completion(max_rounds=400), _md(m)))
    return out


def _physarum() -> SeedRuns:
    from epc.models.trail_network import PhysarumModel, PhysarumParams
    out: SeedRuns = []
    for s in range(N_SEEDS):
        m = PhysarumModel(PhysarumParams(seed=s))
        out.append((m.simulate(), _md(m)))
    return out


def _pheromone() -> SeedRuns:
    from epc.models.territoriality import PheromoneRepulsionModel, PheromoneRepulsionParams
    out: SeedRuns = []
    for s in range(N_SEEDS):
        m = PheromoneRepulsionModel(PheromoneRepulsionParams(seed=s))
        out.append((m.simulate(), _md(m)))
    return out


def _no_reinforcement() -> SeedRuns:
    from epc.models.division_of_labor import NoReinforcementModel
    out: SeedRuns = []
    for s in range(N_SEEDS):
        try:
            m = NoReinforcementModel(seed=s)
        except TypeError:
            m = NoReinforcementModel()
        m.setup()
        h = [m.step() for _ in range(1000)]
        out.append((h, _md(m)))
    return out


MODEL_BUILDERS = {
    "schelling": _bp("build_p1_positives"), "abp": _bp("build_p2_positives"),
    "gray_scott": _bp("build_p3_positives"), "scent_marking_territory": _bp("build_p4_positives"),
    "vicsek": _bp("build_p5_positives"), "dorsogna": _bp("build_p6_positives"),
    "lane_formation": _bp("build_p7_positives"), "nagel_schreckenberg": _bp("build_p8_positives"),
    "kuramoto": _bp("build_p9_positives"), "kuramoto_nonlocal": _bp("build_p10_positives"),
    "lotka_volterra": _bp("build_p11_positives"), "rps_spatial": _bp("build_p12_positives"),
    "greenberg_hastings": _bp("build_p13_positives"), "btw_sandpile": _bp("build_p14_positives"),
    "game_of_life": _bp("build_p15_positives"), "hopfield": _bp("build_p16_positives"),
    "collective_sensing": _bp("build_p17_positives"), "voter": _bp("build_p18_positives"),
    "informed_minority": _bp("build_p19_positives"), "autoinducer_quorum": _bp("build_p20_positives"),
    "hegselmann_krause": _bp("build_p21_positives"), "sir_epidemic": _bp("build_p22_positives"),
    "minority_game": _bp("build_p23_positives"), "proportional_homeostat": _bp("build_p24_positives"),
    "canalized_landscape": _bp("build_p25_positives"), "bistable_double_well": _bp("build_p26_positives"),
    "nowak_may": _bp("build_p27_positives"), "yard_sale": _bp("build_p28_positives"),
    "ant_trail_network": _bp("build_p29_positives"), "autopoiesis": _bp("build_p30_positives"),
    "response_threshold": _bp("build_p32_positives"),
    "boolean_grn": _recog("recog_p16_boolean_grn"), "fraction_threshold": _recog("recog_p20_fraction_threshold"),
    "el_farol": _recog("recog_p23_el_farol"), "multi_basin_grn": _recog("recog_p25_multibasin"),
    "zhang_sequential": _zhang, "zhang_threaded": _zhang, "physarum_network": _physarum,
    "pheromone_repulsion_territory": _pheromone, "no_reinforcement": _no_reinforcement,
}

PRIMARY = {name: list(getattr(MODEL_REGISTRY[name], "primary_patterns", []) or [])
           for name in MODEL_REGISTRY}


def _median_tier(tiers: List[str]) -> str:
    ranks = sorted(TIER_RANK.get(t, 0) for t in tiers)
    return RANK_TIER[ranks[len(ranks) // 2]]


def run():
    t0 = time.time()
    names = sorted(MODEL_REGISTRY.keys())
    assert not (set(names) - set(MODEL_BUILDERS)), f"missing: {set(names)-set(MODEL_BUILDERS)}"
    fns = build_detector_fns()
    pairs = get_compatible_pairs()
    print(f"{len(names)} models x {len(DETECTOR_REGISTRY)} detectors, {len(pairs)} compatible cells, "
          f"{N_SEEDS} seeds each", flush=True)

    hist: Dict[str, Any] = {}
    for n in names:
        tb = time.time()
        try:
            hist[n] = MODEL_BUILDERS[n]()
            print(f"  built {n:<30} {len(hist[n])} seeds ({time.time()-tb:.1f}s)", flush=True)
        except Exception as e:
            hist[n] = ("BUILD_ERROR", repr(e))
            print(f"  BUILD_ERROR {n}: {e}", flush=True)

    cells = []
    for model, det in pairs:
        h = hist[model]
        on_diag = det in PRIMARY.get(model, [])
        if isinstance(h, tuple) and h and h[0] == "BUILD_ERROR":
            cells.append({"model": model, "detector": det, "tier": "build_error",
                          "detected": False, "on_diagonal": on_diag, "err": h[1][:120]})
            continue
        if det not in fns:
            tier = "n/a_nonredundancy" if det == "P31" else "no_detector_fn"
            cells.append({"model": model, "detector": det, "tier": tier,
                          "detected": False, "on_diagonal": on_diag})
            continue
        seed_tiers, err = [], None
        for history, meta in h:
            try:
                res = fns[det](history, meta)
                seed_tiers.append(_verdict(res) if _detected(res) else "none")
            except Exception as e:
                seed_tiers.append("input_rejected")
                err = repr(e)[:120]
        tier = _median_tier(seed_tiers)
        rec = {"model": model, "detector": det, "tier": tier, "seed_tiers": seed_tiers,
               "detected": tier in ("screening", "confirmation", "definitive"), "on_diagonal": on_diag}
        if err:
            rec["err"] = err
        cells.append(rec)

    out = {"n_models": len(names), "n_detectors": len(DETECTOR_REGISTRY),
           "n_compatible_cells": len(pairs), "n_seeds": N_SEEDS, "tier_rule": "median",
           "cells": cells}
    json.dump(out, open("analysis/outputs/transfer_matrix_audit.json", "w"), indent=1)

    tiers = Counter(c["tier"] for c in cells)
    diag = [c for c in cells if c.get("on_diagonal")]
    diag_def = [c for c in diag if c["tier"] == "definitive"]
    off_fire = [c for c in cells if c.get("detected") and not c.get("on_diagonal")]
    print(f"\n=== AUDIT COMPLETE ({time.time()-t0:.0f}s) ===", flush=True)
    print("median-tier histogram over all compatible cells:", dict(tiers), flush=True)
    print(f"on-diagonal cells: {len(diag)}; DEFINITIVE: {len(diag_def)}; "
          f"fired(any tier): {sum(1 for c in diag if c['detected'])}", flush=True)
    for c in sorted(diag, key=lambda c: c["detector"]):
        tag = "" if c["detected"] else "  <-- not fired"
        print(f"   {c['model']:<30} x {c['detector']:<4} -> {c['tier']:<14} {c.get('seed_tiers','')}{tag}", flush=True)
    print(f"\noff-diagonal firings (co-occurrences): {len(off_fire)}", flush=True)
    for c in off_fire:
        print(f"   {c['model']} x {c['detector']} -> {c['tier']} {c.get('seed_tiers','')}", flush=True)
    special = [c for c in cells if c["tier"] in ("build_error", "n/a_nonredundancy", "no_detector_fn")
               or "input_rejected" in c.get("seed_tiers", [])]
    print(f"\nnon-tiered / special cells: {len(special)}", flush=True)
    for c in special:
        print(f"   {c['model']} x {c['detector']} -> {c['tier']}: {c.get('err','')}", flush=True)
    print("wrote analysis/outputs/transfer_matrix_audit.json", flush=True)


if __name__ == "__main__":
    run()
