"""P27 Nowak-May spatial prisoner's dilemma: literature reproduction v2 (group_A.md, P27).

Target: Szabo & Fath 2007 Fig. 26 reading, cooperator frequency 0.316 (+/- 0.005 reading)
for 9/5 < b < 2, Moore + self, random start rho_C(0) = 0.5.

Gating runs: NowakMayModel(L, L, b=1.9, init_mode='random', init_coop_fraction=0.5,
boundary='periodic', seed=s), L = 1000 (gating), 500 and 200 (finite size), seeds 0-9,
1500 generations; burn-in generations 1-500, measure f_C over 501-1500.
Pass (all three): |mean f_C(1000) - 0.316| <= 0.010; |mean f_C(500) - mean f_C(1000)| <= 0.005;
at most 1 of 10 seeds absorbing at L = 1000. The break-point invariance check (L = 200,
seed 0, 300 generations, b = 1.85 and 1.95 bit-identical to b = 1.9 at every generation)
must pass before the gate is evaluated.
Non-gating: b = 1.8 at L = 1000, seeds 0-9; comparison of the gating mean with 12 ln 2 - 8.

Subcommands
  invariance                 break-point invariance check (uses NowakMayModel.step itself)
  run --L L --b B --seeds A-B
  analyze                    writes analysis/outputs/repro_v2/P27.json
"""
from __future__ import annotations

import argparse
import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common  # noqa: E402
import nm_fast  # noqa: E402

from epc.models.nowak_may import NowakMayModel  # noqa: E402

TARGET = 0.316
TOL_1 = 0.010
TOL_2 = 0.005
MAX_ABSORBING = 1
ANALYTIC_NM = 12 * np.log(2) - 8          # 0.31777, attributed to Nowak & May 1992 (non-gating)
GENS = 1500
BURN = 500
L_GATE, L_FS, L_SMALL = 1000, 500, 200
B_GATE, B_BREAK = 1.9, 1.8
INV_L, INV_SEED, INV_GENS = 200, 0, 300
RAW = common.RAW_DIR
INV_FILE = RAW / "P27_breakpoint_invariance.json"


def raw_path(L: int, b: float) -> Path:
    return RAW / f"P27_L{L}_b{b:g}.jsonl"


def make(L: int, b: float, seed: int) -> NowakMayModel:
    return NowakMayModel(rows=L, cols=L, b=b, init_mode="random", init_coop_fraction=0.5,
                         boundary="periodic", seed=seed)


def cmd_invariance() -> None:
    t0 = time.time()
    L, seed, gens = INV_L, INV_SEED, INV_GENS
    models = {b: make(L, b, seed) for b in (1.9, 1.85, 1.95)}
    for m in models.values():
        m.setup()
    first_diff = {1.85: None, 1.95: None}
    same0 = all(np.array_equal(models[1.9].grid, models[b].grid) for b in (1.85, 1.95))
    for g in range(1, gens + 1):
        for m in models.values():
            m.step()
        for b in (1.85, 1.95):
            if first_diff[b] is None and not np.array_equal(models[1.9].grid, models[b].grid):
                first_diff[b] = g
    ok = same0 and all(v is None for v in first_diff.values())
    res = {"check": "break-point invariance (implementation check, pre-registered)",
           "L": L, "seed": seed, "generations": gens, "b_reference": 1.9, "b_compared": [1.85, 1.95],
           "model_step_used": "epc.models.nowak_may.NowakMayModel.step",
           "initial_grids_identical": bool(same0),
           "first_generation_with_difference": {str(k): v for k, v in first_diff.items()},
           "pass": bool(ok), "runtime_s": time.time() - t0, "git": common.git_info(),
           "timestamp": common.now_iso()}
    common.write_json(INV_FILE, res)
    print("invariance", "PASS" if ok else "FAIL", first_diff)
    if not ok:
        raise SystemExit("break-point invariance check failed: implementation bug, stop before gating runs")


def run_one(L: int, b: float, seed: int) -> dict:
    t0 = time.time()
    model = make(L, b, seed)
    model.setup()
    n_c = [int(np.count_nonzero(model.grid == 0))]
    for _ in range(GENS):
        n_c.append(nm_fast.generation(model))
    n_c_arr = np.array(n_c)
    fc = n_c_arr / (L * L)
    absorbed = np.nonzero((n_c_arr == 0) | (n_c_arr == L * L))[0]
    win = fc[BURN + 1:GENS + 1]                     # generations 501..1500
    return {"L": L, "b": b, "seed": seed, "generations": GENS, "burn_in_generations": BURN,
            "fC_time_mean_501_1500": float(win.mean()),
            "fC_temporal_sd_501_1500": float(win.std(ddof=1)),
            "fC_final": float(fc[-1]), "fC_initial": float(fc[0]),
            "absorbing": bool(absorbed.size > 0),
            "first_absorbing_generation": int(absorbed[0]) if absorbed.size else None,
            "n_coop_series": n_c,                    # n_C after generations 0..1500
            "runtime_s": time.time() - t0}


def cmd_run(L: int, b: float, seeds: list[int]) -> None:
    git = common.git_info()
    out = raw_path(L, b)
    done = {r["seed"] for r in common.read_jsonl(out)}
    for s in seeds:
        if s in done:
            continue
        rec = run_one(L, b, s)
        rec.update(git_commit=git["commit"], tracked_tree_clean=git["tracked_tree_clean"],
                   timestamp=common.now_iso())
        common.append_jsonl(out, rec)
        print(f"P27 L={L} b={b} seed={s} {rec['runtime_s']:.1f}s", flush=True)


def summarize(recs: list[dict]) -> dict:
    recs = sorted(recs, key=lambda r: r["seed"])
    keep = [r for r in recs if not r["absorbing"]]
    m = np.array([r["fC_time_mean_501_1500"] for r in keep])
    return {"n_seeds": len(recs), "seeds": [r["seed"] for r in recs],
            "n_absorbing": len(recs) - len(keep),
            "absorbing_seeds": [r["seed"] for r in recs if r["absorbing"]],
            "mean_fC": float(m.mean()) if len(m) else None,
            "se_fC": float(m.std(ddof=1) / np.sqrt(len(m))) if len(m) > 1 else None,
            "sd_across_seeds": float(m.std(ddof=1)) if len(m) > 1 else None,
            "per_seed": [{"seed": r["seed"], "fC_time_mean_501_1500": r["fC_time_mean_501_1500"],
                          "fC_temporal_sd_501_1500": r["fC_temporal_sd_501_1500"], "fC_final": r["fC_final"],
                          "absorbing": r["absorbing"], "first_absorbing_generation": r["first_absorbing_generation"],
                          "runtime_s": r["runtime_s"]} for r in recs],
            "runtime_s": float(sum(r["runtime_s"] for r in recs))}


def cmd_analyze() -> None:
    inv = json.loads(INV_FILE.read_text())
    sets = {}
    allrecs = []
    for (L, b) in [(L_GATE, B_GATE), (L_FS, B_GATE), (L_SMALL, B_GATE), (L_GATE, B_BREAK)]:
        recs = common.read_jsonl(raw_path(L, b))
        allrecs += recs
        if recs:
            sets[f"L{L}_b{b:g}"] = {"L": L, "b": b, **summarize(recs)}
    out = {"pattern": "P27", "title": "Nowak-May spatial prisoner's dilemma, spatial reciprocity",
           "spec": "docs/validation_rebuild/repro_v2/group_A.md#p27",
           "target": {"value": TARGET, "source": "Szabo & Fath 2007 Fig. 26 (reading, +/- 0.005)",
                      "analytic_attributed_to_Nowak_May_1992": float(ANALYTIC_NM)},
           "settings": {"model": "epc.models.nowak_may.NowakMayModel", "init_mode": "random",
                        "init_coop_fraction": 0.5, "boundary": "periodic", "generations": GENS,
                        "burn_in": "generations 1-500", "measurement": "generations 501-1500",
                        "estimator": "per-seed time mean of f_C; ensemble = mean over non-absorbing seeds, SE = SD/sqrt(n)",
                        "absorbing": "f_C reaches 0 or 1 at any generation; counted and excluded from the mean"},
           "breakpoint_invariance_check": inv, "runs": sets}
    gate_ok = bool(inv["pass"])
    s1000, s500 = sets.get(f"L{L_GATE}_b{B_GATE:g}"), sets.get(f"L{L_FS}_b{B_GATE:g}")
    rules = {}
    if s1000 and s500:
        d1 = abs(s1000["mean_fC"] - TARGET)
        d2 = abs(s500["mean_fC"] - s1000["mean_fC"])
        rules = {
            "rule1": {"text": "|mean f_C(L=1000) - 0.316| <= 0.010", "value": s1000["mean_fC"],
                      "abs_diff": d1, "pass": bool(d1 <= TOL_1)},
            "rule2": {"text": "|mean f_C(L=500) - mean f_C(L=1000)| <= 0.005", "L500": s500["mean_fC"],
                      "L1000": s1000["mean_fC"], "abs_diff": d2, "pass": bool(d2 <= TOL_2)},
            "rule3": {"text": "at most 1 of the 10 seeds at L=1000 absorbing", "n_absorbing": s1000["n_absorbing"],
                      "n_seeds": s1000["n_seeds"], "pass": bool(s1000["n_absorbing"] <= MAX_ABSORBING
                                                                and s1000["n_seeds"] == 10)},
        }
    out["pass_rules"] = rules
    out["gate_evaluable"] = gate_ok
    out["primary_pass"] = bool(gate_ok and rules and all(r["pass"] for r in rules.values()))
    out["reported_only"] = {
        f"L{L_SMALL}_b{B_GATE:g}": sets.get(f"L{L_SMALL}_b{B_GATE:g}"),
        f"b{B_BREAK:g}_L{L_GATE}_breakpoint": sets.get(f"L{L_GATE}_b{B_BREAK:g}"),
        "gating_mean_minus_12ln2_minus_8": (s1000["mean_fC"] - float(ANALYTIC_NM)) if s1000 else None,
    }
    out["deviations"] = common.load_deviations("P27")
    out["seeds"] = "0-9 for every L and b"
    out["code_commits_of_runs"] = sorted({r["git_commit"] for r in allrecs})
    out["all_runs_from_clean_tree"] = all(r["tracked_tree_clean"] for r in allrecs)
    out["analysis_git"] = common.git_info()
    out["bit_identity_checks"] = "analysis/outputs/repro_v2/bit_identity_checks.json"
    out["total_compute_runtime_s"] = float(sum(r["runtime_s"] for r in allrecs) + inv["runtime_s"])
    out["fC_series_file"] = "analysis/outputs/repro_v2/raw/P27_L*_b*.jsonl (n_coop_series, generations 0-1500)"
    out["timestamp"] = common.now_iso()
    common.write_json(common.OUT_DIR / "P27.json", out)
    print("P27 primary", "PASS" if out["primary_pass"] else "FAIL")


def main() -> None:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    sub.add_parser("invariance")
    r = sub.add_parser("run"); r.add_argument("--L", type=int, required=True)
    r.add_argument("--b", type=float, required=True); r.add_argument("--seeds", required=True)
    sub.add_parser("analyze")
    a = ap.parse_args()
    if a.cmd == "invariance":
        cmd_invariance()
    elif a.cmd == "run":
        lo, hi = a.seeds.split("-")
        cmd_run(a.L, a.b, list(range(int(lo), int(hi) + 1)))
    else:
        cmd_analyze()


if __name__ == "__main__":
    main()
