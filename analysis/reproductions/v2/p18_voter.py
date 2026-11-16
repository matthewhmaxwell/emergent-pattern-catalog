"""P18 2-D voter model: literature reproduction v2 (group_A.md, section P18).

Targets: exact infinite-lattice interface density rho(t) of the n.n. voter model
(Frachebourg & Krapivsky 1996 Eqs. 9-11, evaluated by group_A_support/voter_exact_crw.py)
and the asymptotic slope d(1/rho)/d ln t = 2/pi (Dornic et al. 2001 Eq. 1).

Gating runs: VoterModel(512, 512, init_mode='random', neighborhood='von_neumann',
boundary='periodic', seed=s), seeds 0-19 (extended to 0-39 if the precision requirement
SE(3000) <= 0.01 rho_exact(3000) is not met), 3000 sweeps; wall_density recorded at
t = 10, 30, 100, 300, 500, 1000, 2000, 3000.
Finite-size check: von Neumann L = 256, seeds 0-9. Secondary: Moore L = 256, seeds 0-9.
T2: |rho_bar(t) - rho_exact(t)| <= 0.03 rho_exact(t) at t = 100, 300, 1000, 3000 (L = 512).
T1: OLS slope of 1/rho_bar(t) on ln t over t = 300, 500, 1000, 2000, 3000 in [0.586, 0.688].
Finite size: rho_bar(3000) at L = 256 and 512 differ by <= 2% (relative to L = 512).

Subcommands
  run --L L --nb {von_neumann,moore} --seeds A-B
  precision           checks SE(3000) at L=512 vN; runs seeds 20-39 if the requirement fails
  theory              runs voter_exact_crw.py and caches its values
  analyze             writes analysis/outputs/repro_v2/P18.json
"""
from __future__ import annotations

import argparse
import contextlib
import io
import json
import runpy
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))
import common  # noqa: E402
import voter_fast  # noqa: E402

from epc.models import voter as voter_mod  # noqa: E402
from epc.models.voter import VoterModel  # noqa: E402

REC_T = [10, 30, 100, 300, 500, 1000, 2000, 3000]
T_GATE = [100, 300, 1000, 3000]
T_FIT = [300, 500, 1000, 2000, 3000]
T_MAX = 3000
SLOPE_WINDOW = (0.586, 0.688)
TOL_T2 = 0.03
TOL_FS = 0.02
PREC = 0.01
N_BOOT = 2000
BOOT_SEED = 20261117
L_GATE, L_FS = 512, 256
N_GATE_SEEDS, N_MAX_SEEDS = 20, 40
RAW = common.RAW_DIR
THEORY_FILE = RAW / "P18_theory_voter_exact_crw.json"


def raw_path(L: int, nb: str) -> Path:
    return RAW / f"P18_L{L}_{nb}.jsonl"


def run_one(L: int, nb: str, seed: int) -> dict:
    t0 = time.time()
    model = VoterModel(rows=L, cols=L, init_mode="random", neighborhood=nb, boundary="periodic", seed=seed)
    model.setup()
    rec = {"0": voter_mod._wall_density(model.grid, model._n_offsets)}
    mag = {"0": float(2.0 * np.mean(model.grid == 1) - 1.0)}
    want = set(REC_T)
    for t in range(1, T_MAX + 1):
        voter_fast.sweep(model)
        if t in want:
            rec[str(t)] = voter_mod._wall_density(model.grid, model._n_offsets)
            mag[str(t)] = float(2.0 * np.mean(model.grid == 1) - 1.0)
    return {"L": L, "neighborhood": nb, "seed": seed, "sweeps": T_MAX, "wall_density": rec,
            "magnetization": mag, "runtime_s": time.time() - t0}


def cmd_run(L: int, nb: str, seeds: list[int]) -> None:
    git = common.git_info()
    out = raw_path(L, nb)
    done = {r["seed"] for r in common.read_jsonl(out)}
    for s in seeds:
        if s in done:
            continue
        rec = run_one(L, nb, s)
        rec.update(git_commit=git["commit"], tracked_tree_clean=git["tracked_tree_clean"],
                   timestamp=common.now_iso())
        common.append_jsonl(out, rec)
        print(f"P18 L={L} {nb} seed={s} {rec['runtime_s']:.1f}s", flush=True)


def cmd_theory() -> dict:
    if THEORY_FILE.exists():
        return json.loads(THEORY_FILE.read_text())
    t0 = time.time()
    buf = io.StringIO()
    with contextlib.redirect_stdout(buf):
        ns = runpy.run_path(str(common.SUPPORT_DIR / "voter_exact_crw.py"))
    res = ns["res"]
    th = {"script": "docs/validation_rebuild/repro_v2/group_A_support/voter_exact_crw.py",
          "stdout": buf.getvalue(), "runtime_s": time.time() - t0,
          "rho_exact_von_neumann": {str(t): float(res[("vonNeumann", t)]) for t in ns["ts"]},
          "rho_exact_moore": {str(t): float(res[("Moore", t)]) for t in ns["ts"]}}
    for name, key in (("von_neumann", "rho_exact_von_neumann"), ("moore", "rho_exact_moore")):
        th[f"exact_slope_{name}_300_3000"] = float(np.polyfit(np.log(T_FIT), [1 / th[key][str(t)] for t in T_FIT], 1)[0])
    g = 0.5772156649
    th["FK_Eq23_reported_only"] = {str(t): float(np.pi / (2 * np.log(t) + np.log(256))) for t in ns["ts"]}
    th["FK_Eq23_plus_2gamma_reported_only"] = {str(t): float(np.pi / (2 * np.log(t) + np.log(256) + 2 * g))
                                               for t in ns["ts"]}
    common.write_json(THEORY_FILE, th)
    print(buf.getvalue())
    return th


def _matrix(recs: list[dict], ts: list[int]) -> np.ndarray:
    recs = sorted(recs, key=lambda r: r["seed"])
    return np.array([[r["wall_density"][str(t)] for t in ts] for r in recs])   # (n_seeds, n_t)


def _slope(rho_bar: np.ndarray) -> float:
    if np.any(rho_bar <= 0):          # consensus in every seed (only possible on tiny smoke lattices)
        return float("nan")
    return float(np.polyfit(np.log(T_FIT), 1.0 / rho_bar, 1)[0])


def cmd_precision() -> None:
    th = cmd_theory()
    recs = [r for r in common.read_jsonl(raw_path(L_GATE, "von_neumann")) if r["seed"] < N_GATE_SEEDS]
    if len(recs) < N_GATE_SEEDS:
        raise SystemExit(f"only {len(recs)} of {N_GATE_SEEDS} gating seeds present")
    v = _matrix(recs, [3000])[:, 0]
    se = float(v.std(ddof=1) / np.sqrt(len(v)))
    req = PREC * th["rho_exact_von_neumann"]["3000"]
    res = {"rule": "SE(t=3000) <= 0.01 rho_exact(3000) with seeds 0-19; otherwise add seeds up to 40",
           "n_seeds": len(v), "se_3000": se, "requirement": req, "met": bool(se <= req),
           "action": "none" if se <= req else "ran seeds 20-39; all rules judged on 40 seeds",
           "timestamp": common.now_iso()}
    common.write_json(RAW / "P18_precision_check.json", res)
    print(res)
    if not res["met"]:
        cmd_run(L_GATE, "von_neumann", list(range(N_GATE_SEEDS, N_MAX_SEEDS)))


def _block(recs: list[dict], exact: dict, ts: list[int]) -> dict:
    M = _matrix(recs, ts)
    bar = M.mean(0)
    se = M.std(0, ddof=1) / np.sqrt(M.shape[0])
    rows = {}
    for j, t in enumerate(ts):
        e = exact.get(str(t))
        rows[str(t)] = {"rho_bar": float(bar[j]), "se": float(se[j]), "sd": float(M[:, j].std(ddof=1)),
                        "rho_exact": e, "rel_diff": (float(bar[j] / e - 1) if e else None),
                        "per_seed": M[:, j].tolist()}
    return {"n_seeds": M.shape[0], "seeds": sorted(r["seed"] for r in recs), "by_t": rows,
            "runtime_s": float(sum(r["runtime_s"] for r in recs))}


def _slope_with_ci(recs: list[dict]) -> dict:
    M = _matrix(recs, T_FIT)
    s = _slope(M.mean(0))
    rng = np.random.default_rng(BOOT_SEED)
    n = M.shape[0]
    bs = np.array([_slope(M[rng.integers(0, n, n)].mean(0)) for _ in range(N_BOOT)])
    bs = bs[np.isfinite(bs)] if np.any(np.isfinite(bs)) else np.array([np.nan])
    return {"slope": s, "fit_t": T_FIT, "fit": "OLS of 1/rho_bar(t) on ln t",
            "bootstrap": {"n_resamples": N_BOOT, "rng_seed": BOOT_SEED, "resampling": "seeds with replacement",
                          "ci95": [float(np.percentile(bs, 2.5)), float(np.percentile(bs, 97.5))],
                          "sd": float(bs.std(ddof=1))}}


def cmd_analyze() -> None:
    th = cmd_theory()
    ex_vn, ex_mo = th["rho_exact_von_neumann"], th["rho_exact_moore"]
    prec = json.loads((RAW / "P18_precision_check.json").read_text())
    vn512 = common.read_jsonl(raw_path(L_GATE, "von_neumann"))
    n_gate = N_GATE_SEEDS if prec["met"] else N_MAX_SEEDS
    vn512 = [r for r in vn512 if r["seed"] < n_gate]
    if len(vn512) != n_gate:
        raise SystemExit(f"expected {n_gate} gating seeds, found {len(vn512)}")
    vn256 = common.read_jsonl(raw_path(L_FS, "von_neumann"))
    mo256 = common.read_jsonl(raw_path(L_FS, "moore"))
    g = _block(vn512, ex_vn, [0] + REC_T)
    t2_rows = {str(t): {"rho_bar": g["by_t"][str(t)]["rho_bar"], "rho_exact": ex_vn[str(t)],
                        "rel_diff": g["by_t"][str(t)]["rel_diff"],
                        "pass": bool(abs(g["by_t"][str(t)]["rho_bar"] - ex_vn[str(t)]) <= TOL_T2 * ex_vn[str(t)])}
               for t in T_GATE}
    t2_pass = all(r["pass"] for r in t2_rows.values())
    sl = _slope_with_ci(vn512)
    t1_pass = SLOPE_WINDOW[0] <= sl["slope"] <= SLOPE_WINDOW[1]
    out = {"pattern": "P18", "title": "2-D voter model, consensus and coarsening",
           "spec": "docs/validation_rebuild/repro_v2/group_A.md#p18", "theory": th,
           "settings": {"model": "epc.models.voter.VoterModel", "init_mode": "random", "boundary": "periodic",
                        "gating": "L=512, neighborhood='von_neumann'", "sweeps": T_MAX,
                        "observable": "voter._wall_density (fraction of disagreeing n.n. bonds)",
                        "record_t": REC_T},
           "precision_check": prec, "gating_L512_von_neumann": g}
    out["T2"] = {"rule": "|rho_bar(t) - rho_exact(t)| <= 0.03 rho_exact(t) for t in {100, 300, 1000, 3000}, "
                         "von Neumann, L = 512", "by_t": t2_rows, "pass": bool(t2_pass)}
    out["T1"] = {"rule": "OLS slope of 1/rho_bar on ln t over t = 300..3000 in [0.586, 0.688] (2/pi +/- 8%)",
                 "asymptotic_2_over_pi": float(2 / np.pi),
                 "exact_finite_t_slope": th["exact_slope_von_neumann_300_3000"], **sl, "pass": bool(t1_pass)}
    out["primary_rules"] = ["T2", "T1"]
    out["primary_pass"] = bool(t2_pass and t1_pass)
    if vn256:
        b256 = _block(vn256, ex_vn, [0] + REC_T)
        r512, r256 = g["by_t"]["3000"]["rho_bar"], b256["by_t"]["3000"]["rho_bar"]
        rel = abs(r256 - r512) / r512 if r512 > 0 else float("inf")
        out["finite_size_check"] = {"rule": "|rho_bar_256(3000) - rho_bar_512(3000)| <= 0.02 rho_bar_512(3000)",
                                    "reported_only": True, "L256_von_neumann": b256,
                                    "rho_bar_512_3000": r512, "rho_bar_256_3000": r256, "rel_diff": rel,
                                    "pass": bool(rel <= TOL_FS),
                                    "slope_L256": _slope_with_ci(vn256)}
    if mo256:
        bm = _block(mo256, ex_mo, [0] + REC_T)
        out["moore_secondary"] = {"reported_only": True, "L256_moore": bm, "slope": _slope_with_ci(mo256),
                                  "exact_finite_t_slope": th["exact_slope_moore_300_3000"],
                                  "asymptotic_4_over_3pi_derived": float(4 / (3 * np.pi))}
    allrecs = vn512 + vn256 + mo256
    out["deviations"] = common.load_deviations("P18")
    out["seeds"] = {"L512_von_neumann": f"0-{n_gate - 1}", "L256_von_neumann": "0-9", "L256_moore": "0-9"}
    out["code_commits_of_runs"] = sorted({r["git_commit"] for r in allrecs})
    out["all_runs_from_clean_tree"] = all(r["tracked_tree_clean"] for r in allrecs)
    out["analysis_git"] = common.git_info()
    out["bit_identity_checks"] = "analysis/outputs/repro_v2/bit_identity_checks.json"
    out["total_compute_runtime_s"] = float(sum(r["runtime_s"] for r in allrecs))
    out["timestamp"] = common.now_iso()
    common.write_json(common.OUT_DIR / "P18.json", out)
    print("P18 T2", "PASS" if t2_pass else "FAIL", "T1", "PASS" if t1_pass else "FAIL", "slope %.4f" % sl["slope"])


def main() -> None:
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run"); r.add_argument("--L", type=int, required=True)
    r.add_argument("--nb", required=True); r.add_argument("--seeds", required=True)
    sub.add_parser("precision"); sub.add_parser("theory"); sub.add_parser("analyze")
    a = ap.parse_args()
    if a.cmd == "run":
        lo, hi = a.seeds.split("-")
        cmd_run(a.L, a.nb, list(range(int(lo), int(hi) + 1)))
    elif a.cmd == "precision":
        cmd_precision()
    elif a.cmd == "theory":
        cmd_theory()
    else:
        cmd_analyze()


if __name__ == "__main__":
    main()
