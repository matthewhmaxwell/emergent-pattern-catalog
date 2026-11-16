"""Repro v2, P12: Reichenbach, Mobilia & Frey (2007) Nature 448:1046, Fig. 2b.

Frozen spec: docs/validation_rebuild/repro_v2/group_C.md, P12 sections 3 (Protocol)
and 4 (Tolerance and pass rule). Nothing here changes a target, tolerance, size,
seed count or rule; where the spec leaves a detail open, the choice is fixed in
this file before any full run and listed in IMPLEMENTATION_CHOICES.

Run (one process per call; resumable; results streamed one JSON line per realisation):
    PYTHONPATH=. python analysis/repro_v2/p12_extinction.py run --variant primary --L 40
    PYTHONPATH=. python analysis/repro_v2/p12_extinction.py run --variant primary --L 100 [--shard 0/2]
    PYTHONPATH=. python analysis/repro_v2/p12_extinction.py run --variant D2 --L 100
    PYTHONPATH=. python analysis/repro_v2/p12_extinction.py queue --jobs primary:100,D2:40,D2:100   (n workers)
    PYTHONPATH=. python analysis/repro_v2/p12_extinction.py analyze
"""
from __future__ import annotations

import argparse
import glob
import json
import math
import os
import subprocess
import sys
import time

import numpy as np

sys.path.insert(0, ".")
from epc.models.rps_spatial import RPSSpatialModel  # noqa: E402

RAW_DIR = "analysis/outputs/repro_v2/raw/p12"
OUT_JSON = "analysis/outputs/repro_v2/P12.json"

L_VALUES = [40, 100]                                   # index(L) = 0, 1 (optional L = 200 would be 2)
M_EXPONENTS = list(range(-1, 6))                       # M = 1e-4 * 2^(k/2), k = -1..5
M_GRID = [1e-4 * 2.0 ** (k / 2.0) for k in M_EXPONENTS]
M_TABLE = ["7.07e-5", "1.00e-4", "1.41e-4", "2.00e-4", "2.83e-4", "4.00e-4", "5.66e-4"]
P_FIG = {
    40: [0.427, 0.599, 0.776, 0.908, 0.945, 0.976, 0.987],
    100: [0.247, 0.407, 0.507, 0.611, 0.734, 0.911, 0.959],
    200: [0.104, 0.167, 0.210, 0.309, 0.534, 0.830, 0.964],
}
N_REAL = 500
SIGMA = MU = 1.0

IMPLEMENTATION_CHOICES = [
    "M grid uses the exact values 1e-4 * 2^(k/2), k = -1..5 (the spec's table prints them rounded to 3 digits).",
    "Seed of realisation r at (L, M) = 120000 + 1000*index(L) + 10*index(M) + r, r = 0..499, index(L) = 0 for L = 40 "
    "and 1 for L = 100, index(M) = 0..6 in ascending M, exactly as the spec's formula. Note: with r up to 499 the "
    "formula's seed ranges overlap ACROSS M values within one L (e.g. (index(M)=1, r=0) and (index(M)=0, r=10) share "
    "seed 120010), so neighbouring M points share initial grids and random streams (common random numbers); within "
    "one (L, M) point all 500 seeds are distinct, and L = 40 and L = 100 never share seeds. The spec's remark 'so "
    "streams never overlap' holds across L only. The formula is applied as written.",
    "G_max = N(2 + M N/2) is not an integer; each realisation runs floor(G_max) generations. Absorption is checked "
    "at the end of every generation (RPSSpatialModel.run), so 'absorbed by G_max' == absorbed at a generation "
    "<= floor(G_max).",
    "Each realisation: RPSSpatialModel(rows=L, cols=L, mobility=M, selection_rate=1, reproduction_rate=1, "
    "neighborhood='von_neumann', init_mode='random', seed, kernel='numba'); run(n_steps=floor(G_max), "
    "record_every=floor(G_max)+1). kernel='numba' is the compiled loop checked bit-identical to the reference loop "
    "(tests/test_rps_numba_kernel.py).",
    "Diagnostic D1 (G = N) is read from the primary runs: absorbed by generation N, i.e. absorption generation <= N "
    "(the primary runs use the same seeds and M and run at least 2N generations).",
    "Diagnostic D2 = our code at mobility M/2 run to floor(N(2 + (M/2)N/2)) generations, same seeds as the primary "
    "point at M. The diagnostics are evaluated with the same three conditions as the primary rule and are "
    "reported only.",
    "Optional L = 200 (spec: 'run only if the budget allows', 30-45 CPU-h) is not run.",
    "Wilson 95% interval with z = 1.959964.",
]


def seed_for(iL: int, iM: int, r: int) -> int:
    return 120000 + 1000 * iL + 10 * iM + r


def g_max(L: int, M: float) -> float:
    N = L * L
    return N * (SIGMA + MU + M * N / 2.0)


def git_commit() -> str:
    try:
        return subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip()
    except Exception:
        return "unknown"


def git_dirty() -> bool:
    try:
        out = subprocess.check_output(["git", "status", "--porcelain", "--untracked-files=no"], text=True)
        return bool(out.strip())
    except Exception:
        return True


def run_one(L: int, M_run: float, G: int, seed: int) -> dict:
    t0 = time.time()
    m = RPSSpatialModel(rows=L, cols=L, mobility=M_run, selection_rate=SIGMA, reproduction_rate=MU,
                        neighborhood="von_neumann", init_mode="random", seed=seed, kernel="numba")
    m.run(n_steps=G, record_every=G + 1)
    absorbed = bool(m.is_converged())
    grid = m._grid
    counts = [int((grid == s).sum()) for s in (1, 2, 3)]
    return {
        "absorbed": absorbed,
        "generations": int(m.step_count),
        "final_counts_ABC": counts,
        "runtime_s": round(time.time() - t0, 3),
    }


def task_list(variant: str, L: int):
    iL = L_VALUES.index(L)
    tasks = []
    for r in range(N_REAL):
        for iM, M in enumerate(M_GRID):
            seed = seed_for(iL, iM, r)
            if variant == "primary":
                M_run = M
            elif variant == "D2":
                M_run = M / 2.0
            else:
                raise ValueError(variant)
            G = int(math.floor(g_max(L, M_run)))
            tasks.append({"variant": variant, "L": L, "iM": iM, "M": M, "M_run": M_run, "r": r,
                          "seed": seed, "G_max": G, "G_max_exact": g_max(L, M_run)})
    return tasks


def done_keys(variant: str, L: int):
    keys = set()
    for fn in glob.glob(os.path.join(RAW_DIR, f"{variant}_L{L}_*.jsonl")):
        with open(fn) as fh:
            for line in fh:
                line = line.strip()
                if not line:
                    continue
                d = json.loads(line)
                keys.add((d["iM"], d["r"]))
    return keys


def cmd_run(args):
    os.makedirs(RAW_DIR, exist_ok=True)
    shard_i, shard_n = (int(x) for x in args.shard.split("/"))
    tasks = task_list(args.variant, args.L)
    done = done_keys(args.variant, args.L)
    # Shards split the tasks still to do (start shards of one batch together).
    todo = [t for t in tasks if (t["iM"], t["r"]) not in done][shard_i::shard_n]
    if args.limit:
        todo = todo[: args.limit]
    out_fn = os.path.join(RAW_DIR, f"{args.variant}_L{args.L}_{time.strftime('%Y%m%dT%H%M%S')}_shard{shard_i}of{shard_n}.jsonl")
    commit = git_commit()
    print(f"[p12] {args.variant} L={args.L} shard {args.shard}: {len(todo)} tasks to run, {len(done)} of {len(tasks)} done "
          f"(commit {commit[:10]}, dirty={git_dirty()})", flush=True)
    t_start = time.time()
    with open(out_fn, "a") as fh:
        for k, t in enumerate(todo):
            res = run_one(t["L"], t["M_run"], t["G_max"], t["seed"])
            rec = {**t, **res, "commit": commit}
            fh.write(json.dumps(rec) + "\n")
            fh.flush()
            if (k + 1) % 50 == 0 or k + 1 == len(todo):
                el = time.time() - t_start
                print(f"[p12] {k + 1}/{len(todo)} done, {el / 3600:.2f} h elapsed, "
                      f"ETA {el / (k + 1) * (len(todo) - k - 1) / 3600:.2f} h", flush=True)


CLAIM_DIR = os.path.join(RAW_DIR, "claims")


def cmd_queue(args):
    """Work queue for several workers: tasks are claimed with O_EXCL lock files, so any number of workers
    started on the same job list split the remaining tasks dynamically. The per-realisation computation
    (run_one) is the same as in `run`. Clear CLAIM_DIR only when no worker is running."""
    os.makedirs(RAW_DIR, exist_ok=True)
    os.makedirs(CLAIM_DIR, exist_ok=True)
    commit = git_commit()
    wid = f"{os.getpid()}"
    for job in args.jobs.split(","):
        variant, L = job.split(":")
        L = int(L)
        tasks = task_list(variant, L)
        done = done_keys(variant, L)
        todo = [t for t in tasks if (t["iM"], t["r"]) not in done]
        out_fn = os.path.join(RAW_DIR, f"{variant}_L{L}_{time.strftime('%Y%m%dT%H%M%S')}_q{wid}.jsonl")
        print(f"[p12q {wid}] job {variant} L={L}: {len(todo)} tasks not done (commit {commit[:10]}, "
              f"dirty={git_dirty()})", flush=True)
        n_run = 0
        t_start = time.time()
        with open(out_fn, "a") as fh:
            for t in todo:
                claim = os.path.join(CLAIM_DIR, f"{variant}_L{L}_{t['iM']}_{t['r']}")
                try:
                    os.close(os.open(claim, os.O_CREAT | os.O_EXCL | os.O_WRONLY))
                except FileExistsError:
                    continue
                res = run_one(t["L"], t["M_run"], t["G_max"], t["seed"])
                fh.write(json.dumps({**t, **res, "commit": commit}) + "\n")
                fh.flush()
                n_run += 1
                if n_run % 50 == 0:
                    left = len(todo) - len(done_keys(variant, L) & {(x["iM"], x["r"]) for x in todo})
                    print(f"[p12q {wid}] {variant} L={L}: {n_run} run by this worker, ~{left} left, "
                          f"{(time.time() - t_start) / 3600:.2f} h", flush=True)
        print(f"[p12q {wid}] job {variant} L={L} finished by this worker ({n_run} run)", flush=True)
    print(f"[p12q {wid}] QUEUE_DONE", flush=True)


# ---------------------------------------------------------------- analysis

def wilson(k: int, n: int, z: float = 1.959964):
    if n == 0:
        return (float("nan"), float("nan"))
    p = k / n
    den = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / den
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / den
    return (centre - half, centre + half)


def load_records(variant: str, L: int):
    recs = {}
    for fn in sorted(glob.glob(os.path.join(RAW_DIR, f"{variant}_L{L}_*.jsonl"))):
        with open(fn) as fh:
            for line in fh:
                line = line.strip()
                if line:
                    d = json.loads(line)
                    recs[(d["iM"], d["r"])] = d
    return recs


def estimate(recs, L: int, cutoff=None):
    """P_ext per M. cutoff=None: absorbed by the run's own G_max; cutoff=G: absorbed at generation <= G."""
    rows = []
    for iM, M in enumerate(M_GRID):
        rs = [recs[(iM, r)] for r in range(N_REAL) if (iM, r) in recs]
        n = len(rs)
        if cutoff is None:
            k = sum(1 for d in rs if d["absorbed"])
        else:
            k = sum(1 for d in rs if d["absorbed"] and d["generations"] <= cutoff)
        p = k / n if n else float("nan")
        lo, hi = wilson(k, n)
        gens = [d["generations"] for d in rs if d["absorbed"]]
        rows.append({
            "M": M, "M_table": M_TABLE[iM], "n": n, "k_absorbed": k, "P_hat": p, "wilson95": [lo, hi],
            "P_fig": P_FIG[L][iM], "deviation": p - P_FIG[L][iM],
            "abs_dev_le_0.10": bool(abs(p - P_FIG[L][iM]) <= 0.10),
            "G_max": (rs[0]["G_max"] if (cutoff is None and rs) else cutoff),
            "absorption_generation_quartiles": (np.percentile(gens, [25, 50, 75]).tolist() if gens else None),
            "absorption_generation_hist": _hist(gens, rs[0]["G_max"] if rs else 1),
        })
    devs = [row["deviation"] for row in rows]
    return {
        "points": rows,
        "all_within_0.10": all(row["abs_dev_le_0.10"] for row in rows),
        "max_abs_deviation": max(abs(d) for d in devs),
        "mean_signed_deviation": float(np.mean(devs)),
        "mean_signed_within_0.05": bool(abs(float(np.mean(devs))) <= 0.05),
        "complete": all(row["n"] == N_REAL for row in rows),
    }


def _hist(gens, gmax):
    if not gens:
        return {"edges": [], "counts": []}
    edges = np.linspace(0, gmax, 11)
    c, e = np.histogram(gens, bins=edges)
    return {"edges": [float(x) for x in e], "counts": [int(x) for x in c]}


def rule(est40, est100):
    c1 = est100["all_within_0.10"] and est100["mean_signed_within_0.05"]
    c2 = est40["all_within_0.10"] and est40["mean_signed_within_0.05"]
    sign = []
    for iM, M in enumerate(M_GRID):
        if M <= 2.83e-4 + 1e-9:
            p40 = est40["points"][iM]["P_hat"]
            p100 = est100["points"][iM]["P_hat"]
            sign.append({"M": M, "P_L100": p100, "P_L40": p40, "L100_lt_L40": bool(p100 < p40)})
    c3 = all(s["L100_lt_L40"] for s in sign)
    return {
        "condition_1_L100": bool(c1),
        "condition_2_L40": bool(c2),
        "condition_3_sharpening": bool(c3),
        "sharpening_points": sign,
        "PASS": bool(c1 and c2 and c3),
    }


def cmd_analyze(args):
    meta = {"commit_at_analysis": git_commit()}
    prim = {L: load_records("primary", L) for L in L_VALUES}
    est = {L: estimate(prim[L], L) for L in L_VALUES}
    primary_rule = rule(est[40], est[100])
    # ordering at all 7 points (reported)
    ordering = [{"M": M, "P_L40": est[40]["points"][i]["P_hat"], "P_L100": est[100]["points"][i]["P_hat"]}
                for i, M in enumerate(M_GRID)]
    # Mc statement: L = 100 points at M >= 4e-4 within tolerance (L = 200 not run)
    mc_points = [est[100]["points"][i] for i, M in enumerate(M_GRID) if M >= 4e-4 - 1e-12]
    mc_repro = all(p["abs_dev_le_0.10"] for p in mc_points)
    # D1
    d1 = {L: estimate(prim[L], L, cutoff=L * L) for L in L_VALUES}
    d1_rule = rule(d1[40], d1[100])
    # D2
    d2recs = {L: load_records("D2", L) for L in L_VALUES}
    d2 = {L: (estimate(d2recs[L], L) if d2recs[L] else None) for L in L_VALUES}
    d2_rule = rule(d2[40], d2[100]) if (d2[40] and d2[100]) else None
    runtimes = {}
    commits = set()
    for name, rr in [("primary", prim), ("D2", d2recs)]:
        for L in L_VALUES:
            vals = [d["runtime_s"] for d in rr[L].values()]
            commits |= {d["commit"] for d in rr[L].values()}
            runtimes[f"{name}_L{L}"] = {"n_runs": len(vals), "cpu_seconds_total": round(float(np.sum(vals)), 1)}
    import numba
    out = {
        "pattern": "P12",
        "source": "Reichenbach, Mobilia & Frey 2007, Nature 448:1046, Fig. 2b (digitized in group_C.md)",
        "spec": "docs/validation_rebuild/repro_v2/group_C.md (P12 sections 3-4), frozen at 39688e3",
        "protocol": {
            "sigma": SIGMA, "mu": MU, "neighborhood": "von_neumann", "init": "random (each of 4 states w.p. 1/4)",
            "L_values": L_VALUES, "M_grid": M_GRID, "realizations_per_point": N_REAL,
            "G_max": "floor(N(2 + M N/2)) generations of N elementary steps",
            "seed_formula": "120000 + 1000*index(L) + 10*index(M) + r",
        },
        "implementation_choices": IMPLEMENTATION_CHOICES,
        "primary": {"L40": est[40], "L100": est[100]},
        "primary_rule": primary_rule,
        "PASS": primary_rule["PASS"],
        "reported_only": {
            "ordering_all_M": ordering,
            "Mc_statement_reproduced_L100_only": bool(mc_repro),
            "Mc_note": "L = 200 not run; the statement is judged on the L = 100 points at M >= 4e-4 only.",
            "D1_G_equals_N": {"L40": d1[40], "L100": d1[100], "rule": d1_rule},
            "D2_M_half": ({"L40": d2[40], "L100": d2[100], "rule": d2_rule} if d2_rule else "not run"),
            "optional_L200": "not run (budget)",
        },
        "runtime": runtimes,
        "code_commits_of_runs": sorted(commits),
        "numba_version": numba.__version__,
        "numpy_version": np.__version__,
        **meta,
    }
    os.makedirs(os.path.dirname(OUT_JSON), exist_ok=True)
    with open(OUT_JSON, "w") as fh:
        json.dump(out, fh, indent=1)
    print(json.dumps({"PASS": out["PASS"], "rule": primary_rule,
                      "L40": [(p["M_table"], round(p["P_hat"], 3), p["P_fig"]) for p in est[40]["points"]],
                      "L100": [(p["M_table"], round(p["P_hat"], 3), p["P_fig"]) for p in est[100]["points"]],
                      "mean_dev": (est[40]["mean_signed_deviation"], est[100]["mean_signed_deviation"])},
                     indent=1))


def main():
    ap = argparse.ArgumentParser()
    sub = ap.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("--variant", choices=["primary", "D2"], required=True)
    r.add_argument("--L", type=int, choices=L_VALUES, required=True)
    r.add_argument("--shard", default="0/1")
    r.add_argument("--limit", type=int, default=0)
    q = sub.add_parser("queue")
    q.add_argument("--jobs", required=True, help="comma list of variant:L, e.g. primary:100,D2:40,D2:100")
    sub.add_parser("analyze")
    args = ap.parse_args()
    if args.cmd == "run":
        cmd_run(args)
    elif args.cmd == "queue":
        cmd_queue(args)
    else:
        cmd_analyze(args)


if __name__ == "__main__":
    main()
