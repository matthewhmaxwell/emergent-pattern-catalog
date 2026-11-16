"""Repro v2, P28: yard-sale model with redistribution.

P28-A (primary): Boghosian et al. 2017 Physica A 476:15, Eq. (61): G_inf = 1/(1 + 2.3682 tau^0.5943).
P28-B (reported): Boghosian 2014 PRE 89:042804, slope of 1/a(t) = 2.04046e-5.

Frozen spec: docs/validation_rebuild/repro_v2/group_C.md, P28 sections 3 (Protocol) and 4 (Tolerance and pass
rule). Nothing here changes a target, tolerance, size, seed count or rule. Details the spec leaves open are fixed
in IMPLEMENTATION_CHOICES before any full run, and the spec's conditional steps (stationarity rerun, finite-f
extrapolation) are applied automatically in the order the spec gives them.

The model is epc.models.yard_sale.YardSale unchanged (parameters only). Calling convention from the spec:
always step(R) with R = redistribute_every = N/2, so each tax follows its own R transactions.

Usage: PYTHONPATH=. python analysis/repro_v2/p28_yard_sale.py all
Per-run results stream to analysis/outputs/repro_v2/raw/p28_runs.jsonl (resumable).
"""
from __future__ import annotations

import argparse
import json
import math
import os
import subprocess
import sys
import time

import numpy as np
from scipy.optimize import brentq

sys.path.insert(0, ".")
from epc.models.yard_sale import YardSale  # noqa: E402

RAW = "analysis/outputs/repro_v2/raw/p28_runs.jsonl"
OUT_JSON = "analysis/outputs/repro_v2/P28.json"

# P28-A
N_A = 5000
R_A = N_A // 2
TAUS = [0.2, 0.5, 1.0]
SEEDS_A = list(range(2801, 2811))
MEAS_TBAR = 50.0
TOL_A = 0.02
HALF_TOL = 0.005
FCHECK_TOL = 0.01

# P28-B
N_B = 50000
F_B = 0.25
W0_B = 100.0
SEEDS_B = list(range(2821, 2831))
T_END_B = 80_000_000
CHECKPOINTS_B = [k * T_END_B // 16 for k in range(1, 17)]
BATCH_B = 1_000_000
SLOPE_SRC = 2.04046e-5
INTERCEPT_SRC = 3.93264

IMPLEMENTATION_CHOICES = [
    "P28-A model: YardSale(n_agents=5000, f, lambda_save=0, chi=tau*f^2, redistribute_every=2500, w0=1.0, "
    "init_mode='equal', seed); advanced only by step(2500), so each tax (w <- (1-chi) w + chi mean) follows its own "
    "2500 transactions. Time per step(2500) = 2 f^2 * 2500 / 5000 = f^2 in tau-bar units.",
    "Burn-in max(20/tau, 20) tau-bar = round(t_burn / f^2) blocks; measurement 50 tau-bar = round(50 / f^2) blocks; "
    "the Gini of the state returned by each step(2500) (i.e. after that block's tax) is sampled every block.",
    "Per-seed estimate = time average of the sampled Gini over the measurement window; reported value = mean over "
    "the 10 seeds, SE = sd/sqrt(10).",
    "Stationarity rule, applied to every P28-A run set (each tau, f, init): first-half and second-half time "
    "averages are formed per seed; the rule compares their 10-seed means. If they differ by more than 0.005, all 10 "
    "seeds of that set are rerun once with the burn-in doubled, and the rerun replaces the set (no further "
    "doubling). Per-seed half differences are also reported.",
    "Finite-f check: tau = 0.5 at f = 0.025 (chi = 3.125e-4, R = 2500), seeds 2801-2810, same tau-bar burn-in and "
    "window. If |G(0.05) - G(0.025)| > 0.01 (seed means), f = 0.025 is run at all three tau and the gated value at "
    "each tau is the linear extrapolation in f^2 to 0: G0 = G(0.025) + (G(0.025) - G(0.05)) * 0.025^2 / "
    "(0.05^2 - 0.025^2).",
    "Reported, not gated: tau = 0.5, f = 0.05, init_mode='exponential', seeds 2801-2810.",
    "P28-B model: YardSale(n_agents=50000, f=0.25, lambda_save=0, chi=0, w0=100, init_mode='equal', seed), seeds "
    "2821-2830, advanced in step(1e6) batches; the 16 equally spaced times are t = k * 5e6, k = 1..16 (up to 8e7).",
    "P28-B estimator: A(w) = fraction of agents with wealth > w, evaluated at every agent's wealth; least-squares "
    "fit A = alpha + beta ln w over agents with 0.2 <= A <= 0.8; ln w_min = (1 - alpha)/beta; a from Eq. (37) "
    "w_min = (W/2a) csch(N/2a) exp(-N/2a) with W = 5e6 and N = 5e4, solved in log form "
    "ln w_min = ln W - ln a - N/a - ln(1 - exp(-N/a)) (brentq on ln a); per seed, least squares of 1/a on t "
    "(t in transactions) gives slope and intercept; the rule uses the mean slope over the 10 seeds.",
]


def git_commit():
    try:
        return subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip()
    except Exception:
        return "unknown"


def g_fit(tau):
    return 1.0 / (1.0 + 2.3682 * tau ** 0.5943)


# ------------------------------------------------------------------ P28-A

def run_A(tau, f, seed, init="equal", burn_mult=1, meas_tbar=MEAS_TBAR, burn_tbar=None):
    t0 = time.time()
    chi = tau * f * f
    m = YardSale(n_agents=N_A, f=f, lambda_save=0.0, chi=chi, redistribute_every=R_A, w0=1.0,
                 init_mode=init, seed=seed)
    st = m.setup()
    g_init = st["gini"]
    tb_block = 2.0 * f * f * R_A / N_A
    t_burn = max(20.0 / tau, 20.0) * burn_mult if burn_tbar is None else burn_tbar
    n_burn = int(round(t_burn / tb_block))
    n_meas = int(round(meas_tbar / tb_block))
    for _ in range(n_burn):
        m.step(R_A)
    g = np.empty(n_meas)
    for k in range(n_meas):
        g[k] = m.step(R_A)["gini"]
    h = n_meas // 2
    dec = max(1, n_meas // 100)
    return {
        "kind": "A", "tau": tau, "f": f, "chi": chi, "seed": seed, "init": init, "burn_mult": burn_mult,
        "t_burn": t_burn, "n_burn_blocks": n_burn, "n_meas_blocks": n_meas, "transactions": (n_burn + n_meas) * R_A,
        "G_init": g_init, "G_mean": float(g.mean()), "G_first_half": float(g[:h].mean()),
        "G_second_half": float(g[h:].mean()), "G_sd": float(g.std()),
        "G_first_sample": float(g[0]), "G_series_decimated": g[::dec].tolist(), "decimation": dec,
        "total_wealth_end": float(m.wealth.sum()),
        "runtime_s": round(time.time() - t0, 1),
    }


def keyA(tau, f, seed, init, burn_mult):
    return ("A", tau, f, seed, init, burn_mult)


# ------------------------------------------------------------------ P28-B

def solve_a(ln_wmin, W, N):
    def F(ln_a):
        a = math.exp(ln_a)
        y = N / a
        return math.log(W) - ln_a - y - math.log(-math.expm1(-y)) - ln_wmin
    lo, hi = -20.0, 60.0
    if not (F(lo) < 0 < F(hi)):
        return float("nan")
    return math.exp(brentq(F, lo, hi, xtol=1e-14, rtol=1e-14, maxiter=500))


def fit_wmin(w):
    N = len(w)
    sw = np.sort(w)
    A = (N - np.searchsorted(sw, sw, side="right")) / N
    sel = (A >= 0.2) & (A <= 0.8)
    lw = np.log(sw[sel])
    beta, alpha = np.polyfit(lw, A[sel], 1)
    ln_wmin = (1.0 - alpha) / beta
    return {"alpha": float(alpha), "beta": float(beta), "ln_wmin": float(ln_wmin), "n_fit": int(sel.sum()),
            "ln_w_poorest": float(np.log(sw[0])), "ln_w_median": float(np.log(np.median(sw))),
            "ln_w_richest": float(np.log(sw[-1]))}


def run_B(seed):
    t0 = time.time()
    m = YardSale(n_agents=N_B, f=F_B, lambda_save=0.0, chi=0.0, w0=W0_B, init_mode="equal", seed=seed)
    m.setup()
    W = N_B * W0_B
    rows = []
    done = 0
    for tc in CHECKPOINTS_B:
        while done < tc:
            n = min(BATCH_B, tc - done)
            m.step(n)
            done += n
        fw = fit_wmin(m.wealth)
        a = solve_a(fw["ln_wmin"], W, N_B)
        rows.append({"t": tc, **fw, "a": a, "inv_a": (1.0 / a if a == a else float("nan")),
                     "total_wealth": float(m.wealth.sum())})
    t = np.array([r["t"] for r in rows], float)
    ia = np.array([r["inv_a"] for r in rows], float)
    ok = np.isfinite(ia)
    if ok.sum() >= 2:
        slope, intercept = np.polyfit(t[ok], ia[ok], 1)
    else:
        slope = intercept = float("nan")
    return {"kind": "B", "seed": seed, "checkpoints": rows, "slope": float(slope), "intercept": float(intercept),
            "n_points_fit": int(ok.sum()), "runtime_s": round(time.time() - t0, 1)}


# ------------------------------------------------------------------ bookkeeping

def load_raw():
    recs = {}
    if os.path.exists(RAW):
        with open(RAW) as fh:
            for line in fh:
                line = line.strip()
                if line:
                    d = json.loads(line)
                    if d["kind"] == "A":
                        recs[keyA(d["tau"], d["f"], d["seed"], d["init"], d["burn_mult"])] = d
                    else:
                        recs[("B", d["seed"])] = d
    return recs


def append(rec, commit):
    rec["commit"] = commit
    os.makedirs(os.path.dirname(RAW), exist_ok=True)
    with open(RAW, "a") as fh:
        fh.write(json.dumps(rec) + "\n")


def ensure_A(tau, f, init, burn_mult, commit):
    recs = load_raw()
    for s in SEEDS_A:
        if keyA(tau, f, s, init, burn_mult) not in recs:
            r = run_A(tau, f, s, init=init, burn_mult=burn_mult)
            append(r, commit)
            print(f"[p28] A tau={tau} f={f} init={init} burn x{burn_mult} seed={s}: G={r['G_mean']:.4f} "
                  f"halves {r['G_first_half']:.4f}/{r['G_second_half']:.4f} ({r['runtime_s']} s)", flush=True)
    return load_raw()


def set_stats(recs, tau, f, init, burn_mult):
    rs = [recs[keyA(tau, f, s, init, burn_mult)] for s in SEEDS_A]
    g = np.array([r["G_mean"] for r in rs])
    h1 = np.array([r["G_first_half"] for r in rs])
    h2 = np.array([r["G_second_half"] for r in rs])
    return {
        "tau": tau, "f": f, "init": init, "burn_mult": burn_mult, "t_burn": rs[0]["t_burn"],
        "G_mean": float(g.mean()), "G_se": float(g.std(ddof=1) / math.sqrt(len(g))), "G_per_seed": g.tolist(),
        "half1_mean": float(h1.mean()), "half2_mean": float(h2.mean()),
        "half_diff": float(h2.mean() - h1.mean()), "half_rule_ok": bool(abs(h2.mean() - h1.mean()) <= HALF_TOL),
        "per_seed_half_diff": (h2 - h1).tolist(),
        "cpu_seconds": float(sum(r["runtime_s"] for r in rs)),
        "transactions_per_seed": rs[0]["transactions"],
    }


def ensure_stationary_set(tau, f, init, commit, log):
    """Run a P28-A set; apply the stationarity rule (double the burn-in once)."""
    recs = ensure_A(tau, f, init, 1, commit)
    st = set_stats(recs, tau, f, init, 1)
    entry = {"first": st}
    if not st["half_rule_ok"]:
        print(f"[p28] stationarity rule triggered for tau={tau} f={f} init={init}: halves differ by "
              f"{st['half_diff']:.4f}; rerunning with doubled burn-in", flush=True)
        recs = ensure_A(tau, f, init, 2, commit)
        st2 = set_stats(recs, tau, f, init, 2)
        entry["rerun_burn_x2"] = st2
        entry["used"] = st2
    else:
        entry["used"] = st
    log.append(entry)
    return entry["used"]


def run_all(commit):
    log = []
    used = {}
    for tau in TAUS:
        used[(tau, 0.05)] = ensure_stationary_set(tau, 0.05, "equal", commit, log)
    used[(0.5, 0.025)] = ensure_stationary_set(0.5, 0.025, "equal", commit, log)
    dG = abs(used[(0.5, 0.05)]["G_mean"] - used[(0.5, 0.025)]["G_mean"])
    trig = dG > FCHECK_TOL
    print(f"[p28] finite-f check: |G(0.05)-G(0.025)| = {dG:.4f}; extrapolation triggered = {trig}", flush=True)
    if trig:
        for tau in (0.2, 1.0):
            used[(tau, 0.025)] = ensure_stationary_set(tau, 0.025, "equal", commit, log)
    exp_set = ensure_stationary_set(0.5, 0.05, "exponential", commit, log)
    # P28-B
    recs = load_raw()
    for s in SEEDS_B:
        if ("B", s) not in recs:
            r = run_B(s)
            append(r, commit)
            print(f"[p28] B seed={s}: slope={r['slope']:.4e} intercept={r['intercept']:.4g} ({r['runtime_s']} s)",
                  flush=True)
    return analyze(commit)


def analyze(commit):
    recs = load_raw()
    log, used = [], {}
    # rebuild the decisions deterministically from the raw records (same rules as run_all)
    def used_set(tau, f, init):
        st = set_stats(recs, tau, f, init, 1)
        entry = {"first": st}
        if not st["half_rule_ok"]:
            st2 = set_stats(recs, tau, f, init, 2)
            entry["rerun_burn_x2"] = st2
            entry["used"] = st2
        else:
            entry["used"] = st
        log.append(entry)
        return entry["used"]
    for tau in TAUS:
        used[(tau, 0.05)] = used_set(tau, 0.05, "equal")
    used[(0.5, 0.025)] = used_set(0.5, 0.025, "equal")
    dG = used[(0.5, 0.05)]["G_mean"] - used[(0.5, 0.025)]["G_mean"]
    trig = abs(dG) > FCHECK_TOL
    if trig:
        for tau in (0.2, 1.0):
            used[(tau, 0.025)] = used_set(tau, 0.025, "equal")
    exp_set = used_set(0.5, 0.05, "exponential")
    rows = []
    for tau in TAUS:
        g05 = used[(tau, 0.05)]
        if trig:
            g025 = used[(tau, 0.025)]
            k = 0.025 ** 2 / (0.05 ** 2 - 0.025 ** 2)
            G = g025["G_mean"] + (g025["G_mean"] - g05["G_mean"]) * k
            se = math.sqrt(((1 + k) * g025["G_se"]) ** 2 + (k * g05["G_se"]) ** 2)
            src = "f^2 -> 0 extrapolation from f = 0.05 and 0.025"
        else:
            G, se, src = g05["G_mean"], g05["G_se"], "f = 0.05"
        tgt = g_fit(tau)
        rows.append({"tau": tau, "G_target_Eq61": tgt, "G_sim": G, "G_sim_se": se, "gated_value_from": src,
                     "deviation": G - tgt, "within_0.02": bool(abs(G - tgt) <= TOL_A),
                     "G_target_at_2tau_reference": g_fit(2 * tau)})
    pass_A = all(r["within_0.02"] for r in rows)
    # P28-B
    Bs = [recs[("B", s)] for s in SEEDS_B if ("B", s) in recs]
    slopes = np.array([b["slope"] for b in Bs]) if Bs else np.array([])
    B = None
    if len(Bs) == len(SEEDS_B):
        ms = float(np.mean(slopes))
        B = {
            "slope_target": SLOPE_SRC, "band": [1.73e-5, 2.35e-5],
            "mean_slope": ms, "slope_se": float(np.std(slopes, ddof=1) / math.sqrt(len(slopes))),
            "slopes_per_seed": slopes.tolist(),
            "mean_intercept": float(np.mean([b["intercept"] for b in Bs])),
            "intercept_source_reported_only": INTERCEPT_SRC,
            "PASS_B": bool(1.73e-5 <= ms <= 2.35e-5),
            "per_seed": [{"seed": b["seed"], "slope": b["slope"], "intercept": b["intercept"],
                          "checkpoints": b["checkpoints"], "runtime_s": b["runtime_s"]} for b in Bs],
        }
    cpu = {}
    for k, r in recs.items():
        name = f"A_tau{r['tau']}_f{r['f']}_{r['init']}_burnx{r['burn_mult']}" if r["kind"] == "A" else "B"
        cpu[name] = round(cpu.get(name, 0.0) + r["runtime_s"], 1)
    out = {
        "pattern": "P28",
        "sources": {
            "P28-A": "Boghosian, Devitt-Lee, Johnson, Marcq & Wang 2017, Physica A 476:15, Eq. (61)",
            "P28-B": "Boghosian 2014, PRE 89:042804, Fig. 5 fit 1/a(t) = 3.93264 + 0.0000204046 t",
        },
        "spec": "docs/validation_rebuild/repro_v2/group_C.md (P28 sections 3-4), frozen at 39688e3",
        "implementation_choices": IMPLEMENTATION_CHOICES,
        "P28_A": {
            "rows": rows,
            "stationarity_log": log,
            "finite_f_check": {"G_f0.05": used[(0.5, 0.05)]["G_mean"], "G_f0.025": used[(0.5, 0.025)]["G_mean"],
                               "diff": dG, "threshold": FCHECK_TOL, "extrapolation_triggered": bool(trig)},
            "PASS_A": bool(pass_A),
        },
        "P28_B": B if B else "incomplete",
        "reported_only": {
            "exponential_init_tau0.5": {"G_mean": exp_set["G_mean"], "G_se": exp_set["G_se"],
                                        "equal_init_G_mean": used[(0.5, 0.05)]["G_mean"],
                                        "G_target": g_fit(0.5)},
            "P28_B": "reported separately; does not gate the P28 verdict",
        },
        "PASS": bool(pass_A),
        "runtime_cpu_seconds": cpu,
        "code_commits_of_runs": sorted({r["commit"] for r in recs.values()}),
        "commit_at_analysis": commit,
        "numpy_version": np.__version__,
    }
    os.makedirs(os.path.dirname(OUT_JSON), exist_ok=True)
    with open(OUT_JSON, "w") as fh:
        json.dump(out, fh, indent=1)
    print(json.dumps({"PASS_A": pass_A, "rows": rows, "finite_f": out["P28_A"]["finite_f_check"],
                      "B_mean_slope": (B or {}).get("mean_slope"), "PASS_B": (B or {}).get("PASS_B")}, indent=1))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["all", "analyze", "smoke"])
    args = ap.parse_args()
    commit = git_commit()
    if args.cmd == "smoke":
        # Timing only: 0.125 + 0.125 tau-bar from the equal start (protocol: 40 + 50), seed outside the protocol.
        r = run_A(0.5, 0.05, 999, meas_tbar=0.125, burn_tbar=0.125)
        r.pop("G_series_decimated")
        print(json.dumps({k: r[k] for k in ("n_burn_blocks", "n_meas_blocks", "transactions", "runtime_s",
                                            "total_wealth_end")}, indent=1))
        print("us per transaction:", r["runtime_s"] / r["transactions"] * 1e6)
        # Estimator check on synthetic data (no model output): P_c(w) = a/w on [w_min, w_max] with Eqs. (35)-(37).
        W, N = 5e6, 50000
        for a_true in (300.0, 3000.0, 30000.0):
            y = N / a_true
            ln_wmin = math.log(W) - math.log(a_true) - y - math.log(-math.expm1(-y))
            ln_wmax = ln_wmin + y
            lw = np.random.default_rng(1).uniform(ln_wmin, ln_wmax, size=N)
            fw = fit_wmin(np.exp(lw))
            print(f"synthetic a={a_true}: Eq37 ln w_min={ln_wmin:.4f}, fitted {fw['ln_wmin']:.4f}, "
                  f"solved a={solve_a(fw['ln_wmin'], W, N):.2f}, solved a from exact w_min={solve_a(ln_wmin, W, N):.6f}")
        return
    if args.cmd == "all":
        run_all(commit)
    else:
        analyze(commit)


if __name__ == "__main__":
    main()
