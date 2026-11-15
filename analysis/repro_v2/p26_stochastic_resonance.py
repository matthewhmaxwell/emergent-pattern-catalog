"""Repro v2, P26: Gammaitoni, Haenggi, Jung & Marchesoni (1998) RMP 70:223, Fig. 7 (A = 0.1, Omega = 0.1).

Frozen spec: docs/validation_rebuild/repro_v2/group_C.md, P26 sections 3 (Protocol) and 4 (Tolerance and
pass rule). Nothing here changes a target, tolerance, size, seed count or rule; details the spec leaves open
are fixed in IMPLEMENTATION_CHOICES before any full run. The script runs the phases in the spec's order and
applies the spec's conditional steps (dt rule, frame check) automatically:

    main      a = b = 1 frame, dt = 0.01, 15 D values x 10 seeds x 400 trials
    dtcheck   D in {0.10, 0.15, 0.30} at dt = 0.005, same seeds
    dtrerun   only if the dt rule triggers: dt = 0.005 over the D window [0.12, 0.18] (+ D = 0.20)
    frame     canonical frame (a = 4, b = 1, A = 0.8, Omega = 0.4, D x 16, dt = 0.0025), D in {0.10, 0.15, 0.30}
    analyze   writes analysis/outputs/repro_v2/P26.json

Usage: PYTHONPATH=. python analysis/repro_v2/p26_stochastic_resonance.py all
Per-run results stream to analysis/outputs/repro_v2/raw/p26_runs.jsonl (resumable).
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

sys.path.insert(0, ".")
from epc.models.stochastic_resonance import DoubleWellParams, double_well_fourier_driver  # noqa: E402

RAW = "analysis/outputs/repro_v2/raw/p26_runs.jsonl"
OUT_JSON = "analysis/outputs/repro_v2/P26.json"

A_SRC = 0.1
OMEGA_SRC = 0.1
D_GRID = [0.06, 0.08, 0.10, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.20, 0.25, 0.30, 0.40, 0.50]
D_WINDOW = [0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18]
D_CHECK = [0.10, 0.15, 0.30]
SEEDS = list(range(1001, 1011))
N_TRIALS = 400
BURN_PERIODS = 10
MEAS_PERIODS = 50
ANNULUS = 5
BOOT_N = 1000
BOOT_SEED = 2626

TARGET = {"eta_max": 12.23, "D_max": 0.150, "eta_020": 10.8, "eta_030": 6.5}

IMPLEMENTATION_CHOICES = [
    "Signal frequency f = Omega/(2 pi) computed exactly (spec prints f = 0.0159155 for a = b = 1 and 0.063662 for "
    "the canonical frame); the projection frequency is Omega = 2 pi f.",
    "Step counts: n_burn = round(10 periods / dt), n_meas = round(50 periods / dt) (62,832 and 314,159 at dt = 0.01, "
    "as the spec states; 125,664 and 628,319 at dt = 0.005; the canonical frame at dt = 0.0025 has the same step "
    "counts as a = b = 1 at dt = 0.01). T = n_meas * dt.",
    "Random numbers: run (frame, dt, D, seed) uses np.random.default_rng(seed) and draws standard normals step by "
    "step for the 400 trials (row = step). 'Same seeds' for the dt and frame checks means the same seed numbers. "
    "Because the canonical frame at dt = 0.0025 is an exact rescaling of the a = b = 1 Euler map at dt = 0.01, the "
    "frame check with the same seeds reproduces the main runs up to rounding: it tests the scaling and the code "
    "path, not an independent statistical replicate. At dt = 0.005 the increments are different draws, so the dt "
    "comparison includes seed noise.",
    "Trial k starts in the well -1 (k even) or +1 (k odd), as BistableDoubleWell.simulate().",
    "Estimator exactly as spec: c_k = (2/T) sum x(t) exp(-i Omega t) dt over the window; per seed M = mean_k c_k, "
    "s2 = sum_k |c_k - M|^2 / (n - 1); eta = (|M|^2 - s2/n) / A^2. Mean and SE (sd/sqrt(10)) over the 10 seeds.",
    "dt rule: rel_k = |eta_dt0.005(D_k) - eta_dt0.01(D_k)| / eta_dt0.01(D_k) for D_k in {0.10, 0.15, 0.30} "
    "(seed means). If any rel_k > 0.02, every gated value is taken from dt = 0.005: dt = 0.005 is run over the "
    "window D in [0.12, 0.18] and also at D = 0.20 (needed for condition 3); D = 0.30 at dt = 0.005 comes from the "
    "check itself.",
    "Frame check: 16 * eta_canonical vs eta (a = b = 1, dt = 0.01) at each check D; agreement means "
    "|diff| <= 2 * sqrt(SE1^2 + SE2^2). A failure is recorded and the gating evaluation is withheld "
    "(spec: investigate before gating).",
    "Peak: np.polyfit(D, seed-mean eta, 2) on the 7 window points; D_max = -b/(2a), eta_max = c - b^2/(4a) "
    "(a < 0 required). Bootstrap: 1000 resamples of the 10 seeds with replacement (same resample at every D), "
    "rng seed 2626; SE = sd of the bootstrap vertex values.",
    "SNR (reported only): per trial X_k(w_m) = sum x exp(-i w_m t) dt at the 10 bins w_m = 2 pi m / T, "
    "m = m0 +- 1..5 (m0 = 50, the bin of Omega); background S_N = mean over bins of "
    "sum_k |X_k - mean_k X_k|^2 / ((n - 1) T); SNR = pi * (|M|^2 - s2/n) / S_N (RMP Eq. 2.13 convention with "
    "S(w) = int e^{i w tau} <x x> d tau). Reported: grid D of local maxima of the seed-mean SNR and a 3-point "
    "quadratic vertex around the interior local maximum, against Delta V / 2 = 0.125.",
]


def git_commit():
    try:
        return subprocess.check_output(["git", "rev-parse", "HEAD"], text=True).strip()
    except Exception:
        return "unknown"


def frame_params(frame: str, dt: float):
    if frame == "unit":
        a, b, A, Om, Dscale = 1.0, 1.0, A_SRC, OMEGA_SRC, 1.0
    elif frame == "canonical":
        a, b, A, Om, Dscale = 4.0, 1.0, 8.0 * A_SRC, 4.0 * OMEGA_SRC, 16.0
    else:
        raise ValueError(frame)
    p = DoubleWellParams(a=a, b=b, signal_amplitude=A, signal_frequency=Om / (2.0 * math.pi), dt=dt)
    return p, Om, Dscale


def run_one(frame: str, dt: float, D_src: float, seed: int, n_trials: int = N_TRIALS,
            burn_periods: float = BURN_PERIODS, meas_periods: float = MEAS_PERIODS):
    p, Om, Dscale = frame_params(frame, dt)
    D = D_src * Dscale
    period = 2.0 * math.pi / Om
    n_burn = int(round(burn_periods * period / dt))
    n_meas = int(round(meas_periods * period / dt))
    T = n_meas * dt
    m0 = int(round(Om * T / (2.0 * math.pi)))
    bins = [m for m in range(m0 - ANNULUS, m0 + ANNULUS + 1) if m != m0]
    omegas = [Om] + [2.0 * math.pi * m / T for m in bins]
    t0 = time.time()
    out = double_well_fourier_driver(p, D, n_trials=n_trials, n_burn=n_burn, n_meas=n_meas, omegas=omegas,
                                     rng=np.random.default_rng(seed))
    X = out["X"]
    c = (2.0 / T) * X[0]
    n = n_trials
    M = c.mean()
    s2 = float(np.sum(np.abs(c - M) ** 2) / (n - 1))
    A = p.signal_amplitude
    amp2 = float(abs(M) ** 2 - s2 / n)
    eta = amp2 / A ** 2
    Xb = X[1:]
    resid = Xb - Xb.mean(axis=1, keepdims=True)
    SN = float(np.mean(np.sum(np.abs(resid) ** 2, axis=1) / (n - 1)) / T)
    snr = math.pi * amp2 / SN
    return {
        "frame": frame, "dt": dt, "D_src": D_src, "D_model": D, "seed": seed, "n_trials": n,
        "n_burn": n_burn, "n_meas": n_meas, "T": T,
        "M_re": float(M.real), "M_im": float(M.imag), "s2": s2, "eta_hat": eta,
        "S_N": SN, "snr": snr, "annulus_bins": bins, "bin_omega": m0,
        "frac_x_positive_end": float(np.mean(out["x_final"] > 0)),
        "runtime_s": round(time.time() - t0, 2),
    }


def load_raw():
    recs = {}
    if os.path.exists(RAW):
        with open(RAW) as fh:
            for line in fh:
                line = line.strip()
                if line:
                    d = json.loads(line)
                    recs[(d["frame"], d["dt"], d["D_src"], d["seed"])] = d
    return recs


def ensure(jobs, commit):
    recs = load_raw()
    todo = [j for j in jobs if j not in recs]
    print(f"[p26] {len(todo)} of {len(jobs)} runs to do", flush=True)
    os.makedirs(os.path.dirname(RAW), exist_ok=True)
    with open(RAW, "a") as fh:
        for k, (frame, dt, D, seed) in enumerate(todo):
            r = run_one(frame, dt, D, seed)
            r["commit"] = commit
            fh.write(json.dumps(r) + "\n")
            fh.flush()
            print(f"[p26] {frame} dt={dt} D={D} seed={seed} eta={r['eta_hat']:.3f} ({r['runtime_s']} s) "
                  f"[{k + 1}/{len(todo)}]", flush=True)
    return load_raw()


def seed_stats(recs, frame, dt, D, key="eta_hat", factor=1.0):
    v = np.array([recs[(frame, dt, D, s)][key] for s in SEEDS]) * factor
    return float(v.mean()), float(v.std(ddof=1) / math.sqrt(len(v))), v


def dt_rule(recs):
    rows, trig = [], False
    for D in D_CHECK:
        e1, s1, _ = seed_stats(recs, "unit", 0.01, D)
        e2, s2, _ = seed_stats(recs, "unit", 0.005, D)
        rel = abs(e2 - e1) / e1
        trig |= rel > 0.02
        rows.append({"D": D, "eta_dt0.01": e1, "se_dt0.01": s1, "eta_dt0.005": e2, "se_dt0.005": s2,
                     "rel_diff": rel, "exceeds_0.02": bool(rel > 0.02)})
    return rows, bool(trig)


def quad_vertex(Ds, etas):
    a, b, c = np.polyfit(Ds, etas, 2)
    if a >= 0:
        return float("nan"), float("nan"), (a, b, c)
    return float(-b / (2 * a)), float(c - b * b / (4 * a)), (a, b, c)


def analyze(recs, commit, t_wall):
    rows_dt, trig = dt_rule(recs)
    gate_dt = 0.005 if trig else 0.01
    # curve (main)
    curve = []
    for D in D_GRID:
        e, se, v = seed_stats(recs, "unit", 0.01, D)
        sn, snse, _ = seed_stats(recs, "unit", 0.01, D, key="snr")
        curve.append({"D": D, "eta_mean": e, "eta_se": se, "eta_per_seed": v.tolist(), "snr_mean": sn, "snr_se": snse})
    curve_gate = [c for c in curve] if gate_dt == 0.01 else []
    if gate_dt == 0.005:
        for D in D_WINDOW + [0.20, 0.30]:
            e, se, v = seed_stats(recs, "unit", 0.005, D)
            curve_gate.append({"D": D, "eta_mean": e, "eta_se": se, "eta_per_seed": v.tolist()})
    by_D = {c["D"]: c for c in curve_gate}
    # peak
    Ds = np.array(D_WINDOW)
    E = np.array([by_D[D]["eta_per_seed"] for D in D_WINDOW])        # (7, 10)
    Dmax, emax, coef = quad_vertex(Ds, E.mean(axis=1))
    rng = np.random.default_rng(BOOT_SEED)
    bD, bE = [], []
    for _ in range(BOOT_N):
        idx = rng.integers(0, len(SEEDS), size=len(SEEDS))
        d_, e_, _ = quad_vertex(Ds, E[:, idx].mean(axis=1))
        bD.append(d_)
        bE.append(e_)
    bD, bE = np.array(bD), np.array(bE)
    peak = {
        "D_max": Dmax, "eta_max": emax, "quadratic_coefficients": [float(x) for x in coef],
        "bootstrap_n": BOOT_N, "bootstrap_seed": BOOT_SEED,
        "D_max_se": float(np.nanstd(bD, ddof=1)), "eta_max_se": float(np.nanstd(bE, ddof=1)),
        "D_max_boot_2.5_97.5": np.nanpercentile(bD, [2.5, 97.5]).tolist(),
        "eta_max_boot_2.5_97.5": np.nanpercentile(bE, [2.5, 97.5]).tolist(),
        "bootstrap_nonconcave_fits": int(np.isnan(bD).sum()),
        "dt_used": gate_dt,
    }
    e20, e30 = by_D[0.20]["eta_mean"], by_D[0.30]["eta_mean"]
    conds = {
        "1_eta_max_in_[11.62,12.84]": bool(12.23 * 0.95 <= emax <= 12.23 * 1.05),
        "2_|D_max-0.150|<=0.015": bool(abs(Dmax - 0.150) <= 0.015),
        "3_eta(0.20)_in_10.8(1+-0.08)": bool(10.8 * 0.92 <= e20 <= 10.8 * 1.08),
        "4_eta(0.30)_in_6.5(1+-0.08)": bool(6.5 * 0.92 <= e30 <= 6.5 * 1.08),
    }
    # frame check
    frame_rows, frame_ok = [], True
    for D in D_CHECK:
        e1, s1, _ = seed_stats(recs, "unit", 0.01, D)
        e2, s2, _ = seed_stats(recs, "canonical", 0.0025, D, factor=16.0)
        ok = abs(e2 - e1) <= 2 * math.sqrt(s1 ** 2 + s2 ** 2)
        frame_ok &= ok
        frame_rows.append({"D": D, "eta_unit": e1, "se_unit": s1, "16_eta_canonical": e2, "se_canonical_x16": s2,
                           "abs_diff": abs(e2 - e1), "within_2_combined_SE": bool(ok)})
    passed = all(conds.values()) if frame_ok else None
    # SNR maximum (reported)
    snr = np.array([c["snr_mean"] for c in curve])
    locmax = [D_GRID[i] for i in range(1, len(D_GRID) - 1) if snr[i] > snr[i - 1] and snr[i] > snr[i + 1]]
    snr_vertex = None
    if locmax:
        i = D_GRID.index(locmax[0])
        dv, ev, _ = quad_vertex(np.array(D_GRID[i - 1:i + 2]), snr[i - 1:i + 2])
        snr_vertex = {"D": dv, "snr": ev}
    runtime = {}
    for (fr, dt, D, s), r in recs.items():
        k = f"{fr}_dt{dt}"
        runtime[k] = runtime.get(k, 0.0) + r["runtime_s"]
    out = {
        "pattern": "P26",
        "source": "Gammaitoni et al. 1998 RMP 70:223, Fig. 7 (Floquet, A = 0.1, Omega = 0.1), digitized in group_C.md",
        "spec": "docs/validation_rebuild/repro_v2/group_C.md (P26 sections 3-4), frozen at 39688e3",
        "targets": TARGET,
        "protocol": {"a": 1, "b": 1, "A": A_SRC, "Omega": OMEGA_SRC, "dt": 0.01, "D_grid": D_GRID, "seeds": SEEDS,
                     "trials_per_seed": N_TRIALS, "burn_in_periods": BURN_PERIODS, "measurement_periods": MEAS_PERIODS},
        "implementation_choices": IMPLEMENTATION_CHOICES,
        "eta_curve_dt0.01": curve,
        "dt_check": {"rows": rows_dt, "triggered": trig, "gated_values_from_dt": gate_dt},
        "gated_curve": curve_gate if gate_dt == 0.005 else "same as eta_curve_dt0.01",
        "peak": peak,
        "eta_0.20": {"mean": e20, "se": by_D[0.20]["eta_se"]},
        "eta_0.30": {"mean": e30, "se": by_D[0.30]["eta_se"]},
        "conditions": conds,
        "frame_check": {"rows": frame_rows, "passed": bool(frame_ok)},
        "PASS": passed,
        "reported_only": {
            "snr_curve": [{"D": c["D"], "snr_mean": c["snr_mean"], "snr_se": c["snr_se"]} for c in curve],
            "snr_local_maxima_grid": locmax,
            "snr_local_max_vertex": snr_vertex,
            "two_state_prediction_DeltaV_over_2": 0.125,
        },
        "runtime_cpu_seconds": {k: round(v, 1) for k, v in runtime.items()},
        "wall_seconds_this_invocation": round(t_wall, 1),
        "code_commits_of_runs": sorted({r["commit"] for r in recs.values()}),
        "commit_at_analysis": commit,
        "numpy_version": np.__version__,
    }
    os.makedirs(os.path.dirname(OUT_JSON), exist_ok=True)
    with open(OUT_JSON, "w") as fh:
        json.dump(out, fh, indent=1)
    print(json.dumps({"PASS": passed, "conditions": conds, "peak": {k: peak[k] for k in ("D_max", "eta_max", "D_max_se", "eta_max_se")},
                      "eta_0.20": e20, "eta_0.30": e30, "dt_check": rows_dt, "frame_ok": frame_ok,
                      "snr_locmax": locmax}, indent=1))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("cmd", choices=["all", "analyze", "smoke"])
    args = ap.parse_args()
    commit = git_commit()
    t0 = time.time()
    if args.cmd == "smoke":
        r = run_one("unit", 0.01, 0.15, 999, n_trials=20, burn_periods=0.5, meas_periods=1)
        print(json.dumps(r, indent=1))
        return
    if args.cmd == "all":
        ensure([("unit", 0.01, D, s) for D in D_GRID for s in SEEDS], commit)
        recs = ensure([("unit", 0.005, D, s) for D in D_CHECK for s in SEEDS], commit)
        rows, trig = dt_rule(recs)
        print(f"[p26] dt rule: {rows} triggered={trig}", flush=True)
        if trig:
            ensure([("unit", 0.005, D, s) for D in D_WINDOW + [0.20] for s in SEEDS], commit)
        ensure([("canonical", 0.0025, D, s) for D in D_CHECK for s in SEEDS], commit)
    analyze(load_raw(), commit, time.time() - t0)


if __name__ == "__main__":
    main()
