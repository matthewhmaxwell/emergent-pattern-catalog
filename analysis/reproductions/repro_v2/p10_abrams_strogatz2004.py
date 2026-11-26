"""repro_v2 P10: Abrams & Strogatz 2004 chimera, profile and existence boundary (group_B.md, section P10).

All times below are paper time units (t_paper = t_ours * N / 2pi). dt_ours = dt_paper * 2pi / N (option
``dt_paper`` of KuramotoNonlocal); IC = init_mode "as2004" (6 exp(-30 x^2) U[-1/2, 1/2]).

Stage A (T10.1): N=256, A=0.995, beta=0.18, seeds 0-9, dt_paper=0.025, T=5000, window [2500, 5000].
  r_bar: time mean of |<e^{i theta}>| sampled every 1.0; Omega_i = (theta_i(5000)-theta_i(2500))/2500;
  locked set = modal 0.001-wide bin of Omega plus oscillators within 0.002 of it; Delta = -median(Omega_locked);
  f_drift = 1 - n_locked/N; R(x) = |Z0 + A (Zc cos x + Zs sin x)| on a 1024-point grid every 10.0,
  R_min/R_max averaged over snapshots. Chimera: r_bar < 0.95 and 0.3 < f_drift < 0.95.
  PASS: medians over chimera seeds within r 0.02, Delta 0.010, f_drift 0.03, R_min 0.03, R_max 0.03 of the
  analytic stable branch.
Stage B (T10.2): from each Stage A chimera: beta 0.20, 0.22, 0.24 (hold 500 each), then 0.2400..0.2800 step
  0.0025 (hold 1000 each). Collapse: r > 0.95 over the last 200 of a hold and f_drift < 0.02.
  PASS: median beta_collapse in [0.2500, 0.2675] and no seed holds a chimera at beta >= 0.2700.
Stage C (T10.3, reported): Stage A at beta=0.05.
dt check: two seeds at dt_paper=0.0125; Stage A estimators must change by < half their tolerance.
Secondary: Stage A + B at N=128 (10 seeds).
"""
from __future__ import annotations

import time
from typing import Any

import numpy as np

from common import RAW, git_info, read_json  # noqa: E402

from epc.models.kuramoto_nonlocal import KuramotoNonlocal

A_KERNEL = 0.995
DTP = 0.025
DTP_CHECK = 0.0125
T_A, WIN_A = 5000, 2500
SEEDS = list(range(10))
TARGETS = {  # frozen, group_B.md P10 table (stable branch of Eqs. 12-13, A=0.995)
    0.18: {"r_bar": 0.6287, "Delta": 0.7320, "f_drift": 0.7008, "R_min": 0.4732, "R_max": 0.8068},
    0.05: {"r_bar": 0.5760, "Delta": 0.7160, "f_drift": 0.8526, "R_min": 0.4214, "R_max": 0.7326},
}
TOLS = {"r_bar": 0.02, "Delta": 0.010, "f_drift": 0.03, "R_min": 0.03, "R_max": 0.03}
BETA_SN = 0.2635
RAMP = [(0.20, 500), (0.22, 500), (0.24, 500)] + [(round(0.24 + 0.0025 * k, 4), 1000) for k in range(17)]
COLLAPSE_LAST = 200
BIN_W, NEIGH = 0.001, 0.002
NGRID = 1024
T10_2_BAND = (0.2500, 0.2675)
T10_2_NO_CHIMERA_FROM = 0.2700
FALLBACK_SEED0 = 100
FALLBACK_MAX_ATTEMPTS = 20

RAWDIR = RAW / "P10"


# ---------------------------------------------------------------- integration
class Runner:
    """Drives KuramotoNonlocal.run() in chunks of whole paper-time units; theta stays unwrapped."""

    def __init__(self, N: int, beta: float, seed: int, dt_paper: float, init_mode: str = "as2004",
                 theta: np.ndarray | None = None):
        dt_ours = dt_paper * 2 * np.pi / N
        self.spu = int(round(1.0 / dt_paper))  # steps per paper time unit
        self.m = KuramotoNonlocal(N=N, A=A_KERNEL, beta=beta, seed=seed, init_mode=init_mode,
                                  dt_paper=dt_paper, record_dt=self.spu * dt_ours)
        assert abs(self.m.params.dt - dt_ours) < 1e-18
        if theta is not None:
            self.m.theta = np.array(theta, dtype=float, copy=True)
        self.steps = 0

    def advance(self, units: int) -> list[float]:
        """Integrate `units` paper-time units; return r sampled at the end of each unit."""
        hist = self.m.run(n_frames=units)
        assert hist[-1]["step"] == units * self.spu
        self.steps += units * self.spu
        return [h["r"] for h in hist]

    @property
    def theta(self) -> np.ndarray:
        return self.m.theta.copy()


def locked_set(Om: np.ndarray) -> dict[str, Any]:
    b = np.floor(Om / BIN_W).astype(np.int64)
    vals, counts = np.unique(b, return_counts=True)
    k = int(np.argmax(counts))
    lo, hi = vals[k] * BIN_W, (vals[k] + 1) * BIN_W
    locked = (Om >= lo - NEIGH) & (Om < hi + NEIGH)
    return {"mask": locked, "modal_bin": [float(lo), float(hi)], "modal_count": int(counts[k]),
            "n_tied_modal_bins": int(np.sum(counts == counts[k]))}


def freq_stats(theta_a: np.ndarray, theta_b: np.ndarray, span: float) -> dict[str, Any]:
    Om = (theta_b - theta_a) / span
    ls = locked_set(Om)
    mask = ls.pop("mask")
    N = len(Om)
    return {"Delta": float(-np.median(Om[mask])), "f_drift": float(1 - mask.sum() / N),
            "n_locked": int(mask.sum()), **ls, "Omega": Om}


def R_profile(theta: np.ndarray, x: np.ndarray) -> tuple[float, float, float]:
    """Min, max of R(x) = |Z0 + A (Zc cos x + Zs sin x)| on a 1024-point grid; also the argmin location."""
    e = np.exp(1j * theta)
    Z0, Zc, Zs = e.mean(), (np.cos(x) * e).mean(), (np.sin(x) * e).mean()
    xg = 2 * np.pi * np.arange(NGRID) / NGRID
    R = np.abs(Z0 + A_KERNEL * (Zc * np.cos(xg) + Zs * np.sin(xg)))
    return float(R.min()), float(R.max()), float(xg[np.argmin(R)])


def stage_a(N: int, beta: float, seed: int, dt_paper: float = DTP, init_mode: str = "as2004",
            theta0: np.ndarray | None = None, T: int | None = None, win0: int | None = None,
            save_key: str | None = None) -> dict[str, Any]:
    T = T_A if T is None else T
    win0 = WIN_A if win0 is None else win0
    t0, c0 = time.time(), time.process_time()
    run = Runner(N, beta, seed, dt_paper, init_mode=init_mode, theta=theta0)
    x = run.m.x
    r_coarse, t = [], 0
    while t < win0:
        rs = run.advance(10); t += 10
        r_coarse.append(rs[-1])
    th_w0 = run.theta
    r_win, snaps = [], []
    while t < T:
        rs = run.advance(10); t += 10
        r_win.extend(rs)
        r_coarse.append(rs[-1])
        snaps.append(R_profile(run.m.theta, x))
    th_end = run.theta
    fs = freq_stats(th_w0, th_end, T - win0)
    Om = fs.pop("Omega")
    r_win = np.array(r_win)
    snaps = np.array(snaps)
    est = {"r_bar": float(r_win.mean()), "Delta": fs["Delta"], "f_drift": fs["f_drift"],
           "R_min": float(snaps[:, 0].mean()), "R_max": float(snaps[:, 1].mean())}
    chim = bool(est["r_bar"] < 0.95 and 0.3 < est["f_drift"] < 0.95)
    out = {"N": N, "beta": beta, "seed": seed, "A": A_KERNEL, "dt_paper": dt_paper,
           "dt_ours": run.m.params.dt, "init_mode": init_mode, "T_paper": T, "window": [win0, T],
           "estimators": est, "chimera": chim, "freq": fs,
           "r_window_sd": float(r_win.std()), "r_window_min": float(r_win.min()), "r_window_max": float(r_win.max()),
           "R_min_snap_sd": float(snaps[:, 0].std()), "R_max_snap_sd": float(snaps[:, 1].std()),
           "incoherent_center_track_every10": snaps[:, 2].tolist(),
           "r_every10_full_run": r_coarse, "steps": run.steps,
           "Omega_sorted_by_position": Om.tolist()}
    if save_key:
        RAWDIR.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(RAWDIR / f"{save_key}_theta.npz", theta_end=th_end, theta_w0=th_w0)
    out["wall_s"] = time.time() - t0
    out["cpu_s"] = time.process_time() - c0
    return out


def hold(N: int, beta: float, seed: int, theta: np.ndarray, units: int, dt_paper: float = DTP,
         last: int = COLLAPSE_LAST) -> tuple[dict[str, Any], np.ndarray]:
    run = Runner(N, beta, seed, dt_paper, theta=theta)
    r_all = run.advance(units - last)
    th_a = run.theta
    r_last = np.array(run.advance(last))
    th_b = run.theta
    fs = freq_stats(th_a, th_b, last)
    fs.pop("Omega")
    rec = {"beta": beta, "hold": units, "r_mean_last200": float(r_last.mean()),
           "r_min_last200": float(r_last.min()), "f_drift_last200": fs["f_drift"], "Delta_last200": fs["Delta"],
           "n_locked_last200": fs["n_locked"],
           "collapse": bool(r_last.min() > 0.95 and fs["f_drift"] < 0.02),
           "chimera": bool(r_last.mean() < 0.95 and 0.3 < fs["f_drift"] < 0.95),
           "r_every10": (r_all + r_last.tolist())[9::10]}
    return rec, th_b


def stage_b(N: int, seed: int, theta_start: np.ndarray, dt_paper: float = DTP) -> dict[str, Any]:
    """Ramp per RAMP. After the first collapse one further hold is run as a confirmation and the ramp stops
    (pre-declared reading: full sync is absorbing; beta_collapse only needs the first collapse)."""
    t0, c0 = time.time(), time.process_time()
    theta = theta_start
    holds: list[dict[str, Any]] = []
    i_collapse = None
    for i, (beta, units) in enumerate(RAMP):
        rec, theta = hold(N, beta, seed, theta, units, dt_paper)
        rec["ramp_index"] = i
        holds.append(rec)
        if i_collapse is None:
            if rec["collapse"]:
                i_collapse = i
            continue
        break  # this was the confirmation hold
    beta_collapse = RAMP[i_collapse][0] if i_collapse is not None else float("inf")
    confirm = holds[i_collapse + 1]["collapse"] if (i_collapse is not None and len(holds) > i_collapse + 1) else None
    ge27 = [h for h in holds if h["beta"] >= T10_2_NO_CHIMERA_FROM - 1e-9]
    return {"N": N, "seed": seed, "dt_paper": dt_paper, "ramp": RAMP, "holds": holds,
            "beta_collapse": beta_collapse, "collapsed_within_ramp": i_collapse is not None,
            "collapse_ramp_index": i_collapse,
            "stopped_after_confirmation_hold": i_collapse is not None and len(holds) < len(RAMP),
            "confirmation_hold_still_collapsed": confirm,
            "chimera_at_hold_with_beta_ge_0.27": any(h["chimera"] for h in ge27),
            "wall_s": time.time() - t0, "cpu_s": time.process_time() - c0}


def fallback_realization(k: int, N: int = 256) -> dict[str, Any]:
    """Pre-declared fallback (only run if < 5 of 10 paper-IC seeds are chimeras): asymmetric-Gaussian IC
    (asym_amp=2.0, asym_sigma=0.5, seed 100+k) at beta=0.05 held 500, ramp 0.06..0.18 by 0.01 (500 each),
    then the Stage A protocol at beta=0.18 (T=5000, window [2500, 5000]) from the ramped state."""
    t0, c0 = time.time(), time.process_time()
    seed = FALLBACK_SEED0 + k
    run = Runner(N, 0.05, seed, DTP, init_mode="asymmetric_gaussian")
    run.advance(500)
    theta = run.theta
    ramp_recs = []
    for j in range(1, 14):
        beta = round(0.05 + 0.01 * j, 2)
        rec, theta = hold(N, beta, seed, theta, 500)
        ramp_recs.append(rec)
    a = stage_a(N, 0.18, seed, theta0=theta, save_key=f"P10_F_N{N}_k{k}")
    a["origin"] = "fallback_ramp"
    a["ramp"] = ramp_recs
    a["wall_s"] = time.time() - t0
    a["cpu_s"] = time.process_time() - c0
    return a


# ---------------------------------------------------------------- tasks
def a_id(N: int, beta: float, seed: int, dtp: float = DTP) -> str:
    return f"P10_A_N{N}_b{beta:g}_s{seed}_dtp{dtp:g}"


def b_id(N: int, seed: int) -> str:
    return f"P10_B_N{N}_s{seed}"


def _res(task_id: str) -> dict[str, Any] | None:
    f = RAWDIR / f"{task_id}.json"
    return read_json(f)["result"] if f.exists() else None


def stage_a_results(N: int, beta: float = 0.18, dtp: float = DTP) -> list[dict[str, Any] | None]:
    return [_res(a_id(N, beta, s, dtp)) for s in SEEDS]


def tasks_stage_a(N: int, beta: float = 0.18) -> list[dict[str, Any]]:
    return [dict(id=a_id(N, beta, s), kind="A", N=N, beta=beta, seed=s, dtp=DTP) for s in SEEDS]


def tasks_stage_b(N: int) -> tuple[list[dict[str, Any]], int]:
    t, blocked = [], 0
    for s in SEEDS:
        a = _res(a_id(N, 0.18, s))
        if a is None:
            blocked += 1
        elif a["chimera"]:
            t.append(dict(id=b_id(N, s), kind="B", N=N, seed=s, src=a_id(N, 0.18, s)))
    return t, blocked


def chimera_seeds(N: int) -> list[int] | None:
    ra = stage_a_results(N)
    if any(r is None for r in ra):
        return None
    return [r["seed"] for r in ra if r["chimera"]]


def tasks_dtcheck() -> tuple[list[dict[str, Any]], int]:
    cs = chimera_seeds(256)
    if cs is None:
        return [], 1
    return [dict(id=a_id(256, 0.18, s, DTP_CHECK), kind="A", N=256, beta=0.18, seed=s, dtp=DTP_CHECK)
            for s in cs[:2]], 0


def fallback_need() -> int | None:
    cs = chimera_seeds(256)
    if cs is None:
        return None
    return 0 if len(cs) >= 5 else 10 - len(cs)


def fallback_used_ks() -> list[int]:
    """Lowest-k fallback attempts that ended in a chimera, up to the number needed."""
    need = fallback_need() or 0
    ks = []
    for k in range(FALLBACK_MAX_ATTEMPTS):
        r = _res(f"P10_F_N256_k{k}")
        if r is not None and r["chimera"] and len(ks) < need:
            ks.append(k)
    return ks


def tasks_fallback() -> tuple[list[dict[str, Any]], int]:
    need = fallback_need()
    if need is None:
        return [], 1
    t: list[dict[str, Any]] = []
    have = pending = 0
    for k in range(FALLBACK_MAX_ATTEMPTS):
        if have >= need:
            break
        r = _res(f"P10_F_N256_k{k}")
        if r is None:
            if have + pending < need:
                t.append(dict(id=f"P10_F_N256_k{k}", kind="F", k=k))
                pending += 1
        elif r["chimera"]:
            have += 1
            if _res(f"P10_BF_N256_k{k}") is None:
                t.append(dict(id=f"P10_BF_N256_k{k}", kind="BF", k=k))
    return t, 0


def run_task(task: dict[str, Any]) -> dict[str, Any]:
    k = task["kind"]
    if k == "A":
        return stage_a(task["N"], task["beta"], task["seed"], dt_paper=task["dtp"], save_key=task["id"])
    if k == "B":
        th = np.load(RAWDIR / f"{task['src']}_theta.npz")["theta_end"]
        return stage_b(task["N"], task["seed"], th)
    if k == "F":
        return fallback_realization(task["k"])
    if k == "BF":
        th = np.load(RAWDIR / f"P10_F_N256_k{task['k']}_theta.npz")["theta_end"]
        return stage_b(256, FALLBACK_SEED0 + task["k"], th)
    raise ValueError(k)


# ---------------------------------------------------------------- aggregation
def _profile_rule(results: list[dict[str, Any]], beta: float) -> dict[str, Any]:
    chim = [r for r in results if r["chimera"]]
    out: dict[str, Any] = {"n_runs": len(results), "n_chimera": len(chim),
                           "basin_fraction": len(chim) / len(results) if results else None,
                           "per_seed": [{"seed": r["seed"], "origin": r.get("origin", "paper_IC"),
                                         "chimera": r["chimera"], **r["estimators"],
                                         "n_locked": r["freq"]["n_locked"],
                                         "modal_bin": r["freq"]["modal_bin"]} for r in results],
                           "targets": TARGETS[beta], "tolerances": TOLS}
    if chim:
        med = {k: float(np.median([r["estimators"][k] for r in chim])) for k in TOLS}
        out["median_over_chimera_seeds"] = med
        out["diff"] = {k: med[k] - TARGETS[beta][k] for k in TOLS}
        out["within"] = {k: abs(out["diff"][k]) <= TOLS[k] for k in TOLS}
        out["all_five_within"] = all(out["within"].values())
    else:
        out["all_five_within"] = False
    return out


def aggregate() -> dict[str, Any]:
    res: dict[str, Any] = {}
    allraw = [read_json(f) for f in sorted(RAWDIR.glob("P10_*.json"))]
    res["n_runs"] = len(allraw)
    res["runtime_cpu_h_total"] = sum(r["result"]["cpu_s"] for r in allraw) / 3600
    res["runtime_wall_h_total_summed_over_runs"] = sum(r["result"]["wall_s"] for r in allraw) / 3600
    res["commits"] = sorted({r["git"]["commit"] for r in allraw if r.get("git")})
    res["any_code_dirty"] = any(r["git"].get("code_dirty") for r in allraw if r.get("git"))

    ra = [r for r in stage_a_results(256) if r is not None]
    fb_all = [_res(f"P10_F_N256_k{k}") for k in range(FALLBACK_MAX_ATTEMPTS)]
    fb_all = [r for r in fb_all if r is not None]
    used_ks = fallback_used_ks()
    fb = [_res(f"P10_F_N256_k{k}") for k in used_ks]
    a_all = ra + fb
    res["fallback"] = {"needed_realizations": fallback_need(), "attempts_run": len(fb_all),
                       "attempts_chimera": [r["seed"] for r in fb_all if r["chimera"]],
                       "used_seeds": [100 + k for k in used_ks]}
    complete_a = len(ra) == len(SEEDS)
    t101 = _profile_rule(a_all, 0.18)
    t101["paper_IC_basin_fraction"] = (sum(r["chimera"] for r in ra) / len(ra)) if ra else None
    t101["fallback_used"] = bool(fb)
    t101["complete"] = complete_a
    t101["outcome"] = ("PASS" if t101["all_five_within"] else "FAIL") if complete_a else "INCOMPLETE"
    t101["rule"] = ("medians over chimera seeds: |r_bar-0.6287|<=0.02, |Delta-0.7320|<=0.010, "
                    "|f_drift-0.7008|<=0.03, |R_min-0.4732|<=0.03, |R_max-0.8068|<=0.03")
    res["T10_1"] = t101

    # dt check
    rows = []
    cs = chimera_seeds(256)
    for s in (cs or [])[:2]:
        r1, r2 = _res(a_id(256, 0.18, s)), _res(a_id(256, 0.18, s, DTP_CHECK))
        if r2 is None:
            continue
        d = {k: r2["estimators"][k] - r1["estimators"][k] for k in TOLS}
        rows.append({"seed": s, "chimera_at_dt0.0125": r2["chimera"], "est_dt0.025": r1["estimators"],
                     "est_dt0.0125": r2["estimators"], "change": d,
                     "lt_half_tol": {k: abs(d[k]) < TOLS[k] / 2 for k in TOLS}})
    res["dt_check"] = {"rule": "two seeds at dt_paper=0.0125: every Stage A estimator changes by < half its "
                               "tolerance (seeds = the two lowest-numbered Stage A chimera seeds)",
                       "rows": rows,
                       "passed": (all(all(r["lt_half_tol"].values()) and r["chimera_at_dt0.0125"] for r in rows)
                                  if len(rows) == 2 else None)}

    def boundary(N: int, with_fallback: bool) -> dict[str, Any]:
        bs = []
        for s in SEEDS:
            r = _res(b_id(N, s))
            if r is not None:
                r["origin"] = "paper_IC"
                bs.append(r)
        if with_fallback:
            for k in fallback_used_ks():
                r = _res(f"P10_BF_N{N}_k{k}")
                if r is not None:
                    r["origin"] = "fallback_ramp"
                    bs.append(r)
        csN = chimera_seeds(N)
        expected = (len(csN) if csN is not None else None)
        out = {"per_seed": [{"seed": b["seed"], "origin": b["origin"], "beta_collapse": b["beta_collapse"],
                             "chimera_at_hold_with_beta_ge_0.27": b["chimera_at_hold_with_beta_ge_0.27"],
                             "confirmation_hold_still_collapsed": b["confirmation_hold_still_collapsed"],
                             "holds": [{k: h[k] for k in ("beta", "r_mean_last200", "r_min_last200",
                                                          "f_drift_last200", "collapse", "chimera")}
                                       for h in b["holds"]]} for b in bs],
               "n_ramps": len(bs), "n_paper_IC_chimera_seeds": expected}
        if bs:
            bc = np.array([b["beta_collapse"] for b in bs], dtype=float)
            med = float(np.median(bc))
            out["median_beta_collapse"] = med
            out["beta_collapse_values"] = bc.tolist()
            out["median_in_band"] = bool(T10_2_BAND[0] <= med <= T10_2_BAND[1])
            no_chim = not any(b["chimera_at_hold_with_beta_ge_0.27"] for b in bs)
            all_by_027 = all(b["beta_collapse"] <= T10_2_NO_CHIMERA_FROM + 1e-9 for b in bs)
            out["no_hold_at_beta_ge_0.27_classified_chimera"] = no_chim
            out["every_seed_collapsed_by_beta_0.27"] = all_by_027
            out["no_chimera_at_beta_ge_0.27"] = no_chim and all_by_027
        return out

    b256 = boundary(256, True)
    complete_b = complete_a and b256["n_ramps"] >= sum(r["chimera"] for r in a_all) and b256["n_ramps"] > 0
    b256["complete"] = complete_b
    b256["rule"] = "median beta_collapse in [0.2500, 0.2675] and no seed holds a chimera at beta >= 0.2700 (N=256)"
    b256["beta_SN_continuum"] = BETA_SN
    b256["outcome"] = (("PASS" if b256.get("median_in_band") and b256.get("no_chimera_at_beta_ge_0.27") else "FAIL")
                       if complete_b else "INCOMPLETE")
    res["T10_2"] = b256

    rc = [_res(a_id(256, 0.05, s)) for s in SEEDS]
    rc = [r for r in rc if r is not None]
    t103 = _profile_rule(rc, 0.05)
    t103["complete"] = len(rc) == len(SEEDS)
    t103["outcome"] = ("PASS" if t103["all_five_within"] else "FAIL") if t103["complete"] else "INCOMPLETE"
    t103["rule"] = "as T10.1 with the beta=0.05 analytic row (reported, not gating)"
    res["T10_3"] = t103

    r128 = [r for r in stage_a_results(128) if r is not None]
    sec = {"stage_A": _profile_rule(r128, 0.18) if r128 else None, "stage_B": boundary(128, False),
           "note": "secondary, reported: finite-N trend of beta_collapse (expected to approach 0.2635 from below)"}
    res["secondary_N128"] = sec
    res["optional_N512"] = "not run"
    return res
