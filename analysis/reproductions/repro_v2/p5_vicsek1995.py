"""repro_v2 P5: Vicsek et al. 1995 Fig. 2a, v_a(eta) at fixed density (group_B.md, section P5).

Protocol (frozen, group_B.md P5 sections 3-4):
  sizes N=400/L=10 (primary), N=100/L=5, N=40/L=3.1; optional N=4000/L=31.6 at eta in {2.0, 2.5, 3.0}
  eta in {0.5, ..., 4.5}; v=0.03, r=1, dt=1; random IC
  primary runs: position_update="old_velocity" (Vicsek Eq. 1); N=400 also with the catalog order
  ("new_velocity") as a sensitivity result
  seeds 1000-1009 (1000-1004 for N=4000)
  burn-in 2e4 steps; measurement 8e4 steps, phi sampled every 10 steps; if first- and second-half
  means differ by > 0.02, the run's measurement is extended once to 2e5 steps and the extended value used
  estimator: per-seed time mean of phi; v_sim = mean over seeds; SE = sd/sqrt(n); 95% t-interval
  PASS (primary): |v_sim - v_pub| <= 0.04 at all 9 eta for N=400.

Measurement streams phi from step() (run() would store every step).
"""
from __future__ import annotations

import time
from pathlib import Path
from typing import Any

import numpy as np

from common import RAW, OUT, SUPPORT, git_info, read_json, write_json, t_interval  # noqa: E402

from epc.models.vicsek import VicsekModel, polarization

ETAS = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5]
BOX = {400: 10.0, 100: 5.0, 40: 3.1, 4000: 31.6}
TARGETS = {  # frozen, group_B.md P5 section 1 table
    400: [0.981, 0.920, 0.830, 0.709, 0.566, 0.380, 0.133, 0.081, 0.062],
    100: [0.983, 0.932, 0.848, 0.741, 0.621, 0.461, 0.278, 0.157, 0.121],
    40: [0.985, 0.936, 0.863, 0.767, 0.654, 0.515, 0.380, 0.260, 0.199],
}
ETAS_4000 = [2.0, 2.5, 3.0]
TARGETS_4000 = {2.0: 0.645, 2.5: 0.463, 3.0: 0.142}
SEEDS = list(range(1000, 1010))
SEEDS_4000 = list(range(1000, 1005))
BURN = 20_000
MEAS = 80_000
MEAS_EXT = 200_000
SAMPLE_EVERY = 10
HALF_TOL = 0.02
TOL = 0.04
TOL_4000 = 0.06
SPEED, RADIUS, DT = 0.03, 1.0, 1.0

RAWDIR = RAW / "P5"


def task_id(N: int, eta: float, seed: int, update: str) -> str:
    return f"P5_N{N}_eta{eta:.1f}_s{seed}_{update}"


def all_tasks(include_optional: bool = True) -> list[dict[str, Any]]:
    t = []
    for N in (400, 100, 40):
        for update in (["old_velocity", "new_velocity"] if N == 400 else ["old_velocity"]):
            for eta in ETAS:
                for s in SEEDS:
                    t.append(dict(id=task_id(N, eta, s, update), N=N, eta=eta, seed=s, update=update))
    if include_optional:
        for eta in ETAS_4000:
            for s in SEEDS_4000:
                t.append(dict(id=task_id(4000, eta, s, "old_velocity"), N=4000, eta=eta, seed=s,
                              update="old_velocity"))
    return t


def run_cell(N: int, eta: float, seed: int, update: str,
             burn: int | None = None, meas: int | None = None, meas_ext: int | None = None,
             save_series: bool = True) -> dict[str, Any]:
    burn = BURN if burn is None else burn
    meas = MEAS if meas is None else meas
    meas_ext = MEAS_EXT if meas_ext is None else meas_ext
    t0, c0 = time.time(), time.process_time()
    m = VicsekModel(n_particles=N, box_size=BOX[N], speed=SPEED, noise=eta,
                    interaction_radius=RADIUS, dt=DT, init_mode="random", seed=seed,
                    position_update=update, wrap_guard=True)
    m.setup()
    for _ in range(burn):
        m.step()
    phis: list[float] = []

    def measure(n_steps: int) -> None:
        for _ in range(n_steps // SAMPLE_EVERY):
            for _ in range(SAMPLE_EVERY - 1):
                m.step()
            st = m.step()
            phis.append(float(polarization(st)))

    measure(meas)
    n0 = len(phis)
    h1, h2 = float(np.mean(phis[: n0 // 2])), float(np.mean(phis[n0 // 2:]))
    extended = abs(h1 - h2) > HALF_TOL
    out: dict[str, Any] = {
        "N": N, "L": BOX[N], "eta": eta, "seed": seed, "position_update": update,
        "burn_in_steps": burn, "measurement_steps_initial": meas, "sample_every": SAMPLE_EVERY,
        "half_means_initial": [h1, h2], "half_diff_initial": h1 - h2, "extended": extended,
    }
    if extended:
        measure(meas_ext - meas)
        n1 = len(phis)
        out["half_means_extended"] = [float(np.mean(phis[: n1 // 2])), float(np.mean(phis[n1 // 2:]))]
    ph = np.asarray(phis)
    out["measurement_steps_used"] = len(ph) * SAMPLE_EVERY
    out["n_samples"] = len(ph)
    out["phi_mean"] = float(ph.mean())
    out["phi_sd_time"] = float(ph.std())
    out["block_means_1e4_steps"] = [float(b.mean()) for b in np.array_split(ph, max(1, len(ph) // 1000))]
    out["wall_s"] = time.time() - t0
    out["cpu_s"] = time.process_time() - c0
    if save_series:
        RAWDIR.mkdir(parents=True, exist_ok=True)
        np.savez_compressed(RAWDIR / (task_id(N, eta, seed, update) + "_phi.npz"),
                            phi=ph.astype(np.float32))
    return out


def run_task(task: dict[str, Any]) -> dict[str, Any]:
    return run_cell(task["N"], task["eta"], task["seed"], task["update"])


# ---------------------------------------------------------------- aggregation
def check_targets_against_digitization() -> dict[str, Any]:
    """Re-read the frozen digitization JSON and confirm the targets (3 d.p.)."""
    d = read_json(SUPPORT / "p5_vicsek1995_fig2a_digitized.json")
    keys = {400: "N=400 (x)", 100: "N=100 (+)", 40: "N=40  (square)", 4000: "N=4000 (triangle)"}
    res: dict[str, Any] = {}
    ok_all = True
    for N, key in keys.items():
        pts = d[key]["pts"]
        etas = ETAS if N != 4000 else ETAS_4000
        tg = TARGETS[N] if N != 4000 else [TARGETS_4000[e] for e in ETAS_4000]
        rows = []
        for eta, t in zip(etas, tg):
            near = [p for p in pts if abs(p[0] - eta) < 0.03]
            v = round(near[0][1], 3) if near else None
            ok = v is not None and abs(v - t) < 5e-4 + 1e-12
            ok_all &= ok
            rows.append({"eta": eta, "digitized": near[0] if near else None, "frozen_target": t, "match": ok})
        res[f"N={N}"] = rows
    res["all_match"] = ok_all
    return res


def straightest_line_beta(etas: np.ndarray, v: np.ndarray) -> dict[str, Any]:
    """Vicsek's procedure (descriptive): choose eta_c so ln v vs ln((eta_c-eta)/eta_c) is straightest
    (max R^2 of an OLS line over points with eta < eta_c); beta = slope."""
    best = None
    for ec in np.arange(2.6, 4.6, 0.001):
        sel = etas < ec - 1e-9
        if sel.sum() < 3:
            continue
        xx = np.log((ec - etas[sel]) / ec)
        yy = np.log(v[sel])
        A = np.vstack([xx, np.ones_like(xx)]).T
        coef, *_ = np.linalg.lstsq(A, yy, rcond=None)
        resid = yy - A @ coef
        r2 = 1 - np.sum(resid ** 2) / np.sum((yy - yy.mean()) ** 2)
        if best is None or r2 > best["r2"]:
            best = {"eta_c": float(ec), "beta": float(coef[0]), "r2": float(r2), "n_points": int(sel.sum()),
                    "etas_used": etas[sel].tolist()}
    return best or {}


def aggregate(raw_dir: Path | None = None) -> dict[str, Any]:
    raw_dir = RAWDIR if raw_dir is None else raw_dir
    res: dict[str, Any] = {}
    cells: dict[tuple, list] = {}
    runs = []
    for t in all_tasks(include_optional=True):
        f = raw_dir / (t["id"] + ".json")
        if f.exists():
            r = read_json(f)
            runs.append(r)
            cells.setdefault((t["N"], t["update"], t["eta"]), []).append(r)
    res["n_runs_found"] = len(runs)
    res["runtime_cpu_h_total"] = sum(r["result"]["cpu_s"] for r in runs) / 3600
    res["runtime_wall_h_total_summed_over_runs"] = sum(r["result"]["wall_s"] for r in runs) / 3600
    res["commits"] = sorted({r["git"]["commit"] for r in runs if r.get("git")})
    res["any_code_dirty"] = any(r["git"].get("code_dirty") for r in runs if r.get("git"))

    def series(N: int, update: str, etas: list[float], targets: list[float] | None, tol: float,
               seeds_expected: int) -> dict[str, Any]:
        rows = []
        complete = True
        for i, eta in enumerate(etas):
            rr = sorted(cells.get((N, update, eta), []), key=lambda r: r["result"]["seed"])
            vals = np.array([r["result"]["phi_mean"] for r in rr])
            if len(vals) < seeds_expected:
                complete = False
            ti = t_interval(vals) if len(vals) else {}
            row = {"eta": eta, "seeds": [r["result"]["seed"] for r in rr],
                   "seed_means": vals.tolist(),
                   "extended_seeds": [r["result"]["seed"] for r in rr if r["result"]["extended"]],
                   "half_diffs_initial": [r["result"]["half_diff_initial"] for r in rr],
                   **{f"v_sim_{k}" if k == "mean" else k: v for k, v in ti.items()}}
            if targets is not None:
                row["v_pub"] = targets[i]
                if len(vals):
                    row["diff"] = row["v_sim_mean"] - targets[i]
                    row["within_tol"] = bool(abs(row["diff"]) <= tol)
            rows.append(row)
        out = {"N": N, "L": BOX[N], "position_update": update, "rows": rows, "complete": complete}
        if targets is not None and complete:
            out["max_abs_diff"] = max(abs(r["diff"]) for r in rows)
            out["all_within_tol"] = all(r["within_tol"] for r in rows)
        return out

    res["N400_old_velocity"] = series(400, "old_velocity", ETAS, TARGETS[400], TOL, 10)
    res["N100_old_velocity"] = series(100, "old_velocity", ETAS, TARGETS[100], TOL, 10)
    res["N40_old_velocity"] = series(40, "old_velocity", ETAS, TARGETS[40], TOL, 10)
    res["N400_new_velocity_sensitivity"] = series(400, "new_velocity", ETAS, TARGETS[400], TOL, 10)
    res["N4000_old_velocity_optional"] = series(4000, "old_velocity", ETAS_4000,
                                                [TARGETS_4000[e] for e in ETAS_4000], TOL_4000, 5)

    rules: dict[str, Any] = {}
    p = res["N400_old_velocity"]
    rules["primary_N400"] = {
        "rule": "|v_sim - v_pub| <= 0.04 at all 9 eta, N=400 L=10, old_velocity",
        "evaluated": p["complete"],
        "outcome": ("PASS" if p.get("all_within_tol") else "FAIL") if p["complete"] else "INCOMPLETE",
        "max_abs_diff": p.get("max_abs_diff"),
        "failing_etas": [r["eta"] for r in p["rows"] if r.get("within_tol") is False],
    }
    for N in (100, 40):
        s = res[f"N{N}_old_velocity"]
        rules[f"secondary_N{N}_reported"] = {
            "rule": f"|v_sim - v_pub| <= 0.04 at all 9 eta, N={N}",
            "outcome": ("PASS" if s.get("all_within_tol") else "FAIL") if s["complete"] else "INCOMPLETE",
            "max_abs_diff": s.get("max_abs_diff"),
            "failing_etas": [r["eta"] for r in s["rows"] if r.get("within_tol") is False],
        }
    order_rows = []
    if all(res[f"N{N}_old_velocity"]["complete"] for N in (40, 100, 400)):
        for i, eta in enumerate(ETAS):
            if eta < 2.5:
                continue
            v40 = res["N40_old_velocity"]["rows"][i]["v_sim_mean"]
            v100 = res["N100_old_velocity"]["rows"][i]["v_sim_mean"]
            v400 = res["N400_old_velocity"]["rows"][i]["v_sim_mean"]
            order_rows.append({"eta": eta, "v40": v40, "v100": v100, "v400": v400, "ordered": v40 > v100 > v400})
        rules["ordering_check_reported"] = {
            "rule": "v(40) > v(100) > v(400) at every eta >= 2.5 (old_velocity)",
            "outcome": "PASS" if all(r["ordered"] for r in order_rows) else "FAIL", "rows": order_rows}
    else:
        rules["ordering_check_reported"] = {"outcome": "INCOMPLETE"}
    s = res["N400_new_velocity_sensitivity"]
    rules["sensitivity_N400_new_velocity_reported"] = {
        "rule": "same 0.04 rule applied descriptively to the catalog's update order (not pass/fail)",
        "within_0.04_at_all_eta": s.get("all_within_tol"), "max_abs_diff": s.get("max_abs_diff"),
        "failing_etas": [r["eta"] for r in s["rows"] if r.get("within_tol") is False]}
    s = res["N4000_old_velocity_optional"]
    rules["optional_N4000_reported"] = {
        "rule": "|v_sim - v_pub| <= 0.06 at eta in {2.0, 2.5, 3.0}, 5 seeds (reported, not pass/fail)",
        "run": any(len(r["seeds"]) for r in s["rows"]),
        "within_0.06_at_all_eta": s.get("all_within_tol"), "max_abs_diff": s.get("max_abs_diff")}
    if p["complete"]:
        etas = np.array(ETAS)
        v = np.array([r["v_sim_mean"] for r in p["rows"]])
        rules["beta_descriptive"] = {
            "procedure": "eta_c scanned 2.600-4.599 step 0.001; OLS of ln v_sim on ln((eta_c-eta)/eta_c) over "
                         "grid points with eta < eta_c (>= 3 points); eta_c chosen at maximum R^2 (straightest "
                         "line); beta = slope. N=400 old_velocity means. Descriptive only; published "
                         "beta = 0.45 +/- 0.07.",
            **straightest_line_beta(etas, v)}
    res["rules"] = rules
    res["target_check_against_digitization"] = check_targets_against_digitization()
    return res
