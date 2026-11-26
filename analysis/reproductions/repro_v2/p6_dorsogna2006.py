"""repro_v2 P6: D'Orsogna et al. 2006 PRL, ring radius Eq. (4) and energy exponent (group_B.md, section P6).

T6.1 (primary): C_a=1, l_a=1, C_r=0.5, l_r=0.5, alpha=1, beta=0.5; N in {100,200,300,400}; seeds 2000-2009;
  RK4 dt=0.01, T=200, window t in [150, 200] sampled every 0.1; init_radius=1.0 (uniform in [-1,1]^2).
  dt check: seed 2000 at N=100 and N=400 with dt=0.005; R must change by < 0.1%, else all N rerun at dt=0.005.
  Gates (time means over window): radial CV < 0.05; |L| > 0.95 (open-space AngularMomentumMetric);
  mean speed within 2% of sqrt(2); max |x_i - x_cm| < 2R; first/second-half R differ < 0.5% (else extend once
  to T=400).
  PASS: >= 8/10 seeds pass all gates at every N AND |R_sim/R_Eq4 - 1| <= 0.03 at every N (R_sim = median over
  valid seeds).
T6.2 (reported): |R_sim/R_exact - 1| <= 0.01 at every N.
T6.3 (reported): |U_total| exponent; catastrophic (C_a=0.5, C_r=1, l_a=2, l_r=0.5, beta=0.5, alpha=1) in [1.9, 2.1];
  H-stable (l_r=1.5, alpha=0.5) in [0.9, 1.1]; N in {50,100,200,400}; 5 seeds; dt=0.01, T=100.
"""
from __future__ import annotations

import re
import subprocess
import sys
import time
from typing import Any

import numpy as np

from common import RAW, SUPPORT, git_info, read_json, t_interval  # noqa: E402

from epc.models.dorsogna_spp import DOrsognaSPPModel
from epc.metrics.collective_motion import AngularMomentumMetric

RING = dict(C_a=1.0, l_a=1.0, C_r=0.5, l_r=0.5, alpha=1.0, beta=0.5)
NS = [100, 200, 300, 400]
SEEDS = list(range(2000, 2010))
DT = 0.01
DT_CHECK = 0.005
T_END, T_WIN0 = 200.0, 150.0
T_EXT_END, T_EXT_WIN0 = 400.0, 350.0
SAMPLE_DT = 0.1
INIT_RADIUS = 1.0
V0 = np.sqrt(2.0)
R_EQ4 = {100: 0.17426, 200: 0.11503, 300: 0.09133, 400: 0.07783}     # frozen (group_B.md P6 table)
R_EXACT = {100: 0.17630, 200: 0.11567, 300: 0.09166, 400: 0.07804}   # frozen
TOL_EQ4, TOL_EXACT = 0.03, 0.01
MIN_VALID = 8

CAT = dict(C_a=0.5, C_r=1.0, l_a=2.0, l_r=0.5, beta=0.5, alpha=1.0)
HST = dict(C_a=0.5, C_r=1.0, l_a=2.0, l_r=1.5, beta=0.5, alpha=0.5)
NS_U = [50, 100, 200, 400]
SEEDS_U = list(range(2000, 2005))
T_U = 100.0

RAWDIR = RAW / "P6"


def ring_id(N: int, seed: int, dt: float) -> str:
    return f"P6_T61_N{N}_s{seed}_dt{dt:g}"


def energy_id(setname: str, N: int, seed: int) -> str:
    return f"P6_T63_{setname}_N{N}_s{seed}"


def snapshot(m: DOrsognaSPPModel) -> dict[str, float]:
    x, v = m.positions, m.velocities
    cm = x.mean(axis=0)
    d = np.linalg.norm(x - cm, axis=1)
    R = float(d.mean())
    L = float(AngularMomentumMetric.compute_instant({"positions": x, "velocities": v}, box_size=None))
    sp = np.linalg.norm(v, axis=1)
    pol = float(np.linalg.norm((v / sp[:, None]).mean(axis=0)))
    return {"R": R, "cv": float(d.std() / R), "L": L, "speed": float(sp.mean()), "maxd": float(d.max()),
            "polarization": pol}


def run_ring(N: int, seed: int, dt: float = DT) -> dict[str, Any]:
    t0, c0 = time.time(), time.process_time()
    m = DOrsognaSPPModel(n_particles=N, dt=dt, init_mode="random", init_radius=INIT_RADIUS, seed=seed, **RING)
    m.setup()
    per = int(round(SAMPLE_DT / dt))
    steps_done = 0

    def advance_to(t: float) -> None:
        nonlocal steps_done
        target = int(round(t / dt))
        while steps_done < target:
            m.step()
            steps_done += 1

    def window(t0w: float, t1w: float) -> list[dict[str, float]]:
        nonlocal steps_done
        advance_to(t0w)
        snaps = []
        n = int(round((t1w - t0w) / SAMPLE_DT))
        for _ in range(n):
            for _ in range(per):
                m.step()
                steps_done += 1
            snaps.append(snapshot(m))
        return snaps

    def summarize(snaps: list[dict[str, float]], win: list[float]) -> dict[str, Any]:
        arr = {k: np.array([s[k] for s in snaps]) for k in snaps[0]}
        n = len(snaps)
        R = float(arr["R"].mean())
        h1, h2 = float(arr["R"][: n // 2].mean()), float(arr["R"][n // 2:].mean())
        s = {"window": win, "n_samples": n, "R": R, "R_half_means": [h1, h2],
             "R_half_reldiff": abs(h2 / h1 - 1.0),
             "cv": float(arr["cv"].mean()), "L": float(arr["L"].mean()), "speed": float(arr["speed"].mean()),
             "maxd": float(arr["maxd"].mean()), "polarization": float(arr["polarization"].mean()),
             "R_series": arr["R"][::10].tolist()}
        s["gates"] = {
            "radial_cv_lt_0.05": s["cv"] < 0.05,
            "absL_gt_0.95": abs(s["L"]) > 0.95,
            "speed_within_2pct_sqrt2": abs(s["speed"] / V0 - 1.0) < 0.02,
            "maxd_lt_2R": s["maxd"] < 2 * R,
            "converged_halves_lt_0.5pct": s["R_half_reldiff"] < 0.005,
        }
        return s

    first = summarize(window(T_WIN0, T_END), [T_WIN0, T_END])
    out: dict[str, Any] = {"N": N, "seed": seed, "dt": dt, "init_radius": INIT_RADIUS, "params": RING,
                           "window_T200": first}
    used = first
    if not first["gates"]["converged_halves_lt_0.5pct"]:
        used = summarize(window(T_EXT_WIN0, T_EXT_END), [T_EXT_WIN0, T_EXT_END])
        out["window_T400_extension"] = used
    out["extended"] = used is not first
    out["used"] = used
    out["valid"] = all(used["gates"].values())
    out["state_label"] = label_state(used)
    RAWDIR.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(RAWDIR / (ring_id(N, seed, dt) + "_final.npz"), x=m.positions, v=m.velocities)
    out["final_time"] = steps_done * dt
    out["wall_s"] = time.time() - t0
    out["cpu_s"] = time.process_time() - c0
    return out


def label_state(s: dict[str, Any]) -> str:
    if all(s["gates"].values()):
        return "co-rotating ring"
    if s["polarization"] > 0.9:
        return "translating flock/clump"
    if s["cv"] < 0.05 and abs(s["L"]) <= 0.95:
        return "ring, not co-rotating (mixed rotation)"
    if abs(s["L"]) > 0.95 and s["cv"] >= 0.05:
        return "co-rotating but not a thin ring (mill/annulus)"
    if not s["gates"]["maxd_lt_2R"]:
        return "escaped particles"
    return "other"


def run_energy(setname: str, N: int, seed: int) -> dict[str, Any]:
    t0, c0 = time.time(), time.process_time()
    P = CAT if setname == "catastrophic" else HST
    m = DOrsognaSPPModel(n_particles=N, dt=DT, init_mode="random", init_radius=INIT_RADIUS, seed=seed, **P)
    m.setup()
    for _ in range(int(round(T_U / DT))):
        m.step()
    x = m.positions
    d = np.linalg.norm(x[:, None, :] - x[None, :, :], axis=2)
    U = P["C_r"] * np.exp(-d / P["l_r"]) - P["C_a"] * np.exp(-d / P["l_a"])
    np.fill_diagonal(U, 0.0)
    Utot = float(U.sum())  # sum_i U(x_i), U(x_i) = sum_{j != i}: pair sum, double counted
    cm = x.mean(axis=0)
    rr = np.linalg.norm(x - cm, axis=1)
    sp = np.linalg.norm(m.velocities, axis=1)
    L = float(AngularMomentumMetric.compute_instant({"positions": x, "velocities": m.velocities}, box_size=None))
    return {"set": setname, "params": P, "N": N, "seed": seed, "dt": DT, "T": T_U, "init_radius": INIT_RADIUS,
            "U_total": Utot, "abs_U_total": abs(Utot), "mean_radius": float(rr.mean()), "max_radius": float(rr.max()),
            "mean_speed": float(sp.mean()), "L": L, "finite": bool(np.isfinite(x).all()),
            "wall_s": time.time() - t0, "cpu_s": time.process_time() - c0}


def ring_tasks_dtcheck() -> list[dict[str, Any]]:
    t = []
    for N in (400, 100):
        for dt in (DT, DT_CHECK):
            t.append(dict(id=ring_id(N, 2000, dt), kind="ring", N=N, seed=2000, dt=dt))
    return t


def dt_check_outcome() -> dict[str, Any] | None:
    """None until all four dt-check runs exist."""
    rows = []
    for N in (100, 400):
        f1, f2 = RAWDIR / (ring_id(N, 2000, DT) + ".json"), RAWDIR / (ring_id(N, 2000, DT_CHECK) + ".json")
        if not (f1.exists() and f2.exists()):
            return None
        r1, r2 = read_json(f1)["result"], read_json(f2)["result"]
        rel = abs(r2["used"]["R"] / r1["used"]["R"] - 1.0)
        rows.append({"N": N, "seed": 2000, "R_dt0.01": r1["used"]["R"], "R_dt0.005": r2["used"]["R"],
                     "rel_change": rel, "lt_0.1pct": rel < 0.001})
    ok = all(r["lt_0.1pct"] for r in rows)
    return {"rows": rows, "passed": ok, "dt_for_T61": DT if ok else DT_CHECK}


def ring_tasks_main(dt_main: float) -> list[dict[str, Any]]:
    t = []
    for N in (400, 300, 200, 100):
        for s in SEEDS:
            t.append(dict(id=ring_id(N, s, dt_main), kind="ring", N=N, seed=s, dt=dt_main))
    return t


def energy_tasks() -> list[dict[str, Any]]:
    t = []
    for setname in ("catastrophic", "hstable"):
        for N in sorted(NS_U, reverse=True):
            for s in SEEDS_U:
                t.append(dict(id=energy_id(setname, N, s), kind="energy", set=setname, N=N, seed=s))
    return t


def run_task(task: dict[str, Any]) -> dict[str, Any]:
    if task["kind"] == "ring":
        return run_ring(task["N"], task["seed"], task["dt"])
    return run_energy(task["set"], task["N"], task["seed"])


# ---------------------------------------------------------------- aggregation
def regenerate_theory() -> dict[str, Any]:
    out = subprocess.run([sys.executable, str(SUPPORT / "p6_solve_dorsogna2006_eq4.py")],
                         capture_output=True, text=True, timeout=3600).stdout
    regen = {}
    for line in out.splitlines():
        mm = re.match(r"\s*(\d+)\s+\[([\d.]+)\]\s+\[([\d.]+)\]", line)
        if mm and int(mm.group(1)) in R_EQ4:
            regen[int(mm.group(1))] = (float(mm.group(2)), float(mm.group(3)))
    match = all(abs(regen[N][0] - R_EQ4[N]) < 6e-6 and abs(regen[N][1] - R_EXACT[N]) < 6e-6 for N in R_EQ4)
    return {"script": "group_B_support/p6_solve_dorsogna2006_eq4.py",
            "regenerated": {N: {"R_eq4": a, "R_exact": b} for N, (a, b) in regen.items()},
            "matches_frozen": match}


def aggregate() -> dict[str, Any]:
    res: dict[str, Any] = {}
    dtc = dt_check_outcome()
    res["dt_check"] = dtc
    dt_main = dtc["dt_for_T61"] if dtc else None
    t61 = {}
    if dt_main is not None:
        rows = []
        all_n_ok = True
        complete = True
        for N in NS:
            rr = []
            for s in SEEDS:
                f = RAWDIR / (ring_id(N, s, dt_main) + ".json")
                if f.exists():
                    rr.append(read_json(f)["result"])
            if len(rr) < len(SEEDS):
                complete = False
            valid = [r for r in rr if r["valid"]]
            Rv = np.array([r["used"]["R"] for r in valid])
            row = {"N": N, "n_seeds_run": len(rr), "n_valid": len(valid),
                   "per_seed": [{"seed": r["seed"], "valid": r["valid"], "state": r["state_label"],
                                 "extended_to_T400": r["extended"], "R": r["used"]["R"], "cv": r["used"]["cv"],
                                 "L": r["used"]["L"], "speed": r["used"]["speed"], "maxd": r["used"]["maxd"],
                                 "polarization": r["used"]["polarization"],
                                 "R_half_reldiff": r["used"]["R_half_reldiff"], "gates": r["used"]["gates"]}
                                for r in rr],
                   "R_Eq4": R_EQ4[N], "R_exact": R_EXACT[N]}
            if len(valid):
                row["R_sim_median"] = float(np.median(Rv))
                row["R_valid_range"] = [float(Rv.min()), float(Rv.max())]
                row["ratio_to_Eq4_minus_1"] = row["R_sim_median"] / R_EQ4[N] - 1
                row["ratio_to_exact_minus_1"] = row["R_sim_median"] / R_EXACT[N] - 1
                row["within_3pct_Eq4"] = abs(row["ratio_to_Eq4_minus_1"]) <= TOL_EQ4
                row["within_1pct_exact"] = abs(row["ratio_to_exact_minus_1"]) <= TOL_EXACT
            row["ge_8_valid"] = len(valid) >= MIN_VALID
            all_n_ok &= row["ge_8_valid"] and bool(row.get("within_3pct_Eq4"))
            rows.append(row)
        t61 = {"dt": dt_main, "rows": rows, "complete": complete,
               "outcome": ("PASS" if all_n_ok else "FAIL") if complete else "INCOMPLETE"}
        t62_ok = all(r.get("within_1pct_exact") for r in rows)
        res["T6_2"] = {"rule": "|R_sim/R_exact - 1| <= 0.01 at every N (reported)",
                       "outcome": ("PASS" if t62_ok else "FAIL") if complete else "INCOMPLETE",
                       "values": {r["N"]: r.get("ratio_to_exact_minus_1") for r in rows}}
    res["T6_1"] = {"rule": ">= 8 of 10 seeds pass all gates at every N and |R_sim/R_Eq4 - 1| <= 0.03 at every N",
                   **t61}
    # T6.3
    t63: dict[str, Any] = {}
    complete3 = True
    for setname, lo, hi in (("catastrophic", 1.9, 2.1), ("hstable", 0.9, 1.1)):
        pts = []
        for N in NS_U:
            us = []
            for s in SEEDS_U:
                f = RAWDIR / (energy_id(setname, N, s) + ".json")
                if f.exists():
                    us.append(read_json(f)["result"])
            if len(us) < len(SEEDS_U):
                complete3 = False
            pts.append({"N": N, "abs_U_per_seed": [u["abs_U_total"] for u in us],
                        "median_abs_U": float(np.median([u["abs_U_total"] for u in us])) if us else None,
                        "mean_radius_per_seed": [u["mean_radius"] for u in us],
                        "L_per_seed": [u["L"] for u in us], "finite": all(u["finite"] for u in us)})
        entry = {"params": CAT if setname == "catastrophic" else HST, "points": pts, "band": [lo, hi]}
        if all(p["median_abs_U"] for p in pts):
            xs = np.log([p["N"] for p in pts]); ys = np.log([p["median_abs_U"] for p in pts])
            slope = float(np.polyfit(xs, ys, 1)[0])
            entry["exponent"] = slope
            entry["in_band"] = lo <= slope <= hi
        t63[setname] = entry
    ok3 = all(t63[k].get("in_band") for k in t63)
    res["T6_3"] = {"rule": "exponent in [1.9, 2.1] (catastrophic, alpha=1) and [0.9, 1.1] (H-stable, alpha=0.5); "
                           "reported", "outcome": ("PASS" if ok3 else "FAIL") if complete3 else "INCOMPLETE", **t63}
    allraw = [read_json(f) for f in sorted(RAWDIR.glob("P6_*.json"))]
    res["runtime_cpu_h_total"] = sum(r["result"]["cpu_s"] for r in allraw) / 3600
    res["runtime_wall_h_total_summed_over_runs"] = sum(r["result"]["wall_s"] for r in allraw) / 3600
    res["n_runs"] = len(allraw)
    res["commits"] = sorted({r["git"]["commit"] for r in allraw if r.get("git")})
    res["any_code_dirty"] = any(r["git"].get("code_dirty") for r in allraw if r.get("git"))
    res["theory_regeneration"] = regenerate_theory()
    return res
