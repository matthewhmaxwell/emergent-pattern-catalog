"""Work queue for the repro_v2 group B runs (P5, P6, P10).

Usage (from the worktree root, as the non-root user):
  python analysis/reproductions/repro_v2/group_b_queue.py worker      # run tasks until none are left
  python analysis/reproductions/repro_v2/group_b_queue.py status
  python analysis/reproductions/repro_v2/group_b_queue.py smoke       # tiny sizes; informs no pass rule
  python analysis/reproductions/repro_v2/group_b_queue.py aggregate   # write analysis/outputs/repro_v2/P*.json

At most two workers are started (shared machine). Tasks are claimed with an O_EXCL lock file; each finished task
writes raw/<pattern>/<task>.json (ignored by git). The task list is rebuilt after every task, so tasks whose
inputs depend on earlier results (P6 runs after the dt check, P10 Stage B after Stage A) appear when ready.
"""
from __future__ import annotations

import json
import os
import socket
import sys
import time
from pathlib import Path
from typing import Any

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import common  # noqa: E402
from common import RAW, OUT, git_info, read_json, write_json  # noqa: E402
import p5_vicsek1995 as p5  # noqa: E402
import p6_dorsogna2006 as p6  # noqa: E402
import p10_abrams_strogatz2004 as p10  # noqa: E402

MODULES = {"P5": p5, "P6": p6, "P10": p10}
SKIP_OPTIONAL_FILE = RAW / "SKIP_OPTIONAL"   # create to drop optional tasks (P5 N=4000) from the queue


def build_tasks() -> tuple[list[tuple[str, dict[str, Any]]], int]:
    """Ordered task list (priority order) and the number of tasks still blocked on earlier results."""
    T: list[tuple[str, dict[str, Any]]] = []
    blocked = 0
    T += [("P6", t) for t in p6.ring_tasks_dtcheck()]
    T += [("P10", t) for t in p10.tasks_stage_a(256, 0.18)]
    tb, b = p10.tasks_stage_b(256); T += [("P10", t) for t in tb]; blocked += b
    dtc = p6.dt_check_outcome()
    if dtc is None:
        blocked += 1
    else:
        T += [("P6", t) for t in p6.ring_tasks_main(dtc["dt_for_T61"])]
    p5t = p5.all_tasks(include_optional=not SKIP_OPTIONAL_FILE.exists())
    T += [("P5", t) for t in p5t if t["N"] == 400 and t["update"] == "old_velocity"]
    tf, b = p10.tasks_fallback(); T += [("P10", t) for t in tf]; blocked += b
    td, b = p10.tasks_dtcheck(); T += [("P10", t) for t in td]; blocked += b
    T += [("P10", t) for t in p10.tasks_stage_a(256, 0.05)]
    T += [("P5", t) for t in p5t if t["N"] in (100, 40)]
    T += [("P5", t) for t in p5t if t["N"] == 400 and t["update"] == "new_velocity"]
    T += [("P6", t) for t in p6.energy_tasks()]
    T += [("P10", t) for t in p10.tasks_stage_a(128, 0.18)]
    tb, b = p10.tasks_stage_b(128); T += [("P10", t) for t in tb]; blocked += b
    T += [("P5", t) for t in p5t if t["N"] == 4000]
    return T, blocked


def paths(pat: str, tid: str) -> tuple[Path, Path]:
    d = RAW / pat
    d.mkdir(parents=True, exist_ok=True)
    return d / f"{tid}.json", d / f"{tid}.lock"


def claim(pat: str, tid: str) -> bool:
    res, lock = paths(pat, tid)
    if res.exists():
        return False
    try:
        fd = os.open(lock, os.O_CREAT | os.O_EXCL | os.O_WRONLY)
    except FileExistsError:
        return False
    with os.fdopen(fd, "w") as fh:
        fh.write(json.dumps({"pid": os.getpid(), "host": socket.gethostname(), "t": time.time()}))
    if res.exists():  # finished between the check and the lock
        os.remove(lock)
        return False
    return True


def worker() -> None:
    name = f"pid{os.getpid()}"
    while True:
        tasks, blocked = build_tasks()
        pending = [(p, t) for p, t in tasks if not paths(p, t["id"])[0].exists()]
        ran = False
        for pat, task in pending:
            if not claim(pat, task["id"]):
                continue
            res_path, lock = paths(pat, task["id"])
            started = time.time()
            print(f"[{time.strftime('%H:%M:%S')}] {name} start {task['id']}", flush=True)
            try:
                result = MODULES[pat].run_task(task)
            except Exception as exc:
                import traceback
                write_json(res_path.with_suffix(".error.json"),
                           {"task": task, "error": repr(exc), "traceback": traceback.format_exc()})
                os.remove(lock)
                print(f"[{time.strftime('%H:%M:%S')}] {name} ERROR {task['id']}: {exc!r}", flush=True)
                raise
            write_json(res_path, {"task": task, "git": git_info(), "worker": name,
                                  "started": started, "finished": time.time(), "result": result})
            os.remove(lock)
            print(f"[{time.strftime('%H:%M:%S')}] {name} done  {task['id']} "
                  f"({(time.time() - started) / 60:.1f} min)", flush=True)
            ran = True
            break
        if ran:
            continue
        in_progress = [t for p, t in pending if paths(p, t["id"])[1].exists()]
        if not pending and blocked == 0:
            print(f"[{time.strftime('%H:%M:%S')}] {name}: queue empty, exiting", flush=True)
            return
        time.sleep(60)


def status() -> None:
    tasks, blocked = build_tasks()
    done = sum(paths(p, t["id"])[0].exists() for p, t in tasks)
    running = [t["id"] for p, t in tasks if paths(p, t["id"])[1].exists()]
    by = {}
    for p, t in tasks:
        d = by.setdefault(p, [0, 0])
        d[0] += 1
        d[1] += paths(p, t["id"])[0].exists()
    print(f"tasks listed {len(tasks)}, done {done}, blocked groups {blocked}")
    for p, (n, dn) in by.items():
        print(f"  {p}: {dn}/{n}")
    print("running:", running)
    errs = list(RAW.glob("*/*.error.json"))
    if errs:
        print("errors:", [e.name for e in errs])


DECL = OUT / "declarations_group_B.json"
DEVS = OUT / "deviations_group_B.json"


def aggregate() -> None:
    decl = read_json(DECL) if DECL.exists() else {}
    devs = read_json(DEVS) if DEVS.exists() else {}
    for pat, mod in MODULES.items():
        d = decl.get(pat, {})
        out = {"pattern": pat, "spec": "docs/validation_rebuild/repro_v2/group_B.md",
               "plan": "docs/validation_rebuild/repro_v2/REPRO_V2_PLAN.md (frozen at 39688e3)",
               "implementation_readings_declared_before_full_runs": d.get("readings", []),
               "deviations": d.get("deviations", []) + devs.get(pat, []),
               "not_run": d.get("not_run", []),
               "aggregated_at": time.strftime("%Y-%m-%d %H:%M:%S %Z"),
               "aggregate_code": git_info(),
               **mod.aggregate()}
        write_json(OUT / f"{pat}.json", out)
        print("wrote", OUT / f"{pat}.json")


def smoke() -> None:
    """Tiny sizes and lengths only: exercises every code path, informs no pass rule."""
    import numpy as np
    t = time.time()
    r = p5.run_cell(40, 2.0, 1, "old_velocity", burn=50, meas=200, meas_ext=400, save_series=False)
    print("P5 smoke", {k: r[k] for k in ("phi_mean", "extended", "n_samples")})
    m = p5.VicsekModel(n_particles=5, box_size=3.1, seed=0, wrap_guard=True)
    m.setup(); m.positions[0, 0] = -1e-17; m.positions %= 3.1
    m.positions[m.positions >= 3.1] = 0.0; m.step()
    print("P5 guard ok")
    old = (p6.T_END, p6.T_WIN0, p6.T_EXT_END, p6.T_EXT_WIN0)
    p6.T_END, p6.T_WIN0, p6.T_EXT_END, p6.T_EXT_WIN0 = 2.0, 1.0, 3.0, 2.5
    r = p6.run_ring(20, 1, 0.01)
    print("P6 smoke ring", r["valid"], r["state_label"], r["extended"], round(r["used"]["R"], 4))
    p6.T_END, p6.T_WIN0, p6.T_EXT_END, p6.T_EXT_WIN0 = old
    oldT = p6.T_U; p6.T_U = 0.5
    r = p6.run_energy("catastrophic", 20, 1); print("P6 smoke energy", r["U_total"])
    p6.T_U = oldT
    r = p10.stage_a(32, 0.18, 0, T=40, win0=20)
    print("P10 smoke A", r["estimators"], r["chimera"])
    oldramp = p10.RAMP
    p10.RAMP = [(0.20, 220), (0.27, 220), (0.28, 220)]
    th = np.zeros(32) + 0.01 * np.arange(32)
    r = p10.stage_b(32, 0, th)
    print("P10 smoke B", r["beta_collapse"], r["confirmation_hold_still_collapsed"], len(r["holds"]))
    p10.RAMP = oldramp
    print(f"smoke total {time.time() - t:.1f} s")


if __name__ == "__main__":
    cmd = sys.argv[1] if len(sys.argv) > 1 else "status"
    {"worker": worker, "status": status, "aggregate": aggregate, "smoke": smoke}[cmd]()
