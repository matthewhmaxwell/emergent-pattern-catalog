"""Runs the planted-pipeline verification jobs (census.planted_pipeline.verify's job list, its _vworker unchanged) under the
plan's compute cap, which the frozen driver does not enforce (PLANTED_PIPELINE_PLAN step 5: "A condition that cannot
finish within 24 CPU-hours per model is recorded as not run"; PLANTED_PIPELINE_DEVIATIONS notes: applied per condition).
Written 2026-11-04, after the gate runs started and before any verification output existed. Rule (PLANTED_PIPELINE_DEVIATIONS P6):

  projected CPU-h P = flag-stage seconds (rows/NXX.json; 3 seeds) x 10/3 seeds x m / 3600
      conditions: total_s (real + knock-out), m = 1 census, 4 area4, 3 len3; sub-programs: sim_s (real program), m = 1
  P > 36 (1.5 x cap)          -> condition recorded as not run, no attempt
  otherwise                   -> run in its own process with RLIMIT_CPU = 24 CPU-h; killed by the limit -> not run
  sub-programs of a model with area4 or len3 not run -> not run: verify.classify_one returns ARTIFACT (or NOT RUNNABLE)
      before it reads sub-programs, so they cannot change the class

A job that is not run writes no part file, so census.planted_pipeline.assemble records it as "not run".

  python -m census.planted_verify_capped plan          -> verify/capped_plan.json (no simulation)
  python -m census.planted_verify_capped run [--workers 6]
"""
import argparse, glob, itertools, json, os, resource, signal, subprocess, sys, time
from concurrent.futures import ThreadPoolExecutor, as_completed

from census import planted_pipeline as P

CAP_H = 24.0
SKIP_OVER_H = 36.0
MULT = {"census": 1, "area4": 4, "len3": 3}


def jobs():
    """census.planted_pipeline.verify's job list, built the same way."""
    out = []
    for n in P.IN_SCOPE:
        if not glob.glob(f"{P.OUT}/cards/pipeline__{n}__*.json"): continue
        if not any(P.final_verdicts(P.J(c)["class_id"].replace(":", "__")) for c in glob.glob(f"{P.OUT}/cards/pipeline__{n}__*.json")): continue
        inter = P.interaction_rules(P.module(n))
        out += [(n, "cond", c) for c in ("census", "area4", "len3")]
        if len(inter) >= 2:
            out += [(n, "sub", sub) for k in range(1, len(inter)) for sub in itertools.combinations(inter, k)]
    return out


def tag(job):
    n, kind, arg = job
    return arg if kind == "cond" else "sub_" + "-".join(map(str, arg))


def plan():
    rows, skipped_cond = [], set()
    js = jobs()
    for n, kind, arg in js:
        r = P.J(f"{P.OUT}/rows/{n}.json")
        secs = r["total_s"] if kind == "cond" else r["sim_s"]
        p = secs * (P.VSEEDS / P.N_SEEDS) * (MULT[arg] if kind == "cond" else 1) / 3600
        d = "not run: projected over 1.5 x cap" if kind == "cond" and p > SKIP_OVER_H else "run (24 CPU-h limit)"
        if d.startswith("not run") and arg in ("area4", "len3"): skipped_cond.add(n)
        rows.append({"card": n, "kind": kind, "arg": list(arg) if kind == "sub" else arg, "tag": tag((n, kind, arg)),
                     "flag_stage_s": secs, "projected_cpu_h": round(p, 2), "decision": d})
    for x in rows:
        if x["kind"] == "sub" and x["card"] in skipped_cond:
            x["decision"] = "not run: area4 or len3 not run, so verify.classify_one cannot read sub-programs"
    res = {"written": time.strftime("%Y-%m-%d %H:%M:%S"), "cap_cpu_h": CAP_H, "skip_over_cpu_h": SKIP_OVER_H,
           "n_jobs_driver": len(js), "jobs": rows,
           "run_cpu_h_projected": round(sum(x["projected_cpu_h"] for x in rows if x["decision"].startswith("run")), 1)}
    P.dump(res, f"{P.OUT}/verify/capped_plan.json")
    for x in rows: print(f"{x['card']} {x['tag']:>10} {x['projected_cpu_h']:7.2f} h  {x['decision']}")
    print(len(js), "driver jobs;", sum(x["decision"].startswith("run") for x in rows), "to run; projected",
          res["run_cpu_h_projected"], "CPU-h")
    return res


def _limit():
    s = int(CAP_H * 3600)
    resource.setrlimit(resource.RLIMIT_CPU, (s, s + 60))
    os.nice(10)


CHILD = ("import json, sys; from census import planted_pipeline as P; j = json.loads(sys.argv[1]); "
         "j = (j[0], j[1], tuple(j[2]) if j[1] == 'sub' else j[2]); print(P._vworker(j)[1], flush=True)")


def one(x):
    job = [x["card"], x["kind"], x["arg"]]
    fn = f"{P.OUT}/verify/parts/{x['card']}__{x['tag']}.json"
    if os.path.exists(fn): return x, "cached", None
    t = time.time()
    pr = subprocess.run([sys.executable, "-c", CHILD, json.dumps(job)], preexec_fn=_limit, capture_output=True, text=True,
                        env={**os.environ, "PYTHONPATH": os.getcwd()})
    if pr.returncode != 0:
        killed = pr.returncode in (-signal.SIGXCPU, -signal.SIGKILL)
        return x, ("KILLED by the 24 CPU-h limit" if killed else f"FAILED rc={pr.returncode}: {pr.stderr[-800:]}"), round(time.time() - t)
    return x, pr.stdout.strip(), round(time.time() - t)


def run(workers):
    pl = plan(); todo = [x for x in pl["jobs"] if x["decision"].startswith("run")]
    todo.sort(key=lambda x: -x["projected_cpu_h"])                         # longest first
    log = P.J(f"{P.OUT}/verify/capped_log.json", [])
    with ThreadPoolExecutor(workers) as ex:
        for fut in as_completed([ex.submit(one, x) for x in todo]):
            x, msg, wall = fut.result()
            print("verify", x["card"], x["tag"], msg, f"{wall}s" if wall else "", flush=True)
            log.append({"card": x["card"], "tag": x["tag"], "result": msg, "wall_s": wall, "at": time.strftime("%Y-%m-%d %H:%M:%S")})
            P.dump(log, f"{P.OUT}/verify/capped_log.json")
    print("VERIFY DONE", time.strftime("%Y-%m-%d %H:%M:%S"), flush=True)


if __name__ == "__main__":
    ap = argparse.ArgumentParser(); ap.add_argument("cmd", choices=["plan", "run"]); ap.add_argument("--workers", type=int, default=6)
    a = ap.parse_args()
    plan() if a.cmd == "plan" else run(a.workers)
