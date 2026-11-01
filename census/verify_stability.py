# post-hoc diagnostic: does each representative (and its knock-out) diverge numerically at census settings? classes unchanged
import json, multiprocessing as mp, os, sys
os.environ.setdefault("OMP_NUM_THREADS", "1")
def one(it):
    from census import verify as V, sim, knockout as KO
    p = V.grammar(it["grammar"]).decode(it["bits"]); kp = KO.knockout_program(p); s = V.fresh_seed(it["bits"])
    o = sim.run(p, seed0=s, S=2)["meta"]; k = sim.run(kp, seed0=s, S=2)["meta"] if kp is not None else {}
    return {"pid": it["pid"], "unstable": bool(o.get("unstable")), "steps_run": o.get("steps_run"),
            "ko_unstable": bool(k.get("unstable")), "ko_steps_run": k.get("steps_run")}
if __name__ == "__main__":
    inv = json.load(open("census/runs/verify/inventory.json"))
    with mp.get_context("fork").Pool(6, maxtasksperchild=1) as pool:
        res = list(pool.imap_unordered(one, inv))
    json.dump(sorted(res, key=lambda r: r["pid"]), open("census/runs/verify/stability.json", "w"), indent=1)
    bad = [r for r in res if r["unstable"] or r["ko_unstable"]]
    print("unstable real", sum(r["unstable"] for r in res), "unstable knock-out", sum(r["ko_unstable"] for r in res))
    for r in sorted(bad, key=lambda r: r["pid"]): print(r)
    print("STABILITY DONE")
