"""Census runner (pilot draft v0): enumerate -> simulate -> filter -> one JSON line per program. Resumable.

  python -m census.runner --out census/runs/<tag> --max-bits L [--min-bits l] [--workers 4] [--limit N] [--null]

Every canonical G1 program of length in [min_bits, max_bits] gets a stable index (enumeration order) and a stable
seed (hash of its bit string), so any program can be re-simulated exactly later — raw trajectories are not stored.
Worker w processes indices with idx % workers == w and appends to results_w.jsonl; on restart, finished indices are
skipped. Run inside tmux as matthewhmaxwell, niced.
"""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")                            # one thread per worker: workers are the parallelism
import argparse, hashlib, json, sys, time, traceback
from multiprocessing import Process

from census import grammar_g1 as G


def seed_of(bits):
    return int(hashlib.sha1(bits.encode()).hexdigest()[:8], 16)


def compact(ev):
    """keep what the map, triage and digest need; drop bulky per-seed detail."""
    out = {}
    for v, e in ev.items():
        out[v] = {"verdict": e["verdict"], "pattern": e["pattern"], "em_seeds": e["emergent_seeds"],
                  "seeds": [{k: s.get(k) for k in ("emergent", "em_score", "em_kind", "is_complex", "C", "psi",
                                                   "collapsed", "verdict", "pattern", "tier", "fired", "confirmed", "fp")}
                            for s in e["seeds"]]}
    return out


def worker(w, nw, progs, outdir, null, battery):
    from census import sim, filter as FL
    import glob
    path = os.path.join(outdir, f"results_{w}.jsonl"); done = set()
    for f in glob.glob(os.path.join(outdir, "results_*.jsonl")):          # all workers' files: resume is safe
        for line in open(f):                                              # even if the worker count changes
            try: done.add(json.loads(line)["idx"])
            except Exception: pass
    with open(path, "a") as fh:
        for idx, item in enumerate(progs):
            if item is None or idx % nw != w or idx in done: continue
            bits, p = item
            row = {"idx": idx, "bits": bits, "len": len(bits), "prog": G.describe(p), "layers": p.layers,
                   "n_rules": len(p.rules), "null": null}
            try:
                t0 = time.time(); out = sim.run(p, seed0=seed_of(bits), null_shuffle=null); t1 = time.time()
                ev = FL.evaluate(out, p, battery=battery); t2 = time.time()
                row.update({"sim_s": round(t1 - t0, 2), "filter_s": round(t2 - t1, 2), "meta": out["meta"],
                            "views": compact(ev)})
                flags = [e["verdict"] for e in ev.values()]
                row["status"] = ("UNCLASSIFIED" if "UNCLASSIFIED" in flags else "EMERGENT" if "EMERGENT" in flags
                                 else "MATCH" if "MATCH" in flags else "NONE")
            except Exception as e:
                row.update({"status": "ERROR", "error": f"{type(e).__name__}: {e}", "trace": traceback.format_exc()[-800:]})
            fh.write(json.dumps(row, default=str) + "\n"); fh.flush()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True); ap.add_argument("--max-bits", type=int, required=True)
    ap.add_argument("--min-bits", type=int, default=0); ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--limit", type=int, default=0); ap.add_argument("--null", action="store_true")
    ap.add_argument("--battery", action="store_true", help="per-program battery (reference mode; default = cluster-first)")
    ap.add_argument("--sample", type=int, default=0, help="fixed-seed random subsample of N programs (indices kept)")
    a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
    progs = [(b, p) for b, p in G.enumerate_programs(a.max_bits) if len(b) >= a.min_bits]
    progs.sort(key=lambda x: (len(x[0]), x[0]))                      # shortest first, deterministic
    if a.limit: progs = progs[:a.limit]
    if a.sample and a.sample < len(progs):
        import random
        keep = set(random.Random(20260927).sample(range(len(progs)), a.sample))
        progs = [pp if i in keep else None for i, pp in enumerate(progs)]
    json.dump({"n_programs": sum(x is not None for x in progs), "max_bits": a.max_bits, "min_bits": a.min_bits, "workers": a.workers,
               "null": a.null, "battery": a.battery, "started": time.strftime("%Y-%m-%d %H:%M:%S")},
              open(os.path.join(a.out, "manifest.json"), "w"), indent=1)
    print(f"{sum(x is not None for x in progs)} programs, {a.workers} workers -> {a.out}", flush=True)
    ps = [Process(target=worker, args=(w, a.workers, progs, a.out, a.null, a.battery)) for w in range(a.workers)]
    for p in ps: p.start()
    for p in ps: p.join()
    print("ALL DONE", time.strftime("%Y-%m-%d %H:%M:%S"), flush=True)


if __name__ == "__main__":
    main()
