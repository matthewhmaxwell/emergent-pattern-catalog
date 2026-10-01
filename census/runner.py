"""Census runner: enumerate -> simulate -> screen -> interaction knock-out -> textbook measures. Resumable.

  python -m census.runner --out census/runs/<tag> --max-bits L [--min-bits l] [--workers 5] [--limit N] [--sample N]

This is the pipeline that passed the validation gate (confirm3, 2026-09-28): a view is FLAGGED when the screen says
emergent in >= 2 of 3 seeds AND the interaction knock-out says the pattern depends on interaction (census.knockout).
For flagged views the per-seed fingerprint and every applicable textbook measure are stored, so naming and grouping
(census.triage) run afterwards from the stored rows — the reference library can grow without re-simulating.

Every canonical G1 program of length in [min_bits, max_bits] gets a stable index (enumeration order) and a stable
seed (hash of its bit string), so any program can be re-simulated exactly later — raw trajectories are not stored.
Worker w processes indices with idx % workers == w and appends to results_w.jsonl; on restart, finished indices are
skipped (all workers' files are read, so the worker count may change). Run inside tmux as matthewhmaxwell, niced.
"""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
    os.environ.setdefault(_v, "1")                            # one thread per worker: workers are the parallelism
import argparse, hashlib, json, sys, time, traceback
from multiprocessing import Process

from census import grammar_g1 as G


def seed_of(bits):
    return int(hashlib.sha1(bits.encode()).hexdigest()[:8], 16)


def program_list(max_bits, cache_dir):
    """every canonical G1 program of length <= max_bits as bit strings, shortest first (then lexicographic).
    Enumeration is slow (hours at 24 bits), so the list is cached on disk and reused by every start / resume."""
    path = os.path.join(cache_dir, f"programs_g1_{max_bits}.txt")
    if os.path.exists(path):
        return [l.strip() for l in open(path) if l.strip()]
    bits = sorted((b for b, _ in G.enumerate_programs(max_bits)), key=lambda b: (len(b), b))
    os.makedirs(cache_dir, exist_ok=True); tmp = path + ".tmp"
    with open(tmp, "w") as fh: fh.write("\n".join(bits) + "\n")
    os.replace(tmp, path)
    return bits


def evaluate(p, bits):
    """one program -> (views dict, status, timings). The validated flag = emergent (>= 2 seeds) AND interaction-driven."""
    from census.knockout import run_with_knockout, is_interaction_free
    from census.reflib import applicable_checks
    t0 = time.time(); out, real, driven, ko_out = run_with_knockout(p, seed_of(bits)); t1 = time.time()
    views = {}; fperr = 0
    for v, seeds in real.items():
        n_em = sum(1 for x in seeds if x["emergent"])
        flagged = sum(1 for s, x in enumerate(seeds) if x["emergent"] and driven[v][s]) >= 2
        row = {"flagged": bool(flagged), "em_seeds": n_em, "driven": bool(driven[v][0]) if seeds else False, "seeds": []}
        checks = applicable_checks(p, v) if flagged else {}
        for s, x in enumerate(seeds):
            fperr += "fp_error" in x["fp"]
            sd = {"emergent": x["emergent"], "evidence": round(float(x["evidence"]), 4)}
            if flagged and x["emergent"]:
                sd["fp"] = {k: round(float(val), 5) for k, val in x["fp"].items()}; sd["checks"] = {}
                for c, fn in checks.items():
                    try: sd["checks"][c] = bool(fn(out, s)[0])
                    except Exception: sd["checks"][c] = False
            row["seeds"].append(sd)
        views[v] = row
    status = "FLAGGED" if any(x["flagged"] for x in views.values()) else "NONE"
    return views, status, {"sim_s": round(t1 - t0, 2), "post_s": round(time.time() - t1, 2), "knockout_run": ko_out is not None,
                           "interaction_free": is_interaction_free(p), "fp_errors": fperr,
                           "meta": {k: out["meta"][k] for k in ("steps_run", "frozen_at", "absorbed_at", "unstable")}}


def worker(w, nw, progs, outdir):
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
                   "n_rules": len(p.rules), "templates": sorted({r.tmpl for r in p.rules})}
            try:
                views, status, info = evaluate(p, bits); row.update(info); row["views"] = views; row["status"] = status
            except Exception as e:
                row.update({"status": "ERROR", "error": f"{type(e).__name__}: {e}", "trace": traceback.format_exc()[-800:]})
            fh.write(json.dumps(row, default=str) + "\n"); fh.flush()


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", required=True); ap.add_argument("--max-bits", type=int, required=True)
    ap.add_argument("--min-bits", type=int, default=0); ap.add_argument("--workers", type=int, default=5)
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--cache-dir", default="census/programs", help="where enumerated program lists are cached")
    ap.add_argument("--sample", type=int, default=0, help="fixed-seed random subsample of N programs (indices kept)")
    a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
    progs = [(b, G.decode(b)) for b in program_list(a.max_bits, a.cache_dir) if len(b) >= a.min_bits]
    if a.limit: progs = progs[:a.limit]
    if a.sample and a.sample < len(progs):
        import random
        keep = set(random.Random(20260927).sample(range(len(progs)), a.sample))
        progs = [pp if i in keep else None for i, pp in enumerate(progs)]
    json.dump({"n_programs": sum(x is not None for x in progs), "max_bits": a.max_bits, "min_bits": a.min_bits, "workers": a.workers,
               "pipeline": "validated (screen + knock-out + textbook measures)", "sample": a.sample,
               "started": time.strftime("%Y-%m-%d %H:%M:%S")},
              open(os.path.join(a.out, "manifest.json"), "w"), indent=1)
    print(f"{sum(x is not None for x in progs)} programs, {a.workers} workers -> {a.out}", flush=True)
    ps = [Process(target=worker, args=(w, a.workers, progs, a.out)) for w in range(a.workers)]
    for p in ps: p.start()
    for p in ps: p.join()
    print("ALL DONE", time.strftime("%Y-%m-%d %H:%M:%S"), flush=True)


if __name__ == "__main__":
    main()
