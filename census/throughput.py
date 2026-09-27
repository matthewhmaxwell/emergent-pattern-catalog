"""Pilot T6/T7: throughput + null base rate on a random sample of G1 programs (pilot only; not census data).

  python -m census.throughput SAMPLE MAXBITS [--null]
Samples SAMPLE programs uniformly from all canonical G1 programs of length <= MAXBITS, runs sim + screen (+ known
filter on emergent views), and reports per-substrate timing and verdict rates. --null runs every sampled program in
the mean-field shuffle null instead."""
import sys, time, json, random, collections
from census import grammar_g1 as G, sim, filter as FL

n, maxbits, null = int(sys.argv[1]), int(sys.argv[2]), "--null" in sys.argv
progs = [p for _, p in G.enumerate_programs(maxbits)]
random.seed(12345); sample = random.sample(progs, min(n, len(progs)))
print(f"{len(progs)} canonical programs <= {maxbits} bits; sampling {len(sample)}; null={null}", flush=True)
rows = []
for i, p in enumerate(sample):
    t0 = time.time(); out = sim.run(p, seed0=i, null_shuffle=null); t1 = time.time()
    ev = FL.evaluate(out, p); t2 = time.time()
    verd = {v: (e["verdict"], e["pattern"], e["emergent_seeds"]) for v, e in ev.items()}
    rows.append({"prog": G.describe(p), "bits": len(G.encode(p)), "layers": p.layers, "sim_s": t1 - t0,
                 "filter_s": t2 - t1, "steps": out["meta"]["steps_run"], "frozen": out["meta"]["frozen_at"],
                 "unstable": out["meta"]["unstable"], "views": verd})
    print(f"[{i+1}/{len(sample)}] {t1-t0:5.1f}s + {t2-t1:5.1f}s  {verd}  {G.describe(p)}", flush=True)
by = collections.defaultdict(list)
for r in rows: by[r["layers"]].append(r)
summ = {}
for L, rs in sorted(by.items()):
    vs = [v for r in rs for v in r["views"].values()]
    summ[L] = {"n": len(rs), "mean_sim_s": sum(r["sim_s"] for r in rs) / len(rs),
               "mean_filter_s": sum(r["filter_s"] for r in rs) / len(rs),
               "frozen_frac": sum(r["frozen"] is not None for r in rs) / len(rs),
               "view_verdicts": dict(collections.Counter(v[0] for v in vs)),
               "matches": dict(collections.Counter(v[1] for v in vs if v[0] == "MATCH"))}
print(json.dumps(summ, indent=1))
tag = "null" if null else "census"
json.dump({"summary": summ, "rows": rows}, open(f"census/pilot_throughput_{tag}_{maxbits}b.json", "w"), indent=1, default=str)
