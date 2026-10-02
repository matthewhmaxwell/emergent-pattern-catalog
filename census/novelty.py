"""Novelty arm (pilot draft v0; DESIGN §5, ~35% of compute): CVT-MAP-Elites over COUPLED G1 programs.

  python -m census.novelty --out census/runs/<tag> --evals N [--workers 5] [--centroids 512] [--seed-from RUN]

Search space: G1 programs that couple at least two ingredients — two or more layers joined by at least one coupling
rule (agents<->cells, agents<->fields, cells<->fields), or a single layer with an adaptation rule (GAME) plus another
mechanism. Niches: CVT-MAP-Elites (Vassiliades et al. 2017) — centroids = k-means over the fingerprint descriptors of
random coupled programs (one CVT per layer combination, since descriptors differ by combination). Each niche keeps
its elite: highest quality = max emergence score over the program's views, ties -> shorter program. Variation:
parameter change, template swap, add rule (<= 5), drop rule, header tweak — always re-validated and re-coupled.
Every evaluation goes through the VALIDATED pipeline (screen + knock-out + textbook measures); only FLAGGED programs can
become elites. Outputs results_0.jsonl (runner-format rows, so census.triage names them like census programs) +
elites.json.
"""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import argparse, json, random, time
from multiprocessing import Pool
import numpy as np
from census import grammar_g1 as G

COUPLING = {"A.TURNCELL", "A.WRITECELL", "A.DEPOSIT", "A.CLIMB", "C.EMIT", "C.SENSE"}
ADAPT = {"C.GAME", "N.GAME"}
MULTI = ["CA", "CF", "AF", "CAF"]
MAX_RULES = 5


def coupled(p):
    names = {r.tmpl for r in p.rules}
    if len(p.layers) >= 2: return bool(names & COUPLING)
    return bool(names & ADAPT) and len(p.rules) >= 2


def _rand_rule(p, rng, prefer=None):
    vt = [t for t in G.valid_templates(p) if prefer is None or t in prefer] or G.valid_templates(p)
    for _ in range(50):
        t = rng.choice(vt)
        prm = tuple((n, rng.choice(G._param_options(k, p))[0]) for n, k in G.T[t]["params"])
        if G._ok_types(t, prm): return G.Rule(t, prm)
    return None


def _rand_header(rng):
    combo = rng.choice(MULTI + ["C", "N"])
    k = {L: rng.randint(1, 4) for L in "CAN" if L in combo}
    agent = (rng.choice([v for v, _ in G.hier_codes(G.GRIDS["speed"])]), rng.choice([v for v, _ in G.hier_codes(G.GRIDS["angle"])])) if "A" in combo else ()
    nf = rng.choice([1, 2]) if "F" in combo else 0
    D = tuple(rng.choice([v for v, _ in G.hier_codes(G.GRIDS["diff"])]) for _ in range(nf))
    return G.make(combo, k, nf, D, agent)


def _repair(p, rng):
    """drop invalid rules, make sure the program is coupled; None if impossible."""
    vt = set(G.valid_templates(p))
    rules = [r for r in p.rules if r.tmpl in vt and G._ok_types(r.tmpl, r.params)
             and all(v in [x for x, _ in G._param_options(k, p)] for (n, v), (_, k) in zip(r.params, G.T[r.tmpl]["params"]))]
    q = G.make(p.layers, p.k, p.nf, p.D, p.agent, rules)
    for _ in range(5):
        if q.rules and coupled(q): return q
        need = COUPLING if len(q.layers) >= 2 else (ADAPT if not ({r.tmpl for r in q.rules} & ADAPT) else None)
        r = _rand_rule(q, rng, need)
        if r is None or len(q.rules) >= MAX_RULES: return None
        q = G.make(q.layers, q.k, q.nf, q.D, q.agent, list(q.rules) + [r])
    return q if q.rules and coupled(q) else None


def random_program(rng):
    for _ in range(100):
        q = _repair(_rand_header(rng), rng)
        if q is not None: return q
    raise RuntimeError("could not build a coupled program")


def mutate(p, rng):
    for _ in range(20):
        op = rng.choice(["param", "param", "swap", "add", "drop", "header"])
        rules = list(p.rules); q = p
        if op == "param" and rules:
            i = rng.randrange(len(rules)); r = rules[i]; ps = G.T[r.tmpl]["params"]
            if not ps: continue
            j = rng.randrange(len(ps)); n, k = ps[j]
            prm = list(r.params); prm[j] = (n, rng.choice(G._param_options(k, p))[0]); rules[i] = G.Rule(r.tmpl, tuple(prm))
            q = G.make(p.layers, p.k, p.nf, p.D, p.agent, rules)
        elif op == "swap" and rules:
            r = _rand_rule(p, rng)
            if r is None: continue
            rules[rng.randrange(len(rules))] = r; q = G.make(p.layers, p.k, p.nf, p.D, p.agent, rules)
        elif op == "add" and len(rules) < MAX_RULES:
            r = _rand_rule(p, rng)
            if r is None: continue
            rules.insert(rng.randint(0, len(rules)), r); q = G.make(p.layers, p.k, p.nf, p.D, p.agent, rules)
        elif op == "drop" and len(rules) > 1:
            rules.pop(rng.randrange(len(rules))); q = G.make(p.layers, p.k, p.nf, p.D, p.agent, rules)
        elif op == "header":
            k = dict(p.k)
            if k and rng.random() < 0.5:
                L = rng.choice(list(k)); k[L] = rng.randint(1, 4)
            agent = p.agent
            if agent and rng.random() < 0.5:
                agent = (rng.choice([v for v, _ in G.hier_codes(G.GRIDS["speed"])]), rng.choice([v for v, _ in G.hier_codes(G.GRIDS["angle"])]))
            q = G.make(p.layers, k, p.nf, p.D, agent, rules)
        else:
            continue
        q = _repair(q, rng)
        if q is not None and G.encode(q) != G.encode(p): return q
    return random_program(rng)


def evaluate(bits):
    """one program through the VALIDATED pipeline (census.runner.evaluate). Row format = runner row (so census.triage
    names novelty-arm programs exactly like census programs) + quality + per-view mean fingerprint."""
    from census.runner import evaluate as run_eval
    p = G.decode(bits); t = time.time(); views, status, info = run_eval(p, bits); fpm = {}; q = 0.0
    for v, e in views.items():
        if not e["flagged"]: continue
        ss = [s for s in e["seeds"] if s.get("emergent") and "fp" in s]; fps = [s["fp"] for s in ss]
        keys = sorted(set().union(*fps)) if fps else []
        fpm[v] = {k: float(np.mean([f.get(k, 0.0) for f in fps])) for k in keys}
        q = max(q, float(np.mean([s["evidence"] for s in ss])) if ss else 0.0)
    return {"bits": bits, "len": len(bits), "prog": G.describe(p), "layers": p.layers, "n_rules": len(p.rules),
            "templates": sorted({r.tmpl for r in p.rules}), "grammar": "g1", "status": status, "views": views,
            "quality": q, "fp_mean": fpm, "fp_errors": info["fp_errors"], "seconds": round(time.time() - t, 2)}


def descriptor(res, keys):
    """fixed-order descriptor: per-view fingerprint features (0 where a view is absent / not emergent)."""
    return np.array([res["fp_mean"].get(v, {}).get(f, 0.0) for v, f in keys], float)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--out", required=True); ap.add_argument("--evals", type=int, default=2000)
    ap.add_argument("--workers", type=int, default=5); ap.add_argument("--centroids", type=int, default=256)
    ap.add_argument("--init", type=int, default=200); ap.add_argument("--batch", type=int, default=20)
    ap.add_argument("--seed", type=int, default=0)
    ap.add_argument("--patience", type=int, default=20000, help="stop after this many evaluations with no NEW niche")
    a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
    rng = random.Random(a.seed); arch = open(os.path.join(a.out, "results_0.jsonl"), "a"); counter = [0]
    json.dump({"grammar": "g1", "n_programs": a.evals, "pipeline": "novelty arm (CVT-MAP-Elites) on the validated pipeline",
               "started": time.strftime("%Y-%m-%d %H:%M:%S")}, open(os.path.join(a.out, "manifest.json"), "w"), indent=1)

    def write(r):
        r["idx"] = counter[0]; counter[0] += 1; arch.write(json.dumps(r, default=str) + "\n")
    with Pool(a.workers) as pool:
        init = [G.encode(random_program(rng)) for _ in range(a.init)]
        res = [r for r in pool.imap_unordered(evaluate, init)]
        for r in res: write(r)
        arch.flush()
        # descriptor keys per layer combination from the initial random sample; CVT centroids by k-means
        from scipy.cluster.vq import kmeans2
        cvt = {}
        for combo in sorted({r["layers"] for r in res}):
            rs = [r for r in res if r["layers"] == combo and r["status"] == "FLAGGED"]
            keys = sorted({(v, f) for r in rs for v, d in r["fp_mean"].items() for f in d})
            if not keys or len(rs) < 4: continue
            X = np.array([descriptor(r, keys) for r in rs]); mu, sd = X.mean(0), X.std(0) + 1e-9
            k = min(a.centroids, max(2, len(rs) // 2))
            C, _ = kmeans2((X - mu) / sd, k, minit="++", seed=a.seed)
            cvt[combo] = {"keys": keys, "mu": mu, "sd": sd, "C": C}
        elites = {}

        def place(r):
            c = cvt.get(r["layers"])
            if c is None or r["status"] != "FLAGGED" or r["quality"] <= 0: return False
            z = (descriptor(r, c["keys"]) - c["mu"]) / c["sd"]; cell = (r["layers"], int(((c["C"] - z) ** 2).sum(1).argmin()))
            cur = elites.get(cell)
            if cur is None or (r["quality"], -r["len"]) > (cur["quality"], -cur["len"]):
                r["new_niche"] = cur is None; elites[cell] = r; return True
            return False
        for r in res: place(r)
        n = len(res); since_new = 0
        while n < a.evals and since_new < a.patience:
            if elites:
                parents = [G.decode(rng.choice(list(elites.values()))["bits"]) for _ in range(a.batch)]
                kids = [G.encode(mutate(p, rng)) for p in parents]
            else:
                kids = [G.encode(random_program(rng)) for _ in range(a.batch)]
            for r in pool.imap_unordered(evaluate, kids):
                r["new_elite"] = place(r); write(r); n += 1
                since_new = 0 if r.get("new_niche") else since_new + 1
            arch.flush()
            json.dump({f"{k[0]}:{k[1]}": {"quality": v["quality"], "len": v["len"], "prog": v["prog"], "bits": v["bits"]}
                       for k, v in elites.items()}, open(os.path.join(a.out, "elites.json"), "w"), indent=1)
            print(f"{n} evals, {len(elites)} niches filled, {since_new} since last new niche", flush=True)


if __name__ == "__main__":
    main()
