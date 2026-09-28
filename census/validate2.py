"""VALIDATION GATE v2 (owner requirement 2026-09-27: error rates must be established before scaling, because even 1%
wrong over ~1M programs is a disaster). Criteria fixed here BEFORE the first run; nothing scales until PASS.

  python -m census.validate2 [--workers 4] [--out census/validation/v2] [--n-null 300] [--n-test 300]

What is tested = the census pipeline exactly: simulate (3 seeds, 1000 steps) -> screen -> fingerprint -> null region
-> naming by the validated battery on a LONG re-run (5000 steps) of positives.

Changes from v1 (v1 FAILED; see VALIDATION v1):
- NON-CIRCULAR NEGATIVES. v1 put the scrambled runs in the null set and then scored them as negatives (they passed
  by construction). v2 builds the null region only from (a) interaction-free programs and (b) scrambled runs of
  NULL-SET programs, and scores a DISJOINT held-out set of scrambled runs of other programs.
- Scrambled runs are negatives only for their SPATIAL views (grid, agents, network adjacency, fields): scrambling
  destroys spatial structure but not every non-spatial collective effect (e.g. mean-field consensus, all-to-all sync),
  so phase and transient views of scrambled runs are not scored.
- SIZE: >= 300 held-out negatives, so 0 false flags bounds the false-flag rate at <= 1% (rule of three, 95%).
- Positives are parameter FAMILIES of each known phenomenon inside the regime where it is known to occur.

PASS criteria (all must hold):
  N1  0 of the held-out negatives flagged (flagged = emergent in >= 2 of 3 seeds AND outside the null region).
  N2  0 of the verified-disordered regime negatives flagged.
  P1  >= 90% of positive variants flagged.
  P3  0 wrong names (any positive named with a catalog pattern other than its expected one).
  P2  (reported, target >= 80%) positives with an expected pattern named with it (>= 2 of 3 seeds, >= confirmation).
"""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import argparse, collections, json, random, time
from multiprocessing import Pool
import numpy as np
from census import grammar_g1 as G
from census.triage import FAMILY, is_null_program

R = G.Rule
SPATIAL_VIEWS = {"C", "A", "F0", "F1"}                    # N (static graph) views carry node states only -> spatial
                                                          # on a graph; scored too
SCORED_SCRAMBLED = SPATIAL_VIEWS | {"N"}


def positive_families():
    """{family: ([programs], expected pattern or None)} — each variant inside the known regime."""
    fam = {}
    fam["voter consensus"] = ([G.make("C", {"C": k}, rules=[R("C.COPY", (("p", p),))]) for k in (2, 3) for p in (0.3, 0.1, 0.03)], "P18")
    fam["Schelling segregation"] = ([G.make("C", {"C": k}, rules=[R("C.SCHELL", (("th", t),))]) for k in (3, 4) for t in (0.3, 0.4, 0.5, 0.6)], "P1")
    fam["majority coarsening"] = ([G.make("C", {"C": k}, rules=[R("C.MAJ", (("p", p),))]) for k in (2, 3) for p in (1.0, 0.3, 0.1)], None)
    fam["cyclic-CA spirals/waves"] = ([G.make("C", {"C": k}, rules=[R("C.CYCLE", (("th", t), ("p", 1.0)))]) for k in (3, 4) for t in (1, 2, 3)], None)
    fam["BTW sandpile"] = ([G.make("C", {"C": 4}, rules=[R("C.SAND", (("p", p),))]) for p in (0.0003, 0.001, 0.003)], "P14")
    fam["spatial PD (Nowak-May)"] = ([G.make("C", {"C": 2}, rules=[R("C.GAME", (("g", (1.9, 0.0)),))])], "P27")
    fam["network sync (Kuramoto)"] = ([G.make("N", {"N": 1}, rules=[R("N.KURA", (("K", K),))]) for K in (1.0, 0.5, 0.3)], "P9")
    fam["lattice phase waves (Kuramoto)"] = ([G.make("C", {"C": 1}, rules=[R("C.KURA", (("K", K),))]) for K in (1.0, 0.5, 0.3)], None)
    fam["co-evolving voter fragmentation"] = ([G.make("N", {"N": 2}, rules=[R("N.COPY", (("p", 1.0),)), R("N.REWIRE", (("p", p),))]) for p in (1.0, 0.5)], "P34")
    fam["flocking (Vicsek)"] = ([G.make("A", {"A": 1}, agent=(v, e), rules=[R("A.ALIGN", (("r", r),))])
                                 for v in (1.0, 0.3) for e in (0.1, 0.3) for r in (1.0, 2.0)], "P5")
    fam["network voter consensus"] = ([G.make("N", {"N": 2}, rules=[R("N.COPY", (("p", p),))]) for p in (1.0, 0.3)], None)
    return fam


REGIME_NEG = {
    "near-free walkers": G.make("A", {"A": 1}, agent=(1.0, 3.14159), rules=[R("A.REPEL", (("b", "any"), ("r", 0.5)))]),
    "sync at tiny coupling (network)": G.make("N", {"N": 1}, rules=[R("N.KURA", (("K", 0.0003),))]),
    "neighbour swaps (mixing)": G.make("C", {"C": 2}, rules=[R("C.SWAP", (("p", 1.0),))]),
    "neighbour swaps, 3 types": G.make("C", {"C": 3}, rules=[R("C.SWAP", (("p", 0.3),))]),
}


def _eval(job):
    """job = (tag, role, bits, shuffle, name_it) -> per-view flags + fingerprints (+ long-run battery names)."""
    tag, role, bits, shuffle, name_it = job
    from census import sim, filter as FL
    from census.runner import seed_of
    p = G.decode(bits); t = time.time()
    out = sim.run(p, seed0=seed_of(bits), null_shuffle=shuffle); ev = FL.evaluate(out, p, battery=False); views = {}
    for v, e in ev.items():
        fps = [sd.get("fp") or {} for sd in e["seeds"] if sd.get("emergent")]
        keys = sorted(set().union(*fps)) if fps else []
        views[v] = {"emergent_seeds": e["emergent_seeds"], "fp": {k: float(np.mean([f.get(k, 0.0) for f in fps])) for k in keys},
                    "em": [sd.get("em_score") for sd in e["seeds"]]}
    names = {}
    if name_it:
        from epc.phase2a.panel import _detected, _verdict
        longo = sim.run(p, seed0=seed_of(bits), steps=5000, rec_scale=5); B, FULL = FL._load()
        for s in range(3):
            for vn, h, md in FL.views(longo, p, s):
                h2 = FL._thin(h, FL.CONFIRM_FRAMES); hits = []
                for pid, fn in FULL.items():
                    try:
                        res = fn(h2, md)
                        if _detected(res) and FL._TIER.get(_verdict(res), 0) >= 2: hits.append((FL._TIER[_verdict(res)], pid))
                    except Exception:
                        pass
                names.setdefault(vn, []).append(max(hits)[1] if hits else None)
    return {"tag": tag, "role": role, "shuffle": shuffle, "prog": G.describe(p), "layers": p.layers, "views": views,
            "names": names, "seconds": round(time.time() - t, 1)}


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--out", default="census/validation/v2"); ap.add_argument("--n-null", type=int, default=300)
    ap.add_argument("--n-test", type=int, default=300); ap.add_argument("--maxbits", type=int, default=18)
    a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
    progs = [p for _, p in G.enumerate_programs(a.maxbits)]
    free = [p for p in progs if is_null_program(G.describe(p))]
    inter = [p for p in progs if not is_null_program(G.describe(p))]
    rng = random.Random(20260928); rng.shuffle(inter)
    null_set, test_set = inter[:a.n_null], inter[a.n_null:a.n_null + a.n_test]
    jobs = [(f"free:{i}", "null-free", G.encode(p), False, False) for i, p in enumerate(free)]
    jobs += [(f"nullscr:{i}", "null-scrambled", G.encode(p), True, False) for i, p in enumerate(null_set)]
    jobs += [(f"testscr:{i}", "neg-scrambled", G.encode(p), True, False) for i, p in enumerate(test_set)]
    jobs += [(n, "neg-regime", G.encode(p), False, False) for n, p in REGIME_NEG.items()]
    fams = positive_families()
    jobs += [(f"{f}#{i}", "pos", G.encode(p), False, True) for f, (ps, _) in fams.items() for i, p in enumerate(ps)]
    print(f"{len(free)} interaction-free, {len(null_set)} null-set scrambled, {len(test_set)} held-out scrambled, "
          f"{len(REGIME_NEG)} regime negatives, {sum(len(ps) for ps, _ in fams.values())} positive variants", flush=True)
    res = []
    with Pool(a.workers) as pool:
        for r in pool.imap_unordered(_eval, jobs):
            res.append(r)
            if len(res) % 50 == 0: print(f"{len(res)}/{len(jobs)}", flush=True)
    json.dump(res, open(os.path.join(a.out, "results.json"), "w"), indent=1, default=str)

    # null region per family from interaction-free + null-set scrambled only (never from scored cases)
    nulls = collections.defaultdict(list); scored = collections.defaultdict(list)
    for r in res:
        for v, d in r["views"].items():
            if d["emergent_seeds"] < 2 or not d["fp"]: continue
            fam = FAMILY.get(v, v)
            if r["role"] in ("null-free", "null-scrambled"): nulls[fam].append(d["fp"])
            else: scored[fam].append((r, v, d["fp"]))
    trivial = set(); Rs = {}
    for fam, items in scored.items():
        allf = nulls.get(fam, []) + [f for _, _, f in items]; keys = sorted(set().union(*allf))
        X = np.array([[f.get(k, 0.0) for k in keys] for f in allf])
        med = np.median(X, 0); mad = np.median(np.abs(X - med), 0) * 1.4826; mad[mad < 1e-9] = 1.0
        Z = np.clip((X - med) / mad, -8, 8); Nz, Sz = Z[:len(nulls.get(fam, []))], Z[len(nulls.get(fam, [])):]
        if len(Nz) >= 2:
            dnn = np.sqrt(((Nz[:, None] - Nz[None]) ** 2).sum(-1)); np.fill_diagonal(dnn, np.inf); Rr = float(np.percentile(dnn.min(1), 95))
            dmin = np.sqrt(((Sz[:, None] - Nz[None]) ** 2).sum(-1)).min(1)
        else:
            Rr, dmin = None, np.full(len(Sz), np.inf)
        Rs[fam] = {"null_views": len(Nz), "R": Rr}
        for i, (r, v, _) in enumerate(items):
            if Rr is not None and dmin[i] <= Rr: trivial.add((r["tag"], v))

    def flagged(r, only=None):
        return [v for v, d in r["views"].items() if d["emergent_seeds"] >= 2 and (r["tag"], v) not in trivial
                and (only is None or v in only)]
    negs = [r for r in res if r["role"] == "neg-scrambled"]; regs = [r for r in res if r["role"] == "neg-regime"]
    pos = [r for r in res if r["role"] == "pos"]
    n1 = [r for r in negs if flagged(r, SCORED_SCRAMBLED)]; n2 = [r for r in regs if flagged(r)]
    p1 = [r for r in pos if flagged(r)]
    exp_of = {f: e for f, (_, e) in fams.items()}
    p2n = p2d = 0; wrong = []
    fam_rows = collections.defaultdict(lambda: [0, 0, 0, 0])            # variants, flagged, named-right, expected
    for r in pos:
        f = r["tag"].split("#")[0]; e = exp_of[f]; fr = fam_rows[f]; fr[0] += 1; fr[1] += bool(flagged(r))
        got = {n for v, ns in r["names"].items() for n in set(ns) if n and ns.count(n) >= 2}
        if e: p2d += 1; fr[3] += 1
        if e and e in got: p2n += 1; fr[2] += 1
        bad = got - ({e} if e else set())
        if bad: wrong.append((r["tag"], sorted(bad)))
    crit = {f"N1 held-out scrambled negatives flagged = {len(n1)} / {len(negs)} (need 0)": len(n1) == 0,
            f"N2 regime negatives flagged = {len(n2)} / {len(regs)} (need 0)": len(n2) == 0,
            f"P1 positives flagged = {len(p1)} / {len(pos)} (need >= 90%)": len(p1) >= 0.9 * len(pos),
            f"P3 wrong names = {len(wrong)} (need 0)": len(wrong) == 0}
    L = ["# Validation gate v2 — " + time.strftime("%Y-%m-%d %H:%M"), "",
         "Null region (built only from interaction-free programs + scrambled runs of the NULL set): " + json.dumps(Rs), "",
         "## Positive families", "", "| family | variants | flagged | named correctly | expected name |", "|---|---|---|---|---|"]
    for f, (nv, nf, nr, ne) in fam_rows.items():
        L.append(f"| {f} | {nv} | {nf} | {nr if ne else '-'}{'/' + str(ne) if ne else ''} | {exp_of[f] or '— (no detector)'} |")
    L += ["", "## Criteria", ""] + [f"- {'PASS' if v else 'FAIL'} — {k}" for k, v in crit.items()]
    L += [f"- (reported) P2 named correctly = {p2n} / {p2d}"]
    if n1: L += ["", "### Held-out negatives flagged"] + [f"- `{r['prog'][:140]}` views {flagged(r, SCORED_SCRAMBLED)}" for r in n1[:30]]
    if n2: L += ["", "### Regime negatives flagged"] + [f"- {r['tag']} views {flagged(r)}" for r in n2]
    if wrong: L += ["", "### Wrong names"] + [f"- {t}: {b}" for t, b in wrong]
    miss = [r for r in pos if not flagged(r)]
    if miss: L += ["", "### Positives not flagged"] + [f"- {r['tag']}: `{r['prog'][:120]}`" for r in miss]
    ub = 3.0 / max(len(negs), 1)
    L += ["", f"**GATE: {'PASS' if all(crit.values()) else 'FAIL'}** — false-flag rate upper bound (95%, rule of three) "
          f"{'<= %.1f%%' % (100 * ub) if not n1 else 'n/a (flags observed)'}"]
    open(os.path.join(a.out, "VALIDATION.md"), "w").write("\n".join(L) + "\n"); print("\n".join(L))


if __name__ == "__main__":
    main()
