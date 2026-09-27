"""VALIDATION GATE (owner requirement 2026-09-27): before any census-scale run, the full pipeline must handle known
positives and known negatives correctly. Criteria are fixed here BEFORE the first run; nothing scales until PASS.

  python -m census.validate [--workers 5] [--out census/validation/<tag>]

Pipeline under test = exactly the census pipeline: simulate (3 seeds) -> screen -> fingerprint -> null region
(interaction-free programs + structure-scrambled runs) -> validated 36-detector battery (full frames) for naming.

POSITIVES (known emergent behaviour). `expect` = the catalog pattern that should be named, or None when the catalog
has no detector for it (then only "flagged" is required).
NEGATIVES (should NOT be flagged): interaction-free rules, known models in their disordered/dead regime, and every
positive re-run with its structure scrambled after each step (null_shuffle).

PASS criteria (all must hold):
  N1  0 negatives flagged (flagged = emergent in >= 2 of 3 seeds AND outside the null region).
  P1  >= 90% of positives flagged.
  P2  >= 80% of positives that have an expected pattern are named with it (>= 2 of 3 seeds, >= confirmation tier).
  P3  0 wrong names: no positive or negative is named with a catalog pattern other than its expected one.
"""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import argparse, collections, json, time
from multiprocessing import Pool
import numpy as np
from census import grammar_g1 as G

R = G.Rule
POS = {  # name: (program, expected pattern or None)
    "voter (p=0.3)":        (G.make("C", {"C": 2}, rules=[R("C.COPY", (("p", 0.3),))]), "P18"),
    "Schelling":            (G.make("C", {"C": 3}, rules=[R("C.SCHELL", (("th", 0.5),))]), "P1"),
    "majority coarsening":  (G.make("C", {"C": 2}, rules=[R("C.MAJ", (("p", 0.1),))]), None),
    "cyclic-CA spirals":    (G.make("C", {"C": 3}, rules=[R("C.CYCLE", (("th", 3), ("p", 1.0)))]), None),
    "excitable waves (GH)": (G.make("C", {"C": 4}, rules=[R("C.CYCLE", (("th", 1), ("p", 1.0)))]), "P13"),
    "BTW sandpile":         (G.make("C", {"C": 4}, rules=[R("C.SAND", (("p", 0.0003),))]), "P14"),
    "spatial PD (Nowak-May)": (G.make("C", {"C": 2}, rules=[R("C.GAME", (("g", (1.9, 0.0)),))]), "P27"),
    "lattice sync (Kuramoto)": (G.make("C", {"C": 1}, rules=[R("C.KURA", (("K", 1.0),))]), "P9"),
    "network sync (Kuramoto)": (G.make("N", {"N": 1}, rules=[R("N.KURA", (("K", 1.0),))]), "P9"),
    "co-evolving voter":    (G.make("N", {"N": 2}, rules=[R("N.COPY", (("p", 1.0),)), R("N.REWIRE", (("p", 1.0),))]), "P34"),
    "flocking (Vicsek)":    (G.make("A", {"A": 1}, agent=(1.0, 0.3), rules=[R("A.ALIGN", (("r", 3.0),))]), "P5"),
    "MIPS (quorum slow-down)": (G.make("A", {"A": 1}, agent=(1.0, 1.0), rules=[R("A.QUORUM", (("th", 3), ("r", 3.0)))]), "P2"),
    "Langton ants":         (G.make("CA", {"C": 2, "A": 1}, agent=(1.0, 0.0), rules=[R("A.TURNCELL"), R("A.WRITECELL")]), None),
}
NEG = {  # name: program (interaction-free or disordered regime)
    "random flips (lattice)":  G.make("C", {"C": 2}, rules=[R("C.FLIP", (("a", 0), ("b", 1), ("p", 0.1))), R("C.FLIP", (("a", 1), ("b", 0), ("p", 0.1)))]),
    "random flips (network)":  G.make("N", {"N": 2}, rules=[R("N.FLIP", (("a", 0), ("b", 1), ("p", 0.1))), R("N.FLIP", (("a", 1), ("b", 0), ("p", 0.1)))]),
    "neighbour swaps (mixing)": G.make("C", {"C": 2}, rules=[R("C.SWAP", (("p", 1.0),))]),
    "diffusion + decay only":  G.make("F", nf=1, D=(0.2,), rules=[R("F.DECAY", (("f", 0), ("d", 0.1)))]),
    "flocking at max noise":   G.make("A", {"A": 1}, agent=(1.0, 3.14159), rules=[R("A.ALIGN", (("r", 1.0),))]),
    "sync at tiny coupling (network)": G.make("N", {"N": 1}, rules=[R("N.KURA", (("K", 0.0003),))]),
    "near-free walkers":       G.make("A", {"A": 1}, agent=(1.0, 3.14159), rules=[R("A.REPEL", (("b", "any"), ("r", 0.5)))]),
}
INTERACTION_FREE = {"random flips (lattice)", "random flips (network)", "diffusion + decay only"}


def _run(job):
    """job = (name, role, bits, shuffle) -> per-view: flagged seeds, fingerprints, battery names per seed."""
    name, role, bits, shuffle = job
    from census import sim, filter as FL
    from census.runner import seed_of
    from epc.phase2a.panel import _detected, _verdict
    p = G.decode(bits); t = time.time(); out = sim.run(p, seed0=seed_of(bits), null_shuffle=shuffle)
    ev = FL.evaluate(out, p, battery=False); B, FULL = FL._load(); views = {}
    for v, e in ev.items():
        names = []
        for s, sd in enumerate(e["seeds"]):
            nm = None
            if sd.get("emergent"):
                h, md = {n: (h, md) for n, h, md in FL.views(out, p, s)}[v]; h = FL._thin(h, FL.CONFIRM_FRAMES); hits = []
                for pid, fn in FULL.items():
                    try:
                        res = fn(h, md)
                        if _detected(res) and FL._TIER.get(_verdict(res), 0) >= 2: hits.append((FL._TIER[_verdict(res)], pid))
                    except Exception:
                        pass
                nm = max(hits)[1] if hits else None
            names.append(nm)
        fps = [sd.get("fp") or {} for sd in e["seeds"] if sd.get("emergent")]
        keys = sorted(set().union(*fps)) if fps else []
        views[v] = {"emergent_seeds": e["emergent_seeds"], "names": names,
                    "em": [sd.get("em_score") for sd in e["seeds"]], "kinds": [sd.get("em_kind") for sd in e["seeds"]],
                    "fp": {k: float(np.mean([f.get(k, 0.0) for f in fps])) for k in keys}}
    return {"name": name, "role": role, "shuffle": shuffle, "prog": G.describe(p), "views": views,
            "seconds": round(time.time() - t, 1)}


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--workers", type=int, default=5)
    ap.add_argument("--out", default="census/validation/v1"); a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
    jobs = [(n, "pos", G.encode(p), False) for n, (p, _) in POS.items()]
    jobs += [(n, "neg", G.encode(p), False) for n, p in NEG.items()]
    jobs += [(n + " [scrambled]", "neg", G.encode(p), True) for n, (p, _) in POS.items()]
    with Pool(a.workers) as pool:
        res = []
        for r in pool.imap_unordered(_run, jobs):
            res.append(r); print(f"done {r['name']} ({r['seconds']}s)", flush=True)
    json.dump(res, open(os.path.join(a.out, "results.json"), "w"), indent=1, default=str)

    # null region per family: interaction-free negatives + scrambled runs (same rule as census triage)
    from census.triage import FAMILY
    fam_views = collections.defaultdict(list)
    for r in res:
        for v, d in r["views"].items():
            if d["emergent_seeds"] >= 2 and d["fp"]: fam_views[FAMILY.get(v, v)].append((r, v, d))
    trivial = set()
    for fam, items in fam_views.items():
        keys = sorted(set().union(*[d["fp"] for _, _, d in items]))
        X = np.array([[d["fp"].get(k, 0.0) for k in keys] for _, _, d in items])
        med = np.median(X, 0); mad = np.median(np.abs(X - med), 0) * 1.4826; mad[mad < 1e-9] = 1.0
        Z = np.clip((X - med) / mad, -8, 8)
        isnull = np.array([r["shuffle"] or r["name"] in INTERACTION_FREE for r, _, _ in items])
        if isnull.sum() == 0: continue
        Nz = Z[isnull]
        if len(Nz) >= 2:
            dnn = np.sqrt(((Nz[:, None] - Nz[None]) ** 2).sum(-1)); np.fill_diagonal(dnn, np.inf); Rr = np.percentile(dnn.min(1), 95)
        else:
            Rr = 0.0
        d = np.sqrt(((Z[:, None] - Nz[None]) ** 2).sum(-1)).min(1)
        for i, (r, v, _) in enumerate(items):
            if isnull[i] or d[i] <= Rr: trivial.add((r["name"], v))

    def outcome(r):
        flagged_views = [v for v, d in r["views"].items() if d["emergent_seeds"] >= 2 and (r["name"], v) not in trivial]
        named = collections.Counter(n for v in r["views"].values() for n in v["names"] if n)
        names2 = {n for v in r["views"].values() for n in set(v["names"]) if n and v["names"].count(n) >= 2}
        return flagged_views, names2, named
    L = ["# Validation gate — " + time.strftime("%Y-%m-%d %H:%M"), "", "| case | role | flagged views | named (>=2 seeds) | expected | ok |", "|---|---|---|---|---|---|"]
    n1 = p1 = p2 = p2d = p3 = 0; npos = sum(1 for r in res if r["role"] == "pos")
    for r in sorted(res, key=lambda r: (r["role"] != "pos", r["name"])):
        fl, nm, _ = outcome(r)
        if r["role"] == "pos":
            exp = POS[r["name"]][1]; ok_flag = bool(fl); p1 += ok_flag
            if exp: p2d += 1; p2 += exp in nm
            wrong = nm - ({exp} if exp else set()); p3 += len(wrong)
            ok = ok_flag and (exp is None or exp in nm) and not wrong
        else:
            exp = "none"; wrong = nm; p3 += len(wrong); n1 += bool(fl); ok = not fl and not wrong
        L.append(f"| {r['name']} | {r['role']} | {', '.join(fl) or '-'} | {', '.join(sorted(nm)) or '-'} | {exp} | {'PASS' if ok else 'FAIL'} |")
    crit = {"N1 negatives flagged == 0": n1 == 0, f"P1 positives flagged >= 90% ({p1}/{npos})": p1 >= 0.9 * npos,
            f"P2 named correctly >= 80% ({p2}/{p2d})": p2 >= 0.8 * p2d, f"P3 wrong names == 0 ({p3})": p3 == 0}
    L += ["", "## Criteria", ""] + [f"- {'PASS' if v else 'FAIL'} — {k}" for k, v in crit.items()]
    L += ["", f"**GATE: {'PASS' if all(crit.values()) else 'FAIL'}** (negatives flagged: {n1})"]
    open(os.path.join(a.out, "VALIDATION.md"), "w").write("\n".join(L) + "\n"); print("\n".join(L))


if __name__ == "__main__":
    main()
