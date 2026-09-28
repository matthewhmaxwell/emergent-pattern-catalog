"""VALIDATION GATE (round 5) — screen + null region + REFERENCE-LIBRARY naming. Criteria fixed BEFORE the run.

  python -m census.validate3 --lib census/reflib/v1/library.json [--workers 4] [--out census/validation/v3]

Positives = the textbook-VERIFIED library examples (census.reflib). Naming accuracy is measured LEAVE-ONE-VARIANT-OUT:
each parameter setting (all its seeds) is removed from the library and then named by the rest, so no example is
ever named by itself or by its own seeds.
Negatives (non-circular): the null region is built ONLY from interaction-free programs (v3 lesson: scrambled runs
still flock / sync / reach consensus in mean field, so they must not define "noise"). Scored: held-out scrambled runs
(grid/field views only) + REGIME negatives for agents, networks and oscillators, each verified disordered by its
textbook measure (polarization, order parameter r, majority share).

NAMING RULE (round 5, stricter): per view family, robust-z fingerprints; k = 5 nearest library examples; name =
class c only if ALL 5 are class c AND the nearest class-c example is within R_c (90th percentile of class-c
nearest-other-variant distances). Classes with < 3 parameter variants in the library are never auto-named (they
cannot be validated leave-one-variant-out). Otherwise UNNAMED (-> literature check).

PASS criteria (all must hold):
  N1  0 of >= 300 HELD-OUT interaction-free programs flagged (disjoint from the 400 that define the null region).
      Round-4 lesson: scrambled runs are NOT negatives — they still show population-level collective dynamics
      (type-fraction oscillations, reactions, mean-field alignment/sync/consensus), which is real emergence, just
      not spatial. They are dropped from scoring and reported separately as a diagnostic.
  N2  0 textbook-verified regime negatives flagged; verification now needs BOTH global AND local disorder
      (round-4 lesson: noisy aligners with low global polarization still formed local aligned clumps).
  N3  0 negatives NAMED (a flagged negative that would also get a catalog name).
  P1  >= 90% of verified library examples flagged (emergent in the screen AND outside the null region).
  P3  0 WRONG names in leave-one-variant-out (named, but with another class).
  P2  (reported; target >= 80%) flagged verified examples named correctly.
"""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import argparse, collections, json, random, time
from multiprocessing import Pool
import numpy as np
from census import grammar_g1 as G
from census.triage import FAMILY, is_null_program
from census.validate2 import _eval

SCORED_SCRAMBLED = {"C", "F0", "F1"}
R_ = G.Rule
REGIME_NEG = {}                                        # textbook-verified disordered programs (checked in _disorder)
for v_ in (1.0, 0.3):
    for r_ in (0.5, 1.0):
        REGIME_NEG[f"walkers repel r={r_} v={v_}"] = G.make("A", {"A": 1}, agent=(v_, 3.14159), rules=[R_("A.REPEL", (("b", "any"), ("r", r_)))])
        REGIME_NEG[f"align at max noise r={r_} v={v_}"] = G.make("A", {"A": 1}, agent=(v_, 3.14159), rules=[R_("A.ALIGN", (("r", r_),))])
for K_ in (0.0003, 0.001, 0.002):
    REGIME_NEG[f"network oscillators K={K_}"] = G.make("N", {"N": 1}, rules=[R_("N.KURA", (("K", K_),))])
    REGIME_NEG[f"lattice oscillators K={K_}"] = G.make("C", {"C": 1}, rules=[R_("C.KURA", (("K", K_),))])
REGIME_NEG["neighbour swaps"] = G.make("C", {"C": 2}, rules=[R_("C.SWAP", (("p", 1.0),))])
REGIME_NEG["neighbour swaps, 3 types"] = G.make("C", {"C": 3}, rules=[R_("C.SWAP", (("p", 0.3),))])


def disorder_ok(tag, res_row):
    """textbook disorder check, from the run's own fingerprint values (seed means)."""
    fp = {k: v for d in res_row["views"].values() for k, v in d["fp"].items()}
    if tag.startswith(("walkers", "align")):
        return fp.get("a_polar", 0.0) < 0.2 and fp.get("a_local_align", 1.0) < 0.05 and fp.get("a_nn_ratio", 1.0) > 0.8
    if tag.startswith("network oscillators"): return fp.get("p_r", 0.0) < 0.3 and fp.get("p_r_std", 1.0) < 0.1
    if tag.startswith("lattice oscillators"): return fp.get("p_local_r", 0.0) < 0.6
    return True

K, VOTES, RPCT, MIN_VARIANTS = 5, 5, 90, 3


def zspace(fps_lists):
    keys = sorted(set().union(*[f for lst in fps_lists for f in lst]))
    X = np.array([[f.get(k, 0.0) for k in keys] for lst in fps_lists for f in lst], float)
    med = np.median(X, 0); mad = np.median(np.abs(X - med), 0) * 1.4826; mad[mad < 1e-9] = 1.0
    return keys, med, mad


def Z(fps, keys, med, mad):
    X = np.array([[f.get(k, 0.0) for k in keys] for f in fps], float).reshape(len(fps), len(keys))
    return np.clip((X - med) / mad, -8, 8)


class Namer:
    def __init__(self, lib_rows, keys, med, mad):
        self.keys, self.med, self.mad = keys, med, mad
        self.cls = np.array([r["class"] for r in lib_rows]); self.var = np.array([f"{r['class']}#{r['variant']}" for r in lib_rows])
        self.Z = Z([r["fp"] for r in lib_rows], keys, med, mad); self.R = {}
        for c in set(self.cls):                       # R_c: nearest OTHER-variant same-class distance, 95th pct
            idx = np.flatnonzero(self.cls == c); d = []
            for i in idx:
                o = idx[self.var[idx] != self.var[i]]
                if len(o): d.append(np.sqrt(((self.Z[o] - self.Z[i]) ** 2).sum(1)).min())
            self.R[c] = float(np.percentile(d, RPCT)) if d else 0.0
        self.nameable = {c for c in set(self.cls) if len(set(self.var[self.cls == c])) >= MIN_VARIANTS}

    def name(self, fp, exclude_variant=None):
        z = Z([fp], self.keys, self.med, self.mad)[0]; m = np.ones(len(self.cls), bool)
        if exclude_variant is not None: m &= self.var != exclude_variant
        if m.sum() < K: return None, None
        d = np.sqrt(((self.Z[m] - z) ** 2).sum(1)); cl = self.cls[m]; nn = np.argsort(d)[:K]
        top, cnt = collections.Counter(cl[nn]).most_common(1)[0]
        dc = d[cl == top].min()
        return (top if cnt >= VOTES and dc <= self.R[top] and top in self.nameable else None), float(dc)


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--lib", required=True); ap.add_argument("--workers", type=int, default=4)
    ap.add_argument("--out", default="census/validation/v3"); ap.add_argument("--n-null", type=int, default=300)
    ap.add_argument("--n-test", type=int, default=400); ap.add_argument("--maxbits", type=int, default=18)
    a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
    lib = json.load(open(a.lib)); ex = [r for r in lib["examples"] if r["verified"]]
    lib_bits = {r["bits"] for r in lib["examples"]}
    progs = [p for _, p in G.enumerate_programs(a.maxbits)]
    free = [p for p in progs if is_null_program(G.describe(p))]
    inter = [p for p in progs if not is_null_program(G.describe(p)) and G.encode(p) not in lib_bits]
    rng = random.Random(20260928); rng.shuffle(inter); rng.shuffle(free)
    free, free_test = free[:400], free[400:400 + a.n_test]   # null region / held-out scored (disjoint)
    null_set, test_set = [], inter[:150]                  # scrambled: diagnostic only (not scored)
    jobs = [(f"free:{i}", "null-free", G.encode(p), False, False) for i, p in enumerate(free)]
    jobs += [(f"freetest:{i}", "neg-free", G.encode(p), False, False) for i, p in enumerate(free_test)]
    jobs += [(f"testscr:{i}", "diag-scrambled", G.encode(p), True, False) for i, p in enumerate(test_set)]
    jobs += [(n, "neg-regime", G.encode(p), False, False) for n, p in REGIME_NEG.items()]
    print(f"negatives: {len(free)} interaction-free + {len(null_set)} null-set scrambled (null region); scored: "
          f"{len(test_set)} held-out scrambled + {len(REGIME_NEG)} regime | positives: {len(ex)} verified library examples", flush=True)
    res = []
    with Pool(a.workers) as pool:
        for r in pool.imap_unordered(_eval, jobs):
            res.append(r)
            if len(res) % 100 == 0: print(f"{len(res)}/{len(jobs)}", flush=True)
    json.dump(res, open(os.path.join(a.out, "negatives.json"), "w"), indent=1, default=str)

    # per family: null views (never scored), scored negative views, library examples
    nullfp = collections.defaultdict(list); negv = collections.defaultdict(list); diag = collections.defaultdict(list)
    for r in res:
        for v, d in r["views"].items():
            if d["emergent_seeds"] < 2 or not d["fp"]: continue
            fam = FAMILY.get(v, v)
            if r["role"] in ("null-free", "null-scrambled"): nullfp[fam].append(d["fp"])
            elif (r["role"] == "neg-regime" and disorder_ok(r["tag"], r)) or r["role"] == "neg-free":
                negv[fam].append((r, v, d["fp"]))
            elif r["role"] == "diag-scrambled": diag[fam].append((r, v))
    libfam = collections.defaultdict(list)
    for e in ex: libfam[FAMILY.get(e["view"], e["view"])].append(e)
    out = {"n1": [], "n2": [], "n3": [], "pos": []}; Rnull = {}
    for fam in set(libfam) | set(negv):
        keys, med, mad = zspace([nullfp.get(fam, []), [e["fp"] for e in libfam.get(fam, [])], [f for *_, f in negv.get(fam, [])]])
        Nz = Z(nullfp[fam], keys, med, mad) if nullfp.get(fam) else None
        if Nz is not None and len(Nz) >= 2:
            dnn = np.sqrt(((Nz[:, None] - Nz[None]) ** 2).sum(-1)); np.fill_diagonal(dnn, np.inf); Rn = float(np.percentile(dnn.min(1), 95))
        else:
            Rn = None
        Rnull[fam] = {"null_views": 0 if Nz is None else len(Nz), "R": Rn}
        in_null = lambda fp: Rn is not None and float(np.sqrt(((Nz - Z([fp], keys, med, mad)[0]) ** 2).sum(1)).min()) <= Rn
        namer = Namer(libfam[fam], keys, med, mad) if len(libfam.get(fam, [])) >= K + 1 else None
        for r, v, fp in negv.get(fam, []):
            if in_null(fp): continue
            (out["n2"] if r["role"] == "neg-regime" else out["n1"]).append((r["tag"], v, r["prog"][:120]))
            if namer:
                nm, _ = namer.name(fp)
                if nm: out["n3"].append((r["tag"], v, nm))
        for e in libfam.get(fam, []):
            flagged = e["screened"] and not in_null(e["fp"])
            nm, d = namer.name(e["fp"], exclude_variant=f"{e['class']}#{e['variant']}") if (namer and flagged) else (None, None)
            out["pos"].append({"class": e["class"], "variant": e["variant"], "seed": e["seed"], "flagged": flagged,
                               "named": nm, "correct": nm == e["class"]})
    pos = out["pos"]; flagged = [p for p in pos if p["flagged"]]; named = [p for p in flagged if p["named"]]
    wrong = [p for p in named if not p["correct"]]; right = [p for p in named if p["correct"]]
    n_neg = sum(1 for r in res if r["role"] == "neg-free")
    reg_ok = [r["tag"] for r in res if r["role"] == "neg-regime" and disorder_ok(r["tag"], r)]
    reg_drop = [r["tag"] for r in res if r["role"] == "neg-regime" and not disorder_ok(r["tag"], r)]
    crit = {f"N1 held-out interaction-free negatives flagged = {len(out['n1'])} / {n_neg}": len(out["n1"]) == 0,
            f"N2 verified regime negatives flagged = {len(out['n2'])} / {len(reg_ok)} (dropped as not disordered: {reg_drop})": len(out["n2"]) == 0,
            f"N3 negatives named = {len(out['n3'])}": len(out["n3"]) == 0,
            f"P1 verified examples flagged = {len(flagged)} / {len(pos)} (need >= 90%)": len(flagged) >= 0.9 * len(pos),
            f"P3 wrong names = {len(wrong)} / {len(named)} named": len(wrong) == 0}
    tab = collections.defaultdict(lambda: [0, 0, 0, 0])
    for p in pos: t = tab[p["class"]]; t[0] += 1; t[1] += p["flagged"]; t[2] += p["correct"]; t[3] += bool(p["named"]) and not p["correct"]
    conf = collections.Counter((p["class"], p["named"]) for p in wrong)
    L = ["# Validation gate v3 — reference-library naming — " + time.strftime("%Y-%m-%d %H:%M"), "",
         f"Null region per family: {json.dumps(Rnull)}", "",
         "| behaviour | verified examples | flagged | named correctly | named WRONG |", "|---|---|---|---|---|"]
    L += [f"| {c} | {n} | {f} | {rr} | {w} |" for c, (n, f, rr, w) in sorted(tab.items())]
    L += ["", "## Criteria", ""] + [f"- {'PASS' if v else 'FAIL'} — {k}" for k, v in crit.items()]
    L += [f"- (reported) P2 flagged examples named correctly = {len(right)} / {len(flagged)} "
          f"({100 * len(right) / max(len(flagged), 1):.0f}%; target >= 80%); unnamed (-> literature) = {len(flagged) - len(named)}"]
    if wrong: L += ["", "### Wrong names (true -> named)"] + [f"- {a} -> {b}: {n}" for (a, b), n in conf.most_common()]
    L += ["", f"Diagnostic (not scored): scrambled runs with a flagged view, by family: "
          + json.dumps({f: len({r['tag'] for r, _ in x}) for f, x in diag.items()})]
    if out["n1"]: L += ["", "### Held-out negatives flagged"] + [f"- {t} [{v}] `{p}`" for t, v, p in out["n1"][:30]]
    if out["n2"]: L += ["", "### Regime negatives flagged"] + [f"- {t} [{v}]" for t, v, _ in out["n2"]]
    if out["n3"]: L += ["", "### Negatives named"] + [f"- {t} [{v}] -> {n}" for t, v, n in out["n3"]]
    ub = 3.0 / max(n_neg, 1)
    L += ["", f"**GATE: {'PASS' if all(crit.values()) else 'FAIL'}**" + (f" — false-flag rate <= {100 * ub:.1f}% (95%)" if not out["n1"] else "")]
    open(os.path.join(a.out, "VALIDATION.md"), "w").write("\n".join(L) + "\n"); print("\n".join(L))
    json.dump(out, open(os.path.join(a.out, "scored.json"), "w"), indent=1, default=str)


if __name__ == "__main__":
    main()
