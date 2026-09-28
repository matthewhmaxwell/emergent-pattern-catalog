"""VALIDATION GATE — round 9 (owner-approved rule change 2026-09-28: MULTI-LABEL naming; a name is wrong only if
the run does not actually show that behaviour). Criteria fixed BEFORE the run; nothing scales until PASS.

NAMES PER RUN: PRIMARY = nearest-example namer + its textbook check (below); ALSO-SHOWS = every other behaviour whose
textbook measure passes (measure defined for the program: the game measure needs a game rule, the segregation
measure needs vacancies). Status: NAMED (has a primary), KNOWN-BEHAVIOUR-ATYPICAL-LOOK (only also-shows -> review
queue, so a new twist on a known behaviour is not hidden), UNNAMED (-> literature check).
Because every name is now verified by a textbook measure, naming correctness rests on those measures. They are
tested where it matters: a named behaviour must DISAPPEAR when the run's interactions are knocked out (C1). (Measures
passing on unflagged negatives — e.g. one-way conversion looks like "consensus" — are reported only: such runs are
never flagged, so never named.)

  python -m census.validate4 --lib census/reflib/v4/library.json [--workers 5] [--out census/validation/v6]

FLAG (per view, per seed) = census screen says emergent AND the interaction knock-out says interaction-driven
(census.knockout). A view is flagged if >= 2 of 3 seeds are.

NAMING = OPEN-SET namer: per view family, robust-z fingerprints (scale floored at 0.5 x std); each nameable class
(>= 3 parameter variants) gets a shrunk-covariance Gaussian (shrinkage 0.5 toward the diagonal, variance floor 0.1); its typicality threshold tau_c = 95th percentile of
its own examples' leave-one-variant-out Mahalanobis distances (nested: an example's own variant never sets its
threshold). A view is named c only if it is typical for EXACTLY ONE class; typical for none or several -> UNNAMED
(-> literature check).

TESTS
  P1  recall: >= 90% of textbook-verified library examples flagged.
  P3  names the run does not show (primary or also-shows failing its textbook measure): 0 (implementation check).
  C1  every name given (primary or also-shows) disappears under the interaction knock-out: the knock-out run fails
      that behaviour's textbook measure. 0 exceptions.
  (reported) P4 novelty-risk indicator, leave-one-CLASS-out: withheld-class examples that still get a PRIMARY name
      (listed, with whether the named behaviour is genuinely co-present).
  N1  0 flagged among verified-trivial INTERACTING programs: a real interaction rule swamped by one-way conversion
      (state verified uniform at the end), plus regime negatives verified disordered globally AND locally.
  N2  0 flagged among 400 held-out interaction-free programs (guaranteed by the knock-out design — reported as a
      consistency check of the implementation, not as evidence).
  F0  0 fingerprint errors among all evaluated views (integrity: round 6 found the agent fingerprint had failed
      silently in every earlier round — a numpy 2.x API change — so errors now fail the gate outright).
  (reported) named-correctly rate among flagged verified examples (target >= 80%).
"""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import argparse, collections, json, random, time
from multiprocessing import Pool
import numpy as np
from census import grammar_g1 as G
from census.triage import FAMILY
from census.knockout import is_interaction_free

R = G.Rule
SHRINK, PCT, MIN_VARIANTS, VAR_FLOOR = 0.5, 95, 3, 0.1


# ------------------------------------------------------------------ negatives
def trivial_interacting():
    """a real interaction rule swamped by a one-way conversion (FLIP 0->1 every step), k = 2 cells / nodes."""
    out = {}
    base_c = G.make("C", {"C": 2}); base_n = G.make("N", {"N": 2})
    for base, flip in ((base_c, "C.FLIP"), (base_n, "N.FLIP")):
        for t in G.valid_templates(base):
            if t in ("C.FLIP", "N.FLIP") or t.startswith(("C.KURA", "N.KURA")): continue
            prm = tuple((n, G._param_options(k, base)[0][0]) for n, k in G.T[t]["params"])
            if not G._ok_types(t, prm): continue
            out[f"swamped {t}"] = G.make(base.layers, base.k, rules=[R(flip, (("a", 0), ("b", 1), ("p", 1.0))), R(t, prm)])
    return out


def regime_negatives():
    out = {}
    for v in (1.0, 0.3):
        for r in (0.5, 1.0):
            out[f"walkers repel r={r} v={v}"] = G.make("A", {"A": 1}, agent=(v, 3.14159), rules=[R("A.REPEL", (("b", "any"), ("r", r)))])
            out[f"align at max noise r={r} v={v}"] = G.make("A", {"A": 1}, agent=(v, 3.14159), rules=[R("A.ALIGN", (("r", r),))])
    for K in (0.0003, 0.001, 0.002):
        out[f"network oscillators K={K}"] = G.make("N", {"N": 1}, rules=[R("N.KURA", (("K", K),))])
        out[f"lattice oscillators K={K}"] = G.make("C", {"C": 1}, rules=[R("C.KURA", (("K", K),))])
    out["neighbour swaps"] = G.make("C", {"C": 2}, rules=[R("C.SWAP", (("p", 1.0),))])
    out["neighbour swaps, 3 types"] = G.make("C", {"C": 3}, rules=[R("C.SWAP", (("p", 0.3),))])
    return out


def verified_negative(tag, out, real):
    """textbook verification on the raw run (all seeds): trivial / disordered, else the case is dropped."""
    fps = [x["fp"] for seeds in real.values() for x in seeds]
    avg = lambda k, d: float(np.mean([f.get(k, d) for f in fps if k in f])) if any(k in f for f in fps) else d
    if tag.startswith("swamped"):
        st = out["C"] if "C" in out else out["nt"]
        return bool(all((np.unique(st[s, -1]).size == 1) for s in range(st.shape[0])))
    if tag.startswith(("walkers", "align")):
        H = out["head"][:, -50:]; pol = float(np.abs(np.exp(1j * H).mean(-1)).mean())
        return pol < 0.2 and avg("a_local_align", 0.0) < 0.05
    if tag.startswith("network oscillators"):
        r = float(np.abs(np.exp(1j * out["Nph"][:, -20:]).mean(-1)).mean()); return r < 0.3
    if tag.startswith("lattice oscillators"):
        z = np.exp(1j * out["Cph"][:, -1]); loc = sum(np.roll(np.roll(z, dy, -2), dx, -1) for dy, dx in ((0, 1), (1, 0), (0, -1), (-1, 0))) / 4
        return float(np.abs(loc).mean()) < 0.6
    if tag.startswith("neighbour swaps"):
        C = out["C"][:, -1].astype(float); C = C - C.mean((1, 2), keepdims=True); v = (C ** 2).mean((1, 2))
        nb = sum(np.roll(np.roll(C, dy, -2), dx, -1) for dy, dx in ((0, 1), (1, 0), (0, -1), (-1, 0))) / 4
        return bool(np.all((C * nb).mean((1, 2)) / np.maximum(v, 1e-9) < 0.1))
    return True


SEED_OFFSET = 0


def _neg(job):
    tag, role, bits = job
    from census.runner import seed_of
    from census.knockout import run_with_knockout
    p = G.decode(bits); out, real, driven, _ = run_with_knockout(p, seed_of(bits) + SEED_OFFSET)
    flagged = [v for v, seeds in real.items() if sum(1 for s, x in enumerate(seeds) if x["emergent"] and driven[v][s]) >= 2]
    ok = verified_negative(tag, out, real) if role != "free" else True
    fp = {v: {k: float(np.mean([x["fp"].get(k, 0.0) for x in seeds])) for k in set().union(*[x["fp"] for x in seeds])} for v, seeds in real.items()}
    fperr = sum(1 for seeds in real.values() for x in seeds if "fp_error" in x["fp"])
    from census.reflib import applicable_checks
    passed = []                                      # textbook measures passing on >= 2 of 3 seeds, any view
    for v in real:
        for c, fn in applicable_checks(p, v).items():
            n_ok = 0
            for s in range(len(real[v])):
                try: n_ok += bool(fn(out, s)[0])
                except Exception: pass
            if n_ok >= 2: passed.append(c)
    return {"tag": tag, "role": role, "prog": G.describe(p), "verified": ok, "flagged": flagged, "fp": fp,
            "fp_errors": fperr, "checks_passed": passed}


# ------------------------------------------------------------------ open-set namer
def _robust(X):
    """robust centre/scale with a floor: a feature that is near-constant in most rows (MAD ~ 0) is scaled by its
    standard deviation instead, so it cannot blow up into +/-8 z-units (round-7 degeneracy)."""
    med = np.median(X, 0); mad = np.median(np.abs(X - med), 0) * 1.4826
    scale = np.maximum(mad, 0.5 * X.std(0)); scale[scale < 1e-6] = 1.0
    return med, scale


def _fit(Zc):
    mu = Zc.mean(0); S = np.atleast_2d(np.cov(Zc, rowvar=False)) if len(Zc) > 1 else np.eye(Zc.shape[1])
    S = (1 - SHRINK) * S + SHRINK * np.diag(np.diag(S)) + VAR_FLOOR * np.eye(len(mu))   # variance floor (z-units)
    return mu, np.linalg.inv(S)


def _d2(z, mu, ic): d = z - mu; return float(d @ ic @ d)


class NearestNamer:
    K, VOTES, RPCT, MARGIN = 5, 5, 90, 1.5

    def __init__(self, rows, keys, med, mad):
        self.keys, self.med, self.mad = keys, med, mad
        self.Z = self._z([r["fp"] for r in rows]); self.cls = np.array([r["class"] for r in rows])
        self.var = np.array([f"{r['class']}#{r['variant']}" for r in rows])
        self.classes = sorted(c for c in set(self.cls) if len(set(self.var[self.cls == c])) >= MIN_VARIANTS)

    def _z(self, fps):
        X = np.array([[f.get(k, 0.0) for k in self.keys] for f in fps], float).reshape(len(fps), len(self.keys))
        return np.clip((X - self.med) / self.mad, -8, 8)

    def name(self, fp, checks, exclude_variant=None, exclude_class=None):
        z = self._z([fp])[0]; m = np.ones(len(self.cls), bool)
        if exclude_variant: m &= self.var != exclude_variant
        if exclude_class: m &= self.cls != exclude_class
        if m.sum() < self.K: return None
        Zr, cr, vr = self.Z[m], self.cls[m], self.var[m]
        d = np.sqrt(((Zr - z) ** 2).sum(1)); o = np.argsort(d)[:self.K]
        top, c = collections.Counter(cr[o]).most_common(1)[0]
        if top not in self.classes or top == exclude_class or c < self.VOTES: return None
        idx = np.flatnonzero(cr == top); nd = []
        for j in idx:
            oo = idx[vr[idx] != vr[j]]
            if len(oo): nd.append(np.sqrt(((Zr[oo] - Zr[j]) ** 2).sum(1)).min())
        R = np.percentile(nd, self.RPCT) if nd else 0.0
        dc = d[cr == top].min(); dother = d[cr != top].min() if (cr != top).any() else np.inf
        if dc > R or dother < self.MARGIN * dc: return None
        return top if checks.get(top, False) else None                       # name-then-verify


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("--lib", required=True); ap.add_argument("--workers", type=int, default=5)
    ap.add_argument("--out", default="census/validation/v6")
    ap.add_argument("--neg-cache", default=None, help="reuse negatives.json from a run with identical sim + fingerprint code")
    ap.add_argument("--seed-offset", type=int, default=0, help="fresh negative runs for confirmation")
    ap.add_argument("--free-start", type=int, default=0, help="fresh held-out interaction-free programs for confirmation")
    a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
    global SEED_OFFSET; SEED_OFFSET = a.seed_offset
    lib = json.load(open(a.lib)); ex = [r for r in lib["examples"] if r["verified"]]
    progs = [p for _, p in G.enumerate_programs(18)]
    free = [p for p in progs if is_interaction_free(p)]; random.Random(20260928).shuffle(free)
    jobs = [(t, "trivial", G.encode(p)) for t, p in trivial_interacting().items()]
    jobs += [(t, "regime", G.encode(p)) for t, p in regime_negatives().items()]
    jobs += [(f"free:{i}", "free", G.encode(p)) for i, p in enumerate(free[a.free_start:a.free_start + 400])]
    print(f"negatives: {sum(j[1] == 'trivial' for j in jobs)} swamped-interaction + {sum(j[1] == 'regime' for j in jobs)} regime "
          f"+ 400 interaction-free | positives: {len(ex)} verified library examples", flush=True)
    if a.neg_cache:
        neg = json.load(open(a.neg_cache))
    else:
        with Pool(a.workers) as pool: neg = list(pool.imap_unordered(_neg, jobs))
    json.dump(neg, open(os.path.join(a.out, "negatives.json"), "w"), indent=1, default=str)

    # ---- positives: flag, LOVO naming, LOCO open-set
    fam_rows = collections.defaultdict(list)
    for e in ex: fam_rows[FAMILY.get(e["view"], e["view"])].append(e)
    pos = []; loco = []
    for fam, rows in fam_rows.items():
        X = np.array([[r["fp"].get(k, 0.0) for k in sorted(set().union(*[r["fp"] for r in rows]))] for r in rows])
        keys = sorted(set().union(*[r["fp"] for r in rows])); med, mad = _robust(X)
        namer = NearestNamer(rows, keys, med, mad)
        from census.reflib import PRECONDITION
        for r in rows:
            prog = G.decode(r["bits"]); fl = bool(r["screened"] and r["driven"]); var = f"{r['class']}#{r['variant']}"
            ck = {c: ok for c, ok in r.get("checks", {}).items() if PRECONDITION.get(c, lambda q: True)(prog)}
            nm = namer.name(r["fp"], ck, exclude_variant=var) if fl else None
            also = sorted(c for c, ok in ck.items() if ok and c != nm) if fl else []
            bad = [x for x in ([nm] if nm else []) + also if not ck.get(x, False)]
            ko_survive = [x for x in ([nm] if nm else []) + also if r.get("ko_checks", {}).get(x, False)]
            pos.append({"class": r["class"], "flagged": fl, "named": nm, "also": also, "not_shown": bad,
                        "ko_survive": ko_survive, "nameable": r["class"] in namer.classes})
            if fl:
                nm2 = namer.name(r["fp"], ck, exclude_class=r["class"])
                loco.append({"class": r["class"], "named_as": nm2, "copresent": bool(nm2 and ck.get(nm2, False))})
    flagged = [p for p in pos if p["flagged"]]
    notshown = [p for p in flagged if p["not_shown"]]
    right = [p for p in flagged if p["named"] == p["class"]]
    other_primary = [p for p in flagged if p["named"] and p["named"] != p["class"]]
    unk = [x for x in loco if x["named_as"]]
    status = collections.Counter("NAMED" if p["named"] else ("KNOWN-BEHAVIOUR-ATYPICAL-LOOK" if p["also"] else "UNNAMED") for p in flagged)
    c1 = [(p["class"], p["ko_survive"]) for p in flagged if p["ko_survive"]]
    neg_measures = [(n["tag"], n["checks_passed"]) for n in neg if (n["role"] == "free" or n["verified"]) and n.get("checks_passed")]
    trv = [n for n in neg if n["role"] == "trivial" and n["verified"]]; reg = [n for n in neg if n["role"] == "regime" and n["verified"]]
    fre = [n for n in neg if n["role"] == "free"]
    n1 = [n for n in trv + reg if n["flagged"]]; n2 = [n for n in fre if n["flagged"]]
    dropped = [n["tag"] for n in neg if n["role"] in ("trivial", "regime") and not n["verified"]]
    fperr = sum(n.get("fp_errors", 0) for n in neg) + sum(1 for e in lib["examples"] if "fp_error" in e["fp"])
    crit = {f"F0 fingerprint errors = {fperr}": fperr == 0,
            f"P1 verified examples flagged = {len(flagged)} / {len(pos)} (need >= 90%)": len(flagged) >= 0.9 * len(pos),
            f"P3 names the run does not show = {len(notshown)} (implementation check)": len(notshown) == 0,
            f"C1 names whose behaviour survives the interaction knock-out = {len(c1)} / {len(flagged)}": len(c1) == 0,
            f"N1 verified trivial/disordered interacting negatives flagged = {len(n1)} / {len(trv) + len(reg)}": len(n1) == 0,
            f"N2 interaction-free flagged = {len(n2)} / {len(fre)} (consistency check)": len(n2) == 0}
    tab = collections.defaultdict(lambda: [0, 0, 0, 0])
    for p in pos: t = tab[p["class"]]; t[0] += 1; t[1] += p["flagged"]; t[2] += p["named"] == p["class"]; t[3] += bool(p["named"]) and p["named"] != p["class"]
    L = ["# Validation gate — round 9 (multi-label, verified names) — " + time.strftime("%Y-%m-%d %H:%M"), "",
         "| behaviour | verified | flagged | primary = own behaviour | primary = another (co-present) behaviour | withheld -> primary name (LOCO) |",
         "|---|---|---|---|---|---|"]
    lc = collections.Counter(x["class"] for x in unk)
    L += [f"| {c} | {n} | {f} | {rr} | {w} | {lc.get(c, 0)} |" for c, (n, f, rr, w) in sorted(tab.items())]
    L += ["", f"Status of flagged verified examples: {dict(status)}",
          "Also-shows pairs (behaviour -> also shows): " + json.dumps(dict(collections.Counter((p['class'], a) for p in flagged for a in p['also']).most_common(12)), default=str)]
    L += ["", "## Criteria", ""] + [f"- {'PASS' if v else 'FAIL'} — {k}" for k, v in crit.items()]
    L += [f"- (reported) primary name = own behaviour: {len(right)} / {len(flagged)} ({100 * len(right) / max(len(flagged), 1):.0f}%; target >= 80%)",
          f"- (reported) P4 novelty-risk: withheld-class examples still given a primary name = {len(unk)} / {len(loco)} "
          f"(named behaviour genuinely co-present in {sum(x['copresent'] for x in unk)})",
          f"- negatives dropped as NOT verified trivial/disordered (not scored): {dropped}"]
    if other_primary: L += ["", "### Primary name = another behaviour (verified co-present)"] + [f"- {a} -> {b}: {n}" for (a, b), n in collections.Counter((p['class'], p['named']) for p in other_primary).most_common()]
    if notshown: L += ["", "### Names not shown by the run (implementation bug!)"] + [f"- {p['class']}: {p['not_shown']}" for p in notshown[:20]]
    if c1: L += ["", "### Names surviving the knock-out"] + [f"- {a}: {b}" for (a, b) in collections.Counter((c, tuple(x)) for c, x in c1).most_common()]
    L += ["", f"(info) textbook measures passing on UNFLAGGED negatives (never named): {len(neg_measures)}"]
    if unk: L += ["", "### Withheld behaviour -> primary name (LOCO)"] + [f"- {a} -> {b}: {n}" for (a, b), n in collections.Counter((x['class'], x['named_as']) for x in unk).most_common()]
    if n1: L += ["", "### Negatives flagged"] + [f"- {n['tag']} {n['flagged']}" for n in n1]
    if n2: L += ["", "### Interaction-free flagged (implementation bug!)"] + [f"- {n['prog'][:120]}" for n in n2]
    L += ["", f"**GATE: {'PASS' if all(crit.values()) else 'FAIL'}**"]
    open(os.path.join(a.out, "VALIDATION.md"), "w").write("\n".join(L) + "\n"); print("\n".join(L))


if __name__ == "__main__":
    main()
