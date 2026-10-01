"""DEV-SET exploration (not a gate): compare namer designs on the saved library (LOVO wrong names, LOCO unknown->named,
named-correct rate). The chosen design is then fixed and tested on FRESH data (confirmation run)."""
import json, sys, collections, numpy as np
from census.namer import FAMILY
from census.validate4 import _robust
lib = [r for r in json.load(open(sys.argv[1]))["examples"] if r["verified"] and r["screened"] and r["driven"]]
fam = collections.defaultdict(list)
for r in lib: fam[FAMILY.get(r["view"], r["view"])].append(r)

def evaluate(namefn, label):
    tot = wrong = right = unk = n = 0; conf = collections.Counter()
    for f, rows in fam.items():
        keys = sorted(set().union(*[r["fp"] for r in rows]))
        X = np.array([[r["fp"].get(k, 0.0) for k in keys] for r in rows]); med, sc = _robust(X)
        Z = np.clip((X - med) / sc, -8, 8); cls = np.array([r["class"] for r in rows]); var = np.array([f"{r['class']}#{r['variant']}" for r in rows])
        nameable = {c for c in set(cls) if len(set(var[cls == c])) >= 3}
        for i in range(len(rows)):
            n += 1
            m = var != var[i]; nm = namefn(Z[i], Z[m], cls[m], var[m], nameable)
            if nm == cls[i]: right += 1
            elif nm: wrong += 1; conf[(cls[i], nm)] += 1
            m2 = cls != cls[i]; nm2 = namefn(Z[i], Z[m2], cls[m2], var[m2], nameable - {cls[i]})
            if nm2: unk += 1; conf[("LOCO " + cls[i], nm2)] += 1
    print(f"{label:48s} correct {right}/{n} ({100*right/n:.0f}%) | wrong {wrong} | unknown->named {unk} | top: {conf.most_common(4)}")

def knn(K, votes, pct, margin):
    def f(z, Zr, cr, vr, nameable):
        if len(Zr) == 0: return None
        d = np.sqrt(((Zr - z) ** 2).sum(1)); o = np.argsort(d)[:K]; top, c = collections.Counter(cr[o]).most_common(1)[0]
        if top not in nameable or c < votes: return None
        # R_top: within-class nearest OTHER-variant distances (from the reference rows only)
        idx = np.flatnonzero(cr == top); nd = []
        for j in idx:
            oo = idx[vr[idx] != vr[j]]
            if len(oo): nd.append(np.sqrt(((Zr[oo] - Zr[j]) ** 2).sum(1)).min())
        R = np.percentile(nd, pct) if nd else 0
        dc = d[cr == top].min(); dother = d[cr != top].min() if (cr != top).any() else np.inf
        return top if dc <= R and dother >= margin * dc else None
    return f

for K, votes, pct, margin in [(5, 4, 95, 1.0), (5, 5, 95, 1.0), (5, 5, 90, 1.0), (5, 5, 90, 1.5), (5, 5, 80, 1.5), (7, 7, 90, 1.5), (5, 5, 90, 2.0)]:
    evaluate(knn(K, votes, pct, margin), f"kNN k={K} votes={votes} R@{pct}pct margin={margin}")
