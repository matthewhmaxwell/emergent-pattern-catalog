"""Check the network node-state screen: pure noise must not fire; real network behaviours must."""
import collections, numpy as np
from census import grammar_g1 as G, sim, filter as FL
R = G.Rule
CASES = {
 "NOISE two-way flips 0.3": G.make("N", {"N": 2}, rules=[R("N.FLIP", (("a", 0), ("b", 1), ("p", 0.3))), R("N.FLIP", (("a", 1), ("b", 0), ("p", 0.3)))]),
 "NOISE two-way flips 0.03, 3 types": G.make("N", {"N": 3}, rules=[R("N.FLIP", (("a", 0), ("b", 1), ("p", 0.03))), R("N.FLIP", (("a", 1), ("b", 2), ("p", 0.03))), R("N.FLIP", (("a", 2), ("b", 0), ("p", 0.03)))]),
 "NOISE-DROWNED copy 0.03 + flips 0.3": G.make("N", {"N": 2}, rules=[R("N.COPY", (("p", 0.03),)), R("N.FLIP", (("a", 0), ("b", 1), ("p", 0.3))), R("N.FLIP", (("a", 1), ("b", 0), ("p", 0.3)))]),
 "POS network voter": G.make("N", {"N": 2}, rules=[R("N.COPY", (("p", 1.0),))]),
 "POS network voter 3 types slow": G.make("N", {"N": 3}, rules=[R("N.COPY", (("p", 0.1),))]),
 "POS network majority": G.make("N", {"N": 2}, rules=[R("N.MAJ", (("p", 0.1),))]),
 "POS network cyclic": G.make("N", {"N": 3}, rules=[R("N.CYCLE", (("th", 1), ("p", 1.0)))]),
}
for name, p in CASES.items():
    res = collections.Counter(); ex = None
    for sd in range(7):
        out = sim.run(p, seed0=2000 + sd)
        for s in range(3):
            vn, hist, md = FL.views(out, p, s)[0]; th = FL._screen_hist(vn, hist, "adj" in out)
            sc = FL.screen(th, adj=out["adj0"][s], network=True); res[sc["emergent"]] += 1; ex = sc
    print(f"{name:38s} emergent in {res[True]:2d} / 21 seed-runs | last: kind={ex['em_kind']} score={ex['em_score']} z={ex.get('agree_z')} gain={ex.get('consensus_gain')} osc={ex.get('oscillation')}", flush=True)
