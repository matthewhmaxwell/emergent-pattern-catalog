"""Spot-check the knock-out rule on known positives and the negatives that failed earlier (2 seed sets each)."""
from census import grammar_g1 as G
from census.knockout import run_with_knockout
R = G.Rule
one = lambda k: [R("N.FLIP", (("a", t), ("b", 1), ("p", 1.0))) for t in range(k) if t != 1]
CASES = [("POS flocking", G.make("A", {"A": 1}, agent=(1.0, 0.3), rules=[R("A.ALIGN", (("r", 3.0),))])),
         ("POS voter", G.make("C", {"C": 2}, rules=[R("C.COPY", (("p", 0.3),))])),
         ("POS network voter", G.make("N", {"N": 2}, rules=[R("N.COPY", (("p", 1.0),))])),
         ("POS net sync", G.make("N", {"N": 1}, rules=[R("N.KURA", (("K", 1.0),))])),
         ("POS Schelling", G.make("C", {"C": 3}, rules=[R("C.SCHELL", (("th", 0.5),))])),
         ("POS co-evolving", G.make("N", {"N": 2}, rules=[R("N.COPY", (("p", 1.0),)), R("N.REWIRE", (("p", 1.0),))])),
         ("NEG noise-drowned net", G.make("N", {"N": 2}, rules=[R("N.COPY", (("p", 0.01),)), R("N.FLIP", (("a", 0), ("b", 1), ("p", 0.5))), R("N.FLIP", (("a", 1), ("b", 0), ("p", 0.5)))])),
         ("NEG swamped rewire k=3", G.make("N", {"N": 3}, rules=[R("N.REWIRE", (("p", 1.0),))] + one(3))),
         ("NEG align max noise r=0.5 v=0.3", G.make("A", {"A": 1}, agent=(0.3, 3.14159), rules=[R("A.ALIGN", (("r", 0.5),))]))]
for name, p in CASES:
    res = []
    for sd in (12345, 999, 31337):
        out, real, drv, ko = run_with_knockout(p, sd)
        flagged = [v for v, s in real.items() if sum(x["emergent"] for x in s) >= 2 and drv[v][0]]
        res.append(",".join(flagged) or "-")
    print(f"{name:34s} flagged views per seed set: {res}", flush=True)
