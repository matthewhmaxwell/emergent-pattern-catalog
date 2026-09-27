"""Pilot T7: simulator + filter cost per rule template (single-rule programs, 3 seeds, 1000 steps)."""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import time
from census import grammar_g1 as G, sim, filter as FL
R = G.Rule
P = {
 "C.COPY": G.make("C", {"C": 2}, rules=[R("C.COPY", (("p", 1.0),))]),
 "C.MAJ": G.make("C", {"C": 3}, rules=[R("C.MAJ", (("p", 1.0),))]),
 "C.CONTAG": G.make("C", {"C": 2}, rules=[R("C.CONTAG", (("a", 0), ("b", 1), ("th", 3), ("p", 1.0)))]),
 "C.CYCLE": G.make("C", {"C": 3}, rules=[R("C.CYCLE", (("th", 3), ("p", 1.0)))]),
 "C.SWAP": G.make("C", {"C": 2}, rules=[R("C.SWAP", (("p", 1.0),))]),
 "C.SCHELL": G.make("C", {"C": 3}, rules=[R("C.SCHELL", (("th", 0.5),))]),
 "C.SAND": G.make("C", {"C": 4}, rules=[R("C.SAND", (("p", 0.001),))]),
 "C.GAME": G.make("C", {"C": 2}, rules=[R("C.GAME", (("g", (1.9, 0.0)),))]),
 "C.KURA": G.make("C", {"C": 1}, rules=[R("C.KURA", (("K", 1.0),))]),
 "A.ALIGN": G.make("A", {"A": 1}, agent=(1.0, 0.3), rules=[R("A.ALIGN", (("r", 3.0),))]),
 "A.ATTRACT r8": G.make("A", {"A": 1}, agent=(1.0, 0.3), rules=[R("A.ATTRACT", (("b", "any"), ("r", 8.0)))]),
 "A.CYCLE": G.make("A", {"A": 3}, agent=(1.0, 0.3), rules=[R("A.CYCLE", (("th", 1), ("r", 1.0)))]),
 "N.COPY": G.make("N", {"N": 2}, rules=[R("N.COPY", (("p", 1.0),))]),
 "N.MAJ": G.make("N", {"N": 3}, rules=[R("N.MAJ", (("p", 1.0),))]),
 "N.KURA": G.make("N", {"N": 1}, rules=[R("N.KURA", (("K", 1.0),))]),
 "N.REWIRE": G.make("N", {"N": 2}, rules=[R("N.REWIRE", (("p", 1.0),))]),
 "N.GAME": G.make("N", {"N": 2}, rules=[R("N.GAME", (("g", (1.9, 0.0)),))]),
 "F.AUTOCAT": G.make("F", nf=2, D=(0.2, 0.1), rules=[R("F.AUTOCAT", (("f", 0), ("g", 1), ("p", 1.0)))]),
}
for n, p in P.items():
    t = time.time(); out = sim.run(p, seed0=5); ts = time.time() - t
    t = time.time(); FL.evaluate(out, p, battery=False); tf = time.time() - t
    print(f"{n:14s} sim {ts:5.2f}s  filter {tf:5.2f}s", flush=True)
