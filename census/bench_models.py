"""Hand encodings of textbook models in G1 (recovery benchmark, DESIGN §7). Pilot draft."""
from . import grammar_g1 as G
R = G.Rule
bench = {
 "voter":            G.make("C", {"C": 2}, rules=[R("C.COPY", (("p", 1.0),))]),
 "majority coarsen": G.make("C", {"C": 2}, rules=[R("C.MAJ", (("p", 1.0),))]),
 "Schelling":        G.make("C", {"C": 3}, rules=[R("C.SCHELL", (("th", 0.5),))]),
 "cyclic CA spirals":G.make("C", {"C": 3}, rules=[R("C.CYCLE", (("th", 3), ("p", 1.0)))]),
 "BTW sandpile":     G.make("C", {"C": 4}, rules=[R("C.SAND", (("p", 0.001),))]),
 "Nowak-May PD":     G.make("C", {"C": 2}, rules=[R("C.GAME", (("g", (1.9, 0.0)),))]),
 "lattice Kuramoto": G.make("C", {"C": 1}, rules=[R("C.KURA", (("K", 0.1),))]),
 "Vicsek":           G.make("A", {"A": 1}, agent=(0.3, 0.3), rules=[R("A.ALIGN", (("r", 1.0),))]),
 "MIPS (quorum)":    G.make("A", {"A": 1}, agent=(1.0, 1.0), rules=[R("A.QUORUM", (("th", 3), ("r", 1.0)))]),
 "net Kuramoto":     G.make("N", {"N": 1}, rules=[R("N.KURA", (("K", 0.1),))]),
 "co-evolving voter":G.make("N", {"N": 2}, rules=[R("N.COPY", (("p", 1.0),)), R("N.REWIRE", (("p", 1.0),))]),
 "Greenberg-Hastings":G.make("C", {"C": 3}, rules=[R("C.CYCLE", (("th", 1), ("p", 1.0)))]),
 "Langton ant(s)":   G.make("CA", {"C": 2, "A": 1}, agent=(1.0, 0.0), rules=[R("A.TURNCELL"), R("A.WRITECELL")]),
 "Gray-Scott":       G.make("F", nf=2, D=(0.2, 0.1), rules=[R("F.AUTOCAT", (("f", 0), ("g", 1), ("p", 1.0))),
                                  R("F.FEED", (("f", 0), ("F", 0.03))), R("F.DECAY", (("f", 1), ("d", 0.1)))]),
 "chemotaxis (KS)":  G.make("AF", {"A": 1}, nf=1, D=(0.2,), agent=(1.0, 0.3), rules=[R("A.DEPOSIT", (("f", 0), ("p", 1.0))), R("A.CLIMB", (("f", 0),)), R("F.DECAY", (("f", 0), ("d", 0.1)))]),
 "RPS agents":       G.make("A", {"A": 3}, agent=(0.3, 1.0), rules=[R("A.CYCLE", (("th", 1), ("r", 1.0)))]),
}
BENCH = bench

if __name__ == "__main__":
  for name, p in bench.items():
    b = G.encode(p); cb = G.canonical_bits(p); q = G.decode(b)
    assert q.rules == p.rules
    print(f"{name:20s} {len(b):3d} bits (canonical {len(cb)})  {G.describe(p)}")
