import sys, time, collections
pass
from census import grammar_g1 as G1, grammar_g2 as G2
R2 = G2.Rule2
def P(layers, k=None, nf=0, D=(), agent=(), rules=()): return G1.make(layers, k, nf, D, agent, rules)
bench = {
 "voter":        P("C", {"C": 2}, rules=[R2("C", "any", "ALWAYS", (), "COPY_RAND", (), 1.0)]),
 "voter p=0.3":  P("C", {"C": 2}, rules=[R2("C", "any", "ALWAYS", (), "COPY_RAND", (), 0.3)]),
 "majority":     P("C", {"C": 2}, rules=[R2("C", "any", "ALWAYS", (), "COPY_MAJ", (), 1.0)]),
 "Schelling":    P("C", {"C": 3}, rules=[R2("C", "any", "SAME_LT", (("th", 0.5),), "MOVE_EMPTY", (), 1.0)]),
 "cyclic CA":    P("C", {"C": 3}, rules=[R2("C", "any", "NEXT_GE", (("th", 3),), "SET_NEXT", (), 1.0)]),
 "GH":           P("C", {"C": 3}, rules=[R2("C", "any", "NEXT_GE", (("th", 1),), "SET_NEXT", (), 1.0)]),
 "sandpile":     P("C", {"C": 4}, rules=[R2("C", "any", "ALWAYS", (), "ADD_GRAIN", (), 0.001)]),
 "Nowak-May":    P("C", {"C": 2}, rules=[R2("C", "any", "ALWAYS", (), "IMITATE_BEST", (("g", (1.9, 0.0)),), 1.0)]),
 "lat Kuramoto": P("C", {"C": 1}, rules=[R2("C", "any", "ALWAYS", (), "PHASE_COUPLE", (("K", 0.1),), 1.0)]),
 "Vicsek":       P("A", {"A": 1}, agent=(0.3, 0.3), rules=[R2("A", "any", "ALWAYS", (), "HEAD_MEAN", (("b", "any"), ("r", 1.0)), 1.0)]),
 "MIPS":         P("A", {"A": 1}, agent=(1.0, 1.0), rules=[R2("A", "any", "CNT_GE", (("b", "any"), ("th", 3), ("r", 1.0)), "SLOW", (), 1.0)]),
 "net Kuramoto": P("N", {"N": 1}, rules=[R2("N", "any", "ALWAYS", (), "PHASE_COUPLE", (("K", 0.1),), 1.0)]),
 "coevo voter":  P("N", {"N": 2}, rules=[R2("N", "any", "ALWAYS", (), "COPY_RAND", (), 1.0), R2("N", "any", "ALWAYS", (), "REWIRE_SAME", (), 1.0)]),
 "Langton ants": P("CA", {"C": 2, "A": 1}, agent=(1.0, 0.0), rules=[R2("A", "any", "CELL_EVEN", (), "TURN_L", (), 1.0), R2("A", "any", "CELL_ODD", (), "TURN_R", (), 1.0), R2("A", "any", "ALWAYS", (), "CELL_NEXT", (), 1.0)]),
 "Gray-Scott":   P("F", nf=2, D=(0.2, 0.1), rules=[R2("F", "any", "ALWAYS", (), "AUTOCAT", (("f", 0), ("g", 1), ("p", 1.0)), 1.0), R2("F", "any", "ALWAYS", (), "RELAX1", (("f", 0), ("F", 0.03)), 1.0), R2("F", "any", "ALWAYS", (), "SCALE", (("f", 1), ("d", 0.1)), 1.0)]),
}
for n, p in bench.items():
    b = G2.encode(p); cb = G2.canonical_bits(p); print(f"{n:14s} G2 {len(b):3d} bits (canonical {len(cb)})")
t = time.time(); cnt = collections.Counter()
for bits, p in G2.enumerate_programs(int(sys.argv[1]) if len(sys.argv) > 1 else 13): cnt[len(bits)] += 1
print("G2 counts:", sorted(cnt.items()), "total", sum(cnt.values()), "%.1fs" % (time.time() - t))
