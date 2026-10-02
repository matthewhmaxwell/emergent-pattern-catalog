"""Pre-freeze check: known behaviours hand-encoded in G2 must be flagged and named like their G1 versions
(validated pipeline: screen + knock-out + textbook measures + library namer). Also 3 G2 negatives."""
import os, sys
os.environ.setdefault("OMP_NUM_THREADS", "1")
from multiprocessing import Pool
from census import grammar_g1 as G1, grammar_g2 as G2

R2 = G2.Rule2
def P(layers, k=None, nf=0, D=(), agent=(), rules=()): return G1.make(layers, k, nf, D, agent, rules)
CASES = {  # name: (program, expected behaviour or None for "must NOT be flagged")
 "voter": (P("C", {"C": 2}, rules=[R2("C", "any", "ALWAYS", (), "COPY_RAND", (), 0.3)]), "domain coarsening"),
 "majority": (P("C", {"C": 2}, rules=[R2("C", "any", "ALWAYS", (), "COPY_MAJ", (), 0.1)]), "domain coarsening"),
 "cyclic CA": (P("C", {"C": 3}, rules=[R2("C", "any", "NEXT_GE", (("th", 3),), "SET_NEXT", (), 1.0)]), "cyclic waves and spirals"),
 "lattice oscillators": (P("C", {"C": 1}, rules=[R2("C", "any", "ALWAYS", (), "PHASE_COUPLE", (("K", 1.0),), 1.0)]), "lattice phase locking"),
 "network oscillators": (P("N", {"N": 1}, rules=[R2("N", "any", "ALWAYS", (), "PHASE_COUPLE", (("K", 1.0),), 1.0)]), "network synchronization"),
 "network voter": (P("N", {"N": 2}, rules=[R2("N", "any", "ALWAYS", (), "COPY_RAND", (), 1.0)]), "network consensus"),
 "co-evolving voter": (P("N", {"N": 2}, rules=[R2("N", "any", "ALWAYS", (), "COPY_RAND", (), 1.0), R2("N", "any", "ALWAYS", (), "REWIRE_SAME", (), 1.0)]), "network fragmentation"),
 "flocking": (P("A", {"A": 1}, agent=(1.0, 0.3), rules=[R2("A", "any", "ALWAYS", (), "HEAD_MEAN", (("b", "any"), ("r", 3.0)), 1.0)]), "flocking"),
 "NEG random flips": (P("C", {"C": 2}, rules=[R2("C", 0, "ALWAYS", (), "SET", (("b", 1),), 0.1), R2("C", 1, "ALWAYS", (), "SET", (("b", 0),), 0.1)]), None),
 "NEG uncoupled turners": (P("A", {"A": 1}, agent=(1.0, 0.3), rules=[R2("A", "any", "ALWAYS", (), "TURN_L", (), 0.1)]), None),
 "NEG swamped copy": (P("C", {"C": 2}, rules=[R2("C", 0, "ALWAYS", (), "SET", (("b", 1),), 1.0), R2("C", "any", "ALWAYS", (), "COPY_RAND", (), 1.0)]), None),
}

def one(item):
    name, bits, lib = item
    from census.runner import evaluate
    from census.namer import CensusNamer
    p = G2.decode(bits); views, status, info = evaluate(p, bits); nm = CensusNamer(lib); out = {}
    for v, e in views.items():
        if e["flagged"] and v != "C.aval": r = nm.name_view(v, e["seeds"]); out[v] = (r["behaviour"], r["shown"], r["status"])
    return name, status, out, info["fp_errors"]

if __name__ == "__main__":
    lib = sys.argv[1]; ok = True
    with Pool(4) as pool:
        for name, status, out, fperr in pool.imap(one, [(n, G2.encode(p), lib) for n, (p, _) in CASES.items()]):
            exp = CASES[name][1]; shown = {b for _, sh, _ in out.values() for b in sh}
            good = (status == "NONE") if exp is None else (exp in shown)
            ok &= good and fperr == 0
            print(f"{'PASS' if good else 'FAIL'}  {name:22s} expect {str(exp):26s} status {status:8s} {out}", flush=True)
    print("ALL PASS" if ok else "SOME FAILED")
