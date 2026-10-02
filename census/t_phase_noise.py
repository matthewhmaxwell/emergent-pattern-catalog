"""Check the PHASE-view screen: uncoupled / far-too-weakly-coupled oscillators must not fire; locked ones must.

  python -m census.t_phase_noise [n_seedsets]

final3 (2026-10-02) flagged `C kC=2 | C.KURA(K=0.002)` (C.phase) although an identical-physics case passed in
final2: this measures how often each screen channel fires per seed-run on phase fields with no order.
"""
import collections, sys
from multiprocessing import Pool
from census import grammar_g1 as G, sim, filter as FL
R = G.Rule
CASES = {
 "NULL lattice K=0": G.make("C", {"C": 1}, rules=[R("C.KURA", (("K", 0.0),))]),
 "NULL lattice K=0.002 k=1": G.make("C", {"C": 1}, rules=[R("C.KURA", (("K", 0.002),))]),
 "NULL lattice K=0.002 k=2": G.make("C", {"C": 2}, rules=[R("C.KURA", (("K", 0.002),))]),
 "NULL network K=0": G.make("N", {"N": 1}, rules=[R("N.KURA", (("K", 0.0),))]),
 "NULL network K=0.002": G.make("N", {"N": 1}, rules=[R("N.KURA", (("K", 0.002),))]),
 "POS lattice K=1.0": G.make("C", {"C": 1}, rules=[R("C.KURA", (("K", 1.0),))]),
 "POS lattice K=0.1": G.make("C", {"C": 1}, rules=[R("C.KURA", (("K", 0.1),))]),
 "POS network K=1.0": G.make("N", {"N": 1}, rules=[R("N.KURA", (("K", 1.0),))]),
}


def one(job):
    name, sd = job; p = CASES[name]; out = sim.run(p, seed0=7000 + 31 * sd); rows = []
    for s in range(3):
        for vn, hist, md in FL.views(out, p, s):
            if not vn.endswith("phase"): continue
            sc = FL.screen(FL._screen_hist(vn, hist, "adj" in out), complexity=False)   # as in the pipeline
            rows.append((name, sc["emergent"], sc["em_score"] >= 0.5, sc["is_complex"], sc["consensus_gain"] >= 0.2,
                         sc["em_kind"], sc["em_score"]))
    return rows


if __name__ == "__main__":
    n = int(sys.argv[1]) if len(sys.argv) > 1 else 10
    with Pool(5) as pool: res = [r for rows in pool.imap_unordered(one, [(c, sd) for c in CASES for sd in range(n)]) for r in rows]
    for c in CASES:
        rs = [r for r in res if r[0] == c]; kinds = collections.Counter(r[5] for r in rs if r[2])
        print(f"{c:28s} emergent {sum(r[1] for r in rs):3d} / {len(rs)} | generic>=0.5: {sum(r[2] for r in rs):3d}  complex: {sum(r[3] for r in rs):3d}"
              f"  gain: {sum(r[4] for r in rs):3d} | generic kinds when fired: {dict(kinds)} | scores {sorted(round(r[6], 2) for r in rs)[-4:]}", flush=True)
