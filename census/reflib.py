"""Reference library (owner-approved 2026-09-28): a better known-pattern detector, built from our own catalog.

Each known phenomenon is simulated many times IN THE CENSUS WORLD (same substrate, same run length, same screen and
fingerprint as the census). Every example must pass its TEXTBOOK CHECK (the classic order parameter for that
phenomenon, computed from the raw state) before it enters the library — so no mislabelled example gets in. A census
view is then NAMED by its nearest library examples; its error rate is measured on held-out parameter settings
(leave-one-variant-out) and on negatives (census.validate3).

  python -m census.reflib build [--out census/reflib/v1] [--workers 4] [--seedsets 2]

Classes map to catalog entries where the catalog has the phenomenon ("P18" ...) or are marked KNOWN-UNCATALOGUED
(with the textbook reference) when the phenomenon is well known but not in the 37-pattern catalog.
"""
import os
for _v in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
    os.environ.setdefault(_v, "1")
import argparse, json, sys, time
from multiprocessing import Pool
import numpy as np
from census import grammar_g1 as G

R = G.Rule


def _interface(g): return float(((g != np.roll(g, 1, 0)).mean() + (g != np.roll(g, 1, 1)).mean()) / 2)


def _moran(g):
    x = g.astype(float); x = x - x.mean(); v = (x ** 2).mean()
    if v < 1e-12: return 1.0                                 # uniform = fully ordered (consensus)
    nb = sum(np.roll(np.roll(x, dy, 0), dx, 1) for dy, dx in ((0, 1), (1, 0), (0, -1), (-1, 0))) / 4
    return float((x * nb).mean() / v)


def _same_share(g):
    occ = g != 0; same = sum((np.roll(np.roll(g, dy, 0), dx, 1) == g) & np.roll(np.roll(occ, dy, 0), dx, 1)
                             for dy in (-1, 0, 1) for dx in (-1, 0, 1) if (dy, dx) != (0, 0))
    tot = sum(np.roll(np.roll(occ, dy, 0), dx, 1) for dy in (-1, 0, 1) for dx in (-1, 0, 1) if (dy, dx) != (0, 0))
    m = occ & (tot > 0); return float((same[m] / tot[m]).mean()) if m.any() else 0.0


# ------------------------------------------------------------------ textbook checks (per seed s, raw output `o`)
def chk_coarsening(o, s):
    C = o["C"][s]; e, l = _interface(C[5]), np.mean([_interface(g) for g in C[-25:]])
    return l < 0.7 * e and _moran(C[-1]) >= 0.3, {"interface_early": e, "interface_late": l, "moran_late": _moran(C[-1])}


def chk_schelling(o, s):
    C = o["C"][s]; a, b = _same_share(C[0]), _same_share(C[-1])
    return b - a >= 0.10, {"same_share_start": a, "same_share_end": b}


def chk_cyclic(o, s):
    C = o["C"][s]; act = float(np.mean([(C[-i] != C[-i - 1]).mean() for i in range(1, 11)])); m = _moran(C[-1])
    k = int(C.max()) + 1; ch = C[-20:][1:] != C[-20:][:-1]
    cyc = float(((C[-20:][1:] == (C[-20:][:-1] + 1) % max(k, 2)) & ch).sum() / max(ch.sum(), 1)) - 1.0 / max(k - 1, 1)
    return act >= 0.05 and m >= 0.1 and k >= 3 and cyc >= 0.2, {"late_activity": act, "moran_late": m, "cyclic_excess": cyc}


def chk_game(o, s):
    C = o["C"][s]; c = float((C[-1] == 0).mean()); act = float((C[-1] != C[-6]).mean())
    return 0.05 < c < 0.95 and act > 0.01, {"coop_fraction": c, "activity": act}


def chk_net_sync(o, s):
    r = float(np.abs(np.exp(1j * o["Nph"][s, -20:]).mean(-1)).mean()); return r >= 0.6, {"r_late": r}


def chk_lat_sync(o, s):
    z = np.exp(1j * o["Cph"][s, -1]); loc = sum(np.roll(np.roll(z, dy, 0), dx, 1) for dy, dx in ((0, 1), (1, 0), (0, -1), (-1, 0))) / 4
    lr = float(np.abs(loc).mean()); return lr >= 0.8, {"local_r": lr}


def chk_flock(o, s):
    pol = float(np.abs(np.exp(1j * o["head"][s, -50:]).mean(-1)).mean()); return pol >= 0.5, {"polarization": pol}


def chk_coevo(o, s):
    import networkx as nx
    A = o["adj"][s, -1]; g = nx.from_numpy_array(A.astype(int)); giant = max(len(c) for c in nx.connected_components(g)) / A.shape[0]
    return giant <= 0.9, {"giant_fraction": giant}


def chk_net_consensus(o, s):
    nt = o["nt"][s]; m = float(np.bincount(nt[-1], minlength=4).max() / nt.shape[1]); m0 = float(np.bincount(nt[0], minlength=4).max() / nt.shape[1])
    return m >= 0.8 and m - m0 >= 0.2, {"majority_start": m0, "majority_end": m}


def classes():
    """{class: (catalog mapping, view, [programs], check)} — variants inside the regime where each is known to occur."""
    c = {}
    c["voter coarsening"] = ("P18", "C", [G.make("C", {"C": k}, rules=[R("C.COPY", (("p", p),))]) for k in (2, 3, 4) for p in (1.0, 0.3, 0.1, 0.03)], chk_coarsening)
    c["majority-rule coarsening"] = ("KNOWN-UNCATALOGUED (zero-temperature Glauber / majority-vote coarsening; Bray 1994 review)", "C",
                                     [G.make("C", {"C": k}, rules=[R("C.MAJ", (("p", p),))]) for k in (2, 3, 4) for p in (1.0, 0.3, 0.1, 0.03)], chk_coarsening)
    c["Schelling segregation"] = ("P1", "C", [G.make("C", {"C": k}, rules=[R("C.SCHELL", (("th", t),))]) for k in (3, 4) for t in (0.3, 0.4, 0.5, 0.6, 0.7)], chk_schelling)
    c["cyclic-CA waves/spirals"] = ("KNOWN-UNCATALOGUED (cyclic cellular automata; Fisch, Gravner & Griffeath 1991)", "C",
                                    [G.make("C", {"C": k}, rules=[R("C.CYCLE", (("th", t), ("p", p)))]) for k in (3, 4) for t in (1, 2, 3) for p in (1.0, 0.3)], chk_cyclic)
    c["spatial PD chaos (Nowak-May)"] = ("P27", "C", [G.make("C", {"C": 2}, rules=[R("C.GAME", (("g", (1.9, 0.0)),))])], chk_game)
    c["network synchronization"] = ("P9", "N.phase", [G.make("N", {"N": k}, rules=[R("N.KURA", (("K", K),))]) for k in (1, 2) for K in (1.0, 0.5, 0.3)], chk_net_sync)
    c["lattice phase locking"] = ("P9", "C.phase", [G.make("C", {"C": k}, rules=[R("C.KURA", (("K", K),))]) for k in (1, 2) for K in (1.0, 0.5, 0.3, 0.1)], chk_lat_sync)
    c["co-evolving network fragmentation"] = ("P34", "N", [G.make("N", {"N": k}, rules=[R("N.COPY", (("p", p),)), R("N.REWIRE", (("p", q),))])
                                                          for k in (2, 3) for p in (1.0, 0.3) for q in (1.0, 0.5)], chk_coevo)
    c["network voter consensus"] = ("P18", "N", [G.make("N", {"N": k}, rules=[R("N.COPY", (("p", p),))]) for k in (2, 3) for p in (1.0, 0.3, 0.1)], chk_net_consensus)
    c["flocking"] = ("P5", "A", [G.make("A", {"A": 1}, agent=(v, e), rules=[R("A.ALIGN", (("r", r),))])
                                 for v in (1.0, 0.3) for e in (0.1, 0.3, 1.0) for r in (1.0, 2.0, 3.0)], chk_flock)
    return c


def _build_one(job):
    cls, vi, bits, view, seedset = job
    from census.runner import seed_of
    from census.knockout import run_with_knockout
    p = G.decode(bits); out, real, driven = run_with_knockout(p, seed_of(bits) + 7919 * seedset)
    allc = classes(); check = allc[cls][3]; rows = []
    same_view = {c: v[3] for c, v in allc.items() if v[1] == view}
    for s, rs in enumerate(real.get(view, [])):
        ok, measures = check(out, s)
        checks = {}
        for c, fn in same_view.items():
            try: checks[c] = bool(fn(out, s)[0])
            except Exception: checks[c] = False
        rows.append({"class": cls, "variant": vi, "bits": bits, "prog": G.describe(p), "view": view, "seedset": seedset,
                     "seed": s, "verified": bool(ok), "measures": {k: round(float(v), 4) for k, v in measures.items()},
                     "screened": rs["emergent"], "driven": bool(driven[view][s]), "em_score": rs["em_score"], "fp": rs["fp"],
                     "checks": checks})
    return rows


def main():
    ap = argparse.ArgumentParser(); ap.add_argument("cmd"); ap.add_argument("--out", default="census/reflib/v1")
    ap.add_argument("--workers", type=int, default=4); ap.add_argument("--seedsets", type=int, default=2)
    ap.add_argument("--seedset-start", type=int, default=0, help="fresh data for confirmation runs: start at an unused seed set")
    a = ap.parse_args(); os.makedirs(a.out, exist_ok=True)
    jobs = [(cls, vi, G.encode(p), view, ss) for cls, (_, view, progs, _) in classes().items()
            for vi, p in enumerate(progs) for ss in range(a.seedset_start, a.seedset_start + a.seedsets)]
    print(f"{len(jobs)} simulations ({len(classes())} classes)", flush=True)
    rows = []
    with Pool(a.workers) as pool:
        for rr in pool.imap_unordered(_build_one, jobs):
            rows += rr
            if len(rows) % 60 == 0: print(f"{len(rows)} examples", flush=True)
    json.dump({"classes": {k: {"catalog": v[0], "view": v[1]} for k, v in classes().items()}, "examples": rows},
              open(os.path.join(a.out, "library.json"), "w"), indent=1)
    import collections
    t = collections.defaultdict(lambda: [0, 0, 0, 0])
    for r in rows: x = t[r["class"]]; x[0] += 1; x[1] += r["verified"]; x[2] += r["verified"] and r["screened"]; x[3] += r["verified"] and r["screened"] and r["driven"]
    print("class: examples / textbook-verified / + flagged by screen / + survives interaction knock-out")
    for k, (n, v, vs, vd) in t.items(): print(f"  {k:36s} {n:4d} {v:4d} {vs:4d} {vd:4d}")


if __name__ == "__main__":
    main()
