"""Pilot: time screen components + fingerprint per view type."""
import time
from census import grammar_g1 as G, sim, filter as FL
from census.fingerprint import fingerprint
from epc.phase2a.emergence import generic_emergence
from epc.phase2a.novelty_tripwire import model_free_complexity
R = G.Rule
progs = {"grid voter": G.make("C", {"C": 2}, rules=[R("C.COPY", (("p", 0.3),))]),
         "lattice phase": G.make("C", {"C": 1}, rules=[R("C.KURA", (("K", 1.0),))]),
         "network voter": G.make("N", {"N": 2}, rules=[R("N.COPY", (("p", 1.0),))]),
         "network phase": G.make("N", {"N": 1}, rules=[R("N.KURA", (("K", 1.0),))]),
         "agents vicsek": G.make("A", {"A": 1}, agent=(1.0, 0.3), rules=[R("A.ALIGN", (("r", 3.0),))]),
         "field": G.make("F", nf=1, D=(0.2,), rules=[R("F.FEED", (("f", 0), ("F", 0.1)))])}
for name, p in progs.items():
    out = sim.run(p, seed0=3)
    for vn, h, md in FL.views(out, p, 0):
        th = FL._thin(h, FL.SCREEN_FRAMES)
        t = time.time(); em = generic_emergence(th, seed=0); t1 = time.time() - t
        t = time.time(); mf = model_free_complexity(th); t2 = time.time() - t
        t = time.time(); fp = fingerprint(vn, th); t3 = time.time() - t
        print(f"{name:14s} view {vn:8s} frames {len(th):3d} | generic {t1:6.2f}s  model-free {t2:6.2f}s  fingerprint {t3:5.2f}s | em {em.get('score'):.2f} {em.get('kind')} complex={mf.get('is_complex')}", flush=True)
