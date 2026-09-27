"""Pilot diagnostic: run each expected detector directly on a census view (fast settings) and print its verdict."""
import sys, traceback
from census.bench_models import BENCH
from census import sim, filter as FL

CASES = [("cyclic CA spirals", "C", "P12"), ("Greenberg-Hastings", "C", "P13"), ("Vicsek", "A", "P5"),
         ("BTW sandpile", "C.aval", "P14"), ("Nowak-May PD", "C", "P27"), ("lattice Kuramoto", "C.phase", "P9"),
         ("co-evolving voter", "N", "P34"), ("MIPS (quorum)", "A", "P2"), ("RPS agents", "A", "P12"),
         ("net Kuramoto", "N.phase", "P9")]
F = FL._fast_fns()
for name, view, pid in CASES:
    if sys.argv[1:] and name not in sys.argv[1:]: continue
    p = BENCH[name]; out = sim.run(p, seed0=1)
    V = {n: (h, md) for n, h, md in FL.views(out, p, 0)}; h, md = V[view]; h = FL._thin(h, 125)
    try:
        r = F[pid](h, md)
        pm = getattr(r, "primary_metric", None)
        if isinstance(pm, dict): pm = {k: (round(v, 3) if isinstance(v, float) else v) for k, v in pm.items()}
        print(f"{name} [{pid}] detected={getattr(r, 'detected', None)} tier={getattr(r, 'tier', None)} "
              f"primary={str(pm)[:260]} warnings={str(getattr(r, 'warnings', ''))[:220]} notes={str(getattr(r, 'notes', ''))[:160]}",
              flush=True)
    except Exception as e:
        tb = traceback.format_exc().strip().splitlines()
        print(f"{name} [{pid}] EXCEPTION {type(e).__name__}: {str(e)[:200]} @ {tb[-3].strip()[:160]}", flush=True)
