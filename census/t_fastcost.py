"""Pilot: cost + firing of the FAST battery (reduced nulls) as an extra screen channel, per view type."""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import time
from census import grammar_g1 as G, sim, filter as FL
from epc.phase2a.panel import _detected, _verdict
from census.validate import POS, NEG
F = FL._fast_fns()
cases = [(n, p, False) for n, (p, _) in POS.items()] + [(n, p, False) for n, p in NEG.items()] + [(n + " [scr]", p, True) for n, (p, _) in list(POS.items())[:6]]
for name, p, shuf in cases:
    out = sim.run(p, seed0=11, null_shuffle=shuf)
    for vn, h, md in FL.views(out, p, 0):
        if vn == "C.aval": continue
        hh = FL._screen_hist(vn, h, "adj" in out); t = time.time(); fired = []
        for pid, fn in F.items():
            try:
                r = fn(hh, md)
                if _detected(r): fired.append(f"{pid}:{_verdict(r)}")
            except Exception:
                pass
        sc = FL.screen(hh)
        print(f"{name:34s} {vn:8s} fast-battery {time.time() - t:5.1f}s fired={fired} | generic em={sc['em_score']} complex={sc['is_complex']}", flush=True)
