"""Pilot: does a LONG re-run (5000 steps, ~501 frames) let the validated detectors name the known positives?"""
import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
import time, sys
from multiprocessing import Pool
import numpy as np
from census import grammar_g1 as G

def one(item):
    from census import sim, filter as FL
    from epc.phase2a.panel import _detected, _verdict
    name, bits = item; p = G.decode(bits); t = time.time()
    out = sim.run(p, seed0=11, steps=5000, rec_scale=5, S=1); B, FULL = FL._load(); lines = []
    extra = ""
    if "head" in out: extra = f" polarization(late)={np.abs(np.exp(1j * out['head'][0, -50:]).mean(-1)).mean():.2f}"
    if "Nph" in out: extra = f" r(late)={np.abs(np.exp(1j * out['Nph'][0, -20:]).mean(-1)).mean():.2f}"
    if "Cph" in out: extra = f" r(late)={np.abs(np.exp(1j * out['Cph'][0, -20:]).reshape(20, -1).mean(-1)).mean():.2f}"
    for vn, h, md in FL.views(out, p, 0):
        h2 = FL._thin(h, FL.CONFIRM_FRAMES); hits = []
        for pid, fn in FULL.items():
            try:
                r = fn(h2, md)
                if _detected(r): hits.append(f"{pid}:{_verdict(r)}")
            except Exception:
                pass
        lines.append(f"{vn}={hits}")
    return f"{name:32s} {time.time() - t:6.0f}s{extra} | " + " ".join(lines)

if __name__ == "__main__":
    from census.validate import POS, NEG
    items = [(n, G.encode(p)) for n, (p, _) in POS.items()] + [(n, G.encode(NEG[n])) for n in ("flocking at max noise",)]
    with Pool(int(sys.argv[1]) if len(sys.argv) > 1 else 4) as pool:
        for line in pool.imap_unordered(one, items): print(line, flush=True)
