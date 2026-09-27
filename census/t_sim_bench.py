"""Pilot sanity run: every benchmark model, short run, timing + one textbook order parameter."""
import time, numpy as np
from census.bench_models import BENCH
from census import sim

def order(name, o):
    if "C" in o and name in ("voter", "majority coarsen", "Nowak-May PD"):
        C = o["C"][:, -1]; return "majority-type frac %.2f | interface frac %.2f" % (
            np.mean([max(np.bincount(c.ravel(), minlength=2)) / c.size for c in C]),
            np.mean(C != np.roll(C, 1, -1)))
    if name == "Schelling":
        C = o["C"]; occ = C != 0
        same = lambda X: np.mean([(np.roll(X, 1, -1) == X)[X != 0].mean() for X in X])
        return "same-neighbour frac start %.2f -> end %.2f" % (same(C[:, 0]), same(C[:, -1]))
    if name in ("cyclic CA spirals", "Greenberg-Hastings"):
        C = o["C"]; ch = (C[:, -1] != C[:, -2]).mean(); return "active frac (last frame) %.2f" % ch
    if name == "BTW sandpile":
        a = o["aval"][:, -500:].ravel(); a = a[a > 0]
        return "avalanches %d, max %d, mean %.1f" % (len(a), a.max() if len(a) else 0, a.mean() if len(a) else 0)
    if "Kuramoto" in name:
        ph = o["Cph"] if "Cph" in o else o["Nph"]; R = lambda x: np.abs(np.exp(1j * x).reshape(x.shape[0], -1).mean(-1)).mean()
        return "order R start %.2f -> end %.2f" % (R(ph[:, 0]), R(ph[:, -1]))
    if name in ("Vicsek", "MIPS (quorum)", "RPS agents", "chemotaxis (KS)"):
        h = o["head"]; pol = lambda x: np.abs(np.exp(1j * x).mean(-1)).mean()
        pos = o["pos"][:, -1]; H = np.mean([np.histogram2d(p[:, 0], p[:, 1], bins=16)[0].var() for p in pos])
        at = o["at"][:, -1]; ty = np.mean([np.bincount(a, minlength=3).max() / a.size for a in at])
        return "polarization %.2f->%.2f | density var %.1f | top type frac %.2f" % (pol(h[:, 0]), pol(h[:, -1]), H, ty)
    if name == "co-evolving voter":
        return "majority frac %.2f | edges %d" % (np.mean([np.bincount(n, minlength=2).max() / n.size for n in o["nt"][:, -1]]), o["adj"][:, -1].sum() // 2 // o["adj"].shape[0])
    if name == "Langton ant(s)":
        return "cells flipped vs start %.2f" % (o["C"][:, -1] != o["C"][:, 0]).mean()
    if name == "Gray-Scott":
        F = o["F"][:, -1]; return "v mean %.3f std %.3f | u std %.3f" % (F[:, 1].mean(), F[:, 1].std(), F[:, 0].std())
    return ""

for name, p in BENCH.items():
    t = time.time(); o = sim.run(p, seed0=1, steps=1000); dt = time.time() - t
    m = o["meta"]; print(f"{name:20s} {dt:6.2f}s  steps={m['steps_run']:4d} frozen={m['frozen_at']} unstable={m['unstable']} | {order(name, o)}", flush=True)
