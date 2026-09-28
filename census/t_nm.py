"""Pilot check: Nowak-May spatial PD dynamics (cooperator fraction over time) vs the literature."""
import numpy as np
from census import grammar_g1 as G, sim
R = G.Rule
for g in [(1.9, 0.0), (1.3, 0.0)]:
    p = G.make("C", {"C": 2}, rules=[R("C.GAME", (("g", g),))]); out = sim.run(p, seed0=3, steps=600)
    C = out["C"]; coop = (C == 0).mean((2, 3))                     # (S, frames); type 0 = cooperate
    ch = (C[:, 1:] != C[:, :-1]).mean((2, 3))
    print(f"T={g[0]} S={g[1]}: coop fraction at frames 0,5,25,100,300 = {np.round(coop[0, [0, 5, 25, 100, 300]], 3)} | "
          f"late activity {ch[:, -50:].mean():.3f} | seeds final coop {np.round(coop[:, -1], 2)}")
