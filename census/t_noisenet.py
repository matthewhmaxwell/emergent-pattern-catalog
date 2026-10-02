"""Diagnostic: why is the noise-drowned network flagged? real vs knock-out features, per seed."""
import numpy as np
from census import grammar_g1 as G, sim
from census.knockout import per_seed_view_results, knockout_program, fp_distance
R = G.Rule
p = G.make("N", {"N": 2}, rules=[R("N.COPY", (("p", 0.01),)), R("N.FLIP", (("a", 0), ("b", 1), ("p", 0.5))), R("N.FLIP", (("a", 1), ("b", 0), ("p", 0.5)))])
for seed0 in (12345, 999):
    out = sim.run(p, seed0=seed0); real = per_seed_view_results(out, p)
    kp = knockout_program(p); ko = per_seed_view_results(sim.run(kp, seed0=seed0), kp)
    print("seed0", seed0)
    for s in range(3):
        print("  seed", s, "real emergent", real["N"][s]["emergent"], real["N"][s]["em_score"], "| ko emergent", ko["N"][s]["emergent"], ko["N"][s]["em_score"])
    for k in sorted(real["N"][0]["fp"]):
        rv = [round(real["N"][s]["fp"].get(k, 0), 3) for s in range(3)]; kv = [round(ko["N"][s]["fp"].get(k, 0), 3) for s in range(3)]
        print(f"    {k:20s} real {rv}  ko {kv}")
    x = out["nt"][0, -1]; a = out["adj0"][0].astype(float)
    print("  neighbour same-type share", round(float((a * (x[:, None] == x[None, :])).sum() / a.sum()), 3), "(0.5 = no correlation)")
