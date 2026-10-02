"""Diagnostic: reproduce the flagged noise-drowned network negative of the final gate (same seed) and show why."""
import numpy as np
from census import grammar_g1 as G, sim
from census.runner import seed_of
from census.knockout import per_seed_view_results, knockout_program, interaction_driven, ORDER
R = G.Rule
p = G.make("N", {"N": 2}, rules=[R("N.COPY", (("p", 0.03),)), R("N.FLIP", (("a", 0), ("b", 1), ("p", 0.3))), R("N.FLIP", (("a", 1), ("b", 0), ("p", 0.3)))])
seed0 = seed_of(G.encode(p)) + 161803
real = per_seed_view_results(sim.run(p, seed0=seed0), p); kp = knockout_program(p); ko = per_seed_view_results(sim.run(kp, seed0=seed0), kp)
print("driven:", interaction_driven(real, ko))
for s in range(3):
    print(f" seed {s}: real emergent={real['N'][s]['emergent']} evidence={real['N'][s]['evidence']:.3f} | ko emergent={ko['N'][s]['emergent']} evidence={ko['N'][s]['evidence']:.3f}")
for f in ORDER:
    if f in real["N"][0]["fp"]:
        print(f"   {f:20s} real {[round(real['N'][s]['fp'][f], 3) for s in range(3)]}  ko {[round(ko['N'][s]['fp'][f], 3) for s in range(3)]}")
