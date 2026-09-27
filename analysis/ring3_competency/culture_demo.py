"""#21 cumulative-culture ratchet: known recipe prefix accumulates across generations (cumulative) vs flat when
each generation starts from scratch (no inheritance). Animated line chart."""
import sys, numpy as np
ROOT = "/home/matthewhmaxwell/emergent-pattern-catalog"
sys.path.insert(0, ROOT); sys.path.insert(0, ROOT + "/analysis/ring3_competency")
from gallery.mk_ring3_demo import ASSETS
import render_diagram as RD
import openculture_ppo as OC

gens = 60
tc, _ = OC.run(gens, inherit=True, seed=0)
tn, _ = OC.run(gens, inherit=False, seed=0)
idx = np.linspace(0, gens - 1, 16).astype(int)
xs = [int(i) for i in idx]
series = [("cumulative culture", [tc[i] for i in idx], "#38a169"),
          ("no inheritance", [tn[i] for i in idx], "#888888")]
spr = RD.line_reveal(xs, series, "Cumulative culture (#21): recipe builds across generations",
                     ASSETS + "/ring3_culture", "generation", "known recipe length", ymax=OC.L * 1.08)
print("culture cum", round(tc[-1], 1), "noinh", round(tn[-1], 1), "of", OC.L, spr)
