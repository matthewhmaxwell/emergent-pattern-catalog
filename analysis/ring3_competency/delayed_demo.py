"""#5 delayed gratification: a reward value rises then vanishes (unknown peak height); the agent must cash out
near the peak (waiting too long loses it). Line: reward over time with the agent's cash-out marked."""
import sys, random
ROOT = "/home/matthewhmaxwell/emergent-pattern-catalog"
sys.path.insert(0, ROOT); sys.path.insert(0, ROOT + "/analysis/ring3_competency")
from gallery.mk_ring3_demo import ASSETS
import render_diagram as RD
import evolve_delayed as ED                      # evolves on import -> ED.best (ns=12 stateful)
best = ED.best; NS = 12
pick = None
for k in range(500):
    r = random.Random(k); c, h = ED.curve(r); s = 0; cash = None
    for t, v in enumerate(c):
        b = min(9, int(v * 9 + 1e-9)); a, s = best[b * NS + s]; s %= NS
        if a == 1:
            cash = t; break
    if cash is not None and cash >= 2 and 5 <= len(c) <= 28:
        payoff = c[cash] / h
        if payoff >= 0.85:
            pick = (k, list(c), cash, payoff); break
if pick is None:
    r = random.Random(0); c, h = ED.curve(r); pick = (0, list(c), max(1, len(c) // 2), 0.5)
k, c, cash, payoff = pick
xs = list(range(len(c)))
print("delayed seed", k, "cash", cash, "payoff", round(payoff, 2),
      RD.line_reveal(xs, [("reward available", c, "#dd6b20")],
                     f"Delayed gratification (#5): cash out near the peak (payoff {payoff:.2f})",
                     ASSETS + "/ring3_delayed", "time", "reward value", vmark=cash, vlabel="cash out", ymax=max(c) * 1.18))
