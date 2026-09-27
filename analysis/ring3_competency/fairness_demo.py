"""#13 fair alternation: two symmetric agents play an iterated claim/yield game; the fair solution is TURN-TAKING
(one claims each round, alternating). Timeline: rows=2 agents, cols=rounds; C=claim (gold), Y=yield (gray)."""
import sys, numpy as np, torch
ROOT = "/home/matthewhmaxwell/emergent-pattern-catalog"
sys.path.insert(0, ROOT); sys.path.insert(0, ROOT + "/analysis/ring3_competency")
from gallery.mk_ring3_demo import ASSETS
import render_diagram as RD
import fairness_ppo as W

net = W.Pol(); net.load_state_dict(torch.load(ROOT + "/analysis/ring3_competency/fairness_net.pt", map_location="cpu")); net.eval()
# MUST sample (not argmax): symmetric agents on identical obs can only break symmetry stochastically.
best = None
for seed in range(600):
    O, Aa, SY, LP, V, R, st = W.rollout(net, 1, seed, greedy=False)
    A0 = Aa[0]
    one_claim = float(((A0[0] + A0[1]) == 1).mean())            # exactly one claims each round
    flips = float(((A0[:, 1:] != A0[:, :-1]).mean()))           # agents alternate
    score = one_claim + 0.3 * flips
    if best is None or score > best[0]:
        best = (score, seed, A0, one_claim)
    if one_claim >= 0.92:
        break
_, seed, Aarr, oc = best
cols = W.ROUNDS
steps = []
for k in range(1, cols + 1):
    fr = []
    for r in range(k):
        for i in range(2):
            claim = Aarr[i, r] == 1
            fr.append((r, 1 - i, "#d69e2e" if claim else "#e2e8f0", "C" if claim else "Y", "#333"))
    steps.append(fr)
leg = [("claim", "#d69e2e"), ("yield", "#e2e8f0")]
print("fairness seed", seed, "score", round(float(best[0]), 2),
      RD.timeline(steps, "Fair alternation (#13): agents take turns claiming", ASSETS + "/ring3_fairness", leg))
