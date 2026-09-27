"""#14 compositional count-communication: a Speaker watches a K-bit stream, COUNTS the items, and emits a 3-slot
message (one symbol can't label 6 counts, so it must use a multi-slot compositional code); a Listener decodes the
count. Diagram: top row = stream bits (green=item), bottom row = the 3-slot message; title shows count -> decoded."""
import sys, numpy as np, torch
ROOT = "/home/matthewhmaxwell/emergent-pattern-catalog"
sys.path.insert(0, ROOT); sys.path.insert(0, ROOT + "/analysis/ring3_competency")
from gallery.mk_ring3_demo import ASSETS
import render_diagram as RD
import comp_ppo as W

d = torch.load(ROOT + "/analysis/ring3_competency/comp_net.pt", map_location="cpu")
spk = W.Speaker(); spk.load_state_dict(d["spk"]); spk.eval()
lis = W.Listener(); lis.load_state_dict(d["lis"]); lis.eval()
pick = None
for seed in range(400):
    ep = W.episode(spk, lis, 1, seed, greedy=True)
    bits = ep["stream"][0].numpy()[:, 1].astype(int); count = int(bits.sum())
    sym = ep["sym"][0].numpy(); pred = int(ep["pred"].numpy()[0])
    if pred == count and 2 <= count <= 4:
        pick = (seed, bits, sym, count, pred); break
if pick is None:
    ep = W.episode(spk, lis, 1, 0, greedy=True); bits = ep["stream"][0].numpy()[:, 1].astype(int)
    pick = (0, bits, ep["sym"][0].numpy(), int(bits.sum()), int(ep["pred"].numpy()[0]))
seed, bits, sym, count, pred = pick
SYMCOL = ["#4363d8", "#dd6b20"]
full_stream = [(j, 1, ("#38a169" if bits[j] else "#e2e8f0"), str(int(bits[j])), ("white" if bits[j] else "#333")) for j in range(W.K)]
steps = []
for k in range(1, W.K + 1):
    steps.append(full_stream[:k])
for m in range(1, W.M + 1):
    steps.append(full_stream + [(j, 0, SYMCOL[int(sym[j])], str(int(sym[j])), "white") for j in range(m)])
title = f"Compositional comms (#14): count {count} -> 3-slot message -> decoded {pred}" + (" OK" if pred == count else "")
leg = [("stream item", "#38a169"), ("no item", "#e2e8f0"), ("message slot", "#dd6b20")]
print("comp seed", seed, "count", count, "pred", pred,
      RD.timeline(steps, title, ASSETS + "/ring3_comp", leg))
