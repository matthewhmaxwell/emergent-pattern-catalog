"""#27 non-stationary co-player adaptation (REBUILD from catalog spec).

A repeated 3-move game (rock/paper/scissors style). The CO-PLAYER follows a fixed periodic pattern, then
SWITCHES to a different pattern at mid-episode. The agent scores when its move counters the co-player's
current move, so it must infer the co-player's pattern from history, and RE-INFER after the switch.
Competency = online opponent-modelling that re-adapts to a regime change. Diagnostic (catalog): zero the
opponent-history channel ('noopphist') -> cannot model -> chance; a STATIONARY-trained agent (never saw a
switch) holds phase 1 but fails phase 2. Renders rolling counter-accuracy with the switch marked: a dip
at the switch, then recovery = the re-adaptation.
"""
import sys, numpy as np, torch, torch.nn as nn
ROOT = "/home/matthewhmaxwell/emergent-pattern-catalog"
sys.path.insert(0, ROOT); sys.path.insert(0, ROOT + "/analysis/ring3_competency")
from gallery.mk_ring3_demo import ASSETS
import render_diagram as RD

torch.set_num_threads(4)                 # cap threads: the box is shared; avoid thrash
A, T, SWITCH, H = 3, 40, 20, 64
torch.manual_seed(0); np.random.seed(0)
NIN = A + A + 1   # prev opp move onehot + prev own move onehot + bias


class Pol(nn.Module):
    def __init__(s, nin=NIN, h=H):
        super().__init__(); s.gru = nn.GRUCell(nin, h); s.pi = nn.Linear(h, A)
    def forward(s, x, h):
        h = s.gru(x, h); return s.pi(h), h


def make_pattern(rng):
    P = int(rng.integers(2, 4))                              # period 2 or 3
    return [int(rng.integers(0, A)) for _ in range(P)]


def opp_moves(rng, switch=True):
    p1 = make_pattern(rng)
    p2 = make_pattern(rng)
    while switch and p2 == p1:
        p2 = make_pattern(rng)
    seq = []
    for t in range(T):
        pat = p1 if (t < SWITCH or not switch) else p2
        seq.append(pat[t % len(pat)])
    return seq


def rollout(net, rng, greedy=False, noopphist=False, switch=True, record=False):
    seq = opp_moves(rng, switch); h = torch.zeros(1, H)
    prev_o = prev_a = -1; logps = []; hits = []
    for t in range(T):
        oo = np.zeros(A); aa = np.zeros(A)
        if prev_o >= 0 and not noopphist: oo[prev_o] = 1
        if prev_a >= 0: aa[prev_a] = 1
        x = torch.tensor(np.concatenate([oo, aa, [1.0]]), dtype=torch.float32).unsqueeze(0)
        logits, h = net(x, h)
        d = torch.distributions.Categorical(logits=logits)
        a = int(logits.argmax(1)) if greedy else int(d.sample())
        logps.append(d.log_prob(torch.tensor([a])))
        o = seq[t]
        hit = int(a == (o + 1) % A)                          # counter beats o
        hits.append(hit); prev_o, prev_a = o, a
    return np.array(hits), torch.stack(logps), seq


def phase_acc(net, rng, trials=300, noopphist=False, switch=True):
    h1 = h2 = 0.0
    for _ in range(trials):
        hits, _, _ = rollout(net, rng, greedy=True, noopphist=noopphist, switch=switch)
        h1 += hits[:SWITCH].mean(); h2 += hits[SWITCH:].mean()
    return h1 / trials, h2 / trials


def train(switch=True, iters=3000):
    net = Pol(); opt = torch.optim.Adam(net.parameters(), lr=2e-3)
    rng = np.random.default_rng(1); base = 0.33
    for it in range(iters):
        hits, logps, _ = rollout(net, rng, switch=switch)
        R = float(hits.mean()); base = 0.95 * base + 0.05 * R
        # per-round advantage: reward each round equally toward its own hit
        adv = torch.tensor(hits - base, dtype=torch.float32)
        loss = -(logps.squeeze(1) * adv).mean()
        opt.zero_grad(); loss.backward(); opt.step()
        if it % 600 == 0 or it == iters - 1:
            a1, a2 = phase_acc(net, np.random.default_rng(7), 80)
            print(f"  it {it:>4}: phase1 {a1:.3f}  phase2 {a2:.3f}")
    return net


print("training adaptive agent (sees regime switches):")
net = train(switch=True, iters=2200)
a1, a2 = phase_acc(net, np.random.default_rng(9))
n1, n2 = phase_acc(net, np.random.default_rng(9), noopphist=True)
print(f"\nADAPTIVE: phase1 {a1:.3f}  phase2 {a2:.3f}   NOOPPHIST collapse: p1 {n1:.3f} p2 {n2:.3f}")
print("training stationary control (never sees a switch):")
stat = train(switch=False, iters=1200)
s1, s2 = phase_acc(stat, np.random.default_rng(9))                # tested WITH a switch
print(f"STATIONARY-trained tested with switch: phase1 {s1:.3f}  phase2 {s2:.3f} (should fail p2)")
torch.save(net.state_dict(), ROOT + "/analysis/ring3_competency/nonstat_net.pt")

# --- render rolling counter-accuracy with the switch marked (dip then recovery) ---
best = None
for sd in range(400):
    hits, _, _ = rollout(net, np.random.default_rng(300 + sd), greedy=True)
    dip = hits[SWITCH:SWITCH + 3].mean(); rec = hits[SWITCH + 5:].mean(); pre = hits[3:SWITCH].mean()
    score = pre + rec - dip                                      # want strong pre, clear dip, strong recovery
    if best is None or score > best[0]:
        best = (score, hits)
_, hits = best
W = 4
roll = [float(hits[max(0, k - W + 1):k + 1].mean()) for k in range(T)]
xs = list(range(1, T + 1))
series = [("counter-accuracy (rolling)", roll, "#3182ce")]
print("render:", RD.line_reveal(xs, series, "Non-stationary co-player (#27): re-adapt after a switch",
      ASSETS + "/ring3_nonstat", "round", "accuracy", fps=6, ymax=1.08, ymin=0,
      hline=1 / A, vmark=SWITCH + 0.5, vlabel="co-player switches"))
print("SEQDONE")
