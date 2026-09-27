"""#23 group-size regulation (REBUILD from catalog spec).

n agents each choose active/inactive each round; the SHARED reward peaks when the ACTIVE COUNT equals a
target K*(t) that switches every few rounds across [1, n-1]. Competency = agents use IDENTITY to take
stable complementary roles whose sum tracks K*, and RE-FORM the assignment when K* switches. Diagnostic
(from catalog): drop the id channel ('noid') -> all agents identical -> they can only match p*N in
expectation, never the exact moving target -> collapse. This script trains a shared recurrent policy,
prints the full-vs-noid gap as the interventional check, then renders target K*(t) vs actual count.
"""
import sys, numpy as np, torch, torch.nn as nn
ROOT = "/home/matthewhmaxwell/emergent-pattern-catalog"
sys.path.insert(0, ROOT); sys.path.insert(0, ROOT + "/analysis/ring3_competency")
from gallery.mk_ring3_demo import ASSETS
import render_diagram as RD

N, T, SWITCH, H = 6, 24, 4, 64
torch.manual_seed(0); np.random.seed(0)
NIN = N + 3   # onehot id (N) + K*/N + last_count/N + bias


class Pol(nn.Module):
    def __init__(s, nin=NIN, h=H):
        super().__init__(); s.gru = nn.GRUCell(nin, h); s.pi = nn.Linear(h, 2)
    def forward(s, x, h):
        h = s.gru(x, h); return s.pi(h), h


def make_targets(rng):
    ks, k = [], int(rng.integers(1, N))
    for t in range(T):
        if t and t % SWITCH == 0:
            k = int(rng.integers(1, N))
        ks.append(k)
    return ks


def rollout(net, rng, greedy=False, noid=False):
    ks = make_targets(rng); h = torch.zeros(N, H); last = 0.0
    logps, counts, rews = [], [], []
    idm = np.zeros((N, N)) if noid else np.eye(N)
    for t in range(T):
        x = np.concatenate([idm, np.full((N, 1), ks[t] / N),
                            np.full((N, 1), last / N), np.ones((N, 1))], 1)
        logits, h = net(torch.tensor(x, dtype=torch.float32), h)
        d = torch.distributions.Categorical(logits=logits)
        a = logits.argmax(1) if greedy else d.sample()
        logps.append(d.log_prob(a))
        cnt = int(a.sum().item()); last = cnt; counts.append(cnt)
        rews.append(1.0 - abs(cnt - ks[t]) / N)
    return ks, counts, rews, torch.stack(logps)   # logps: (T, N)


def exact_match(net, rng, trials=200, noid=False):
    hit = 0
    for _ in range(trials):
        ks, counts, _, _ = rollout(net, rng, greedy=True, noid=noid)
        hit += np.mean([c == k for c, k in zip(counts, ks)])
    return hit / trials


net = Pol(); opt = torch.optim.Adam(net.parameters(), lr=3e-3)
rng = np.random.default_rng(1); base = 0.5
for it in range(2000):
    ks, counts, rews, logps = rollout(net, rng)
    R = torch.tensor(rews, dtype=torch.float32)          # (T,)
    base = 0.95 * base + 0.05 * R.mean().item()
    adv = (R - base).unsqueeze(1)                         # (T,1) broadcast over agents
    loss = -(logps * adv).mean()
    opt.zero_grad(); loss.backward(); opt.step()
    if it % 400 == 0 or it == 1999:
        print(f"  it {it:>4}: exact-match {exact_match(net, np.random.default_rng(7), 60):.3f}")

full = exact_match(net, np.random.default_rng(7))
noid = exact_match(net, np.random.default_rng(7), noid=True)
print(f"\nFULL exact-match regulation: {full:.3f}   NOID (drop identity) collapse: {noid:.3f}")
torch.save(net.state_dict(), ROOT + "/analysis/ring3_competency/group_reg_net.pt")

# --- render a clean episode: target K*(t) vs actual active count, tracking ---
rng = np.random.default_rng(3)
best = None
for s in range(400):
    ks, counts, _, _ = rollout(net, np.random.default_rng(100 + s), greedy=True)
    err = np.mean([abs(c - k) for c, k in zip(counts, ks)])
    if best is None or err < best[0]:
        best = (err, ks, counts)
    if err == 0:
        break
_, ks, counts = best
xs = list(range(1, T + 1))
series = [("target K*", [float(k) for k in ks], "#e53e3e"),
          ("active count", [float(c) for c in counts], "#3182ce")]
styles = [dict(ls="--", lw=3.2, marker="", zorder=2, drawstyle="steps-mid"),   # target: dashed step underneath
          dict(ls="-", lw=1.6, marker="o", ms=4.0, zorder=4)]                   # count: dots ride on top
print("render:", RD.line_reveal(xs, series, "Group-size regulation (#23): count tracks moving target",
                                 ASSETS + "/ring3_group_reg", "round", "# active agents",
                                 fps=6, ymax=N + 0.5, ymin=0, styles=styles))
print("SEQDONE")
