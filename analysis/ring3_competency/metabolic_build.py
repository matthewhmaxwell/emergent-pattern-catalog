"""#25 metabolic-commons foraging (REBUILD from catalog spec).

An agent survives by reading its OWN internal energy and harvesting a REGROWING commons SUSTAINABLY. A
1-D row of food cells regrows only while stock remains (a cell driven to 0 is dead and never recovers --
the tragedy-of-the-commons pressure). Each step costs metabolism; eating restores energy. Competency =
homeostatic regulation: eat when energy is low, but leave enough stock that the commons regrows, so the
agent outlives a greedy strip-miner. Diagnostic (catalog): zero the own-energy channel ('novel') ->
cannot time intake -> collapse; memory-wipe -> collapse. Renders the sustainable run: food row staying
alive (green) + the energy bar holding, vs the greedy baseline that crashes the commons.
"""
import sys, numpy as np, torch, torch.nn as nn
ROOT = "/home/matthewhmaxwell/emergent-pattern-catalog"
sys.path.insert(0, ROOT); sys.path.insert(0, ROOT + "/analysis/ring3_competency")
from gallery.mk_ring3_demo import ASSETS
import render_extra as RE

L, FMAX, REGROW, METAB, EAT, E0, EMAX, T, H = 6, 4, 3, 3.0, 7.0, 20.0, 30.0, 80, 64
torch.manual_seed(0); np.random.seed(0)
NIN = 1 + 3 + L   # energy + food(left,here,right) + position onehot


class Pol(nn.Module):
    def __init__(s, nin=NIN, h=H):
        super().__init__(); s.gru = nn.GRUCell(nin, h); s.pi = nn.Linear(h, 3)  # left, right, eat/stay
    def forward(s, x, h):
        h = s.gru(x, h); return s.pi(h), h


class Commons:
    """1-D regrowing commons. Cell at 0 is dead. Regrows +1 every REGROW steps while level>0."""
    def __init__(s, rng):
        s.f = np.full(L, FMAX, float); s.p = int(rng.integers(0, L)); s.e = E0; s.t = 0
    def obs(s, novel=False):
        fl = s.f[s.p - 1] if s.p > 0 else 0.0
        fr = s.f[s.p + 1] if s.p < L - 1 else 0.0
        pos = np.eye(L)[s.p]
        en = 0.0 if novel else s.e / EMAX
        return np.concatenate([[en], [fl / FMAX, s.f[s.p] / FMAX, fr / FMAX], pos])
    def step(s, a):
        if a == 0: s.p = max(0, s.p - 1)
        elif a == 1: s.p = min(L - 1, s.p + 1)
        else:                                    # eat
            if s.f[s.p] > 0:
                s.f[s.p] -= 1; s.e = min(EMAX, s.e + EAT)
        s.e -= METAB; s.t += 1
        if s.t % REGROW == 0:
            s.f = np.where(s.f > 0, np.minimum(FMAX, s.f + 1), 0.0)
        return s.e > 0


def rollout(net, rng, greedy=False, novel=False, record=False):
    env = Commons(rng); h = torch.zeros(1, H); logps = []; alive = 0; frames = []
    for _ in range(T):
        x = torch.tensor(env.obs(novel), dtype=torch.float32).unsqueeze(0)
        logits, h = net(x, h)
        d = torch.distributions.Categorical(logits=logits)
        a = int(logits.argmax(1)) if greedy else int(d.sample())
        logps.append(d.log_prob(torch.tensor([a])))
        if record:
            frames.append((env.f.copy(), env.p, env.e))
        ok = env.step(a); alive += 1
        if not ok:
            break
    if record:
        frames.append((env.f.copy(), env.p, max(0.0, env.e)))
    return alive / T, torch.stack(logps) if logps else None, frames


def survival(net, rng, trials=200, novel=False):
    return np.mean([rollout(net, rng, greedy=True, novel=novel)[0] for _ in range(trials)])


net = Pol(); opt = torch.optim.Adam(net.parameters(), lr=2e-3)
rng = np.random.default_rng(1); base = 0.3
for it in range(2500):
    frac, logps, _ = rollout(net, rng)
    if logps is None:
        continue
    base = 0.95 * base + 0.05 * frac
    loss = -(logps.sum() * (frac - base))
    opt.zero_grad(); loss.backward(); opt.step()
    if it % 500 == 0 or it == 2499:
        print(f"  it {it:>4}: survival {survival(net, np.random.default_rng(7), 60):.3f}")

surv = survival(net, np.random.default_rng(7))
nov = survival(net, np.random.default_rng(7), novel=True)
# greedy strip-miner baseline: always eat -> crashes the commons
class Greedy:
    def rollout(s, rng):
        env = Commons(rng); alive = 0
        for _ in range(T):
            if not env.step(2): break
            alive += 1
        return alive / T
gr = np.mean([Greedy().rollout(np.random.default_rng(k)) for k in range(200)])
print(f"\nSUSTAINABLE survival: {surv:.3f}   NOVEL (blind to own energy) collapse: {nov:.3f}   "
      f"greedy strip-miner: {gr:.3f}")
torch.save(net.state_dict(), ROOT + "/analysis/ring3_competency/metabolic_net.pt")

# --- render best sustainable run: food row (green by level) + agent + energy bar ---
best = None
for sd in range(300):
    frac, _, frames = rollout(net, np.random.default_rng(200 + sd), greedy=True, record=True)
    if best is None or frac > best[0]:
        best = (frac, frames)
    if frac >= 0.98:
        break
_, frames = best
frames = frames[::2] if len(frames) > 34 else frames        # subsample for a snappy clip


def green(level):
    t = level / FMAX
    if t <= 0: return "#e2e8f0"                              # dead/empty cell = gray
    r = int(0xba + (0x22 - 0xba) * t); g = int(0xe0 + (0x8a - 0xe0) * t); b = int(0xba + (0x3f - 0xba) * t)
    return f"#{r:02x}{g:02x}{b:02x}"


anim = []
for (f, p, e) in frames:
    markers = []
    for i in range(L):                                       # food row at y=1
        markers.append((i, 1, green(f[i]), "s", True, 320))
    ebar = int(round((e / EMAX) * 6))                        # energy bar (0-6 cells) at x=L+1
    for j in range(6):
        col = "#ed8936" if j < ebar else "#edf2f7"
        markers.append((L + 1, j * 0.5, col, "s", True, 160))
    anim.append({"markers": markers, "agents": [(p, 1, "#3182ce")]})
leg = [("food (level)", "s", "#4a9c5d", "#4a9c5d"), ("agent", "o", "#3182ce", "white"),
       ("energy", "s", "#ed8936", "#ed8936")]
print("render:", RE.render_markers(anim, ASSETS + "/ring3_metabolic",
      "Metabolic commons (#25): harvest sustainably, stay fed",
      leg, extent=(-0.6, L + 1.8, -0.6, 3.2), fps=7))
print("SEQDONE")
