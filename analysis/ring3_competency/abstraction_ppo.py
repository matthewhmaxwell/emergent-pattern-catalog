"""Track-1 probe (audit gap #1): ABSTRACTION / systematic generalization. The agent must navigate to the
goal whose COLOUR matches a cue colour. Colours are random feature VECTORS (not one-hot), so 'match'
requires computing similarity(cue, goal) -- a RELATION. It is trained on a SUBSET of colours and tested on
HELD-OUT colours. If it learned the relation (compare cue to each goal), it generalizes to unseen colours
(held-out accuracy high = abstraction). If it memorized colour->goal instances, held-out collapses to
chance (= memorization, the systematic-generalization failure nets often show). The train->held-out gap IS
the competency; both outcomes are informative (genuine abstraction, or a clean negative).

Held-out colours appear as DISTRACTORS during training (their input dims are seen) but never as the
cue-match, so the test isolates whether the EQUALITY relation transfers.

Run on VPS epc-venv: python abstraction_ppo.py [--iters 800]
"""
import numpy as np, sys, json, os, torch, torch.nn as nn

N = 7; C = 8; D = 8; G = 3; T = 16; H = 96; NACT = 5
TRAIN = list(range(6)); HELD = [6, 7]                 # train colours vs held-out colours
ODIM = D + G * (D + 2)
DIRS = np.array([(0, 0), (1, 0), (-1, 0), (0, 1), (0, -1)])
COLVEC = np.random.default_rng(12345).standard_normal((C, D)).astype(np.float32)   # fixed colour vectors
torch.set_num_threads(8)


class Pol(nn.Module):
    def __init__(self):
        super().__init__()
        self.body = nn.Sequential(nn.Linear(ODIM, H), nn.Tanh(), nn.Linear(H, H), nn.Tanh())
        self.pi = nn.Linear(H, NACT); self.v = nn.Linear(H, 1)

    def forward(self, x):
        z = self.body(x); return self.pi(z), self.v(z).squeeze(-1)


class VecAbs:
    def __init__(self, B, seed, colours):
        self.B, self.colours = B, colours; self.rng = np.random.default_rng(seed); self.reset()

    def reset(self):
        B = self.B
        self.cue = np.array([self.rng.choice(self.colours) for _ in range(B)])     # cue colour (from given set)
        self.gcol = np.zeros((B, G), int); self.gpos = np.zeros((B, G, 2), int); self.ap = np.zeros((B, 2), int)
        for b in range(B):
            others = [c for c in range(C) if c != self.cue[b]]
            self.rng.shuffle(others)
            cols = [self.cue[b]] + others[:G - 1]; self.rng.shuffle(cols)        # one matching + distractors
            self.gcol[b] = cols
            cells = self.rng.permutation(N * N)[:G + 1]
            for g in range(G): self.gpos[b, g] = (cells[g] % N, cells[g] // N)
            self.ap[b] = (cells[G] % N, cells[G] // N)
        self.match = np.array([list(self.gcol[b]).index(self.cue[b]) for b in range(B)])
        self.done = np.zeros(B, bool); self.win = np.zeros(B, bool)
        return self.obs()

    def obs(self):
        B = self.B; out = np.zeros((B, ODIM), np.float32); o = 0
        out[:, o:o + D] = COLVEC[self.cue]; o += D
        for g in range(G):
            out[:, o:o + D] = COLVEC[self.gcol[:, g]]; o += D
            out[:, o:o + 2] = np.sign(self.gpos[:, g] - self.ap); o += 2
        return out

    def step(self, a):
        B = self.B; rew = np.full(B, -0.02, np.float32)
        # shaping: encourage closing on the NEAREST goal (does not reveal which is the match)
        d_before = np.min(np.abs(self.gpos - self.ap[:, None, :]).sum(2), 1)
        for b in range(B):
            if self.done[b]: rew[b] = 0.0; continue
            np_ = self.ap[b] + DIRS[a[b]]
            if 0 <= np_[0] < N and 0 <= np_[1] < N: self.ap[b] = np_
            for g in range(G):
                if tuple(self.ap[b]) == tuple(self.gpos[b, g]):
                    win = (g == self.match[b]); rew[b] += 1.0 if win else -0.5
                    self.done[b] = True; self.win[b] = win; break
        d_after = np.min(np.abs(self.gpos - self.ap[:, None, :]).sum(2), 1)
        rew += 0.03 * (d_before - d_after) * (~self.done)
        return rew

    def acc(self): return float(self.win.mean())


def gae(R, V, TERM, ACT, gamma=0.98, lam=0.95):
    B, Tt = R.shape; adv = np.zeros((B, Tt), np.float32); last = np.zeros(B, np.float32)
    for t in reversed(range(Tt)):
        nextv = V[:, t + 1] if t + 1 < Tt else np.zeros(B, np.float32)
        nonterm = 1.0 - TERM[:, t]
        delta = R[:, t] + gamma * nextv * nonterm - V[:, t]
        last = delta + gamma * lam * nonterm * last; adv[:, t] = last * ACT[:, t]
    return adv, adv + V


def rollout(net, env, greedy=False):
    B = env.B; obs = env.obs()
    O = np.zeros((B, T, ODIM), np.float32); A = np.zeros((B, T), int); LP = np.zeros((B, T), np.float32)
    V = np.zeros((B, T), np.float32); R = np.zeros((B, T), np.float32); TERM = np.zeros((B, T), np.float32); ACT = np.zeros((B, T), np.float32)
    for t in range(T):
        with torch.no_grad():
            lg, v = net(torch.from_numpy(obs))
        d = torch.distributions.Categorical(logits=lg); a = lg.argmax(1) if greedy else d.sample()
        O[:, t] = obs; A[:, t] = a.numpy(); LP[:, t] = d.log_prob(a).numpy(); V[:, t] = v.numpy()
        ACT[:, t] = (~env.done).astype(np.float32)
        before = env.done.copy(); r = env.step(a.numpy()); R[:, t] = r
        TERM[:, t] = (env.done & ~before).astype(np.float32); obs = env.obs()
    return O, A, LP, V, R, TERM, ACT


def train(iters, seed=0, B=256):
    torch.manual_seed(seed); net = Pol(); opt = torch.optim.Adam(net.parameters(), lr=3e-3)
    for it in range(iters):
        env = VecAbs(B, 1000 + it, TRAIN)
        O, A, LP, V, R, TERM, ACT = rollout(net, env)
        adv, ret = gae(R, V, TERM, ACT); m = ACT.sum()
        adv = (adv - (adv * ACT).sum() / m) / (np.sqrt(((adv - (adv * ACT).sum() / m) ** 2 * ACT).sum() / m) + 1e-8)
        Ot, At, LPt = torch.from_numpy(O), torch.from_numpy(A), torch.from_numpy(LP)
        advt, rett, actt = torch.from_numpy(adv), torch.from_numpy(ret), torch.from_numpy(ACT)
        for _ in range(4):
            lg, v = net(Ot.reshape(-1, ODIM)); lg = lg.reshape(B, T, NACT); v = v.reshape(B, T)
            d = torch.distributions.Categorical(logits=lg)
            ratio = torch.exp(d.log_prob(At) - LPt)
            s1 = ratio * advt; s2 = torch.clamp(ratio, 0.8, 1.2) * advt
            loss = -(torch.min(s1, s2) * actt).sum() / actt.sum() + 0.5 * (((v - rett) ** 2) * actt).sum() / actt.sum() - 0.01 * (d.entropy() * actt).sum() / actt.sum()
            opt.zero_grad(); loss.backward(); nn.utils.clip_grad_norm_(net.parameters(), 0.5); opt.step()
        if it % 100 == 0 or it == iters - 1:
            print(f"  iter {it:>3}: train-acc {_acc(net, 7000+it, TRAIN):.2f}  held-out-acc {_acc(net, 7700+it, HELD):.2f}", flush=True)
    return net


def _acc(net, seed, colours, B=500):
    env = VecAbs(B, seed, colours); rollout(net, env, greedy=True); return env.acc()


if __name__ == "__main__":
    a = sys.argv
    iters = int(a[a.index("--iters") + 1]) if "--iters" in a else 800
    here = __file__.rsplit("/", 1)[0]
    print(f"ABSTRACTION probe: go to the goal whose colour MATCHES the cue (random colour vectors); train {TRAIN} test {HELD}", flush=True)
    net = train(iters)
    tr = np.mean([_acc(net, 9000 + k, TRAIN) for k in range(3)])
    ho = np.mean([_acc(net, 9100 + k, HELD) for k in range(3)])
    chance = 1.0 / G
    print(f"\nFINAL: train-acc {tr:.2f} | HELD-OUT-acc {ho:.2f} (chance ~{chance:.2f})", flush=True)
    # 3-way band: memorization = held-out AT chance (pure per-instance lookup); genuine = held-out clears
    # the strong bar; PARTIAL = held-out well above chance but below the bar (genuine-but-incomplete transfer).
    strong = tr >= 0.7 and ho >= 0.6
    partial = tr >= 0.7 and ho >= chance + 0.15 and not strong
    if strong:      verdict = "GENUINE ABSTRACTION (relational generalization)"
    elif partial:   verdict = "PARTIAL ABSTRACTION (genuine but incomplete transfer)"
    elif tr >= 0.7: verdict = "MEMORIZATION (no generalization; held-out at chance)"
    else:           verdict = "under-solved"
    print(f"VERDICT: {verdict} (train {tr:.2f}, held-out {ho:.2f}, gap {tr-ho:+.2f})", flush=True)
    torch.save(net.state_dict(), f"{here}/abstraction_net.pt")
    json.dump({"train_acc": round(float(tr), 2), "heldout_acc": round(float(ho), 2),
               "chance": round(float(chance), 2), "strong_abstraction": bool(strong),
               "partial_abstraction": bool(partial)}, open(f"{here}/abstraction_result.json", "w"))
