"""Track-1 probe (audit gap #3, the LAST untouched audit cell): IMITATION / SOCIAL LEARNING. The learner must
reach a CORRECT goal that is HIDDEN from it (not in its observation). A DEMONSTRATOR agent (scripted, knows
the goal) walks greedily toward the correct goal. The learner sees the demonstrator's position + last move
and the G goal positions, but NOT which goal is correct. The only way to score is to READ the demonstrator's
behaviour to identify the goal, then go there = goal-emulation / observational learning.

Diagnostics (the demand):
  - FROZEN demonstrator (sits at centre, reveals nothing) -> learner has no goal info -> should fall to chance
    (1/G). This is the core social-dependence test: without the other agent's behaviour the goal is unfindable.
  - MISLEADING demonstrator (walks to a WRONG goal) -> a learner that genuinely COPIES the demonstrator's
    choice is led to the wrong goal -> accuracy collapses. Confirms the signal used is the demonstrator's
    choice, not some other cue.

Run on VPS epc-venv: python imitation_ppo.py [--iters 800]
"""
import numpy as np, sys, json, os, torch, torch.nn as nn

N = 7; G = 3; T = 16; H = 96; NACT = 5
ODIM = 2 + G * 2 + 2 + NACT             # self xy(2), each goal xy(G*2), demonstrator xy(2), demo last-move(5)
DIRS = np.array([(0, 0), (1, 0), (-1, 0), (0, 1), (0, -1)])
torch.set_num_threads(8)


class Pol(nn.Module):
    def __init__(self):
        super().__init__()
        self.body = nn.Sequential(nn.Linear(ODIM, H), nn.Tanh(), nn.Linear(H, H), nn.Tanh())
        self.pi = nn.Linear(H, NACT); self.v = nn.Linear(H, 1)

    def forward(self, x):
        z = self.body(x); return self.pi(z), self.v(z).squeeze(-1)


def greedy_step(pos, target):
    """One Manhattan-greedy move index toward target (x first, then y)."""
    d = target - pos
    if d[0] != 0: return 1 if d[0] > 0 else 2
    if d[1] != 0: return 3 if d[1] > 0 else 4
    return 0


class VecImit:
    def __init__(self, B, seed, mode="normal"):
        self.B, self.mode = B, mode; self.rng = np.random.default_rng(seed); self.reset()

    def reset(self):
        B = self.B
        self.gpos = np.zeros((B, G, 2), int); self.ap = np.zeros((B, 2), int); self.dp = np.zeros((B, 2), int)
        self.correct = np.zeros(B, int); self.demo_goal = np.zeros(B, int); self.dmove = np.zeros(B, int)
        for b in range(B):
            cells = self.rng.permutation(N * N)[:G + 2]
            for g in range(G): self.gpos[b, g] = (cells[g] % N, cells[g] // N)
            self.ap[b] = (cells[G] % N, cells[G] // N)
            self.dp[b] = (N // 2, N // 2)                       # demonstrator starts at centre
            self.correct[b] = self.rng.integers(G)
            if self.mode == "mislead":                          # demonstrator walks to a WRONG goal
                wrong = [g for g in range(G) if g != self.correct[b]]
                self.demo_goal[b] = wrong[self.rng.integers(len(wrong))]
            else:
                self.demo_goal[b] = self.correct[b]
        self.done = np.zeros(B, bool); self.win = np.zeros(B, bool)
        return self.obs()

    def obs(self):
        B = self.B; out = np.zeros((B, ODIM), np.float32); o = 0
        out[:, o:o + 2] = self.ap / (N - 1); o += 2
        for g in range(G):
            out[:, o:o + 2] = self.gpos[:, g] / (N - 1); o += 2
        out[:, o:o + 2] = self.dp / (N - 1); o += 2
        out[np.arange(B), o + self.dmove] = 1.0; o += NACT
        return out

    def step(self, a):
        B = self.B; rew = np.full(B, -0.02, np.float32)
        d_before = np.min(np.abs(self.gpos - self.ap[:, None, :]).sum(2), 1)
        # move demonstrator (unless frozen)
        if self.mode != "frozen":
            for b in range(B):
                mv = greedy_step(self.dp[b], self.gpos[b, self.demo_goal[b]])
                self.dmove[b] = mv; npd = self.dp[b] + DIRS[mv]
                if 0 <= npd[0] < N and 0 <= npd[1] < N: self.dp[b] = npd
        # move learner
        for b in range(B):
            if self.done[b]: rew[b] = 0.0; continue
            np_ = self.ap[b] + DIRS[a[b]]
            if 0 <= np_[0] < N and 0 <= np_[1] < N: self.ap[b] = np_
            for g in range(G):
                if tuple(self.ap[b]) == tuple(self.gpos[b, g]):
                    win = (g == self.correct[b]); rew[b] += 1.0 if win else -0.3
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
        before = env.done.copy(); R[:, t] = env.step(a.numpy())
        TERM[:, t] = (env.done & ~before).astype(np.float32); obs = env.obs()
    return O, A, LP, V, R, TERM, ACT


def train(iters, seed=0, B=256):
    torch.manual_seed(seed); net = Pol(); opt = torch.optim.Adam(net.parameters(), lr=3e-3)
    for it in range(iters):
        env = VecImit(B, 1000 + it)
        O, A, LP, V, R, TERM, ACT = rollout(net, env)
        adv, ret = gae(R, V, TERM, ACT); m = ACT.sum()
        mu = (adv * ACT).sum() / m; adv = (adv - mu) / (np.sqrt(((adv - mu) ** 2 * ACT).sum() / m) + 1e-8)
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
            print(f"  iter {it:>3}: normal-acc {_acc(net, 7000+it,'normal'):.2f}  frozen-acc {_acc(net, 7700+it,'frozen'):.2f}", flush=True)
    return net


def _acc(net, seed, mode, B=500):
    env = VecImit(B, seed, mode); rollout(net, env, greedy=True); return env.acc()


if __name__ == "__main__":
    a = sys.argv
    iters = int(a[a.index("--iters") + 1]) if "--iters" in a else 800
    here = __file__.rsplit("/", 1)[0]
    print(f"IMITATION probe: reach the HIDDEN correct goal by observing a demonstrator; {G} goals, chance {1/G:.2f}", flush=True)
    net = train(iters)
    nm = np.mean([_acc(net, 9000 + k, "normal") for k in range(3)])
    fr = np.mean([_acc(net, 9100 + k, "frozen") for k in range(3)])
    ms = np.mean([_acc(net, 9200 + k, "mislead") for k in range(3)])
    chance = 1 / G
    print(f"\nFINAL: normal {nm:.2f} | frozen-demo {fr:.2f} | misleading-demo {ms:.2f} (chance {chance:.2f})", flush=True)
    genuine = nm >= 0.7 and fr <= chance + 0.12          # strong when demonstrator present, ~chance when frozen
    print(f"VERDICT: {'GENUINE SOCIAL LEARNING (imitation/emulation)' if genuine else 'partial/none'} "
          f"(needs-demonstrator={nm-fr:+.2f}; misled-by-demo={nm-ms:+.2f})", flush=True)
    torch.save(net.state_dict(), f"{here}/imitation_net.pt")
    json.dump({"normal": round(float(nm), 2), "frozen": round(float(fr), 2), "mislead": round(float(ms), 2),
               "chance": round(float(chance), 2), "genuine": bool(genuine)}, open(f"{here}/imitation_result.json", "w"))
