"""Phase D3: symmetry-breaking under contest -- the probe built to DEFEAT the shortcuts that kept D1/D2
landing on known competencies. Two ANONYMOUS agents (no agent-id in the observation -> they cannot
pre-split roles by identity) share a policy. One goal is HIGH-value (+2), the others LOW (+0.5). Team
reward is shared, so the optimum is to SPLIT -- one agent takes high, the other a low -- but if BOTH go
high they COLLIDE (-1 each). With anonymity + a shared policy + identical obs structure, the agents must
BREAK SYMMETRY to decide who takes the contested goal: either via the channel (negotiate/claim) or via an
emergent convention on relative geometry. Neither a shared deterministic rule (no id) nor pure passive
observation (symmetric) trivially solves it.

Metrics: capture = fraction of episodes where the high goal is taken by EXACTLY one agent (contest
resolved); collision = fraction where BOTH pile onto high (resolution failed). Discriminators:
channel-scramble (collapse => resolution uses COMMUNICATION) and blind-partner (collapse => uses
relative geometry / the partner). The agent-observer then NAMES the symmetry-breaking mechanism and asks
whether the literature already has a clean name for it.

Run on VPS epc-venv: python contest_ppo.py [--iters 800]
"""
import numpy as np, sys, json, os, torch, torch.nn as nn

N = 5; G = 3; K = 4; MAXT = 18; H = 64
ODIM = G * 3 + 2 + K              # per-goal: dir(2)+is_high(1); + other rel-dir(2) + other last-symbol(K). NO id.
DIRS = np.array([(0, 0), (1, 0), (-1, 0), (0, 1), (0, -1)])
torch.set_num_threads(8)


class Agent(nn.Module):
    def __init__(self):
        super().__init__()
        self.net = nn.Sequential(nn.Linear(ODIM, H), nn.Tanh(), nn.Linear(H, H), nn.Tanh())
        self.move = nn.Linear(H, 5); self.sym = nn.Linear(H, K); self.v = nn.Linear(H, 1)

    def forward(self, x):
        z = self.net(x); return self.move(z), self.sym(z), self.v(z).squeeze(-1)


class VecContest:
    def __init__(self, B, seed):
        self.B = B; self.rng = np.random.default_rng(seed); self.reset()

    def reset(self):
        B = self.B
        self.goal = np.zeros((B, G, 2), int); self.p = np.zeros((B, 2, 2), int)
        for b in range(B):
            idx = self.rng.permutation(N * N)[:G + 2]
            self.goal[b, :, 0], self.goal[b, :, 1] = idx[:G] % N, idx[:G] // N
            self.p[b, 0] = (idx[G] % N, idx[G] // N); self.p[b, 1] = (idx[G + 1] % N, idx[G + 1] // N)
        self.high = self.rng.integers(0, G, size=B)            # which goal is high-value (random index)
        self.last_sym = np.zeros((B, 2), int)
        return None

    def on_goal(self, i):
        eq = (self.goal == self.p[:, i, None, :]).all(2)
        return np.where(eq.any(1), eq.argmax(1), -1)

    def obs(self, i, scramble=False, blind=False):
        B = self.B; out = np.zeros((B, ODIM), np.float32); o = 0
        for g in range(G):
            out[:, o:o + 2] = np.sign(self.goal[:, g, :] - self.p[:, i, :]); o += 2
            out[:, o] = (self.high == g).astype(np.float32); o += 1   # is this goal the high-value one
        oth = 1 - i
        out[:, o:o + 2] = 0.0 if blind else np.sign(self.p[:, oth, :] - self.p[:, i, :]); o += 2
        sym = self.rng.integers(0, K, size=B) if scramble else self.last_sym[:, oth]
        if not blind:
            out[np.arange(B), o + sym] = 1.0
        o += K
        return out

    def step(self, mv):
        B = self.B
        for i in range(2):
            np_ = self.p[:, i, :] + DIRS[mv[:, i]]
            inb = (np_[:, 0] >= 0) & (np_[:, 0] < N) & (np_[:, 1] >= 0) & (np_[:, 1] < N)
            self.p[:, i, :] = np.where(inb[:, None], np_, self.p[:, i, :])
        r = np.full((B, 2), -0.01, np.float32)
        for i in range(2):
            r[:, i] += 0.02 * (self.on_goal(i) >= 0)
        return r

    def terminal(self):                                        # shared team reward
        B = self.B; g0, g1 = self.on_goal(0), self.on_goal(1); team = np.zeros(B, np.float32)
        both_high = (g0 == self.high) & (g1 == self.high)
        one_high = ((g0 == self.high) ^ (g1 == self.high))
        team += 2.0 * one_high                                  # high captured by exactly one
        team += -2.0 * both_high                                # collision (both -1)
        # low-goal credit: an agent alone-ish on a low goal
        for i, gi in enumerate((g0, g1)):
            on_low = (gi >= 0) & (gi != self.high)
            team += 0.5 * on_low
        return np.stack([team, team], 1)

    def stats(self):
        g0, g1 = self.on_goal(0), self.on_goal(1)
        capture = float((((g0 == self.high) ^ (g1 == self.high))).mean())
        collision = float(((g0 == self.high) & (g1 == self.high)).mean())
        return capture, collision


def gae(R, V, gamma=0.99, lam=0.95):
    A, Tt = R.shape; adv = np.zeros((A, Tt), np.float32); last = np.zeros(A, np.float32)
    for t in reversed(range(Tt)):
        nextv = V[:, t + 1] if t + 1 < Tt else np.zeros(A, np.float32)
        delta = R[:, t] + gamma * nextv - V[:, t]
        last = delta + gamma * lam * last; adv[:, t] = last
    return adv, adv + V


def rollout(net, env, greedy=False, scramble=False, blind=False):
    B = env.B; env.reset()
    O = np.zeros((B, 2, MAXT, ODIM), np.float32); MV = np.zeros((B, 2, MAXT), int); SY = np.zeros((B, 2, MAXT), int)
    LP = np.zeros((B, 2, MAXT), np.float32); V = np.zeros((B, 2, MAXT), np.float32); R = np.zeros((B, 2, MAXT), np.float32)
    for t in range(MAXT):
        for i in range(2):
            ob = env.obs(i, scramble=scramble, blind=blind); O[:, i, t] = ob
            with torch.no_grad():
                ml, sl, v = net(torch.from_numpy(ob))
            md, sd = torch.distributions.Categorical(logits=ml), torch.distributions.Categorical(logits=sl)
            mv = ml.argmax(1) if greedy else md.sample(); sy = sl.argmax(1) if greedy else sd.sample()
            MV[:, i, t] = mv.numpy(); SY[:, i, t] = sy.numpy()
            LP[:, i, t] = (md.log_prob(mv) + sd.log_prob(sy)).numpy(); V[:, i, t] = v.numpy()
        r = env.step(np.stack([MV[:, 0, t], MV[:, 1, t]], 1))
        env.last_sym = np.stack([SY[:, 0, t], SY[:, 1, t]], 1); R[:, :, t] = r
    R[:, :, -1] += env.terminal()
    return O, MV, SY, LP, V, R


def train(iters, seed=0, B=320):
    torch.manual_seed(seed); net = Agent(); opt = torch.optim.Adam(net.parameters(), lr=3e-3)
    for it in range(iters):
        env = VecContest(B, 1000 + seed * 99999 + it)
        O, MV, SY, LP, V, R = rollout(net, env)
        O2 = O.reshape(2 * B, MAXT, ODIM); MV2, SY2, LP2 = MV.reshape(2 * B, MAXT), SY.reshape(2 * B, MAXT), LP.reshape(2 * B, MAXT)
        V2, R2 = V.reshape(2 * B, MAXT), R.reshape(2 * B, MAXT)
        adv, ret = gae(R2, V2); adv = (adv - adv.mean()) / (adv.std() + 1e-8)
        Ot, MVt, SYt, LPt = torch.from_numpy(O2), torch.from_numpy(MV2), torch.from_numpy(SY2), torch.from_numpy(LP2)
        advt, rett = torch.from_numpy(adv), torch.from_numpy(ret)
        for _ in range(4):
            ml, sl, v = net(Ot.reshape(-1, ODIM))
            ml = ml.reshape(2 * B, MAXT, 5); sl = sl.reshape(2 * B, MAXT, K); v = v.reshape(2 * B, MAXT)
            md, sd = torch.distributions.Categorical(logits=ml), torch.distributions.Categorical(logits=sl)
            ratio = torch.exp(md.log_prob(MVt) + sd.log_prob(SYt) - LPt)
            s1 = ratio * advt; s2 = torch.clamp(ratio, 0.8, 1.2) * advt
            loss = -torch.min(s1, s2).mean() + 0.5 * ((v - rett) ** 2).mean() - 0.01 * (md.entropy() + sd.entropy()).mean()
            opt.zero_grad(); loss.backward(); nn.utils.clip_grad_norm_(net.parameters(), 0.5); opt.step()
        if it % 100 == 0 or it == iters - 1:
            cap, col = _eval(net, 7000 + it)
            print(f"  iter {it:>3}: capture {cap:.2f}  collision {col:.2f}", flush=True)
    return net


def _eval(net, seed, scramble=False, blind=False, B=500):
    env = VecContest(B, seed); rollout(net, env, greedy=True, scramble=scramble, blind=blind); return env.stats()


if __name__ == "__main__":
    a = sys.argv
    iters = int(a[a.index("--iters") + 1]) if "--iters" in a else 800
    here = __file__.rsplit("/", 1)[0]
    print(f"PHASE D3 symmetry-breaking under contest: anonymous agents, 1 high(+2)/2 low(+0.5) goals", flush=True)
    net = train(iters)
    cap, col = np.mean([_eval(net, 9000 + k) for k in range(3)], 0)
    capS, colS = np.mean([_eval(net, 9100 + k, scramble=True) for k in range(3)], 0)
    capB, colB = np.mean([_eval(net, 9200 + k, blind=True) for k in range(3)], 0)
    print(f"\nFINAL  capture/collision: normal {cap:.2f}/{col:.2f} | chan-scramble {capS:.2f}/{colS:.2f} | blind-partner {capB:.2f}/{colB:.2f}", flush=True)
    print(f"  (random-greedy floor would be ~collision-heavy; genuine symmetry-breaking = high capture, low collision)", flush=True)
    uses_comms = (cap - capS) >= 0.2; uses_partner = (cap - capB) >= 0.2
    print(f"VERDICT: symmetry-breaking {'ACHIEVED' if cap >= 0.7 and col <= 0.2 else 'partial/none'} "
          f"(capture {cap:.2f}, collision {col:.2f}); mechanism: comms={uses_comms} partner-geometry={uses_partner}", flush=True)
    torch.save(net.state_dict(), f"{here}/contest_net.pt")
    json.dump({"capture": round(float(cap), 3), "collision": round(float(col), 3),
               "scramble_capture": round(float(capS), 3), "blind_capture": round(float(capB), 3),
               "uses_comms": bool(uses_comms), "uses_partner": bool(uses_partner)}, open(f"{here}/contest_result.json", "w"))
