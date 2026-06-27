"""Phase A: does a GRADIENT learner (recurrent PPO) solve the exact-counting task that evolution
(OpenAI-ES and SNES) could not (exact-k ~0.10-0.17)? The agent must collect EXACTLY k target objects
(supply > k, so it must STOP, not exhaust the supply) then reach the goal -- with NO oracle count
(it must internalize the tally in the GRU hidden state). If PPO solves it, the ceiling was the
LEARNING METHOD (evolution + local optima), and we now hold a capable learner for the off-map hunt.
If PPO also fails, the ceiling is deeper than expected.

Recurrent PPO: batched numpy env (B parallel), GRU policy+value, fixed-horizon rollouts with
done-masking, full-sequence BPTT, GAE, clipped surrogate, entropy bonus.

Run (on VPS epc-venv): python counting_ppo.py [--k 3] [--supply 5] [--iters 300]
"""
import numpy as np, sys, torch, torch.nn as nn

N = 9; OBS = 5; NACT = 5; MAXT = 60
DIRS = np.array([(0, 0), (1, 0), (-1, 0), (0, 1), (0, -1)])
torch.set_num_threads(8)


class VecCount:
    def __init__(self, B, k, supply, seed):
        self.B, self.k, self.S = B, k, supply
        self.rng = np.random.default_rng(seed)
        self.reset()

    def reset(self):
        B, S = self.B, self.S
        self.obj = np.zeros((B, S, 2), int); self.goal = np.zeros((B, 2), int); self.pos = np.zeros((B, 2), int)
        for b in range(B):
            idx = self.rng.permutation(N * N)[:S + 2]
            xs, ys = idx % N, idx // N
            self.obj[b, :, 0], self.obj[b, :, 1] = xs[:S], ys[:S]
            self.goal[b] = (xs[S], ys[S]); self.pos[b] = (xs[S + 1], ys[S + 1])
        self.alive = np.ones((B, S), bool); self.count = np.zeros(B, int)
        self.done = np.zeros(B, bool); self.just = np.zeros(B)
        return self.obs()

    def obs(self):
        B = self.B
        d = np.abs(self.obj - self.pos[:, None, :]).sum(2).astype(float)
        d[~self.alive] = 1e9
        nearest = self.obj[np.arange(B), d.argmin(1)]
        has = self.alive.any(1)
        dobj = np.sign(nearest - self.pos) * has[:, None]
        dgoal = np.sign(self.goal - self.pos)
        return np.concatenate([dobj, dgoal, self.just[:, None]], 1).astype(np.float32)

    def step(self, act):
        B = self.B
        active = ~self.done
        mv = DIRS[act]
        newp = self.pos + mv
        inb = (newp[:, 0] >= 0) & (newp[:, 0] < N) & (newp[:, 1] >= 0) & (newp[:, 1] < N)
        self.pos = np.where((inb & active)[:, None], newp, self.pos)
        rew = np.where(active, -0.02, 0.0)
        self.just = np.zeros(B)
        on = (self.obj == self.pos[:, None, :]).all(2) & self.alive
        for b in np.where(on.any(1) & active)[0]:
            j = np.where(on[b])[0][0]
            self.alive[b, j] = False; self.count[b] += 1; self.just[b] = 1.0
            rew[b] += 0.5 if self.count[b] <= self.k else -1.0
        atgoal = (self.pos == self.goal).all(1) & active
        term = np.zeros(B, bool)
        for b in np.where(atgoal)[0]:
            rew[b] += 2.0 if self.count[b] == self.k else -0.5
            self.done[b] = True; term[b] = True
        return self.obs(), rew.astype(np.float32), term, active


class Policy(nn.Module):
    def __init__(self, H=64):
        super().__init__()
        self.gru = nn.GRU(OBS, H, batch_first=True)
        self.pi = nn.Linear(H, NACT); self.v = nn.Linear(H, 1)

    def forward(self, x, h=None):
        out, h = self.gru(x, h)
        return self.pi(out), self.v(out).squeeze(-1), h


def rollout(net, env, k, supply, seed):
    env.__init__(env.B, k, supply, seed)
    B = env.B
    obs = env.reset()
    O = np.zeros((B, MAXT, OBS), np.float32); A = np.zeros((B, MAXT), int)
    LP = np.zeros((B, MAXT), np.float32); V = np.zeros((B, MAXT), np.float32)
    R = np.zeros((B, MAXT), np.float32); TERM = np.zeros((B, MAXT), np.float32); ACT = np.zeros((B, MAXT), np.float32)
    h = None
    for t in range(MAXT):
        ot = torch.from_numpy(obs)[:, None, :]
        with torch.no_grad():
            logit, val, h = net(ot, h)
        logit = logit[:, 0]; dist = torch.distributions.Categorical(logits=logit)
        a = dist.sample()
        O[:, t] = obs; A[:, t] = a.numpy(); LP[:, t] = dist.log_prob(a).numpy(); V[:, t] = val[:, 0].numpy()
        active = ~env.done.copy()
        obs, r, term, _ = env.step(a.numpy())
        R[:, t] = r; TERM[:, t] = term; ACT[:, t] = active
    return O, A, LP, V, R, TERM, ACT


def gae(R, V, TERM, ACT, gamma=0.99, lam=0.95):
    B, T = R.shape
    adv = np.zeros((B, T), np.float32); last = np.zeros(B, np.float32)
    for t in reversed(range(T)):
        nextv = V[:, t + 1] if t + 1 < T else np.zeros(B, np.float32)
        nonterm = 1.0 - TERM[:, t]
        delta = R[:, t] + gamma * nextv * nonterm - V[:, t]
        last = delta + gamma * lam * nonterm * last
        adv[:, t] = last * ACT[:, t]
    return adv, adv + V


def evaluate(net, k, supply, seed, B=400):
    env = VecCount(B, k, supply, seed); obs = env.reset(); h = None
    for t in range(MAXT):
        with torch.no_grad():
            logit, _, h = net(torch.from_numpy(obs)[:, None, :], h)
        a = logit[:, 0].argmax(1).numpy()
        obs, r, term, _ = env.step(a)
    exact = ((env.done) & (env.count == k)).mean()
    return float(exact), float(env.count.mean())


if __name__ == "__main__":
    a = sys.argv
    k = int(a[a.index("--k") + 1]) if "--k" in a else 3
    supply = int(a[a.index("--supply") + 1]) if "--supply" in a else 5
    iters = int(a[a.index("--iters") + 1]) if "--iters" in a else 300
    B = 256
    net = Policy(64); opt = torch.optim.Adam(net.parameters(), lr=3e-3)
    env = VecCount(B, k, supply, 0)
    print(f"recurrent PPO on exact-counting: collect EXACTLY k={k} (supply={supply}), then goal. iters={iters}")
    for it in range(iters):
        O, A, LP, V, R, TERM, ACT = rollout(net, env, k, supply, 1000 + it)
        adv, ret = gae(R, V, TERM, ACT)
        m = ACT.sum()
        adv = (adv - (adv * ACT).sum() / m) / (np.sqrt(((adv - (adv * ACT).sum() / m) ** 2 * ACT).sum() / m) + 1e-8)
        Ot = torch.from_numpy(O); At = torch.from_numpy(A); LPt = torch.from_numpy(LP)
        advt = torch.from_numpy(adv); rett = torch.from_numpy(ret); actt = torch.from_numpy(ACT)
        for _ in range(4):
            logit, val, _ = net(Ot)
            dist = torch.distributions.Categorical(logits=logit)
            nlp = dist.log_prob(At)
            ratio = torch.exp(nlp - LPt)
            s1 = ratio * advt; s2 = torch.clamp(ratio, 0.8, 1.2) * advt
            pol = -(torch.min(s1, s2) * actt).sum() / actt.sum()
            vf = (((val - rett) ** 2) * actt).sum() / actt.sum()
            ent = (dist.entropy() * actt).sum() / actt.sum()
            loss = pol + 0.5 * vf - 0.01 * ent
            opt.zero_grad(); loss.backward(); nn.utils.clip_grad_norm_(net.parameters(), 0.5); opt.step()
        if it % 30 == 0 or it == iters - 1:
            ex, mc = evaluate(net, k, supply, 7)
            print(f"  iter {it:>3}: exact-k {ex:.2f}  mean count {mc:.2f} (target {k})")
    ex, mc = evaluate(net, k, supply, 77)
    ex2, mc2 = evaluate(net, k, supply * 2, 78)   # does it STOP at k under doubled supply? (the counting test)
    print(f"\nFINAL: exact-k {ex:.2f}  mean count {mc:.2f} (target {k})")
    print(f"  doubled-supply: mean count {mc2:.2f} exact-k {ex2:.2f}  (a real counter stops at {k} regardless of supply)")
    print("VERDICT:", "PPO COUNTS (exact-k>=0.6 and stops at k under doubled supply) -> learning method WAS the ceiling; capable learner in hand"
          if (ex >= 0.6 and mc2 < k + 1.0) else f"PPO did NOT cleanly count (exact-k {ex:.2f}) -> ceiling deeper than the learner")
    torch.save(net.state_dict(), __file__.rsplit("/", 1)[0] + "/counting_ppo_net.pt")
