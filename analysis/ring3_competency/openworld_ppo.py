"""Phase B: the open-ended off-map hunt with a CAPABLE learner. Reuses the open-ended generator
(reward rule sampled from openworld's count/exact/order/collect grammar -- nothing baked in) but trains
recurrent PPO (validated in Phase A) instead of the ceilinged evolutionary search. For each sampled
environment: train a GRU-PPO agent + a memoryless PPO baseline over domain-randomized layouts; the MCC
filter keeps environments PPO solves but the memoryless baseline cannot (i.e. that genuinely demand
internal state). Demanding environments + their agents go to the agent-observer (Layer 2) to NAME the
competency and classify it vs the map {navigation, memory, delayed-grat, regulation, +combinations,
+counting/accumulation} or flag OFF-MAP. Honest bound held: literature gate before any novelty claim.

Run on VPS epc-venv: python openworld_ppo.py [--envs 8] [--iters 200] [--seed 11] [--mode open]
"""
import numpy as np, sys, json, os, torch, torch.nn as nn
from openworld import sample_predicate, descr

N = 9; T = 3; PER = 5; MAXT = 70; H = 64
OBS = T * 2 + 2 + T; NACT = 5
DIRS = np.array([(0, 0), (1, 0), (-1, 0), (0, 1), (0, -1)])
torch.set_num_threads(8)


def pred_frac_vec(lits, counts, first):
    """counts (B,T) int, first (B,T) int (-1 = untouched). Returns (B,) fraction of literals satisfied."""
    B = counts.shape[0]; ok = np.zeros(B)
    for lit in lits:
        if lit[0] == "count":
            ok += counts[:, lit[1]] >= lit[2]
        elif lit[0] == "exact":
            ok += counts[:, lit[1]] == lit[2]
        elif lit[0] == "order":
            a, b = lit[1], lit[2]
            ok += (first[:, a] >= 0) & (first[:, b] >= 0) & (first[:, a] < first[:, b])
        else:  # collect
            sub = lit[1]; good = np.ones(B, bool)
            for t in sub: good &= counts[:, t] >= 1
            ok += good
    return ok / len(lits)


class VecOpen:
    def __init__(self, B, lits, seed):
        self.B, self.lits = B, lits
        self.rng = np.random.default_rng(seed)
        self.S = T * PER
        self.otype = np.tile(np.repeat(np.arange(T), PER), (B, 1))   # (B,S)
        self.reset()

    def reset(self):
        B, S = self.B, self.S
        self.obj = np.zeros((B, S, 2), int); self.goal = np.zeros((B, 2), int); self.pos = np.zeros((B, 2), int)
        for b in range(B):
            idx = self.rng.permutation(N * N)[:S + 2]
            self.obj[b, :, 0], self.obj[b, :, 1] = idx[:S] % N, idx[:S] // N
            self.goal[b] = (idx[S] % N, idx[S] // N); self.pos[b] = (idx[S + 1] % N, idx[S + 1] // N)
        self.alive = np.ones((B, S), bool); self.counts = np.zeros((B, T), int); self.first = -np.ones((B, T), int)
        self.done = np.zeros(B, bool); self.just = -np.ones(B, int); self.step_i = 0
        self.frac = pred_frac_vec(self.lits, self.counts, self.first)
        return self.obs()

    def obs(self):
        B = self.B; out = np.zeros((B, OBS), np.float32); o = 0
        for t in range(T):
            mask = (self.otype == t) & self.alive
            d = np.abs(self.obj - self.pos[:, None, :]).sum(2).astype(float)
            d[~mask] = 1e9
            j = d.argmin(1); nearest = self.obj[np.arange(B), j]; has = mask.any(1)
            out[:, o:o + 2] = np.sign(nearest - self.pos) * has[:, None]; o += 2
        out[:, o:o + 2] = np.sign(self.goal - self.pos); o += 2
        for b in range(B):
            if self.just[b] >= 0: out[b, o + self.just[b]] = 1.0
        return out

    def step(self, act):
        B = self.B; active = ~self.done
        newp = self.pos + DIRS[act]
        inb = (newp[:, 0] >= 0) & (newp[:, 0] < N) & (newp[:, 1] >= 0) & (newp[:, 1] < N)
        self.pos = np.where((inb & active)[:, None], newp, self.pos)
        self.just = -np.ones(B, int)
        on = (self.obj == self.pos[:, None, :]).all(2) & self.alive
        for b in np.where(on.any(1) & active)[0]:
            j = np.where(on[b])[0][0]; t = self.otype[b, j]
            self.alive[b, j] = False; self.counts[b, t] += 1
            if self.first[b, t] < 0: self.first[b, t] = self.step_i
            self.just[b] = t
        newfrac = pred_frac_vec(self.lits, self.counts, self.first)
        rew = np.where(active, -0.01 + 2.0 * (newfrac - self.frac), 0.0)
        self.frac = newfrac
        atgoal = (self.pos == self.goal).all(1) & active
        term = np.zeros(B, bool)
        for b in np.where(atgoal)[0]:
            rew[b] += 3.0 if self.frac[b] >= 0.999 else -0.3
            self.done[b] = True; term[b] = True
        self.step_i += 1
        return self.obs(), rew.astype(np.float32), term, active


class Policy(nn.Module):
    def __init__(self):
        super().__init__()
        self.gru = nn.GRU(OBS, H, batch_first=True)
        self.pi = nn.Linear(H, NACT); self.v = nn.Linear(H, 1)

    def forward(self, x, h=None):
        out, h = self.gru(x, h)
        return self.pi(out), self.v(out).squeeze(-1), h


class MLPPolicy(nn.Module):
    """TRUE memoryless baseline: feedforward, per-timestep, no recurrence -> rollout and update are
    trivially consistent (the previous GRU-with-recurrent=False baseline was broken: memoryless
    rollout but recurrent update)."""
    def __init__(self):
        super().__init__()
        self.body = nn.Sequential(nn.Linear(OBS, H), nn.Tanh(), nn.Linear(H, H), nn.Tanh())
        self.pi = nn.Linear(H, NACT); self.v = nn.Linear(H, 1)

    def forward(self, x, h=None):
        z = self.body(x)
        return self.pi(z), self.v(z).squeeze(-1), None


def run_episode_batch(net, env):
    B = env.B; obs = env.reset()
    O = np.zeros((B, MAXT, OBS), np.float32); A = np.zeros((B, MAXT), int)
    LP = np.zeros((B, MAXT), np.float32); V = np.zeros((B, MAXT), np.float32)
    R = np.zeros((B, MAXT), np.float32); TERM = np.zeros((B, MAXT), np.float32); ACT = np.zeros((B, MAXT), np.float32)
    h = None
    for t in range(MAXT):
        ot = torch.from_numpy(obs)[:, None, :]
        with torch.no_grad():
            logit, val, h = net(ot, h)
        dist = torch.distributions.Categorical(logits=logit[:, 0]); a = dist.sample()
        O[:, t] = obs; A[:, t] = a.numpy(); LP[:, t] = dist.log_prob(a).numpy(); V[:, t] = val[:, 0].numpy()
        ACT[:, t] = (~env.done).astype(np.float32)
        obs, r, term, _ = env.step(a.numpy()); R[:, t] = r; TERM[:, t] = term
    return O, A, LP, V, R, TERM, ACT


def gae(R, V, TERM, ACT, gamma=0.99, lam=0.95):
    B, Tt = R.shape; adv = np.zeros((B, Tt), np.float32); last = np.zeros(B, np.float32)
    for t in reversed(range(Tt)):
        nextv = V[:, t + 1] if t + 1 < Tt else np.zeros(B, np.float32)
        nonterm = 1.0 - TERM[:, t]
        delta = R[:, t] + gamma * nextv * nonterm - V[:, t]
        last = delta + gamma * lam * nonterm * last
        adv[:, t] = last * ACT[:, t]
    return adv, adv + V


def success_rate(net, lits, seed, B=400):
    env = VecOpen(B, lits, seed); obs = env.reset(); h = None
    for t in range(MAXT):
        with torch.no_grad():
            logit, _, h = net(torch.from_numpy(obs)[:, None, :], h)
        a = logit[:, 0].argmax(1).numpy(); obs, r, term, _ = env.step(a)
    return float(((env.done) & (env.frac >= 0.999)).mean())


def train(lits, net_cls, iters, seed, B=256):
    torch.manual_seed(seed); net = net_cls(); opt = torch.optim.Adam(net.parameters(), lr=3e-3)
    env = VecOpen(B, lits, seed)
    for it in range(iters):
        env.__init__(B, lits, 1000 + seed * 9999 + it)
        O, A, LP, V, R, TERM, ACT = run_episode_batch(net, env)
        adv, ret = gae(R, V, TERM, ACT); m = ACT.sum()
        mean = (adv * ACT).sum() / m; std = np.sqrt(((adv - mean) ** 2 * ACT).sum() / m) + 1e-8
        adv = (adv - mean) / std
        Ot, At, LPt = torch.from_numpy(O), torch.from_numpy(A), torch.from_numpy(LP)
        advt, rett, actt = torch.from_numpy(adv), torch.from_numpy(ret), torch.from_numpy(ACT)
        for _ in range(4):
            logit, val, _ = net(Ot)
            dist = torch.distributions.Categorical(logits=logit)
            ratio = torch.exp(dist.log_prob(At) - LPt)
            s1 = ratio * advt; s2 = torch.clamp(ratio, 0.8, 1.2) * advt
            loss = -(torch.min(s1, s2) * actt).sum() / actt.sum() \
                + 0.5 * (((val - rett) ** 2) * actt).sum() / actt.sum() \
                - 0.01 * (dist.entropy() * actt).sum() / actt.sum()
            opt.zero_grad(); loss.backward(); nn.utils.clip_grad_norm_(net.parameters(), 0.5); opt.step()
    return net


if __name__ == "__main__":
    a = sys.argv
    nenv = int(a[a.index("--envs") + 1]) if "--envs" in a else 8
    iters = int(a[a.index("--iters") + 1]) if "--iters" in a else 200
    mode = a[a.index("--mode") + 1] if "--mode" in a else "open"
    rmeta = np.random.default_rng(int(a[a.index("--seed") + 1]) if "--seed" in a else 11)
    here = __file__.rsplit("/", 1)[0]; out = f"{here}/openworld_ppo_sweep.json"
    # RESUMABLE: reload completed envs so a crash/reboot loses <=1 env; checkpoint after each env.
    prior = {r["env"]: r for r in json.load(open(out))} if os.path.exists(out) else {}
    results = []
    print(f"PHASE B open-ended PPO hunt: {nenv} envs, iters={iters}, mode={mode}, PER={PER} "
          f"(resuming: {sorted(prior)} cached)", flush=True)
    for e in range(nenv):
        lits = sample_predicate(rmeta, mode=mode)              # advance rng every env to keep predicates aligned
        if e in prior:
            results.append(prior[e]); print(f"env {e}: (cached) {prior[e]['pred']}", flush=True); continue
        ppo = train(lits, Policy, iters, 10 + e)
        mem = train(lits, MLPPolicy, iters, 100 + e)
        ps = success_rate(ppo, lits, 500 + e); ms = success_rate(mem, lits, 600 + e)
        dem = ps >= 0.6 and (ps - ms) >= 0.25
        if dem: torch.save(ppo.state_dict(), f"{here}/openworld_ppo_env{e}.pt")
        results.append({"env": e, "pred": descr(lits), "lits": lits, "ppo": round(ps, 2), "mem": round(ms, 2),
                        "gap": round(ps - ms, 2), "demanding": dem})
        json.dump(results, open(out, "w"), indent=1)           # checkpoint after EACH env
        print(f"env {e}: {descr(lits):44s} PPO {ps:.2f}  memoryless {ms:.2f}  gap {ps - ms:+.2f}  {'DEMANDING' if dem else ''}", flush=True)
    dem = [x for x in results if x["demanding"]]
    print(f"\n{len(dem)}/{nenv} environments DEMAND competency (PPO>=0.6 AND gap>=0.25): {[x['env'] for x in dem]}", flush=True)
    json.dump(results, open(out, "w"), indent=1)
    print("saved openworld_ppo_sweep.json", flush=True)
