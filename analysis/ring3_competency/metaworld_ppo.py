"""Phase C: the genuine off-map probe. Enriches the generator past obviously-known object-level
competencies with rules whose competency is a RUNTIME / META capability:
  - cued   : a start-cue (seen only at t=0) names the good type -> conditional rule-selection + memory
  - switch : the good type flips after each good collection -> track-and-adapt
  - infer  : the good type is HIDDEN; the agent must INFER it from reward feedback within the episode,
             then exploit it -> in-episode rule inference (meta-learning / RL^2). The strongest off-map
             candidate: the competency is not "execute rule X" but "figure out the rule, then exploit it".

Each episode: collecting a good-type object = +1, bad = -1, step = -0.02; fixed horizon. The agent
observes direction-to-nearest per type, last-touched type, and LAST-REWARD-SIGN (the feedback channel
that makes inference possible), plus a cue channel (only populated for `cued`). A memoryless MLP control
cannot infer/track/remember -> ~chance precision; a GRU can. MCC keeps rules the GRU solves but the
memoryless baseline cannot. Demanding rules go to the agent-observer (name blind -> debunk -> classify
vs the map {commitment, storage, accumulation, reactive...} or flag OFF-MAP). Honest bound: literature
gate before any novelty claim (meta-learning is KNOWN -> Tier-2 at best, but the richest face yet).

Run on VPS epc-venv: python metaworld_ppo.py [--iters 300] [--rules cued,infer,switch]
"""
import numpy as np, sys, json, os, torch, torch.nn as nn

N = 9; T = 3; PER = 4; MAXT = 45; H = 64
OBS = T * 2 + T + 1 + T          # dirs/type, last-touched onehot, last-reward-sign, cue
NACT = 5
DIRS = np.array([(0, 0), (1, 0), (-1, 0), (0, 1), (0, -1)])
torch.set_num_threads(8)


class VecMeta:
    def __init__(self, B, rule, seed, ablate_fb=False, flip_at=None):
        self.B, self.rule = B, rule
        self.ablate_fb = ablate_fb      # zero the reward-feedback channel (debunk: kills genuine inference)
        self.flip_at = flip_at          # resample the good type at this step (debunk: tests re-adaptation)
        self.rng = np.random.default_rng(seed); self.S = T * PER
        self.otype = np.tile(np.repeat(np.arange(T), PER), (B, 1))
        self.reset()

    def reset(self):
        B, S = self.B, self.S
        self.obj = np.zeros((B, S, 2), int); self.pos = np.zeros((B, 2), int)
        for b in range(B):
            idx = self.rng.permutation(N * N)[:S + 1]
            self.obj[b, :, 0], self.obj[b, :, 1] = idx[:S] % N, idx[:S] // N
            self.pos[b] = (idx[S] % N, idx[S] // N)
        self.alive = np.ones((B, S), bool)
        self.g = self.rng.integers(0, T, size=B)                # hidden good type per env
        self.last_t = -np.ones(B, int); self.last_r = np.zeros(B)
        self.good = np.zeros(B, int); self.bad = np.zeros(B, int); self.step_i = 0
        return self.obs()

    def obs(self):
        B = self.B; out = np.zeros((B, OBS), np.float32); o = 0
        for t in range(T):
            m = (self.otype == t) & self.alive
            d = np.abs(self.obj - self.pos[:, None, :]).sum(2).astype(float); d[~m] = 1e9
            nr = self.obj[np.arange(B), d.argmin(1)]; has = m.any(1)
            out[:, o:o + 2] = np.sign(nr - self.pos) * has[:, None]; o += 2
        for b in range(B):
            if self.last_t[b] >= 0: out[b, o + self.last_t[b]] = 1.0
        o += T
        out[:, o] = 0.0 if self.ablate_fb else self.last_r; o += 1
        if self.rule == "cued" and self.step_i == 0:
            for b in range(B): out[b, o + self.g[b]] = 1.0      # cue: good type, t=0 only
        o += T
        return out

    def step(self, act):
        B = self.B
        if self.flip_at is not None and self.step_i == self.flip_at:
            self.g = (self.g + 1) % T                           # flip the good type mid-episode
        newp = self.pos + DIRS[act]
        inb = (newp[:, 0] >= 0) & (newp[:, 0] < N) & (newp[:, 1] >= 0) & (newp[:, 1] < N)
        self.pos = np.where(inb[:, None], newp, self.pos)
        rew = np.full(B, -0.01); self.last_t = -np.ones(B, int); self.last_r = np.zeros(B)
        on = (self.obj == self.pos[:, None, :]).all(2) & self.alive
        for b in np.where(on.any(1))[0]:
            j = np.where(on[b])[0][0]; t = self.otype[b, j]; self.alive[b, j] = False
            good = (t == self.g[b])
            # bad penalty kept SMALL (-0.3) so exploring-then-exploiting beats abstaining; the
            # feedback SIGN stays full (+1/-1) so inference is still possible from last_r.
            rew[b] += 1.0 if good else -0.3
            self.last_t[b] = t; self.last_r[b] = 1.0 if good else -1.0
            if good: self.good[b] += 1
            else: self.bad[b] += 1
            if self.rule == "switch" and good: self.g[b] = (self.g[b] + 1) % T
        self.step_i += 1
        return self.obs(), rew.astype(np.float32), np.zeros(B, bool), np.ones(B, bool)


class Policy(nn.Module):
    def __init__(self):
        super().__init__(); self.gru = nn.GRU(OBS, H, batch_first=True)
        self.pi = nn.Linear(H, NACT); self.v = nn.Linear(H, 1)

    def forward(self, x, h=None):
        out, h = self.gru(x, h); return self.pi(out), self.v(out).squeeze(-1), h


class MLPPolicy(nn.Module):                                     # true memoryless control
    def __init__(self):
        super().__init__()
        self.body = nn.Sequential(nn.Linear(OBS, H), nn.Tanh(), nn.Linear(H, H), nn.Tanh())
        self.pi = nn.Linear(H, NACT); self.v = nn.Linear(H, 1)

    def forward(self, x, h=None):
        z = self.body(x); return self.pi(z), self.v(z).squeeze(-1), None


def rollout(net, env):
    B = env.B; obs = env.reset()
    O = np.zeros((B, MAXT, OBS), np.float32); A = np.zeros((B, MAXT), int)
    LP = np.zeros((B, MAXT), np.float32); V = np.zeros((B, MAXT), np.float32); R = np.zeros((B, MAXT), np.float32)
    h = None
    for t in range(MAXT):
        with torch.no_grad():
            logit, val, h = net(torch.from_numpy(obs)[:, None, :], h)
        dist = torch.distributions.Categorical(logits=logit[:, 0]); a = dist.sample()
        O[:, t] = obs; A[:, t] = a.numpy(); LP[:, t] = dist.log_prob(a).numpy(); V[:, t] = val[:, 0].numpy()
        obs, r, _, _ = env.step(a.numpy()); R[:, t] = r
    return O, A, LP, V, R


def gae(R, V, gamma=0.99, lam=0.95):
    B, Tt = R.shape; adv = np.zeros((B, Tt), np.float32); last = np.zeros(B, np.float32)
    for t in reversed(range(Tt)):
        nextv = V[:, t + 1] if t + 1 < Tt else np.zeros(B, np.float32)
        delta = R[:, t] + gamma * nextv - V[:, t]
        last = delta + gamma * lam * last; adv[:, t] = last
    return adv, adv + V


def evaluate(net, rule, seed, B=400):
    env = VecMeta(B, rule, seed); obs = env.reset(); h = None; net_r = np.zeros(B)
    for t in range(MAXT):
        with torch.no_grad():
            logit, _, h = net(torch.from_numpy(obs)[:, None, :], h)
        a = logit[:, 0].argmax(1).numpy(); obs, r, _, _ = env.step(a); net_r += r
    coll = env.good + env.bad
    prec = float((env.good[coll > 0] / coll[coll > 0]).mean()) if (coll > 0).any() else 0.0
    return prec, float(net_r.mean())


def train(rule, net_cls, iters, seed, B=256, progress=False):
    torch.manual_seed(seed); net = net_cls(); opt = torch.optim.Adam(net.parameters(), lr=3e-3)
    for it in range(iters):
        if progress and (it % 150 == 0 or it == iters - 1):
            pp, pr = evaluate(net, rule, 777)
            print(f"    [{rule}] iter {it:>3}: prec {pp:.2f} netR {pr:+.1f}", flush=True)
        env = VecMeta(B, rule, 1000 + seed * 9999 + it)
        O, A, LP, V, R = rollout(net, env)
        adv, ret = gae(R, V); adv = (adv - adv.mean()) / (adv.std() + 1e-8)
        Ot, At, LPt = torch.from_numpy(O), torch.from_numpy(A), torch.from_numpy(LP)
        advt, rett = torch.from_numpy(adv), torch.from_numpy(ret)
        for _ in range(4):
            logit, val, _ = net(Ot); dist = torch.distributions.Categorical(logits=logit)
            ratio = torch.exp(dist.log_prob(At) - LPt)
            s1 = ratio * advt; s2 = torch.clamp(ratio, 0.8, 1.2) * advt
            loss = -torch.min(s1, s2).mean() + 0.5 * ((val - rett) ** 2).mean() - 0.01 * dist.entropy().mean()
            opt.zero_grad(); loss.backward(); nn.utils.clip_grad_norm_(net.parameters(), 0.5); opt.step()
    return net


if __name__ == "__main__":
    a = sys.argv
    iters = int(a[a.index("--iters") + 1]) if "--iters" in a else 300
    rules = (a[a.index("--rules") + 1].split(",")) if "--rules" in a else ["cued", "infer", "switch"]
    here = __file__.rsplit("/", 1)[0]; out = f"{here}/metaworld_sweep.json"
    prior = {r["rule"]: r for r in json.load(open(out))} if os.path.exists(out) else {}
    results = []
    print(f"PHASE C meta/conditional hunt: rules={rules}, iters={iters} (resuming: {sorted(prior)})", flush=True)
    for i, rule in enumerate(rules):
        if rule in prior:
            results.append(prior[rule]); print(f"{rule}: (cached)", flush=True); continue
        ppo = train(rule, Policy, iters, 10 + i, progress=True)
        mem = train(rule, MLPPolicy, iters, 200 + i)
        pp, pr = evaluate(ppo, rule, 500); mp, mr = evaluate(mem, rule, 600)
        dem = pp >= 0.6 and (pp - mp) >= 0.2
        if dem: torch.save(ppo.state_dict(), f"{here}/metaworld_{rule}.pt")
        results.append({"rule": rule, "ppo_prec": round(pp, 2), "ppo_netR": round(pr, 1),
                        "mem_prec": round(mp, 2), "mem_netR": round(mr, 1), "gap": round(pp - mp, 2), "demanding": dem})
        json.dump(results, open(out, "w"), indent=1)
        print(f"{rule:7s}: PPO prec {pp:.2f} netR {pr:+.1f} | memoryless prec {mp:.2f} netR {mr:+.1f} | gap {pp-mp:+.2f} {'DEMANDING' if dem else ''}", flush=True)
    dem = [x["rule"] for x in results if x["demanding"]]
    print(f"\nDEMANDING (PPO prec>=0.6 AND gap>=0.2): {dem}", flush=True)
    json.dump(results, open(out, "w"), indent=1); print("saved metaworld_sweep.json", flush=True)
