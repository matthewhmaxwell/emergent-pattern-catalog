"""Phase-2.5 substrate widening: ANTI-coordination / task-allocation world.

The openhunt/paramworld family only rewards CONVERGENCE (all agents meet on one cell = rendezvous). Every
catalogued competency lives inside that pro-coordination assumption. This world expresses the behavior that
family structurally CANNOT: agents must reach DISTINCT targets - a division-of-labor / task-allocation problem.

Task: there are n agents and n target cells. Reward is maximized when EACH target is covered by EXACTLY ONE
agent (a perfect matching). Two agents on the same target, or a target left uncovered, scores worse. Solving
this requires the agents to SPLIT UP and specialize - anti-coordination, not rendezvous.

Two capability regimes (the heterogeneity knob):
  - homogeneous (hetero=False): all agents identical; must break symmetry to allocate (who takes which target).
  - heterogeneous (hetero=True): each target has a "type" and each agent a matching type-affinity, so an
    efficient allocation is by capability (role specialization), not just arbitrary symmetry-breaking.

Ablation channels the battery can fingerprint (introspected from __init__):
  - blind   : zero others' relative positions (can't see partners -> can't avoid doubling up)
  - noid    : zero the agent's own identity one-hot (can't tell which agent it is -> can't hold a stable role)
  - mute    : scramble the discrete signal (if a channel is used to negotiate the allocation)
  - memwipe : reset recurrent state each step (if allocation is integrated over time)

Metric: alloc_score() = mean fraction of targets covered by exactly one agent at episode end (1.0 = perfect
division of labor; lower = collisions or gaps).
"""
import numpy as np, torch, torch.nn as nn

N = 9; K = 3; T = 30; H = 96; NMOVE = 5; NTYPE = 2
DIRS = np.array([(0, 0), (1, 0), (-1, 0), (0, 1), (0, -1)])
torch.set_num_threads(8)


def torus_rel(a, b):
    return ((b - a + N // 2) % N) - N // 2


class AllocWorld:
    def __init__(self, B, seed, cfg=None, blind=False, mute=False, memwipe=False, noid=False):
        self.B = B; self.cfg = cfg
        self.blind = blind; self.mute = mute; self.memwipe = memwipe; self.noid = noid
        self.n = cfg["n_ag"]; self.has_obs = cfg["has_obs"]; self.has_sig = cfg["has_sig"]
        self.hetero = cfg.get("hetero", False); self.reward_shared = cfg.get("reward_shared", True)
        self.asym_init = cfg.get("asym_init", False)  # agents start in distinct spatial bands (symmetry-break)
        self._t = 0
        # obs per agent: own id one-hot(n) + [rel to each target (2 each)] + [target type one-hot(NTYPE) each if hetero]
        #                + own type one-hot(NTYPE if hetero) + [each other agent rel (2 each)] + [others' sig (K each)]
        self.odim = (self.n                                   # own identity
                     + self.n * 2                             # rel to each of n targets
                     + (self.n * NTYPE if self.hetero else 0) # each target's type
                     + (NTYPE if self.hetero else 0)          # own type
                     + (self.n - 1) * 2                       # rel to each other agent
                     + (self.n - 1) * K)                      # others' signals
        self.rng = np.random.default_rng(seed); self.reset()

    def reset(self):
        B, n = self.B, self.n
        if self.asym_init:
            # each agent starts in a distinct vertical band => positional symmetry between agents is broken,
            # but the TASK (reach distinct targets) and REWARD are unchanged; agents must still DISCOVER to split.
            self.ap = np.zeros((B, n, 2), int)
            band = N / n
            for i in range(n):
                self.ap[:, i, 0] = self.rng.integers(int(i*band), max(int(i*band)+1, int((i+1)*band)), size=B) % N
                self.ap[:, i, 1] = self.rng.integers(0, N, size=B)
        else:
            self.ap = self.rng.integers(0, N, size=(B, n, 2))
        self.tg = self.rng.integers(0, N, size=(B, n, 2))            # n target cells
        self.ttype = self.rng.integers(0, NTYPE, size=(B, n))        # each target's type
        self.atype = self.rng.integers(0, NTYPE, size=(B, n))        # each agent's type-affinity
        self.lsig = np.zeros((B, n), int)
        return self.obs()

    def obs(self):
        B, n = self.B, self.n
        out = np.zeros((B, n, self.odim), np.float32)
        for i in range(n):
            o = 0
            # own identity one-hot (the noid channel zeros this)
            if not self.noid:
                out[:, i, i] = 1.0
            o += n
            # relative position to each target
            for k in range(n):
                if self.blind or not self.has_obs:
                    rel = np.zeros((B, 2), np.float32)
                else:
                    rel = torus_rel(self.ap[:, i], self.tg[:, k]) / (N // 2)
                out[:, i, o:o+2] = rel; o += 2
            # target types + own type (heterogeneous regime only)
            if self.hetero:
                for k in range(n):
                    out[np.arange(B), i, o + self.ttype[:, k]] = 1.0; o += NTYPE
                out[np.arange(B), i, o + self.atype[:, i]] = 1.0; o += NTYPE
            # relative position to each OTHER agent (to avoid doubling up)
            for j in range(n):
                if j == i:
                    continue
                if self.blind or not self.has_obs:
                    rel = np.zeros((B, 2), np.float32)
                else:
                    rel = torus_rel(self.ap[:, i], self.ap[:, j]) / (N // 2)
                out[:, i, o:o+2] = rel; o += 2
            # others' signals
            for j in range(n):
                if j == i:
                    continue
                if self.has_sig:
                    sj = self.rng.integers(0, K, size=B) if self.mute else self.lsig[:, j]
                    out[np.arange(B), i, o + sj] = 1.0
                o += K
        return out

    def step(self, mv, sg):
        B, n = self.B, self.n
        self._t += 1
        for i in range(n):
            self.ap[:, i] = (self.ap[:, i] + DIRS[mv[:, i]]) % N
        self.lsig = sg.copy()
        # per-step shaping is off; reward is computed at episode end via alloc_score-style credit
        rew = np.full((B, n), -0.01, np.float32)
        # count coverage each step: reward agents that sit alone on a target (encourages splitting)
        for i in range(n):
            on_t = np.zeros(B, bool)
            for k in range(n):
                on = np.all(self.ap[:, i] == self.tg[:, k], axis=1)
                # bonus if this agent covers target k AND type matches (hetero) / any (homo)
                if self.hetero:
                    match = (self.atype[:, i] == self.ttype[:, k])
                    on_t |= (on & match)
                else:
                    on_t |= on
            # collision penalty: two agents on same cell
            same = np.zeros(B, bool)
            for j in range(n):
                if j == i: continue
                same |= np.all(self.ap[:, i] == self.ap[:, j], axis=1)
            rew[:, i] += 0.1 * on_t.astype(np.float32) - 0.05 * same.astype(np.float32)
        if self.reward_shared:
            rew[:] = rew.mean(axis=1, keepdims=True)
        return rew

    def alloc_score(self):
        """Fraction of targets covered by exactly one (type-matching, if hetero) agent, at current state."""
        B, n = self.B, self.n
        cov = np.zeros((B, n), int)  # how many agents cover each target
        for k in range(n):
            for i in range(n):
                on = np.all(self.ap[:, i] == self.tg[:, k], axis=1)
                if self.hetero:
                    on = on & (self.atype[:, i] == self.ttype[:, k])
                cov[:, k] += on.astype(int)
        exactly_one = (cov == 1).sum(axis=1) / n
        return float(exactly_one.mean())


ID_DIM = 8  # width of the learned identity embedding on the skip path

class Pol(nn.Module):
    def __init__(self, odim, n=None, id_skip=False):
        super().__init__(); self.gru = nn.GRU(odim, H, batch_first=True)
        # id_skip: route a learned per-agent identity embedding AROUND the GRU straight to the heads, so the
        # fixed role token is available at decision time regardless of what the (churn-prone) recurrent state does.
        # The identity one-hot is the first n dims of obs; noid ablation zeros it, so the skip goes dead under noid.
        self.id_skip = id_skip; self.n = n
        head_in = H + (ID_DIM if id_skip else 0)
        if id_skip:
            self.id_emb = nn.Linear(n, ID_DIM, bias=False)
        self.mv = nn.Linear(head_in, NMOVE); self.sg = nn.Linear(head_in, K); self.v = nn.Linear(head_in, 1)

    def forward(self, x, h=None):
        z, h = self.gru(x, h)
        if self.id_skip:
            idvec = self.id_emb(x[..., :self.n])  # constant-per-agent learned vector, memory-independent
            z = torch.cat([z, idvec], dim=-1)
        return self.mv(z), self.sg(z), self.v(z).squeeze(-1), h


def gae(R, V, gamma=0.97, lam=0.95):
    B, Tt = R.shape; adv = np.zeros((B, Tt), np.float32); last = np.zeros(B, np.float32)
    for t in reversed(range(Tt)):
        nextv = V[:, t+1] if t+1 < Tt else np.zeros(B, np.float32)
        delta = R[:, t] + gamma*nextv - V[:, t]; last = delta + gamma*lam*last; adv[:, t] = last
    return adv, adv + V


def rollout(net, env, greedy=False):
    B, n, od = env.B, env.n, env.odim; obs = env.obs(); h = [None]*n
    O = np.zeros((B, n, T, od), np.float32); MV = np.zeros((B, n, T), int); SG = np.zeros((B, n, T), int)
    LP = np.zeros((B, n, T), np.float32); V = np.zeros((B, n, T), np.float32); R = np.zeros((B, n, T), np.float32)
    memwipe = getattr(env, "memwipe", False)
    for t in range(T):
        mv = np.zeros((B, n), int); sg = np.zeros((B, n), int)
        for i in range(n):
            hin = None if memwipe else h[i]
            with torch.no_grad():
                ml, sl, v, h[i] = net(torch.from_numpy(obs[:, i])[:, None, :], hin)
            md = torch.distributions.Categorical(logits=ml[:, 0]); sd = torch.distributions.Categorical(logits=sl[:, 0])
            m = ml[:, 0].argmax(1) if greedy else md.sample(); s = sl[:, 0].argmax(1) if greedy else sd.sample()
            O[:, i, t] = obs[:, i]; MV[:, i, t] = m.numpy(); SG[:, i, t] = s.numpy()
            LP[:, i, t] = (md.log_prob(m) + sd.log_prob(s)).numpy(); V[:, i, t] = v[:, 0].numpy()
            mv[:, i] = m.numpy(); sg[:, i] = s.numpy()
        r = env.step(mv, sg); R[:, :, t] = r; obs = env.obs()
    return O, MV, SG, LP, V, R


def train(cfg, iters, seed=0, B=192):
    torch.manual_seed(seed); od = AllocWorld(1, 0, cfg).odim
    net = Pol(od); opt = torch.optim.Adam(net.parameters(), lr=3e-3); n = cfg["n_ag"]
    for it in range(iters):
        env = AllocWorld(B, 1000 + it, cfg)
        O, MV, SG, LP, V, R = rollout(net, env)
        O2 = O.reshape(B*n, T, od); MV2, SG2 = MV.reshape(B*n, T), SG.reshape(B*n, T)
        LP2, V2, R2 = LP.reshape(B*n, T), V.reshape(B*n, T), R.reshape(B*n, T)
        adv, ret = gae(R2, V2); adv = (adv - adv.mean()) / (adv.std() + 1e-8)
        Ot, MVt, SGt, LPt = torch.from_numpy(O2), torch.from_numpy(MV2), torch.from_numpy(SG2), torch.from_numpy(LP2)
        advt, rett = torch.from_numpy(adv), torch.from_numpy(ret)
        for _ in range(4):
            ml, sl, v, _ = net(Ot)
            md = torch.distributions.Categorical(logits=ml); sd = torch.distributions.Categorical(logits=sl)
            ratio = torch.exp(md.log_prob(MVt) + sd.log_prob(SGt) - LPt)
            s1 = ratio*advt; s2 = torch.clamp(ratio, 0.8, 1.2)*advt
            loss = -torch.min(s1, s2).mean() + 0.5*((v-rett)**2).mean() - 0.01*(md.entropy()+sd.entropy()).mean()
            opt.zero_grad(); loss.backward(); nn.utils.clip_grad_norm_(net.parameters(), 0.5); opt.step()
    return net


def eval_alloc(net, cfg, seed, B=400, **kw):
    env = AllocWorld(B, seed, cfg, **kw)
    # roll out greedily to a terminal state and score the final allocation, averaged over the last few steps
    O, MV, SG, LP, V, R = rollout(net, env, greedy=True)
    return env.alloc_score()
