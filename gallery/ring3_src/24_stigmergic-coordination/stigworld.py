"""Phase-2.5 substrate widening (axis 5): STIGMERGIC / environment-mediated coordination.

The three role-partitioning competencies (division of labor, size regulation, turn-taking) all fingerprint into
one saturated region: identity + observation + memory. To find a NEW mechanism the search must load DIFFERENT
channels. Stigmergy is coordination through MARKS LEFT IN THE ENVIRONMENT (like ant pheromone trails): agents do
NOT see each other; they deposit onto and sense a shared field, and the field itself carries the coordination
state. The hypothesis is that this INVERTS the role-partitioning fingerprint -- the environment becomes the memory,
so a `nomark` ablation (field never persists) is load-bearing while identity (noid) is irrelevant and recurrent
memory (memwipe) matters less.

Task: n agents on an NxN torus must COVER a set of target cells (visit each at least once within the episode).
Agents cannot observe other agents at all. They observe only: their own local view of the shared PHEROMONE field
(a decaying scalar deposited by anyone who has been near) and whether a target is at/near their cell. To cover
efficiently and avoid all piling onto the same targets, they must read and write the field: cells already visited
carry pheromone, steering others toward unvisited targets. That division of exploration is mediated ENTIRELY by
the environment, not by seeing peers or by identity.

Ablation channels (battery introspects __init__):
  - nomark  : the pheromone field is zeroed every step (no environmental persistence) -> the coordination medium
              is destroyed. THE key channel; if stigmergic, this collapses coverage.
  - blind   : zero the local field + target view -> agents can't sense anything -> random walk
  - noid    : zero own identity one-hot -> should be HARMLESS if truly stigmergic (agents interchangeable)
  - memwipe : reset recurrent state each step -> should matter LESS than for role-partitioning (state is in the
              field, not the GRU)

Metric: cover_score() = fraction of target cells visited by episode end, averaged over batch. High = efficient
stigmergic coverage; near-chance = no environmental coordination.
"""
import numpy as np, torch, torch.nn as nn

N = 7; T = 24; H = 96
DIRS = np.array([[0, 0], [0, 1], [0, -1], [1, 0], [-1, 0]])   # stay, E, W, S, N
NMOVE = 5
DECAY = 0.85            # pheromone decay per step
DEPOSIT = 1.0
torch.set_num_threads(6)


class StigWorld:
    def __init__(self, B, seed, cfg=None, blind=False, memwipe=False, noid=False, nocost=False, nomark=False):
        self.B = B; self.cfg = cfg or {}
        self.blind = blind; self.memwipe = memwipe; self.noid = noid; self.nomark = nomark
        self.n = self.cfg.get("n_ag", 4); self.n_tgt = self.cfg.get("n_tgt", 12)
        self.decay = self.cfg.get("decay", DECAY)
        global N
        N = self.cfg.get("grid", 11)                                  # set grid size for this run (module-global, all instances share)
        # obs per agent: own id(n) + local 3x3 pheromone patch(9) + local 3x3 target patch(9) + own pos (2) + progress(1)
        self.odim = self.n + 9 + 9 + 2 + 1
        self.rng = np.random.default_rng(seed); self.reset()

    def reset(self):
        B, n = self.B, self.n
        self.ap = self.rng.integers(0, N, size=(B, n, 2))
        self.field = np.zeros((B, N, N), np.float32)                 # shared pheromone
        # target cells
        self.tgt = np.zeros((B, N, N), np.float32)
        self.visited = np.zeros((B, N, N), bool)
        for b in range(B):
            idx = self.rng.choice(N*N, size=self.n_tgt, replace=False)
            for c in idx:
                self.tgt[b, c//N, c%N] = 1.0
        self._t = 0
        return self.obs()

    def _patch(self, grid, pos):
        """3x3 toroidal patch of `grid` around each agent position. grid[B,N,N], pos[B,2] -> [B,9]."""
        B = self.B; out = np.zeros((B, 9), np.float32); o = 0
        for dr in (-1, 0, 1):
            for dc in (-1, 0, 1):
                r = (pos[:, 0] + dr) % N; c = (pos[:, 1] + dc) % N
                out[:, o] = grid[np.arange(B), r, c]; o += 1
        return out

    def obs(self):
        B, n = self.B, self.n
        out = np.zeros((B, n, self.odim), np.float32)
        fld = np.zeros_like(self.field) if (self.blind or self.nomark) else self.field
        tg = np.zeros_like(self.tgt) if self.blind else self.tgt
        for i in range(n):
            o = 0
            if not self.noid: out[:, i, i] = 1.0
            o += n
            out[:, i, o:o+9] = self._patch(fld, self.ap[:, i]); o += 9
            out[:, i, o:o+9] = self._patch(tg, self.ap[:, i]); o += 9
            out[:, i, o:o+2] = self.ap[:, i] / N; o += 2
            out[:, i, o] = self._t / T; o += 1
        return out

    def step(self, mv, sg=None):
        B, n = self.B, self.n
        self._t += 1
        for i in range(n):
            self.ap[:, i] = (self.ap[:, i] + DIRS[mv[:, i]]) % N
        # decay + deposit pheromone at agent cells
        self.field *= self.decay
        for i in range(n):
            self.field[np.arange(B), self.ap[:, i, 0], self.ap[:, i, 1]] += DEPOSIT
        if self.nomark:
            self.field[:] = 0.0                                       # environmental memory destroyed
        # mark visited targets; reward = NEW targets covered this step (shared)
        new = np.zeros(B, np.float32)
        for i in range(n):
            r, c = self.ap[:, i, 0], self.ap[:, i, 1]
            hit = (self.tgt[np.arange(B), r, c] > 0) & (~self.visited[np.arange(B), r, c])
            self.visited[np.arange(B), r, c] |= hit
            new += hit.astype(np.float32)
        r = new / self.n_tgt
        rew = np.repeat(r[:, None], n, axis=1).astype(np.float32)
        return rew

    def cover_score(self):
        covered = (self.visited & (self.tgt > 0)).reshape(self.B, -1).sum(1)
        return covered / self.n_tgt


class Pol(nn.Module):
    def __init__(self, odim):
        super().__init__(); self.gru = nn.GRU(odim, H, batch_first=True)
        self.a = nn.Linear(H, NMOVE); self.v = nn.Linear(H, 1)

    def forward(self, x, h=None):
        z, h = self.gru(x, h); return self.a(z), self.v(z).squeeze(-1), h


def gae(Rw, V, gamma=0.96, lam=0.95):
    B, Tt = Rw.shape; adv = np.zeros((B, Tt), np.float32); last = np.zeros(B, np.float32)
    for t in reversed(range(Tt)):
        nextv = V[:, t+1] if t+1 < Tt else np.zeros(B, np.float32)
        delta = Rw[:, t] + gamma*nextv - V[:, t]; last = delta + gamma*lam*last; adv[:, t] = last
    return adv, adv + V


def rollout(net, env, greedy=False):
    B, n, od = env.B, env.n, env.odim; obs = env.obs(); h = [None]*n
    O = np.zeros((B, n, T, od), np.float32); A = np.zeros((B, n, T), int)
    LP = np.zeros((B, n, T), np.float32); V = np.zeros((B, n, T), np.float32); Rw = np.zeros((B, n, T), np.float32)
    memwipe = getattr(env, "memwipe", False)
    for t in range(T):
        mv = np.zeros((B, n), int)
        for i in range(n):
            hin = None if memwipe else h[i]
            with torch.no_grad():
                al, v, h[i] = net(torch.from_numpy(obs[:, i])[:, None, :], hin)
            d = torch.distributions.Categorical(logits=al[:, 0])
            a_ = al[:, 0].argmax(1) if greedy else d.sample()
            O[:, i, t] = obs[:, i]; A[:, i, t] = a_.numpy(); LP[:, i, t] = d.log_prob(a_).numpy(); V[:, i, t] = v[:, 0].numpy()
            mv[:, i] = a_.numpy()
        r = env.step(mv); Rw[:, :, t] = r; obs = env.obs()
    return O, A, LP, V, Rw


def train(cfg, iters, seed=0, B=192):
    torch.manual_seed(seed); od = StigWorld(1, 0, cfg).odim; n = cfg["n_ag"]
    net = Pol(od); opt = torch.optim.Adam(net.parameters(), lr=3e-3)
    for it in range(iters):
        env = StigWorld(B, 1000 + it, cfg)
        O, A, LP, V, Rw = rollout(net, env)
        O2 = O.reshape(B*n, T, od); A2 = A.reshape(B*n, T)
        LP2, V2, R2 = LP.reshape(B*n, T), V.reshape(B*n, T), Rw.reshape(B*n, T)
        adv, ret = gae(R2, V2); adv = (adv - adv.mean()) / (adv.std() + 1e-8)
        Ot, At, LPt, advt, rett = map(torch.from_numpy, (O2, A2, LP2, adv, ret))
        for _ in range(4):
            al, v, _ = net(Ot); d = torch.distributions.Categorical(logits=al)
            ratio = torch.exp(d.log_prob(At) - LPt)
            s1 = ratio*advt; s2 = torch.clamp(ratio, 0.8, 1.2)*advt
            loss = -torch.min(s1, s2).mean() + 0.5*((v-rett)**2).mean() - 0.01*d.entropy().mean()
            opt.zero_grad(); loss.backward(); nn.utils.clip_grad_norm_(net.parameters(), 0.5); opt.step()
    return net


def eval_cover(net, cfg, seed, B=500, **kw):
    env = StigWorld(B, seed, cfg, **kw); rollout(net, env, greedy=True)
    return float(env.cover_score().mean())
