"""EPC Ring-3 pre-registration v1 (osf.io/up7m9) — environments, trainer and evaluation for P1–P4.

Trainer = dial2_ppo.py at the frozen commit ac27a8d, generalized to four tasks with IDENTICAL PPO settings:
GRU hidden 96, batch 192, Adam lr 3e-3, 4 PPO epochs/iter, clip 0.2, value coef 0.5, entropy coef 0.02,
grad-clip 0.5, gamma 0.97, lambda 0.95, advantage normalization, team reward, parameter sharing, 1500 iters.
Env-specific additions only: action-head masks (P1: talk step uses only the signal head, act step only the move head;
P3: single agent, no signal head), ablations, fair baselines, scripted oracle ceilings, cross-play.
Evaluation samples both action heads (DEVIATIONS.md D1). Written before any data; hashed before the first run.

  python prereg_envs.py smoke                                   # code check + chance/ceiling levels (no task data)
  python prereg_envs.py train P1 full 1                         # task, variant, seed
  python prereg_envs.py train P3 full 1 --rho 0.5
  python prereg_envs.py crossplay                               # P4 self/cross-play over saved nets
"""
import numpy as np, sys, json, os, torch, torch.nn as nn

torch.set_num_threads(int(os.environ.get("PREREG_THREADS", "2")))
H, K, N = 96, 8, 7
DIRS = np.array([(0, 0), (1, 0), (-1, 0), (0, 1), (0, -1)])
HERE = os.path.dirname(os.path.abspath(__file__))
OUT = os.path.join(HERE, "prereg_runs")
ITERS = int(os.environ.get("PREREG_ITERS", "1500"))
EVAL_B, EVAL_REPS, EVAL_SEED0 = 400, 5, 900000          # 5 x 400 = 2000 evaluation episodes, disjoint seeds


class Pol(nn.Module):
    def __init__(self, odim, nmove, k):
        super().__init__(); self.gru = nn.GRU(odim, H, batch_first=True)
        self.mv = nn.Linear(H, nmove); self.sg = nn.Linear(H, max(k, 1)); self.v = nn.Linear(H, 1)

    def forward(self, x, h=None):
        z, h = self.gru(x, h); return self.mv(z), self.sg(z), self.v(z).squeeze(-1), h


class RandPol(nn.Module):                                    # uniform random policy (chance level)
    def __init__(self, nmove, k): super().__init__(); self.nm, self.k = nmove, max(k, 1)

    def forward(self, x, h=None):
        B, Tq = x.shape[:2]
        return torch.zeros(B, Tq, self.nm), torch.zeros(B, Tq, self.k), torch.zeros(B, Tq), None


def _onehot(o, rows, col0, idx):
    o[rows, col0 + idx] = 1.0


# ------------------------------------------------------------------ P1: symmetric anti-coordination (one-shot)
class P1Env:
    n, nmove, k, T = 2, 2, K, 2
    odim = 2 + K + K                                          # phase(2), own symbol, partner symbol

    def __init__(self, B, seed, mute=False, channel=True, **_):
        self.B, self.mute, self.channel = B, mute, channel; self.rng = np.random.default_rng(seed)
        self.t = 0; self.sym = np.zeros((B, 2), int); self.succ = np.zeros(B)

    def obs(self):
        B = self.B; o = np.zeros((B, 2, self.odim), np.float32); r = np.arange(B)
        if self.t == 0:
            o[:, :, 0] = 1.0
        else:
            o[:, :, 1] = 1.0
            for i in range(2):
                oi = o[:, i]; _onehot(oi, r, 2, self.sym[:, i])
                if self.channel:
                    pj = self.rng.integers(0, K, size=B) if self.mute else self.sym[:, 1 - i]
                    _onehot(oi, r, 2 + K, pj)
        return o

    def step(self, mv, sg):
        rew = np.zeros(self.B, np.float32)
        if self.t == 0:
            self.sym = sg.copy()
        else:
            d = (mv[:, 0] != mv[:, 1]).astype(np.float32); rew += d; self.succ = d
        self.t += 1; return rew

    def mask(self, t): return (t == 1, t == 0)               # (move/act head active, signal head active)

    def score(self): return float(self.succ.mean())


# ------------------------------------------------------------------ P2: asymmetry (alpha=1) + partner observation
class P2Env:
    n, nmove, k, T = 2, 5, K, 24
    odim = 2 + 2 + 2 + 1 + K                                  # own pos, partner pos, target(speaker), spk flag, signal

    def __init__(self, B, seed, mute=False, blind=False, channel=True, partner_obs=True, **_):
        self.B, self.mute, self.blind, self.channel, self.pobs = B, mute, blind, channel, partner_obs
        self.rng = np.random.default_rng(seed)
        self.ap = self.rng.integers(0, N, size=(B, 2, 2)); self.tgt = self.rng.integers(0, N, size=(B, 2))
        self.speaker = self.rng.integers(0, 2, size=B); self.lsig = np.zeros((B, 2), int); self.hit = np.zeros(B, int)

    def obs(self):
        B = self.B; o = np.zeros((B, 2, self.odim), np.float32); r = np.arange(B)
        for i in range(2):
            j = 1 - i
            o[:, i, 0:2] = self.ap[:, i] / (N - 1)
            if self.pobs:
                pp = self.rng.integers(0, N, size=(B, 2)) if self.blind else self.ap[:, j]
                o[:, i, 2:4] = pp / (N - 1)
            sp = (self.speaker == i)
            o[sp, i, 4:6] = self.tgt[sp] / (N - 1); o[:, i, 6] = sp.astype(np.float32)
            if self.channel:
                sj = self.rng.integers(0, K, size=B) if self.mute else self.lsig[:, j]
                oi = o[:, i]; _onehot(oi, r, 7, sj)
        return o

    def step(self, mv, sg):
        B = self.B; rew = np.full(B, -0.01, np.float32)
        for i in range(2): self.ap[:, i] = np.clip(self.ap[:, i] + DIRS[mv[:, i]], 0, N - 1)
        self.lsig = sg.copy()
        both = np.all(self.ap[:, 0] == self.tgt, axis=1) & np.all(self.ap[:, 1] == self.tgt, axis=1)
        hb = np.where(both)[0]
        if len(hb):
            rew[hb] += 1.0; self.hit[hb] += 1
            self.tgt[hb] = self.rng.integers(0, N, size=(len(hb), 2))
            self.speaker[hb] = self.rng.integers(0, 2, size=len(hb))
            self.ap[hb] = self.rng.integers(0, N, size=(len(hb), 2, 2))
        return rew

    def mask(self, t): return (True, True)

    def score(self): return float(self.hit.mean())


# ------------------------------------------------------------------ P3: memory-perception dial (single agent)
class P3Env:
    n, nmove, k, T = 1, 5, 1, 24
    odim = 2 + 2 + 2 + 2 + 2                                  # own pos, goal A, goal B, cue (t=0 only), landmark

    def __init__(self, B, seed, rho=1.0, **_):
        self.B = B; self.rng = np.random.default_rng(seed); self.t = 0
        cells = np.array([self.rng.choice(N * N, 3, replace=False) for _ in range(B)])
        self.pos = np.stack([cells[:, 0] // N, cells[:, 0] % N], 1)
        self.g = np.stack([np.stack([cells[:, 1] // N, cells[:, 1] % N], 1),
                           np.stack([cells[:, 2] // N, cells[:, 2] % N], 1)], 1)       # (B, 2 goals, 2)
        self.correct = self.rng.integers(0, 2, size=B); self.lm = self.rng.random(B) < rho
        self.done = np.zeros(B, bool); self.right = np.zeros(B, bool)

    def obs(self):
        B = self.B; o = np.zeros((B, 1, self.odim), np.float32); r = np.arange(B)
        o[:, 0, 0:2] = self.pos / (N - 1); o[:, 0, 2:4] = self.g[:, 0] / (N - 1); o[:, 0, 4:6] = self.g[:, 1] / (N - 1)
        if self.t == 0: o[r, 0, 6 + self.correct] = 1.0
        o[r[self.lm], 0, 8 + self.correct[self.lm]] = 1.0
        return o

    def step(self, mv, sg):
        B = self.B; act = ~self.done; rew = np.where(act, -0.01, 0.0).astype(np.float32)
        self.pos[act] = np.clip(self.pos[act] + DIRS[mv[act, 0]], 0, N - 1)
        for gi in range(2):
            at = act & np.all(self.pos == self.g[:, gi], axis=1)
            ok = at & (self.correct == gi); bad = at & (self.correct != gi)
            rew[ok] += 1.0; rew[bad] -= 0.3; self.right |= ok; self.done |= at
        self.t += 1; return rew

    def mask(self, t): return (True, False)

    def score(self): return float(self.right.mean())


# ------------------------------------------------------------------ P4: blind rendezvous (rung-2 world), cross-play
class P4Env:
    n, nmove, k, T = 2, 5, K, 24
    odim = 2 + K                                              # own pos, partner signal

    def __init__(self, B, seed, mute=False, channel=True, **_):
        self.B, self.mute, self.channel = B, mute, channel; self.rng = np.random.default_rng(seed)
        self.ap = self._distinct(B); self.lsig = np.zeros((B, 2), int); self.hit = np.zeros(B, int)
        self.meets = np.zeros((N, N), int)                    # exploratory: where rendezvous happens

    def _distinct(self, m):
        c = np.array([self.rng.choice(N * N, 2, replace=False) for _ in range(m)])
        return np.stack([np.stack([c[:, 0] // N, c[:, 0] % N], 1), np.stack([c[:, 1] // N, c[:, 1] % N], 1)], 1)

    def obs(self):
        B = self.B; o = np.zeros((B, 2, self.odim), np.float32); r = np.arange(B)
        for i in range(2):
            o[:, i, 0:2] = self.ap[:, i] / (N - 1)
            if self.channel:
                sj = self.rng.integers(0, K, size=B) if self.mute else self.lsig[:, 1 - i]
                oi = o[:, i]; _onehot(oi, r, 2, sj)
        return o

    def step(self, mv, sg):
        B = self.B; rew = np.full(B, -0.01, np.float32)
        for i in range(2): self.ap[:, i] = np.clip(self.ap[:, i] + DIRS[mv[:, i]], 0, N - 1)
        self.lsig = sg.copy()
        hb = np.where(np.all(self.ap[:, 0] == self.ap[:, 1], axis=1))[0]
        if len(hb):
            np.add.at(self.meets, (self.ap[hb, 0, 0], self.ap[hb, 0, 1]), 1)
            rew[hb] += 1.0; self.hit[hb] += 1; self.ap[hb] = self._distinct(len(hb))
        return rew

    def mask(self, t): return (True, True)

    def score(self): return float(self.hit.mean())


ENVS = {"P1": P1Env, "P2": P2Env, "P3": P3Env, "P4": P4Env}
# training variants (fair baselines = trained WITHOUT the channel)
VARIANTS = {"P1": {"full": {}, "nochannel": {"channel": False}},
            "P2": {"full": {}, "noobs": {"partner_obs": False}, "nochannel": {"channel": False}},
            "P3": {"full": {}, "memoryless": {"_memoryless": True}},
            "P4": {"full": {}, "nochannel": {"channel": False}}}
# evaluation-time ablations applied to a trained model
ABLATIONS = {"P1": {"mute": {"mute": True}}, "P2": {"blind": {"blind": True}, "mute": {"mute": True}},
             "P3": {"memwipe": {"_memwipe": True}}, "P4": {"mute": {"mute": True}}}


# ------------------------------------------------------------------ rollout / train / evaluate (dial2_ppo.py core)
def rollout(nets, env, memoryless=False):
    B, n, T = env.B, env.n, env.T; obs = env.obs(); h = [None] * n
    O = np.zeros((B, n, T, env.odim), np.float32); MV = np.zeros((B, n, T), int); SG = np.zeros((B, n, T), int)
    LP = np.zeros((B, n, T), np.float32); V = np.zeros((B, n, T), np.float32); R = np.zeros((B, T), np.float32)
    MM = np.zeros(T, np.float32); MS = np.zeros(T, np.float32)
    for t in range(T):
        am, asg = env.mask(t); MM[t], MS[t] = float(am), float(asg)
        mv = np.zeros((B, n), int); sg = np.zeros((B, n), int)
        for i in range(n):
            with torch.no_grad():
                ml, sl, v, hn = nets[i](torch.from_numpy(obs[:, i])[:, None, :], None if memoryless else h[i])
            h[i] = hn
            md = torch.distributions.Categorical(logits=ml[:, 0]); sd = torch.distributions.Categorical(logits=sl[:, 0])
            m, s = md.sample(), sd.sample()                     # sampled for both heads (D1)
            O[:, i, t] = obs[:, i]; MV[:, i, t] = m.numpy(); SG[:, i, t] = s.numpy()
            LP[:, i, t] = (md.log_prob(m) * am + sd.log_prob(s) * asg).numpy(); V[:, i, t] = v[:, 0].numpy()
            mv[:, i] = m.numpy(); sg[:, i] = s.numpy()
        R[:, t] = env.step(mv, sg); obs = env.obs()
    return O, MV, SG, LP, V, R, MM, MS


def gae(R, V, gamma=0.97, lam=0.95):
    B, Tt = R.shape; adv = np.zeros((B, Tt), np.float32); last = np.zeros(B, np.float32)
    for t in reversed(range(Tt)):
        nextv = V[:, t + 1] if t + 1 < Tt else np.zeros(B, np.float32)
        delta = R[:, t] + gamma * nextv - V[:, t]; last = delta + gamma * lam * last; adv[:, t] = last
    return adv, adv + V


def make(task, **kw):
    E = ENVS[task]; envkw = {k: v for k, v in kw.items() if not k.startswith("_")}
    return lambda B, seed, **abl: E(B, seed, **{**envkw, **{k: v for k, v in abl.items() if not k.startswith("_")}})


def train(task, variant, seed, rho=None, iters=ITERS, B=192, log=True):
    E = ENVS[task]; vkw = dict(VARIANTS[task][variant]); memoryless = vkw.pop("_memoryless", False)
    if rho is not None: vkw["rho"] = rho
    mk = make(task, **vkw); torch.manual_seed(seed); net = Pol(E.odim, E.nmove, E.k)
    opt = torch.optim.Adam(net.parameters(), lr=3e-3); n, T = E.n, E.T
    for it in range(iters):
        env = mk(B, 1000 + seed * 99991 + it)
        O, MV, SG, LP, V, R, MM, MS = rollout([net] * n, env, memoryless=memoryless)
        Rag = np.repeat(R[:, None, :], n, axis=1)
        O2 = O.reshape(B * n, T, E.odim); MV2, SG2 = MV.reshape(B * n, T), SG.reshape(B * n, T)
        LP2, V2, R2 = LP.reshape(B * n, T), V.reshape(B * n, T), Rag.reshape(B * n, T)
        adv, ret = gae(R2, V2); adv = (adv - adv.mean()) / (adv.std() + 1e-8)
        Ot = torch.from_numpy(O2); MVt, SGt, LPt = torch.from_numpy(MV2), torch.from_numpy(SG2), torch.from_numpy(LP2)
        advt, rett = torch.from_numpy(adv), torch.from_numpy(ret); MMt, MSt = torch.from_numpy(MM), torch.from_numpy(MS)
        for _ in range(4):
            if memoryless:                                   # no recurrence: every step processed from h=0
                ml, sl, v, _ = net(Ot.reshape(B * n * T, 1, E.odim))
                ml, sl, v = ml.reshape(B * n, T, -1), sl.reshape(B * n, T, -1), v.reshape(B * n, T)
            else:
                ml, sl, v, _ = net(Ot)
            md = torch.distributions.Categorical(logits=ml); sd = torch.distributions.Categorical(logits=sl)
            ratio = torch.exp(md.log_prob(MVt) * MMt + sd.log_prob(SGt) * MSt - LPt)
            s1 = ratio * advt; s2 = torch.clamp(ratio, 0.8, 1.2) * advt
            ent = md.entropy() * MMt + sd.entropy() * MSt
            loss = -torch.min(s1, s2).mean() + 0.5 * ((v - rett) ** 2).mean() - 0.02 * ent.mean()
            opt.zero_grad(); loss.backward(); nn.utils.clip_grad_norm_(net.parameters(), 0.5); opt.step()
        if log and (it % 150 == 0 or it == iters - 1):
            print(f"  [{task}/{variant}{'' if rho is None else f'/rho={rho}'} s{seed}] iter {it:>4}: "
                  f"S {evaluate(task, [net] * n, vkw, memoryless=memoryless, reps=1, B=200):.3f}", flush=True)
    return net, vkw, memoryless


def evaluate(task, nets, envkw, memoryless=False, reps=EVAL_REPS, B=EVAL_B, **abl):
    mk = make(task, **envkw); memwipe = abl.pop("_memwipe", False); out = []
    for r in range(reps):
        env = mk(B, EVAL_SEED0 + r, **abl); rollout(nets, env, memoryless=memoryless or memwipe); out.append(env.score())
    return float(np.mean(out))


# ------------------------------------------------------------------ scripted oracle ceilings + chance (smoke)
def _toward(p, g):
    d = g - p; m = np.zeros(len(p), int)
    m[d[:, 0] > 0] = 1; m[d[:, 0] < 0] = 2
    zr = d[:, 0] == 0; m[zr & (d[:, 1] > 0)] = 3; m[zr & (d[:, 1] < 0)] = 4
    return m


def _bfs_move(p, g, block):
    """first move of a shortest grid path p->g that never enters `block` (the wrong goal)."""
    from collections import deque
    start, goal, blk = tuple(int(x) for x in p), tuple(int(x) for x in g), tuple(int(x) for x in block)
    if start == goal: return 0
    prev = {start: None}; q = deque([start])
    while q:
        c = q.popleft()
        if c == goal: break
        for k in range(1, 5):
            nb = (c[0] + DIRS[k][0], c[1] + DIRS[k][1])
            if 0 <= nb[0] < N and 0 <= nb[1] < N and nb != blk and nb not in prev:
                prev[nb] = (c, k); q.append(nb)
    c = goal
    while prev[c][0] != start: c = prev[c][0]
    return prev[c][1]


def oracle_score(task, rho=1.0, reps=EVAL_REPS, B=EVAL_B):
    out = []
    for r in range(reps):
        rng = np.random.default_rng(EVAL_SEED0 + 77 + r)
        if task == "P1":
            env = P1Env(B, EVAL_SEED0 + r); s = rng.integers(0, K, size=(B, 2)); env.step(np.zeros((B, 2), int), s)
            a0 = np.where(s[:, 0] > s[:, 1], 0, np.where(s[:, 0] < s[:, 1], 1, rng.integers(0, 2, size=B)))
            a1 = np.where(s[:, 1] > s[:, 0], 0, np.where(s[:, 1] < s[:, 0], 1, rng.integers(0, 2, size=B)))
            env.step(np.stack([a0, a1], 1), np.zeros((B, 2), int))
        elif task == "P2":
            env = P2Env(B, EVAL_SEED0 + r)
            for _ in range(env.T):
                spk = env.speaker; m = np.zeros((B, 2), int)
                for i in range(2):
                    goal = np.where((spk == i)[:, None], env.tgt, env.ap[:, 1 - i]); m[:, i] = _toward(env.ap[:, i], goal)
                env.step(m, np.zeros((B, 2), int))
        elif task == "P3":
            env = P3Env(B, EVAL_SEED0 + r, rho=rho); ix = np.arange(B)
            for _ in range(env.T):
                g, w = env.g[ix, env.correct], env.g[ix, 1 - env.correct]
                m = np.array([0 if env.done[b] else _bfs_move(env.pos[b], g[b], w[b]) for b in ix])
                env.step(m[:, None], np.zeros((B, 1), int))
        else:
            env = P4Env(B, EVAL_SEED0 + r); c = np.array([N // 2, N // 2])
            for _ in range(env.T):
                m = np.stack([_toward(env.ap[:, i], np.tile(c, (B, 1))) for i in range(2)], 1); env.step(m, np.zeros((B, 2), int))
        out.append(env.score())
    return float(np.mean(out))


def smoke():
    os.makedirs(OUT, exist_ok=True); res = {}
    for task, E in ENVS.items():
        rp = RandPol(E.nmove, E.k)
        chance = evaluate(task, [rp] * E.n, {"rho": 1.0} if task == "P3" else {})
        ceil = oracle_score(task)
        train(task, "full", 0, rho=1.0 if task == "P3" else None, iters=50, log=False)   # execution check only
        res[task] = {"chance_random_policy": round(chance, 4), "ceiling_oracle": round(ceil, 4)}
        print(f"{task}: random-policy chance {chance:.4f} | scripted-oracle ceiling {ceil:.4f} | 50-iter train OK", flush=True)
    res["P3"]["chance_registered"] = 0.50                     # PREREG P3: 'random choice between goals = 0.50'
    json.dump(res, open(os.path.join(OUT, "smoke.json"), "w"), indent=1); print("SMOKE DONE", flush=True)


def run_train(task, variant, seed, rho=None):
    os.makedirs(OUT, exist_ok=True)
    tag = f"{task}_{variant}{'' if rho is None else f'_rho{rho}'}_s{seed}"
    net, vkw, memoryless = train(task, variant, seed, rho=rho)
    torch.save(net.state_dict(), os.path.join(OUT, tag + ".pt"))
    E = ENVS[task]; row = {"task": task, "variant": variant, "seed": seed, "rho": rho,
                           "S": evaluate(task, [net] * E.n, vkw, memoryless=memoryless)}
    if variant == "full":
        for name, abl in ABLATIONS[task].items():
            row["S_" + name] = evaluate(task, [net] * E.n, vkw, **dict(abl))
    json.dump(row, open(os.path.join(OUT, tag + ".json"), "w"), indent=1)
    print("DONE", tag, json.dumps(row), flush=True)


def crossplay():
    E = P4Env; res = {}
    for variant in ("full", "nochannel"):
        vkw = VARIANTS["P4"][variant]; nets = {}
        for s in (1, 2, 3):
            p = os.path.join(OUT, f"P4_{variant}_s{s}.pt")
            if os.path.exists(p):
                net = Pol(E.odim, E.nmove, E.k); net.load_state_dict(torch.load(p)); nets[s] = net
        pairs = {f"{a}x{b}": evaluate("P4", [nets[a], nets[b]], vkw) for a in nets for b in nets}
        modal = {}                                           # EXPLORATORY (not scored): modal self-play meeting cell
        for s, net in nets.items():
            env = make("P4", **vkw)(EVAL_B, EVAL_SEED0 + 500 + s); rollout([net, net], env)
            flat = int(env.meets.argmax()); tot = int(env.meets.sum())
            modal[str(s)] = {"cell": [flat // N, flat % N], "share": round(env.meets.max() / max(tot, 1), 3)}
        res[variant] = {"pairs": pairs, "exploratory_modal_meeting_cell": modal,
                        "self_play_mean": float(np.mean([v for k, v in pairs.items() if k[0] == k[-1]])),
                        "cross_play_mean": float(np.mean([v for k, v in pairs.items() if k[0] != k[-1]]))}
    json.dump(res, open(os.path.join(OUT, "P4_crossplay.json"), "w"), indent=1); print(json.dumps(res, indent=1))


if __name__ == "__main__":
    a = sys.argv
    if a[1] == "smoke":
        smoke()
    elif a[1] == "train":
        rho = float(a[a.index("--rho") + 1]) if "--rho" in a else None
        run_train(a[2], a[3], int(a[4]), rho=rho)
    elif a[1] == "crossplay":
        crossplay()
