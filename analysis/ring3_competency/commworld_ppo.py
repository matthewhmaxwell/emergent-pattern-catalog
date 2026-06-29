"""Phase D1: the richer world (multi-agent) -- validate the substrate + reach the first competency that
CANNOT exist single-agent: emergent COMMUNICATION. A referential signaling game embedded in a gridworld:
  - 3 colored goals at random cells; one color is the TARGET (sampled per episode, hidden from the listener).
  - SPEAKER observes the target color and emits a symbol (1 of K) -- it has no other action.
  - LISTENER observes the speaker's symbol + the direction to each colored goal, and navigates. Reward
    (shared): +1 if the listener reaches the TARGET-color goal, -0.3 wrong goal, -0.02/step.
There is no built-in symbol->color mapping; speaker and listener must CO-ADAPT a convention (Lewis
signaling). A single agent cannot do this. Joint PPO trains both.

Verification (the multi-agent analogue of feedback-ablation): SCRAMBLE the channel (feed the listener a
random symbol). If success collapses to chance (1/G), the symbol genuinely CARRIED the target info =
genuine communication, not a shortcut. A no-channel listener is the MCC baseline.

Run on VPS epc-venv: python commworld_ppo.py [--iters 400]
"""
import numpy as np, sys, json, os, torch, torch.nn as nn

N = 7; G = 3; K = 3; MAXT = 25; H = 64
LOBS = K + G * 2                 # listener obs: symbol one-hot + dir-to-each-color-goal
DIRS = np.array([(0, 0), (1, 0), (-1, 0), (0, 1), (0, -1)])
torch.set_num_threads(8)


class Speaker(nn.Module):        # target one-hot -> symbol
    def __init__(self):
        super().__init__(); self.net = nn.Sequential(nn.Linear(G, H), nn.Tanh())
        self.pi = nn.Linear(H, K); self.v = nn.Linear(H, 1)

    def forward(self, x):
        z = self.net(x); return self.pi(z), self.v(z).squeeze(-1)


class Listener(nn.Module):       # symbol + goal-dirs -> move
    def __init__(self):
        super().__init__(); self.net = nn.Sequential(nn.Linear(LOBS, H), nn.Tanh(), nn.Linear(H, H), nn.Tanh())
        self.pi = nn.Linear(H, 5); self.v = nn.Linear(H, 1)

    def forward(self, x):
        z = self.net(x); return self.pi(z), self.v(z).squeeze(-1)


class Vec:
    def __init__(self, B, seed):
        self.B = B; self.rng = np.random.default_rng(seed); self.reset()

    def reset(self):
        B = self.B
        self.goal = np.zeros((B, G, 2), int); self.pos = np.zeros((B, 2), int)
        for b in range(B):
            idx = self.rng.permutation(N * N)[:G + 1]
            self.goal[b, :, 0], self.goal[b, :, 1] = idx[:G] % N, idx[:G] // N
            self.pos[b] = (idx[G] % N, idx[G] // N)
        self.target = self.rng.integers(0, G, size=B)         # hidden target color
        self.done = np.zeros(B, bool); self.reached = -np.ones(B, int)
        return self.target

    def lobs(self, sym):                                       # sym: (B,) symbol ints
        B = self.B; out = np.zeros((B, LOBS), np.float32)
        out[np.arange(B), sym] = 1.0                           # symbol one-hot
        o = K
        for g in range(G):
            out[:, o:o + 2] = np.sign(self.goal[:, g, :] - self.pos); o += 2
        return out

    def step(self, move):
        B = self.B; active = ~self.done
        newp = self.pos + DIRS[move]
        inb = (newp[:, 0] >= 0) & (newp[:, 0] < N) & (newp[:, 1] >= 0) & (newp[:, 1] < N)
        self.pos = np.where((inb & active)[:, None], newp, self.pos)
        rew = np.where(active, -0.02, 0.0); term = np.zeros(B, bool)
        for b in np.where(active)[0]:
            for g in range(G):
                if tuple(self.pos[b]) == tuple(self.goal[b, g]):
                    rew[b] += 1.0 if g == self.target[b] else -0.3
                    self.done[b] = True; term[b] = True; self.reached[b] = g
                    break
        return rew.astype(np.float32), term, active


def gae(R, V, TERM, ACT, gamma=0.99, lam=0.95):
    B, Tt = R.shape; adv = np.zeros((B, Tt), np.float32); last = np.zeros(B, np.float32)
    for t in reversed(range(Tt)):
        nextv = V[:, t + 1] if t + 1 < Tt else np.zeros(B, np.float32)
        nonterm = 1.0 - TERM[:, t]
        delta = R[:, t] + gamma * nextv * nonterm - V[:, t]
        last = delta + gamma * lam * nonterm * last; adv[:, t] = last * ACT[:, t]
    return adv, adv + V


def episode(spk, lis, env, ablate=False, greedy=False):
    B = env.B; tgt = env.reset()
    th = torch.from_numpy(np.eye(G, dtype=np.float32)[tgt])
    slog, sval = spk(th); sdist = torch.distributions.Categorical(logits=slog)
    sym = (slog.argmax(1) if greedy else sdist.sample())
    slp = sdist.log_prob(sym); symn = sym.numpy()
    fed = env.rng.integers(0, K, size=B) if ablate else symn   # scramble channel for the debunk
    O = np.zeros((B, MAXT, LOBS), np.float32); A = np.zeros((B, MAXT), int)
    LP = np.zeros((B, MAXT), np.float32); V = np.zeros((B, MAXT), np.float32)
    R = np.zeros((B, MAXT), np.float32); TERM = np.zeros((B, MAXT), np.float32); ACT = np.zeros((B, MAXT), np.float32)
    for t in range(MAXT):
        obs = env.lobs(fed); ot = torch.from_numpy(obs)
        with torch.no_grad():
            llog, lval = lis(ot)
        ldist = torch.distributions.Categorical(logits=llog)
        a = llog.argmax(1) if greedy else ldist.sample()
        O[:, t] = obs; A[:, t] = a.numpy(); LP[:, t] = ldist.log_prob(a).numpy(); V[:, t] = lval.numpy()
        ACT[:, t] = (~env.done).astype(np.float32)
        r, term, _ = env.step(a.numpy()); R[:, t] = r; TERM[:, t] = term
    success = float(((env.reached == env.target)).mean())
    return dict(th=th, sym=sym, slp=slp.detach(), sval=sval.detach(), R=R, O=O, A=A, LP=LP, V=V,
                TERM=TERM, ACT=ACT, success=success)


def train(iters=400, B=256, seed=0):
    torch.manual_seed(seed)
    spk, lis = Speaker(), Listener()
    opt = torch.optim.Adam(list(spk.parameters()) + list(lis.parameters()), lr=3e-3)
    env = Vec(B, seed)
    for it in range(iters):
        env = Vec(B, 1000 + seed * 9999 + it)
        ep = episode(spk, lis, env)
        Rtot = ep["R"].sum(1)                                  # episode return per env (shared)
        adv, ret = gae(ep["R"], ep["V"], ep["TERM"], ep["ACT"])
        m = ep["ACT"].sum()
        advn = (adv - (adv * ep["ACT"]).sum() / m) / (np.sqrt(((adv - (adv * ep["ACT"]).sum() / m) ** 2 * ep["ACT"]).sum() / m) + 1e-8)
        Ot, At, LPt = torch.from_numpy(ep["O"]), torch.from_numpy(ep["A"]), torch.from_numpy(ep["LP"])
        advt, rett, actt = torch.from_numpy(advn), torch.from_numpy(ret), torch.from_numpy(ep["ACT"])
        sadv = torch.from_numpy(Rtot.astype(np.float32)) - ep["sval"]
        sadvn = (sadv - sadv.mean()) / (sadv.std() + 1e-8)
        for _ in range(4):
            # listener
            llog, lval = lis(Ot.reshape(-1, LOBS)); llog = llog.reshape(B, MAXT, 5); lval = lval.reshape(B, MAXT)
            ld = torch.distributions.Categorical(logits=llog)
            lratio = torch.exp(ld.log_prob(At) - LPt)
            l1 = lratio * advt; l2 = torch.clamp(lratio, 0.8, 1.2) * advt
            lloss = -(torch.min(l1, l2) * actt).sum() / actt.sum() + 0.5 * (((lval - rett) ** 2) * actt).sum() / actt.sum() - 0.01 * (ld.entropy() * actt).sum() / actt.sum()
            # speaker (bandit: one symbol per episode)
            slog, sval = spk(ep["th"]); sd = torch.distributions.Categorical(logits=slog)
            sratio = torch.exp(sd.log_prob(ep["sym"]) - ep["slp"])
            s1 = sratio * sadvn; s2 = torch.clamp(sratio, 0.8, 1.2) * sadvn
            sloss = -torch.min(s1, s2).mean() + 0.5 * ((sval - torch.from_numpy(Rtot.astype(np.float32))) ** 2).mean() - 0.02 * sd.entropy().mean()
            opt.zero_grad(); (lloss + sloss).backward(); nn.utils.clip_grad_norm_(list(spk.parameters()) + list(lis.parameters()), 0.5); opt.step()
        if it % 80 == 0 or it == iters - 1:
            ev = np.mean([episode(spk, lis, Vec(400, 7000 + it + k), greedy=True)["success"] for k in range(2)])
            print(f"  iter {it:>3}: success {ev:.2f}", flush=True)
    return spk, lis


if __name__ == "__main__":
    a = sys.argv
    iters = int(a[a.index("--iters") + 1]) if "--iters" in a else 400
    here = __file__.rsplit("/", 1)[0]
    print(f"PHASE D1 emergent-communication (referential signaling): G={G} colors, K={K} symbols, chance={1/G:.2f}", flush=True)
    spk, lis = train(iters)
    normal = np.mean([episode(spk, lis, Vec(500, 9000 + k), greedy=True)["success"] for k in range(3)])
    ablate = np.mean([episode(spk, lis, Vec(500, 9100 + k), ablate=True, greedy=True)["success"] for k in range(3)])
    print(f"\nFINAL: success {normal:.2f}  |  channel-SCRAMBLED success {ablate:.2f}  (chance {1/G:.2f})", flush=True)
    print("VERDICT:", "EMERGENT COMMUNICATION (success>=0.6 AND scrambled collapses toward chance) -> genuine multi-agent competency"
          if (normal >= 0.6 and ablate <= 0.45) else f"no clean comms (normal {normal:.2f}, scrambled {ablate:.2f})", flush=True)
    torch.save({"spk": spk.state_dict(), "lis": lis.state_dict()}, f"{here}/commworld_net.pt")
    json.dump({"normal": round(float(normal), 3), "scrambled": round(float(ablate), 3), "chance": round(1/G, 3)},
              open(f"{here}/commworld_result.json", "w"))
