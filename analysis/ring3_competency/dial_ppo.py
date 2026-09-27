"""Cost-ladder DIAL experiment (paper §8: cost-hierarchy / communication-niche, causal upgrade).

The coordination cascade (ring3_coordination_cascade.md) established the communication niche across
8 SEPARATE worlds. Track-2 asks for the single missing "causal-ladder experiment": ONE environment
with a continuously dialable information structure, showing the signal channel is AVAILABLE THE WHOLE
TIME but DECLINED under symmetry and recruited exactly as information asymmetry increases.

Design (extends commhunt_ppo.py, same proven joint-PPO base):
  - 2 agents on an N x N grid, G goals, one is the hidden TARGET. Referential: team reward iff the
    MOVER (agent 1) ends on the target.
  - agent 0 (SPEAKER) always sees the target. agent 1 (MOVER) sees the target with prob (1 - alpha),
    per episode (mask `sees1`). alpha is the ASYMMETRY DIAL: alpha=0 -> fully symmetric (both always
    know the target, a shared referent / focal point); alpha=1 -> fully asymmetric (only the speaker
    knows). A K-symbol signal channel (speaker's last symbol, seen by the mover) is present at EVERY
    alpha. The speaker is never rewarded for reaching the target itself, so its position is
    uninformative -> the channel is the only medium for the asymmetric case.
  - Ablations (the fingerprint battery, run per alpha): normal; channel-SCRAMBLE (mover gets a random
    symbol) -> collapse => the policy USES the channel; blind-PARTNER (zero the partner rel-pos) ->
    tests whether partner OBSERVATION is used instead.

Predicted (communication niche): channel-reliance = success(normal) - success(scramble) is ~0 at
alpha=0 (channel present but DECLINED under symmetry) and rises monotonically toward alpha=1; the
speaker->target mutual information I(symbol; target) rises from ~0 in lockstep. That is the
available-but-declined-channel-under-symmetry result, as a single causal dial.

Run (VPS epc-venv, ~8 min/config):
  PYTHONDONTWRITEBYTECODE=1 /home/matthewhmaxwell/epc-venv/bin/python dial_ppo.py --iters 500 \
      --alphas 0,0.25,0.5,0.75,1 --seeds 2
"""
import numpy as np, sys, json, os, torch, torch.nn as nn

N = 5; G = 3; K = 4; MAXT = 18; H = 64
ODIM = G * 2 + 2 + K + 1 + G     # goal-dirs, other rel-dir, other last-symbol, agent-id, (visible) target
DIRS = np.array([(0, 0), (1, 0), (-1, 0), (0, 1), (0, -1)])
torch.set_num_threads(int(os.environ.get("DIAL_THREADS", "4")))  # polite: share the box with other jobs


class Agent(nn.Module):
    def __init__(self):
        super().__init__()
        self.net = nn.Sequential(nn.Linear(ODIM, H), nn.Tanh(), nn.Linear(H, H), nn.Tanh())
        self.move = nn.Linear(H, 5); self.sym = nn.Linear(H, K); self.v = nn.Linear(H, 1)

    def forward(self, x):
        z = self.net(x); return self.move(z), self.sym(z), self.v(z).squeeze(-1)


class DialEnv:
    """Referential rendezvous with an ASYMMETRY dial alpha. agent0=speaker (always sees target),
    agent1=mover (sees target iff a per-episode Bernoulli(1-alpha) coin)."""
    def __init__(self, B, alpha, seed):
        self.B, self.alpha = B, alpha; self.rng = np.random.default_rng(seed); self.reset()

    def reset(self):
        B = self.B
        self.goal = np.zeros((B, G, 2), int); self.p = np.zeros((B, 2, 2), int)
        for b in range(B):
            idx = self.rng.permutation(N * N)[:G + 2]
            self.goal[b, :, 0], self.goal[b, :, 1] = idx[:G] % N, idx[:G] // N
            self.p[b, 0] = (idx[G] % N, idx[G] // N); self.p[b, 1] = (idx[G + 1] % N, idx[G + 1] // N)
        self.tgt = self.rng.integers(0, G, size=B)
        self.sees1 = self.rng.random(B) >= self.alpha         # mover sees target iff True
        self.last_sym = np.zeros((B, 2), int)
        return None

    def on_goal(self, i):
        eq = (self.goal == self.p[:, i, None, :]).all(2)
        return np.where(eq.any(1), eq.argmax(1), -1)

    def obs(self, i, scramble=False, blind=False):
        B = self.B; out = np.zeros((B, ODIM), np.float32); o = 0
        for g in range(G):
            out[:, o:o + 2] = np.sign(self.goal[:, g, :] - self.p[:, i, :]); o += 2
        oth = 1 - i
        rel = np.sign(self.p[:, oth, :] - self.p[:, i, :])
        out[:, o:o + 2] = 0.0 if blind else rel; o += 2
        sym = self.rng.integers(0, K, size=B) if scramble else self.last_sym[:, oth]
        out[np.arange(B), o + sym] = 1.0; o += K            # partner's last symbol (present at all alpha)
        out[:, o] = i; o += 1                                # agent id
        if i == 0:                                           # speaker always sees target
            out[np.arange(B), o + self.tgt] = 1.0
        else:                                                # mover sees target only where sees1
            vis = self.sees1
            out[np.arange(B)[vis], o + self.tgt[vis]] = 1.0
        o += G
        return out

    def step(self, mv):
        B = self.B
        for i in range(2):
            np_ = self.p[:, i, :] + DIRS[mv[:, i]]
            inb = (np_[:, 0] >= 0) & (np_[:, 0] < N) & (np_[:, 1] >= 0) & (np_[:, 1] < N)
            self.p[:, i, :] = np.where(inb[:, None], np_, self.p[:, i, :])
        r = np.full((B, 2), -0.01, np.float32)
        r[:, 1] += 0.03 * (self.on_goal(1) == self.tgt)      # shape the mover toward the TARGET goal (env-side signal, not leaked to obs)
        return r

    def terminal(self):
        win = (self.on_goal(1) == self.tgt).astype(np.float32)
        r = np.zeros((self.B, 2), np.float32); r[:, 0] += win; r[:, 1] += win
        return r

    def success(self, mask=None):
        ok = (self.on_goal(1) == self.tgt)
        return float(ok.mean()) if mask is None else float(ok[mask].mean()) if mask.any() else float("nan")


def gae(R, V, gamma=0.99, lam=0.95):
    A, Tt = R.shape; adv = np.zeros((A, Tt), np.float32); last = np.zeros(A, np.float32)
    for t in reversed(range(Tt)):
        nextv = V[:, t + 1] if t + 1 < Tt else np.zeros(A, np.float32)
        delta = R[:, t] + gamma * nextv - V[:, t]
        last = delta + gamma * lam * last; adv[:, t] = last
    return adv, adv + V


def rollout(net, env, greedy=False, scramble=False, blind=False, log_sym=False):
    B = env.B; env.reset()
    O = np.zeros((B, 2, MAXT, ODIM), np.float32); MV = np.zeros((B, 2, MAXT), int); SY = np.zeros((B, 2, MAXT), int)
    LP = np.zeros((B, 2, MAXT), np.float32); V = np.zeros((B, 2, MAXT), np.float32); R = np.zeros((B, 2, MAXT), np.float32)
    for t in range(MAXT):
        for i in range(2):
            ob = env.obs(i, scramble=scramble, blind=blind); O[:, i, t] = ob
            with torch.no_grad():
                ml, sl, v = net(torch.from_numpy(ob))
            md, sd = torch.distributions.Categorical(logits=ml), torch.distributions.Categorical(logits=sl)
            mv = ml.argmax(1) if greedy else md.sample()
            sy = sl.argmax(1) if greedy else sd.sample()
            MV[:, i, t] = mv.numpy(); SY[:, i, t] = sy.numpy()
            LP[:, i, t] = (md.log_prob(mv) + sd.log_prob(sy)).numpy(); V[:, i, t] = v.numpy()
        r = env.step(np.stack([MV[:, 0, t], MV[:, 1, t]], 1))
        env.last_sym = np.stack([SY[:, 0, t], SY[:, 1, t]], 1)
        R[:, :, t] = r
    R[:, :, -1] += env.terminal()
    if log_sym:
        return O, MV, SY, LP, V, R, SY[:, 0, 0].copy(), env.tgt.copy(), env.sees1.copy()  # speaker 1st symbol
    return O, MV, SY, LP, V, R


def train(alpha, iters, seed=0, B=320):
    torch.manual_seed(seed); net = Agent(); opt = torch.optim.Adam(net.parameters(), lr=3e-3)
    for it in range(iters):
        env = DialEnv(B, alpha, 1000 + seed * 99999 + it)
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
            nlp = md.log_prob(MVt) + sd.log_prob(SYt); ratio = torch.exp(nlp - LPt)
            s1 = ratio * advt; s2 = torch.clamp(ratio, 0.8, 1.2) * advt
            loss = -torch.min(s1, s2).mean() + 0.5 * ((v - rett) ** 2).mean() - 0.01 * (md.entropy() + sd.entropy()).mean()
            opt.zero_grad(); loss.backward(); nn.utils.clip_grad_norm_(net.parameters(), 0.5); opt.step()
        if it % 100 == 0 or it == iters - 1:
            print(f"  [alpha={alpha:.2f} seed{seed}] iter {it:>3}: success {_succ(net, alpha, 7000 + it):.2f}", flush=True)
    return net


def _succ(net, alpha, seed, scramble=False, blind=False, B=600, subset=None):
    env = DialEnv(B, alpha, seed); _ = rollout(net, env, greedy=True, scramble=scramble, blind=blind)
    if subset == "asym":  return env.success(~env.sees1)
    if subset == "sym":   return env.success(env.sees1)
    return env.success()


def _mi_sym_tgt(net, alpha, seed, B=1500):
    """Empirical mutual information I(speaker-first-symbol ; target) over ASYMMETRIC episodes (bits)."""
    env = DialEnv(B, alpha, seed); *_, sym0, tgt, sees1 = rollout(net, env, greedy=True, log_sym=True)
    m = ~sees1
    if m.sum() < 50: return float("nan")
    s, g = sym0[m], tgt[m]; joint = np.zeros((K, G))
    for a, b in zip(s, g): joint[a, b] += 1
    joint /= joint.sum(); ps = joint.sum(1, keepdims=True); pg = joint.sum(0, keepdims=True)
    nz = joint > 0
    return float((joint[nz] * np.log2(joint[nz] / (ps @ pg)[nz])).sum())


if __name__ == "__main__":
    a = sys.argv
    iters = int(a[a.index("--iters") + 1]) if "--iters" in a else 500
    alphas = [float(x) for x in (a[a.index("--alphas") + 1].split(",") if "--alphas" in a else ["0", "0.25", "0.5", "0.75", "1"])]
    seeds = int(a[a.index("--seeds") + 1]) if "--seeds" in a else 2
    here = __file__.rsplit("/", 1)[0]; out = f"{here}/dial_sweep.json"
    prior = {f'{r["alpha"]:.2f}': r for r in json.load(open(out))} if os.path.exists(out) else {}
    print(f"COST-LADDER DIAL: alphas={alphas}, iters={iters}, seeds={seeds} (resuming {sorted(prior)})", flush=True)
    results = []
    for alpha in alphas:
        key = f"{alpha:.2f}"
        if key in prior:
            results.append(prior[key]); print(f"alpha={key}: (cached)", flush=True); continue
        norm = scr = bld = asym = sym = mi = 0.0
        for s in range(seeds):
            net = train(alpha, iters, seed=s)
            norm += np.mean([_succ(net, alpha, 9000 + k) for k in range(3)]) / seeds
            scr  += np.mean([_succ(net, alpha, 9100 + k, scramble=True) for k in range(3)]) / seeds
            bld  += np.mean([_succ(net, alpha, 9200 + k, blind=True) for k in range(3)]) / seeds
            asym += _succ(net, alpha, 9300, subset="asym") / seeds
            sym  += (_succ(net, alpha, 9300, subset="sym") if alpha < 1 else float("nan")) / seeds
            mi   += _mi_sym_tgt(net, alpha, 9400) / seeds
            torch.save(net.state_dict(), f"{here}/dial_a{key}_s{s}.pt")
        chan_reliance = norm - scr; obs_reliance = norm - bld
        row = {"alpha": alpha, "success": round(norm, 3), "scramble": round(scr, 3), "blind": round(bld, 3),
               "channel_reliance": round(chan_reliance, 3), "obs_reliance": round(obs_reliance, 3),
               "succ_asym": round(asym, 3), "succ_sym": round(sym, 3), "mi_sym_tgt_bits": round(mi, 3)}
        results.append(row); json.dump(results, open(out, "w"), indent=1)
        print(f"alpha={key}: succ {norm:.2f} | scramble {scr:.2f} | blind {bld:.2f} "
              f"| CHANNEL-RELIANCE {chan_reliance:+.2f} | obs-reliance {obs_reliance:+.2f} "
              f"| MI(sym;tgt) {mi:.2f} bits", flush=True)
    json.dump(results, open(out, "w"), indent=1)
    print("saved dial_sweep.json", flush=True)
