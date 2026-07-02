"""OPEN-ENDED N-agent hunt -- RUNG 7 (a DIFFERENT competency KIND: collective AGREEMENT under conflict).

Rungs 0-6 were all spatial/temporal COORDINATION. Rung 7 tests a different kind -- social choice / consensus
under CONFLICTING preferences -- to check whether it, too, reduces to the closed media set or introduces a new
primitive. 3 agents, M shared options. Each agent has a PRIVATE preferred option (random). Each of T rounds all
pick an option; reward +1 to ALL iff all three pick the SAME option (consensus), plus a small +0.3 private
bonus to an agent if the consensus option is its own preference (so each prefers its own to win, but consensus
matters more -- a 3-way Battle-of-the-Sexes). Media: mutual last-pick OBSERVATION + a signal CHANNEL.

Fingerprint:
  - BLIND (zero others' last picks) / MUTE (scramble others' signals).
  - Expected: converge on a FOCAL option (cheapest; both ablations invariant) -- consensus reduces to a focal
    point, like rung 0. If instead it rides observation (join the plurality) or the channel, that still sits IN
    the closed set. A NEW medium would be something outside {focal-point, observation, communication}.

Run on VPS epc-venv: python openhunt9_ppo.py [--iters 1000]
"""
import numpy as np, sys, json, torch, torch.nn as nn

N_AG = 3; M = 3; K = 4; T = 12; H = 96
ODIM = M + 1 + (N_AG - 1) * M + (N_AG - 1) * K   # my-pref one-hot(M), round-frac(1), others' last picks(M each), others' last signals(K each)
torch.set_num_threads(8)


class Pol(nn.Module):
    def __init__(self):
        super().__init__(); self.gru = nn.GRU(ODIM, H, batch_first=True)
        self.pick = nn.Linear(H, M); self.sg = nn.Linear(H, K); self.v = nn.Linear(H, 1)

    def forward(self, x, h=None):
        z, h = self.gru(x, h); return self.pick(z), self.sg(z), self.v(z).squeeze(-1), h


class VecCons:
    def __init__(self, B, seed, blind=False, mute=False):
        self.B, self.blind, self.mute = B, blind, mute
        self.rng = np.random.default_rng(seed); self.reset()

    def reset(self):
        B = self.B
        self.pref = self.rng.integers(0, M, size=(B, N_AG))          # private preferences
        self.lpick = np.full((B, N_AG), -1, int); self.lsig = np.zeros((B, N_AG), int)
        self.t = 0; self.cons = np.zeros(B, int)
        return self.obs()

    def obs(self):
        B = self.B; out = np.zeros((B, N_AG, ODIM), np.float32)
        for i in range(N_AG):
            o = 0; out[np.arange(B), i, o + self.pref[:, i]] = 1.0; o += M
            out[:, i, o] = self.t / T; o += 1
            for j in range(N_AG):
                if j == i: continue
                if not self.blind:
                    pk = self.lpick[:, j]; valid = pk >= 0
                    out[np.arange(B)[valid], i, o + pk[valid]] = 1.0
                o += M
            for j in range(N_AG):
                if j == i: continue
                s = self.rng.integers(0, K, size=B) if self.mute else self.lsig[:, j]
                out[np.arange(B), i, o + s] = 1.0; o += K
        return out

    def step(self, pick, sg):
        B = self.B; rew = np.zeros(B, np.float32)
        same = (pick[:, 0] == pick[:, 1]) & (pick[:, 1] == pick[:, 2])
        rew += same.astype(np.float32)
        for i in range(N_AG):
            rew += 0.3 * (same & (pick[:, i] == self.pref[:, i])).astype(np.float32) / N_AG   # small private bonus
        self.cons += same.astype(int); self.lpick = pick.copy(); self.lsig = sg.copy(); self.t += 1
        return rew

    def consensus(self): return float(self.cons.mean())


def gae(R, V, gamma=0.9, lam=0.95):
    B, Tt = R.shape; adv = np.zeros((B, Tt), np.float32); last = np.zeros(B, np.float32)
    for t in reversed(range(Tt)):
        nextv = V[:, t+1] if t+1 < Tt else np.zeros(B, np.float32)
        delta = R[:, t] + gamma*nextv - V[:, t]; last = delta + gamma*lam*last; adv[:, t] = last
    return adv, adv + V


def rollout(net, env, greedy=False):
    B, A = env.B, N_AG; obs = env.obs(); h = [None]*A
    O = np.zeros((B, A, T, ODIM), np.float32); PK = np.zeros((B, A, T), int); SG = np.zeros((B, A, T), int)
    LP = np.zeros((B, A, T), np.float32); V = np.zeros((B, A, T), np.float32); R = np.zeros((B, T), np.float32)
    for t in range(T):
        pk = np.zeros((B, A), int); sg = np.zeros((B, A), int)
        for i in range(A):
            with torch.no_grad():
                pl, sl, v, h[i] = net(torch.from_numpy(obs[:, i])[:, None, :], h[i])
            pd = torch.distributions.Categorical(logits=pl[:, 0]); sd = torch.distributions.Categorical(logits=sl[:, 0])
            p = pl[:, 0].argmax(1) if greedy else pd.sample(); s = sl[:, 0].argmax(1) if greedy else sd.sample()
            O[:, i, t] = obs[:, i]; PK[:, i, t] = p.numpy(); SG[:, i, t] = s.numpy()
            LP[:, i, t] = (pd.log_prob(p) + sd.log_prob(s)).numpy(); V[:, i, t] = v[:, 0].numpy()
            pk[:, i] = p.numpy(); sg[:, i] = s.numpy()
        R[:, t] = env.step(pk, sg); obs = env.obs()
    return O, PK, SG, LP, V, R


def train(iters, seed=0, B=256):
    torch.manual_seed(seed); net = Pol(); opt = torch.optim.Adam(net.parameters(), lr=3e-3)
    for it in range(iters):
        env = VecCons(B, 1000 + it)
        O, PK, SG, LP, V, R = rollout(net, env)
        Rag = np.repeat(R[:, None, :], N_AG, axis=1)
        O2 = O.reshape(B*N_AG, T, ODIM); PK2, SG2 = PK.reshape(B*N_AG, T), SG.reshape(B*N_AG, T)
        LP2, V2, R2 = LP.reshape(B*N_AG, T), V.reshape(B*N_AG, T), Rag.reshape(B*N_AG, T)
        adv, ret = gae(R2, V2); adv = (adv - adv.mean()) / (adv.std() + 1e-8)
        Ot, PKt, SGt, LPt = torch.from_numpy(O2), torch.from_numpy(PK2), torch.from_numpy(SG2), torch.from_numpy(LP2)
        advt, rett = torch.from_numpy(adv), torch.from_numpy(ret)
        for _ in range(4):
            pl, sl, v, _ = net(Ot)
            pd = torch.distributions.Categorical(logits=pl); sd = torch.distributions.Categorical(logits=sl)
            ratio = torch.exp(pd.log_prob(PKt) + sd.log_prob(SGt) - LPt)
            s1 = ratio*advt; s2 = torch.clamp(ratio, 0.8, 1.2)*advt
            loss = -torch.min(s1, s2).mean() + 0.5*((v-rett)**2).mean() - 0.02*(pd.entropy()+sd.entropy()).mean()
            opt.zero_grad(); loss.backward(); nn.utils.clip_grad_norm_(net.parameters(), 0.5); opt.step()
        if it % 100 == 0 or it == iters - 1:
            print(f"  iter {it:>3}: consensus-rounds {_cons(net, 7000+it):.2f}", flush=True)
    return net


def _cons(net, seed, B=400, **kw):
    env = VecCons(B, seed, **kw); rollout(net, env, greedy=False); return env.consensus()


if __name__ == "__main__":
    a = sys.argv
    iters = int(a[a.index("--iters")+1]) if "--iters" in a else 1000
    here = __file__.rsplit("/", 1)[0]
    print(f"OPEN-ENDED HUNT rung 7: 3-agent CONSENSUS under conflicting preferences ({M} options, {T} rounds)", flush=True)
    net = train(iters)
    norm = np.mean([_cons(net, 9000+k) for k in range(3)])
    blind = np.mean([_cons(net, 9200+k, blind=True) for k in range(3)])
    mute = np.mean([_cons(net, 9300+k, mute=True) for k in range(3)])
    print(f"\nFINGERPRINT (consensus-rounds/{T}): normal {norm:.2f} | blind {blind:.2f} | mute {mute:.2f}", flush=True)
    ob = norm - blind; ch = norm - mute
    print(f"  needs-OBSERVATION = {ob:+.2f}; needs-CHANNEL = {ch:+.2f}", flush=True)
    mech = ("focal point (both ablations invariant)" if ob < 0.8 and ch < 0.8 else
            "observation (join the plurality)" if ob >= ch else "channel (announce/negotiate)")
    print(f"RUNG-7 MECHANISM: {mech} -- {'IN the closed set' if True else ''} (consensus reduces to a known medium)", flush=True)
    torch.save(net.state_dict(), f"{here}/openhunt9_net.pt")
    json.dump({"normal": round(float(norm),2), "blind": round(float(blind),2), "mute": round(float(mute),2),
               "mechanism": mech}, open(f"{here}/openhunt9_result.json","w"))
