"""OPEN-ENDED N-agent hunt (SCALE, new-to-science attempt) -- RUNG 0 (substrate + fingerprint).

Discovery engine: build a rich multi-agent world -> FINGERPRINT which coordination mechanism the agents
actually use (full ablation battery) -> (later rungs) REMOVE that mechanism and re-run, climbing the cascade.
A rung whose emergent mechanism resists literature-naming is the new-to-science candidate.

RUNG 0 world: N_AG agents forage HEAVY food -- a food item is collected only when >=2 agents occupy its cell
SIMULTANEOUSLY (a solo agent gets nothing -> genuine cooperation is required). Collected food gives a SHARED
team reward and respawns. Agents have BOTH an observation channel (they see others' positions) AND a signal
channel (K symbols). Fingerprint the rung-0 mechanism with ablations:
  - SOLO (1 agent) -> ~0 food: confirms cooperation is genuinely required (heavy food).
  - BLIND (zero others' positions) collapse => coordination rides OBSERVATION (mutual-position rendezvous).
  - MUTE (scramble others' signals) collapse => coordination rides the CHANNEL.
  - MEMWIPE (reset hidden each step) collapse => coordination needs recurrent memory.
Whatever survives / collapses names the rung-0 mechanism (expected: mutual-observation rendezvous = KNOWN).

Run on VPS epc-venv: python openhunt_ppo.py [--iters 800]
"""
import numpy as np, sys, json, torch, torch.nn as nn

N = 9; N_AG = 3; F = 2; K = 3; T = 30; H = 96; NMOVE = 5
ODIM = 2 + F * 2 + (N_AG - 1) * 2 + (N_AG - 1) * K       # own xy, each food rel, others rel, others last-signal
DIRS = np.array([(0, 0), (1, 0), (-1, 0), (0, 1), (0, -1)])
torch.set_num_threads(8)


class Pol(nn.Module):
    def __init__(self):
        super().__init__(); self.gru = nn.GRU(ODIM, H, batch_first=True)
        self.mv = nn.Linear(H, NMOVE); self.sg = nn.Linear(H, K); self.v = nn.Linear(H, 1)

    def forward(self, x, h=None):
        z, h = self.gru(x, h); return self.mv(z), self.sg(z), self.v(z).squeeze(-1), h


def clip(p): return np.clip(p, 0, N - 1)


class VecHunt:
    def __init__(self, B, seed, n_ag=N_AG, blind=False, mute=False):
        self.B, self.n_ag, self.blind, self.mute = B, n_ag, blind, mute
        self.rng = np.random.default_rng(seed); self.reset()

    def reset(self):
        B, A = self.B, self.n_ag
        self.ap = self.rng.integers(0, N, size=(B, A, 2))
        self.fp = self.rng.integers(0, N, size=(B, F, 2))
        self.lsig = np.zeros((B, A), int)
        self.collected = np.zeros(B, int)
        return self.obs()

    def obs(self):
        B, A = self.B, self.n_ag; out = np.zeros((B, A, ODIM), np.float32)
        for i in range(A):
            o = 0; out[:, i, o:o+2] = self.ap[:, i] / (N - 1); o += 2
            for f in range(F): out[:, i, o:o+2] = (self.fp[:, f] - self.ap[:, i]) / (N - 1); o += 2
            for j in range(A):
                if j == i: continue
                rel = np.zeros((B, 2), np.float32) if self.blind else (self.ap[:, j] - self.ap[:, i]) / (N - 1)
                out[:, i, o:o+2] = rel; o += 2
            for j in range(A):
                if j == i: continue
                sj = self.rng.integers(0, K, size=B) if self.mute else self.lsig[:, j]
                out[np.arange(B), i, o + sj] = 1.0; o += K
        return out

    def step(self, mv, sg):
        B, A = self.B, self.n_ag; rew = np.full(B, -0.01, np.float32)
        for i in range(A): self.ap[:, i] = clip(self.ap[:, i] + DIRS[mv[:, i]])
        self.lsig = sg.copy()
        for b in range(B):
            for f in range(F):
                onit = sum(1 for i in range(A) if tuple(self.ap[b, i]) == tuple(self.fp[b, f]))
                if onit >= 2:                                     # HEAVY: needs >=2 agents on the cell
                    rew[b] += 1.0; self.collected[b] += 1
                    self.fp[b, f] = self.rng.integers(0, N, size=2)
        return rew

    def food(self): return float(self.collected.mean())


def gae(R, V, gamma=0.97, lam=0.95):
    B, Tt = R.shape; adv = np.zeros((B, Tt), np.float32); last = np.zeros(B, np.float32)
    for t in reversed(range(Tt)):
        nextv = V[:, t+1] if t+1 < Tt else np.zeros(B, np.float32)
        delta = R[:, t] + gamma*nextv - V[:, t]; last = delta + gamma*lam*last; adv[:, t] = last
    return adv, adv + V


def rollout(net, env, greedy=False):
    B, A = env.B, env.n_ag; obs = env.obs(); h = [None]*A
    O = np.zeros((B, A, T, ODIM), np.float32); MV = np.zeros((B, A, T), int); SG = np.zeros((B, A, T), int)
    LP = np.zeros((B, A, T), np.float32); V = np.zeros((B, A, T), np.float32); R = np.zeros((B, T), np.float32)
    for t in range(T):
        mv = np.zeros((B, A), int); sg = np.zeros((B, A), int)
        for i in range(A):
            with torch.no_grad():
                ml, sl, v, h[i] = net(torch.from_numpy(obs[:, i])[:, None, :], h[i])
            md = torch.distributions.Categorical(logits=ml[:, 0]); sd = torch.distributions.Categorical(logits=sl[:, 0])
            m = ml[:, 0].argmax(1) if greedy else md.sample(); s = sl[:, 0].argmax(1) if greedy else sd.sample()
            O[:, i, t] = obs[:, i]; MV[:, i, t] = m.numpy(); SG[:, i, t] = s.numpy()
            LP[:, i, t] = (md.log_prob(m) + sd.log_prob(s)).numpy(); V[:, i, t] = v[:, 0].numpy()
            mv[:, i] = m.numpy(); sg[:, i] = s.numpy()
        R[:, t] = env.step(mv, sg); obs = env.obs()
    return O, MV, SG, LP, V, R


def train(iters, seed=0, B=192):
    torch.manual_seed(seed); net = Pol(); opt = torch.optim.Adam(net.parameters(), lr=3e-3)
    for it in range(iters):
        env = VecHunt(B, 1000 + it)
        O, MV, SG, LP, V, R = rollout(net, env)
        A = env.n_ag; Rag = np.repeat(R[:, None, :], A, axis=1)          # shared team reward -> every agent
        O2 = O.reshape(B*A, T, ODIM); MV2, SG2 = MV.reshape(B*A, T), SG.reshape(B*A, T)
        LP2, V2, R2 = LP.reshape(B*A, T), V.reshape(B*A, T), Rag.reshape(B*A, T)
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
        if it % 100 == 0 or it == iters - 1:
            print(f"  iter {it:>3}: food/episode {_food(net, 7000+it):.2f}", flush=True)
    return net


def _food(net, seed, B=400, **kw):
    env = VecHunt(B, seed, **kw); rollout(net, env, greedy=False); return env.food()


if __name__ == "__main__":
    a = sys.argv
    iters = int(a[a.index("--iters")+1]) if "--iters" in a else 800
    here = __file__.rsplit("/", 1)[0]
    print(f"OPEN-ENDED HUNT rung 0: {N_AG} agents, HEAVY food (needs 2 on a cell), obs+signal channels", flush=True)
    net = train(iters)
    norm = np.mean([_food(net, 9000+k) for k in range(3)])
    solo = np.mean([_food(net, 9100+k, n_ag=1) for k in range(3)])
    blind = np.mean([_food(net, 9200+k, blind=True) for k in range(3)])
    mute = np.mean([_food(net, 9300+k, mute=True) for k in range(3)])
    print(f"\nFINGERPRINT (food/episode): normal {norm:.2f} | solo(1 agent) {solo:.2f} | blind {blind:.2f} | mute {mute:.2f}", flush=True)
    uses_obs = norm - blind; uses_sig = norm - mute; needs_coop = norm - solo
    print(f"  needs-cooperation (heavy) = {needs_coop:+.2f}; rides-OBSERVATION = {uses_obs:+.2f}; rides-CHANNEL = {uses_sig:+.2f}", flush=True)
    mech = "observation-rendezvous" if uses_obs >= 0.5 and uses_obs > uses_sig else ("channel" if uses_sig >= 0.5 else "unclear")
    print(f"RUNG-0 MECHANISM: {mech} (expected known: mutual-position rendezvous). Next rung: remove it, re-run.", flush=True)
    torch.save(net.state_dict(), f"{here}/openhunt_net.pt")
    json.dump({"normal": round(float(norm),2), "solo": round(float(solo),2), "blind": round(float(blind),2),
               "mute": round(float(mute),2), "mechanism": mech}, open(f"{here}/openhunt_result.json","w"))
