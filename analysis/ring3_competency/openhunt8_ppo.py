"""OPEN-ENDED N-agent hunt -- RUNG 6 (3+ agents: INDIRECT / TRANSITIVE coordination -- a NEW medium?).

The completeness claim (from the spatial+temporal cascades) is that all coordination reduces to a CLOSED set of
media {focal-point, observation, communication, social-memory}. 3+ agents is the regime that could break it,
because it unlocks INDIRECT coordination the 2-agent primitives cannot express. Cleanest instance: a RELAY.

World: 3 agents in a chain. SOURCE (agent 0) sees a target cell. SEEKER (agent 2) must reach it. RELAY (agent 1)
sits between: the source's signal is heard ONLY by the relay; the relay's signal ONLY by the seeker; there is NO
direct source->seeker link. The relay NEVER sees the target -- so to succeed it must FORWARD information it
cannot itself use (transitive communication). Seeker moves; source + relay are stationary.

Fingerprint (a genuine relay chain requires BOTH hops):
  - MUTE-A (scramble source->relay) -> collapse (relay gets noise).
  - MUTE-B (scramble relay->seeker) -> collapse (seeker gets noise).
Question: is relayed coordination a NEW primitive medium, or just COMMUNICATION composed across a hop (which
keeps the closed set intact -> completeness holds -> still 0 new-to-science)? Expected KNOWN (routing / gossip /
multi-hop signaling); the point is to test whether 3 agents produce anything OUTSIDE the closed set.

Run on VPS epc-venv: python openhunt8_ppo.py [--iters 1500]
"""
import numpy as np, sys, json, torch, torch.nn as nn

N = 7; N_AG = 3; K = 8; T = 22; H = 96; NMOVE = 5
ODIM = 3 + 2 + 2 + K              # role one-hot(3), own pos(2; seeker), target(2; source), predecessor's last signal(K)
DIRS = np.array([(0, 0), (1, 0), (-1, 0), (0, 1), (0, -1)])
torch.set_num_threads(8)


class Pol(nn.Module):
    def __init__(self):
        super().__init__(); self.gru = nn.GRU(ODIM, H, batch_first=True)
        self.mv = nn.Linear(H, NMOVE); self.sg = nn.Linear(H, K); self.v = nn.Linear(H, 1)

    def forward(self, x, h=None):
        z, h = self.gru(x, h); return self.mv(z), self.sg(z), self.v(z).squeeze(-1), h


class VecRelay:
    def __init__(self, B, seed, muteA=False, muteB=False):
        self.B, self.muteA, self.muteB = B, muteA, muteB
        self.rng = np.random.default_rng(seed); self.reset()

    def reset(self):
        B = self.B
        self.seek = self.rng.integers(0, N, size=(B, 2))
        self.tgt = self.rng.integers(0, N, size=(B, 2))
        self.lsig = np.zeros((B, N_AG), int); self.hit = np.zeros(B, int)
        return self.obs()

    def obs(self):
        B = self.B; out = np.zeros((B, N_AG, ODIM), np.float32)
        for i in range(N_AG):
            o = 0; out[:, i, o + i] = 1.0; o += 3                      # role
            if i == 2: out[:, i, o:o+2] = self.seek / (N - 1)          # seeker position
            o += 2
            if i == 0: out[:, i, o:o+2] = self.tgt / (N - 1)           # source sees the target
            o += 2
            # predecessor's last signal: agent i hears agent i-1 (chain A->B->C); agent 0 hears nothing
            if i == 0:
                pass
            else:
                s = self.lsig[:, i - 1].copy()
                if i == 1 and self.muteA: s = self.rng.integers(0, K, size=B)   # source->relay
                if i == 2 and self.muteB: s = self.rng.integers(0, K, size=B)   # relay->seeker
                out[np.arange(B), i, o + s] = 1.0
            o += K
        return out

    def step(self, mv, sg):
        B = self.B; rew = np.full(B, -0.01, np.float32)
        self.seek = np.clip(self.seek + DIRS[mv[:, 2]], 0, N - 1)      # only the seeker moves
        self.lsig = sg.copy()
        at = np.all(self.seek == self.tgt, axis=1)
        for b in range(B):
            if at[b]:
                rew[b] += 1.0; self.hit[b] += 1
                self.tgt[b] = self.rng.integers(0, N, size=2); self.seek[b] = self.rng.integers(0, N, size=2)
        return rew

    def hits(self): return float(self.hit.mean())


def gae(R, V, gamma=0.97, lam=0.95):
    B, Tt = R.shape; adv = np.zeros((B, Tt), np.float32); last = np.zeros(B, np.float32)
    for t in reversed(range(Tt)):
        nextv = V[:, t+1] if t+1 < Tt else np.zeros(B, np.float32)
        delta = R[:, t] + gamma*nextv - V[:, t]; last = delta + gamma*lam*last; adv[:, t] = last
    return adv, adv + V


def rollout(net, env, greedy=False):
    B, A = env.B, N_AG; obs = env.obs(); h = [None]*A
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
        env = VecRelay(B, 1000 + it)
        O, MV, SG, LP, V, R = rollout(net, env)
        Rag = np.repeat(R[:, None, :], N_AG, axis=1)
        O2 = O.reshape(B*N_AG, T, ODIM); MV2, SG2 = MV.reshape(B*N_AG, T), SG.reshape(B*N_AG, T)
        LP2, V2, R2 = LP.reshape(B*N_AG, T), V.reshape(B*N_AG, T), Rag.reshape(B*N_AG, T)
        adv, ret = gae(R2, V2); adv = (adv - adv.mean()) / (adv.std() + 1e-8)
        Ot, MVt, SGt, LPt = torch.from_numpy(O2), torch.from_numpy(MV2), torch.from_numpy(SG2), torch.from_numpy(LP2)
        advt, rett = torch.from_numpy(adv), torch.from_numpy(ret)
        for _ in range(4):
            ml, sl, v, _ = net(Ot)
            md = torch.distributions.Categorical(logits=ml); sd = torch.distributions.Categorical(logits=sl)
            ratio = torch.exp(md.log_prob(MVt) + sd.log_prob(SGt) - LPt)
            s1 = ratio*advt; s2 = torch.clamp(ratio, 0.8, 1.2)*advt
            loss = -torch.min(s1, s2).mean() + 0.5*((v-rett)**2).mean() - 0.02*(md.entropy()+sd.entropy()).mean()
            opt.zero_grad(); loss.backward(); nn.utils.clip_grad_norm_(net.parameters(), 0.5); opt.step()
        if it % 100 == 0 or it == iters - 1:
            print(f"  iter {it:>3}: hits/episode {_hits(net, 7000+it):.2f}", flush=True)
    return net


def _hits(net, seed, B=400, **kw):
    env = VecRelay(B, seed, **kw); rollout(net, env, greedy=False); return env.hits()


if __name__ == "__main__":
    a = sys.argv
    iters = int(a[a.index("--iters")+1]) if "--iters" in a else 1500
    here = __file__.rsplit("/", 1)[0]
    print(f"OPEN-ENDED HUNT rung 6: 3-agent RELAY chain (source->relay->seeker; relay never sees target)", flush=True)
    net = train(iters)
    norm = np.mean([_hits(net, 9000+k) for k in range(3)])
    mA = np.mean([_hits(net, 9300+k, muteA=True) for k in range(3)])
    mB = np.mean([_hits(net, 9400+k, muteB=True) for k in range(3)])
    print(f"\nFINGERPRINT (hits/episode): normal {norm:.2f} | mute-A(src->relay) {mA:.2f} | mute-B(relay->seeker) {mB:.2f}", flush=True)
    needs_A = norm - mA; needs_B = norm - mB
    relay = norm >= 1.0 and needs_A >= 0.4 and needs_B >= 0.4
    print(f"  needs-hop-A = {needs_A:+.2f}; needs-hop-B = {needs_B:+.2f}", flush=True)
    print(f"RUNG-6 MECHANISM: {'RELAY (transitive communication; BOTH hops required) = composed communication, IN the closed set' if relay else 'degenerate/weak'}", flush=True)
    torch.save(net.state_dict(), f"{here}/openhunt8_net.pt")
    json.dump({"normal": round(float(norm),2), "muteA": round(float(mA),2), "muteB": round(float(mB),2),
               "relay": bool(relay)}, open(f"{here}/openhunt8_result.json","w"))
