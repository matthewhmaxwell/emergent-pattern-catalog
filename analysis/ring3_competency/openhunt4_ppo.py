"""OPEN-ENDED N-agent hunt -- RUNG 3 (the cascade's top: introduce information ASYMMETRY -> is COMMUNICATION
FINALLY forced?).

Cascade: rung 0 used an environmental FOCAL POINT (food); rung 1 (focal point removed) -> OBSERVATION
(pursuit); rung 2 (observation removed) -> a COORDINATE FOCAL POINT (a fixed convention cell in absolute
coords) -- STILL no channel. The lesson: SYMMETRIC coordination never forces communication, because a shared
focal point is always cheaper. Communication should be forced ONLY by information ASYMMETRY.

RUNG 3 world: a private TARGET cell is shown to ONLY ONE agent (the SPEAKER); BOTH must meet AT the target;
the LISTENER sees its own position + the speaker's SIGNAL but NOT the target and NOT the speaker's position.
No focal point resolves the target (it is random each round), so the speaker MUST communicate it. Shared
policy, asymmetric input (the listener's target field is zeroed) -> the net must act as speaker or listener
by role. Reward +1 when both agents occupy the target; then a new random target + respawn.

Fingerprint:
  - MUTE (scramble the speaker's signal) -> COLLAPSE => the channel is now REQUIRED = communication FINALLY
    forced (contrast rungs 0-2 where mute was invariant). This is the cascade's top rung.
Expected KNOWN: emergent referential signaling of a location (Lewis game). The point is the CASCADE: comms is
forced iff information is asymmetric -> a coordination-media cost hierarchy with a communication niche.

Run on VPS epc-venv: python openhunt4_ppo.py [--iters 1500]
"""
import numpy as np, sys, json, torch, torch.nn as nn

N = 7; K = 8; T = 24; H = 96; NMOVE = 5
ODIM = 2 + 2 + 1 + K              # own pos(2), target(2; zero for listener), am-I-speaker flag(1), other signal(K)
DIRS = np.array([(0, 0), (1, 0), (-1, 0), (0, 1), (0, -1)])
torch.set_num_threads(8)


class Pol(nn.Module):
    def __init__(self):
        super().__init__(); self.gru = nn.GRU(ODIM, H, batch_first=True)
        self.mv = nn.Linear(H, NMOVE); self.sg = nn.Linear(H, K); self.v = nn.Linear(H, 1)

    def forward(self, x, h=None):
        z, h = self.gru(x, h); return self.mv(z), self.sg(z), self.v(z).squeeze(-1), h


class VecRef:
    def __init__(self, B, seed, mute=False):
        self.B, self.mute = B, mute; self.rng = np.random.default_rng(seed); self.reset()

    def reset(self):
        B = self.B
        self.ap = self.rng.integers(0, N, size=(B, 2, 2))
        self.tgt = self.rng.integers(0, N, size=(B, 2))
        self.speaker = self.rng.integers(0, 2, size=B)          # which agent sees the target
        self.lsig = np.zeros((B, 2), int); self.hit = np.zeros(B, int)
        return self.obs()

    def obs(self):
        B = self.B; out = np.zeros((B, 2, ODIM), np.float32)
        for i in range(2):
            j = 1 - i; o = 0
            out[:, i, o:o+2] = self.ap[:, i] / (N - 1); o += 2
            is_sp = (self.speaker == i)
            out[is_sp, i, o:o+2] = self.tgt[is_sp] / (N - 1); o += 2   # only the speaker sees the target
            out[:, i, o] = is_sp.astype(np.float32); o += 1
            sj = self.rng.integers(0, K, size=B) if self.mute else self.lsig[:, j]
            out[np.arange(B), i, o + sj] = 1.0
        return out

    def step(self, mv, sg):
        B = self.B; rew = np.full(B, -0.01, np.float32)
        for i in range(2): self.ap[:, i] = np.clip(self.ap[:, i] + DIRS[mv[:, i]], 0, N - 1)
        self.lsig = sg.copy()
        both = np.all(self.ap[:, 0] == self.tgt, axis=1) & np.all(self.ap[:, 1] == self.tgt, axis=1)
        for b in range(B):
            if both[b]:
                rew[b] += 1.0; self.hit[b] += 1
                self.tgt[b] = self.rng.integers(0, N, size=2); self.speaker[b] = self.rng.integers(0, 2)
                self.ap[b] = self.rng.integers(0, N, size=(2, 2))
        return rew

    def hits(self): return float(self.hit.mean())


def gae(R, V, gamma=0.97, lam=0.95):
    B, Tt = R.shape; adv = np.zeros((B, Tt), np.float32); last = np.zeros(B, np.float32)
    for t in reversed(range(Tt)):
        nextv = V[:, t+1] if t+1 < Tt else np.zeros(B, np.float32)
        delta = R[:, t] + gamma*nextv - V[:, t]; last = delta + gamma*lam*last; adv[:, t] = last
    return adv, adv + V


def rollout(net, env, greedy=False):
    B = env.B; obs = env.obs(); h = [None, None]
    O = np.zeros((B, 2, T, ODIM), np.float32); MV = np.zeros((B, 2, T), int); SG = np.zeros((B, 2, T), int)
    LP = np.zeros((B, 2, T), np.float32); V = np.zeros((B, 2, T), np.float32); R = np.zeros((B, T), np.float32)
    for t in range(T):
        mv = np.zeros((B, 2), int); sg = np.zeros((B, 2), int)
        for i in range(2):
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
        env = VecRef(B, 1000 + it)
        O, MV, SG, LP, V, R = rollout(net, env)
        Rag = np.repeat(R[:, None, :], 2, axis=1)
        O2 = O.reshape(B*2, T, ODIM); MV2, SG2 = MV.reshape(B*2, T), SG.reshape(B*2, T)
        LP2, V2, R2 = LP.reshape(B*2, T), V.reshape(B*2, T), Rag.reshape(B*2, T)
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
    env = VecRef(B, seed, **kw); rollout(net, env, greedy=False); return env.hits()


if __name__ == "__main__":
    a = sys.argv
    iters = int(a[a.index("--iters")+1]) if "--iters" in a else 1500
    here = __file__.rsplit("/", 1)[0]
    print(f"OPEN-ENDED HUNT rung 3: asymmetric-info rendezvous (only speaker sees target) -- is COMMUNICATION forced?", flush=True)
    net = train(iters)
    norm = np.mean([_hits(net, 9000+k) for k in range(3)])
    mute = np.mean([_hits(net, 9300+k, mute=True) for k in range(3)])
    print(f"\nFINGERPRINT (hits/episode): normal {norm:.2f} | mute {mute:.2f}", flush=True)
    uses_sig = norm - mute
    forced = uses_sig >= 0.4 and norm >= 1.0
    print(f"  rides-CHANNEL = {uses_sig:+.2f}", flush=True)
    print(f"RUNG-3 MECHANISM: {'COMMUNICATION FORCED (emergent referential location-signaling under asymmetry)' if forced else 'weak/none'}", flush=True)
    torch.save(net.state_dict(), f"{here}/openhunt4_net.pt")
    json.dump({"normal": round(float(norm),2), "mute": round(float(mute),2), "channel_forced": bool(forced)},
              open(f"{here}/openhunt4_result.json","w"))
