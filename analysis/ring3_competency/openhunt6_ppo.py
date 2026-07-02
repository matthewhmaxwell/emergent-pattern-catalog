"""OPEN-ENDED N-agent hunt -- RUNG 5 (change DOMAIN: TEMPORAL coordination -- does the cost hierarchy generalize
off the spatial axis, and is TIMING a NEW medium or does it reduce to known primitives?).

The spatial cascade (rungs 0-4) showed a coordination-media cost hierarchy. Rung 5 tests whether the SAME
hierarchy governs a different coordination domain -- TIME -- and whether temporal coordination needs a NEW
primitive medium. If it reduces to {temporal focal point (clock/count) < observation < communication}, that is
strong evidence the primitive-media set is CLOSED -- which would EXPLAIN the 0-new-to-science bound (every
multi-agent competency = minimal composition of a closed set of catalogued media; a NEW medium is what would be
required for new-to-science, and none appears).

World: 2 agents must FIRE SIMULTANEOUSLY. Each step an agent picks wait/fire (+ a signal). A hard REFRACTORY
(cooldown C after firing) forbids spamming -> they must agree a PHASE. Reward +1 when both fire the same step;
-0.1 when exactly one fires (wasted). Media available: a shared CLOCK (t/T), mutual last-action OBSERVATION,
and a signal CHANNEL. Fingerprint which one they use, then (later rungs) remove it and climb.
  - NO-CLOCK (zero the clock) -> if it survives, a recurrent agent COUNTS internally (temporal focal point via
    counting) rather than needing the explicit clock.
  - BLIND (zero other's last action); MUTE (scramble other's signal).
Expected: temporal FOCAL POINT (explicit clock or internal count) -- NOT a new medium => primitive set closed.

Run on VPS epc-venv: python openhunt6_ppo.py [--iters 1000]
"""
import numpy as np, sys, json, torch, torch.nn as nn

T = 24; C = 4; K = 4; H = 96                          # C = refractory cooldown after firing
ODIM = 1 + 1 + 1 + K                                  # clock, own cooldown, other's last action, other's last signal
torch.set_num_threads(8)


class Pol(nn.Module):
    def __init__(self):
        super().__init__(); self.gru = nn.GRU(ODIM, H, batch_first=True)
        self.act = nn.Linear(H, 2); self.sg = nn.Linear(H, K); self.v = nn.Linear(H, 1)

    def forward(self, x, h=None):
        z, h = self.gru(x, h); return self.act(z), self.sg(z), self.v(z).squeeze(-1), h


class VecSync:
    def __init__(self, B, seed, noclock=False, blind=False, mute=False):
        self.B, self.noclock, self.blind, self.mute = B, noclock, blind, mute
        self.rng = np.random.default_rng(seed); self.reset()

    def reset(self):
        B = self.B
        self.cd = np.zeros((B, 2), int); self.lact = np.zeros((B, 2), int); self.lsig = np.zeros((B, 2), int)
        self.t = 0; self.sync = np.zeros(B, int)
        return self.obs()

    def obs(self):
        B = self.B; out = np.zeros((B, 2, ODIM), np.float32)
        for i in range(2):
            j = 1 - i; o = 0
            out[:, i, o] = 0.0 if self.noclock else self.t / T; o += 1
            out[:, i, o] = self.cd[:, i] / C; o += 1
            out[:, i, o] = 0.0 if self.blind else self.lact[:, j]; o += 1
            sj = self.rng.integers(0, K, size=B) if self.mute else self.lsig[:, j]
            out[np.arange(B), i, o + sj] = 1.0
        return out

    def step(self, act, sg):
        B = self.B; rew = np.zeros(B, np.float32)
        can = (self.cd == 0)
        fire = (act == 1) & can                                    # can only fire when off cooldown
        both = fire[:, 0] & fire[:, 1]; one = fire[:, 0] ^ fire[:, 1]
        rew += 1.0 * both - 0.1 * one
        self.sync += both.astype(int)
        for i in range(2):
            self.cd[:, i] = np.where(fire[:, i], C, np.maximum(self.cd[:, i] - 1, 0))
        self.lact = fire.astype(int); self.lsig = sg.copy(); self.t += 1
        return rew

    def synced(self): return float(self.sync.mean())


def gae(R, V, gamma=0.95, lam=0.95):
    B, Tt = R.shape; adv = np.zeros((B, Tt), np.float32); last = np.zeros(B, np.float32)
    for t in reversed(range(Tt)):
        nextv = V[:, t+1] if t+1 < Tt else np.zeros(B, np.float32)
        delta = R[:, t] + gamma*nextv - V[:, t]; last = delta + gamma*lam*last; adv[:, t] = last
    return adv, adv + V


def rollout(net, env, greedy=False):
    B = env.B; obs = env.obs(); h = [None, None]
    O = np.zeros((B, 2, T, ODIM), np.float32); A = np.zeros((B, 2, T), int); SG = np.zeros((B, 2, T), int)
    LP = np.zeros((B, 2, T), np.float32); V = np.zeros((B, 2, T), np.float32); R = np.zeros((B, T), np.float32)
    for t in range(T):
        act = np.zeros((B, 2), int); sg = np.zeros((B, 2), int)
        for i in range(2):
            with torch.no_grad():
                al, sl, v, h[i] = net(torch.from_numpy(obs[:, i])[:, None, :], h[i])
            ad = torch.distributions.Categorical(logits=al[:, 0]); sd = torch.distributions.Categorical(logits=sl[:, 0])
            aa = al[:, 0].argmax(1) if greedy else ad.sample(); s = sl[:, 0].argmax(1) if greedy else sd.sample()
            O[:, i, t] = obs[:, i]; A[:, i, t] = aa.numpy(); SG[:, i, t] = s.numpy()
            LP[:, i, t] = (ad.log_prob(aa) + sd.log_prob(s)).numpy(); V[:, i, t] = v[:, 0].numpy()
            act[:, i] = aa.numpy(); sg[:, i] = s.numpy()
        R[:, t] = env.step(act, sg); obs = env.obs()
    return O, A, SG, LP, V, R


def train(iters, seed=0, B=256):
    torch.manual_seed(seed); net = Pol(); opt = torch.optim.Adam(net.parameters(), lr=3e-3)
    for it in range(iters):
        env = VecSync(B, 1000 + it)
        O, A, SG, LP, V, R = rollout(net, env)
        Rag = np.repeat(R[:, None, :], 2, axis=1)
        O2 = O.reshape(B*2, T, ODIM); A2, SG2 = A.reshape(B*2, T), SG.reshape(B*2, T)
        LP2, V2, R2 = LP.reshape(B*2, T), V.reshape(B*2, T), Rag.reshape(B*2, T)
        adv, ret = gae(R2, V2); adv = (adv - adv.mean()) / (adv.std() + 1e-8)
        Ot, At, SGt, LPt = torch.from_numpy(O2), torch.from_numpy(A2), torch.from_numpy(SG2), torch.from_numpy(LP2)
        advt, rett = torch.from_numpy(adv), torch.from_numpy(ret)
        for _ in range(4):
            al, sl, v, _ = net(Ot)
            ad = torch.distributions.Categorical(logits=al); sd = torch.distributions.Categorical(logits=sl)
            ratio = torch.exp(ad.log_prob(At) + sd.log_prob(SGt) - LPt)
            s1 = ratio*advt; s2 = torch.clamp(ratio, 0.8, 1.2)*advt
            loss = -torch.min(s1, s2).mean() + 0.5*((v-rett)**2).mean() - 0.02*(ad.entropy()+sd.entropy()).mean()
            opt.zero_grad(); loss.backward(); nn.utils.clip_grad_norm_(net.parameters(), 0.5); opt.step()
        if it % 100 == 0 or it == iters - 1:
            print(f"  iter {it:>3}: syncs/episode {_sync(net, 7000+it):.2f}", flush=True)
    return net


def _sync(net, seed, B=400, **kw):
    env = VecSync(B, seed, **kw); rollout(net, env, greedy=False); return env.synced()


if __name__ == "__main__":
    a = sys.argv
    iters = int(a[a.index("--iters")+1]) if "--iters" in a else 1000
    here = __file__.rsplit("/", 1)[0]
    print(f"OPEN-ENDED HUNT rung 5 (TEMPORAL): 2 agents fire SIMULTANEOUSLY, refractory C={C}; clock+obs+channel", flush=True)
    net = train(iters)
    norm = np.mean([_sync(net, 9000+k) for k in range(3)])
    nock = np.mean([_sync(net, 9100+k, noclock=True) for k in range(3)])
    blind = np.mean([_sync(net, 9200+k, blind=True) for k in range(3)])
    mute = np.mean([_sync(net, 9300+k, mute=True) for k in range(3)])
    print(f"\nFINGERPRINT (syncs/episode): normal {norm:.2f} | no-clock {nock:.2f} | blind {blind:.2f} | mute {mute:.2f}", flush=True)
    print(f"  needs-CLOCK = {norm-nock:+.2f}; needs-OBSERVATION = {norm-blind:+.2f}; needs-CHANNEL = {norm-mute:+.2f}", flush=True)
    if norm - nock >= 0.8: mech = "explicit temporal focal point (clock)"
    elif nock >= 2.0: mech = "internal-count temporal focal point (memory; clock not needed)"
    elif norm - blind >= 0.8: mech = "observation (reactive follow)"
    elif norm - mute >= 0.8: mech = "communication"
    else: mech = "unclear"
    print(f"RUNG-5 MECHANISM: {mech} -- temporal coordination reduces to KNOWN primitives (no new medium)" , flush=True)
    torch.save(net.state_dict(), f"{here}/openhunt6_net.pt")
    json.dump({"normal": round(float(norm),2), "noclock": round(float(nock),2), "blind": round(float(blind),2),
               "mute": round(float(mute),2), "mechanism": mech}, open(f"{here}/openhunt6_result.json","w"))
