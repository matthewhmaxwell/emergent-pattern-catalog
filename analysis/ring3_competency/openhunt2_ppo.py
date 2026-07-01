"""OPEN-ENDED N-agent hunt -- RUNG 1 (climb the cascade: REMOVE the focal point rung 0 exploited).

Rung 0 finding: with food + obs + signal, agents coordinated the 2-agent rendezvous through the shared FOOD
as a Schelling FOCAL POINT -- neither observation nor the channel was used (both ablation-invariant). To climb
the cascade we delete every free focal point and see what mechanism is forced next.

RUNG 1 world: pure 2-agent RENDEZVOUS on a TOROIDAL grid (wraps -> no geometric centre) with EGOCENTRIC-only
observation (each agent sees ONLY the other's RELATIVE position + the other's last signal -- NO absolute
coordinates, NO food, NO landmark). Every free focal point (environmental / geometric / coordinate / index)
is gone. Reward +1 when the two agents share a cell (then both respawn). Fingerprint:
  - BLIND (zero the other's relative position) -> if it collapses, rendezvous RIDES OBSERVATION (home in on
    the partner). The only thing left is the signal channel.
  - MUTE (scramble the other's signal) -> if invariant, the channel is NOT used.
Expected (known): observation-based pursuit rendezvous (Alpern). The interesting question the cascade asks:
with all focal points removed, is COMMUNICATION finally forced, or does observation alone still suffice?

Run on VPS epc-venv: python openhunt2_ppo.py [--iters 800]
"""
import numpy as np, sys, json, torch, torch.nn as nn

N = 9; N_AG = 2; K = 3; T = 30; H = 96; NMOVE = 5
ODIM = 2 + K                                          # other's RELATIVE (toroidal) position, other's last signal
DIRS = np.array([(0, 0), (1, 0), (-1, 0), (0, 1), (0, -1)])
torch.set_num_threads(8)


class Pol(nn.Module):
    def __init__(self):
        super().__init__(); self.gru = nn.GRU(ODIM, H, batch_first=True)
        self.mv = nn.Linear(H, NMOVE); self.sg = nn.Linear(H, K); self.v = nn.Linear(H, 1)

    def forward(self, x, h=None):
        z, h = self.gru(x, h); return self.mv(z), self.sg(z), self.v(z).squeeze(-1), h


def torus_rel(a, b):                                  # shortest signed vector from a to b on the torus, in [-N/2, N/2]
    return ((b - a + N // 2) % N) - N // 2


class VecRdv:
    def __init__(self, B, seed, blind=False, mute=False):
        self.B, self.blind, self.mute = B, blind, mute
        self.rng = np.random.default_rng(seed); self.reset()

    def reset(self):
        B = self.B
        self.ap = self.rng.integers(0, N, size=(B, 2, 2))
        self.lsig = np.zeros((B, 2), int); self.met = np.zeros(B, int)
        return self.obs()

    def obs(self):
        B = self.B; out = np.zeros((B, 2, ODIM), np.float32)
        for i in range(2):
            j = 1 - i; o = 0
            rel = np.zeros((B, 2), np.float32) if self.blind else torus_rel(self.ap[:, i], self.ap[:, j]) / (N // 2)
            out[:, i, o:o+2] = rel; o += 2
            sj = self.rng.integers(0, K, size=B) if self.mute else self.lsig[:, j]
            out[np.arange(B), i, o + sj] = 1.0
        return out

    def step(self, mv, sg):
        B = self.B; rew = np.full(B, -0.01, np.float32)
        for i in range(2): self.ap[:, i] = (self.ap[:, i] + DIRS[mv[:, i]]) % N     # toroidal wrap
        self.lsig = sg.copy()
        same = np.all(self.ap[:, 0] == self.ap[:, 1], axis=1)
        for b in range(B):
            if same[b]:
                rew[b] += 1.0; self.met[b] += 1
                self.ap[b] = self.rng.integers(0, N, size=(2, 2))                    # respawn both
        return rew

    def meets(self): return float(self.met.mean())


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
        env = VecRdv(B, 1000 + it)
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
            loss = -torch.min(s1, s2).mean() + 0.5*((v-rett)**2).mean() - 0.01*(md.entropy()+sd.entropy()).mean()
            opt.zero_grad(); loss.backward(); nn.utils.clip_grad_norm_(net.parameters(), 0.5); opt.step()
        if it % 100 == 0 or it == iters - 1:
            print(f"  iter {it:>3}: meets/episode {_meets(net, 7000+it):.2f}", flush=True)
    return net


def _meets(net, seed, B=400, **kw):
    env = VecRdv(B, seed, **kw); rollout(net, env, greedy=False); return env.meets()


if __name__ == "__main__":
    a = sys.argv
    iters = int(a[a.index("--iters")+1]) if "--iters" in a else 800
    here = __file__.rsplit("/", 1)[0]
    print(f"OPEN-ENDED HUNT rung 1: 2-agent TOROIDAL rendezvous, EGOCENTRIC obs (no focal point)", flush=True)
    net = train(iters)
    norm = np.mean([_meets(net, 9000+k) for k in range(3)])
    blind = np.mean([_meets(net, 9200+k, blind=True) for k in range(3)])
    mute = np.mean([_meets(net, 9300+k, mute=True) for k in range(3)])
    print(f"\nFINGERPRINT (meets/episode): normal {norm:.2f} | blind {blind:.2f} | mute {mute:.2f}", flush=True)
    uses_obs = norm - blind; uses_sig = norm - mute
    print(f"  rides-OBSERVATION = {uses_obs:+.2f}; rides-CHANNEL = {uses_sig:+.2f}", flush=True)
    mech = ("observation-rendezvous" if uses_obs >= 0.5 and uses_obs > uses_sig
            else "channel-rendezvous" if uses_sig >= 0.5 else "focal-point-survived(unclear)")
    print(f"RUNG-1 MECHANISM: {mech}", flush=True)
    torch.save(net.state_dict(), f"{here}/openhunt2_net.pt")
    json.dump({"normal": round(float(norm),2), "blind": round(float(blind),2), "mute": round(float(mute),2),
               "mechanism": mech}, open(f"{here}/openhunt2_result.json","w"))
