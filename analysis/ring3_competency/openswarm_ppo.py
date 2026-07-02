"""OPEN-ENDED hunt -- NEW TERRITORY (user: climb into a materially different regime). SWARM-scale COLLECTIVE
DECISION under distributed noisy evidence (25 agents) -- collective computation, not pairwise coordination.

A hidden binary state s (act-is-GOOD=1 / act-is-BAD=0), random per episode. Each of N agents gets a PRIVATE
NOISY signal about s (accuracy p). Each of T rounds every agent picks ACT or ABSTAIN. If the fraction acting
>= quorum Q: all agents get +1 if s=1, -1 if s=0; if < Q: 0. So the swarm should collectively reach quorum
IFF the (distributed, noisy) evidence says s=1 -- a collective decision aggregating evidence no single agent has.

Questions the swarm regime raises:
  1. Does the swarm achieve COLLECTIVE accuracy > individual signal accuracy p (wisdom of crowds / Condorcet)?
  2. Is that a NEW coordination medium, or does it reduce to INDEPENDENT voting + the environment's vote-count
     (no medium: each agent votes its signal, the quorum rule aggregates)? Fingerprint: BLIND (remove the
     observed crowd fraction) -> if collective accuracy SURVIVES, aggregation is done by independent voting, not
     by any coordination medium (completeness holds). If it COLLAPSES, the swarm needs to OBSERVE the crowd
     (quorum/recruitment). Also HERD-check: does observing the crowd HELP or HURT vs independent voting?

Run on VPS epc-venv: python openswarm_ppo.py [--iters 800]
"""
import numpy as np, sys, json, torch, torch.nn as nn

N_AG = 25; T = 12; Q = 0.5; P = 0.70; H = 96
ODIM = 1 + 1 + 1                        # own private signal, observed crowd-fraction-acting last round, round-frac
torch.set_num_threads(8)


class Pol(nn.Module):
    def __init__(self):
        super().__init__(); self.gru = nn.GRU(ODIM, H, batch_first=True)
        self.act = nn.Linear(H, 2); self.v = nn.Linear(H, 1)

    def forward(self, x, h=None):
        z, h = self.gru(x, h); return self.act(z), self.v(z).squeeze(-1), h


class VecSwarm:
    def __init__(self, B, seed, blind=False):
        self.B, self.blind = B, blind; self.rng = np.random.default_rng(seed); self.reset()

    def reset(self):
        B = self.B
        self.s = self.rng.integers(0, 2, size=B)                          # hidden state per episode
        flip = self.rng.random((B, N_AG)) > P
        self.sig = np.where(flip, 1 - self.s[:, None], self.s[:, None])   # private noisy signals (accuracy P)
        self.frac = np.zeros(B, np.float32); self.t = 0; self.last_frac = np.zeros(B, np.float32)
        return self.obs()

    def obs(self):
        B = self.B; out = np.zeros((B, N_AG, ODIM), np.float32)
        out[:, :, 0] = self.sig
        out[:, :, 1] = 0.0 if self.blind else self.frac[:, None]
        out[:, :, 2] = self.t / T
        return out

    def step(self, act):
        B = self.B; self.frac = act.mean(1).astype(np.float32); self.last_frac = self.frac.copy()
        reached = self.frac >= Q
        r = np.where(reached, np.where(self.s == 1, 1.0, -1.0), 0.0).astype(np.float32)
        self.t += 1
        return np.repeat(r[:, None], N_AG, axis=1)                        # shared collective reward

    def decision_correct(self):
        reached = self.last_frac >= Q
        return ((reached & (self.s == 1)) | (~reached & (self.s == 0))).astype(np.float32)


def gae(R, V, gamma=0.9, lam=0.95):
    B, Tt = R.shape; adv = np.zeros((B, Tt), np.float32); last = np.zeros(B, np.float32)
    for t in reversed(range(Tt)):
        nextv = V[:, t+1] if t+1 < Tt else np.zeros(B, np.float32)
        delta = R[:, t] + gamma*nextv - V[:, t]; last = delta + gamma*lam*last; adv[:, t] = last
    return adv, adv + V


def rollout(net, env, greedy=False):
    B, A = env.B, N_AG; obs = env.obs(); h = None
    O = np.zeros((B, A, T, ODIM), np.float32); AC = np.zeros((B, A, T), int)
    LP = np.zeros((B, A, T), np.float32); V = np.zeros((B, A, T), np.float32); R = np.zeros((B, T), np.float32)
    for t in range(T):
        flat = torch.from_numpy(obs.reshape(B*A, 1, ODIM))
        with torch.no_grad():
            al, v, h = net(flat, h)
        al = al[:, 0].reshape(B, A, 2); v = v[:, 0].reshape(B, A)
        d = torch.distributions.Categorical(logits=al); a = al.argmax(-1) if greedy else d.sample()
        O[:, :, t] = obs; AC[:, :, t] = a.numpy(); LP[:, :, t] = d.log_prob(a).numpy(); V[:, :, t] = v.numpy()
        R[:, t] = env.step(a.numpy())[:, 0]; obs = env.obs()
    return O, AC, LP, V, R


def train(iters, seed=0, B=96, blind=False):
    torch.manual_seed(seed); net = Pol(); opt = torch.optim.Adam(net.parameters(), lr=3e-3)
    for it in range(iters):
        env = VecSwarm(B, 1000 + it, blind=blind)
        O, AC, LP, V, R = rollout(net, env)
        Rag = np.repeat(R[:, None, :], N_AG, axis=1)
        O2 = O.reshape(B*N_AG, T, ODIM); AC2 = AC.reshape(B*N_AG, T)
        LP2, V2, R2 = LP.reshape(B*N_AG, T), V.reshape(B*N_AG, T), Rag.reshape(B*N_AG, T)
        adv, ret = gae(R2, V2); adv = (adv - adv.mean()) / (adv.std() + 1e-8)
        Ot, ACt, LPt = torch.from_numpy(O2), torch.from_numpy(AC2), torch.from_numpy(LP2)
        advt, rett = torch.from_numpy(adv), torch.from_numpy(ret)
        for _ in range(4):
            al, v, _ = net(Ot)
            d = torch.distributions.Categorical(logits=al)
            ratio = torch.exp(d.log_prob(ACt) - LPt); s1 = ratio*advt; s2 = torch.clamp(ratio, 0.8, 1.2)*advt
            loss = -torch.min(s1, s2).mean() + 0.5*((v-rett)**2).mean() - 0.01*d.entropy().mean()
            opt.zero_grad(); loss.backward(); nn.utils.clip_grad_norm_(net.parameters(), 0.5); opt.step()
        if it % 100 == 0 or it == iters - 1:
            print(f"  iter {it:>3}: collective-accuracy {_acc(net, 7000+it):.2f}", flush=True)
    return net


def _acc(net, seed, B=300, **kw):
    env = VecSwarm(B, seed, **kw); rollout(net, env, greedy=False); return float(env.decision_correct().mean())


if __name__ == "__main__":
    a = sys.argv
    iters = int(a[a.index("--iters")+1]) if "--iters" in a else 800
    here = __file__.rsplit("/", 1)[0]
    print(f"OPEN-ENDED HUNT swarm: {N_AG}-agent collective decision under noisy evidence (individual acc P={P})", flush=True)
    print("== training WITH crowd observation ==", flush=True); net = train(iters)
    print("== training BLIND (fair independent-voting / Condorcet baseline) ==", flush=True); netb = train(iters, blind=True)
    norm = np.mean([_acc(net, 9000+k) for k in range(4)])
    abl = np.mean([_acc(net, 9200+k, blind=True) for k in range(4)])           # trained-then-ablated (degenerate)
    fair = np.mean([_acc(netb, 9200+k, blind=True) for k in range(4)])          # FAIR blind-trained baseline
    print(f"\nFINGERPRINT (collective decision accuracy): with-obs {norm:.2f} | blind-TRAINED(Condorcet) {fair:.2f} | trained-then-ablated {abl:.2f} | individual P={P:.2f}", flush=True)
    wisdom = norm - P; obs_gain = norm - fair; condorcet_gain = fair - P
    print(f"  wisdom-of-crowds (with-obs - individual) = {wisdom:+.2f}; Condorcet gain (blind-trained - individual) = {condorcet_gain:+.2f}; crowd-obs adds = {obs_gain:+.2f}", flush=True)
    if condorcet_gain >= 0.1 and obs_gain >= 0.05:
        mech = "QUORUM/RECRUITMENT beats INDEPENDENT VOTING: both aggregate distributed evidence (Condorcet), but crowd-observation (quorum) adds a further gain -- an AGGREGATE-observation MODE that only appears at swarm scale"
    elif condorcet_gain >= 0.1:
        mech = "INDEPENDENT VOTING (Condorcet) suffices -- wisdom-of-crowds without a coordination medium; crowd-obs adds little"
    else:
        mech = "crowd-observation (quorum) is REQUIRED for wisdom-of-crowds -- independent voting insufficient here"
    print(f"SWARM MECHANISM: {mech}", flush=True)
    torch.save(net.state_dict(), f"{here}/openswarm_net.pt")
    json.dump({"with_obs": round(float(norm),2), "blind_trained_condorcet": round(float(fair),2),
               "trained_then_ablated": round(float(abl),2), "individual": P, "wisdom": round(float(wisdom),2),
               "condorcet_gain": round(float(condorcet_gain),2), "crowd_obs_adds": round(float(obs_gain),2),
               "mechanism": mech}, open(f"{here}/openswarm_result.json","w"))
