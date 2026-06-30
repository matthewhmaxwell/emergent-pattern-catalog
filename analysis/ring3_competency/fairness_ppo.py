"""Phase D4: forced negotiation -- remove the geometric asymmetry that D3 used, so the ONLY way to break
symmetry is the channel and/or memory. Abstract ITERATED contest: two FULLY SYMMETRIC agents (shared
recurrent policy, identical observations) play ROUNDS rounds. Each round each agent picks CLAIM or YIELD
(+ emits a symbol). Per-agent payoff: one claims -> claimer +2, yielder +0.5; both claim -> -1 each
(collision); both yield -> 0. Each agent is SELFISH (own reward), so both want to claim, but both-claim
is mutually bad: an iterated Hawk-Dove. Symmetric agents on identical obs can only reliably get
exactly-one-claim by establishing a TURN ORDER -- which needs the channel and/or memory of past rounds.
Over many rounds the natural symmetric-cooperative solution is ALTERNATION (take turns) = an emergent
FAIRNESS norm. We hunt for it and ask whether the literature has a clean name for the specific mechanism.

Metrics: efficiency = fraction of rounds with exactly-one-claim; fairness = 1 - 2*|agent0's share of the
claims - 0.5| (1 = perfect alternation/fair, 0 = one agent dominates); alt = lag-1 anti-correlation of an
agent's own claim (1 = strict alternation). Discriminators: channel-scramble (random symbols -> collapse
=> turn-order rides the CHANNEL) and memory-ablate (reset hidden each round -> collapse => rides MEMORY).

Run on VPS epc-venv: python fairness_ppo.py [--iters 700]
"""
import numpy as np, sys, json, os, torch, torch.nn as nn

ROUNDS = 14; K = 4; H = 64
ODIM = 2 + 2 + K + 1             # own last action(2), other last action(2), other last symbol(K), first-flag(1)
torch.set_num_threads(8)


class Pol(nn.Module):
    def __init__(self):
        super().__init__(); self.gru = nn.GRU(ODIM, H, batch_first=True)
        self.act = nn.Linear(H, 2); self.sym = nn.Linear(H, K); self.v = nn.Linear(H, 1)

    def forward(self, x, h=None):
        out, h = self.gru(x, h); return self.act(out), self.sym(out), self.v(out).squeeze(-1), h


def rollout(net, B, seed, greedy=False, scramble=False, memwipe=False):
    rng = np.random.default_rng(seed)
    la = np.zeros((B, 2), int); lsym = np.zeros((B, 2), int)     # last action / symbol per agent
    O = np.zeros((B, 2, ROUNDS, ODIM), np.float32); A = np.zeros((B, 2, ROUNDS), int); SY = np.zeros((B, 2, ROUNDS), int)
    LP = np.zeros((B, 2, ROUNDS), np.float32); V = np.zeros((B, 2, ROUNDS), np.float32); R = np.zeros((B, 2, ROUNDS), np.float32)
    claims = np.zeros((B, 2), int)
    h = [None, None]
    for t in range(ROUNDS):
        acts = np.zeros((B, 2), int)
        for i in range(2):
            oth = 1 - i; ob = np.zeros((B, ODIM), np.float32); o = 0
            ob[np.arange(B), o + la[:, i]] = 1.0; o += 2
            ob[np.arange(B), o + la[:, oth]] = 1.0; o += 2
            sy = rng.integers(0, K, size=B) if scramble else lsym[:, oth]
            ob[np.arange(B), o + sy] = 1.0; o += K
            ob[:, o] = float(t == 0)
            O[:, i, t] = ob
            hin = None if memwipe else h[i]
            with torch.no_grad():
                al, sl, v, h[i] = net(torch.from_numpy(ob)[:, None, :], hin)
            ad, sd = torch.distributions.Categorical(logits=al[:, 0]), torch.distributions.Categorical(logits=sl[:, 0])
            a = al[:, 0].argmax(1) if greedy else ad.sample(); s = sl[:, 0].argmax(1) if greedy else sd.sample()
            A[:, i, t] = a.numpy(); SY[:, i, t] = s.numpy()
            LP[:, i, t] = (ad.log_prob(a) + sd.log_prob(s)).numpy(); V[:, i, t] = v[:, 0].numpy()
            acts[:, i] = a.numpy()
        a0, a1 = acts[:, 0], acts[:, 1]
        one0 = (a0 == 1) & (a1 == 0); one1 = (a1 == 1) & (a0 == 0)
        both = (a0 == 1) & (a1 == 1)
        R[:, 0, t] = 2.0 * one0 + 0.5 * one1 - 1.0 * both
        R[:, 1, t] = 2.0 * one1 + 0.5 * one0 - 1.0 * both
        claims[:, 0] += a0; claims[:, 1] += a1
        la = acts.copy(); lsym = np.stack([SY[:, 0, t], SY[:, 1, t]], 1)
    eff = float((((A[:, 0] + A[:, 1]) == 1).mean()))            # exactly one claim
    tot = claims.sum(1) + 1e-9; share0 = claims[:, 0] / tot
    fairness = float((1 - 2 * np.abs(share0 - 0.5)).mean())
    alt = float(np.mean([(A[:, i, 1:] != A[:, i, :-1]).mean() for i in range(2)]))   # how often an agent flips
    return O, A, SY, LP, V, R, dict(eff=eff, fairness=fairness, alt=alt)


def gae(R, V, gamma=0.95, lam=0.95):
    Aa, Tt = R.shape; adv = np.zeros((Aa, Tt), np.float32); last = np.zeros(Aa, np.float32)
    for t in reversed(range(Tt)):
        nextv = V[:, t + 1] if t + 1 < Tt else np.zeros(Aa, np.float32)
        delta = R[:, t] + gamma * nextv - V[:, t]
        last = delta + gamma * lam * last; adv[:, t] = last
    return adv, adv + V


def train(iters, seed=0, B=512):
    torch.manual_seed(seed); net = Pol(); opt = torch.optim.Adam(net.parameters(), lr=3e-3)
    for it in range(iters):
        O, A, SY, LP, V, R, m = rollout(net, B, 1000 + seed * 99999 + it)
        O2 = O.reshape(2 * B, ROUNDS, ODIM); A2, SY2, LP2 = A.reshape(2 * B, ROUNDS), SY.reshape(2 * B, ROUNDS), LP.reshape(2 * B, ROUNDS)
        V2, R2 = V.reshape(2 * B, ROUNDS), R.reshape(2 * B, ROUNDS)
        adv, ret = gae(R2, V2); adv = (adv - adv.mean()) / (adv.std() + 1e-8)
        Ot, At, SYt, LPt = torch.from_numpy(O2), torch.from_numpy(A2), torch.from_numpy(SY2), torch.from_numpy(LP2)
        advt, rett = torch.from_numpy(adv), torch.from_numpy(ret)
        for _ in range(4):
            al, sl, v, _ = net(Ot)
            ad, sd = torch.distributions.Categorical(logits=al), torch.distributions.Categorical(logits=sl)
            ratio = torch.exp(ad.log_prob(At) + sd.log_prob(SYt) - LPt)
            s1 = ratio * advt; s2 = torch.clamp(ratio, 0.8, 1.2) * advt
            loss = -torch.min(s1, s2).mean() + 0.5 * ((v - rett) ** 2).mean() - 0.02 * (ad.entropy() + sd.entropy()).mean()
            opt.zero_grad(); loss.backward(); nn.utils.clip_grad_norm_(net.parameters(), 0.5); opt.step()
        if it % 100 == 0 or it == iters - 1:
            _, _, _, _, _, _, ev = rollout(net, 800, 7000 + it, greedy=False)
            print(f"  iter {it:>3}: efficiency {ev['eff']:.2f}  fairness {ev['fairness']:.2f}  flip-rate {ev['alt']:.2f}", flush=True)
    return net


if __name__ == "__main__":
    a = sys.argv
    iters = int(a[a.index("--iters") + 1]) if "--iters" in a else 700
    here = __file__.rsplit("/", 1)[0]
    print(f"PHASE D4 forced negotiation / fairness: iterated symmetric contest, {ROUNDS} rounds", flush=True)
    net = train(iters)
    def ev(**kw):
        outs = [rollout(net, 800, 9000 + k, greedy=False, **kw)[6] for k in range(3)]
        return {k: round(float(np.mean([o[k] for o in outs])), 2) for k in ("eff", "fairness", "alt")}
    norm = ev(); scr = ev(scramble=True); mem = ev(memwipe=True)
    print(f"\nFINAL: normal {norm} | channel-scramble {scr} | memory-wiped {mem}", flush=True)
    uses_comms = norm["eff"] - scr["eff"] >= 0.15; uses_mem = norm["eff"] - mem["eff"] >= 0.15
    fairnorm = norm["eff"] >= 0.7 and norm["fairness"] >= 0.7
    print(f"VERDICT: {'EMERGENT FAIR ALTERNATION' if fairnorm else 'partial/none'} "
          f"(eff {norm['eff']}, fairness {norm['fairness']}); mechanism comms={uses_comms} memory={uses_mem}", flush=True)
    torch.save(net.state_dict(), f"{here}/fairness_net.pt")
    json.dump({"normal": norm, "scramble": scr, "memwipe": mem, "uses_comms": bool(uses_comms),
               "uses_memory": bool(uses_mem), "fair_norm": bool(fairnorm)}, open(f"{here}/fairness_result.json", "w"))
