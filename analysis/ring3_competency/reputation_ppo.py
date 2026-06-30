"""Track-1 probe (audit gap #2): REPUTATION / indirect reciprocity -- a 4th coordination MEDIUM (social
memory), and the mechanism that could stabilize the honest/deceptive behaviour minimal signaling could
not. A focal agent (recurrent) plays T rounds, each against ONE of M partners (chosen at random, so
partners RECUR). Each partner has a fixed hidden TYPE (cooperator or defector). Each round the focal sees
WHICH partner it faces (an id) and must HELP (cooperate, costly) or PASS. Help a cooperator -> +1; help a
defector -> -1 (exploited); pass -> 0. After the round the partner's type is REVEALED (third-party
observation). To do well the focal must build a REPUTATION map (partner-id -> type) in memory from the
revealed stream and condition: help known cooperators, pass known defectors.

Diagnostics: (a) ANONYMIZE -- hide the partner id -> the focal cannot tell who it faces -> reputations are
useless -> payoff collapses to ~0; (b) MEMORY-WIPE -- reset the GRU each round -> cannot accumulate
reputations -> collapse. Genuine reputation = payoff(normal) >> payoff(anonymized) AND >> payoff(memoryless).

Run on VPS epc-venv: python reputation_ppo.py [--iters 600]
"""
import numpy as np, sys, json, os, torch, torch.nn as nn

M = 5; T = 30; H = 64; ODIM = M + M + 1; NACT = 2   # cur-partner id(M), last-partner id(M), last-type(+1/-1/0)
torch.set_num_threads(8)


class Pol(nn.Module):
    def __init__(self):
        super().__init__(); self.gru = nn.GRU(ODIM, H, batch_first=True)
        self.pi = nn.Linear(H, NACT); self.v = nn.Linear(H, 1)

    def forward(self, x, h=None):
        out, h = self.gru(x, h); return self.pi(out), self.v(out).squeeze(-1), h


def rollout(net, B, seed, greedy=False, anon=False, memwipe=False):
    rng = np.random.default_rng(seed)
    ptype = (rng.random((B, M)) < 0.5)                       # each partner: cooperator(True)/defector
    O = np.zeros((B, T, ODIM), np.float32); A = np.zeros((B, T), int); LP = np.zeros((B, T), np.float32)
    V = np.zeros((B, T), np.float32); R = np.zeros((B, T), np.float32)
    last_pid = -np.ones(B, int); last_type = np.zeros(B)
    coop_help = 0; coop_n = 0; def_pass = 0; def_n = 0
    h = None
    for t in range(T):
        pid = rng.integers(0, M, size=B)
        ob = np.zeros((B, ODIM), np.float32)
        if not anon:
            ob[np.arange(B), pid] = 1.0
            for b in range(B):
                if last_pid[b] >= 0: ob[b, M + last_pid[b]] = 1.0
            ob[:, 2 * M] = last_type
        O[:, t] = ob
        hin = None if memwipe else h
        with torch.no_grad():
            lg, v, h = net(torch.from_numpy(ob)[:, None, :], hin)
        d = torch.distributions.Categorical(logits=lg[:, 0]); a = lg[:, 0].argmax(1) if greedy else d.sample()
        A[:, t] = a.numpy(); LP[:, t] = d.log_prob(a).numpy(); V[:, t] = v[:, 0].numpy()
        an = a.numpy(); coop = ptype[np.arange(B), pid]
        r = np.where(an == 1, np.where(coop, 1.0, -1.0), 0.0).astype(np.float32)
        R[:, t] = r
        coop_help += ((an == 1) & coop).sum(); coop_n += coop.sum()
        def_pass += ((an == 0) & ~coop).sum(); def_n += (~coop).sum()
        last_pid = pid; last_type = np.where(coop, 1.0, -1.0)
    avg_r = float(R.sum(1).mean() / T)
    discr = float(0.5 * (coop_help / max(coop_n, 1)) + 0.5 * (def_pass / max(def_n, 1)))   # help coops + pass defectors
    return O, A, LP, V, R, dict(avg_r=avg_r, discr=discr)


def gae(R, V, gamma=0.95, lam=0.95):
    B, Tt = R.shape; adv = np.zeros((B, Tt), np.float32); last = np.zeros(B, np.float32)
    for t in reversed(range(Tt)):
        nextv = V[:, t + 1] if t + 1 < Tt else np.zeros(B, np.float32)
        delta = R[:, t] + gamma * nextv - V[:, t]
        last = delta + gamma * lam * last; adv[:, t] = last
    return adv, adv + V


def train(iters, seed=0, B=384):
    torch.manual_seed(seed); net = Pol(); opt = torch.optim.Adam(net.parameters(), lr=3e-3)
    for it in range(iters):
        O, A, LP, V, R, m = rollout(net, B, 1000 + it)
        adv, ret = gae(R, V); adv = (adv - adv.mean()) / (adv.std() + 1e-8)
        Ot, At, LPt = torch.from_numpy(O), torch.from_numpy(A), torch.from_numpy(LP)
        advt, rett = torch.from_numpy(adv), torch.from_numpy(ret)
        for _ in range(4):
            lg, v, _ = net(Ot); d = torch.distributions.Categorical(logits=lg)
            ratio = torch.exp(d.log_prob(At) - LPt)
            s1 = ratio * advt; s2 = torch.clamp(ratio, 0.8, 1.2) * advt
            loss = -torch.min(s1, s2).mean() + 0.5 * ((v - rett) ** 2).mean() - 0.01 * d.entropy().mean()
            opt.zero_grad(); loss.backward(); nn.utils.clip_grad_norm_(net.parameters(), 0.5); opt.step()
        if it % 100 == 0 or it == iters - 1:
            st = rollout(net, 800, 7000 + it, greedy=True)[5]
            print(f"  iter {it:>3}: avg-reward {st['avg_r']:+.2f}  discrimination {st['discr']:.2f}", flush=True)
    return net


def ev(net, seed, **kw): return rollout(net, 1000, seed, greedy=True, **kw)[5]


if __name__ == "__main__":
    a = sys.argv
    iters = int(a[a.index("--iters") + 1]) if "--iters" in a else 600
    here = __file__.rsplit("/", 1)[0]
    print(f"REPUTATION probe: indirect reciprocity over {M} recurring partners, {T} rounds", flush=True)
    net = train(iters)
    nm = ev(net, 9000); an = ev(net, 9100, anon=True); mw = ev(net, 9200, memwipe=True)
    print(f"\nFINAL avg-reward/discrimination: normal {nm['avg_r']:+.2f}/{nm['discr']:.2f} | anonymized {an['avg_r']:+.2f}/{an['discr']:.2f} | memory-wiped {mw['avg_r']:+.2f}/{mw['discr']:.2f}", flush=True)
    rep = nm['discr'] >= 0.7 and (nm['discr'] - an['discr']) >= 0.15 and (nm['discr'] - mw['discr']) >= 0.15
    print(f"VERDICT: {'GENUINE REPUTATION (social memory)' if rep else 'partial/none'} "
          f"(discrim {nm['discr']:.2f}; needs id={nm['discr']-an['discr']:+.2f}, needs memory={nm['discr']-mw['discr']:+.2f})", flush=True)
    torch.save(net.state_dict(), f"{here}/reputation_net.pt")
    json.dump({"discrim": round(nm['discr'], 2), "discrim_anon": round(an['discr'], 2), "discrim_memwipe": round(mw['discr'], 2),
               "avg_r": round(nm['avg_r'], 2), "reputation": bool(rep)}, open(f"{here}/reputation_result.json", "w"))
