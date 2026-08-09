"""Full-scale Transformer confirmation: does attention hold BOTH compositional channels at full budget?
Reduced-scale (1500it, 3 seeds) gave mean nobody 0.45, freezeblock 0.31 - BOTH loaded, unlike every GRU variant.
This run: 8000it, 3 seeds, causal Transformer. If both necessities hold >=0.3 -> first compositional competency."""
import numpy as np, torch, torch.nn as nn, json, time, sys
torch.set_num_threads(4)                                   # shared box - avoid thread thrash
sys.path.insert(0, ".")
import adapthetero as E
od = E.AdaptHetero(1, 0).odim
SEEDS = [0, 1, 2]; ITERS = 8000; B = 256

class TransformerPolicy(nn.Module):
    def __init__(s, odim=od, d=64, nhead=4, layers=2, maxT=80):
        super().__init__(); s.recurrent = True
        s.inp = nn.Linear(odim, d); s.pos = nn.Parameter(torch.zeros(1, maxT, d))
        layer = nn.TransformerEncoderLayer(d_model=d, nhead=nhead, dim_feedforward=2 * d, batch_first=True, dropout=0.0)
        s.enc = nn.TransformerEncoder(layer, num_layers=layers)
        s.pi = nn.Linear(d, E.NACT); s.v = nn.Linear(d, 1)
    def forward(s, x, h=None):
        B, Tq, _ = x.shape; z = s.inp(x) + s.pos[:, :Tq]
        mask = torch.triu(torch.full((Tq, Tq), float("-inf")), diagonal=1)
        z = s.enc(z, mask=mask); return s.pi(z), s.v(z).squeeze(-1), None

def rollout(net, env):
    B = env.B; T = E.T
    OH = np.zeros((B, T, od), np.float32); A = np.zeros((B, T), np.int64)
    LP = np.zeros((B, T), np.float32); V = np.zeros((B, T), np.float32); R = np.zeros((B, T), np.float32)
    o = env.obs()[:, 0, :]; prev = np.zeros(B, np.float32)
    for t in range(T):
        OH[:, t] = o
        with torch.no_grad(): lg, v, _ = net(torch.from_numpy(OH[:, :t + 1]))
        d = torch.distributions.Categorical(logits=lg[:, -1]); a = d.sample()
        A[:, t] = a.numpy(); LP[:, t] = d.log_prob(a).numpy(); V[:, t] = v[:, -1].numpy()
        env.step(a.numpy()[:, None]); o = env.obs()[:, 0, :]
        cur = env._surv.copy(); R[:, t] = cur - prev; prev = cur
    return OH, A, LP, V, R

def score(net, seed, Bn=300, **abl):
    env = E.AdaptHetero(Bn, seed, n_ag=1, **abl); T = E.T
    OH = np.zeros((Bn, T, od), np.float32); o = env.obs()[:, 0, :]
    for t in range(T):
        OH[:, t] = o
        with torch.no_grad(): lg, _, _ = net(torch.from_numpy(OH[:, :t + 1]))
        env.step(lg[:, -1].argmax(-1).numpy()[:, None]); o = env.obs()[:, 0, :]
    return env.score()

def nec(net, ch, n=8):
    base = np.mean([score(net, 3000 + k) for k in range(n)])
    ab = np.mean([score(net, 3000 + k, **{ch: True}) for k in range(n)])
    return float(1 - ab / max(base, 1e-6))

res = {}
for sd in SEEDS:
    torch.manual_seed(sd); np.random.seed(sd)
    net = TransformerPolicy(); opt = torch.optim.Adam(net.parameters(), lr=3e-4); t0 = time.time()
    for it in range(ITERS):
        env = E.AdaptHetero(B, sd * 99991 + it)
        OH, A, LP, V, R = rollout(net, env); T = E.T
        adv = np.zeros_like(R); last = np.zeros(B, np.float32)
        for t in reversed(range(T)):
            nv = V[:, t + 1] if t + 1 < T else np.zeros(B, np.float32)
            delta = R[:, t] + 0.99 * nv - V[:, t]; last = delta + 0.95 * 0.99 * last; adv[:, t] = last
        ret = adv + V; adv = (adv - adv.mean()) / (adv.std() + 1e-8)
        Ot = torch.from_numpy(OH); At = torch.from_numpy(A); LPt = torch.from_numpy(LP)
        advt = torch.from_numpy(adv.astype(np.float32)); rett = torch.from_numpy(ret.astype(np.float32))
        for _ in range(3):
            lg, v, _ = net(Ot); d = torch.distributions.Categorical(logits=lg)
            nlp = d.log_prob(At); ratio = torch.exp(nlp - LPt)
            s1 = ratio * advt; s2 = torch.clamp(ratio, 0.8, 1.2) * advt
            loss = -torch.min(s1, s2).mean() + 0.5 * ((v - rett) ** 2).mean() - 0.01 * d.entropy().mean()
            opt.zero_grad(); loss.backward(); nn.utils.clip_grad_norm_(net.parameters(), 0.5); opt.step()
        if it % 500 == 0 or it == ITERS - 1:
            sc = np.mean([score(net, 999 + k) for k in range(3)])
            print(f"  s{sd}[{it}] surv {sc:.3f} nobody {nec(net,'nobody',3):.3f} freezeblock {nec(net,'freezeblock',3):.3f}", flush=True)
    torch.save(net.state_dict(), f"adaptheteroTF_s{sd}.pt")
    res[f"seed{sd}"] = {"survival": round(float(np.mean([score(net, 999 + k) for k in range(3)])), 4),
                        "nobody_nec": round(nec(net, "nobody"), 3), "freezeblock_nec": round(nec(net, "freezeblock"), 3),
                        "sec": round(time.time() - t0)}
    json.dump(res, open("adaptheteroTF_result.json", "w"), indent=1)
    print(f"seed {sd} DONE {res[f'seed{sd}']}", flush=True)
print("ALL DONE", flush=True)
