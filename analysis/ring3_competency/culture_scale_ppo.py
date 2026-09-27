"""Multi-scale / temporal-axis experiment (paper §8.2): the across-generation scale as a characterized axis.

openculture_ppo.py demonstrated the cumulative-culture ratchet at ONE recipe length (L=12): cumulative ~11.7/12
vs no-inheritance ~0.4/12. Track-2 asks to move the temporal-scale axis from asserted to DEMONSTRATED. This
sweeps recipe length L and shows:
  (1) the across-generation buffer is a genuine competency SCALE: cumulative culture reliably reaches ~L at EVERY
      length, while the buffer-ablation (no-inheritance) stays flat and low -- so the ratchet advantage grows
      with task complexity L (an unbounded scale vs a bounded one).
  (2) the buffer-ablation IS the instrument's collapse-ablation, applied at the CROSS-GENERATIONAL scale:
      remove the inherited artifact -> the competency collapses. Same interventional method, new temporal scale.
This is the temporal analog of the cost-ladder dial: there, the cheapest coordination MEDIUM is recruited only
when asymmetry forces it; here, the higher temporal SCALE (across-generation accumulation) is the only thing
that reaches recipes a single lifetime cannot -- and its necessity is proven by the collapse when ablated.

Honest bound: mechanism = imitation (#19) + frontier innovation + a persistent artifact (#15), across
generations -> cumulative cultural evolution (Tomasello / Boyd-Richerson / Kirby). 0 new-to-science; the
contribution is the interventional characterization of the axis, NOT the mechanism. Differentiate from Watson &
Szathmary "same mechanism across timescales": the accumulation competency has NO within-lifetime analog for an
unshown recipe -- the scale itself is load-bearing.

Run (VPS epc-venv, cheap ~minutes): DIAL_THREADS=4 PYTHONDONTWRITEBYTECODE=1 nice -n 10 \
  /home/matthewhmaxwell/epc-venv/bin/python culture_scale_ppo.py --Ls 2,4,6,8,12,16,20 --seeds 3
"""
import numpy as np, sys, json, os, torch, torch.nn as nn

M = 4; H = 96; B = 256
torch.set_num_threads(int(os.environ.get("DIAL_THREADS", "4")))


class Pol(nn.Module):
    def __init__(self, idim):
        super().__init__()
        self.body = nn.Sequential(nn.Linear(idim, H), nn.Tanh(), nn.Linear(H, H), nn.Tanh())
        self.pi = nn.Linear(H, M)

    def forward(self, x): return self.pi(self.body(x))


def make_obs(recipe, p, L):
    b = recipe.shape[0]; ID = L + 1 + M; out = np.zeros((b, L, ID), np.float32)
    pos = np.arange(L)
    out[:, pos, pos] = 1.0                                   # position one-hot
    known = pos[None, :] < p[:, None]                        # (b, L)
    out[:, :, L] = known.astype(np.float32)
    bb, ii = np.where(known)
    out[bb, ii, L + 1 + recipe[bb, ii]] = 1.0                # inherited cultural symbol where known
    return out.reshape(b * L, ID)


def prefix_len(sym, recipe, L):
    correct = (sym == recipe); pl = np.zeros(sym.shape[0], int)
    for b in range(sym.shape[0]):
        k = 0
        while k < L and correct[b, k]: k += 1
        pl[b] = k
    return pl


def run(gens, inherit, L, seed=0):
    torch.manual_seed(seed); rng = np.random.default_rng(seed)
    net = Pol(L + 1 + M); opt = torch.optim.Adam(net.parameters(), lr=3e-3)
    recipe = rng.integers(0, M, size=(B, L))
    p = np.zeros(B, int)
    for g in range(gens):
        obs = make_obs(recipe, p if inherit else np.zeros(B, int), L)
        logits = net(torch.from_numpy(obs)); d = torch.distributions.Categorical(logits=logits)
        a = d.sample(); sym = a.numpy().reshape(B, L)
        pl = prefix_len(sym, recipe, L)
        rew = np.zeros((B, L), np.float32)
        for b in range(B): rew[b, :pl[b]] = 1.0
        R = torch.from_numpy(rew.reshape(B * L)); adv = R - R.mean()
        loss = -(d.log_prob(a) * adv).mean() - 0.03 * d.entropy().mean()
        opt.zero_grad(); loss.backward(); opt.step()
        if inherit: p = np.maximum(p, pl)
    if inherit:
        return float(p.mean())
    # no-inheritance: greedy final prefix a single lifetime achieves (buffer ablated)
    fin = net(torch.from_numpy(make_obs(recipe, np.zeros(B, int), L))).argmax(1).numpy().reshape(B, L)
    return float(prefix_len(fin, recipe, L).mean())


if __name__ == "__main__":
    a = sys.argv
    Ls = [int(x) for x in (a[a.index("--Ls") + 1].split(",") if "--Ls" in a else ["2", "4", "6", "8", "12", "16", "20"])]
    seeds = int(a[a.index("--seeds") + 1]) if "--seeds" in a else 3
    here = __file__.rsplit("/", 1)[0]; out = f"{here}/culture_scale_sweep.json"
    prior = {r["L"]: r for r in json.load(open(out))} if os.path.exists(out) else {}
    print(f"TEMPORAL-SCALE sweep: Ls={Ls} seeds={seeds} (resuming {sorted(prior)})", flush=True)
    results = []
    for L in Ls:
        if L in prior:
            results.append(prior[L]); print(f"L={L}: (cached)", flush=True); continue
        gens = max(80, 25 * L)
        cum = np.mean([run(gens, True, L, seed=s) for s in range(seeds)])
        noi = np.mean([run(gens, False, L, seed=100 + s) for s in range(seeds)])
        row = {"L": L, "gens": gens, "cumulative": round(cum, 2), "cumulative_frac": round(cum / L, 3),
               "noinherit": round(noi, 2), "noinherit_frac": round(noi / L, 3), "ratchet_gap": round(cum - noi, 2)}
        results.append(row); json.dump(results, open(out, "w"), indent=1)
        print(f"L={L:>2}: cumulative {cum:.2f}/{L} ({cum/L:.0%}) | no-inherit {noi:.2f}/{L} ({noi/L:.0%}) "
              f"| RATCHET GAP {cum - noi:+.2f}", flush=True)
    json.dump(results, open(out, "w"), indent=1); print("SCALEDONE saved culture_scale_sweep.json", flush=True)
