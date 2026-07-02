"""OPEN-ENDED hunt -- NEW TERRITORY (cumulative-culture shot): a CROSS-GENERATIONAL RATCHET.

Every catalogued competency + the whole coordination cascade is WITHIN-lifetime. Cumulative culture is the one
genuinely untested regime: the collective reaches solutions NO SINGLE GENERATION could, by inheriting and
building on accumulated culture (Tomasello / Boyd-Richerson / Kirby iterated learning).

World: a hidden RECIPE (length-L sequence over M symbols), one per parallel 'lineage'. Culture = the known-
correct PREFIX length p (persists across generations = the inherited artifact). Each generation the agents see
the prefix (positions < p, with their correct symbols) and must output the whole length-L sequence; reward =
length of the correct leading prefix. When the frontier symbol (position p) is discovered, the culture advances
to p+1. A persistent policy learns the GENERAL cultural strategy (copy the inherited prefix + explore the
frontier); the recipe is never memorized because B parallel lineages each have a different random recipe.

The RATCHET test (three conditions):
  - CUMULATIVE (buffer accumulates) -> p should climb toward L.
  - NO-INHERITANCE (buffer always empty, p=0) -> must discover the whole recipe in ONE lifetime -> EXPONENTIALLY
    hard -> plateaus at a low prefix length. The gap = the ratchet.
Mechanistically this REDUCES to imitation (#19, copy the prefix) + innovation (frontier exploration) + a
persistent artifact (#15) -- all catalogued -> KNOWN (cumulative cultural evolution) -> 0 new-to-science. But it
is the FIRST cross-generational (not within-episode) competency, and the RATCHET is a genuine system-level
property. Honest question: does a real ratchet emerge, and is anything about it OUTSIDE the closed set?

Run on VPS epc-venv: python openculture_ppo.py [--gens 60]
"""
import numpy as np, sys, json, torch, torch.nn as nn

L = 12; M = 4; H = 96; B = 256
ODIM = L + 1 + M                                     # position one-hot(L), known-flag(1), cultural symbol one-hot(M)
torch.set_num_threads(8)


class Pol(nn.Module):
    def __init__(self):
        super().__init__()
        self.body = nn.Sequential(nn.Linear(ODIM, H), nn.Tanh(), nn.Linear(H, H), nn.Tanh())
        self.pi = nn.Linear(H, M)

    def forward(self, x): return self.pi(self.body(x))


def make_obs(recipe, p):
    """obs for every (lineage b, position i): position, known-flag, cultural symbol if known."""
    b = recipe.shape[0]; out = np.zeros((b, L, ODIM), np.float32)
    for i in range(L):
        out[:, i, i] = 1.0
        known = i < p
        out[:, i, L] = known.astype(np.float32)
        for bb in range(b):
            if known[bb]: out[bb, i, L + 1 + recipe[bb, i]] = 1.0
    return out.reshape(b * L, ODIM)


def prefix_len(sym, recipe):
    correct = (sym == recipe)                                    # (B, L)
    pl = np.zeros(sym.shape[0], int)
    for b in range(sym.shape[0]):
        k = 0
        while k < L and correct[b, k]: k += 1
        pl[b] = k
    return pl


def run(gens, inherit, seed=0):
    torch.manual_seed(seed); rng = np.random.default_rng(seed)
    net = Pol(); opt = torch.optim.Adam(net.parameters(), lr=3e-3)
    recipe = rng.integers(0, M, size=(B, L))                     # one random recipe per lineage
    p = np.zeros(B, int)                                         # culture prefix length per lineage
    traj = []
    for g in range(gens):
        obs = make_obs(recipe, p if inherit else np.zeros(B, int))
        logits = net(torch.from_numpy(obs))
        d = torch.distributions.Categorical(logits=logits)
        a = d.sample(); sym = a.numpy().reshape(B, L)
        pl = prefix_len(sym, recipe)
        # per-position reward: +1 if position i is within the correct prefix
        rew = np.zeros((B, L), np.float32)
        for b in range(B): rew[b, :pl[b]] = 1.0
        R = torch.from_numpy(rew.reshape(B * L))
        adv = R - R.mean()
        loss = -(d.log_prob(a) * adv).mean() - 0.03 * d.entropy().mean()
        opt.zero_grad(); loss.backward(); opt.step()
        if inherit:
            p = np.maximum(p, pl)                                # culture accumulates the discovered prefix
        traj.append(float((p if inherit else pl).mean()))
    final = float((p if inherit else prefix_len(net(torch.from_numpy(make_obs(recipe, np.zeros(B, int)))).argmax(1).numpy().reshape(B, L), recipe)).mean())
    return traj, (float(p.mean()) if inherit else traj[-1])


if __name__ == "__main__":
    a = sys.argv
    gens = int(a[a.index("--gens") + 1]) if "--gens" in a else 60
    here = __file__.rsplit("/", 1)[0]
    print(f"CUMULATIVE-CULTURE ratchet: recipe length L={L}, M={M} symbols, {gens} generations", flush=True)
    tc, fc = run(gens, inherit=True, seed=0)
    tn, fn = run(gens, inherit=False, seed=0)
    print(f"  CUMULATIVE culture prefix (mean/{L}): gen0 {tc[0]:.2f} -> gen{gens-1} {tc[-1]:.2f}", flush=True)
    print(f"  NO-INHERITANCE prefix   (mean/{L}): gen0 {tn[0]:.2f} -> gen{gens-1} {tn[-1]:.2f}", flush=True)
    ratchet = tc[-1] - tn[-1]
    print(f"\nFINAL: cumulative {tc[-1]:.2f} | no-inheritance {tn[-1]:.2f} (of L={L})", flush=True)
    print(f"  RATCHET gap (cumulative - no-inheritance) = {ratchet:+.2f}", flush=True)
    genuine = tc[-1] >= 0.7 * L and ratchet >= 0.3 * L
    print(f"VERDICT: {'GENUINE CROSS-GENERATIONAL RATCHET (cumulative culture) -- reaches solutions no single generation could' if genuine else 'no clear ratchet'} "
          f"(cumulative reaches {tc[-1]/L:.0%} of L; single-generation plateaus at {tn[-1]/L:.0%})", flush=True)
    json.dump({"cumulative_final": round(tc[-1], 2), "noinherit_final": round(tn[-1], 2), "L": L,
               "ratchet_gap": round(ratchet, 2), "genuine": bool(genuine),
               "cumulative_traj": [round(x, 2) for x in tc]}, open(f"{here}/openculture_result.json", "w"))
