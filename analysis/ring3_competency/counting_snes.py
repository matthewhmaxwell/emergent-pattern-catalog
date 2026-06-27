"""Decisive cheap test (user fork: 'stronger optimizer first'): is the counting ceiling the OPTIMIZER
or deeper? Swap vanilla OpenAI-ES for SNES (Separable Natural Evolution Strategies -- per-coordinate
adaptive variance + natural gradient + rank-based fitness shaping; the CMA-family choice that scales
to ~1000 params where full CMA-ES is impractical) and enlarge the RNN (H=32). Same exact-counting task
(collect EXACTLY k, then goal; supply > k; dense |count-k| shaping). If SNES makes a genuine
threshold-stopping counter emerge where OpenAI-ES failed (0.10), the ceiling was the optimizer (cheap
lift). If SNES also fails, that is decisive evidence the ceiling is deeper -> substrate escalation.

Run: python3 counting_snes.py [--k 3] [--supply 5] [--iters 400] [--H 32]
"""
import numpy as np, sys

a = sys.argv
H = int(a[a.index("--H") + 1]) if "--H" in a else 32
N = 9; OBS = 2 + 2 + 1; NACT = 5
DIRS = [(0, 0), (1, 0), (-1, 0), (0, 1), (0, -1)]


def nparams(): return H * OBS + H * H + H + NACT * H + NACT


def unpack(th):
    i = 0
    Wx = th[i:i + H * OBS].reshape(H, OBS); i += H * OBS
    Wh = th[i:i + H * H].reshape(H, H); i += H * H
    b = th[i:i + H]; i += H
    Wo = th[i:i + NACT * H].reshape(NACT, H); i += NACT * H
    bo = th[i:i + NACT]
    return Wx, Wh, b, Wo, bo


def layout(r, supply):
    cells = [(x, y) for x in range(N) for y in range(N)]
    idx = r.permutation(len(cells))
    return [[cells[idx[j]], True] for j in range(supply)], cells[idx[supply]], cells[idx[supply + 1]]


def sgn(p, q): return (q[0] > p[0]) - (q[0] < p[0]), (q[1] > p[1]) - (q[1] < p[1])


def episode(th, r, k, supply, STEPS=60):
    Wx, Wh, b, Wo, bo = unpack(th)
    objs, goal, start = layout(r, supply)
    p = list(start); h = np.zeros(H); count = 0; just = 0.0; reached = 0
    for _ in range(STEPS):
        best, bd = None, 99
        for op, al in objs:
            if al:
                d = abs(op[0] - p[0]) + abs(op[1] - p[1])
                if d < bd: bd, best = d, op
        x = np.zeros(OBS)
        if best is not None:
            x[0], x[1] = sgn(p, best)
        x[2], x[3] = sgn(p, goal); x[4] = just
        h = np.tanh(Wx @ x + Wh @ h + b)
        ac = int(np.argmax(Wo @ h + bo))
        np_ = (p[0] + DIRS[ac][0], p[1] + DIRS[ac][1])
        if 0 <= np_[0] < N and 0 <= np_[1] < N: p = list(np_)
        just = 0.0
        for ob in objs:
            if ob[1] and tuple(p) == ob[0]:
                ob[1] = False; count += 1; just = 1.0
        if tuple(p) == goal:
            reached = 1; break
    return 0.5 * reached + 0.5 * (1 - min(1.0, abs(count - k) / k)), int(reached and count == k), count


def fitness(th, seed, k, supply, trials=10):
    return np.mean([episode(th, np.random.default_rng(seed * 97 + j), k, supply)[0] for j in range(trials)])


def exact_rate(th, seed, k, supply, trials=60):
    return np.mean([episode(th, np.random.default_rng(seed * 131 + j), k, supply)[1] for j in range(trials)])


def mean_count(th, seed, k, supply, trials=60):
    return np.mean([episode(th, np.random.default_rng(seed * 53 + j), k, supply)[2] for j in range(trials)])


def snes(k, supply, iters=400, seed=1):
    r = np.random.default_rng(seed)
    d = nparams(); lam = 4 + int(3 * np.log(d)); lam = max(lam, 40)
    mu = r.standard_normal(d) * 0.4; sigma = np.ones(d) * 0.5
    eta_mu, eta_sig = 1.0, (3 + np.log(d)) / (5 * np.sqrt(d))
    rk = np.arange(1, lam + 1)
    util = np.maximum(0, np.log(lam / 2 + 1) - np.log(rk)); util = util / util.sum() - 1.0 / lam
    for it in range(iters):
        s = r.standard_normal((lam, d))
        z = mu + sigma * s
        F = np.array([fitness(z[i], 3000 + it, k, supply) for i in range(lam)])
        order = np.argsort(-F)                                  # best first
        u = np.empty(lam); u[order] = util                     # utility by rank
        mu = mu + eta_mu * sigma * (u[:, None] * s).sum(0)
        sigma = sigma * np.exp(eta_sig / 2 * (u[:, None] * (s ** 2 - 1)).sum(0))
    return mu


def openai_es(k, supply, iters=400, npop=50, sg=0.1, lr=0.05, seed=2):
    r = np.random.default_rng(seed); th = r.standard_normal(nparams()) * 0.4
    for it in range(iters):
        eps = r.standard_normal((npop, th.size))
        F = np.array([fitness(th + sg * e, 3000 + it, k, supply) for e in eps])
        A = (F - F.mean()) / (F.std() + 1e-8)
        th = th + lr / (npop * sg) * (eps.T @ A)
    return th


if __name__ == "__main__":
    k = int(a[a.index("--k") + 1]) if "--k" in a else 3
    supply = int(a[a.index("--supply") + 1]) if "--supply" in a else 5
    iters = int(a[a.index("--iters") + 1]) if "--iters" in a else 400
    print(f"counting (exact k={k}, supply={supply}), RNN H={H}, iters={iters} -- SNES vs OpenAI-ES")
    for name, fn in [("SNES      ", snes), ("OpenAI-ES ", openai_es)]:
        th = fn(k, supply, iters=iters)
        er = float(exact_rate(th, 7, k, supply))
        print(f"  {name}: exact-k {er:.2f}   mean count {mean_count(th,7,k,supply):.2f} (target {k})   "
              f"double-supply mean count {mean_count(th,11,k,supply*2):.2f}")
    print("\nVERDICT: a genuine counter = exact-k >= 0.6 AND mean count ~k AND stays ~k under doubled supply.")
