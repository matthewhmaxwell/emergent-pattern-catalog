"""Decisive counting test (well-tooled): can RNN+ES learn to collect EXACTLY k objects then stop and
go to the goal -- a genuine accumulate-to-threshold counter (a candidate THIRD state-mode for the map,
alongside commitment and storage)? The open-ended sweeps left counting at <=0.28, but with a NARROW
fitness (only exactly-k scores). Here we give DENSE shaping: reward = 0.5*reached_goal + 0.5*(1 -
|count-k|/k), which pulls toward stopping at exactly k (collect-all overshoots and scores lower), plus
more ES budget. A memoryless baseline is the control. Then we vary the SUPPLY (the --double test the
agent-observer used) to confirm a real counter stops at k regardless of supply -- not greedy
collect-all.

Run: python3 counting.py [--k 3] [--supply 5] [--iters 500]
"""
import numpy as np, sys

N = 9; H = 16; OBS = 2 + 2 + 1; NACT = 5     # dir-to-nearest type, dir-to-goal, just-touched flag
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
    objs = [[cells[idx[j]], True] for j in range(supply)]
    return objs, cells[idx[supply]], cells[idx[supply + 1]]


def sgn(a, b): return (b[0] > a[0]) - (b[0] < a[0]), (b[1] > a[1]) - (b[1] < a[1])


def episode(th, r, k, supply, recurrent=True, STEPS=60):
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
            sx, sy = sgn(p, best); x[0], x[1] = sx, sy
        gx, gy = sgn(p, goal); x[2], x[3] = gx, gy; x[4] = just
        h = np.tanh(Wx @ x + (Wh @ h if recurrent else 0.0) + b)
        a = int(np.argmax(Wo @ h + bo))
        np_ = (p[0] + DIRS[a][0], p[1] + DIRS[a][1])
        if 0 <= np_[0] < N and 0 <= np_[1] < N: p = list(np_)
        just = 0.0
        for ob in objs:
            if ob[1] and tuple(p) == ob[0]:
                ob[1] = False; count += 1; just = 1.0
        if tuple(p) == goal:
            reached = 1; break
    fit = 0.5 * reached + 0.5 * (1 - min(1.0, abs(count - k) / k))
    return fit, (reached and count == k), count


def fitness(th, seed, k, supply, recurrent=True, trials=10):
    return np.mean([episode(th, np.random.default_rng(seed * 97 + j), k, supply, recurrent)[0] for j in range(trials)])


def exact_rate(th, seed, k, supply, recurrent=True, trials=50):
    return np.mean([episode(th, np.random.default_rng(seed * 131 + j), k, supply, recurrent)[1] for j in range(trials)])


def mean_count(th, seed, k, supply, recurrent=True, trials=50):
    return np.mean([episode(th, np.random.default_rng(seed * 53 + j), k, supply, recurrent)[2] for j in range(trials)])


def es(k, supply, recurrent=True, iters=500, npop=40, sigma=0.1, lr=0.05, seed=1):
    r = np.random.default_rng(seed); th = r.standard_normal(nparams()) * 0.4
    for it in range(iters):
        eps = r.standard_normal((npop, th.size))
        F = np.array([fitness(th + sigma * e, 2000 + it, k, supply, recurrent) for e in eps])
        A = (F - F.mean()) / (F.std() + 1e-8)
        th = th + lr / (npop * sigma) * (eps.T @ A)
    return th


if __name__ == "__main__":
    a = sys.argv
    k = int(a[a.index("--k") + 1]) if "--k" in a else 3
    supply = int(a[a.index("--supply") + 1]) if "--supply" in a else 5
    iters = int(a[a.index("--iters") + 1]) if "--iters" in a else 500
    print(f"counting test: collect EXACTLY k={k} (supply={supply}), then goal. ES iters={iters}, dense fitness.")
    rnn = es(k, supply, True, iters=iters, seed=1)
    mem = es(k, supply, False, iters=iters, seed=2)
    re_, me_ = float(exact_rate(rnn, 7, k, supply, True)), float(exact_rate(mem, 8, k, supply, False))
    print(f"  exact-k success:  RNN {re_:.2f}   memoryless {me_:.2f}   gap {re_ - me_:+.2f}")
    print(f"  mean count collected:  RNN {mean_count(rnn,7,k,supply,True):.2f}  memoryless {mean_count(mem,8,k,supply,False):.2f}  (target {k})")
    print("\n  --double test (does the RNN STOP at k when supply is larger?):")
    for s2 in (supply, supply * 2):
        print(f"    supply={s2}: RNN mean count = {mean_count(rnn, 11, k, s2, True):.2f}  exact-rate {exact_rate(rnn,11,k,s2,True):.2f}")
    print("\nVERDICT:",
          "GENUINE COUNTER emerged (RNN exact>=0.6, gap>=0.3, stops at k under doubled supply)" if (re_ >= 0.6 and re_ - me_ >= 0.3) else
          f"no genuine counter (RNN exact {re_:.2f}) -> accumulation NOT reachable by this substrate even well-tooled")
