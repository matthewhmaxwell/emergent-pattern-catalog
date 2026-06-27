"""Substrate escalation, step 1: richer agent (a small recurrent network) + stronger search
(evolution strategies / OpenAI-ES on continuous weights), tested on the SAME prerequisite task the
FSM+mutation under-solved in rung 2 (visit A before goal G counts; once you leave A a reactive agent
oscillates, so it needs memory -> the RNN's hidden state). If the RNN+ES cleanly solves it where the
FSM did not, the rung-2 failure was a tooling ceiling, and the escalation is the right substrate for
the open-ended off-map hunt.
"""
import numpy as np

N = 9
DIRS = [(0, 0), (1, 0), (-1, 0), (0, 1), (0, -1)]
H = 16
rng0 = np.random.default_rng(0)


def init_theta(r):
    Wx = r.standard_normal((H, 4)) * 0.5
    Wh = r.standard_normal((H, H)) * 0.5
    b = np.zeros(H)
    Wo = r.standard_normal((5, H)) * 0.5
    bo = np.zeros(5)
    return np.concatenate([Wx.ravel(), Wh.ravel(), b, Wo.ravel(), bo])


def unpack(th):
    i = 0
    Wx = th[i:i + H * 4].reshape(H, 4); i += H * 4
    Wh = th[i:i + H * H].reshape(H, H); i += H * H
    b = th[i:i + H]; i += H
    Wo = th[i:i + 5 * H].reshape(5, H); i += 5 * H
    bo = th[i:i + 5]
    return Wx, Wh, b, Wo, bo


def sgn(a, b):
    return np.array([(b[0] > a[0]) - (b[0] < a[0]), (b[1] > a[1]) - (b[1] < a[1])], float)


def episode(th, start, A, G, T=50):
    Wx, Wh, b, Wo, bo = unpack(th)
    p = list(start); h = np.zeros(H); done = 0
    for _ in range(T):
        x = np.concatenate([sgn(p, A), sgn(p, G)])
        h = np.tanh(Wx @ x + Wh @ h + b)
        a = int(np.argmax(Wo @ h + bo))
        np_ = (p[0] + DIRS[a][0], p[1] + DIRS[a][1])
        if 0 <= np_[0] < N and 0 <= np_[1] < N:
            p = list(np_)
        if done == 0 and tuple(p) == A:
            done = 1
        if done == 1 and tuple(p) == G:
            return 2.0
    return float(done)


def cfg(r):
    c = [(x, y) for x in range(N) for y in range(N)]; r.shuffle(c)
    return c[0], c[1], c[2]


def fitness(th, r, trials=12):
    return np.mean([episode(th, *cfg(r)) for _ in range(trials)])


def es(iters=200, npop=40, sigma=0.1, lr=0.05, seed=1):
    r = np.random.default_rng(seed)
    th = init_theta(r)
    for it in range(iters):
        eps = r.standard_normal((npop, th.size))
        fr = np.random.default_rng(1000 + it)
        F = np.array([fitness(th + sigma * e, np.random.default_rng(2000 + it)) for e in eps])
        A = (F - F.mean()) / (F.std() + 1e-8)
        th = th + lr / (npop * sigma) * (eps.T @ A)
        if it % 40 == 0 or it == iters - 1:
            print(f"  iter {it:>3}: mean prereq progress (0-2) = {fitness(th, np.random.default_rng(7)):.2f}")
    return th


print("escalation step 1: RNN (16 hidden) + evolution strategies on the prerequisite task")
th = es()
final = np.mean([fitness(th, np.random.default_rng(900 + k), trials=40) for k in range(4)])
print(f"\nfinal held-out mean progress (0=no A, 1=A only, 2=A then G): {final:.2f}")
print("VERDICT:", "RNN+ES SOLVES it (>=1.8) -> richer substrate fixes the rung-2 tooling ceiling"
      if final >= 1.8 else f"still under-solved ({final:.2f}) -> deeper ceiling, not just FSM tooling")
