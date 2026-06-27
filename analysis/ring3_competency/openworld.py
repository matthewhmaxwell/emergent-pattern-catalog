"""Open-ended off-map hunt, LAYER 1: an open-ended environment GENERATOR (reward RULE sampled from a
compositional grammar, NOT hand-designed) + the validated RNN+ES learner + an MCC filter that keeps
only environments demanding a competency a mundane (memoryless) baseline cannot supply.

Why this design: the moonshot showed the bottleneck is environment design, and rung-2/escalation
showed the FSM was the wrong substrate. So: sample the reward rule from a grammar (count / order /
collect literals, conjoined) so we do NOT bake in which competency is demanded; train the validated
RNN with evolution strategies; and contrast against a memoryless agent of the same architecture. An
environment is kept ("demanding") only if the RNN succeeds and the memoryless baseline fails -- i.e.
the environment genuinely requires internal state / accumulation. The demanding environments + their
agents are then handed to the agent-observer (Layer 2) to NAME the competency and classify it vs the
map {navigation, memory, delayed-gratification, regulation, combinations} or flag it OFF-MAP.

Run: python3 openworld.py [--envs 8] [--iters 100] [--seed 7]
Writes openworld_sweep.json (+ openworld_<env>.npy weights for the demanding ones).
"""
import numpy as np, json, sys

N = 9; T = 3; PER = 5; STEPS = 55; H = 16   # PER>threshold so "collect-all" != "count": forces a real counter
OBS = T * 2 + 2 + T          # dir-to-nearest-remaining per type, dir-to-goal, last-touched one-hot
NACT = 5
DIRS = [(0, 0), (1, 0), (-1, 0), (0, 1), (0, -1)]


# ---------- open-ended reward-rule grammar ----------
def sample_predicate(r, mode="open"):
    if mode == "count":                                     # focused: single count/exact literal
        t = int(r.integers(0, T)); k = int(r.integers(2, 4))
        return [("exact", t, k)] if r.random() < 0.5 else [("count", t, k)]
    lits = []
    for _ in range(int(r.integers(1, 3))):                  # 1 or 2 conjoined literals
        kind = r.choice(["count", "exact", "order", "collect"])
        if kind == "count":
            lits.append(("count", int(r.integers(0, T)), int(r.integers(2, 4))))
        elif kind == "exact":                               # EXACTLY k -> over-collection FAILS -> demands a counter
            lits.append(("exact", int(r.integers(0, T)), int(r.integers(2, 4))))
        elif kind == "order":
            ab = r.choice(T, size=2, replace=False); lits.append(("order", int(ab[0]), int(ab[1])))
        else:
            m = int(r.integers(2, T + 1))
            lits.append(("collect", tuple(sorted(int(x) for x in r.choice(T, size=m, replace=False)))))
    out = []
    for l in lits:                                          # dedup identical literals
        if l not in out: out.append(l)
    return out


def pred_frac(lits, counts, first):
    ok = 0
    for lit in lits:
        if lit[0] == "count":
            ok += counts[lit[1]] >= lit[2]
        elif lit[0] == "exact":
            ok += counts[lit[1]] == lit[2]                  # exactly k -> over-collection breaks it
        elif lit[0] == "order":
            fa, fb = first[lit[1]], first[lit[2]]
            ok += (fa is not None and fb is not None and fa < fb)
        else:
            ok += all(counts[t] >= 1 for t in lit[1])
    return ok / len(lits)


def descr(lits):
    p = []
    for lit in lits:
        if lit[0] == "count": p.append(f"count(t{lit[1]})>={lit[2]}")
        elif lit[0] == "exact": p.append(f"count(t{lit[1]})=={lit[2]}")
        elif lit[0] == "order": p.append(f"t{lit[1]}-before-t{lit[2]}")
        else: p.append("collect(" + ",".join("t" + str(t) for t in lit[1]) + ")")
    return " AND ".join(p)


# ---------- RNN agent ----------
def nparams(): return H * OBS + H * H + H + NACT * H + NACT


def unpack(th):
    i = 0
    Wx = th[i:i + H * OBS].reshape(H, OBS); i += H * OBS
    Wh = th[i:i + H * H].reshape(H, H); i += H * H
    b = th[i:i + H]; i += H
    Wo = th[i:i + NACT * H].reshape(NACT, H); i += NACT * H
    bo = th[i:i + NACT]
    return Wx, Wh, b, Wo, bo


def layout(r):
    cells = [(x, y) for x in range(N) for y in range(N)]
    idx = r.permutation(len(cells)); k = 0
    objs = []
    for t in range(T):
        for _ in range(PER):
            objs.append([cells[idx[k]], t, True]); k += 1
    return objs, cells[idx[k]], cells[idx[k + 1]]           # objs, goal, start


def sgn(a, b):
    return (b[0] > a[0]) - (b[0] < a[0]), (b[1] > a[1]) - (b[1] < a[1])


def episode(th, r, lits, recurrent=True):
    Wx, Wh, b, Wo, bo = unpack(th)
    objs, goal, start = layout(r)
    p = list(start); h = np.zeros(H); counts = [0] * T; first = [None] * T; last = -1
    success = 0
    for step in range(STEPS):
        x = np.zeros(OBS); o = 0
        for t in range(T):
            best, bd = None, 99
            for op, ot, al in objs:
                if al and ot == t:
                    d = abs(op[0] - p[0]) + abs(op[1] - p[1])
                    if d < bd: bd, best = d, op
            if best is not None:
                sx, sy = sgn(p, best); x[o], x[o + 1] = sx, sy
            o += 2
        sx, sy = sgn(p, goal); x[o], x[o + 1] = sx, sy; o += 2
        if last >= 0: x[o + last] = 1.0
        h = np.tanh(Wx @ x + (Wh @ h if recurrent else 0.0) + b)
        a = int(np.argmax(Wo @ h + bo))
        np_ = (p[0] + DIRS[a][0], p[1] + DIRS[a][1])
        if 0 <= np_[0] < N and 0 <= np_[1] < N: p = list(np_)
        last = -1
        for ob in objs:
            if ob[2] and tuple(p) == ob[0]:
                ob[2] = False; counts[ob[1]] += 1
                if first[ob[1]] is None: first[ob[1]] = step
                last = ob[1]
        if tuple(p) == goal and pred_frac(lits, counts, first) == 1.0:
            success = 1; break
    return 0.5 * pred_frac(lits, counts, first) + 0.5 * success, success


def fitness(th, seed, lits, recurrent=True, trials=8):
    return np.mean([episode(th, np.random.default_rng(seed * 97 + k), lits, recurrent)[0] for k in range(trials)])


def success_rate(th, seed, lits, recurrent=True, trials=40):
    return np.mean([episode(th, np.random.default_rng(seed * 131 + k), lits, recurrent)[1] for k in range(trials)])


def es(lits, recurrent=True, iters=100, npop=30, sigma=0.1, lr=0.05, seed=1):
    r = np.random.default_rng(seed); th = r.standard_normal(nparams()) * 0.4
    for it in range(iters):
        eps = r.standard_normal((npop, th.size))
        F = np.array([fitness(th + sigma * e, 2000 + it, lits, recurrent) for e in eps])
        A = (F - F.mean()) / (F.std() + 1e-8)
        th = th + lr / (npop * sigma) * (eps.T @ A)
    return th


if __name__ == "__main__":
    a = sys.argv
    nenv = int(a[a.index("--envs") + 1]) if "--envs" in a else 8
    iters = int(a[a.index("--iters") + 1]) if "--iters" in a else 100
    rmeta = np.random.default_rng(int(a[a.index("--seed") + 1]) if "--seed" in a else 7)
    mode = a[a.index("--mode") + 1] if "--mode" in a else "open"
    prefix = a[a.index("--prefix") + 1] if "--prefix" in a else "openworld"
    here = __file__.rsplit("/", 1)[0]
    results = []
    print(f"open-ended sweep: {nenv} sampled environments, ES iters={iters}, mode={mode}, PER={PER}")
    for e in range(nenv):
        lits = sample_predicate(rmeta, mode=mode)
        rnn = es(lits, True, iters=iters, seed=10 + e)
        mem = es(lits, False, iters=iters, seed=100 + e)
        rs = float(success_rate(rnn, 500 + e, lits, True))
        ms = float(success_rate(mem, 600 + e, lits, False))
        dem = rs >= 0.6 and (rs - ms) >= 0.3
        if dem:
            np.save(f"{here}/{prefix}_env{e}.npy", rnn)
        results.append({"env": e, "pred": descr(lits), "lits": lits, "rnn": round(rs, 2),
                        "mem": round(ms, 2), "gap": round(rs - ms, 2), "demanding": dem})
        print(f"env {e}: {descr(lits):42s} RNN {rs:.2f}  mundane {ms:.2f}  gap {rs - ms:+.2f}  {'DEMANDING' if dem else ''}")
    dem = [x for x in results if x["demanding"]]
    print(f"\n{len(dem)}/{nenv} environments DEMAND competency (RNN>=0.6 AND gap>=0.3)")
    json.dump(results, open(f"{here}/{prefix}_sweep.json", "w"), indent=1)
    print(f"saved {prefix}_sweep.json" + (f" + weights for envs {[x['env'] for x in dem]}" if dem else ""))
