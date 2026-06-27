"""Third competency face: DELAYED GRATIFICATION (cash out at the peak of a rising-then-vanishing
reward of UNKNOWN height). Question: is this a THIRD state-mode, or does it fold into the two we
have (commitment / storage)?

Design that forces the question: the agent sees only the ABSOLUTE current reward value (a bucket
0-9). Height is randomized per episode, so no fixed threshold works (a low-height peak never
reaches a high threshold; a high threshold set low cashes out early on the rise). The ONLY way to
cash out at the peak under unknown height is to notice the value just turned DOWN -- which requires
remembering the recent value. So: memoryless should cap below optimal; stateful should win by
detecting the turn; and if it wins by tracking the recent value, that is STORAGE mode (folds into
memory), not a new mode.
"""
import random

T = 16
NIN = 10


def curve(r):
    p = r.randint(3, T - 3); h = r.uniform(0.4, 1.0); w = r.uniform(2.0, 4.0)
    return [max(0.0, h * (1 - abs(t - p) / w)) for t in range(T)], h


def heldout_curve(r):   # novel shape family: plateau-top bumps + wider, never trained
    p = r.randint(3, T - 4); h = r.uniform(0.4, 1.0); w = r.uniform(4.0, 6.0)
    return [max(0.0, min(h, h * 1.3 * (1 - abs(t - p) / w))) for t in range(T)], h


def random_rule(r, nstates):
    return [(r.randrange(2), r.randrange(nstates)) for _ in range(NIN * nstates)]


def mutate(rule, r, nstates, rate=0.05):
    c = rule[:]
    for i in range(len(c)):
        if r.random() < rate:
            c[i] = (r.randrange(2), r.randrange(nstates))
    return c


def episode(rule, c, h, nstates):
    s = 0
    for v in c:
        b = min(9, int(v * 9 + 1e-9))
        a, s = rule[b * nstates + s]; s %= nstates
        if a == 1:
            return v / h                                       # normalized payoff (1.0 = took at peak)
    return 0.0


def fitness(rule, r, nstates, trials=24, held=False):
    g = heldout_curve if held else curve
    tot = 0.0
    for _ in range(trials):
        c, h = g(r)
        tot += episode(rule, c, h, nstates)
    return tot / trials


def evolve(nstates, gens=60, pop=120, mu=18):
    r = random.Random(123 + nstates)
    P = [random_rule(random.Random(i + nstates * 7), nstates) for i in range(pop)]
    parents = P[:mu]
    for g in range(gens):
        rg = random.Random(5000 + g)
        scored = sorted(((fitness(ru, rg, nstates), ru) for ru in P), key=lambda x: -x[0])
        parents = [ru for _, ru in scored[:mu]]
        P = parents + [mutate(r.choice(parents), r, nstates) for _ in range(pop - mu - 5)] + \
            [random_rule(r, nstates) for _ in range(5)]
    return parents[0]


for ns in (1, 12):
    best = evolve(ns)
    ho = sum(fitness(best, random.Random(900 + k), ns, trials=80, held=True) for k in range(5)) / 5
    tag = "MEMORYLESS (no state)" if ns == 1 else f"STATEFUL ({ns} states, can track recent value)"
    print(f"{tag}: held-out normalized payoff (novel curves) = {ho:.2f}   (1.0 = always cash out at peak)")
