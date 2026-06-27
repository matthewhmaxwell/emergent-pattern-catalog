"""Face #4: REGULATION / homeostasis -- hold a variable at a setpoint (0) against active random
disturbances, in a system with INERTIA (velocity/momentum). The agent observes only the current
position (error bucket), NOT the velocity. A purely reactive controller (push opposite to current
error) overshoots and oscillates because of momentum; to damp it must estimate the trend, which
needs remembering recent error -> internal state.

Question: (1) does regulation REQUIRE state (sharp, like navigation/memory) or is it reactive
(like delayed gratification)? (2) if it needs state, is the state-use a distinct CONTROL/feedback
mode, or does it fold into storage (hold recent error)?
"""
import random

NIN = 9
ACT = [-1, 0, 1]


def quant(x):
    return min(8, max(0, int((x + 1.5) / 3.0 * 9)))


def random_rule(r, nstates):
    return [(r.randrange(3), r.randrange(nstates)) for _ in range(NIN * nstates)]


def mutate(rule, r, nstates, rate=0.05):
    c = rule[:]
    for i in range(len(c)):
        if r.random() < rate:
            c[i] = (r.randrange(3), r.randrange(nstates))
    return c


def episode(rule, nstates, r, T=40, sigma=0.25, damping=0.85, force=0.18, dist="noise"):
    x = r.uniform(-1, 1); v = 0.0; s = 0; sc = 0.0
    A = r.uniform(0.2, 0.5); per = r.uniform(5, 12); ph = r.uniform(0, 6.28)   # structured disturbance params
    for t in range(T):
        a, s = rule[quant(x) * nstates + s]; s %= nstates
        d = r.gauss(0, sigma) if dist == "noise" else A * __import__("math").sin(2 * 3.14159 * t / per + ph)
        v = damping * v + force * ACT[a] + d
        x = max(-1.5, min(1.5, x + v))
        sc += max(0.0, 1.0 - abs(x))                          # 1 at setpoint, 0 at the rails
    return sc / T


def fitness(rule, nstates, r, trials=16, **kw):
    return sum(episode(rule, nstates, r, **kw) for _ in range(trials)) / trials


def evolve(nstates, dist, gens=60, pop=120, mu=18):
    r = random.Random(321 + nstates + (0 if dist == "noise" else 99))
    P = [random_rule(random.Random(i + nstates * 13), nstates) for i in range(pop)]
    parents = P[:mu]
    for g in range(gens):
        rg = random.Random(8000 + g)
        scored = sorted(((fitness(ru, nstates, rg, dist=dist), ru) for ru in P), key=lambda x: -x[0])
        parents = [ru for _, ru in scored[:mu]]
        P = parents + [mutate(r.choice(parents), r, nstates) for _ in range(pop - mu - 5)] + \
            [random_rule(r, nstates) for _ in range(5)]
    return parents[0]


for dist in ("noise", "sine(structured/predictable)"):
    dk = "noise" if dist == "noise" else "sine"
    print(f"\n--- disturbance = {dist} ---")
    for ns in (1, 8):
        best = evolve(ns, dk)
        ho = sum(fitness(best, ns, random.Random(700 + k), trials=80, T=60, dist=dk) for k in range(5)) / 5
        tag = "MEMORYLESS (position only)" if ns == 1 else f"STATEFUL ({ns} states, can track/anticipate)"
        print(f"  {tag}: held-out regulation score = {ho:.2f}")
