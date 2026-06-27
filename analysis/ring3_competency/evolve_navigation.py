"""Directed search for a DISCOVERED competency (v2: domain-randomized).

v1 over-fit a FIXED training set (train 0.92, held-out fail) -- the same Goodhart that bit seed
2959. Fix: every generation the whole population is scored on FRESH RANDOM conditions (random
starts + random single barriers/traps). There is no fixed set to memorize, so the only way to
keep scoring high across generations is to actually generalize. Then verify on a HELD-OUT bar the
training distribution never produces: specific corner starts + an unseen multi-wall serpentine.
"""
import random
from stateful_funnel import run_fsm, NCODES, NSTATES, RES, N


def random_rule(r):
    return [(r.randrange(5), r.randrange(NSTATES)) for _ in range(NCODES * NSTATES)]


def mutate(rule, r, rate=0.03):
    c = rule[:]
    for i in range(len(c)):
        if r.random() < rate:
            c[i] = (r.randrange(5), r.randrange(NSTATES))
    return c


def dist(p):
    return abs(RES[0] - p[0]) + abs(RES[1] - p[1])


def rand_conditions(r, K=16):
    """Fresh random (start, walls) each generation: open, or a random barrier/trap with a random
    gap. Single-wall obstacles only -- multi-wall mazes are reserved for held-out."""
    conds = []
    for _ in range(K):
        s = (r.randrange(N), r.randrange(N))
        if s == RES:
            s = (0, 0)
        kind = r.random()
        if kind < 0.45:
            w = set()
        else:
            col = r.randrange(3, N - 1)
            gap = r.randrange(0, N - 2)
            w = {(col, y) for y in range(0, N) if not (gap <= y <= gap + 1)}
        conds.append((s, w))
    return conds


def fitness(rule, conds):
    tot = 0.0
    for s, w in conds:
        tr, ok = run_fsm(rule, s, RES, w, steps=150)
        if ok:
            tot += 1.0
        else:
            tot += max(0.0, 1.0 - dist(tr[-1]) / (dist(s) or 1))
    return tot / len(conds)


HELDOUT_STARTS = [(0, 0), (12, 0), (12, 6), (5, 0), (11, 2)]
SERP = {(4, y) for y in range(0, N - 2)} | {(8, y) for y in range(3, N)}
def barrier():
    return {(RES[0] - 2, y) for y in range(0, N - 4)}
def trap():
    return {(RES[0] - 2, y) for y in range(2, N)}


def heldout_signature(rule):
    rch = sum(run_fsm(rule, s, RES, set(), 220)[1] for s in HELDOUT_STARTS) / len(HELDOUT_STARTS)
    rer = sum(run_fsm(rule, s, RES, barrier(), 220)[1] for s in HELDOUT_STARTS) / len(HELDOUT_STARTS)
    awy = sum(run_fsm(rule, s, RES, trap(), 240)[1] for s in [(0, 0), (5, 0), (11, 2)]) / 3
    gen = sum(run_fsm(rule, s, RES, SERP, 320)[1] for s in [(0, 0), (12, 0), (5, 0)]) / 3
    return dict(reach=round(rch, 2), reroute=round(rer, 2), away=round(awy, 2), generalize=round(gen, 2))


POP, MU, GENS = 140, 24, 90
r = random.Random(123)
pop = [random_rule(random.Random(i)) for i in range(POP)]
for g in range(GENS):
    conds = rand_conditions(r)                                 # fresh each generation
    scored = sorted(((fitness(ru, conds), ru) for ru in pop), key=lambda x: -x[0])
    if g % 15 == 0 or g == GENS - 1:
        print(f"gen {g:>2}: best on this gen's random conds {scored[0][0]:.3f}")
    parents = [ru for _, ru in scored[:MU]]
    pop = parents + [mutate(r.choice(parents), r) for _ in range(POP - MU - 6)] + [random_rule(r) for _ in range(6)]

# pick the best GENERALIZER among final survivors by held-out, report it honestly
finals = parents[:MU]
sigs = [(heldout_signature(ru), ru) for ru in finals]
def score_sig(s):
    return s["reach"] + s["reroute"] + s["away"] + s["generalize"]
sig, bestrule = max(sigs, key=lambda x: score_sig(x[0]))
print(f"\nbest final survivor HELD-OUT signature (disjoint corners + unseen serpentine): {sig}")
core = sig["reroute"] >= 0.6 and sig["away"] >= 0.6
broad = sig["reach"] >= 0.8 and core and sig["generalize"] >= 0.6
print("VERDICT:", "DISCOVERED COMPETENCY (clears broad held-out bar)" if broad
      else ("PARTIAL (core re-routing, not full general navigation)" if core
            else "no robust competency yet (improved but does not clear the broad bar)"))

# ---- DEBUNK: hammer the best survivor with fresh NOVEL multi-wall mazes never seen ----
from stateful_funnel import dfs_reaches
rr = random.Random(999)
def rand_maze():
    w = set()
    for _ in range(rr.randrange(2, 4)):                        # 2-3 walls (out of training distribution)
        col = rr.randrange(2, N - 1); gap = rr.randrange(0, N - 2)
        w |= {(col, y) for y in range(N) if not (gap <= y <= gap + 1)}
    return w
solv = reached = 0
for _ in range(60):
    w = rand_maze()
    for _ in range(3):
        s = (rr.randrange(N), rr.randrange(N))
        if s == RES or s in w or not dfs_reaches(s, RES, w):
            continue
        solv += 1
        if run_fsm(bestrule, s, RES, w, 400)[1]:
            reached += 1
print(f"DEBUNK on {solv} SOLVABLE novel multi-wall mazes (never seen): reached {reached}/{solv} = {reached/max(1,solv):.2f}")
print("  ->", "CONFIRMED general navigation competency" if solv and reached/solv >= 0.8
      else "serpentine pass was narrow -- not robustly general")
