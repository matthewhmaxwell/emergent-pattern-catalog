"""Moonshot probe: a MULTI-AFFORDANCE world where diverse competencies can emerge, so the search
is not pre-aimed at one named face. An agent with depleting energy must survive by foraging
regrowing food and avoiding hazards. Different survivors may evolve different strategies (reactive
foraging, hazard-avoidance, patch-rotation to exploit regrowth, ...). We evolve survivors, then the
agent-observer NAMES whatever the best one is competently doing and classifies it against the
4-form map {navigation, memory, delayed-gratification, regulation} -- a strategy that fits none is
a candidate for genuinely new competency.

CLI (for the agent-observer; reads the evolved survivor in survivor.json):
  python3 survival_world.py run [--seed S] [--intervene "TYPE@t"]
  TYPE: movefood (teleport all food), addhazard (ring of hazards round the agent), sap (halve
  energy), wallfood (wall off the nearest food). Prints sampled ASCII frames (A=agent F=food
  H=hazard .=empty) + energy + alive.
"""
import sys, json, random

N = 11
E0, EAT, HAZ, STEP, REGROW, T = 15, 7, 5, 1, 12, 100
DIRS = [(0, 0), (1, 0), (-1, 0), (0, 1), (0, -1)]
NSTATES = 6


def layout(r, nf=8, nh=6):
    cells = [(x, y) for x in range(N) for y in range(N)]
    r.shuffle(cells)
    food = set(cells[:nf]); haz = set(cells[nf:nf + nh])
    start = cells[nf + nh]
    return food, haz, start


def sense(p, food, haz):
    near = min(food, key=lambda f: abs(f[0] - p[0]) + abs(f[1] - p[1])) if food else p
    dx = (near[0] > p[0]) - (near[0] < p[0]); dy = (near[1] > p[1]) - (near[1] < p[1])
    fdir = (dx + 1) * 3 + (dy + 1)
    hb = 0
    for k, (ax, ay) in enumerate(DIRS[1:]):
        if (p[0] + ax, p[1] + ay) in haz:
            hb |= 1 << k
    return fdir * 16 + hb


NIN = 9 * 16


def random_rule(r, nstates=NSTATES):
    return [(r.randrange(5), r.randrange(nstates)) for _ in range(NIN * nstates)]


def mutate(rule, r, rate=0.04):
    c = rule[:]
    for i in range(len(c)):
        if r.random() < rate:
            c[i] = (r.randrange(5), r.randrange(NSTATES))
    return c


def run(rule, food0, haz, start, T=T, trace=False, interv=None, r=None):
    food = set(food0); regrow = {}; p = start; e = E0; s = 0; frames = []; eaten = 0
    for t in range(T):
        if interv and t == interv[0]:
            k = interv[1]
            if k == "movefood" and r:
                food = set(layout(r)[0])
            elif k == "addhazard":
                haz = set(haz) | {(p[0] + a, p[1] + b) for a in (-1, 0, 1) for b in (-1, 0, 1) if (a, b) != (0, 0)
                                  for px, py in [(p[0] + a, p[1] + b)] if 0 <= px < N and 0 <= py < N}
            elif k == "sap":
                e = max(1, e // 2)
            elif k == "wallfood" and food:
                nf = min(food, key=lambda f: abs(f[0] - p[0]) + abs(f[1] - p[1]))
                haz = set(haz) | {(nf[0] + a, nf[1] + b) for a in (-1, 0, 1) for b in (-1, 0, 1)}
        if trace and (t % 12 == 0 or t == T - 1):
            frames.append((t, e, render(p, food, haz)))
        a, s = rule[sense(p, food, haz) * NSTATES + s]; s %= NSTATES
        ax, ay = DIRS[a]; np_ = (p[0] + ax, p[1] + ay)
        if 0 <= np_[0] < N and 0 <= np_[1] < N:
            p = np_
        e -= STEP
        if p in food:
            food.discard(p); regrow[p] = REGROW; e += EAT; eaten += 1
        if p in haz:
            e -= HAZ
        for cell in list(regrow):
            regrow[cell] -= 1
            if regrow[cell] <= 0:
                food.add(cell); del regrow[cell]
        if e <= 0:
            return t, eaten, frames
    return T, eaten, frames


def render(p, food, haz):
    rows = []
    for y in range(N):
        row = ""
        for x in range(N):
            row += "A" if (x, y) == p else "F" if (x, y) in food else "H" if (x, y) in haz else "."
        rows.append(row)
    return "\n".join(rows)


def fitness(rule, r, k=4):
    tot = 0
    for _ in range(k):
        food, haz, start = layout(r)
        surv, eaten, _ = run(rule, food, haz, start)
        tot += surv + 0.5 * eaten
    return tot / k


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "run":
        rule = json.load(open(__file__.rsplit("/", 1)[0] + "/survivor.json"))
        rule = [tuple(x) for x in rule]
        a = sys.argv
        seed = int(a[a.index("--seed") + 1]) if "--seed" in a else 0
        interv = None
        if "--intervene" in a:
            spec = a[a.index("--intervene") + 1]; k, t = spec.split("@"); interv = (int(t), k)
        r = random.Random(seed)
        food, haz, start = layout(r)
        surv, eaten, frames = run(rule, food, haz, start, trace=True, interv=interv, r=random.Random(seed + 1))
        print(f"survived {surv}/{T} steps, ate {eaten} food" + (f", intervention {interv}" if interv else ""))
        for t, e, fr in frames:
            print(f"\n--- t={t} energy={e} ---\n{fr}")
        sys.exit(0)
    # evolve survivors (domain-randomized layouts each gen)
    POP, MU, GENS = 120, 18, 60
    r = random.Random(123)
    pop = [random_rule(random.Random(i)) for i in range(POP)]
    parents = pop[:MU]
    for g in range(GENS):
        rg = random.Random(4000 + g)
        scored = sorted(((fitness(ru, rg), ru) for ru in pop), key=lambda x: -x[0])
        if g % 15 == 0 or g == GENS - 1:
            print(f"gen {g:>2}: best fitness {scored[0][0]:.1f}")
        parents = [ru for _, ru in scored[:MU]]
        pop = parents + [mutate(r.choice(parents), r) for _ in range(POP - MU - 5)] + [random_rule(r) for _ in range(5)]
    best = parents[0]
    json.dump([list(x) for x in best], open(__file__.rsplit("/", 1)[0] + "/survivor.json", "w"))
    # benchmark best vs a random rule on fresh layouts
    rr = random.Random(999)
    bf = sum(run(best, *layout(rr))[0] for _ in range(40)) / 40
    rf = sum(run(random_rule(rr), *layout(rr))[0] for _ in range(40)) / 40
    print(f"\nbest survivor: mean survival {bf:.1f}/{T} vs random rule {rf:.1f}/{T}  (saved survivor.json)")
