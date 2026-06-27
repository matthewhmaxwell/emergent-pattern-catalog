"""Toy rung 2 (cheap off-map check): a PREREQUISITE environment whose demanded competency was not
pre-specified as one of the four map forms. The agent must reach a prerequisite tile A *before* the
goal G counts (stepping on G before A does nothing). Navigation is made easy (the agent senses the
direction to A and to G); the challenge is the INSTRUMENTAL ORDERING: go to A first, then commit to
G and not drift back. A reactive (memoryless) agent oscillates -- once it leaves A it again senses A
elsewhere and is pulled back -- so it cannot reliably finish; a stateful agent remembers "A done"
and commits to G. Evolve a survivor, then the agent-observer NAMES the emergent competency
open-endedly and classifies it vs the map {navigation, memory, delayed-gratification, regulation}
or flags it OFF-MAP (e.g. instrumental/means-ends chaining).

CLI (reads seq_survivor.json):
  python3 sequence_world.py run [--seed S] [--intervene "TYPE@t"]
  TYPE: moveA (relocate the prerequisite), moveG (relocate the goal), teleport (jump the agent).
  Prints frames (A=agent, P=prerequisite(unvisited), p=prereq(visited), G=goal-locked,
  g=goal-open) + whether the prerequisite is done.
"""
import sys, json, random

N = 9
DIRS = [(0, 0), (1, 0), (-1, 0), (0, 1), (0, -1)]
NSTATES = 4
NIN = 81  # dirA(9) x dirG(9)


def d9(p, q):
    return ((q[0] > p[0]) - (q[0] < p[0]) + 1) * 3 + ((q[1] > p[1]) - (q[1] < p[1]) + 1)


def random_rule(r, nstates=NSTATES):
    return [(r.randrange(5), r.randrange(nstates)) for _ in range(NIN * nstates)]


def mutate(rule, r, nstates, rate=0.05):
    c = rule[:]
    for i in range(len(c)):
        if r.random() < rate:
            c[i] = (r.randrange(5), r.randrange(nstates))
    return c


def episode(rule, nstates, start, Apos, Gpos, T=60, trace=False, interv=None, r=None):
    p = start; s = 0; done = 0; frames = []
    for t in range(T):
        if interv and t == interv[0]:
            k = interv[1]
            if k == "moveA" and r:
                Apos = (r.randrange(N), r.randrange(N))
            elif k == "moveG" and r:
                Gpos = (r.randrange(N), r.randrange(N))
            elif k == "teleport" and r:
                p = (r.randrange(N), r.randrange(N))
        if trace and (t % 8 == 0 or t == T - 1):
            frames.append((t, done, render(p, Apos, Gpos, done)))
        inp = d9(p, Apos) * 9 + d9(p, Gpos)
        a, s = rule[inp * nstates + s]; s %= nstates
        ax, ay = DIRS[a]; np_ = (p[0] + ax, p[1] + ay)
        if 0 <= np_[0] < N and 0 <= np_[1] < N:
            p = np_
        if done == 0 and p == Apos:
            done = 1
        if done == 1 and p == Gpos:
            return 2, frames                                   # success: A then G
    return done, frames                                        # 0 (no A) or 1 (A but not G)


def render(p, Apos, Gpos, done):
    rows = []
    for y in range(N):
        row = ""
        for x in range(N):
            row += ("A" if (x, y) == p else
                    ("p" if done else "P") if (x, y) == Apos else
                    ("g" if done else "G") if (x, y) == Gpos else ".")
        rows.append(row)
    return "\n".join(rows)


def cfg(r):
    cells = [(x, y) for x in range(N) for y in range(N)]; r.shuffle(cells)
    return cells[0], cells[1], cells[2]                         # start, A, G


def fit(rule, nstates, r, trials=24):
    return sum(episode(rule, nstates, *cfg(r))[0] for _ in range(trials)) / trials


def evolve(nstates, gens=70, pop=120, mu=18, seed=0):
    r = random.Random(seed)
    P = [random_rule(random.Random(seed * 100 + i), nstates) for i in range(pop)]
    parents = P[:mu]
    for g in range(gens):
        rg = random.Random(6000 + g)
        scored = sorted(((fit(ru, nstates, rg), ru) for ru in P), key=lambda x: -x[0])
        parents = [ru for _, ru in scored[:mu]]
        P = parents + [mutate(r.choice(parents), r, nstates) for _ in range(pop - mu - 4)] + \
            [random_rule(r, nstates) for _ in range(4)]
    return parents[0]


if __name__ == "__main__":
    if len(sys.argv) > 1 and sys.argv[1] == "run":
        rule = [tuple(x) for x in json.load(open(__file__.rsplit("/", 1)[0] + "/seq_survivor.json"))]
        a = sys.argv
        seed = int(a[a.index("--seed") + 1]) if "--seed" in a else 0
        interv = None
        if "--intervene" in a:
            k, t = a[a.index("--intervene") + 1].split("@"); interv = (int(t), k)
        r = random.Random(seed); start, Apos, Gpos = cfg(r)
        prog, frames = episode(rule, NSTATES, start, Apos, Gpos, trace=True, interv=interv, r=random.Random(seed + 1))
        print(f"result: {'SUCCESS (A then G)' if prog == 2 else ('reached A only' if prog == 1 else 'never reached A')}" + (f", intervention {interv}" if interv else ""))
        for t, d, fr in frames:
            print(f"\n--- t={t} prereq_done={d} ---\n{fr}")
        sys.exit(0)
    react = evolve(1, seed=1)
    state = evolve(NSTATES, seed=2)
    rr = random.Random(999)
    rf = sum(episode(react, 1, *cfg(rr))[0] for _ in range(80)) / 80
    sf = sum(episode(state, NSTATES, *cfg(rr))[0] for _ in range(80)) / 80
    json.dump([list(x) for x in state], open(__file__.rsplit("/", 1)[0] + "/seq_survivor.json", "w"))
    print(f"mean progress (0=no A, 1=A only, 2=A then G):  reactive(1) {rf:.2f}   stateful({NSTATES}) {sf:.2f}")
    print("(reactive oscillates A<->G once it leaves A -> caps below 2; stateful remembers A-done -> ~2. saved.)")
