"""Program search with INTERNAL STATE. The memoryless reactive test showed 0/4000 solve the
away-from-goal trap -- but a hand-built reactive navigator also failed it, so that 0 was a
REPRESENTATION WALL, not competency rarity: committing to a goal-defying detour needs memory.

Here the rule is a tiny finite-state controller: (sensory_code, state) -> (action, next_state).
Step 1: a hand-built STATEFUL positive control (Pledge wall-follower) must solve the strict trap
-> proves the test is passable, so any 0 means 'not found', not 'impossible'.
Step 2: re-run the funnel over random FSMs with the HONEST strict test (away-from-goal trap from
hard high-y starts) -> is competency now representable AND findable by random search?
"""
import random
from collections import Counter

N = 13
DIRS = [(1, 0), (0, 1), (-1, 0), (0, -1)]                  # E, N, W, S  (index 0..3)
ACTIONS = [(0, 0)] + DIRS                                  # stay, E, N, W, S
NSTATES = 4
NCODES = 9 * 16
RES = (N - 1, N - 1)
STARTS = [(0, 0), (0, N - 1), (N - 1, 0), (3, 6), (6, 3)]
HARD = [(0, N - 1), (2, N - 3), (4, N - 1)]               # high-y starts: force a real away-from-goal detour
TRAP = {(RES[0] - 2, y) for y in range(2, N)}             # wall x=10, y=2..12; only gap at y=0,1


def _free(p, walls):
    return 0 <= p[0] < N and 0 <= p[1] < N and p not in walls


def _dir_code(p, r):
    dx = (r[0] > p[0]) - (r[0] < p[0]); dy = (r[1] > p[1]) - (r[1] < p[1])
    return (dx + 1) * 3 + (dy + 1)


def _walls_code(p, walls):
    b = 0
    for k, (ax, ay) in enumerate(DIRS):
        if not _free((p[0] + ax, p[1] + ay), walls):
            b |= 1 << k
    return b


# ---------- random finite-state controller ----------
def make_fsm(seed):
    rng = random.Random(seed)
    return [(rng.randrange(5), rng.randrange(NSTATES)) for _ in range(NCODES * NSTATES)]


def run_fsm(rule, start, resource, walls, steps=160):
    p = start; s = 0; trace = [p]
    for _ in range(steps):
        if p == resource:
            return trace, True
        code = (_dir_code(p, resource) * 16 + _walls_code(p, walls)) * NSTATES + s
        a, s = rule[code]
        ax, ay = ACTIONS[a]
        np_ = (p[0] + ax, p[1] + ay)
        if _free(np_, walls):
            p = np_
        trace.append(p)
    return trace, p == resource


# ---------- hand-built STATEFUL positive control: Pledge wall-follower ----------
def pledge(start, resource, walls, steps=400):
    p = start; trace = [p]
    facing = 0; counter = 0; following = False
    def goal_dir(p):
        dx = (resource[0] > p[0]) - (resource[0] < p[0]); dy = (resource[1] > p[1]) - (resource[1] < p[1])
        if abs(resource[0]-p[0]) >= abs(resource[1]-p[1]) and dx: return (dx, 0)
        if dy: return (0, dy)
        return (dx, 0)
    for _ in range(steps):
        if p == resource:
            return trace, True
        if not following:
            gd = goal_dir(p)
            if _free((p[0]+gd[0], p[1]+gd[1]), walls):
                p = (p[0]+gd[0], p[1]+gd[1]); trace.append(p); continue
            following = True; counter = 0
            facing = DIRS.index(gd) if gd in DIRS else 0
        # wall-following: keep wall on the LEFT; try left, straight, right, back
        for turn, dc in ((+1, +1), (0, 0), (-1, -1), (2, 2)):  # left, straight, right, back
            nf = (facing + turn) % 4
            np_ = (p[0] + DIRS[nf][0], p[1] + DIRS[nf][1])
            if _free(np_, walls):
                counter += dc
                facing = nf; p = np_; trace.append(p)
                break
        gd = goal_dir(p)
        if counter == 0 and _free((p[0]+gd[0], p[1]+gd[1]), walls):
            following = False
    return trace, p == resource


# ---------- DFS explorer: guaranteed stateful control (proves solvability) ----------
def dfs_reaches(start, resource, walls):
    stack = [start]; seen = {start}
    while stack:
        c = stack.pop()
        if c == resource:
            return True
        for ax, ay in DIRS:
            nb = (c[0]+ax, c[1]+ay)
            if _free(nb, walls) and nb not in seen:
                seen.add(nb); stack.append(nb)
    return False


def surplus(rule):
    reached = sum(run_fsm(rule, s, RES, set())[1] for s in STARTS)
    moved = sum(len(set(run_fsm(rule, s, RES, set())[0])) > 4 for s in STARTS)
    return reached >= 2 and moved >= 2


def strict_competent(rule):
    return all(run_fsm(rule, s, RES, TRAP, steps=200)[1] for s in HARD)


if __name__ == "__main__":
    print("=== POSITIVE CONTROLS on the strict away-from-goal trap (must reach from hard starts) ===")
    for st in HARD:
        pl = pledge(st, RES, TRAP)[1]
        print(f"  start {st}: Pledge wall-follower reached={pl} | DFS(stateful) reaches={dfs_reaches(st, RES, TRAP)}")
    print("\n=== program search over random FINITE-STATE controllers ===")
    M = 6000
    surv = [s for s in range(M) if surplus(make_fsm(s))]
    comp = [s for s in surv if strict_competent(make_fsm(s))]
    print(f"  surplus survivors: {len(surv)} / {M}")
    print(f"  STRICT-competent (solve away-from-goal trap from ALL hard starts): {len(comp)}  {comp[:10]}")
