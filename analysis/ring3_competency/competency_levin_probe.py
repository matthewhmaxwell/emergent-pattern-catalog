"""Levin's competency test, falsifiable form. Goal-directedness = 'same ends through different
means' under a NOVEL barrier (you cannot tell by passive observation; you must intervene).

Cleanest isolation of competency: hold the system, the goal, and the original skill fixed --
change ONLY whether, when a novel obstacle appears, the system RE-DERIVES its means (adaptive)
or REPLAYS a fixed plan it built for the unobstructed world. If the barrier-test is real,
adaptive reaches the goal a new way; replay fails. Run on two unrelated substrates so it isn't
a sorting artifact. A null with no competency (replay) MUST fail, or the test means nothing.
"""
import random

# ================= substrate 1: SORTING (goal = order; barrier = frozen elements) =================
def sortedness(a):
    n = len(a)
    inv = sum(1 for i in range(n) for j in range(i + 1, n) if a[i] > a[j])
    mx = n * (n - 1) // 2
    return 1 - inv / mx if mx else 1.0


def adaptive_sort(a, frozen):
    a = a[:]; plan = []; changed = True
    while changed:
        changed = False
        for i in range(len(a) - 1):
            if frozen[i] or frozen[i + 1]:
                continue
            if a[i] > a[i + 1]:
                a[i], a[i + 1] = a[i + 1], a[i]; plan.append(i); changed = True
    return a, plan


def replay_sort(a, frozen, plan):
    a = a[:]
    for i in plan:                       # blindly execute the plan built for the unobstructed array
        if frozen[i] or frozen[i + 1]:
            continue
        if a[i] > a[i + 1]:
            a[i], a[i + 1] = a[i + 1], a[i]
    return a


def case_sorting(seed=0, n=30, n_frozen=6):
    r = random.Random(seed)
    arr = list(range(n)); r.shuffle(arr)
    frozen = [False] * n
    for i in r.sample(range(n), n_frozen):
        frozen[i] = True
    _, plan0 = adaptive_sort(arr, [False] * n)             # means derived in the UNobstructed world
    base = sortedness(adaptive_sort(arr, [False] * n)[0])
    adaptive = sortedness(adaptive_sort(arr, frozen)[0])   # re-derives means under barrier
    replay = sortedness(replay_sort(arr, frozen, plan0))   # replays fixed means under barrier
    return base, adaptive, replay


# ================= substrate 2: NAVIGATION (goal = reach target; barrier = wall) =================
def bfs(walls, start, target, N):
    from collections import deque
    q = deque([start]); prev = {start: None}
    while q:
        c = q.popleft()
        if c == target:
            break
        x, y = c
        for nx, ny in ((x+1,y),(x-1,y),(x,y+1),(x,y-1)):
            if 0 <= nx < N and 0 <= ny < N and (nx, ny) not in walls and (nx, ny) not in prev:
                prev[(nx, ny)] = c; q.append((nx, ny))
    if target not in prev:
        return None
    path = []; c = target
    while c is not None:
        path.append(c); c = prev[c]
    return path[::-1]


def replay_path(path, walls, target):
    # execute a fixed path; stop the moment it would step into a wall
    pos = path[0]
    for step in path[1:]:
        if step in walls:
            break
        pos = step
    return pos == target


def case_navigation(N=15):
    start, target = (0, 0), (N - 1, N - 1)
    # a U-shaped wall that blocks the direct route and traps a greedy/fixed plan
    walls = set()
    wx = N // 2
    for y in range(0, N - 2):                              # vertical wall with a gap near the top
        walls.add((wx, y))
    path0 = bfs(set(), start, target, N)                  # plan for the UNobstructed world
    rerouted = bfs(walls, start, target, N)               # re-derives means under the wall
    adaptive_reaches = rerouted is not None
    replay_reaches = replay_path(path0, walls, target)    # replays the no-wall plan under the wall
    return adaptive_reaches, replay_reaches


print("=== Levin barrier-test: does 'same-ends-different-means' separate competency from mechanism? ===\n")
print("SUBSTRATE 1 — sorting (barrier = frozen elements):")
bs, ad, rp = [], [], []
for s in range(5):
    b, a, r = case_sorting(s)
    bs.append(b); ad.append(a); rp.append(r)
import statistics as st
print(f"  unobstructed sortedness:           {st.mean(bs):.2f}")
print(f"  adaptive (re-derives means):       {st.mean(ad):.2f}  <- reaches goal despite barrier")
print(f"  replay   (fixed means):            {st.mean(rp):.2f}  <- fixed plan fails under barrier")
print(f"  competency gap (adaptive - replay): {st.mean(ad) - st.mean(rp):+.2f}\n")

print("SUBSTRATE 2 — navigation (barrier = wall across the route):")
a_reach, r_reach = case_navigation()
print(f"  adaptive (re-routes) reaches target: {a_reach}")
print(f"  replay   (fixed path) reaches target: {r_reach}")
print()

sort_ok = (st.mean(ad) - st.mean(rp)) > 0.1
nav_ok = a_reach and not r_reach
print(f"VERDICT: {'barrier-test SEPARATES competency from mechanism on BOTH substrates' if sort_ok and nav_ok else 'inconclusive'}")
print("  (same system, same goal, same original skill — only means-adaptation differs.)")
