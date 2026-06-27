"""HIDDEN ground-truth systems for the agent-observer test. The scientist must NOT read this.
Ground truth:
  A = genuinely competent re-router (reaches target via different paths around any obstacle).
  B = mechanism: replays a fixed plan built for the empty grid; any novel wall on it -> stuck.
  C = TRAP/illusion: greedy downhill toward target. Converges from ANY start with no walls
      (looks strongly goal-directed) but only makes distance-reducing moves, so a concave wall
      that requires temporarily moving AWAY from the target traps it. No re-routing competency.
"""
from collections import deque

N = 14
TARGET = (N - 1, N - 1)


def _free(x, y, walls):
    return 0 <= x < N and 0 <= y < N and (x, y) not in walls


def sys_A(walls, start):
    walls = set(map(tuple, walls))
    q = deque([start]); prev = {start: None}
    while q:
        c = q.popleft()
        if c == TARGET:
            break
        x, y = c
        for nb in ((x+1, y), (x-1, y), (x, y+1), (x, y-1)):
            if _free(*nb, walls) and nb not in prev:
                prev[nb] = c; q.append(nb)
    if TARGET not in prev:
        return {"path": [list(start)], "reached": False}
    path = []; c = TARGET
    while c is not None:
        path.append(list(c)); c = prev[c]
    return {"path": path[::-1], "reached": True}


def sys_B(walls, start):
    walls = set(map(tuple, walls))
    plan = [tuple(p) for p in sys_A([], (0, 0))["path"]]   # fixed plan: empty-grid path from (0,0)
    pos = plan[0]; trace = [list(pos)]
    for step in plan[1:]:
        if step in walls:
            break
        pos = step; trace.append(list(pos))
    return {"path": trace, "reached": pos == TARGET}


def sys_C(walls, start):
    walls = set(map(tuple, walls))
    pos = tuple(start); trace = [list(pos)]; seen = {pos}
    for _ in range(8 * N * N):
        if pos == TARGET:
            break
        x, y = pos; dx = TARGET[0] - x; dy = TARGET[1] - y
        cands = []
        if abs(dx) >= abs(dy):
            if dx: cands.append((x + (1 if dx > 0 else -1), y))
            if dy: cands.append((x, y + (1 if dy > 0 else -1)))
        else:
            if dy: cands.append((x, y + (1 if dy > 0 else -1)))
            if dx: cands.append((x + (1 if dx > 0 else -1), y))
        moved = False
        for nb in cands:                                   # only distance-REDUCING moves
            if _free(*nb, walls) and nb not in seen:
                pos = nb; trace.append(list(pos)); seen.add(pos); moved = True; break
        if not moved:
            break
    return {"path": trace, "reached": pos == TARGET}


SYSTEMS = {"A": sys_A, "B": sys_B, "C": sys_C}
