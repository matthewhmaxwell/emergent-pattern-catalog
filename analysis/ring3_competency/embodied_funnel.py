"""Program-search funnel, beachhead schema = minimal embodied agents.

General funnel (schema-agnostic in spirit): generate simple ALGORITHMS -> cheap surplus filter
(does it do anything non-trivial?) -> competency probe (here cheap+programmatic: same-ends-
different-means under a novel barrier) -> [later: catalog + platonic]. This is program search
(we vary the RULE), not parameter search.

Schema: a tiny gridworld with one resource. The "program" is a random sensorimotor RULE: a
lookup table from (coarse direction-to-resource x local-wall pattern) -> action. The agent is
NEVER told the resource is a goal; competency = it reliably reaches the resource AND re-routes
to it when the normal path is blocked, despite not being designed to (an unexpected competency
in a minimal model -- the Levin/Zhang case).
"""
import random

N = 13
ACTIONS = [(0, 0), (0, 1), (0, -1), (1, 0), (-1, 0)]      # stay,N,S,E,W
NCODES = 9 * 16                                            # 9 dir-codes x 16 wall-patterns


def _dir_code(p, r):
    dx = (r[0] > p[0]) - (r[0] < p[0]); dy = (r[1] > p[1]) - (r[1] < p[1])
    return (dx + 1) * 3 + (dy + 1)


def _walls_code(p, walls):
    b = 0
    for k, (ax, ay) in enumerate(ACTIONS[1:]):
        nx, ny = p[0] + ax, p[1] + ay
        blocked = not (0 <= nx < N and 0 <= ny < N) or (nx, ny) in walls
        b |= (1 if blocked else 0) << k
    return b


def make_rule(seed):
    rng = random.Random(seed)
    return [rng.randrange(5) for _ in range(NCODES)]


def run(rule, start, resource, walls, steps=90):
    p = start; trace = [p]
    for _ in range(steps):
        if p == resource:
            return trace, True
        code = _dir_code(p, resource) * 16 + _walls_code(p, walls)
        ax, ay = ACTIONS[rule[code]]
        nx, ny = p[0] + ax, p[1] + ay
        if 0 <= nx < N and 0 <= ny < N and (nx, ny) not in walls:
            p = (nx, ny)
        trace.append(p)
    return trace, p == resource


# fixed task geometry; vary only the RULE (program search)
RES = (N - 1, N - 1)
STARTS = [(0, 0), (0, N - 1), (N - 1, 0), (3, 6), (6, 3)]


def surplus(rule):
    # generous: does it reach the resource from at least 2 starts, with non-trivial motion?
    reached = 0; moved = 0
    for s in STARTS:
        tr, ok = run(rule, s, RES, set())
        if ok:
            reached += 1
        if len(set(tr)) > 4:
            moved += 1
    return reached >= 2 and moved >= 2, reached


def competency(rule):
    """For starts it reaches: block a cell on its taken path (novel barrier) and a 'trap' wall;
    competent = still reaches the SAME resource by a DIFFERENT path. Greedy/brittle = fails."""
    rerouted = 0; trials = 0; trapped_ok = 0; trap_trials = 0
    for s in STARTS:
        tr, ok = run(rule, s, RES, set())
        if not ok:
            continue
        mid = tr[len(tr) // 2]
        if mid in (s, RES):
            continue
        trials += 1
        tr2, ok2 = run(rule, s, RES, {mid})           # block a cell on its own path
        if ok2 and tr2 != tr:
            rerouted += 1
        # trap: wall the whole row just shy of the resource except a gap AWAY from the approach
        trap = {(RES[0] - 2, y) for y in range(2, N)}
        trap_trials += 1
        tr3, ok3 = run(rule, s, RES, trap)
        if ok3:
            trapped_ok += 1
    if trials == 0:
        return "no-reach", 0.0
    reroute_rate = rerouted / trials
    trap_rate = trapped_ok / max(1, trap_trials)
    if trap_rate >= 0.5:
        return "COMPETENT", trap_rate                  # solves the away-from-goal trap
    if reroute_rate >= 0.5:
        return "re-routes(weak)", reroute_rate         # adapts to simple blocks, trap unknown
    return "brittle/greedy", reroute_rate


if __name__ == "__main__":
    M = 4000
    surv = []
    for seed in range(M):
        rule = make_rule(seed)
        ok, reached = surplus(rule)
        if ok:
            surv.append((seed, reached))
    print(f"program search over {M} random simple rules (embodied schema)")
    print(f"  surplus survivors (reach resource from >=2 starts, non-trivial): {len(surv)}")
    from collections import Counter
    verdicts = Counter()
    competents = []
    for seed, reached in surv:
        v, score = competency(make_rule(seed))
        verdicts[v] += 1
        if v == "COMPETENT":
            competents.append((seed, round(score, 2)))
    print(f"  competency verdicts among survivors: {dict(verdicts)}")
    print(f"  COMPETENT rules (solve away-from-goal trap): {len(competents)}  e.g. {competents[:8]}")
