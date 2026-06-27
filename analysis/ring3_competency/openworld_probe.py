"""Open-ended off-map hunt, LAYER 2: the behavior-tracing + intervention harness the agent-observer
drives to NAME the competency of a demanding environment's trained RNN agent -- WITHOUT being told the
reward rule. The observer watches neutral behavior (which object-types the agent touches, in what
order, how many, whether/when it reaches the goal), hypothesizes what the agent is competent at,
designs an intervention that a genuine competency survives but mere mechanism fails, predicts both,
runs it, and tries to debunk. Only competency surviving debunking is accepted; then it is classified
vs the map {navigation, memory, delayed-gratification, regulation, combinations} or flagged OFF-MAP.

The reward predicate is used ONLY to score success (hidden from the trace) so the observer must infer
the competency from behavior, not be handed it.

CLI:
  python3 openworld_probe.py run --env E [--seed S] [--remove T] [--double] [--hidegoal K]
  --remove T : start with all type-T objects gone (does the agent still pursue them / adapt?)
  --double   : 2x objects per type (does it stop at the right count, or over-collect?)
  --hidegoal K : goal only becomes steppable after step K (does it collect first, then go?)
"""
import numpy as np, json, sys
import openworld as W


def traced_episode(th, r, lits, remove=None, double=False, hidegoal=0, recurrent=True):
    Wx, Wh, b, Wo, bo = W.unpack(th)
    objs, goal, start = W.layout(r)
    if double:
        extra = [[op, ot, True] for op, ot, al in objs]      # duplicate set at shifted cells
        for ob in extra:
            ob[0] = ((ob[0][0] + 1) % W.N, ob[0][1])
        objs = objs + extra
    if remove is not None:
        objs = [ob for ob in objs if ob[1] != remove]
    p = list(start); h = np.zeros(W.H); counts = [0] * W.T; first = [None] * W.T; last = -1
    touches = []; success = 0; goal_step = None
    for step in range(W.STEPS):
        x = np.zeros(W.OBS); o = 0
        for t in range(W.T):
            best, bd = None, 99
            for op, ot, al in objs:
                if al and ot == t:
                    d = abs(op[0] - p[0]) + abs(op[1] - p[1])
                    if d < bd: bd, best = d, op
            if best is not None:
                sx, sy = W.sgn(p, best); x[o], x[o + 1] = sx, sy
            o += 2
        sx, sy = W.sgn(p, goal); x[o], x[o + 1] = sx, sy; o += 2
        if last >= 0: x[o + last] = 1.0
        h = np.tanh(Wx @ x + (Wh @ h if recurrent else 0.0) + b)
        a = int(np.argmax(Wo @ h + bo))
        np_ = (p[0] + W.DIRS[a][0], p[1] + W.DIRS[a][1])
        if 0 <= np_[0] < W.N and 0 <= np_[1] < W.N: p = list(np_)
        last = -1
        for ob in objs:
            if ob[2] and tuple(p) == ob[0]:
                ob[2] = False; counts[ob[1]] += 1
                if first[ob[1]] is None: first[ob[1]] = step
                last = ob[1]; touches.append((step, ob[1]))
        if tuple(p) == goal and step >= hidegoal:
            goal_step = step
            if lits is not None and W.pred_frac(lits, counts, first) == 1.0:
                success = 1
            break
    return {"touches": touches, "counts": counts, "first": first,
            "goal_step": goal_step, "success": success}


if __name__ == "__main__":
    a = sys.argv
    here = __file__.rsplit("/", 1)[0]
    env = int(a[a.index("--env") + 1])
    seed = int(a[a.index("--seed") + 1]) if "--seed" in a else 0
    remove = int(a[a.index("--remove") + 1]) if "--remove" in a else None
    double = "--double" in a
    hidegoal = int(a[a.index("--hidegoal") + 1]) if "--hidegoal" in a else 0
    prefix = a[a.index("--prefix") + 1] if "--prefix" in a else "openworld"
    sweep = {x["env"]: x for x in json.load(open(f"{here}/{prefix}_sweep.json"))}
    lits = [tuple(l) if not isinstance(l[0], str) else (l[0], *l[1:]) for l in sweep[env]["lits"]]
    # normalize collect tuples (json turns inner tuple into list)
    lits = []
    for l in sweep[env]["lits"]:
        if l[0] == "collect":
            lits.append(("collect", tuple(l[1])))
        else:
            lits.append(tuple(l))
    th = np.load(f"{here}/{prefix}_env{env}.npy")
    res = traced_episode(th, np.random.default_rng(seed), lits, remove=remove, double=double, hidegoal=hidegoal)
    tag = []
    if remove is not None: tag.append(f"removed all type-{remove}")
    if double: tag.append("doubled objects")
    if hidegoal: tag.append(f"goal locked until step {hidegoal}")
    print(f"env {env} seed {seed}" + (f"  [intervention: {', '.join(tag)}]" if tag else "  [baseline]"))
    print("  object-types touched in order:", [t for _, t in res["touches"]] or "none")
    print("  counts per type:", res["counts"])
    print("  reached goal at step:", res["goal_step"] if res["goal_step"] is not None else "never")
    print("  (scored success against the hidden rule:", res["success"], ")")
