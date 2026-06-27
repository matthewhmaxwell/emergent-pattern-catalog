"""Phase B, Layer 2: the behavior-tracing + intervention harness the agent-observer drives to NAME a
demanding PPO environment's competency -- blind to the reward rule. Loads a trained PPO net
(openworld_ppo_env<E>.pt) and reports neutral behavior (object-types touched in order, counts per
type, goal-reach step) + interventions, so the observer infers the competency, debunks mirages, and
classifies vs the map or flags OFF-MAP. The hidden rule is used only to score success.

CLI (on VPS epc-venv): python openworld_ppo_probe.py run --env E [--seed S] [--remove T] [--double] [--hidegoal K]
  --remove T   : start with all type-T objects gone
  --double     : 2x objects per type (the decisive counting / threshold test)
  --hidegoal K : goal not enterable until step K
"""
import numpy as np, sys, json, torch
import openworld_ppo as W


def traced(net, lits, seed, remove=None, double=False, hidegoal=0):
    rng = np.random.default_rng(seed)
    per = W.PER * 2 if double else W.PER
    types = [t for t in range(W.T) if t != remove]
    objs = []  # [x,y,type,alive]
    ncells = per * len(types) + 2
    idx = rng.permutation(W.N * W.N)[:ncells]
    k = 0
    for t in types:
        for _ in range(per):
            objs.append([idx[k] % W.N, idx[k] // W.N, t, True]); k += 1
    goal = (idx[k] % W.N, idx[k] // W.N); start = (idx[k + 1] % W.N, idx[k + 1] // W.N)
    objs = np.array([[o[0], o[1]] for o in objs]); otype = np.array([o[2] for o in objs])
    alive = np.ones(len(objs), bool)
    p = list(start); h = None; counts = np.zeros(W.T, int); first = -np.ones(W.T, int)
    touches = []; just = -1; goal_step = None; success = 0
    for step in range(W.MAXT):
        x = np.zeros(W.OBS, np.float32); o = 0
        for t in range(W.T):
            m = (otype == t) & alive
            if m.any():
                d = np.abs(objs - p).sum(1).astype(float); d[~m] = 1e9
                nr = objs[d.argmin()]; x[o], x[o + 1] = np.sign(nr[0] - p[0]), np.sign(nr[1] - p[1])
            o += 2
        x[o], x[o + 1] = np.sign(goal[0] - p[0]), np.sign(goal[1] - p[1]); o += 2
        if just >= 0: x[o + just] = 1.0
        with torch.no_grad():
            logit, _, h = net(torch.from_numpy(x)[None, None, :], h)
        a = int(logit[0, 0].argmax())
        np_ = (p[0] + W.DIRS[a][0], p[1] + W.DIRS[a][1])
        if 0 <= np_[0] < W.N and 0 <= np_[1] < W.N: p = list(np_)
        just = -1
        for j in range(len(objs)):
            if alive[j] and objs[j][0] == p[0] and objs[j][1] == p[1]:
                alive[j] = False; counts[otype[j]] += 1
                if first[otype[j]] < 0: first[otype[j]] = step
                just = int(otype[j]); touches.append((step, just))
        frac = W.pred_frac_vec(lits, counts[None], first[None])[0]
        if tuple(p) == tuple(goal) and step >= hidegoal:
            goal_step = step; success = int(frac >= 0.999); break
    return touches, counts.tolist(), goal_step, success


if __name__ == "__main__":
    a = sys.argv; here = __file__.rsplit("/", 1)[0]
    env = int(a[a.index("--env") + 1]); seed = int(a[a.index("--seed") + 1]) if "--seed" in a else 0
    remove = int(a[a.index("--remove") + 1]) if "--remove" in a else None
    double = "--double" in a; hidegoal = int(a[a.index("--hidegoal") + 1]) if "--hidegoal" in a else 0
    sweep = {x["env"]: x for x in json.load(open(f"{here}/openworld_ppo_sweep.json"))}
    lits = [("collect", tuple(l[1])) if l[0] == "collect" else tuple(l) for l in sweep[env]["lits"]]
    net = W.Policy(); net.load_state_dict(torch.load(f"{here}/openworld_ppo_env{env}.pt")); net.eval()
    tou, counts, gstep, succ = traced(net, lits, seed, remove, double, hidegoal)
    tag = [s for s, c in [(f"removed type-{remove}", remove is not None), ("doubled objects", double), (f"goal locked<{hidegoal}", hidegoal)] if c]
    print(f"env {env} seed {seed}" + (f"  [intervention: {', '.join(tag)}]" if tag else "  [baseline]"))
    print("  object-types touched in order:", [t for _, t in tou] or "none")
    print("  counts per type:", counts)
    print("  reached goal at step:", gstep if gstep is not None else "never")
    print("  (scored success vs hidden rule:", succ, ")")
