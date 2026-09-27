"""Batch C part 2 — #12 contest-resolution (two anonymous agents split: one takes the high-value goal, the
other a low one, avoiding collision) and #15 stigmergy (two agents forage via a shared trail field, no direct
comms/vision). Both spatial -> render_ring3. Run on VPS epc-venv from repo root."""
import sys, numpy as np, torch
ROOT = "/home/matthewhmaxwell/emergent-pattern-catalog"
sys.path.insert(0, ROOT); sys.path.insert(0, ROOT + "/analysis/ring3_competency")
from gallery.mk_ring3_demo import render_ring3, ASSETS
import contest_ppo as CO
import stigmergy_ppo as ST


def contest():
    net = CO.Agent(); net.load_state_dict(torch.load(ROOT + "/analysis/ring3_competency/contest_net.pt", map_location="cpu")); net.eval()
    chosen = None
    for seed in range(200):
        env = CO.VecContest(1, seed); goals = env.goal[0].astype(float).copy(); high = int(env.high[0])
        gc = (float(goals[high][0]), float(goals[high][1]), 0.55)

        def cap():
            g0, g1 = int(env.on_goal(0)[0]), int(env.on_goal(1)[0])
            return {"ap": env.p[0].astype(float).copy(), "acol": np.arange(2), "tgt": goals.copy(),
                    "tgt_done": [bool(g == g0 or g == g1) for g in range(CO.G)], "goal_circle": gc}

        frames = [cap()]
        for t in range(CO.MAXT):
            mv = np.zeros((1, 2), int); sy = np.zeros((1, 2), int)
            for i in range(2):
                out = net(torch.from_numpy(env.obs(i))); ml = out[0]; sl = out[1] if len(out) > 2 else None
                mv[:, i] = ml.argmax(1).numpy()
                if sl is not None: sy[:, i] = sl.argmax(1).numpy()
            env.step(mv); env.last_sym = sy; frames.append(cap())
        g0, g1 = int(env.on_goal(0)[0]), int(env.on_goal(1)[0])
        if (g0 == high) ^ (g1 == high) and g0 >= 0 and g1 >= 0 and g0 != g1:      # one high, one low, no collision
            chosen = (seed, frames); break
    seed, frames = chosen
    leg = [("agent A", "o", "#2b6cb0", "white"), ("agent B", "o", "#dd6b20", "white"),
           ("goals", "s", "none", "#2f855a"), ("high-value (circled)", "o", "none", "#2f855a")]
    print("contest seed", seed, "frames", len(frames),
          render_ring3(frames, ASSETS + "/ring3_contest", "Contest (#12): agents split; one takes the high goal",
                       leg, (-.5, CO.N - .5, -.5, CO.N - .5), coord="xy"))


def stig():
    net = ST.Agent(); net.load_state_dict(torch.load(ROOT + "/analysis/ring3_competency/stigmergy_net.pt", map_location="cpu")); net.eval()
    chosen = None
    for seed in range(120):
        env = ST.VecStig(1, seed)

        def cap():
            food = np.argwhere(env.food[0]).astype(float)
            return {"ap": env.p[0].astype(float).copy(), "acol": np.arange(2),
                    "field": env.visited[0].astype(float).copy(),
                    "food_red": (food.copy() if len(food) else np.zeros((0, 2)))}

        frames = [cap()]
        for t in range(ST.T):
            acts = np.zeros((1, 2), int)
            for i in range(2):
                out = net(torch.from_numpy(env.obs(i))); acts[:, i] = (out[0] if isinstance(out, tuple) else out).argmax(1).numpy()
            env.step(acts); frames.append(cap())
        if int(env.collected[0]) >= 7:
            chosen = (seed, frames); break
    if chosen is None:
        chosen = (0, frames)
    seed, frames = chosen
    leg = [("agents", "o", "#2b6cb0", "white"), ("trail (visited)", "s", "#dd6b20", "#dd6b20"),
           ("food", "D", "none", "#e53e3e")]
    print("stig seed", seed, "collected", int(env.collected[0]), "frames", len(frames),
          render_ring3(frames, ASSETS + "/ring3_stig", "Stigmergy (#15): agents divide ground via shared trails",
                       leg, (-.5, ST.N - .5, -.5, ST.N - .5), coord="xy"))


if __name__ == "__main__":
    contest(); stig()
