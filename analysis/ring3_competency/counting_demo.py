"""Faithful demo for competency #3 — counting / accumulation. Trained recurrent PPO agent must collect EXACTLY
k=3 objects (supply=5, so it must STOP, not exhaust) then reach the goal; the tally lives in the GRU hidden
state. Renders gallery-style: agent + trail, uncollected objects (they vanish as grabbed), goal square.
Run on VPS epc-venv from repo root."""
import sys, numpy as np, torch
ROOT = "/home/matthewhmaxwell/emergent-pattern-catalog"
sys.path.insert(0, ROOT); sys.path.insert(0, ROOT + "/analysis/ring3_competency")
from gallery.mk_ring3_demo import render_ring3, ASSETS
import counting_ppo as W


def load(ckpt):
    sd = torch.load(ckpt, map_location="cpu"); H = int(sd["gru.weight_hh_l0"].shape[1])
    net = W.Policy(H); net.load_state_dict(sd); net.eval(); return net


def episode(net, seed, k=3, supply=5):
    env = W.VecCount(1, k, supply, seed); obs = env.reset(); h = None

    def cap():
        alive = env.obj[0][env.alive[0]]
        return {"ap": env.pos[0].reshape(1, 2).astype(float).copy(), "acolors": ["#2b6cb0"],
                "food_red": (alive.astype(float).copy() if len(alive) else np.zeros((0, 2))),
                "tgt": env.goal[0].reshape(1, 2).astype(float).copy(),
                "tgt_done": [bool((env.pos[0] == env.goal[0]).all())],
                "hud": "collected %d / %d" % (int(env.count[0]), k)}

    frames = [cap()]
    for t in range(W.MAXT):
        with torch.no_grad():
            lg, v, h = net(torch.from_numpy(obs)[:, None, :], h)
        a = lg[:, 0].argmax(1).numpy(); obs, r, term, _ = env.step(a); frames.append(cap())
        if env.done[0]:
            break
    return frames, int(env.count[0]), bool(env.done[0])


def main():
    net = load(ROOT + "/analysis/ring3_competency/counting_ppo_net.pt")
    pick = None
    for seed in range(300):
        fr, cnt, done = episode(net, seed)
        if done and cnt == 3 and 8 <= len(fr) <= 34:
            pick = (seed, fr); break
    if pick is None:
        for seed in range(300):
            fr, cnt, done = episode(net, seed)
            if done and cnt == 3:
                pick = (seed, fr); break
    seed, fr = pick
    leg = [("agent", "o", "#2b6cb0", "white"), ("objects (5 available)", "D", "none", "#e53e3e"),
           ("goal", "s", "none", "#2f855a"), ("reached", "s", "#38a169", "#2f855a")]
    spr = render_ring3(fr, ASSETS + "/ring3_counting",
                       "Counting (#3): collect exactly 3 of 5, then the goal",
                       leg, (-.5, W.N - .5, -.5, W.N - .5), coord="xy")
    print("counting seed", seed, "frames", len(fr), spr)


if __name__ == "__main__":
    main()
