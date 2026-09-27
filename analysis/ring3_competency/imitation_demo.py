"""Faithful demo for competency #19 — imitation / social learning. A scripted demonstrator walks from the centre
toward the correct goal; the trained learner must read that motion and go to the SAME goal. Renders both agents
+ trails, 3 goals, a circle on the correct goal. Run on VPS epc-venv from repo root."""
import sys, numpy as np, torch
ROOT = "/home/matthewhmaxwell/emergent-pattern-catalog"
sys.path.insert(0, ROOT); sys.path.insert(0, ROOT + "/analysis/ring3_competency")
from gallery.mk_ring3_demo import render_ring3, ASSETS
import imitation_ppo as W


def load(ckpt):
    sd = torch.load(ckpt, map_location="cpu")
    net = W.Pol(); net.load_state_dict(sd); net.eval(); return net


def episode(net, seed):
    env = W.VecImit(1, seed, "normal"); obs = env.obs()
    correct = int(env.correct[0]); gpos = env.gpos[0].astype(float).copy()
    gc = (float(gpos[correct][0]), float(gpos[correct][1]), 0.6)

    def cap():
        return {"ap": np.array([env.ap[0], env.dp[0]], float), "acolors": ["#2b6cb0", "#dd6b20"],
                "tgt": gpos.copy(), "tgt_done": [bool((env.ap[0] == env.gpos[0, g]).all()) for g in range(W.G)],
                "goal_circle": gc}

    frames = [cap()]
    for t in range(W.T):
        with torch.no_grad():
            out = net(torch.from_numpy(obs)); lg = out[0] if isinstance(out, tuple) else out
        a = lg.argmax(1).numpy(); env.step(a); obs = env.obs(); frames.append(cap())
        if env.done[0]:
            break
    return frames, bool(env.win[0])


def main():
    net = load(ROOT + "/analysis/ring3_competency/imitation_net.pt")
    pick = None
    for seed in range(300):
        fr, win = episode(net, seed)
        if win and 6 <= len(fr) <= 20:
            pick = (seed, fr); break
    if pick is None:
        for seed in range(300):
            fr, win = episode(net, seed)
            if win:
                pick = (seed, fr); break
    seed, fr = pick
    leg = [("learner", "o", "#2b6cb0", "white"), ("demonstrator", "o", "#dd6b20", "white"),
           ("goals", "s", "none", "#2f855a"), ("correct goal", "o", "none", "#2f855a")]
    spr = render_ring3(fr, ASSETS + "/ring3_imitation", "Imitation (#19): follow the demonstrator's goal",
                       leg, (-.5, W.N - .5, -.5, W.N - .5), coord="xy")
    print("imitation seed", seed, "frames", len(fr), spr)


if __name__ == "__main__":
    main()
