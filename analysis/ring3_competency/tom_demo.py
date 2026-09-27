"""Faithful demo for competency #20 — predictive intention-reading (theory-of-mind component).
Runs the TRAINED recurrent predictor on the real ToM env: a mover walks toward its hidden goal for K steps then
freezes; all 3 goals are EQUIDISTANT from the stop, so a still frame is symmetric — only the mover's MOTION picks
g*. Renders gallery-style: mover dot + trail (the velocity cue), 3 goal squares, the observer's current pick
filled green, and a circle on the TRUE goal. Run on VPS epc-venv from repo root:
  PYTHONDONTWRITEBYTECODE=1 /home/matthewhmaxwell/epc-venv/bin/python analysis/ring3_competency/tom_demo.py
"""
import sys, numpy as np, torch
ROOT = "/home/matthewhmaxwell/emergent-pattern-catalog"
sys.path.insert(0, ROOT); sys.path.insert(0, ROOT + "/analysis/ring3_competency")
from gallery.mk_ring3_demo import render_ring3, ASSETS
import tom_ppo as W


def load(ckpt):
    sd = torch.load(ckpt, map_location="cpu")
    H = int(sd["gru.weight_hh_l0"].shape[1])          # infer hidden size from the checkpoint
    W.H = H
    net = W.Pol(recurrent=True); net.load_state_dict(sd); net.eval(); return net


def episode(net, seed):
    env = W.VecToM(1, seed, "normal"); obs = env.obs(); h = None
    gpos = env.gpos[0].astype(float).copy(); gstar = int(env.gstar[0])
    gc = (float(gpos[gstar][1]), float(gpos[gstar][0]), 0.6)     # circle on true goal (x=col,y=row; no _xy on goal_circle)

    def cap(choice):
        return {"ap": env.mpos[0].reshape(1, 2).astype(float).copy(), "acolors": ["#dd6b20"],
                "tgt": gpos.copy(), "tgt_done": [i == choice for i in range(W.G)], "goal_circle": gc}

    frames = [cap(-1)]; choice = -1
    for t in range(W.T):
        with torch.no_grad():
            lg, v, h = net(torch.from_numpy(obs)[:, None, :], h)
        choice = int(lg[:, 0].argmax(1).numpy()[0]); env.step(np.array([choice])); obs = env.obs()
        frames.append(cap(choice))
    return frames, (choice == gstar), gstar


def main():
    net = load(ROOT + "/analysis/ring3_competency/tom_net.pt")
    chosen = None
    for seed in range(9000, 9080):
        frames, correct, gstar = episode(net, seed)
        if correct:                                             # first clean, correct episode
            chosen = (seed, frames); break
    if chosen is None:
        seed, frames, _, _ = 9000, *episode(net, 9000)
    else:
        seed, frames = chosen
    leg = [("mover", "o", "#dd6b20", "white"), ("goals", "s", "none", "#2f855a"),
           ("observer's pick", "s", "#38a169", "#2f855a"), ("true goal", "o", "none", "#2f855a")]
    spr = render_ring3(frames, ASSETS + "/ring3_tom",
                       "Intention-reading (#20): the observer picks the goal from the mover's motion, not its position",
                       leg, (-.5, W.N - .5, -.5, W.N - .5), coord="grid")
    print("tom seed", seed, "frames", len(frames), spr)


if __name__ == "__main__":
    main()
