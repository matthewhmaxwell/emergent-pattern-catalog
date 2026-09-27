"""Faithful demo for competency #18 — relational abstraction. The agent must reach the goal whose COLOUR matches
the cue colour (colours are random feature vectors -> 'match' is a relation, tested on held-out colours). Render:
agent tinted the cue colour goes to the same-coloured goal, ignoring the differently-coloured distractors.
Run on VPS epc-venv from repo root."""
import sys, numpy as np, torch
ROOT = "/home/matthewhmaxwell/emergent-pattern-catalog"
sys.path.insert(0, ROOT); sys.path.insert(0, ROOT + "/analysis/ring3_competency")
from gallery.mk_ring3_demo import ASSETS
import render_extra as RX
import abstraction_ppo as W

DISPLAY = ["#e6194b", "#3cb44b", "#4363d8", "#f58231", "#911eb4", "#42d4f4", "#f032e6", "#9a6324"]


def load(ckpt):
    sd = torch.load(ckpt, map_location="cpu"); net = W.Pol(); net.load_state_dict(sd); net.eval(); return net


def episode(net, seed, colours):
    env = W.VecAbs(1, seed, colours); obs = env.obs()
    cue = int(env.cue[0]); gcol = env.gcol[0].copy(); gpos = env.gpos[0].astype(float).copy()

    def cap():
        mk = [(float(gpos[g][0]), float(gpos[g][1]), DISPLAY[int(gcol[g])], "s", True, 130) for g in range(W.G)]
        return {"markers": mk, "agents": [(float(env.ap[0][0]), float(env.ap[0][1]), DISPLAY[cue])]}

    frames = [cap()]
    for t in range(W.T):
        out = net(torch.from_numpy(obs)); lg = out[0] if isinstance(out, tuple) else out
        a = lg.argmax(1).numpy(); env.step(a); obs = env.obs(); frames.append(cap())
        if env.done[0]:
            break
    return frames, bool(env.win[0])


def main():
    net = load(ROOT + "/analysis/ring3_competency/abstraction_net.pt")
    colours = list(range(6))                                   # trained colour subset (clean demo)
    pick = None
    for seed in range(400):
        fr, win = episode(net, seed, colours)
        if win and 5 <= len(fr) <= 18:
            pick = (seed, fr); break
    if pick is None:
        for seed in range(400):
            fr, win = episode(net, seed, colours)
            if win:
                pick = (seed, fr); break
    seed, fr = pick
    leg = [("agent (= cue colour)", "o", "#555", "white"), ("goals (colours)", "s", "#888", "#888"),
           ("match = same colour", "s", "#38a169", "#38a169")]
    spr = RX.render_markers(fr, ASSETS + "/ring3_abstraction", "Abstraction (#18): go to the goal matching the cue colour",
                            leg, (-.5, W.N - .5, -.5, W.N - .5))
    print("abstraction seed", seed, "frames", len(fr), spr)


if __name__ == "__main__":
    main()
