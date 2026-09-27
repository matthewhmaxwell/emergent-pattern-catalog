"""Faithful demos for the meta-learning competencies (share the metaworld env): #9 conditional-selection (cued),
#8 rule-tracking (switch), #7 rule-inference (infer). The trained GRU agent collects the GOOD object type
(cued at t=0 / flips after each good pick / hidden-must-infer-from-feedback). Objects colored by type;
collected picks turn green (good) or red-x (bad); agent + trail. Run on VPS epc-venv from repo root."""
import sys, numpy as np, torch
ROOT = "/home/matthewhmaxwell/emergent-pattern-catalog"
sys.path.insert(0, ROOT); sys.path.insert(0, ROOT + "/analysis/ring3_competency")
from gallery.mk_ring3_demo import ASSETS
import render_extra as RX
import metaworld_ppo as W

TYPECOL = ["#dd6b20", "#805ad5", "#00a3c4"]     # 3 object types (distinct from blue agent / green-red picks)
RULES = {"cued": ("Conditional selection (#9): collect the cued type", "ring3_cued"),
         "switch": ("Rule-tracking (#8): good type flips; agent adapts", "ring3_switch"),
         "infer": ("Rule-inference (#7): infer the good type from feedback", "ring3_infer")}


def load(ckpt):
    sd = torch.load(ckpt, map_location="cpu"); H = int(sd["gru.weight_hh_l0"].shape[1]); W.H = H
    net = W.Policy(); net.load_state_dict(sd); net.eval(); return net


def episode(net, rule, seed):
    env = W.VecMeta(1, rule, seed); obs = env.obs(); h = None
    S = env.S; otype = env.otype[0].copy(); state = ["alive"] * S

    def cap():
        mk = []
        for j in range(S):
            x, y = env.obj[0, j]
            if state[j] == "alive":
                mk.append((float(x), float(y), TYPECOL[int(otype[j])], "D", False, 55))
            elif state[j] == "good":
                mk.append((float(x), float(y), "#38a169", "D", True, 62))
            else:
                mk.append((float(x), float(y), "#e53e3e", "x", True, 62))
        d = {"markers": mk, "agents": [(float(env.pos[0, 0]), float(env.pos[0, 1]), "#2b6cb0")]}
        g = int(env.g[0])                                        # current good type index
        if env.rule == "cued":                                  # #9: target KNOWN from the t=0 cue
            d["hud"] = "cued type"; d["hud_color"] = TYPECOL[g]
        elif env.rule == "switch":                              # #8: target FLIPS after each success
            d["hud"] = "good type (flips on hit)"; d["hud_color"] = TYPECOL[g]
        else:                                                   # #7: target hidden, inferred from feedback
            if int(env.good[0]) >= 1:
                d["hud"] = "good type (inferred)"; d["hud_color"] = TYPECOL[g]
            else:
                d["hud"] = "good type: inferring…"
        return d

    frames = [cap()]
    for t in range(W.MAXT):
        alive_before = env.alive[0].copy()
        with torch.no_grad():
            lg, v, h = net(torch.from_numpy(obs)[:, None, :], h)
        a = lg[:, 0].argmax(1).numpy(); obs, r, term, _ = env.step(a)
        for j in np.where(alive_before & ~env.alive[0])[0]:
            state[j] = "good" if env.last_r[0] > 0 else "bad"
        frames.append(cap())
        if int(env.good[0]) >= 4 or t >= 26:
            break
    return frames, int(env.good[0]), int(env.bad[0])


def main():
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    rules = list(RULES) if which == "all" else [which]
    for rule in rules:
        net = load(ROOT + f"/analysis/ring3_competency/metaworld_{rule}.pt")
        pick = None
        for seed in range(400):
            fr, g, b = episode(net, rule, seed)
            if g >= 4 and b <= (1 if rule != "infer" else 2) and 8 <= len(fr) <= 28:
                pick = (seed, fr, g, b); break
        if pick is None:
            for seed in range(400):
                fr, g, b = episode(net, rule, seed)
                if g >= 3:
                    pick = (seed, fr, g, b); break
        seed, fr, g, b = pick
        title, tag = RULES[rule]
        leg = [("agent", "o", "#2b6cb0", "white"), ("object types", "D", "none", "#805ad5"),
               ("good pick", "D", "#38a169", "#38a169"), ("bad pick", "x", "#e53e3e", "#e53e3e")]
        spr = RX.render_markers(fr, ASSETS + "/" + tag, title, leg, (-.5, W.N - .5, -.5, W.N - .5))
        print(f"{rule} seed {seed} good {g} bad {b} frames {len(fr)}", spr)


if __name__ == "__main__":
    main()
