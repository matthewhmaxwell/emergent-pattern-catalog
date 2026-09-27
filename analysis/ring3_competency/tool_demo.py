"""Faithful demo for competency #16 — tool use / instrumental construction. Trained recurrent agent must PUSH a
movable block into the impassable gap (middle row) to build a bridge, then cross to the goal (unreachable
without modifying the environment). Renders agent + block + trails, the gap wall (filled cells become the
bridge), and the goal. Run on VPS epc-venv from repo root."""
import sys, numpy as np, torch
ROOT = "/home/matthewhmaxwell/emergent-pattern-catalog"
sys.path.insert(0, ROOT); sys.path.insert(0, ROOT + "/analysis/ring3_competency")
from gallery.mk_ring3_demo import render_ring3, ASSETS
import tool_ppo as W


def load(ckpt):
    sd = torch.load(ckpt, map_location="cpu"); H = int(sd["gru.weight_hh_l0"].shape[1]); W.H = H
    net = W.Pol(); net.load_state_dict(sd); net.eval(); return net


def episode(net, seed):
    env = W.VecTool(1, seed); env.reset(); obs = env.obs(); h = None

    def cap():
        ax, ay = env.ap[0]; gx, gy = env.goal[0]
        aps = [[ay, ax]]; cols = ["#2b6cb0"]                       # (row=y, col=x) for coord='grid'
        if not env.used[0]:
            bx, by = env.bp[0]; aps.append([by, bx]); cols.append("#8b5a2b")
        health = (~env.filled[0]).astype(float)                    # 1 = still-blocked gap cell; 0 = bridged (vanishes)
        return {"ap": np.array(aps, float), "acolors": cols, "wall": (health, W.GAP),
                "tgt": np.array([[gy, gx]], float), "tgt_done": [bool(env.done[0])]}

    frames = [cap()]
    for t in range(W.T):
        with torch.no_grad():
            lg, v, h = net(torch.from_numpy(obs)[:, None, :], h)
        a = lg[:, 0].argmax(1).numpy(); env.step(a); obs = env.obs(); frames.append(cap())
        if env.done[0]:
            break
    return frames, bool(env.done[0])


def main():
    net = load(ROOT + "/analysis/ring3_competency/tool_net.pt")
    pick = None
    for seed in range(400):
        fr, done = episode(net, seed)
        if done and 8 <= len(fr) <= 24:
            pick = (seed, fr); break
    if pick is None:
        for seed in range(400):
            fr, done = episode(net, seed)
            if done:
                pick = (seed, fr); break
    seed, fr = pick
    leg = [("agent", "o", "#2b6cb0", "white"), ("pushable block", "o", "#8b5a2b", "white"),
           ("gap (impassable)", "s", "#4a5568", "#2d3748"), ("goal", "s", "none", "#2f855a")]
    spr = render_ring3(fr, ASSETS + "/ring3_tool", "Tool use (#16): push a block to bridge the gap",
                       leg, (-.5, W.N - .5, -.5, W.N - .5), coord="grid")
    print("tool seed", seed, "frames", len(fr), spr)


if __name__ == "__main__":
    main()
