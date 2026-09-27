"""Batch C part 1 — #10 communication (openhunt4 rung-3: only the speaker sees the target; the listener must
reach it via the signal) and #11 role-division (commhunt: two agents must end on DIFFERENT goals). Both 2-agent
grids -> render_ring3. Run on VPS epc-venv from repo root."""
import sys, numpy as np, torch
ROOT = "/home/matthewhmaxwell/emergent-pattern-catalog"
sys.path.insert(0, ROOT); sys.path.insert(0, ROOT + "/analysis/ring3_competency")
from gallery.mk_ring3_demo import render_ring3, ASSETS
import openhunt4_ppo as OH
import commhunt_ppo as CH


def comms():
    sd = torch.load(ROOT + "/analysis/ring3_competency/openhunt4_net.pt", map_location="cpu")
    OH.H = int(sd["gru.weight_hh_l0"].shape[1]); net = OH.Pol(); net.load_state_dict(sd); net.eval()
    chosen = None
    for seed in range(200):
        env = OH.VecRef(1, seed); obs = env.obs(); h = [None, None]
        sp = int(env.speaker[0]); tgt = env.tgt[0].astype(float).copy()
        cols = ["#dd6b20" if i == sp else "#2b6cb0" for i in range(2)]
        cap = lambda: {"ap": env.ap[0].astype(float).copy(), "acolors": cols, "tgt": tgt.reshape(1, 2).copy(),
                       "tgt_done": [bool(np.all(env.ap[0, 0] == env.tgt[0]) and np.all(env.ap[0, 1] == env.tgt[0]))]}
        frames = [cap()]
        for t in range(OH.T):
            mv = np.zeros((1, 2), int); sg = np.zeros((1, 2), int)
            for i in range(2):
                ml, sl, v, h[i] = net(torch.from_numpy(obs[:, i])[:, None, :], h[i])
                mv[:, i] = ml[:, 0].argmax(1).numpy(); sg[:, i] = sl[:, 0].argmax(1).numpy()
            env.step(mv, sg)
            if env.hit[0] >= 1:
                frames.append({"ap": np.array([tgt, tgt]), "acolors": cols, "tgt": tgt.reshape(1, 2).copy(),
                               "tgt_done": [True]})     # clean "both met at target" frame (env respawns on contact)
                break
            obs = env.obs(); frames.append(cap())
        if env.hit[0] >= 1 and 5 <= len(frames) <= 20:
            chosen = (seed, frames); break
    seed, frames = chosen
    leg = [("speaker (sees target)", "o", "#dd6b20", "white"), ("listener", "o", "#2b6cb0", "white"),
           ("target", "s", "none", "#2f855a"), ("met", "s", "#38a169", "#2f855a")]
    print("comms seed", seed, "frames", len(frames),
          render_ring3(frames, ASSETS + "/ring3_comms", "Communication (#10): speaker signals target to listener",
                       leg, (-.5, OH.N - .5, -.5, OH.N - .5), coord="xy"))


def roldiv():
    net = CH.Agent(); net.load_state_dict(torch.load(ROOT + "/analysis/ring3_competency/commhunt_role_div.pt", map_location="cpu")); net.eval()
    chosen = None
    for seed in range(200):
        env = CH.Vec2(1, "role_div", seed); goals = env.goal[0].astype(float).copy()

        def cap():
            g0, g1 = int(env.on_goal(0)[0]), int(env.on_goal(1)[0])
            return {"ap": env.p[0].astype(float).copy(), "acol": np.arange(2), "tgt": goals.reshape(-1, 2).copy(),
                    "tgt_done": [bool(g == g0 or g == g1) for g in range(CH.G)]}

        frames = [cap()]
        for t in range(CH.MAXT):
            mv = np.zeros((1, 2), int); sy = np.zeros((1, 2), int)
            for i in range(2):
                ml, sl, v = net(torch.from_numpy(env.obs(i)))
                mv[:, i] = ml.argmax(1).numpy(); sy[:, i] = sl.argmax(1).numpy()
            env.step(mv); env.last_sym = sy
            frames.append(cap())
        g0, g1 = int(env.on_goal(0)[0]), int(env.on_goal(1)[0])
        if g0 >= 0 and g1 >= 0 and g0 != g1:
            chosen = (seed, frames); break
    seed, frames = chosen
    leg = [("agent A", "o", "#2b6cb0", "white"), ("agent B", "o", "#dd6b20", "white"), ("goals (occupied=green)", "s", "none", "#2f855a")]
    print("roldiv seed", seed, "frames", len(frames),
          render_ring3(frames, ASSETS + "/ring3_roldiv", "Role division (#11): agents settle on different goals",
                       leg, (-.5, CH.N - .5, -.5, CH.N - .5), coord="xy"))


if __name__ == "__main__":
    comms(); roldiv()
