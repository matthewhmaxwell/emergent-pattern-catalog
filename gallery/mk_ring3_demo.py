"""Faithful Ring-3 competency demos: run each TRAINED policy in its real env, capture per-step state, render in
the gallery visual language (white bg, no ticks, small title) + role legend + faint trails + optional overlays
(pheromone field, targets, goal circle, wall-health, red/blue food, blocked sites) -> sprite-sheet + MP4 + GIF.

Usage: python gallery/mk_ring3_demo.py <name|all>   (run from repo root, epc-venv)
"""
import sys, os, numpy as np, torch, torch.nn as nn
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from PIL import Image
ROOT = "/home/matthewhmaxwell/emergent-pattern-catalog"; sys.path.insert(0, ROOT)
from gallery.anim import save_animation


class TFPol(nn.Module):
    """Causal Transformer policy for #30 — inlined so we never import the unguarded trainer module."""
    def __init__(self, odim, d, nhead=4, layers=2, maxT=80, nact=5):
        super().__init__()
        self.inp = nn.Linear(odim, d); self.pos = nn.Parameter(torch.zeros(1, maxT, d))
        layer = nn.TransformerEncoderLayer(d_model=d, nhead=nhead, dim_feedforward=2 * d, batch_first=True, dropout=0.0)
        self.enc = nn.TransformerEncoder(layer, num_layers=layers)
        self.pi = nn.Linear(d, nact); self.v = nn.Linear(d, 1)

    def forward(self, x, h=None):
        B, Tq, _ = x.shape; z = self.inp(x) + self.pos[:, :Tq]
        mask = torch.triu(torch.full((Tq, Tq), float("-inf")), diagonal=1)
        z = self.enc(z, mask=mask); return self.pi(z), self.v(z).squeeze(-1), None

ASSETS = ROOT + "/gallery/assets"; os.makedirs(ASSETS, exist_ok=True)
FW = FH = 300
PAL = ["#2b6cb0", "#dd6b20", "#38a169", "#805ad5", "#d53f8c", "#00a3c4"]  # agent categorical
BODYC = ["#e53e3e", "#3182ce"]  # 0=red, 1=blue (body/food/site colour)


def _pil(fig):
    fig.canvas.draw(); w, h = fig.canvas.get_width_height()
    buf = np.frombuffer(fig.canvas.buffer_rgba(), np.uint8).reshape(h, w, 4)[..., :3]
    plt.close(fig); return Image.fromarray(buf).resize((FW, FH))


def _xy(pt, coord):
    return (pt[1], pt[0]) if coord == "grid" else (pt[0], pt[1])


def render_ring3(frames, out_base, title, legend, extent, coord="grid", field_cmap="YlOrBr"):
    trails = {}
    fmax = max((np.asarray(f["field"]).max() for f in frames if f.get("field") is not None), default=1.0) or 1.0
    pil = []
    for fr in frames:
        fig, ax = plt.subplots(figsize=(3, 3), dpi=100)
        ax.set_xlim(extent[0], extent[1]); ax.set_ylim(extent[2], extent[3])
        ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([]); ax.set_title(title, fontsize=7.5)
        if fr.get("field") is not None:
            ax.imshow(np.asarray(fr["field"], float).T, origin="lower", cmap=field_cmap, vmin=0, vmax=fmax,
                      alpha=0.75, extent=(extent[0], extent[1], extent[2], extent[3]), aspect="equal", zorder=0)
        if fr.get("wall") is not None:  # (health[], wall_row)
            hp, wr = fr["wall"]
            for c, hh in enumerate(np.asarray(hp)):
                if hh > 0:
                    ax.scatter(*_xy((wr, c), coord), s=120, marker="s", c="#4a5568",
                               alpha=0.35 + 0.5 * min(hh, 2) / 2, edgecolors="#2d3748", zorder=1)
        for key, mk in (("food_red", 0), ("food_blue", 1)):
            if fr.get(key) is not None:
                for r, c in np.asarray(fr[key], float).reshape(-1, 2):
                    ax.scatter(*_xy((r, c), coord), s=55, marker="D", facecolors="none",
                               edgecolors=BODYC[mk], linewidths=1.6, zorder=2)
        if fr.get("tgt") is not None:
            done = fr.get("tgt_done")
            for j, (r, c) in enumerate(np.asarray(fr["tgt"], float).reshape(-1, 2)):
                fl = done is not None and done[j]
                ax.scatter(*_xy((r, c), coord), s=70, marker="s",
                           facecolors=("#38a169" if fl else "none"), edgecolors="#2f855a", linewidths=1.3, zorder=2)
        for s in fr.get("sites", []):  # (r,c,coloridx,blocked)
            r, c, ci, bl = s
            ax.scatter(*_xy((r, c), coord), s=240, marker="s", facecolors="none",
                       edgecolors=BODYC[int(ci)], linewidths=2.2, zorder=2)
            if bl:
                ax.scatter(*_xy((r, c), coord), s=210, marker="x", c="#111", linewidths=2.4, zorder=3)
        if fr.get("goal_circle") is not None:
            gx, gy, gr = fr["goal_circle"]
            ax.add_patch(plt.Circle((gx, gy), gr, fill=False, ec="#2f855a", lw=1.8, zorder=2))
        ap = np.asarray(fr["ap"], float).reshape(-1, 2); acol = fr.get("acol", np.zeros(len(ap)))
        cols = fr.get("acolors")  # explicit per-agent hex overrides PAL/acol
        for i, pt in enumerate(ap):
            x, y = _xy(pt, coord); trails.setdefault(i, []).append((x, y))
            xs = [p[0] for p in trails[i]][-16:]; ys = [p[1] for p in trails[i]][-16:]
            col = cols[i] if cols is not None else PAL[int(acol[i]) % len(PAL)]
            ax.plot(xs, ys, "-", color=col, alpha=0.28, lw=1.5, zorder=4)
            ax.scatter([x], [y], s=95, color=col, edgecolors="white", linewidths=1.0, zorder=5)
        if legend:
            hs = [Line2D([0], [0], marker=m, color="none", markerfacecolor=fc, markeredgecolor=ec,
                         markersize=8, label=lab, lw=0) for (lab, m, fc, ec) in legend]
            ax.legend(handles=hs, loc="upper center", bbox_to_anchor=(0.5, -0.01), ncol=min(3, len(legend)),
                      fontsize=6, frameon=False, handletextpad=0.3, columnspacing=0.8)
        fig.tight_layout(pad=0.4); pil.append(_pil(fig))
    spr = save_animation(pil, out_base, fps=8)
    try:
        import imageio.v2 as imageio
        imageio.mimwrite(out_base + ".gif", [np.asarray(p) for p in pil], duration=0.12, loop=0)
    except Exception as e:
        print("gif warn:", e)
    return spr


def rollout_gru(net, env, T, n, capture, sub=1):
    obs = env.obs(); h = [None] * n; frames = [capture(env)]
    for t in range(T):
        mv = np.zeros((env.B, n), int)
        for i in range(n):
            with torch.no_grad():
                al, v, h[i] = net(torch.from_numpy(obs[:, i])[:, None, :], h[i])
            mv[:, i] = al[:, 0].argmax(1).numpy()
        env.step(mv); obs = env.obs()
        if (t + 1) % sub == 0: frames.append(capture(env))
    return frames


def load_pol(W, odim, ckpt):
    """Build W.Pol at the checkpoint's actual GRU hidden size (scale nets use H=192, not the module default)."""
    sd = torch.load(ckpt, map_location="cpu")
    W.H = int(sd["gru.weight_hh_l0"].shape[1])
    net = W.Pol(odim); net.load_state_dict(sd); net.eval(); return net


def _src(name):
    d = ROOT + "/gallery/ring3_src/" + name; sys.path.insert(0, d); return d


# ---------------- producers ----------------
def stigmergy():
    d = _src("24_stigmergic-coordination"); import stigworld as W
    cfg = dict(n_ag=4, grid=11, n_tgt=8, decay=0.85); env = W.StigWorld(1, 0, cfg)
    net = W.Pol(env.obs().shape[-1]); net.load_state_dict(torch.load(d + "/stigmergic_coordination_s0.pt", map_location="cpu")); net.eval()
    Ng = np.asarray(env.field).shape[-1]; tc = np.argwhere(np.asarray(env.tgt)[0] > 0)
    cap = lambda e: {"ap": np.asarray(e.ap)[0].copy(), "acol": np.arange(e.n), "field": np.asarray(e.field)[0].copy(),
                     "tgt": tc.copy(), "tgt_done": np.array([bool(np.asarray(e.visited)[0, int(r), int(c)]) for r, c in tc])}
    fr = rollout_gru(net, env, W.T, env.n, cap)
    leg = [("agents (blind)", "o", "#2b6cb0", "white"), ("pheromone", "s", "#dd6b20", "#dd6b20"),
           ("target", "s", "none", "#2f855a"), ("covered", "s", "#38a169", "#2f855a")]
    print("stigmergy cover", round(float(env.cover_score()[0]), 3),
          render_ring3(fr, ASSETS + "/ring3_stigmergy", "Stigmergic coordination (#24): blind agents divide a torus via pheromone trails", leg, (-.5, Ng-.5, -.5, Ng-.5)))


def morphology():
    d = _src("29_morphology-specialization"); import heteroworld as W
    env = W.HeteroWorld(1, 0, n_ag=4)
    net = load_pol(W, env.obs().shape[-1], d + "/morphology_specialization_s0.pt")
    body = np.asarray(env.body)[0]
    def cap(e):
        fr = {"ap": np.asarray(e.ap)[0].copy(), "acolors": [BODYC[int(b)] for b in body]}
        for a, k in (("red", "food_red"), ("blue", "food_blue")):
            if hasattr(e, a): fr[k] = np.asarray(getattr(e, a))[0].copy()
        return fr
    frames = rollout_gru(net, env, W.T, 4, cap, sub=2)
    leg = [("red-body agent", "o", BODYC[0], "white"), ("blue-body agent", "o", BODYC[1], "white"),
           ("red food", "D", "none", BODYC[0]), ("blue food", "D", "none", BODYC[1])]
    print("morphology score", round(env.score(), 3),
          render_ring3(frames, ASSETS + "/ring3_morphology", "Morphology specialization (#29): each body-type harvests its matching food", leg, (-1, W.N, -1, W.N)))


def momentum():
    d = _src("28_momentum-control"); import momentumworld as W
    env = W.MomentumWorld(1, 0, n_ag=2)
    net = load_pol(W, env.obs().shape[-1], d + "/momentum_control_s0.pt")
    def cap(e):
        ag = list(np.asarray(e.ap)[0]) + [np.asarray(e.pp)[0]]  # 2 agents + puck
        g = np.asarray(e.goal)[0]
        return {"ap": np.array(ag), "acolors": ["#2b6cb0", "#2b6cb0", "#dd6b20"],
                "goal_circle": (float(g[0]), float(g[1]), W.GOAL_R)}
    frames = rollout_gru(net, env, W.T, 2, cap)
    leg = [("thruster agents", "o", "#2b6cb0", "white"), ("puck", "o", "#dd6b20", "white"), ("goal", "o", "none", "#2f855a")]
    print("momentum score", round(env.score(), 3),
          render_ring3(frames, ASSETS + "/ring3_momentum", "Momentum control (#28): agents drive a puck through inertia into the goal", leg, (0, W.BOX, 0, W.BOX), coord="xy"))


def niche():
    d = _src("26_niche-construction"); import terraformworld as W
    env = W.TerraformWorld(1, 0, n_ag=2)
    net = load_pol(W, env.obs().shape[-1], d + "/niche_construction_s0.pt")
    def cap(e):
        fr = {"ap": np.asarray(e.ap)[0].copy(), "acol": np.arange(2), "wall": (np.asarray(e.wall)[0].copy(), W.WALL_ROW)}
        if hasattr(e, "goal"): g = np.asarray(e.goal)[0]; fr["tgt"] = np.array([[g[0], g[1]]])
        return fr
    frames = rollout_gru(net, env, W.T, 2, cap, sub=5)
    leg = [("agents", "o", "#2b6cb0", "white"), ("wall (health)", "s", "#4a5568", "#2d3748"), ("goal", "s", "none", "#2f855a")]
    print("niche score", round(env.score(), 3),
          render_ring3(frames, ASSETS + "/ring3_niche", "Niche construction (#26): agents dig a durable breach to reach alternating goals", leg, (-1, W.N, -1, W.N)))


def division():
    d = _src("22_division-of-labor"); import allocworld as W
    cfg = dict(n_ag=2, has_obs=True, has_sig=False, hetero=False, reward_shared=True, asym_init=False)
    env = W.AllocWorld(1, 1, cfg); odim = env.obs().shape[-1]
    net = load_pol(W, odim, d + "/division_of_labor_s0.pt")
    tg = np.asarray(env.tg)[0].reshape(-1, 2)
    cap = lambda e: {"ap": np.asarray(e.ap)[0].copy(), "acol": np.arange(2), "tgt": tg.copy()}
    obs = env.obs(); h = [None, None]; frames = [cap(env)]
    for t in range(W.T):
        mv = np.zeros((1, 2), int); sg = np.zeros((1, 2), int)
        for i in range(2):
            with torch.no_grad():
                ml, sl, v, h[i] = net(torch.from_numpy(obs[:, i])[:, None, :], h[i])
            mv[:, i] = ml[:, 0].argmax(1).numpy(); sg[:, i] = sl[:, 0].argmax(1).numpy()
        env.step(mv, sg); obs = env.obs(); frames.append(cap(env))
    leg = [("agent A", "o", PAL[0], "white"), ("agent B", "o", PAL[1], "white"), ("target cells", "s", "none", "#2f855a")]
    alloc = float(np.asarray(env.alloc_score()).reshape(-1)[0])
    spr = render_ring3(frames, ASSETS + "/ring3_division", "Division of labor (#22): two agents split to one-per-target", leg, (-1, W.N, -1, W.N))
    print("division alloc", round(alloc, 3), spr)


def compositional():
    d = _src("30_compositional-attention"); import adapthetero as E
    env = E.AdaptHetero(1, 0, n_ag=1); odim = env.obs().shape[-1]
    sd = torch.load(d + "/compositional_attention_s0.pt", map_location="cpu")
    dd = int(sd["inp.weight"].shape[0]); maxT = int(sd["pos"].shape[1])
    net = TFPol(odim, dd, nact=E.NACT, maxT=maxT); net.load_state_dict(sd); net.eval()
    body = int(np.asarray(env.body)[0])
    def cap(e):
        blk = int(np.asarray(e._cur_block())[0]); sc = [int(np.asarray(e._site_color(0))[0]), int(np.asarray(e._site_color(1))[0])]
        sites = [(E.SITES[si][0], E.SITES[si][1], sc[si], si == blk) for si in range(2)]
        return {"ap": np.asarray(e.ap)[0].reshape(1, 2), "acolors": [BODYC[body]], "sites": sites}
    seq = [env.obs()[:, 0, :]]; frames = [cap(env)]
    for t in range(E.T):
        X = torch.from_numpy(np.stack(seq, 1))  # (B, t+1, odim)
        with torch.no_grad():
            lg, v, _ = net(X)
        a = lg[:, -1].argmax(-1).numpy()[:, None]
        env.step(a); seq.append(env.obs()[:, 0, :])
        if (t + 1) % 3 == 0: frames.append(cap(env))
    leg = [("agent (by body)", "o", "#805ad5", "white"), ("site serves red", "s", "none", BODYC[0]),
           ("site serves blue", "s", "none", BODYC[1]), ("blocked", "x", "#111", "#111")]
    print("compositional score", round(env.score(), 3),
          render_ring3(frames, ASSETS + "/ring3_compositional", "Compositional attention (#30): agent tracks block + colour signals at once", leg, (-1, E.G, -1, E.G)))


PRODUCERS = {"stigmergy": stigmergy, "morphology": morphology, "momentum": momentum,
             "niche": niche, "division": division, "compositional": compositional}
if __name__ == "__main__":
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    for name in (PRODUCERS if which == "all" else [which]):
        try:
            PRODUCERS[name]()
        except Exception as ex:
            import traceback; print(f"!! {name} FAILED:", repr(ex)); traceback.print_exc()
    print("done")
