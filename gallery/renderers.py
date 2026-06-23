"""EPC gallery renderers — one per morphology family. Pure functions of a model
history -> an image asset (GIF for dynamics, PNG for distributions/series), so
the same renderers serve the static gallery now and a live backend later.
"""
from __future__ import annotations
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import animation

_FPS = 8
_MAXF = 40  # cap animation frames


def _subsample(hist, key=None):
    frames = hist
    if key is not None:
        frames = [f for f in hist if isinstance(f, dict) and key in f]
    if len(frames) > _MAXF:
        idx = np.linspace(0, len(frames) - 1, _MAXF).astype(int)
        frames = [frames[i] for i in idx]
    return frames


def render_point_cloud(history, out_path, title=""):
    frames = _subsample(history, "positions")
    if len(frames) < 2:
        return None
    pos0 = np.asarray(frames[0]["positions"], float)
    allp = np.concatenate([np.asarray(f["positions"], float) for f in frames])
    lo, hi = allp.min(0) - 1, allp.max(0) + 1
    fig, ax = plt.subplots(figsize=(4, 4), dpi=80)
    has_types = "types" in frames[0]
    sc = ax.scatter(pos0[:, 0], pos0[:, 1], s=14, c="#2b6cb0")
    ax.set_xlim(lo[0], hi[0]); ax.set_ylim(lo[1], hi[1]); ax.set_aspect("equal")
    ax.set_xticks([]); ax.set_yticks([]); ax.set_title(title, fontsize=9)
    def upd(i):
        p = np.asarray(frames[i]["positions"], float)
        sc.set_offsets(p[:, :2])
        if has_types:
            sc.set_array(np.asarray(frames[i]["types"], float))
        return sc,
    a = animation.FuncAnimation(fig, upd, frames=len(frames), blit=False)
    a.save(out_path, writer=animation.PillowWriter(fps=_FPS)); plt.close(fig)
    return out_path


def render_grid_field(history, out_path, title=""):
    key = "field" if (history and "field" in history[-1]) else "grid"
    frames = _subsample(history, key)
    if len(frames) < 2:
        return None
    f0 = np.asarray(frames[0][key], float)
    discrete = key == "grid"
    fig, ax = plt.subplots(figsize=(4, 4), dpi=80)
    im = ax.imshow(f0, cmap=("tab10" if discrete else "viridis"),
                   interpolation="nearest", animated=True)
    ax.set_xticks([]); ax.set_yticks([]); ax.set_title(title, fontsize=9)
    def upd(i):
        im.set_array(np.asarray(frames[i][key], float)); return im,
    a = animation.FuncAnimation(fig, upd, frames=len(frames), blit=False)
    a.save(out_path, writer=animation.PillowWriter(fps=_FPS)); plt.close(fig)
    return out_path


def render_phase_circle(history, out_path, title=""):
    key = "theta" if (history and "theta" in history[-1]) else "phases"
    frames = _subsample(history, key)
    if len(frames) < 2:
        return None
    fig, ax = plt.subplots(figsize=(4, 4), dpi=80)
    th0 = np.asarray(frames[0][key], float)
    pts = ax.scatter(np.cos(th0), np.sin(th0), s=18, c="#2b6cb0")
    arrow = ax.plot([0, 0], [0, 0], "-", color="#e53e3e", lw=2)[0]
    circ = plt.Circle((0, 0), 1, fill=False, color="#ccc"); ax.add_patch(circ)
    ax.set_xlim(-1.2, 1.2); ax.set_ylim(-1.2, 1.2); ax.set_aspect("equal")
    ax.set_xticks([]); ax.set_yticks([]); ax.set_title(title, fontsize=9)
    def upd(i):
        th = np.asarray(frames[i][key], float)
        pts.set_offsets(np.c_[np.cos(th), np.sin(th)])
        z = np.exp(1j * th).mean()
        arrow.set_data([0, z.real], [0, z.imag])
        return pts, arrow
    a = animation.FuncAnimation(fig, upd, frames=len(frames), blit=False)
    a.save(out_path, writer=animation.PillowWriter(fps=_FPS)); plt.close(fig)
    return out_path


def _series(history, key):
    vals = [f[key] for f in history if isinstance(f, dict) and key in f
            and np.isscalar(f[key])]
    return np.asarray(vals, float)


def render_timeseries(history, out_path, title="", key=None):
    if key is None:  # auto-pick the scalar per-frame field with most variation
        cand = {}
        for f in history:
            if isinstance(f, dict):
                for k, v in f.items():
                    if isinstance(v, (int, float)) and not isinstance(v, bool) \
                            and k not in ("step", "round", "trial"):
                        cand.setdefault(k, []).append(float(v))
        key = max(cand, key=lambda k: np.std(cand[k]) if len(cand[k]) > 2 else 0,
                  default=None)
    y = _series(history, key) if key else np.array([])
    fig, ax = plt.subplots(figsize=(5, 3.2), dpi=90)
    if y.size:
        ax.plot(y, color="#2b6cb0", lw=1.3)
    ax.set_title(f"{title}\n{key}", fontsize=9); ax.set_xlabel("step")
    fig.tight_layout(); fig.savefig(out_path); plt.close(fig)
    return out_path


def render_distribution(history, out_path, title="", key=None, mode="hist"):
    last = history[-1] if history else {}
    arr = np.asarray(last.get(key, []), float).ravel() if key else np.array([])
    fig, ax = plt.subplots(figsize=(5, 3.2), dpi=90)
    arr = arr[np.isfinite(arr)]
    if arr.size:
        if mode == "loglog":
            arr = arr[arr > 0]
            vals, cnts = np.unique(arr.astype(int), return_counts=True)
            ax.loglog(vals, cnts, "o", ms=4, color="#2b6cb0")
            ax.set_xlabel("size"); ax.set_ylabel("count")
        elif mode == "lorenz":
            s = np.sort(arr); c = np.cumsum(s) / s.sum()
            x = np.linspace(0, 1, len(c))
            ax.plot(x, c, color="#2b6cb0"); ax.plot([0, 1], [0, 1], "--", color="#ccc")
            ax.set_xlabel("cumulative share of agents"); ax.set_ylabel("cumulative wealth")
        else:
            ax.hist(arr, bins=24, color="#2b6cb0")
    ax.set_title(f"{title}\n{key} ({mode})", fontsize=9)
    fig.tight_layout(); fig.savefig(out_path); plt.close(fig)
    return out_path


def render_network(history, out_path, title=""):
    last = history[-1] if history else {}
    fig, ax = plt.subplots(figsize=(4.2, 4), dpi=85)
    if "edge_weights" in last:
        W = np.asarray(last["edge_weights"], float)
        pos = np.asarray(last.get("node_positions",
              np.c_[np.cos(np.linspace(0, 2*np.pi, len(W), endpoint=False)),
                    np.sin(np.linspace(0, 2*np.pi, len(W), endpoint=False))]), float)
        wmax = W.max() or 1
        for i in range(len(W)):
            for j in range(i+1, len(W)):
                w = max(W[i, j], W[j, i])
                if w > 0.05*wmax:
                    ax.plot([pos[i,0],pos[j,0]],[pos[i,1],pos[j,1]],"-",
                            color="#2b6cb0", lw=0.4+3*w/wmax, alpha=0.7)
        ax.scatter(pos[:,0], pos[:,1], s=60, c="#e53e3e", zorder=3)
        ax.set_aspect("equal")
    elif "task_assignments" in last:
        T = np.array([np.asarray(f["task_assignments"]) for f in history
                      if "task_assignments" in f])
        im = ax.imshow(T.T, aspect="auto", cmap="tab10", interpolation="nearest")
        ax.set_xlabel("step"); ax.set_ylabel("agent")
    ax.set_xticks([]); ax.set_yticks([]); ax.set_title(title, fontsize=9)
    fig.tight_layout(); fig.savefig(out_path); plt.close(fig)
    return out_path


def render_spacetime(history, out_path, title="", key="array"):
    rows = [np.asarray(f[key], float) for f in history if isinstance(f, dict) and key in f]
    if len(rows) < 2:
        return None
    L = min(len(r) for r in rows)
    img = np.array([r[:L] for r in rows])  # time x position
    fig, ax = plt.subplots(figsize=(5, 3.6), dpi=90)
    ax.imshow(img, aspect="auto", cmap="viridis", interpolation="nearest", origin="lower")
    ax.set_xlabel("position"); ax.set_ylabel("step"); ax.set_title(title, fontsize=9)
    fig.tight_layout(); fig.savefig(out_path); plt.close(fig)
    return out_path


def render_director(history, out_path, title=""):
    """Nematic director field: coarse-grained orientation theta in [0,pi) as a cyclic
    color map. Nematic domains + ±1/2 topological defects show as color textures and
    singularities. (P33 active nematic.)"""
    from scipy.ndimage import uniform_filter
    frames = _subsample(history, "theta_field")
    if len(frames) < 2:
        return None
    def cg(th):
        th = np.asarray(th, float)
        c = uniform_filter(np.cos(2 * th), 3, mode="wrap")
        s = uniform_filter(np.sin(2 * th), 3, mode="wrap")
        return (np.arctan2(s, c) / 2.0) % np.pi
    fig, ax = plt.subplots(figsize=(4, 4), dpi=80)
    im = ax.imshow(cg(frames[0]["theta_field"]), cmap="twilight",
                   interpolation="bilinear", vmin=0, vmax=np.pi, animated=True)
    ax.set_xticks([]); ax.set_yticks([]); ax.set_title(title, fontsize=9)
    def upd(i):
        im.set_array(cg(frames[i]["theta_field"])); return im,
    a = animation.FuncAnimation(fig, upd, frames=len(frames), blit=False)
    a.save(out_path, writer=animation.PillowWriter(fps=_FPS)); plt.close(fig)
    return out_path


def render_adjacency_blocks(history, out_path, title=""):
    """Adaptive network: adjacency reordered by node state. Coevolution fragments the
    graph into same-state communities -> two diagonal blocks emerge. (P34.)"""
    frames = _subsample(history, "adjacency")
    if len(frames) < 2:
        return None
    def reorder(f):
        A = np.asarray(f["adjacency"], float)
        if "opinions" in f:
            order = np.argsort(np.asarray(f["opinions"]))
            A = A[np.ix_(order, order)]
        return A
    fig, ax = plt.subplots(figsize=(4, 4), dpi=80)
    im = ax.imshow(reorder(frames[0]), cmap="binary", interpolation="nearest", animated=True)
    ax.set_xticks([]); ax.set_yticks([]); ax.set_title(title, fontsize=9)
    def upd(i):
        im.set_array(reorder(frames[i])); return im,
    a = animation.FuncAnimation(fig, upd, frames=len(frames), blit=False)
    a.save(out_path, writer=animation.PillowWriter(fps=_FPS)); plt.close(fig)
    return out_path


def render_hexatic(history, out_path, title=""):
    """Entropy crystallization: disks colored by local hexatic order |psi6_i|.
    Crystalline (high psi6) grains brighten as the packing self-orders. (P35.)"""
    frames = _subsample(history, "positions")
    if len(frames) < 2:
        return None
    allp = np.concatenate([np.asarray(f["positions"], float) for f in frames])
    lo, hi = allp.min(0) - 0.5, allp.max(0) + 0.5
    L = float(allp.max() - allp.min())
    def psi6(p):
        p = np.asarray(p, float); rij = p[:, None, :] - p[None, :, :]; rij -= L * np.round(rij / L)
        d = np.sqrt((rij ** 2).sum(-1)); np.fill_diagonal(d, 1e9); out = np.empty(len(p))
        for i in range(len(p)):
            nn = np.argsort(d[i])[:6]
            out[i] = abs(np.mean(np.exp(6j * np.arctan2(rij[i, nn, 1], rij[i, nn, 0]))))
        return out
    fig, ax = plt.subplots(figsize=(4, 4), dpi=80)
    p0 = np.asarray(frames[0]["positions"], float)
    sc = ax.scatter(p0[:, 0], p0[:, 1], s=30, c=psi6(p0), cmap="viridis", vmin=0, vmax=1)
    ax.set_xlim(lo[0], hi[0]); ax.set_ylim(lo[1], hi[1]); ax.set_aspect("equal")
    ax.set_xticks([]); ax.set_yticks([]); ax.set_title(title, fontsize=9)
    def upd(i):
        p = np.asarray(frames[i]["positions"], float)
        sc.set_offsets(p[:, :2]); sc.set_array(psi6(p)); return sc,
    a = animation.FuncAnimation(fig, upd, frames=len(frames), blit=False)
    a.save(out_path, writer=animation.PillowWriter(fps=_FPS)); plt.close(fig)
    return out_path


RENDERERS = {
    "director": render_director, "adjacency_blocks": render_adjacency_blocks,
    "hexatic": render_hexatic,
    "point_cloud": render_point_cloud, "grid_field": render_grid_field,
    "phase_circle": render_phase_circle, "timeseries": render_timeseries,
    "distribution": render_distribution, "network": render_network, "spacetime": render_spacetime,
}
