"""EPC gallery animations (v4): balanced 3-phase frame selection so every effect
is SEEN TO TAKE SHAPE — hold on the start, dwell across the transition, hold on
the end. Each producer -> PIL frames; save_animation emits playbar sprite + GIF."""
from __future__ import annotations
import io
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image

FW = FH = 300; MAXF = 36; DPI = 75


def _img(fig):
    buf = io.BytesIO(); fig.savefig(buf, format="png", dpi=DPI, facecolor="white")
    plt.close(fig); buf.seek(0)
    im = Image.open(buf).convert("RGB")
    return im.resize((FW, FH)) if im.size != (FW, FH) else im

def _new(title):
    fig, ax = plt.subplots(figsize=(4, 4), dpi=DPI); ax.set_title(title, fontsize=9); return fig, ax

def _headings(f):
    if "headings" in f: return np.asarray(f["headings"], float).ravel()
    if "velocities" in f:
        v = np.asarray(f["velocities"], float); return np.arctan2(v[:, 1], v[:, 0])
    return None

# ---- structure measures (for transition-aware selection) --------------------
def _moran(a):
    a = a - a.mean(); den = float((a*a).sum())
    if den <= 1e-12: return 0.0
    num = float((a[:-1]*a[1:]).sum() + (a[:, :-1]*a[:, 1:]).sum())
    return (num/(a[:-1].size + a[:, :-1].size))/(den/a.size)

def _s_points(f):
    h = _headings(f)
    if h is not None and h.size: return float(abs(np.exp(1j*h).mean()))
    p = np.asarray(f["positions"], float)[:, :2]; n = len(p)
    if n < 3: return 0.0
    d = np.sqrt(((p[:, None]-p[None])**2).sum(-1)); np.fill_diagonal(d, np.inf)
    ext = max((np.ptp(p[:, 0]) + np.ptp(p[:, 1]))/2, 1e-9)
    return float(ext/(d.min(1).mean()+1e-9))

def _s_grid(f):
    k = "field" if "field" in f else "grid"; return _moran(np.asarray(f[k], float))

def _s_phase(f):
    t = np.asarray(f.get("theta", f.get("phases")), float); return float(abs(np.exp(1j*t).mean()))

def _s_road(f):
    v = np.asarray(f.get("velocities", [1]), float); return float(np.mean(v == 0))

def _s_hist(f, key):
    a = np.asarray(f[key], float).ravel(); h, _ = np.histogram(a, bins=16); p = h/max(h.sum(), 1)
    nz = p[p > 0]; return 1.0 - float(-(nz*np.log(nz)).sum()/np.log(16))

def _s_gini(f, key):
    a = np.sort(np.asarray(f[key], float).ravel()); n = len(a); sa = a.sum()
    return float((2*np.arange(1, n+1)-n-1).dot(a)/(n*sa)) if (sa > 0 and n) else 0.0

def _s_net(f):
    return float(np.asarray(f["edge_weights"], float).sum())

def _s_bars(f):
    if "sortedness" in f: return float(f["sortedness"])
    a = np.asarray(f["array"], float); return float((np.diff(a) >= 0).mean())


def _select(fr, sfn, n=MAXF):
    """Hold START, dwell across the TRANSITION window, hold END."""
    if len(fr) <= 6:
        return fr
    try:
        S = np.array([sfn(f) for f in fr], float); good = np.isfinite(S)
        if good.sum() < 3 or np.nanstd(S[good]) < 1e-9:
            return [fr[i] for i in np.linspace(0, len(fr)-1, n).astype(int)]
        S = np.interp(np.arange(len(S)), np.where(good)[0], S[good])
        smin, smax = S.min(), S.max(); rng = smax - smin
        t_lo = int(np.argmax(S >= smin + 0.1*rng))
        after = np.where(S[t_lo:] >= smin + 0.9*rng)[0]
        t_settle = t_lo + int(after[0]) if after.size else len(fr)-1
        tail = S[t_settle:]
        if tail.size and float(tail.min()) < smin + 0.5*rng:   # transient (wave): extend to last active
            la = np.where(S >= smin + 0.1*rng)[0]; t_hi = int(la[-1])
        else:                                                   # rise-and-stay: dwell on the rise only
            t_hi = t_settle
        if t_hi - t_lo < 2: t_lo, t_hi = 0, len(fr)-1
        N = len(fr); H = 5
        idx = ([0]*H + list(np.linspace(0, t_lo, 3, dtype=int))
               + list(np.linspace(t_lo, t_hi, n-2*H-6, dtype=int))
               + list(np.linspace(t_hi, N-1, 3, dtype=int)) + [N-1]*H)
        return [fr[int(i)] for i in idx]
    except Exception:
        return [fr[i] for i in np.linspace(0, len(fr)-1, n).astype(int)]

def _evenly(fr, n=MAXF):
    return fr if len(fr) <= n else [fr[i] for i in np.linspace(0, len(fr)-1, n).astype(int)]

def _has(history, key):
    return [f for f in history if isinstance(f, dict) and key in f]


# ---- snapshot producers (3-phase) -------------------------------------------
def point_cloud(history, title, **_):
    fr = _select(_has(history, "positions"), _s_points)
    if len(fr) < 2: return []
    allp = np.concatenate([np.asarray(f["positions"], float)[:, :2] for f in fr])
    lo, hi = allp.min(0)-1, allp.max(0)+1
    arrows = _headings(fr[0]) is not None; has_t = "types" in fr[0]; has_b = "bonds" in fr[0]
    out = []
    for f in fr:
        p = np.asarray(f["positions"], float)[:, :2]; fig, ax = _new(title)
        if has_b and not arrows:
            for i, j in np.asarray(f.get("bonds", []), int).reshape(-1, 2):
                if i < len(p) and j < len(p):
                    ax.plot([p[i, 0], p[j, 0]], [p[i, 1], p[j, 1]], "-", color="#e53e3e", lw=0.9, alpha=0.7)
        if arrows:
            h = _headings(f)
            ax.quiver(p[:, 0], p[:, 1], np.cos(h), np.sin(h), (h % (2*np.pi)), cmap="hsv",
                      clim=(0, 2*np.pi), scale=22, scale_units="width", width=0.006, headwidth=3, pivot="mid")
        else:
            ax.scatter(p[:, 0], p[:, 1], s=16, c=(np.asarray(f["types"], float) if has_t else "#2b6cb0"),
                       cmap="tab10" if has_t else None)
        ax.set_xlim(lo[0], hi[0]); ax.set_ylim(lo[1], hi[1]); ax.set_aspect("equal")
        ax.set_xticks([]); ax.set_yticks([]); out.append(_img(fig))
    return out

def grid_field(history, title, **_):
    key = "field" if (history and "field" in history[-1]) else "grid"
    fr = _select(_has(history, key), _s_grid)
    if len(fr) < 2: return []
    allv = np.concatenate([np.asarray(f[key], float).ravel() for f in fr])
    vmin, vmax = float(allv.min()), float(allv.max()); disc = key == "grid"
    out = []
    for f in fr:
        fig, ax = _new(title)
        ax.imshow(np.asarray(f[key], float), cmap=("tab10" if disc else "viridis"),
                  vmin=vmin, vmax=vmax, interpolation="nearest")
        ax.set_xticks([]); ax.set_yticks([]); out.append(_img(fig))
    return out

def phase_circle(history, title, **_):
    key = "theta" if (history and "theta" in history[-1]) else "phases"
    fr = _select(_has(history, key), _s_phase)
    if len(fr) < 2: return []
    out = []
    for f in fr:
        th = np.asarray(f[key], float); fig, ax = _new(title)
        ax.add_patch(plt.Circle((0, 0), 1, fill=False, color="#ccc"))
        ax.scatter(np.cos(th), np.sin(th), s=18, c="#2b6cb0")
        z = np.exp(1j*th).mean(); ax.plot([0, z.real], [0, z.imag], "-", color="#e53e3e", lw=2)
        ax.set_xlim(-1.2, 1.2); ax.set_ylim(-1.2, 1.2); ax.set_aspect("equal")
        ax.set_xticks([]); ax.set_yticks([]); out.append(_img(fig))
    return out

def road_1d(history, title, **_):
    fr = _select(_has(history, "positions"), _s_road)
    if len(fr) < 2: return []
    L = float(np.max([np.asarray(f["positions"]).max() for f in fr]))+1
    vmax = float(np.max([np.asarray(f.get("velocities", [1])).max() for f in fr])) or 1
    out = []
    for f in fr:
        p = np.asarray(f["positions"], float).ravel()
        v = np.asarray(f.get("velocities", np.full_like(p, vmax)), float).ravel()
        fig, ax = _new(title)
        ax.scatter(p, np.zeros_like(p), c=v, cmap="RdYlGn", vmin=0, vmax=vmax, s=20)
        ax.set_xlim(0, L); ax.set_ylim(-1, 1); ax.set_yticks([]); ax.set_xlabel("position (red = stopped)")
        out.append(_img(fig))
    return out

def vector_grid(history, title, key="state", **_):
    hist = history
    if hist and isinstance(hist[0], dict) and "trial" in hist[0]:
        t0 = hist[0]["trial"]; hist = [f for f in hist if f.get("trial") == t0]
    fr = _has(hist, key)
    if len(fr) < 2: return []
    nn = np.asarray(fr[0][key]).size; side = int(round(np.sqrt(nn)))
    fr = fr if len(fr) <= MAXF else [fr[i] for i in np.linspace(0, len(fr)-1, MAXF).astype(int)]
    out = []
    for f in fr:
        a = np.asarray(f[key], float).ravel()[:side*side].reshape(side, side)
        fig, ax = _new(title); ax.imshow(a, cmap="bwr", vmin=-1, vmax=1, interpolation="nearest")
        ax.set_xticks([]); ax.set_yticks([]); out.append(_img(fig))
    return out

def histogram_time(history, title, key="opinions", bins=24, rng=None, **_):
    fr = _select(_has(history, key), lambda f: _s_hist(f, key))
    if len(fr) < 2: return []
    allv = np.concatenate([np.asarray(f[key], float).ravel() for f in fr])
    lo, hi = (rng or (float(allv.min()), float(allv.max()))); ymax = 0
    for f in fr:
        h, _ = np.histogram(np.asarray(f[key], float).ravel(), bins=bins, range=(lo, hi)); ymax = max(ymax, h.max())
    out = []
    for f in fr:
        fig, ax = _new(title)
        ax.hist(np.asarray(f[key], float).ravel(), bins=bins, range=(lo, hi), color="#2b6cb0")
        ax.set_xlim(lo, hi); ax.set_ylim(0, ymax*1.05); ax.set_xlabel(key); out.append(_img(fig))
    return out

def lorenz_time(history, title, key="wealth", **_):
    fr = _select(_has(history, key), lambda f: _s_gini(f, key))
    if len(fr) < 2: return []
    out = []
    for f in fr:
        a = np.sort(np.asarray(f[key], float).ravel()); c = np.cumsum(a)/max(a.sum(), 1e-9)
        x = np.linspace(0, 1, len(c)); fig, ax = _new(title)
        ax.plot([0, 1], [0, 1], "--", color="#ccc"); ax.fill_between(x, c, x, color="#2b6cb0", alpha=.3)
        ax.plot(x, c, color="#2b6cb0", lw=2)
        ax.set_xlim(0, 1); ax.set_ylim(0, 1); ax.set_xlabel("share of agents"); ax.set_ylabel("share of wealth")
        out.append(_img(fig))
    return out

def network_time(history, title, **_):
    fr = _select(_has(history, "edge_weights"), _s_net)
    if len(fr) < 2: return []
    last = fr[-1]; nn = len(np.asarray(last["edge_weights"]))
    pos = np.asarray(last.get("node_positions",
          np.c_[np.cos(np.linspace(0, 2*np.pi, nn, endpoint=False)),
                np.sin(np.linspace(0, 2*np.pi, nn, endpoint=False))]), float)
    wmax = max(np.asarray(last["edge_weights"], float).max(), 1e-9)
    out = []
    for f in fr:
        W = np.asarray(f["edge_weights"], float); fig, ax = _new(title)
        for i in range(nn):
            for j in range(i+1, nn):
                w = max(W[i, j], W[j, i])
                if w > 0.04*wmax:
                    ax.plot([pos[i, 0], pos[j, 0]], [pos[i, 1], pos[j, 1]], "-", color="#2b6cb0",
                            lw=0.3+3.5*w/wmax, alpha=.75)
        ax.scatter(pos[:, 0], pos[:, 1], s=70, c="#e53e3e", zorder=3)
        ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([]); out.append(_img(fig))
    return out

def bars_sort(history, title, key="array", **_):
    fr = _select(_has(history, key), _s_bars)
    if len(fr) < 2: return []
    nn = len(np.asarray(fr[0][key]))
    out = []
    for f in fr:
        a = np.asarray(f[key], float); fig, ax = _new(title)
        ax.bar(range(nn), a, color=plt.cm.viridis(a/max(a.max(), 1)), width=1.0)
        ax.set_xticks([]); ax.set_yticks([]); ax.set_xlabel("position"); out.append(_img(fig))
    return out


# ---- progressive producers (inherently show start->change->end) -------------
def line_grow(history, title, key="x", ref=None, **_):
    ys = np.array([float(f[key]) for f in history if isinstance(f, dict) and key in f and np.isscalar(f[key])])
    if ys.size < 2: return []
    if ys.size > 600: ys = ys[np.linspace(0, ys.size-1, 600).astype(int)]
    refv = None
    if ref:
        rv = [float(f[ref]) for f in history if isinstance(f, dict) and ref in f and np.isscalar(f[ref])]
        refv = float(np.median(rv)) if rv else None
    out = []
    for k in np.linspace(0, len(ys)-1, MAXF).astype(int):
        fig, ax = _new(title)
        if refv is not None: ax.axhline(refv, ls="--", color="#e53e3e", lw=1, label=ref)
        ax.plot(ys[:k+1], color="#2b6cb0", lw=1.4)
        ax.set_xlim(0, len(ys)); ax.set_ylim(ys.min()-.05*abs(ys.min())-1e-6, ys.max()*1.05+1e-6)
        ax.set_xlabel("step"); ax.set_ylabel(key)
        if refv is not None: ax.legend(fontsize=8, loc="upper right")
        out.append(_img(fig))
    return out

def multi_line(history, title, key="distance_to_target", group="trial", **_):
    series = {}
    for f in history:
        if isinstance(f, dict) and key in f and group in f:
            series.setdefault(int(f[group]), []).append(float(f[key]))
    series = {g: v for g, v in list(series.items())[:8] if len(v) > 1}
    if not series: return []
    Lmax = max(len(v) for v in series.values()); ymax = max(max(v) for v in series.values())
    out = []
    for k in np.linspace(0, Lmax-1, MAXF).astype(int):
        fig, ax = _new(title)
        for v in series.values(): ax.plot(v[:k+1], lw=1.2, alpha=.8)
        ax.set_xlim(0, Lmax); ax.set_ylim(0, ymax*1.05); ax.set_xlabel("step"); ax.set_ylabel(key)
        out.append(_img(fig))
    return out

def sweep_trace(history, title, x="density", y="fraction_on", **_):
    fr = _evenly(_has(history, x), 120)
    xs = np.array([float(f[x]) for f in fr]); ys = np.array([float(f[y]) for f in fr])
    if xs.size < 2: return []
    out = []
    for k in np.linspace(0, len(xs)-1, MAXF).astype(int):
        fig, ax = _new(title)
        ax.plot(xs[:k+1], ys[:k+1], "-", color="#2b6cb0", lw=1.5)
        ax.scatter([xs[k]], [ys[k]], c="#e53e3e", s=40, zorder=3)
        ax.set_xlim(xs.min(), xs.max()); ax.set_ylim(-.05, 1.05); ax.set_xlabel(x); ax.set_ylabel(y)
        out.append(_img(fig))
    return out

def dist_accumulate(history, title, key="avalanche_sizes", **_):
    a = np.asarray(history[0][key], float); a = a[a > 0]
    if a.size < 10: return []
    out = []
    for fr in np.linspace(0.05, 1.0, MAXF):
        sub = a[:int(fr*len(a))]; vals, cnts = np.unique(sub.astype(int), return_counts=True)
        fig, ax = _new(title); ax.loglog(vals, cnts, "o", ms=4, color="#2b6cb0")
        ax.set_xlim(1, a.max()*1.1); ax.set_ylim(0.8, None)
        ax.set_xlabel("avalanche size"); ax.set_ylabel("count"); out.append(_img(fig))
    return out

def growing_spacetime(history, title, key="task_assignments", **_):
    rows = [np.asarray(f[key], float) for f in history if isinstance(f, dict) and key in f]
    if len(rows) < 2: return []
    L = min(len(r) for r in rows); img = np.array([r[:L] for r in rows])
    if img.shape[0] > 300: img = img[np.linspace(0, img.shape[0]-1, 300).astype(int)]
    out = []
    for c in np.linspace(2, img.shape[0], MAXF).astype(int):
        fig, ax = _new(title)
        ax.imshow(img[:c].T, aspect="auto", cmap="tab10", interpolation="nearest", origin="lower",
                  extent=[0, img.shape[0], 0, img.shape[1]])
        ax.set_xlim(0, img.shape[0]); ax.set_xlabel("step"); ax.set_ylabel("agent"); out.append(_img(fig))
    return out


PRODUCERS = {
    "point_cloud": point_cloud, "grid_field": grid_field, "phase_circle": phase_circle,
    "road_1d": road_1d, "vector_grid": vector_grid, "histogram_time": histogram_time,
    "lorenz_time": lorenz_time, "line_grow": line_grow, "multi_line": multi_line,
    "sweep_trace": sweep_trace, "dist_accumulate": dist_accumulate,
    "network_time": network_time, "bars_sort": bars_sort, "growing_spacetime": growing_spacetime,
}


def save_animation(frames, out_base, fps=10):
    if len(frames) < 2: return None
    n = len(frames); cols = int(np.ceil(np.sqrt(n))); rows = int(np.ceil(n/cols))
    sheet = Image.new("RGB", (cols*FW, rows*FH), "white")
    for i, im in enumerate(frames):
        sheet.paste(im, ((i % cols)*FW, (i//cols)*FH))
    sheet.save(out_base + "_sprite.png", optimize=True)
    frames[0].save(out_base + ".gif", save_all=True, append_images=frames[1:],
                   duration=int(1000/fps), loop=0, optimize=True)
    return {"frames": n, "cols": cols, "rows": rows, "fw": FW, "fh": FH,
            "gif": out_base.split("/")[-1] + ".gif"}

def render(history, viz, out_base, title, args=None):
    frames = PRODUCERS[viz](history, title, **(args or {}))
    return save_animation(frames, out_base) if frames else None
