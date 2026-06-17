"""EPC gallery animations (v4): balanced 3-phase frame selection so every effect
is SEEN TO TAKE SHAPE — hold on the start, dwell across the transition, hold on
the end. Each producer -> PIL frames; save_animation emits playbar sprite + GIF."""
from __future__ import annotations
import io
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image

FW = FH = 300; DPI = 75
NF = 35           # current frame target (set per-pattern by render)
NF_MIN, NF_MAX = 35, 105
MEASURE = None    # per-pattern per-frame structure measure (set by render via FRAME_MEASURE); overrides the producer default in _select


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


def _select(fr, sfn, n=None):
    """Hold START briefly, dwell across the TRANSITION window, hold END — but
    return DISTINCT frames (no duplicate-padding) and backfill across the run so a
    short transition still yields a smooth, non-repetitive animation. A per-pattern
    MEASURE (registered in FRAME_MEASURE) overrides the producer's default sfn so the
    dwell window tracks the actual canonical observable."""
    sfn = MEASURE or sfn
    if n is None: n = NF
    N = len(fr)
    if N <= 6:
        return fr
    try:
        S = np.array([sfn(f) for f in fr], float); good = np.isfinite(S)
        if good.sum() < 3 or np.nanstd(S[good]) < 1e-9:
            return [fr[i] for i in np.linspace(0, N-1, n).astype(int)]
        S = np.interp(np.arange(N), np.where(good)[0], S[good])
        smin, smax = S.min(), S.max(); rng = smax - smin
        t_lo = int(np.argmax(S >= smin + 0.1*rng))
        after = np.where(S[t_lo:] >= smin + 0.9*rng)[0]
        t_settle = t_lo + int(after[0]) if after.size else N-1
        tail = S[t_settle:]
        if tail.size and float(tail.min()) < smin + 0.5*rng:   # transient (wave): extend to last active
            la = np.where(S >= smin + 0.1*rng)[0]; t_hi = int(la[-1])
        else:                                                   # rise-and-stay: dwell on the rise only
            t_hi = t_settle
        if t_hi - t_lo < 2: t_lo, t_hi = 0, N-1
        # Candidate pool: dense across the dwell window + spread over the whole run.
        pool = sorted(set(
            [int(x) for x in np.linspace(0, t_lo, 4)]
            + [int(x) for x in np.linspace(t_lo, t_hi, n)]
            + [int(x) for x in np.linspace(t_hi, N-1, 8)]
            + [int(x) for x in np.linspace(0, N-1, n)]))
        # Keep frames spaced by EQUAL CHANGE in the structure measure (equal-arc-length):
        # this drops frozen/settled states that render identically and spends the budget
        # where the system is actually changing — no duplicate padding.
        eps = 0.012 * rng
        kept = [pool[0]]
        for i in pool[1:]:
            if abs(S[i] - S[kept[-1]]) >= eps:
                kept.append(i)
        if kept[-1] != pool[-1]:
            kept.append(pool[-1])
        if len(kept) > n:                                      # cap to budget, keep equal-S spread
            kept = [kept[k] for k in np.linspace(0, len(kept)-1, n).astype(int)]
        H = 2                                                  # light orientation holds (was 5)
        idx = [0]*H + kept + [N-1]*H
        return [fr[int(i)] for i in idx]
    except Exception:
        return [fr[i] for i in np.linspace(0, N-1, n).astype(int)]

def _evenly(fr, n=None):
    if n is None: n = NF
    return fr if len(fr) <= n else [fr[i] for i in np.linspace(0, len(fr)-1, n).astype(int)]

def _has(history, key):
    return [f for f in history if isinstance(f, dict) and key in f]


# ---- snapshot producers (3-phase) -------------------------------------------
def point_cloud(history, title, sample="auto", **_):
    base=_has(history, "positions")
    fr = _evenly(base) if sample=="even" else _select(base, _s_points)
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

def grid_field(history, title, sample="auto", **_):
    key = "field" if (history and "field" in history[-1]) else "grid"
    base = _has(history, key)
    fr = _evenly(base) if sample=="even" else _select(base, _s_grid)
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
    fr = fr if len(fr) <= NF else [fr[i] for i in np.linspace(0, len(fr)-1, NF).astype(int)]
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
    for k in np.linspace(0, len(ys)-1, NF).astype(int):
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
    for k in np.linspace(0, Lmax-1, NF).astype(int):
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
    for k in np.linspace(0, len(xs)-1, NF).astype(int):
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
    for fr in np.linspace(0.05, 1.0, NF):
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
    for c in np.linspace(2, img.shape[0], NF).astype(int):
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


def _motion(history, viz, args):
    """Mean normalised frame-to-frame change over the recorded run (0..~0.6)."""
    try:
        a = args or {}
        if viz == "grid_field":
            key = "field" if (history and "field" in history[-1]) else "grid"; fr = _has(history, key)
        elif viz in ("point_cloud", "road_1d"):
            key = "positions"; fr = _has(history, key)
        elif viz == "phase_circle":
            key = "theta" if (history and "theta" in history[-1]) else "phases"; fr = _has(history, key)
        elif viz == "vector_grid":
            key = a.get("key", "state"); fr = _has(history, key)
        elif viz == "network_time":
            key = "edge_weights"; fr = _has(history, key)
        elif viz == "bars_sort":
            key = a.get("key", "array"); fr = _has(history, key)
        elif viz in ("histogram_time", "lorenz_time"):
            key = a.get("key"); fr = _has(history, key)
        else:
            return None
        if len(fr) < 3:
            return 0.0
        fr = [fr[i] for i in np.linspace(0, len(fr)-1, min(60, len(fr))).astype(int)]
        if viz in ("point_cloud", "road_1d"):
            P = [np.asarray(f["positions"], float) for f in fr]
            P = [p[:, :2] if p.ndim == 2 else p.reshape(-1, 1) for p in P]
            ext = max(np.ptp(np.concatenate([p[:, 0] for p in P])), 1e-9)
            return float(np.mean([np.mean(np.abs(P[i+1][:min(len(P[i]), len(P[i+1]))]
                          - P[i][:min(len(P[i]), len(P[i+1]))]))
                          for i in range(len(P)-1)]) / ext)
        if viz == "phase_circle":
            T = [np.asarray(f[key], float) for f in fr]
            return float(np.mean([np.mean(np.abs(((T[i+1]-T[i]+np.pi) % (2*np.pi))-np.pi))
                          for i in range(len(T)-1)]) / np.pi)
        A = [np.asarray(f[key], float).ravel() for f in fr]
        m = min(len(x) for x in A); A = [x[:m] for x in A]
        rng = max(np.ptp(np.concatenate(A)), 1e-9)
        return float(np.mean([np.mean(np.abs(A[i+1]-A[i])) for i in range(len(A)-1)]) / rng)
    except Exception:
        return None


def _adaptive_nframes(history, viz, args):
    m = _motion(history, viz, args)
    if m is None:
        return NF_MIN
    return int(np.clip(round(NF_MIN + (m / 0.30) * (NF_MAX - NF_MIN)), NF_MIN, NF_MAX))


def save_animation(frames, out_base, fps=10):
    if len(frames) < 2: return None
    n = len(frames); cols = int(np.ceil(np.sqrt(n))); rows = int(np.ceil(n/cols))
    sheet = Image.new("RGB", (cols*FW, rows*FH), "white")
    for i, im in enumerate(frames):
        sheet.paste(im, ((i % cols)*FW, (i//cols)*FH))
    sheet.save(out_base + "_sprite.png", optimize=True)
    try:
        import imageio.v2 as imageio
        imageio.mimwrite(out_base + ".mp4", [np.asarray(im) for im in frames], fps=12,
                         codec="libx264", macro_block_size=4, ffmpeg_params=["-pix_fmt", "yuv420p"])
        mp4 = out_base.split("/")[-1] + ".mp4"
    except Exception:
        mp4 = None
    return {"frames": n, "cols": cols, "rows": rows, "fw": FW, "fh": FH, "mp4": mp4}

# Per-pattern PER-FRAME structure measure for _select's dwell window — used where the
# producer's generic default (_s_points/_s_grid/_s_net) tracks the wrong observable
# (QA defect class DC1). Each callable takes one history frame dict -> float. Empty
# entries fall back to the producer default. Populated per-pattern in Pass 4b.
FRAME_MEASURE = {}


def render(history, viz, out_base, title, args=None, measure=None):
    global NF, MEASURE
    MEASURE = measure
    try:
        NF = _adaptive_nframes(history, viz, args or {})
        frames = PRODUCERS[viz](history, title, **(args or {}))
    finally:
        MEASURE = None
    return save_animation(frames, out_base, fps=max(10, min(20, NF // 5))) if frames else None
