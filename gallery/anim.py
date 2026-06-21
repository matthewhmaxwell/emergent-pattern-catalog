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

def grid_field(history, title, sample="auto", legend=None, empty_idx=None, **_):
    key = "field" if (history and "field" in history[-1]) else "grid"
    base = _has(history, key)
    fr = _evenly(base) if sample == "even" else _select(base, _s_grid)
    if len(fr) < 2: return []
    allv = np.concatenate([np.asarray(f[key], float).ravel() for f in fr])
    vmin, vmax = float(allv.min()), float(allv.max()); disc = key == "grid"
    if disc:
        from matplotlib.colors import ListedColormap, BoundaryNorm
        import matplotlib.patches as mpatches
        nlev = int(round(vmax)) + 1
        cols = [plt.cm.tab10(i % 10) for i in range(nlev)]
        if empty_idx is not None and 0 <= empty_idx < nlev:
            cols[empty_idx] = (0.91, 0.91, 0.91, 1.0)   # subordinate light gray for empty/idle (QA DC7)
        lcmap = ListedColormap(cols); norm = BoundaryNorm(np.arange(-0.5, nlev + 0.5, 1), nlev)
    out = []
    for f in fr:
        fig, ax = _new(title)
        if disc:                                        # discrete states -> exact per-value tab10 colors (no continuous mis-sampling)
            ax.imshow(np.asarray(f[key], float), cmap=lcmap, norm=norm, interpolation="nearest")
            if legend:
                handles = [mpatches.Patch(color=cols[i], label=legend[i]) for i in range(min(len(legend), nlev))]
                ax.legend(handles=handles, fontsize=6, loc="upper right", framealpha=.9)
        else:
            ax.imshow(np.asarray(f[key], float), cmap="viridis", vmin=vmin, vmax=vmax, interpolation="nearest")
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
    gwmax = max(max(np.asarray(f["edge_weights"], float).max() for f in fr), 1e-9)  # global, for stable line widths
    out = []
    for f in fr:
        W = np.asarray(f["edge_weights"], float); wmf = max(W.max(), 1e-9)
        fig, ax = _new(title)
        for i in range(nn):
            for j in range(i+1, nn):
                w = max(W[i, j], W[j, i])
                if w > 0.06*wmf:                      # PER-FRAME threshold: shows the early dense MESH and the later PRUNE
                    ax.plot([pos[i, 0], pos[j, 0]], [pos[i, 1], pos[j, 1]], "-", color="#2b6cb0",
                            lw=0.4+3.5*w/gwmax, alpha=.75)
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
    import matplotlib.patches as mpatches
    rows = [np.asarray(f[key], float) for f in history if isinstance(f, dict) and key in f]
    if len(rows) < 2: return []
    L = min(len(r) for r in rows); img = np.array([r[:L] for r in rows])
    if img.shape[0] > 300: img = img[np.linspace(0, img.shape[0]-1, 300).astype(int)]
    from matplotlib.colors import ListedColormap, BoundaryNorm
    T = img.shape[0]; ntask = int(img.max()) + 1; base = plt.cm.tab10
    lcmap = ListedColormap([base(t % 10) for t in range(ntask)])         # discrete: each task id -> its exact tab10 color (matches legend)
    norm = BoundaryNorm(np.arange(-0.5, ntask + 0.5, 1), ntask)
    handles = [mpatches.Patch(color=base(t % 10), label=f"task {t}") for t in range(min(ntask, 6))]
    out = []
    for c in np.linspace(max(6, T // NF + 2), T, NF).astype(int):   # start wide enough; grow the raster into a FIXED frame
        fig, ax = _new(title)
        ax.imshow(img[:c].T, aspect="auto", cmap=lcmap, norm=norm, interpolation="nearest",
                  origin="lower", extent=[0, c, 0, img.shape[1]])        # extent matches accumulated steps -> no first-frame stretch
        ax.set_xlim(0, T); ax.set_ylim(0, img.shape[1]); ax.set_xlabel("step"); ax.set_ylabel("agent")
        ax.legend(handles=handles, fontsize=6, loc="upper right", framealpha=.9)
        out.append(_img(fig))
    return out


def pop_time(history, title, keys=("prey_fraction", "predator_fraction"),
             labels=("prey", "predator"), colors=("#2f9e44", "#e53e3e"), ylabel="population fraction", **_):
    """Two populations over time — shows the anti-correlated Lotka-Volterra cycle (P11)."""
    H = [f for f in history if isinstance(f, dict) and keys[0] in f]
    if len(H) < 2: return []
    ser = {k: np.array([float(f[k]) for f in H]) for k in keys}
    L = len(H)
    if L > 600:
        ix = np.linspace(0, L-1, 600).astype(int); ser = {k: v[ix] for k, v in ser.items()}; L = 600
    ymax = max(v.max() for v in ser.values())*1.12 + 1e-9
    out = []
    for k in np.linspace(0, L-1, NF).astype(int):
        fig, ax = _new(title)
        for key, c, lab in zip(keys, colors, labels):
            ax.plot(ser[key][:k+1], color=c, lw=1.7, label=lab)
        ax.set_xlim(0, L); ax.set_ylim(0, ymax); ax.set_xlabel("step"); ax.set_ylabel(ylabel)
        ax.legend(fontsize=8, loc="upper right"); out.append(_img(fig))
    return out

def noise_response(history, title, xkey="noise_level", out_key="x", sig_key="signal", **_):
    """Performance vs noise — the inverted-U signature of stochastic resonance (P26).
    Aggregates the noise sweep in the run and plots signal-tracking |corr| per noise level."""
    H = [f for f in history if isinstance(f, dict) and xkey in f and out_key in f and sig_key in f]
    if len(H) < 10: return []
    nl = np.array([float(f[xkey]) for f in H]); xo = np.array([float(f[out_key]) for f in H])
    sg = np.array([float(f[sig_key]) for f in H]); levels = np.unique(nl)
    if levels.size < 4: return []
    perf = []
    for lv in levels:
        m = nl == lv
        perf.append(abs(np.corrcoef(xo[m], sg[m])[0, 1]) if (m.sum() > 5 and xo[m].std() > 0 and sg[m].std() > 0) else 0.0)
    perf = np.nan_to_num(np.array(perf)); pk = int(np.argmax(perf))
    out = []
    for k in np.linspace(0, levels.size-1, min(NF, levels.size)).astype(int):
        fig, ax = _new(title)
        ax.plot(levels[:k+1], perf[:k+1], "-o", color="#2b6cb0", ms=4)
        if pk <= k and perf[pk] > 0:
            ax.scatter([levels[pk]], [perf[pk]], c="#e53e3e", s=60, zorder=3, label="optimal noise")
            ax.legend(fontsize=8, loc="upper right")
        ax.set_xlim(levels.min(), levels.max()); ax.set_ylim(0, max(perf.max()*1.15, 0.1))
        ax.set_xlabel("noise level"); ax.set_ylabel("signal tracking |corr|"); out.append(_img(fig))
    return out

def phase_space(history, title, **_):
    """Phase vs spatial position — exposes the chimera (a contiguous coherent arc beside an
    incoherent arc), which a phase circle cannot distinguish from ordinary partial sync (P10)."""
    key = "theta" if (history and isinstance(history[-1], dict) and "theta" in history[-1]) else "phases"
    H = [f for f in history if isinstance(f, dict) and key in f and "positions" in f]
    if len(H) < 2: return []
    fr = _select(H, _s_phase)
    out = []
    for f in fr:
        th = np.asarray(f[key], float) % (2*np.pi); pos = np.asarray(f["positions"], float).ravel()
        o = np.argsort(pos); fig, ax = _new(title)
        ax.scatter(pos[o], th[o], s=14, c=th[o], cmap="hsv", vmin=0, vmax=2*np.pi)
        ax.set_xlabel("oscillator position"); ax.set_ylabel("phase θ")
        ax.set_ylim(-0.2, 2*np.pi+0.2); ax.set_yticks([0, np.pi, 2*np.pi]); ax.set_yticklabels(["0", "π", "2π"])
        out.append(_img(fig))
    return out


def spacetime_slice(history, title, cmap="inferno", **_):
    """Space-time raster of a fixed 1-D slice (middle row) over time — traveling waves show as
    diagonal stripes (propagation), so sustained spiral/target waves are visible (P13)."""
    key = "field" if (history and isinstance(history[-1], dict) and "field" in history[-1]) else "grid"
    fr = _has(history, key)
    if len(fr) < 3: return []
    grids = [np.asarray(f[key], float) for f in fr]
    r = grids[0].shape[0] // 2 if grids[0].ndim == 2 else 0
    sl = np.array([(g[r] if g.ndim == 2 else g.ravel()) for g in grids])      # (T, X)
    T = sl.shape[0]
    if T > 400:
        sl = sl[np.linspace(0, T-1, 400).astype(int)]; T = 400
    vmin, vmax = float(sl.min()), float(sl.max())
    out = []
    for c in np.linspace(max(6, T // NF + 2), T, NF).astype(int):
        fig, ax = _new(title)
        ax.imshow(sl[:c], aspect="auto", cmap=cmap, vmin=vmin, vmax=vmax, interpolation="nearest",
                  origin="lower", extent=[0, sl.shape[1], 0, c])
        ax.set_xlim(0, sl.shape[1]); ax.set_ylim(0, T); ax.set_xlabel("space (lattice slice)"); ax.set_ylabel("time ->")
        out.append(_img(fig))
    return out

def pca_funnel(history, title, key="state", group="trial", **_):
    """2-D PCA of per-trial trajectories — many different initial conditions (blue) all funnel to
    the same final state (red): equifinality, which a distance-to-target line collapses away (P25)."""
    H = [f for f in history if isinstance(f, dict) and key in f and group in f]
    if len(H) < 4: return []
    tr = {}
    for f in H:
        tr.setdefault(int(f[group]), []).append(np.asarray(f[key], float).ravel())
    tr = {g: np.array(v) for g, v in list(tr.items()) if len(v) > 1}
    tr = dict(list(tr.items())[:12])
    if len(tr) < 2: return []
    allp = np.vstack(list(tr.values())); mu = allp.mean(0)
    Vt = np.linalg.svd(allp - mu, full_matrices=False)[2]; comp = Vt[:2]
    proj = {g: (v - mu) @ comp.T for g, v in tr.items()}
    P = np.vstack(list(proj.values())); lo, hi = P.min(0) - 0.5, P.max(0) + 0.5
    Lmax = max(len(p) for p in proj.values())
    out = []
    for k in np.linspace(0, Lmax-1, NF).astype(int):
        fig, ax = _new(title)
        for p in proj.values():
            kk = min(int(k), len(p)-1)
            ax.plot(p[:kk+1, 0], p[:kk+1, 1], lw=1.0, alpha=.65, color="#888")
            ax.scatter([p[0, 0]], [p[0, 1]], c="#2b6cb0", s=18, zorder=3)
            ax.scatter([p[kk, 0]], [p[kk, 1]], c="#e53e3e", s=12, zorder=4)
        ax.set_xlim(lo[0], hi[0]); ax.set_ylim(lo[1], hi[1]); ax.set_xlabel("PC1"); ax.set_ylabel("PC2")
        out.append(_img(fig))
    return out


def leadership(history, title, **_):
    """Heading arrows with the INFORMED minority ringed in black and their mean direction drawn as a
    bold arrow — shows a few informed individuals steering the whole group's travel direction (Couzin, P19)."""
    base = _has(history, "positions"); fr = _select(base, _s_points)
    if len(fr) < 2: return []
    allp = np.concatenate([np.asarray(f["positions"], float)[:, :2] for f in fr])
    lo, hi = allp.min(0) - 1, allp.max(0) + 1; span = float((hi - lo).mean())
    out = []
    for f in fr:
        p = np.asarray(f["positions"], float)[:, :2]; h = _headings(f)
        inf = np.asarray(f.get("informed_mask", np.zeros(len(p), bool))).astype(bool)
        fig, ax = _new(title)
        ax.quiver(p[:, 0], p[:, 1], np.cos(h), np.sin(h), (h % (2*np.pi)), cmap="hsv", clim=(0, 2*np.pi),
                  scale=22, scale_units="width", width=0.005, headwidth=3, pivot="mid", alpha=.7)
        if inf.any():
            ax.scatter(p[inf, 0], p[inf, 1], s=55, facecolors="none", edgecolors="black", linewidths=1.4, zorder=4)
            mh = np.angle(np.exp(1j*h[inf]).mean()); c = p.mean(0); s = span*0.22
            ax.annotate("", xy=(c[0]+s*np.cos(mh), c[1]+s*np.sin(mh)), xytext=(c[0], c[1]),
                        arrowprops=dict(arrowstyle="-|>", color="black", lw=2.2))
        ax.set_xlim(lo[0], hi[0]); ax.set_ylim(lo[1], hi[1]); ax.set_aspect("equal")
        ax.set_xticks([]); ax.set_yticks([]); out.append(_img(fig))
    return out

def minority_game(history, title, key="attendance", **_):
    """Attendance vs the resource capacity (N/2) with the random-choice baseline band — coordinated
    play keeps attendance in a TIGHTER band than random (better-than-chance resource use) (P23)."""
    H = [f for f in history if isinstance(f, dict) and key in f]
    if len(H) < 2: return []
    a = np.array([float(f[key]) for f in H]); N = float(H[0].get("n_agents", a.max()*2 or 2))
    cap = N/2; rstd = np.sqrt(N)/2; L = len(H)
    if L > 600:
        a = a[np.linspace(0, L-1, 600).astype(int)]; L = 600
    out = []
    for k in np.linspace(0, L-1, NF).astype(int):
        fig, ax = _new(title)
        ax.axhspan(cap-2*rstd, cap+2*rstd, color="#cbd5e0", alpha=.5, label="random ±2σ")
        ax.axhline(cap, ls="--", color="#e53e3e", lw=1.2, label="capacity N/2")
        ax.plot(a[:k+1], color="#2b6cb0", lw=1.0)
        ax.set_xlim(0, L); ax.set_ylim(cap-3*rstd, cap+3*rstd); ax.set_xlabel("round"); ax.set_ylabel("attendance")
        ax.legend(fontsize=7, loc="upper right"); out.append(_img(fig))
    return out

def gradient_climb(history, title, **_):
    """Swarm climbing a gradient: agents colored by their NOISY local reading, with the group
    centroid's trail (red) showing steady net translation up-gradient — 'many wrongs' average into
    accurate collective navigation no individual achieves (Berdahl, P17)."""
    H = [f for f in history if isinstance(f, dict) and "positions" in f]
    if len(H) < 3: return []
    fr = H if len(H) <= NF else [H[i] for i in np.linspace(0, len(H)-1, NF).astype(int)]
    allp = np.concatenate([np.asarray(f["positions"], float)[:, :2] for f in fr])
    lo, hi = allp.min(0) - 1, allp.max(0) + 1
    cents = np.array([np.asarray(f["positions"], float)[:, :2].mean(0) for f in fr])
    fs_all = np.concatenate([np.asarray(f.get("field_samples", np.zeros(len(np.asarray(f["positions"])))), float).ravel() for f in fr])
    vmin, vmax = float(fs_all.min()), float(fs_all.max())
    fc = fr[0].get("field_center"); fsig = float(fr[0].get("field_sigma", 4.0))   # draw the gradient source if known
    bg = None
    if fc is not None:
        gx, gy = np.meshgrid(np.linspace(lo[0], hi[0], 80), np.linspace(lo[1], hi[1], 80))
        bg = np.exp(-(((gx - fc[0])**2 + (gy - fc[1])**2) / (2 * fsig**2)))
    out = []
    for k, f in enumerate(fr):
        p = np.asarray(f["positions"], float)[:, :2]
        fsv = np.asarray(f.get("field_samples", np.zeros(len(p))), float).ravel()
        fig, ax = _new(title)
        if bg is not None:                              # gradient field: bright = source the group climbs toward
            ax.imshow(bg, extent=[lo[0], hi[0], lo[1], hi[1]], origin="lower", cmap="Greys", alpha=.45, aspect="auto", zorder=0)
        ax.scatter(p[:, 0], p[:, 1], c=fsv, cmap="viridis", vmin=vmin, vmax=vmax, s=18, zorder=2)
        ax.plot(cents[:k+1, 0], cents[:k+1, 1], "-", color="#e53e3e", lw=2)
        ax.scatter([cents[k, 0]], [cents[k, 1]], c="#e53e3e", s=70, marker="*", zorder=4)
        ax.set_xlim(lo[0], hi[0]); ax.set_ylim(lo[1], hi[1]); ax.set_aspect("equal")
        ax.set_xticks([]); ax.set_yticks([]); out.append(_img(fig))
    return out

def traffic_spacetime(history, title, bins=200, **_):
    """Space-time diagram of a 1-D road (position x time, color = speed) — a traffic jam is a
    band of stopped cars that propagates BACKWARD as a red diagonal stripe (Nagel-Schreckenberg, P8)."""
    H = [f for f in history if isinstance(f, dict) and "positions" in f]
    if len(H) < 3: return []
    L = float(H[0].get("L", max(np.asarray(H[0]["positions"]).max(), 1) + 1))
    vmax = float(max(np.asarray(f.get("velocities", [1])).max() for f in H)) or 1.0
    Hs = H if len(H) <= 400 else [H[i] for i in np.linspace(0, len(H)-1, 400).astype(int)]
    rows = []
    for f in Hs:
        p = np.asarray(f["positions"], float).ravel()
        v = np.asarray(f.get("velocities", np.full_like(p, vmax)), float).ravel()
        b = np.clip((p / L * bins).astype(int), 0, bins-1)
        row = np.full(bins, np.nan)
        for bi, vi in zip(b, v):
            row[bi] = vi if np.isnan(row[bi]) else min(row[bi], vi)   # slowest car in the bin (0 = jam)
        rows.append(row)
    img = np.array(rows); T = img.shape[0]
    out = []
    for c in np.linspace(max(6, T // NF + 2), T, NF).astype(int):
        fig, ax = _new(title)
        ax.imshow(img[:c], aspect="auto", cmap="RdYlGn", vmin=0, vmax=vmax, interpolation="nearest",
                  origin="lower", extent=[0, L, 0, c])
        ax.set_xlim(0, L); ax.set_ylim(0, T); ax.set_xlabel("position on road (red = stopped)"); ax.set_ylabel("time ->")
        out.append(_img(fig))
    return out

def density_field(history, title, bins=40, **_):
    """Coarse-grained DENSITY heatmap of continuous particle positions over time — motility-induced
    phase separation shows as a bright dense droplet emerging from a uniform dilute gas (MIPS, P2)."""
    H = [f for f in history if isinstance(f, dict) and "positions" in f]
    if len(H) < 3: return []
    fr = H if len(H) <= NF else [H[i] for i in np.linspace(0, len(H)-1, NF).astype(int)]
    allp = np.concatenate([np.asarray(f["positions"], float)[:, :2] for f in fr])
    L = float(fr[-1].get("box_size", allp.max() + 1))
    hists = []; hmax = 1e-9
    for f in fr:
        p = np.asarray(f["positions"], float)[:, :2]
        h = np.histogram2d(p[:, 0], p[:, 1], bins=bins, range=[[0, L], [0, L]])[0]
        hists.append(h); hmax = max(hmax, h.max())
    out = []
    for h in hists:
        fig, ax = _new(title)
        ax.imshow(h.T, origin="lower", cmap="inferno", vmin=0, vmax=hmax, interpolation="bilinear", extent=[0, L, 0, L])
        ax.set_xticks([]); ax.set_yticks([]); out.append(_img(fig))
    return out

def occupancy_field(history, title, **_):
    """Per-agent occupancy map — each cell colored by the agent that visits it most; distinct
    non-overlapping color regions are the exclusive home ranges of territoriality (P4)."""
    from matplotlib.colors import ListedColormap, BoundaryNorm
    H = [f for f in history if isinstance(f, dict) and "positions" in f]
    if len(H) < 3: return []
    G = int(H[0].get("grid_size", int(np.asarray(H[0]["positions"]).max()) + 1))
    na = int(np.asarray(H[0]["positions"]).shape[0])
    Hs = H if len(H) <= 300 else [H[i] for i in np.linspace(0, len(H)-1, 300).astype(int)]
    counts = np.zeros((na, G, G)); snaps = []; every = max(1, len(Hs) // NF)
    for t, f in enumerate(Hs):
        P = np.asarray(f["positions"], int)
        for a in range(min(na, len(P))):
            counts[a, P[a, 0] % G, P[a, 1] % G] += 1
        if t % every == 0 or t == len(Hs)-1:
            tot = counts.sum(0); dom = counts.argmax(0).astype(float); dom[tot == 0] = np.nan
            snaps.append(dom.copy())
    lcmap = ListedColormap([plt.cm.tab10(a % 10) for a in range(na)]); norm = BoundaryNorm(np.arange(-0.5, na+0.5, 1), na)
    out = []
    for dom in snaps[:NF]:
        fig, ax = _new(title)
        ax.imshow(dom.T, cmap=lcmap, norm=norm, interpolation="nearest", origin="lower")
        ax.set_xticks([]); ax.set_yticks([]); out.append(_img(fig))
    return out

def regulate_time(history, title, key="x", setpoint="setpoint", perturb="perturbation", **_):
    """State vs set-point over time, with perturbation onsets marked — shows homeostatic regulation
    driving the variable back toward its set-point after a kick (Ashby, P24)."""
    H = [f for f in history if isinstance(f, dict) and key in f]
    if len(H) < 2: return []
    x = np.array([float(f[key]) for f in H]); sp = np.array([float(f.get(setpoint, 0)) for f in H])
    pert = np.array([float(f.get(perturb, 0)) for f in H]); L = len(H)
    if L > 600:
        ix = np.linspace(0, L-1, 600).astype(int); x, sp, pert = x[ix], sp[ix], pert[ix]; L = 600
    pon = np.where(np.abs(np.diff(pert)) > 1e-9)[0]
    ymin = min(x.min(), sp.min()); ymax = max(x.max(), sp.max()); pad = 0.1*(ymax - ymin) + 1e-6
    out = []
    for k in np.linspace(0, L-1, NF).astype(int):
        fig, ax = _new(title)
        ax.plot(sp[:k+1], "--", color="#e53e3e", lw=1.3, label="set-point")
        ax.plot(x[:k+1], color="#2b6cb0", lw=1.6, label="state")
        for po in pon:
            if po <= k: ax.axvline(po, color="#b7791f", lw=0.9, alpha=.6)
        ax.set_xlim(0, L); ax.set_ylim(ymin-pad, ymax+pad); ax.set_xlabel("time"); ax.set_ylabel(key)
        ax.legend(fontsize=8, loc="upper right"); out.append(_img(fig))
    return out


PRODUCERS = {
    "point_cloud": point_cloud, "grid_field": grid_field, "phase_circle": phase_circle,
    "road_1d": road_1d, "vector_grid": vector_grid, "histogram_time": histogram_time,
    "lorenz_time": lorenz_time, "line_grow": line_grow, "multi_line": multi_line,
    "sweep_trace": sweep_trace, "dist_accumulate": dist_accumulate,
    "network_time": network_time, "bars_sort": bars_sort, "growing_spacetime": growing_spacetime,
    "pop_time": pop_time, "noise_response": noise_response, "phase_space": phase_space,
    "spacetime_slice": spacetime_slice, "pca_funnel": pca_funnel,
    "traffic_spacetime": traffic_spacetime, "occupancy_field": occupancy_field, "regulate_time": regulate_time,
    "leadership": leadership, "minority_game": minority_game, "gradient_climb": gradient_climb,
    "density_field": density_field,
}


def _motion(history, viz, args):
    """Mean normalized frame-to-frame change over the recorded run (0..~0.6)."""
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
# (QA defect class DC1). Each callable takes one history frame dict -> float (higher =
# more of the phenomenon present). Empty entries fall back to the producer default.
def _fld(f):
    return np.asarray(f.get("field", f.get("grid")), float)

def _m_turing(f):                       # P3: spectral power in the Turing wavelength band (rises as the maze forms)
    a = _fld(f); a = a - a.mean()
    if a.ndim != 2 or a.size < 16: return 0.0
    F = np.abs(np.fft.fftshift(np.fft.fft2(a))); ny, nx = a.shape
    yy, xx = np.ogrid[:ny, :nx]; r = np.sqrt((yy-ny//2)**2 + (xx-nx//2)**2)
    rmax = max(r.max(), 1e-9); band = (r > 0.12*rmax) & (r < 0.5*rmax)
    return float(F[band].sum() / (F.sum() + 1e-9))

def _m_mill(f):                         # P6: |normalized angular momentum| about the centroid (rises as the vortex forms)
    p = np.asarray(f["positions"], float)[:, :2]; h = _headings(f)
    if h is None or len(p) < 3: return 0.0
    r = p - p.mean(0); v = np.c_[np.cos(h), np.sin(h)]
    L = r[:, 0]*v[:, 1] - r[:, 1]*v[:, 0]
    return float(abs((L / (np.linalg.norm(r, axis=1) + 1e-9)).mean()))

def _m_lanes(f):                        # P7: lane order = mean per-y-strip directional purity (rises as counterflow segregates)
    p = np.asarray(f["positions"], float)[:, :2]; h = _headings(f)
    if h is None or len(p) < 6: return 0.0
    vx = np.cos(h); y = p[:, 1]; ymin, ymax = y.min(), y.max() + 1e-9
    bins = ((y - ymin)/(ymax - ymin)*8).astype(int).clip(0, 7)
    bias = [abs(np.sign(vx[bins == b]).mean()) for b in range(8) if (bins == b).sum() >= 3]
    return float(np.mean(bias)) if bias else 0.0

def _m_prune(f):                        # P29: active-edge count — RISES (mesh) then FALLS (prune), so equal-arc sampling catches both
    W = np.asarray(f["edge_weights"], float); m = max(W.max(), 1e-9)
    return float((W > 0.1*m).sum())

def _m_excited(f):                       # P13: excited-cell count — OSCILLATES with each wave, so equal-arc samples across the wave cycle
    if "excited_count" in f: return float(f["excited_count"])
    a = _fld(f); return float((a > 0).mean()) if a.size else 0.0

FRAME_MEASURE = {"P3": _m_turing, "P6": _m_mill, "P7": _m_lanes, "P29": _m_prune}
# P13 uses even sampling over a short formation-rich run (see gallery_runs._p13): the established
# excitable medium is a period-8 limit cycle (only ~8 distinct global states though every cell moves
# each step), so distinctness lives in the formation transient; the mp4 conveys the sustained sweep.
# spacetime_slice (1-D space-time raster) is reserved for genuinely 1-D systems (e.g. P8 traffic,
# where a backward jam wave reads as a clean diagonal) — a 1-D slice through 2-D spirals does not.


def render(history, viz, out_base, title, args=None, measure=None):
    global NF, MEASURE
    MEASURE = measure
    try:
        NF = _adaptive_nframes(history, viz, args or {})
        frames = PRODUCERS[viz](history, title, **(args or {}))
    finally:
        MEASURE = None
    return save_animation(frames, out_base, fps=max(10, min(20, NF // 5))) if frames else None
