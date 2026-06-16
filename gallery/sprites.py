"""Sprite-sheet renderers (v2): velocity-aware. Motion patterns are drawn as
arrows coloured by heading so ALIGNMENT is visible (a flock = one colour, all
arrows parallel; disorder = rainbow). Frame sampling concentrates on where the
emergent structure is actually changing."""
from __future__ import annotations
import io
import numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image

FW = FH = 300
MAXF = 36
DPI = 75


def _fig_img(fig):
    buf = io.BytesIO(); fig.savefig(buf, format="png", dpi=DPI, facecolor="white")
    plt.close(fig); buf.seek(0)
    im = Image.open(buf).convert("RGB")
    return im.resize((FW, FH)) if im.size != (FW, FH) else im


# ---- structure per frame (kind-appropriate) for sampling + rise -------------
def _headings(f):
    if "headings" in f: return np.asarray(f["headings"], float).ravel()
    if "velocities" in f:
        v = np.asarray(f["velocities"], float); return np.arctan2(v[:, 1], v[:, 0])
    return None

def _struct(kind, f):
    try:
        if kind == "point_cloud":
            h = _headings(f)
            if h is not None and h.size:
                return float(abs(np.exp(1j * h).mean()))   # polar order (alignment)
            p = np.asarray(f["positions"], float)[:, :2]    # else spatial clustering
            n = len(p)
            if n < 3: return 0.0
            d = np.sqrt(((p[:, None] - p[None]) ** 2).sum(-1)); np.fill_diagonal(d, np.inf)
            ext = max((np.ptp(p[:, 0]) + np.ptp(p[:, 1])) / 2, 1e-9)
            return float(ext / (d.min(1).mean() + 1e-9))
        if kind == "grid_field":
            key = "field" if "field" in f else "grid"
            a = np.asarray(f[key], float); a = a - a.mean(); den = (a * a).sum()
            if den <= 1e-12: return 0.0
            num = (a[:-1] * a[1:]).sum() + (a[:, :-1] * a[:, 1:]).sum()
            return float((num / (a[:-1].size + a[:, :-1].size)) / (den / a.size))
        if kind == "phase_circle":
            t = np.asarray(f.get("theta", f.get("phases")), float)
            return float(abs(np.exp(1j * t).mean()))
    except Exception:
        return np.nan
    return np.nan


def _frames(history, key, kind):
    fr = [f for f in history if isinstance(f, dict) and key in f]
    if len(fr) <= MAXF:
        return fr
    try:
        S = np.array([_struct(kind, f) for f in fr], float)
        good = np.isfinite(S)
        if good.sum() >= 3 and np.nanstd(S[good]) > 1e-9:
            S = np.interp(np.arange(len(S)), np.where(good)[0], S[good])
            C = np.concatenate([[0.0], np.cumsum(np.abs(np.diff(S)))])
            C = 0.7 * (C / C[-1]) + 0.3 * (np.arange(len(C)) / (len(C) - 1))
            idx = np.unique(np.searchsorted(C, np.linspace(0, C[-1], MAXF)).clip(0, len(fr) - 1))
            return [fr[i] for i in idx]
    except Exception:
        pass
    return [fr[i] for i in np.linspace(0, len(fr) - 1, MAXF).astype(int)]


def _imgs_point_cloud(history, title):
    fr = _frames(history, "positions", "point_cloud")
    if len(fr) < 2: return []
    allp = np.concatenate([np.asarray(f["positions"], float)[:, :2] for f in fr])
    lo, hi = allp.min(0) - 1, allp.max(0) + 1
    use_arrows = _headings(fr[0]) is not None
    has_t = "types" in fr[0] and not use_arrows
    out = []
    for f in fr:
        p = np.asarray(f["positions"], float)[:, :2]
        fig, ax = plt.subplots(figsize=(4, 4), dpi=DPI)
        if use_arrows:
            h = _headings(f)
            ax.quiver(p[:, 0], p[:, 1], np.cos(h), np.sin(h), (h % (2*np.pi)),
                      cmap="hsv", clim=(0, 2*np.pi), scale=22, scale_units="width",
                      width=0.006, headwidth=3, pivot="mid")
        else:
            ax.scatter(p[:, 0], p[:, 1], s=16,
                       c=(np.asarray(f["types"], float) if has_t else "#2b6cb0"),
                       cmap="tab10" if has_t else None)
        ax.set_xlim(lo[0], hi[0]); ax.set_ylim(lo[1], hi[1]); ax.set_aspect("equal")
        ax.set_xticks([]); ax.set_yticks([]); ax.set_title(title, fontsize=9)
        out.append(_fig_img(fig))
    return out


def _imgs_grid_field(history, title):
    key = "field" if (history and "field" in history[-1]) else "grid"
    fr = _frames(history, key, "grid_field")
    if len(fr) < 2: return []
    allv = np.concatenate([np.asarray(f[key], float).ravel() for f in fr])
    vmin, vmax = float(allv.min()), float(allv.max()); discrete = key == "grid"
    out = []
    for f in fr:
        fig, ax = plt.subplots(figsize=(4, 4), dpi=DPI)
        ax.imshow(np.asarray(f[key], float), cmap=("tab10" if discrete else "viridis"),
                  vmin=vmin, vmax=vmax, interpolation="nearest")
        ax.set_xticks([]); ax.set_yticks([]); ax.set_title(title, fontsize=9)
        out.append(_fig_img(fig))
    return out


def _imgs_phase_circle(history, title):
    key = "theta" if (history and "theta" in history[-1]) else "phases"
    fr = _frames(history, key, "phase_circle")
    if len(fr) < 2: return []
    out = []
    for f in fr:
        th = np.asarray(f[key], float)
        fig, ax = plt.subplots(figsize=(4, 4), dpi=DPI)
        ax.add_patch(plt.Circle((0, 0), 1, fill=False, color="#ccc"))
        ax.scatter(np.cos(th), np.sin(th), s=18, c="#2b6cb0")
        z = np.exp(1j * th).mean()
        ax.plot([0, z.real], [0, z.imag], "-", color="#e53e3e", lw=2)
        ax.set_xlim(-1.2, 1.2); ax.set_ylim(-1.2, 1.2); ax.set_aspect("equal")
        ax.set_xticks([]); ax.set_yticks([]); ax.set_title(title, fontsize=9)
        out.append(_fig_img(fig))
    return out


_KIND = {"point_cloud": _imgs_point_cloud, "grid_field": _imgs_grid_field,
         "phase_circle": _imgs_phase_circle}


def _structure_rise(history, kind):
    key = {"point_cloud": "positions",
           "grid_field": ("field" if (history and "field" in history[-1]) else "grid"),
           "phase_circle": ("theta" if (history and "theta" in history[-1]) else "phases")}[kind]
    S = np.array([_struct(kind, f) for f in history if isinstance(f, dict) and key in f], float)
    S = S[np.isfinite(S)]
    if S.size < 2: return {}
    smin, smax = float(S.min()), float(S.max()); rng = smax - smin
    start_frac = (float(S[0]) - smin) / rng if rng > 1e-9 else 1.0
    return {"s0": round(float(S[0]), 3), "send": round(float(S[-1]), 3),
            "smin": round(smin, 3), "smax": round(smax, 3),
            "start_frac": round(start_frac, 2), "flat": bool(rng < 1e-6 or start_frac > 0.6)}


def render_sprite(history, kind, out_path, title=""):
    imgs = _KIND[kind](history, title)
    if len(imgs) < 2: return None
    n = len(imgs); cols = int(np.ceil(np.sqrt(n))); rows = int(np.ceil(n / cols))
    sheet = Image.new("RGB", (cols * FW, rows * FH), "white")
    for i, im in enumerate(imgs):
        sheet.paste(im, ((i % cols) * FW, (i // cols) * FH))
    sheet.save(out_path, optimize=True)
    return {"frames": n, "cols": cols, "rows": rows, "fw": FW, "fh": FH,
            "rise": _structure_rise(history, kind)}
