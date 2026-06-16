"""Sprite-sheet renderers for the animated morphologies -> one PNG of tiled
frames + metadata, so the gallery playbar can step/play/scrub frame-by-frame."""
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
    buf = io.BytesIO(); fig.savefig(buf, format="png", dpi=DPI,
                                    facecolor="white"); plt.close(fig); buf.seek(0)
    im = Image.open(buf).convert("RGB")
    return im.resize((FW, FH)) if im.size != (FW, FH) else im


def _frames(history, key):
    fr = [f for f in history if isinstance(f, dict) and key in f]
    if len(fr) > MAXF:
        idx = np.linspace(0, len(fr) - 1, MAXF).astype(int)
        fr = [fr[i] for i in idx]
    return fr


def _imgs_point_cloud(history, title):
    fr = _frames(history, "positions")
    if len(fr) < 2: return []
    allp = np.concatenate([np.asarray(f["positions"], float)[:, :2] for f in fr])
    lo, hi = allp.min(0) - 1, allp.max(0) + 1
    has_t = "types" in fr[0]
    out = []
    for f in fr:
        p = np.asarray(f["positions"], float)[:, :2]
        fig, ax = plt.subplots(figsize=(4, 4), dpi=DPI)
        ax.scatter(p[:, 0], p[:, 1], s=14,
                   c=(np.asarray(f["types"], float) if has_t else "#2b6cb0"),
                   cmap="tab10" if has_t else None)
        ax.set_xlim(lo[0], hi[0]); ax.set_ylim(lo[1], hi[1]); ax.set_aspect("equal")
        ax.set_xticks([]); ax.set_yticks([]); ax.set_title(title, fontsize=9)
        out.append(_fig_img(fig))
    return out


def _imgs_grid_field(history, title):
    key = "field" if (history and "field" in history[-1]) else "grid"
    fr = _frames(history, key)
    if len(fr) < 2: return []
    allv = np.concatenate([np.asarray(f[key], float).ravel() for f in fr])
    vmin, vmax = float(allv.min()), float(allv.max())
    discrete = key == "grid"
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
    fr = _frames(history, key)
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


def render_sprite(history, kind, out_path, title=""):
    imgs = _KIND[kind](history, title)
    if len(imgs) < 2:
        return None
    n = len(imgs)
    cols = int(np.ceil(np.sqrt(n))); rows = int(np.ceil(n / cols))
    sheet = Image.new("RGB", (cols * FW, rows * FH), "white")
    for i, im in enumerate(imgs):
        sheet.paste(im, ((i % cols) * FW, (i // cols) * FH))
    sheet.save(out_path, optimize=True)
    return {"frames": n, "cols": cols, "rows": rows, "fw": FW, "fh": FH}
