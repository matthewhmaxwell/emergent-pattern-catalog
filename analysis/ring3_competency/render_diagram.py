"""Non-spatial diagram animators for the competencies that have no grid (culture ratchet, reputation timeline,
fairness turn-taking, compositional message). Same output contract as the spatial renderers: sprite + mp4 + gif."""
import sys, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image
ROOT = "/home/matthewhmaxwell/emergent-pattern-catalog"; sys.path.insert(0, ROOT)
from gallery.anim import save_animation
FW = FH = 300


def _pil(fig):
    fig.canvas.draw(); w, h = fig.canvas.get_width_height()
    buf = np.frombuffer(fig.canvas.buffer_rgba(), np.uint8).reshape(h, w, 4)[..., :3]
    plt.close(fig); return Image.fromarray(buf).resize((FW, FH))


def line_reveal(xs, series, title, out_base, xlabel, ylabel, fps=6, ymax=None, ymin=0, hline=None, vmark=None, vlabel="", styles=None):
    """series: list of (label, yvals, color). Reveal points progressively. Optional setpoint hline + a vertical
    event marker (vmark = x-value, e.g. the cash-out step). Optional per-series `styles`: list of dicts with any
    of {ls, lw, ms, marker, zorder, drawstyle} — lets a reference series be drawn dashed/underneath so an
    overlapping series stays visible (default: solid line + round markers for every series)."""
    ymax = ymax or max(max(yv) for _, yv, _ in series) * 1.12
    pil = []
    for k in range(2, len(xs) + 1):
        fig, ax = plt.subplots(figsize=(3, 3), dpi=100)
        for i, (lab, yv, col) in enumerate(series):
            st = (styles[i] if styles and i < len(styles) and styles[i] else {})
            ax.plot(xs[:k], yv[:k], color=col, label=lab, ls=st.get("ls", "-"),
                    lw=st.get("lw", 1.7), marker=st.get("marker", "o"), ms=st.get("ms", 2.5),
                    zorder=st.get("zorder", 3), drawstyle=st.get("drawstyle", "default"))
        if hline is not None:
            ax.axhline(hline, color="#888", ls="--", lw=1.0)
        if vmark is not None and xs[k - 1] >= vmark:
            ax.axvline(vmark, color="#e53e3e", lw=1.6)
            ax.text(vmark, ymax * 0.92, vlabel, fontsize=6, color="#e53e3e", ha="center")
        ax.set_xlim(xs[0], xs[-1]); ax.set_ylim(ymin, ymax)
        ax.set_title(title, fontsize=7.5); ax.set_xlabel(xlabel, fontsize=7); ax.set_ylabel(ylabel, fontsize=7)
        ax.tick_params(labelsize=6); ax.legend(fontsize=6, frameon=False, loc="upper left")
        fig.tight_layout(pad=0.5); pil.append(_pil(fig))
    pil += [pil[-1]] * 3
    spr = save_animation(pil, out_base, fps=fps)
    try:
        import imageio.v2 as imageio
        imageio.mimwrite(out_base + ".gif", [np.asarray(p) for p in pil], duration=0.16, loop=0)
    except Exception as e:
        print("gif warn:", e)
    return spr


def timeline(steps, title, out_base, legend, fps=4):
    """steps: list of frames; each frame = list of cells (col, row, hexcolor, text, textcolor). Reveal columns."""
    ncol = max((c[0] for fr in steps for c in fr), default=0) + 1
    nrow = max((c[1] for fr in steps for c in fr), default=0) + 1
    pil = []
    for fr in steps:
        fig, ax = plt.subplots(figsize=(3, 3), dpi=100)
        ax.set_xlim(-.5, ncol - .5); ax.set_ylim(-.8, nrow - .5); ax.set_aspect("auto")
        ax.set_xticks([]); ax.set_yticks([]); ax.set_title(title, fontsize=7.2)
        for (col, row, hexc, text, tcol) in fr:
            ax.add_patch(plt.Rectangle((col - .42, row - .42), .84, .84, facecolor=hexc, edgecolor="#444", lw=0.6))
            if text:
                ax.text(col, row, text, ha="center", va="center", fontsize=6.5, color=tcol)
        if legend:
            from matplotlib.patches import Patch
            hs = [Patch(facecolor=fc, edgecolor="#444", label=lab) for (lab, fc) in legend]
            ax.legend(handles=hs, loc="upper center", bbox_to_anchor=(0.5, -0.02), ncol=min(4, len(legend)),
                      fontsize=6, frameon=False, handletextpad=0.4, columnspacing=0.8)
        fig.tight_layout(pad=0.5); pil.append(_pil(fig))
    pil += [pil[-1]] * 3
    spr = save_animation(pil, out_base, fps=fps)
    try:
        import imageio.v2 as imageio
        imageio.mimwrite(out_base + ".gif", [np.asarray(p) for p in pil], duration=0.22, loop=0)
    except Exception as e:
        print("gif warn:", e)
    return spr
