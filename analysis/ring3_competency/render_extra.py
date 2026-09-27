"""Lightweight colored-marker renderer for competencies whose visuals need arbitrary per-marker colors/shapes
(metaworld types, abstraction colours, etc.) — a companion to gallery/mk_ring3_demo.render_ring3. Same visual
language (white bg, no ticks, small title, faint agent trails, bottom legend) -> sprite + mp4 + gif."""
import sys, numpy as np
import matplotlib; matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from PIL import Image
ROOT = "/home/matthewhmaxwell/emergent-pattern-catalog"; sys.path.insert(0, ROOT)
from gallery.anim import save_animation
FW = FH = 300


def _pil(fig):
    fig.canvas.draw(); w, h = fig.canvas.get_width_height()
    buf = np.frombuffer(fig.canvas.buffer_rgba(), np.uint8).reshape(h, w, 4)[..., :3]
    plt.close(fig); return Image.fromarray(buf).resize((FW, FH))


def render_markers(frames, out_base, title, legend, extent, fps=8):
    """frames: list of {"markers":[(x,y,hexcolor,shape,filled_bool,size),...], "agents":[(x,y,hexcolor),...]}."""
    trails = {}; pil = []
    for fr in frames:
        fig, ax = plt.subplots(figsize=(3, 3), dpi=100)
        ax.set_xlim(extent[0], extent[1]); ax.set_ylim(extent[2], extent[3]); ax.set_aspect("equal")
        ax.set_xticks([]); ax.set_yticks([]); ax.set_title(title, fontsize=7.5)
        for (x, y, col, shape, filled, size) in fr.get("markers", []):
            ax.scatter([x], [y], s=size, marker=shape, facecolors=(col if filled else "none"),
                       edgecolors=col, linewidths=1.6, zorder=2)
        for i, (x, y, col) in enumerate(fr.get("agents", [])):
            trails.setdefault(i, []).append((x, y))
            xs = [p[0] for p in trails[i]][-16:]; ys = [p[1] for p in trails[i]][-16:]
            ax.plot(xs, ys, "-", color=col, alpha=0.28, lw=1.5, zorder=4)
            ax.scatter([x], [y], s=95, color=col, edgecolors="white", linewidths=1.0, zorder=5)
        hud = fr.get("hud")
        if hud is not None:
            hc = fr.get("hud_color")
            if hc:
                ax.add_patch(plt.Rectangle((0.035, 0.915), 0.055, 0.055, transform=ax.transAxes,
                             facecolor=hc, edgecolor="#333", lw=0.7, zorder=6, clip_on=False))
            ax.text(0.11 if hc else 0.035, 0.955, hud, transform=ax.transAxes, fontsize=7,
                    va="top", ha="left", color="#222", zorder=6)
        if legend:
            hs = [Line2D([0], [0], marker=mk, color="none", markerfacecolor=fc, markeredgecolor=ec,
                         markersize=8, label=lab, lw=0) for (lab, mk, fc, ec) in legend]
            ax.legend(handles=hs, loc="upper center", bbox_to_anchor=(0.5, -0.01), ncol=min(3, len(legend)),
                      fontsize=6, frameon=False, handletextpad=0.3, columnspacing=0.8)
        fig.tight_layout(pad=0.4); pil.append(_pil(fig))
    spr = save_animation(pil, out_base, fps=fps)
    try:
        import imageio.v2 as imageio
        imageio.mimwrite(out_base + ".gif", [np.asarray(p) for p in pil], duration=0.12, loop=0)
    except Exception as e:
        print("gif warn:", e)
    return spr
