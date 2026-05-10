"""Sprint 27 figure: chimera basin fraction heatmap on (A, β) grid.

Replots Sprint 26's phase diagram with the new multi-seed basin-volume
data overlaid. Shows the smooth gradient that single-seed sampling
missed.
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

OUT = Path("/home/claude/epc/analysis/outputs")

with open(OUT / "p10_phase_boundary_multiseed.json") as f:
    blob = json.load(f)

A_grid = np.array(blob["config"]["A_grid"])
betas = sorted(blob["config"]["beta_bulk"] + blob["config"]["beta_boundary"])
basin = np.zeros((len(A_grid), len(betas)))
for i, A in enumerate(A_grid):
    for j, beta in enumerate(betas):
        key = f"A={A:.2f}_beta={beta:.2f}"
        basin[i, j] = blob["grid"][key]["basin_fraction_chimera"]

fig, ax = plt.subplots(figsize=(9, 5))
im = ax.imshow(basin, origin="lower", aspect="auto",
               extent=[-0.5, len(betas) - 0.5, A_grid[0] - 0.005,
                       A_grid[-1] + 0.005],
               vmin=0.0, vmax=1.0, cmap="RdYlBu_r")
ax.set_xticks(range(len(betas)))
ax.set_xticklabels([f"{b:.2f}" for b in betas])
ax.set_xlabel(r"$\beta$ (phase-lag deviation $\pi/2 - \alpha$)")
ax.set_ylabel(r"$A$ (cosine-kernel amplitude)")
ax.set_title("Sprint 27 Phase 1n: chimera basin fraction (5 seeds per cell)\n"
             f"N={blob['config']['N']}, T={blob['config']['T']}, "
             f"IC={blob['config']['init_mode']}")

# Annotate each cell with k/5
for i, A in enumerate(A_grid):
    for j, beta in enumerate(betas):
        key = f"A={A:.2f}_beta={beta:.2f}"
        cell = blob["grid"][key]
        nc = cell["classifications"].count("chimera")
        ns = cell["classifications"].count("sync")
        text_color = "white" if 0.3 < basin[i, j] < 0.75 else "black"
        ax.text(j, A, f"{nc}/{ns}", ha="center", va="center",
                fontsize=9, color=text_color, fontweight="bold")

cbar = fig.colorbar(im, ax=ax)
cbar.set_label("chimera basin fraction (n_chimera / 5 seeds)")

# Legend
ax.text(0.99, 1.02, "cell labels: n_chimera / n_sync (out of 5)",
        transform=ax.transAxes, ha="right", fontsize=9, style="italic")

# Sprint 26 single-seed boundary annotation
ax.text(0.5, -0.18,
        "Sprint 26 single-seed (seed=42) classified β=0.18, 0.20 as "
        "chimera_only and β=0.22 as sync_only.\n"
        "Sprint 27 multi-seed shows continuous basin-volume gradient: "
        "β=0.18 → 60%, β=0.20 → 40%, β=0.22 → 0–20%.",
        transform=ax.transAxes, ha="center", va="top", fontsize=9)

plt.tight_layout()
out = OUT / "p10_basin_volume_multiseed.png"
plt.savefig(out, dpi=130, bbox_inches="tight")
print(f"Wrote {out}")
