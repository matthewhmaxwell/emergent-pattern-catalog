"""Generate figures for Sprint 26 replication output.

Two-panel figure:
  Left:  (A, beta) phase diagram from Phase 1m (heatmap of r_mean)
  Right: chimera basin fraction vs N at beta=0.18 from Phase 1k
"""
from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

OUT = Path("/home/claude/epc/analysis/outputs")

# ---- Phase 1m: phase diagram ------------------------------------------
pd_npz = np.load(OUT / "p10_phase_diagram.npz", allow_pickle=True)
A_grid = pd_npz["A_grid"]
beta_grid = pd_npz["beta_grid"]
r_mean = pd_npz["r_mean"]
classification = pd_npz["classification"]

# ---- Phase 1k: basin fraction summary ---------------------------------
with open(OUT / "p10_lifetime_replication.json") as f:
    lifetime_blob = json.load(f)
sumN = lifetime_blob["summary_per_N"]
Ns = sorted(int(k) for k in sumN.keys())
basin_fracs = [sumN[str(N)]["chimera_basin_fraction"] for N in Ns]
n_seeds = [sumN[str(N)]["n_seeds"] for N in Ns]
censored_fracs = [sumN[str(N)]["censored_fraction"] for N in Ns]
chim_n = [sumN[str(N)]["chimera_basin_count"] for N in Ns]

# 95% CI on basin fractions via Wilson score
def wilson(k, n, z=1.96):
    if n == 0:
        return 0.0, 0.0
    p = k / n
    denom = 1 + z*z/n
    centre = (p + z*z/(2*n)) / denom
    halfwidth = z * np.sqrt(p*(1-p)/n + z*z/(4*n*n)) / denom
    return centre - halfwidth, centre + halfwidth

ci = [wilson(chim_n[i], n_seeds[i]) for i in range(len(Ns))]
ci_lo = [c[0] for c in ci]
ci_hi = [c[1] for c in ci]
err_lo = [basin_fracs[i] - ci_lo[i] for i in range(len(Ns))]
err_hi = [ci_hi[i] - basin_fracs[i] for i in range(len(Ns))]

# ---- Plot -------------------------------------------------------------
fig, axes = plt.subplots(1, 2, figsize=(13, 5))

# Left: heatmap of r_mean. beta on x-axis, A on y-axis.
ax = axes[0]
im = ax.imshow(r_mean, origin="lower", aspect="auto",
               extent=[beta_grid[0] - 0.01, beta_grid[-1] + 0.01,
                       A_grid[0] - 0.005, A_grid[-1] + 0.005],
               vmin=0.4, vmax=1.0, cmap="viridis")
ax.set_xlabel(r"$\beta$ (phase-lag deviation $\pi/2 - \alpha$)")
ax.set_ylabel(r"$A$ (cosine-kernel amplitude)")
ax.set_title("Phase 1m: (A, β) phase diagram\nN=128, seed=42, T=50, single-seed classification")
cbar = fig.colorbar(im, ax=ax)
cbar.set_label(r"$\langle r_{\rm global} \rangle$ (last 10 frames)")

# Overlay chimera/sync boundary contour (where r drops out of chimera band)
# Chimera = r in (0.4, 0.85). Boundary near r=0.85.
ax.contour(beta_grid, A_grid, r_mean, levels=[0.85], colors="white",
           linewidths=2, linestyles="--")
ax.text(0.10, 0.965, "chimera", color="white", fontsize=11, fontweight="bold")
ax.text(0.26, 0.965, "sync", color="white", fontsize=11, fontweight="bold")
# Mark paper's beta=0.18
ax.axvline(0.18, color="orange", linewidth=1.2, alpha=0.8)
ax.text(0.18, 0.955, "  paper β=0.18", color="orange", fontsize=9)

# Right: basin fraction vs N
ax = axes[1]
ax.errorbar(Ns, basin_fracs, yerr=[err_lo, err_hi], fmt="o-",
            capsize=4, markersize=8, color="C0", label="chimera basin fraction")
for i, N in enumerate(Ns):
    ax.annotate(f"{chim_n[i]}/{n_seeds[i]}",
                (N, basin_fracs[i]),
                textcoords="offset points", xytext=(8, -2),
                fontsize=9, color="C0")
ax.axhline(0.5, color="gray", linewidth=0.5, linestyle=":")
ax.set_xscale("log", base=2)
ax.set_xticks(Ns)
ax.set_xticklabels([str(N) for N in Ns])
ax.set_xlabel(r"$N$ (oscillators)")
ax.set_ylabel("fraction of seeds in chimera basin")
ax.set_title("Phase 1k: chimera basin fraction at β=0.18\n"
             "(all surviving seeds censored at T=100)")
ax.set_ylim(0, 1)
ax.grid(True, alpha=0.3)
# Annotate that lifetimes are right-censored
ax.text(0.02, 0.04,
        "All chimera-basin seeds survived T_max=100\n"
        "(0/47 transitions to sync observed)",
        transform=ax.transAxes, fontsize=9,
        bbox=dict(boxstyle="round,pad=0.4", facecolor="wheat", alpha=0.8))

plt.suptitle("Sprint 26 — P10 β=0.18 paper-parameter replication "
             "(Abrams-Strogatz 2004 Fig. 2)", fontsize=12, y=1.02)
plt.tight_layout()

out = OUT / "p10_phase_diagram.png"
plt.savefig(out, dpi=130, bbox_inches="tight")
print(f"Wrote {out}")
