"""Sprint 26 Phase 1m — (A, β) phase-space boundary at fixed seed.

Reproduces the phase-diagram aspect of Abrams-Strogatz 2004 Fig. 2.

Setup: N=128, IC=asymmetric_gaussian, seed=42, T=50, dt=0.025.
Grid: A in {0.95, 0.96, ..., 1.00} × beta in {0.00, 0.02, ..., 0.30}.
Classification per cell:
    sync       : r_global -> 1.0  (mean over last 10 frames > 0.95)
    chimera    : r_global in (0.4, 0.85), low std
    incoherent : r_global -> 0    (mean < 0.2)
    other      : anything else (transient, drift, undecided)

Output: analysis/outputs/p10_phase_diagram.json (classification grid +
metadata) and .npz (grids for plotting).

This is a coarse-grained boundary map, not a basin-volume measurement
(single seed per cell). Topology (location and shape of chimera region)
is the bar; quantitative basin-volume agreement is a stretch goal.
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from epc.models.kuramoto_nonlocal import KuramotoNonlocal


def classify(r_tail: np.ndarray) -> str:
    r_mean = float(np.mean(r_tail))
    r_std = float(np.std(r_tail))
    if r_mean > 0.95:
        return "sync"
    if r_mean < 0.2:
        return "incoherent"
    if 0.4 < r_mean < 0.85 and r_std < 0.10:
        return "chimera"
    return "other"


def main():
    A_grid = np.array([0.95, 0.96, 0.97, 0.98, 0.99, 1.00])
    beta_grid = np.array([0.00, 0.02, 0.04, 0.06, 0.08, 0.10, 0.12,
                          0.14, 0.16, 0.18, 0.20, 0.22, 0.24, 0.26,
                          0.28, 0.30])

    N = 128
    seed = 42
    T = 50
    n_frames = T  # record_dt = 1.0
    tail = 10  # average r over last 10 frames

    n_A, n_b = len(A_grid), len(beta_grid)
    r_mean_grid = np.zeros((n_A, n_b))
    r_std_grid = np.zeros((n_A, n_b))
    classification = np.empty((n_A, n_b), dtype=object)

    t_start = time.time()
    print(f"Phase 1m: scanning {n_A}x{n_b} = {n_A*n_b} cells "
          f"at N={N}, seed={seed}, T={T}")
    for i, A in enumerate(A_grid):
        for j, beta in enumerate(beta_grid):
            t0 = time.time()
            m = KuramotoNonlocal(N=N, A=float(A), beta=float(beta),
                                 seed=seed, init_mode="asymmetric_gaussian")
            hist = m.run(n_frames=n_frames)
            r_tail = np.array([h["r"] for h in hist[-tail:]])
            r_mean_grid[i, j] = np.mean(r_tail)
            r_std_grid[i, j] = np.std(r_tail)
            cls = classify(r_tail)
            classification[i, j] = cls
            print(f"  A={A:.2f} beta={beta:.2f}: "
                  f"r_mean={r_mean_grid[i,j]:.3f} r_std={r_std_grid[i,j]:.3f} "
                  f"-> {cls:10s}  [{time.time()-t0:.1f}s]")

    elapsed = time.time() - t_start
    print(f"\nTotal: {elapsed:.1f}s ({elapsed/(n_A*n_b):.2f}s/cell)")

    # Counts
    counts = {}
    for cls in ["sync", "chimera", "incoherent", "other"]:
        counts[cls] = int(np.sum(classification == cls))
    print(f"\nClassification counts: {counts}")

    # Identify chimera region bounding box
    chim_mask = (classification == "chimera")
    if chim_mask.any():
        ai, bj = np.where(chim_mask)
        a_lo, a_hi = A_grid[ai.min()], A_grid[ai.max()]
        b_lo, b_hi = beta_grid[bj.min()], beta_grid[bj.max()]
        print(f"Chimera region bounding box: "
              f"A in [{a_lo:.2f}, {a_hi:.2f}], "
              f"beta in [{b_lo:.2f}, {b_hi:.2f}]")
    else:
        a_lo = a_hi = b_lo = b_hi = None
        print("No chimera cells found.")

    out_dir = Path("/home/claude/epc/analysis/outputs")
    out_dir.mkdir(parents=True, exist_ok=True)

    np.savez(out_dir / "p10_phase_diagram.npz",
             A_grid=A_grid, beta_grid=beta_grid,
             r_mean=r_mean_grid, r_std=r_std_grid,
             classification=classification.astype(str))

    with open(out_dir / "p10_phase_diagram.json", "w") as f:
        json.dump({
            "sprint": 26,
            "phase": "1m",
            "purpose": "(A, beta) phase-space boundary, single seed",
            "config": {
                "N": N, "seed": seed, "T": T, "dt": 0.025,
                "record_dt": 1.0, "init_mode": "asymmetric_gaussian",
                "tail_frames_for_classification": tail,
            },
            "A_grid": A_grid.tolist(),
            "beta_grid": beta_grid.tolist(),
            "r_mean": r_mean_grid.tolist(),
            "r_std": r_std_grid.tolist(),
            "classification": classification.tolist(),
            "counts": counts,
            "chimera_bounding_box": {
                "A_lo": float(a_lo) if a_lo is not None else None,
                "A_hi": float(a_hi) if a_hi is not None else None,
                "beta_lo": float(b_lo) if b_lo is not None else None,
                "beta_hi": float(b_hi) if b_hi is not None else None,
            } if chim_mask.any() else None,
            "wall_time_s": elapsed,
        }, f, indent=2)
    print(f"\nWrote {out_dir / 'p10_phase_diagram.json'}")
    print(f"Wrote {out_dir / 'p10_phase_diagram.npz'}")


if __name__ == "__main__":
    main()
