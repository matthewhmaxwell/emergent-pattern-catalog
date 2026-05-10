"""Sprint 27 Phase 1n — multi-seed (A, β) phase boundary refinement.

Sprint 26 Phase 1m showed a sharp boundary at β ≈ 0.21 (every cell at
β ≤ 0.20 chimera, every cell at β ≥ 0.22 sync) but with single-seed
classification per cell. This phase samples 5 seeds per cell across two
strips:

  Boundary strip:  β ∈ {0.18, 0.20, 0.22}, A ∈ {0.95..1.00}, 6×3=18 cells
  Bulk strip:      β ∈ {0.05, 0.10},       A ∈ {0.95..1.00}, 6×2=12 cells

Total: 30 cells × 5 seeds = 150 runs. At ~3 s/run for N=128, T=50
≈ 7-8 min wall. Fits inside one timeout.

For each cell we report:
  - basin fraction (chimera basin / 5)
  - r_mean per seed (last-10-frame average)
  - per-cell classification: chimera-only, sync-only, mixed, other

Output: analysis/outputs/p10_phase_boundary_multiseed.json (and .npz).

This addresses Sprint 26 risk-register item #31: "phase diagram is
single-seed; a more rigorous diagram would sample basin volume at
each (A, β) cell with multiple seeds."
"""
from __future__ import annotations

import json
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from epc.models.kuramoto_nonlocal import KuramotoNonlocal


CHIM_LO = 0.4
CHIM_HI = 0.85
SYNC_THRESH = 0.95


def classify_run(r_tail: np.ndarray) -> str:
    r_mean = float(np.mean(r_tail))
    r_std = float(np.std(r_tail))
    if r_mean > SYNC_THRESH:
        return "sync"
    if r_mean < 0.2:
        return "incoherent"
    if CHIM_LO < r_mean < CHIM_HI and r_std < 0.10:
        return "chimera"
    return "other"


def main():
    A_grid = np.array([0.95, 0.96, 0.97, 0.98, 0.99, 1.00])
    beta_boundary = [0.18, 0.20, 0.22]
    beta_bulk = [0.05, 0.10]
    seeds = [0, 1, 42, 200, 500]
    N = 128
    n_frames = 50  # T = 50
    tail = 10

    cells = []
    for beta in beta_bulk + beta_boundary:
        for A in A_grid:
            cells.append((float(A), float(beta)))

    print(f"Sprint 27 Phase 1n: {len(cells)} cells x {len(seeds)} seeds = "
          f"{len(cells)*len(seeds)} runs at N={N}, T={n_frames}")

    out_dir = Path("/home/claude/epc/analysis/outputs")
    out_dir.mkdir(parents=True, exist_ok=True)
    incr_path = out_dir / "p10_phase_boundary_multiseed_incr.json"

    # Resume from incremental file if it exists
    grid: dict = {}
    if incr_path.exists():
        with open(incr_path) as f:
            grid = json.load(f).get("grid", {})
        print(f"Resuming with {len(grid)} cells already done")

    t_start = time.time()
    for A, beta in cells:
        cell_key = f"A={A:.2f}_beta={beta:.2f}"
        if cell_key in grid:
            print(f"  {cell_key}: SKIP (cached)")
            continue
        cell_key = f"A={A:.2f}_beta={beta:.2f}"
        cell = {"A": A, "beta": beta, "seeds": {}, "classifications": []}
        for seed in seeds:
            t0 = time.time()
            m = KuramotoNonlocal(N=N, A=A, beta=beta, seed=seed,
                                 init_mode="asymmetric_gaussian")
            hist = m.run(n_frames=n_frames)
            r_tail = np.array([h["r"] for h in hist[-tail:]])
            cls = classify_run(r_tail)
            cell["seeds"][str(seed)] = {
                "r_mean": float(np.mean(r_tail)),
                "r_std": float(np.std(r_tail)),
                "classification": cls,
                "wall_s": time.time() - t0,
            }
            cell["classifications"].append(cls)

        # Cell-level summary
        n_chim = cell["classifications"].count("chimera")
        n_sync = cell["classifications"].count("sync")
        n_other = len(seeds) - n_chim - n_sync
        cell["basin_fraction_chimera"] = n_chim / len(seeds)
        cell["basin_fraction_sync"] = n_sync / len(seeds)
        if n_chim == len(seeds):
            cell["cell_class"] = "chimera_only"
        elif n_sync == len(seeds):
            cell["cell_class"] = "sync_only"
        elif n_other > 0:
            cell["cell_class"] = "mixed_with_other"
        else:
            cell["cell_class"] = "mixed_chimera_sync"
        grid[cell_key] = cell
        # Incremental save after each cell so partial runs survive timeouts
        with open(incr_path, "w") as f:
            json.dump({"grid": grid}, f)
        print(f"  {cell_key}: chimera={n_chim}/5, sync={n_sync}/5, "
              f"other={n_other}/5  [{cell['cell_class']}]")

    elapsed = time.time() - t_start
    print(f"\nTotal: {elapsed:.1f}s "
          f"({elapsed/(len(cells)*len(seeds)):.2f}s/run)")

    # Aggregate
    cell_class_counts = {}
    for cell in grid.values():
        cell_class_counts[cell["cell_class"]] = (
            cell_class_counts.get(cell["cell_class"], 0) + 1
        )
    print(f"\nCell classification counts: {cell_class_counts}")

    # Bulk vs boundary summary
    bulk_chim_fracs = [grid[f"A={A:.2f}_beta={b:.2f}"]["basin_fraction_chimera"]
                       for b in beta_bulk for A in A_grid]
    boundary_chim_fracs = [grid[f"A={A:.2f}_beta={b:.2f}"]["basin_fraction_chimera"]
                           for b in beta_boundary for A in A_grid]
    print(f"\nBulk chimera basin fractions ({len(bulk_chim_fracs)} cells, "
          f"β ∈ {beta_bulk}):")
    print(f"  min={min(bulk_chim_fracs):.2f}, mean={np.mean(bulk_chim_fracs):.2f}, "
          f"max={max(bulk_chim_fracs):.2f}")
    print(f"Boundary chimera basin fractions ({len(boundary_chim_fracs)} cells, "
          f"β ∈ {beta_boundary}):")
    print(f"  min={min(boundary_chim_fracs):.2f}, "
          f"mean={np.mean(boundary_chim_fracs):.2f}, "
          f"max={max(boundary_chim_fracs):.2f}")

    out_dir = Path("/home/claude/epc/analysis/outputs")
    out_dir.mkdir(parents=True, exist_ok=True)
    with open(out_dir / "p10_phase_boundary_multiseed.json", "w") as f:
        json.dump({
            "sprint": 27,
            "phase": "1n",
            "purpose": ("multi-seed (A, beta) phase boundary refinement; "
                        "addresses Sprint 26 #31"),
            "config": {
                "N": N,
                "T": n_frames,
                "dt": 0.025,
                "record_dt": 1.0,
                "init_mode": "asymmetric_gaussian",
                "A_grid": A_grid.tolist(),
                "beta_boundary": beta_boundary,
                "beta_bulk": beta_bulk,
                "seeds": seeds,
                "tail_frames_for_classification": tail,
                "chim_band": [CHIM_LO, CHIM_HI],
                "sync_threshold": SYNC_THRESH,
            },
            "grid": grid,
            "cell_class_counts": cell_class_counts,
            "bulk_summary": {
                "min": min(bulk_chim_fracs),
                "mean": float(np.mean(bulk_chim_fracs)),
                "max": max(bulk_chim_fracs),
                "n_cells": len(bulk_chim_fracs),
            },
            "boundary_summary": {
                "min": min(boundary_chim_fracs),
                "mean": float(np.mean(boundary_chim_fracs)),
                "max": max(boundary_chim_fracs),
                "n_cells": len(boundary_chim_fracs),
            },
            "wall_time_s": elapsed,
        }, f, indent=2)
    print(f"\nWrote {out_dir / 'p10_phase_boundary_multiseed.json'}")


if __name__ == "__main__":
    main()
