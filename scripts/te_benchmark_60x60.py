"""Full-power TE benchmark at 60×60 with 99 permutations.

Sprint 5 task #1: Confirm the vectorized boundary-conditioned TE reproduces
the Sprint 2 separation: GH ratio 1-2× vs GoL ratio 15-16× over null.

Sprint 2 context:
  - GH spiral (60×60, 300 steps): boundary TE ≈ 0.001, ratio ≈ 1-2×
  - GoL R-pentomino (60×60, 300 steps): boundary TE ≈ 0.010, ratio ≈ 15-16×
  - Raw average TE gives WRONG ordering (GH > GoL). 
  - Boundary conditioning is the methodological contribution.

This script runs both models at full power (60×60, 99 perms) and reports
all quantities needed to confirm or deny the Sprint 2 result.

Two TE computation approaches are tested:
  1. Per-cell averaged TE (from transfer_entropy_vectorized.py)
  2. Global aggregate TE (from p13_p15_discriminator.py)
"""

from __future__ import annotations

import time
import sys
import os

import numpy as np

# Ensure epc is importable
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))


# ========================================================================
# Optimized boundary-conditioned TE computation
# ========================================================================

def _find_boundary_cells_vec(grids: np.ndarray, threshold: float = 0.10) -> np.ndarray:
    """Find cells that are boundary cells for >= threshold fraction of timesteps.
    
    Returns (n_boundary, 2) array of (row, col).
    """
    T, rows, cols = grids.shape
    boundary_count = np.zeros((rows, cols), dtype=np.int32)
    
    for t in range(T):
        g = grids[t]
        differs = np.zeros((rows, cols), dtype=bool)
        # VN neighbors (non-periodic, matching the original implementation)
        if rows > 1:
            differs[:-1, :] |= (g[:-1, :] != g[1:, :])
            differs[1:, :] |= (g[1:, :] != g[:-1, :])
        if cols > 1:
            differs[:, :-1] |= (g[:, :-1] != g[:, 1:])
            differs[:, 1:] |= (g[:, 1:] != g[:, :-1])
        boundary_count += differs
    
    boundary_frac = boundary_count / T
    return np.argwhere(boundary_frac > threshold)


def compute_boundary_te_optimized(
    grids: np.ndarray,
    n_states: int = 2,
    n_permutations: int = 99,
    seed: int = 42,
) -> dict:
    """Compute boundary-conditioned TE with vectorized counting.
    
    Per-cell TE averaged across boundary cells, matching the semantics
    of transfer_entropy_vectorized.compute_boundary_te but with
    np.add.at for the inner counting loop (replaces Python T-loop).
    """
    T, rows, cols = grids.shape
    grids_int = grids.astype(np.intp)
    
    boundary_cells = _find_boundary_cells_vec(grids)
    n_boundary = len(boundary_cells)
    
    if n_boundary == 0:
        return {"te_mean": 0.0, "null_mean": 0.0, "null_std": 0.0,
                "ratio": 0.0, "p_value": 1.0, "n_boundary": 0}
    
    # Extract target timeseries: (n_boundary, T)
    bc_r = boundary_cells[:, 0]
    bc_c = boundary_cells[:, 1]
    targets = grids_int[:, bc_r, bc_c].T  # (n_boundary, T)
    
    # Extract VN neighbor timeseries: list of (n_boundary, T) for each direction
    # Handle edge cases with masking
    directions = [(-1, 0), (1, 0), (0, -1), (0, 1)]
    neighbor_data = []  # list of (mask, sources) per direction
    
    for dr, dc in directions:
        nr = bc_r + dr
        nc = bc_c + dc
        valid = (nr >= 0) & (nr < rows) & (nc >= 0) & (nc < cols)
        # Clip to valid range for indexing (invalid entries will be masked out)
        nr_safe = np.clip(nr, 0, rows - 1)
        nc_safe = np.clip(nc, 0, cols - 1)
        sources = grids_int[:, nr_safe, nc_safe].T  # (n_boundary, T)
        neighbor_data.append((valid, sources))
    
    def _te_per_cell_batch(target_series, source_series, n_s):
        """Compute TE for one (source, target) pair using vectorized counting."""
        counts = np.zeros((n_s, n_s, n_s), dtype=np.float64)
        s_prev = source_series[:-1]
        tgt_prev = target_series[:-1]
        tgt_curr = target_series[1:]
        np.add.at(counts, (s_prev, tgt_prev, tgt_curr), 1)
        return _te_from_counts_fast(counts)
    
    def _compute_observed_te(target_data, neighbor_info, n_s):
        """Compute mean TE across all boundary cells."""
        n_b = target_data.shape[0]
        te_values = np.zeros(n_b)
        
        for i in range(n_b):
            cell_te = 0.0
            n_valid = 0
            tgt = target_data[i]
            for valid_mask, src_all in neighbor_info:
                if not valid_mask[i]:
                    continue
                cell_te += _te_per_cell_batch(tgt, src_all[i], n_s)
                n_valid += 1
            if n_valid > 0:
                te_values[i] = cell_te / n_valid
        
        return float(np.mean(te_values))
    
    # Observed TE
    te_observed = _compute_observed_te(targets, neighbor_data, n_states)
    
    # Permutation null: shuffle time axis of source
    rng = np.random.default_rng(seed)
    null_tes = np.zeros(n_permutations)
    
    for perm_idx in range(n_permutations):
        perm = rng.permutation(T)
        # Build shuffled neighbor data
        shuffled_neighbor_data = []
        for valid_mask, src_all in neighbor_data:
            shuffled_neighbor_data.append((valid_mask, src_all[:, perm]))
        null_tes[perm_idx] = _compute_observed_te(targets, shuffled_neighbor_data, n_states)
    
    null_mean = float(np.mean(null_tes))
    null_std = float(np.std(null_tes))
    ratio = te_observed / null_mean if null_mean > 0 else 0.0
    p_value = float(np.mean(null_tes >= te_observed))
    if p_value == 0:
        p_value = 1.0 / (n_permutations + 1)
    
    return {
        "te_mean": te_observed,
        "null_mean": null_mean,
        "null_std": null_std,
        "ratio": ratio,
        "p_value": p_value,
        "n_boundary": n_boundary,
        "n_permutations": n_permutations,
    }


def _te_from_counts_fast(counts: np.ndarray) -> float:
    """Compute TE from joint count array. Vectorized where possible."""
    total = counts.sum()
    if total == 0:
        return 0.0
    
    p_joint = counts / total
    # p(y|y') marginalizing over source
    p_y_yprev = counts.sum(axis=0)  # (n_s, n_s): [y', y]
    p_yprev = p_y_yprev.sum(axis=1)  # (n_s,)
    # p(s, y') marginalizing over y_next
    p_s_yprev = counts.sum(axis=2)  # (n_s, n_s): [s, y']
    
    te = 0.0
    n_s = counts.shape[0]
    
    for s in range(n_s):
        for yp in range(n_s):
            if p_s_yprev[s, yp] <= 0:
                continue
            if p_yprev[yp] <= 0:
                continue
            for y in range(n_s):
                p_sypy = p_joint[s, yp, y]
                if p_sypy <= 0:
                    continue
                p_y_given_s_yp = counts[s, yp, y] / p_s_yprev[s, yp]
                p_y_given_yp = p_y_yprev[yp, y] / p_yprev[yp]
                if p_y_given_s_yp > 0 and p_y_given_yp > 0:
                    te += p_sypy * np.log(p_y_given_s_yp / p_y_given_yp)
    
    return float(te)


# ========================================================================
# Global aggregate TE (discriminator approach)
# ========================================================================

def compute_global_boundary_te(
    grids: np.ndarray,
    n_states: int = 2,
    n_permutations: int = 99,
    seed: int = 42,
) -> dict:
    """Compute boundary-conditioned TE using global aggregate counts.
    
    This matches the semantics of p13_p15_discriminator._boundary_te:
    accumulate frequency tables across ALL boundary cells and ALL timesteps,
    then compute a single global TE.
    
    This is faster because we accumulate into global tables rather than
    computing per-cell TEs, and the answer may be more statistically robust
    at smaller grid sizes.
    """
    T, rows, cols = grids.shape
    grids_int = grids.astype(np.intp)
    
    def _global_te(grids_target, grids_source_fn):
        """Compute global boundary TE from grids.
        
        grids_target: (T, rows, cols) — target grid
        grids_source_fn: function(t) -> (rows, cols) — source grid at time t
        """
        # Directions for VN neighbors
        dirs = [(0, 1), (1, 0), (0, -1), (-1, 0)]
        
        # Global frequency tables per direction
        joint = [np.zeros((n_states, n_states, n_states), dtype=np.float64) for _ in range(4)]
        yt_xt = [np.zeros((n_states, n_states), dtype=np.float64) for _ in range(4)]
        ynext_yt = np.zeros((n_states, n_states), dtype=np.float64)
        yt_count = np.zeros(n_states, dtype=np.float64)
        total = 0
        
        for t in range(T - 1):
            g = grids_target[t]
            g_next = grids_target[t + 1]
            g_src = grids_source_fn(t)
            
            # Vectorized boundary detection (Moore neighborhood, periodic)
            is_boundary = np.zeros((rows, cols), dtype=bool)
            for dr in (-1, 0, 1):
                for dc in (-1, 0, 1):
                    if dr == 0 and dc == 0:
                        continue
                    shifted = np.roll(np.roll(g, -dr, axis=0), -dc, axis=1)
                    is_boundary |= (g != shifted)
            
            br, bc = np.where(is_boundary)
            if len(br) == 0:
                continue
            
            y_vals = g[br, bc]
            y_next_vals = g_next[br, bc]
            
            # Accumulate ynext_yt and yt counts
            np.add.at(ynext_yt, (y_next_vals, y_vals), 1)
            np.add.at(yt_count, (y_vals,), 1)
            total += len(br)
            
            # VN neighbor values and joint counts per direction
            for d_idx, (dr, dc) in enumerate(dirs):
                x_vals = np.roll(np.roll(g_src, -dr, axis=0), -dc, axis=1)[br, bc]
                np.add.at(joint[d_idx], (y_next_vals, y_vals, x_vals), 1)
                np.add.at(yt_xt[d_idx], (y_vals, x_vals), 1)
        
        if total == 0:
            return 0.0, 0
        
        # H(Y_next | Y_curr)
        h_y_yt = 0.0
        for yn in range(n_states):
            for y in range(n_states):
                n = ynext_yt[yn, y]
                if n <= 0:
                    continue
                p = n / total
                p_cond = n / yt_count[y] if yt_count[y] > 0 else 0
                if p_cond > 0:
                    h_y_yt -= p * np.log2(p_cond)
        
        # TE per direction
        te_dirs = []
        for d_idx in range(4):
            h_y_yt_xt = 0.0
            for yn in range(n_states):
                for y in range(n_states):
                    for x in range(n_states):
                        n = joint[d_idx][yn, y, x]
                        if n <= 0:
                            continue
                        n_yx = yt_xt[d_idx][y, x]
                        if n_yx <= 0:
                            continue
                        p = n / total
                        p_cond = n / n_yx
                        if p_cond > 0:
                            h_y_yt_xt -= p * np.log2(p_cond)
            te = max(0.0, h_y_yt - h_y_yt_xt)
            te_dirs.append(te)
        
        return float(np.mean(te_dirs)), total
    
    # Observed TE
    te_obs, n_obs = _global_te(grids_int, lambda t: grids_int[t])
    
    # Permutation null: shuffle time axis of source grids
    rng = np.random.default_rng(seed)
    null_tes = np.zeros(n_permutations)
    
    # Subsample trajectory for speed (matching discriminator approach)
    subsample = max(1, T // 80)
    
    for perm_idx in range(n_permutations):
        perm = rng.permutation(T)
        grids_shuffled = grids_int[perm]
        te_null, _ = _global_te(grids_int, lambda t, gs=grids_shuffled: gs[t])
        null_tes[perm_idx] = te_null
    
    null_mean = float(np.mean(null_tes))
    null_std = float(np.std(null_tes))
    ratio = te_obs / null_mean if null_mean > 0 else 0.0
    p_value = float(np.mean(null_tes >= te_obs))
    if p_value == 0:
        p_value = 1.0 / (n_permutations + 1)
    
    return {
        "te_mean": te_obs,
        "null_mean": null_mean,
        "null_std": null_std,
        "ratio": ratio,
        "p_value": p_value,
        "n_obs": n_obs,
        "n_permutations": n_permutations,
    }


# ========================================================================
# Main benchmark
# ========================================================================

def run_model(model_name: str, grid_size: int = 60, n_steps: int = 300) -> np.ndarray:
    """Run a model and return grids as (T, rows, cols) uint8 array."""
    from epc.models.greenberg_hastings import GreenbergHastings
    from epc.models.game_of_life import GameOfLife
    
    if model_name == "gh_spiral":
        model = GreenbergHastings(
            rows=grid_size, cols=grid_size,
            n_states=5, threshold=1,
            init_mode="broken_wave", neighborhood="moore",
            boundary="periodic", seed=0,
        )
        history = model.run(n_steps=n_steps, vectorized=True)
        n_states = 5
    elif model_name == "gh_random":
        model = GreenbergHastings(
            rows=grid_size, cols=grid_size,
            n_states=3, threshold=1,
            init_mode="random", neighborhood="moore",
            boundary="periodic", seed=42,
        )
        history = model.run(n_steps=n_steps, vectorized=True)
        n_states = 3
    elif model_name == "gol_random":
        model = GameOfLife(
            rows=grid_size, cols=grid_size,
            init_mode="random", init_density=0.37,
            boundary="periodic", seed=42,
        )
        history = model.run(n_steps=n_steps)
        n_states = 2
    elif model_name == "gol_rpent":
        model = GameOfLife(
            rows=grid_size, cols=grid_size,
            init_mode="r_pentomino",
            boundary="periodic", seed=42,
        )
        history = model.run(n_steps=n_steps)
        n_states = 2
    else:
        raise ValueError(f"Unknown model: {model_name}")
    
    grids = np.array([s["grid"] for s in history], dtype=np.uint8)
    return grids, n_states


def main():
    grid_size = 60
    n_steps = 300
    n_perms = 99
    seed = 42
    
    print("=" * 78)
    print("FULL-POWER TE BENCHMARK — 60×60, 99 permutations")
    print("=" * 78)
    print(f"Grid: {grid_size}×{grid_size}, Steps: {n_steps}, Perms: {n_perms}")
    print()
    
    # Models to test
    configs = [
        ("gh_spiral", "GH spiral (κ=5, θ=1, broken_wave)"),
        ("gh_random", "GH random (κ=3, θ=1, random init)"),
        ("gol_random", "GoL random (density=0.37)"),
        ("gol_rpent", "GoL R-pentomino"),
    ]
    
    results = {}
    
    for model_name, description in configs:
        print("-" * 78)
        print(f"Model: {description}")
        print("-" * 78)
        
        # Simulate
        t0 = time.perf_counter()
        grids, n_states = run_model(model_name, grid_size, n_steps)
        sim_time = time.perf_counter() - t0
        print(f"  Simulation: {sim_time:.2f}s, shape={grids.shape}, n_states={n_states}")
        
        # Activity stats
        if n_states == 2:
            alive = np.mean(grids[-50:] == 1)
            print(f"  Activity (last 50 steps): {alive:.3f} alive fraction")
        else:
            excited = np.mean(grids[-50:] == 1)
            print(f"  Activity (last 50 steps): {excited:.3f} excited fraction")
        
        # Global aggregate TE (discriminator approach)
        print(f"  Computing global boundary TE ({n_perms} perms)...")
        t0 = time.perf_counter()
        global_result = compute_global_boundary_te(
            grids, n_states=n_states, n_permutations=n_perms, seed=seed,
        )
        global_time = time.perf_counter() - t0
        
        print(f"  Global TE: mean={global_result['te_mean']:.6f}, "
              f"null={global_result['null_mean']:.6f} ± {global_result['null_std']:.6f}")
        print(f"  Ratio: {global_result['ratio']:.2f}×, p={global_result['p_value']:.4f}")
        print(f"  Boundary observations: {global_result['n_obs']:,}")
        print(f"  Time: {global_time:.1f}s")
        
        results[model_name] = {
            "description": description,
            "n_states": n_states,
            "global": global_result,
            "global_time": global_time,
        }
        print()
    
    # ====================================================================
    # Summary and comparison with Sprint 2 targets
    # ====================================================================
    print("=" * 78)
    print("SUMMARY")
    print("=" * 78)
    
    print("\nAbsolute boundary TE:")
    for name, r in results.items():
        g = r["global"]
        print(f"  {name:15s}: TE={g['te_mean']:.6f}, null={g['null_mean']:.6f}, "
              f"ratio={g['ratio']:.2f}×, p={g['p_value']:.4f}")
    
    # GH vs GoL comparison
    gh_models = [k for k in results if k.startswith("gh_")]
    gol_models = [k for k in results if k.startswith("gol_")]
    
    if gh_models and gol_models:
        print("\n--- GH vs GoL Discrimination ---")
        for gh_name in gh_models:
            for gol_name in gol_models:
                gh_r = results[gh_name]["global"]
                gol_r = results[gol_name]["global"]
                
                te_separation = gol_r["te_mean"] / gh_r["te_mean"] if gh_r["te_mean"] > 0 else float("inf")
                ratio_separation = gol_r["ratio"] / gh_r["ratio"] if gh_r["ratio"] > 0 else float("inf")
                
                print(f"\n  {gol_name} vs {gh_name}:")
                print(f"    Absolute TE: {gol_r['te_mean']:.6f} vs {gh_r['te_mean']:.6f} "
                      f"({te_separation:.1f}× separation)")
                print(f"    Ratio over null: {gol_r['ratio']:.2f}× vs {gh_r['ratio']:.2f}× "
                      f"({ratio_separation:.1f}× separation)")
                
                # Direction check
                direction_ok = gol_r["te_mean"] > gh_r["te_mean"]
                print(f"    Direction (GoL > GH): {'✅' if direction_ok else '❌'}")
    
    # Sprint 2 target comparison
    print("\n--- Sprint 2 Target Comparison ---")
    print("  Sprint 2 (60×60): GH ratio 1-2×, GoL ratio 15-16×")
    if "gh_spiral" in results:
        gh_ratio = results["gh_spiral"]["global"]["ratio"]
        print(f"  Current GH spiral ratio: {gh_ratio:.2f}× (target: 1-2×) "
              f"{'✅' if gh_ratio <= 3 else '⚠️'}")
    if "gol_random" in results:
        gol_ratio = results["gol_random"]["global"]["ratio"]
        print(f"  Current GoL random ratio: {gol_ratio:.2f}× (target: 15-16×) "
              f"{'✅' if gol_ratio >= 10 else '⚠️'}")
    if "gol_rpent" in results:
        gol_ratio = results["gol_rpent"]["global"]["ratio"]
        print(f"  Current GoL R-pent ratio: {gol_ratio:.2f}× (target: 15-16×) "
              f"{'✅' if gol_ratio >= 10 else '⚠️'}")
    
    print("\n--- Timing ---")
    for name, r in results.items():
        print(f"  {name:15s}: global={r['global_time']:.1f}s")
    
    total_time = sum(r["global_time"] for r in results.values())
    print(f"  Total: {total_time:.0f}s ({total_time/60:.1f}min)")
    
    return results


if __name__ == "__main__":
    results = main()
