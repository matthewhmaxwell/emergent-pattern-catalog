"""Sprint 26 Phase 1k (chunked) — runs ONE N at a time so we fit in a single
bash call. Pass N as command-line arg. Output appends to a per-N JSON.
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
TAU_CONSECUTIVE = 10
INIT_BURNIN_FRAMES = 20


def chimera_alive(r_traj):
    return (r_traj > CHIM_LO) & (r_traj < CHIM_HI)


def measure_lifetime(r_traj, tau=TAU_CONSECUTIVE):
    alive = chimera_alive(r_traj)
    n = len(alive)
    for t in range(n - tau + 1):
        if not alive[t : t + tau].any():
            return t, False
    return n, True


def main():
    N = int(sys.argv[1])
    n_seeds = int(sys.argv[2]) if len(sys.argv) > 2 else 30
    T_max = int(sys.argv[3]) if len(sys.argv) > 3 else 100  # default 100 frames
    seeds = list(range(n_seeds))
    A = 0.995
    beta = 0.18

    results = []
    counts = {"chimera_seeds": 0, "sync_seeds": 0,
              "censored": 0, "completed": 0, "lifetimes": []}

    t_start = time.time()
    print(f"N={N}, beta={beta}, A={A}, {n_seeds} seeds, T_max={T_max}",
          flush=True)
    for seed in seeds:
        t0 = time.time()
        m = KuramotoNonlocal(N=N, A=A, beta=beta, seed=seed,
                             init_mode="asymmetric_gaussian")
        hist = m.run(n_frames=T_max)
        r_traj = np.array([h["r"] for h in hist])
        r_burnin = float(r_traj[:INIT_BURNIN_FRAMES].mean())
        if r_burnin > 0.95:
            counts["sync_seeds"] += 1
            results.append({"seed": seed, "basin": "sync",
                            "r_burnin": r_burnin, "lifetime": 0,
                            "censored": False})
            print(f"  seed={seed:3d}: SYNC (r_burnin={r_burnin:.3f}) "
                  f"[{time.time()-t0:.1f}s]", flush=True)
            continue

        counts["chimera_seeds"] += 1
        lifetime, censored = measure_lifetime(r_traj)
        counts["lifetimes"].append(lifetime)
        if censored:
            counts["censored"] += 1
        else:
            counts["completed"] += 1
        results.append({
            "seed": seed,
            "basin": "chimera",
            "r_burnin": r_burnin,
            "r_traj": r_traj.tolist(),
            "lifetime": int(lifetime),
            "censored": bool(censored),
        })
        print(f"  seed={seed:3d}: CHIMERA r_burnin={r_burnin:.3f} "
              f"lifetime={lifetime:3d}{' (censored)' if censored else ''} "
              f"[{time.time()-t0:.1f}s]", flush=True)

    elapsed = time.time() - t_start
    chim_n = counts["chimera_seeds"]
    if chim_n > 0:
        lt = np.array(counts["lifetimes"])
        median_lt = float(np.median(lt))
        mean_lt = float(np.mean(lt))
        cens_frac = counts["censored"] / chim_n
    else:
        median_lt = mean_lt = cens_frac = float("nan")

    summary = {
        "N": N,
        "T_max": T_max,
        "n_seeds": n_seeds,
        "chimera_basin_count": chim_n,
        "sync_basin_count": counts["sync_seeds"],
        "chimera_basin_fraction": chim_n / n_seeds,
        "censored_count": counts["censored"],
        "completed_count": counts["completed"],
        "censored_fraction": cens_frac,
        "median_lifetime_frames": median_lt,
        "mean_lifetime_frames": mean_lt,
        "wall_time_s": elapsed,
    }
    print(f"\nN={N} summary: {chim_n}/{n_seeds} chimera basin, "
          f"{counts['censored']} censored at T_max={T_max}, "
          f"median lifetime={median_lt:.1f} ({elapsed:.1f}s wall)",
          flush=True)

    out_dir = Path("/home/claude/epc/analysis/outputs")
    out_dir.mkdir(parents=True, exist_ok=True)
    with open(out_dir / f"p10_lifetime_N{N}.json", "w") as f:
        json.dump({"summary": summary,
                   "results": [{k: v for k, v in r.items() if k != "r_traj"}
                               for r in results]},
                  f, indent=2)
    np_arrays = {f"seed{r['seed']}": np.array(r["r_traj"])
                 for r in results if "r_traj" in r}
    if np_arrays:
        np.savez(out_dir / f"p10_lifetime_N{N}_trajectories.npz", **np_arrays)
    print(f"Wrote {out_dir / f'p10_lifetime_N{N}.json'}")


if __name__ == "__main__":
    main()
