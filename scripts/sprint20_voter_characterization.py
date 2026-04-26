"""Sprint 20 Phase 1: voter model characterization.

Goal: empirically measure, before writing any detector:
  1. Domain-wall density decay ρ_w(t) — expected scaling in 2D voter model
  2. Magnetization |m(t)| trajectory
  3. Moran's I trajectory
  4. Consensus time τ_c(L) — expected ∝ L² ln L
  5. Whether the ln(t)/t vs t^(-1/2) distinction is visible at our scales

Output: /tmp/voter_characterization.json with per-seed time series + scaling
fits, plus console summary.
"""
from __future__ import annotations

import json
import time
from typing import Any

import numpy as np

from epc.models.voter import VoterModel


def characterize_one(L: int, seed: int, n_sweeps: int,
                     record_every: int) -> dict[str, Any]:
    """Run voter model and record time series.

    Returns dict with arrays: t, abs_m, wall, moran, consensus_t (or None).
    """
    m = VoterModel(rows=L, cols=L, seed=seed, init_mode="random")
    m.setup()
    ts: list[int] = [0]
    abs_m: list[float] = [abs(m._state_dict()["magnetization"])]
    wall: list[float] = [m._state_dict()["wall_density"]]
    moran: list[float] = [m._state_dict()["moran_i"]]
    consensus_t: int | None = None

    for t in range(1, n_sweeps + 1):
        s = m.step()
        if t % record_every == 0 or s["consensus_reached"]:
            ts.append(t)
            abs_m.append(s["abs_magnetization"])
            wall.append(s["wall_density"])
            moran.append(s["moran_i"])
        if s["consensus_reached"] and consensus_t is None:
            consensus_t = t
            break

    return {
        "L": L,
        "seed": seed,
        "ts": ts,
        "abs_m": abs_m,
        "wall": wall,
        "moran": moran,
        "consensus_t": consensus_t,
        "n_sweeps_run": ts[-1],
    }


def fit_power_law(ts: np.ndarray, ys: np.ndarray,
                  t_min: int, t_max: int) -> tuple[float, float]:
    """Fit log(y) = a + b*log(t) over [t_min, t_max]. Returns (exponent, intercept)."""
    mask = (ts >= t_min) & (ts <= t_max) & (ys > 0)
    if mask.sum() < 3:
        return (np.nan, np.nan)
    lt = np.log(ts[mask])
    ly = np.log(ys[mask])
    b, a = np.polyfit(lt, ly, 1)
    return (float(b), float(a))


def fit_log_over_t(ts: np.ndarray, ys: np.ndarray,
                   t_min: int, t_max: int) -> float:
    """Fit y = c * ln(t) / t over [t_min, t_max]. Returns c."""
    mask = (ts >= t_min) & (ts <= t_max) & (ys > 0) & (ts > 1)
    if mask.sum() < 3:
        return float("nan")
    tv = ts[mask]
    yv = ys[mask]
    # y = c * ln(t)/t  ⇒ c ≈ <y * t / ln(t)>
    c = float(np.mean(yv * tv / np.log(tv)))
    return c


def main() -> None:
    # Config: L=32 short runs to pin consensus time; L=64, 128 longer for
    # coarsening scaling. Per Sprint 19 finite-size rule, we'll run full
    # {64, 128} plus short {32} for consensus.
    t0 = time.time()
    all_runs: list[dict[str, Any]] = []

    # L=32: consensus time characterization. Should reach consensus typically.
    # Expected τ_c ~ 32² ln 32 ≈ 3550 sweeps.
    print("=== L=32: consensus time (5 seeds, up to 15000 sweeps) ===")
    for seed in (0, 1, 42, 123, 2024):
        run = characterize_one(L=32, seed=seed, n_sweeps=15_000, record_every=25)
        all_runs.append(run)
        tc = run["consensus_t"]
        print(f"  seed={seed}: consensus={'t=' + str(tc) if tc else 'NOT REACHED'}, "
              f"n_recorded={len(run['ts'])}")

    # L=64: coarsening dynamics. 2000 sweeps covers the coarsening regime.
    # Expected τ_c ~ 64² ln 64 ≈ 17000 sweeps.
    print("\n=== L=64: coarsening (5 seeds, 2000 sweeps each) ===")
    for seed in (0, 1, 42, 123, 2024):
        run = characterize_one(L=64, seed=seed, n_sweeps=2000, record_every=10)
        all_runs.append(run)
        tc = run["consensus_t"]
        print(f"  seed={seed}: final |m|={run['abs_m'][-1]:.3f}, "
              f"final wall={run['wall'][-1]:.4f}, "
              f"consensus={'t=' + str(tc) if tc else 'not yet'}")

    # L=128: coarsening at larger scale. 1000 sweeps (far below τ_c).
    print("\n=== L=128: coarsening (3 seeds, 1000 sweeps each) ===")
    for seed in (0, 42, 2024):
        run = characterize_one(L=128, seed=seed, n_sweeps=1000, record_every=10)
        all_runs.append(run)
        print(f"  seed={seed}: final |m|={run['abs_m'][-1]:.3f}, "
              f"final wall={run['wall'][-1]:.4f}")

    # Fits: for L=64 and L=128, measure coarsening scaling
    print("\n=== Wall-density decay fits (fit window: t ∈ [30, 500]) ===")
    for L in (64, 128):
        runs_L = [r for r in all_runs if r["L"] == L]
        exps = []
        for r in runs_L:
            ts = np.array(r["ts"])
            ws = np.array(r["wall"])
            exp, _ = fit_power_law(ts, ws, 30, 500)
            exps.append(exp)
            print(f"  L={L}, seed={r['seed']}: power-law exponent b={exp:+.3f}")
        print(f"  L={L} mean exponent: {np.mean(exps):+.3f} "
              f"(expected naive −0.5; true −1 with ln(t) correction)")

    # τ_c(L) from L=32 runs
    print("\n=== Consensus time (L=32 only reaches consensus in budget) ===")
    tcs_32 = [r["consensus_t"] for r in all_runs
              if r["L"] == 32 and r["consensus_t"] is not None]
    if tcs_32:
        print(f"  L=32: n_consense={len(tcs_32)}, "
              f"τ_c mean={np.mean(tcs_32):.0f}, median={np.median(tcs_32):.0f}")
        print(f"  Expected L² ln L = {32**2 * np.log(32):.0f}")
    else:
        print("  L=32: no seed reached consensus in 15000 sweeps.")

    # Save
    out = {
        "runs": all_runs,
        "elapsed_s": time.time() - t0,
        "config": {
            "L32_seeds": [0, 1, 42, 123, 2024],
            "L32_n_sweeps": 15_000,
            "L64_seeds": [0, 1, 42, 123, 2024],
            "L64_n_sweeps": 2000,
            "L128_seeds": [0, 42, 2024],
            "L128_n_sweeps": 1000,
        },
    }
    with open("/tmp/voter_characterization.json", "w") as f:
        json.dump(out, f, indent=2)
    print(f"\n=== Total elapsed: {time.time() - t0:.1f}s ===")
    print("Saved: /tmp/voter_characterization.json")


if __name__ == "__main__":
    main()
