"""Full-power TE benchmark at 60×60 with 99 permutations.

Sprint 5 task #1: Confirm boundary-conditioned TE reproduces the Sprint 2
GH 1-2× vs GoL 15-16× separation at full scale.

Sprint 2 reference (60×60, 300 steps):
  GH spiral:       ratio vs GH control = 1.0×  → P13
  GH random:       ratio vs GH control = 2.1×  → P13
  GoL random:      ratio vs GH control = 15.1× → P15_candidate
  GoL R-pentomino: ratio vs GH control = 16.1× → P15_candidate

Tests:
  1. Numerical agreement: global TE == discriminator._boundary_te (fast)
  2. Full 60×60 benchmark: direction, ratios, classification, significance (slow)
"""

from __future__ import annotations

import time
import sys
import os

import numpy as np
import pytest

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))


def _run_gh(rows, cols, n_states, n_steps, init_mode="broken_wave", seed=0):
    from epc.models.greenberg_hastings import GreenbergHastings
    m = GreenbergHastings(rows=rows, cols=cols, n_states=n_states, threshold=1,
                          init_mode=init_mode, neighborhood="moore",
                          boundary="periodic", seed=seed)
    h = m.run(n_steps=n_steps, vectorized=True)
    return np.array([s["grid"] for s in h], dtype=np.uint8)


def _run_gol(rows, cols, n_steps, init_mode="random", seed=42):
    from epc.models.game_of_life import GameOfLife
    kw = dict(rows=rows, cols=cols, init_mode=init_mode,
              boundary="periodic", seed=seed)
    if init_mode == "random":
        kw["init_density"] = 0.37
    m = GameOfLife(**kw)
    h = m.run(n_steps=n_steps)
    return np.array([s["grid"] for s in h], dtype=np.uint8)


def test_global_te_matches_discriminator():
    """Global TE must match discriminator._boundary_te exactly."""
    from epc.metrics.transfer_entropy_global import _global_te_pass
    from epc.detectors.p13_p15_discriminator import P13P15Discriminator

    grids = _run_gh(25, 25, n_states=3, n_steps=100, init_mode="random", seed=42)

    disc = P13P15Discriminator()
    disc_history = [{"grid": g, "grid_dims": (25, 25)} for g in grids]
    disc_te, _, _ = disc._boundary_te(disc_history)

    global_te, n_obs = _global_te_pass(grids.astype(np.intp), grids.astype(np.intp), 3)

    print(f"\n  Discriminator: {disc_te:.8f}")
    print(f"  Global opt:    {global_te:.8f}")
    assert abs(global_te - disc_te) < 1e-6, \
        f"Mismatch: {global_te:.8f} vs {disc_te:.8f}"
    print("  ✅ Exact agreement")


@pytest.mark.slow
def test_te_benchmark_60x60():
    """Full-power benchmark confirming Sprint 2 GH/GoL separation."""
    from epc.metrics.transfer_entropy_global import compute_global_boundary_te

    G = 60
    T = 300
    P = 99

    print(f"\n{'='*70}")
    print(f"FULL-POWER TE BENCHMARK — {G}×{G}, {P} perms")
    print(f"{'='*70}")

    t_total = time.perf_counter()

    configs = [
        ("gh_spiral",  _run_gh(G, G, 5, T, "broken_wave", 0),  5),
        ("gh_random",  _run_gh(G, G, 3, T, "random", 42),      3),
        ("gol_random", _run_gol(G, G, T, "random", 42),        2),
        ("gol_rpent",  _run_gol(G, G, T, "r_pentomino", 42),   2),
    ]

    results = {}
    for name, grids, ns in configs:
        print(f"\n  {name}: computing ({ns} states)...")
        t0 = time.perf_counter()
        r = compute_global_boundary_te(grids, n_states=ns, n_permutations=P, seed=42)
        dt = time.perf_counter() - t0
        results[name] = r
        print(f"    TE={r.te_mean:.6f}  null={r.null_mean:.6f}±{r.null_std:.6f}  "
              f"p={r.p_value:.4f}  obs={r.n_boundary_obs:,}  {dt:.0f}s")

    elapsed = time.perf_counter() - t_total

    # ---- GH-control ratios (the Sprint 2 metric) ----
    gh_te = results["gh_spiral"].te_mean
    ratios = {name: r.te_mean / gh_te for name, r in results.items()}

    print(f"\n  {'Model':15s} {'TE':>10s} {'vs GH':>8s} {'Sprint 2':>10s}")
    print(f"  {'-'*48}")
    s2 = {"gh_spiral": "1.0×", "gh_random": "2.1×",
          "gol_random": "15.1×", "gol_rpent": "16.1×"}
    for name, _, _ in configs:
        print(f"  {name:15s} {results[name].te_mean:10.6f} "
              f"{ratios[name]:7.1f}× {s2[name]:>10s}")

    print(f"\n  Total: {elapsed:.0f}s ({elapsed/60:.1f}min)")

    # ---- Assertions ----

    # A. Direction: GoL TE >> GH TE
    for gol in ["gol_random", "gol_rpent"]:
        assert results[gol].te_mean > gh_te * 5, \
            f"{gol} TE not >5× GH ({results[gol].te_mean:.6f} vs {gh_te:.6f})"

    # B. GoL/GH ratio in Sprint 2 range [10, 25]
    for gol in ["gol_random", "gol_rpent"]:
        assert 10 <= ratios[gol] <= 25, \
            f"{gol} ratio {ratios[gol]:.1f}× outside [10, 25]"

    # C. GH ratios both ≤ 3 (P13 classification)
    for gh in ["gh_spiral", "gh_random"]:
        assert ratios[gh] <= 3.0, \
            f"{gh} ratio {ratios[gh]:.1f}× should be ≤3"

    # D. Statistical significance for GoL
    for gol in ["gol_random", "gol_rpent"]:
        assert results[gol].p_value <= 0.01, \
            f"{gol} p={results[gol].p_value:.4f} > 0.01"

    # E. Classification: threshold at ratio > 3
    for name in results:
        r = results[name]
        cls = "P15_candidate" if (ratios[name] > 3.0 and r.p_value <= 0.05) else "P13"
        expected = "P15_candidate" if name.startswith("gol") else "P13"
        assert cls == expected, f"{name}: {cls} != {expected}"

    print(f"\n  ✅ Sprint 2 separation CONFIRMED at 60×60.")
    print(f"     GoL/GH = {ratios['gol_random']:.1f}× and {ratios['gol_rpent']:.1f}× "
          f"(target: 15-16×)")


if __name__ == "__main__":
    print("TEST 1: Numerical agreement")
    test_global_te_matches_discriminator()
    print("\n\nTEST 2: Full-power 60×60 benchmark")
    test_te_benchmark_60x60()
    print("\n\nALL TESTS PASSED")
