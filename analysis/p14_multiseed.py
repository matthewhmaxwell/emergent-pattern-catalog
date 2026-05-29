"""
P14 (BTW sandpile) multi-seed extension — Sprint 55 dim2 closure.

Runs N_SEEDS independent seeds of the canonical BTW critical sandpile,
fits the avalanche-size power-law exponent τ for each, and reports
per-seed values plus aggregate (mean, std, CV).

Parameters: L=32, n_drive=30_000, n_burn=3_000. L=32 is in the same 2D BTW
universality class as the canonical L=64 run (τ ≈ 1.20–1.25; finite-size
cutoff at s_max ~ L² = 1024). n_drive=30_000 gives ~12k non-zero avalanches
per seed (~28s/seed, ~560s for 20 seeds) — sufficient for stable MLE τ
estimates with <1% statistical uncertainty. n_burn=3_000 ensures the system
reaches the critical state before the measurement window.

Primary observable: τ (MLE power-law exponent, xmin=1, Clauset et al. 2009).

Output: analysis/outputs/p14_multiseed.json
"""

from __future__ import annotations

import json
import time
from typing import Any, Dict, List

import numpy as np

from epc.models.btw_sandpile import BTWSandpileParams, run_sandpile
from epc.detectors.p14_soc import fit_power_law


# Parameters — L=32 (same 2D universality class as canonical L=64), n_drive=30_000
# gives ~12k non-zero avalanches per seed (reliable τ estimates) at ~28s/seed
# (~560s total for 20 seeds). n_burn=3_000 ensures critical state before
# measurement window.
L = 32
N_DRIVE = 30_000
N_BURN = 3_000
N_SEEDS = 20

# Seed sequence: offset from 0 to avoid overlap with phase2a panel seeds (0–4)
SEED_BASE = 100


def run_seed(seed: int) -> Dict[str, Any]:
    """Run one BTW sandpile seed and return τ and supporting statistics."""
    params = BTWSandpileParams(
        L=L,
        n_drive=N_DRIVE,
        n_burn=N_BURN,
        seed=seed,
    )
    result = run_sandpile(params)
    sizes = result.avalanche_sizes[result.avalanche_sizes > 0]
    fit = fit_power_law(sizes, xmin=1)
    return {
        "seed": seed,
        "tau": float(fit.tau),
        "tau_logbin": float(fit.tau_logbin),
        "n_nonzero_avalanches": int(len(sizes)),
        "max_avalanche_size": int(sizes.max()),
        "lr_vs_exponential": float(fit.lr_vs_exponential),
        "lr_vs_exponential_p": float(fit.lr_vs_exponential_p),
    }


def main() -> None:
    t0 = time.time()
    print(f"P14 BTW sandpile multi-seed: L={L}, n_drive={N_DRIVE}, "
          f"n_burn={N_BURN}, N_seeds={N_SEEDS}")

    per_seed: List[Dict[str, Any]] = []
    for i in range(N_SEEDS):
        seed = SEED_BASE + i
        t_seed = time.time()
        result = run_seed(seed)
        elapsed = time.time() - t_seed
        per_seed.append(result)
        print(f"  seed {seed:3d}: τ={result['tau']:.4f}  "
              f"n_avs={result['n_nonzero_avalanches']:6d}  "
              f"({elapsed:.1f}s)")

    taus = np.array([r["tau"] for r in per_seed])
    tau_mean = float(np.mean(taus))
    tau_std = float(np.std(taus, ddof=1))
    tau_cv = float(tau_std / tau_mean)

    print(f"\nτ aggregate ({N_SEEDS} seeds):")
    print(f"  mean = {tau_mean:.4f}")
    print(f"  std  = {tau_std:.4f}")
    print(f"  CV   = {tau_cv:.4f}")

    output: Dict[str, Any] = {
        "sprint": 55,
        "description": "P14 BTW sandpile multi-seed τ extension",
        "params": {
            "L": L,
            "n_drive": N_DRIVE,
            "n_burn": N_BURN,
            "n_seeds": N_SEEDS,
            "seed_base": SEED_BASE,
        },
        "primary_observable": "tau (MLE power-law exponent, xmin=1)",
        "per_seed": per_seed,
        "aggregate": {
            "tau_mean": tau_mean,
            "tau_std": tau_std,
            "tau_cv": tau_cv,
            "tau_min": float(taus.min()),
            "tau_max": float(taus.max()),
        },
        "elapsed_seconds": time.time() - t0,
    }

    out_path = "analysis/outputs/p14_multiseed.json"
    with open(out_path, "w") as f:
        json.dump(output, f, indent=2)
    print(f"\nSaved → {out_path}  (elapsed {output['elapsed_seconds']:.1f}s)")


if __name__ == "__main__":
    main()
