"""
P21 (Hegselmann-Krause opinion dynamics) multi-seed extension — Sprint 60 dim2 closure.

Runs N_SEEDS independent seeds at the canonical ε=0.20 (clear fragmentation regime,
well below the ε_c ≈ 0.24–0.27 consensus transition) and reports per-seed converged
cluster count plus aggregate (mean, std, CV).

Primary observable: converged cluster count (histogram bin-merging, merge_gap=0.05).

Parameters match the existing reproduction script (p21_hegselmann2002.py):
  N=100 agents, uniform [0,1] IC, deterministic synchronous HK update,
  convergence_tol=1e-6, T_max=10000.

Output: analysis/outputs/p21_multiseed.json
"""

from __future__ import annotations

import json
import os
import sys
import time
from typing import Any, Dict, List

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), ".."))
from epc.models.hegselmann_krause import HegselmannKrauseModel

# Parameters
N_AGENTS = 100
EPSILON = 0.20  # canonical fragmentation regime (well below ε_c)
T_MAX = 10000
MERGE_GAP = 0.05

N_SEEDS = 20
SEED_BASE = 100


def count_clusters(opinions: np.ndarray, merge_gap: float = MERGE_GAP) -> int:
    """Count distinct opinion clusters with min-gap merging."""
    sorted_op = np.sort(opinions)
    clusters = 1
    for i in range(1, len(sorted_op)):
        if sorted_op[i] - sorted_op[i - 1] > merge_gap:
            clusters += 1
    return clusters


def run_seed(seed: int) -> Dict[str, Any]:
    """Run one HK seed and return cluster count."""
    model = HegselmannKrauseModel(
        n_agents=N_AGENTS,
        epsilon=EPSILON,
        init_mode="uniform",
        convergence_tol=1e-6,
        seed=seed,
    )
    history = model.run(n_steps=T_MAX)
    final_opinions = history[-1]["opinions"]
    n_clusters = count_clusters(final_opinions)
    n_steps_actual = len(history) - 1  # subtract initial state

    return {
        "seed": seed,
        "cluster_count": n_clusters,
        "n_steps_converged": n_steps_actual,
    }


def main() -> None:
    t0 = time.time()
    print(f"P21 HK opinion multi-seed: N={N_AGENTS}, ε={EPSILON}, "
          f"N_seeds={N_SEEDS}")

    per_seed: List[Dict[str, Any]] = []
    for i in range(N_SEEDS):
        seed = SEED_BASE + i
        t_seed = time.time()
        result = run_seed(seed)
        elapsed = time.time() - t_seed
        per_seed.append(result)
        print(f"  seed {seed:3d}: clusters={result['cluster_count']}  "
              f"steps={result['n_steps_converged']}  ({elapsed:.1f}s)")

    counts = np.array([r["cluster_count"] for r in per_seed])
    count_mean = float(np.mean(counts))
    count_std = float(np.std(counts, ddof=1))
    count_cv = float(count_std / count_mean) if count_mean > 0 else float("nan")

    print(f"\ncluster_count aggregate ({N_SEEDS} seeds):")
    print(f"  mean = {count_mean:.2f}")
    print(f"  std  = {count_std:.2f}")
    print(f"  CV   = {count_cv:.4f}")

    output: Dict[str, Any] = {
        "sprint": 60,
        "description": "P21 HK opinion multi-seed cluster_count extension",
        "params": {
            "n_agents": N_AGENTS,
            "epsilon": EPSILON,
            "t_max": T_MAX,
            "merge_gap": MERGE_GAP,
            "n_seeds": N_SEEDS,
            "seed_base": SEED_BASE,
        },
        "primary_observable": "converged cluster count (histogram bin-merge, gap=0.05)",
        "per_seed": per_seed,
        "aggregate": {
            "cluster_count_mean": count_mean,
            "cluster_count_std": count_std,
            "cluster_count_cv": count_cv,
            "cluster_count_min": int(counts.min()),
            "cluster_count_max": int(counts.max()),
            "cluster_count_median": int(np.median(counts)),
        },
        "elapsed_seconds": time.time() - t0,
    }

    out_path = "analysis/outputs/p21_multiseed.json"
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, "w") as f:
        json.dump(output, f, indent=2)
    print(f"\nSaved → {out_path}  (elapsed {output['elapsed_seconds']:.1f}s)")


if __name__ == "__main__":
    main()
