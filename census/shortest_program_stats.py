"""Post hoc agreement statistics for the census table of shortest programs (RESULTS_census_v1.md section 2(a)).

These were computed after the registered analyses and are reported in the catalog paper as such:
  - Spearman rho between the shortest-program lengths in the two grammars over the ten named behaviors,
    with a percentile bootstrap 95% CI over behaviors;
  - rho with one behavior kept per credited catalog pattern (P9 and P18 each account for two behaviors), for every
    choice of the kept behavior, with an exact one-sided permutation p-value.

  python census/shortest_program_stats.py
"""
import itertools
import numpy as np
from scipy.stats import spearmanr

# behavior, catalog pattern credited, shortest program in grammar 1 (bits), in grammar 2 (bits)
TABLE = [
    ("Network synchronization", "P9", 9, 11),
    ("Lattice phase locking", "P9", 9, 11),
    ("Network consensus", "P18", 10, 14),
    ("Network fragmentation", "P34", 11, 14),
    ("Domain coarsening", "P18", 11, 14),
    ("Segregation with vacancies", "P1", 11, 18),
    ("Cyclic waves and spirals", None, 13, 18),
    ("Flocking", "P5", 14, 17),
    ("Power-law avalanches", "P14", 15, 18),
    ("Spatial game chaos", "P27", 15, 16),
]
N_BOOT, SEED = 10_000, 0


def rho(rows):
    return spearmanr([r[2] for r in rows], [r[3] for r in rows]).statistic


def bootstrap_ci(rows, rng):
    vals = []
    for _ in range(N_BOOT):
        s = [rows[i] for i in rng.integers(0, len(rows), len(rows))]
        r = rho(s)
        if np.isfinite(r):
            vals.append(r)
    return np.percentile(vals, [2.5, 97.5]), len(vals)


def perm_p(rows):
    """exact one-sided permutation p: share of all orderings of grammar 2 with rho at least the observed value."""
    g1 = np.array([r[2] for r in rows]); g2 = np.array([r[3] for r in rows])
    obs = spearmanr(g1, g2).statistic
    null = np.array([spearmanr(g1, g2[list(o)]).statistic for o in itertools.permutations(range(len(g2)))])
    return float(np.mean(null >= obs - 1e-12))


def main():
    rng = np.random.default_rng(SEED)
    r = rho(TABLE); (lo, hi), n_ok = bootstrap_ci(TABLE, rng)
    print(f"all ten behaviors: rho = {r:.3f}, bootstrap 95% CI {lo:.2f} to {hi:.2f} ({n_ok} finite resamples)")
    groups = {}
    for row in TABLE:
        groups.setdefault(row[1] or row[0], []).append(row)
    for choice in itertools.product(*groups.values()):
        kept = list(choice)
        print(f"one per pattern (n = {len(kept)}; kept: {', '.join(k[0] for k in kept if len(groups[k[1] or k[0]]) > 1)}): "
              f"rho = {rho(kept):.3f}, exact one-sided permutation p = {perm_p(kept):.3f}")


if __name__ == "__main__":
    main()
