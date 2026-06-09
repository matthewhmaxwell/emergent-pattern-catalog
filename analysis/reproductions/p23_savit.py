"""
P23 dim1 reproduction — Savit curve (Minority Game).

Anchor:
    Savit, R., Manuca, R. & Riolo, R. (1999). Adaptive competition,
    market efficiency, and phase transitions. Physical Review Letters,
    82(10), 2203–2206.

    Canonical result: Scaled fluctuation σ²/N as a function of the control
    parameter α = 2^m / N shows a characteristic minimum (efficient phase)
    below the random-choice baseline σ²/N = 1/4, with σ²/N rising toward
    (and above) the baseline in the symmetric phase (small α).

Reproduction targets:
    1. Interior minimum: σ²/N has a clear interior minimum as a function
       of α (not at the boundary of the scanned range).
    2. Below random baseline: The minimum σ²/N < 0.25.
    3. Symmetric phase: At least one small-α point has σ²/N > 0.25.

Parameters:
    N = 101 (odd, standard)
    m = 1..11 → α = 2^m / 101
    n_rounds = 3000 (burn-in 600, analyze last 2400)
    Seeds: average over 10 seeds for smooth curve
"""

import json
import sys
import os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

from epc.models.minority_game import MinorityGame, MinorityGameParams


def run_reproduction() -> dict:
    """Reproduce the Savit curve: σ²/N vs α = 2^m / N."""
    N = 101
    m_values = list(range(1, 12))
    n_rounds = 3000
    burn_in = 600
    n_seeds = 10

    results_by_m = []

    for m in m_values:
        alpha = 2 ** m / N
        svs = []
        for s in range(n_seeds):
            mg = MinorityGame(MinorityGameParams(
                n_agents=N, memory=m, n_rounds=n_rounds, seed=42 + s,
            ))
            history = mg.simulate()
            att = np.array([r['attendance'] for r in history[burn_in:]])
            sv = float(np.var(att, ddof=0) / N)
            svs.append(sv)

        mean_sv = float(np.mean(svs))
        std_sv = float(np.std(svs, ddof=1))

        results_by_m.append({
            'memory': m,
            'alpha': round(alpha, 6),
            'scaled_variance_mean': round(mean_sv, 6),
            'scaled_variance_std': round(std_sv, 6),
            'n_seeds': n_seeds,
        })

    # Extract the curve
    alphas = [r['alpha'] for r in results_by_m]
    sv_means = [r['scaled_variance_mean'] for r in results_by_m]

    # Find the minimum
    min_idx = int(np.argmin(sv_means))
    min_alpha = alphas[min_idx]
    min_sv = sv_means[min_idx]

    # Tolerance checks
    random_baseline = 0.25
    # 1. Interior minimum (not at first or last point)
    is_interior = 0 < min_idx < len(sv_means) - 1
    # 2. Below random baseline
    below_baseline = min_sv < random_baseline
    # 3. Symmetric phase above baseline
    symmetric_above = any(sv > random_baseline for sv in sv_means[:min_idx + 1])

    passes_tolerance = is_interior and below_baseline and symmetric_above

    output = {
        'reproduction': 'p23_savit_curve',
        'pattern': 'P23',
        'anchor': 'Savit, Manuca & Riolo (1999) — σ²/N vs α = 2^m/N',
        'parameters': {
            'N': N,
            'm_range': [1, 11],
            'n_rounds': n_rounds,
            'burn_in': burn_in,
            'n_seeds': n_seeds,
        },
        'curve': results_by_m,
        'minimum': {
            'alpha': min_alpha,
            'memory': m_values[min_idx],
            'scaled_variance': round(min_sv, 6),
        },
        'random_baseline': random_baseline,
        'tolerance_checks': {
            'interior_minimum': {
                'min_index': min_idx,
                'n_points': len(sv_means),
                'pass': is_interior,
            },
            'below_random_baseline': {
                'min_sv': round(min_sv, 6),
                'threshold': random_baseline,
                'pass': below_baseline,
            },
            'symmetric_phase_above_baseline': {
                'max_sv_small_alpha': round(max(sv_means[:min_idx + 1]), 6),
                'pass': symmetric_above,
            },
        },
        'passes_tolerance': passes_tolerance,
    }

    return output


if __name__ == '__main__':
    result = run_reproduction()
    out_path = os.path.join(
        os.path.dirname(__file__), '..', 'outputs',
        'p23_savit_reproduction.json',
    )
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, 'w') as f:
        json.dump(result, f, indent=2)
    print(json.dumps(result, indent=2))
    print(f"\nWritten to {out_path}")
    print(f"passes_tolerance: {result['passes_tolerance']}")
