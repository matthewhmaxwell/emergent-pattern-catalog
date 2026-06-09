"""
P23 dim2 — multi-seed stability at a canonical α.

Runs ≥20 seeds at α ≈ 0.63 (m=6, N=101) — the efficient phase — and
reports σ²/N mean ± std ± CV. Emits analysis/outputs/p23_multiseed.json.
"""

import json
import sys
import os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

from epc.models.minority_game import MinorityGame, MinorityGameParams
from epc.detectors.p23_anticoordination import P23AnticoordinationDetector


def run_multiseed() -> dict:
    """Run multi-seed campaign at canonical α."""
    N = 101
    m = 6
    n_rounds = 2000
    burn_in = 400
    n_seeds = 25

    svs = []
    ac1s = []
    tiers = []
    confidences = []

    det = P23AnticoordinationDetector(n_permutations=199, seed=42)

    for s in range(n_seeds):
        mg = MinorityGame(MinorityGameParams(
            n_agents=N, memory=m, n_rounds=n_rounds, seed=100 + s,
        ))
        history = mg.simulate()
        att = np.array([r['attendance'] for r in history[burn_in:]], dtype=float)

        sv = float(np.var(att, ddof=0) / N)
        xm = att - np.mean(att)
        c0 = np.sum(xm ** 2)
        ac1 = float(np.sum(xm[:-1] * xm[1:]) / c0) if c0 > 0 else 0.0

        result = det.detect(history, mg.get_metadata())

        svs.append(sv)
        ac1s.append(ac1)
        tiers.append(result.tier.value)
        confidences.append(result.confidence)

    sv_arr = np.array(svs)
    ac1_arr = np.array(ac1s)

    output = {
        'reproduction': 'p23_multiseed',
        'pattern': 'P23',
        'parameters': {
            'N': N,
            'memory': m,
            'alpha': round(2 ** m / N, 6),
            'n_rounds': n_rounds,
            'burn_in': burn_in,
            'n_seeds': n_seeds,
        },
        'scaled_variance': {
            'mean': round(float(np.mean(sv_arr)), 6),
            'std': round(float(np.std(sv_arr, ddof=1)), 6),
            'cv': round(float(np.std(sv_arr, ddof=1) / np.mean(sv_arr)), 4),
            'min': round(float(np.min(sv_arr)), 6),
            'max': round(float(np.max(sv_arr)), 6),
        },
        'lag1_autocorrelation': {
            'mean': round(float(np.mean(ac1_arr)), 6),
            'std': round(float(np.std(ac1_arr, ddof=1)), 6),
        },
        'detection': {
            'tiers': {
                str(k): int(v)
                for k, v in zip(*np.unique(tiers, return_counts=True))
            },
            'confidence_mean': round(float(np.mean(confidences)), 4),
            'confidence_std': round(float(np.std(confidences, ddof=1)), 4),
            'all_detected': all(t != 'screening' or c > 0 for t, c in zip(tiers, confidences)),
        },
        'per_seed': [
            {
                'seed': 100 + s,
                'scaled_variance': round(svs[s], 6),
                'lag1_autocorrelation': round(ac1s[s], 6),
                'tier': tiers[s],
                'confidence': round(confidences[s], 4),
            }
            for s in range(n_seeds)
        ],
    }

    return output


if __name__ == '__main__':
    result = run_multiseed()
    out_path = os.path.join(
        os.path.dirname(__file__), '..', 'outputs',
        'p23_multiseed.json',
    )
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, 'w') as f:
        json.dump(result, f, indent=2)
    print(json.dumps(result, indent=2))
    print(f"\nWritten to {out_path}")
