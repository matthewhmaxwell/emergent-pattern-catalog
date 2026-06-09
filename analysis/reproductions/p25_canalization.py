"""
P25 dim1 reproduction — Canalized restoration / equifinality (Waddington 1957).

Anchor:
    Waddington, C. H. (1957). The Strategy of the Genes. Allen & Unwin.
    Huang, S. et al. (2005). Cell Fates as High-Dimensional Attractor
    States of a Complex Gene Regulatory Network. PRL 94, 128701.

    Canonical result: a system with a canalized basin funnels diverse
    initial conditions to the same target macrostate. The convergence
    variance ratio (Var(finals) / Var(ICs)) is ≪ 1, and the basin
    volume (fraction of ICs that converge) is large (≥ 0.8).

Reproduction targets:
    1. Convergence variance ratio < 0.1.
       The final-state variance across diverse ICs is less than 10%
       of the initial-state variance.
       Tolerance: ratio < 0.10.

    2. Basin volume ≥ 0.8.
       At least 80% of randomly sampled ICs converge to the target.
       Tolerance: basin_volume ≥ 0.80.

    3. Detector reaches at least CONFIRMATION tier.
       Tolerance: tier ∈ {confirmation, definitive}.

Parameters:
    CanalizedLandscape: n_dims=10, basin_strength=2.0, ic_spread=5.0,
    n_ics=20, n_steps=200, dt=0.05, seed=42.
"""

import json
import sys
import os

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', '..'))

from epc.models.canalization import (
    CanalizedLandscape,
    CanalizedLandscapeParams,
)
from epc.detectors.p25_equifinality import (
    P25EquifinalityDetector,
    extract_observation_bundle,
    _convergence_variance_ratio,
    _basin_volume,
)


def run_reproduction() -> dict:
    """Run the P25 canalized restoration reproduction."""

    # ── Canalized system ──
    params = CanalizedLandscapeParams(
        n_dims=10,
        n_ics=20,
        n_steps=200,
        basin_strength=2.0,
        ic_spread=5.0,
        dt=0.05,
        seed=42,
    )
    model = CanalizedLandscape(params)
    history = model.simulate()

    # Extract bundle
    bundle = extract_observation_bundle(history)
    ics = bundle['ics']
    finals = bundle['finals']
    converged = bundle['converged']

    # Compute metrics
    conv_ratio = _convergence_variance_ratio(ics, finals)
    basin_vol = _basin_volume(converged)

    ic_var = float(np.sum(np.var(ics, axis=0)))
    final_var = float(np.sum(np.var(finals, axis=0)))

    # ── Run detector ──
    det = P25EquifinalityDetector(n_permutations=199, seed=42)
    result = det.detect(history, model.get_metadata())

    # ── Tolerance checks ──
    tol1_pass = conv_ratio < 0.10
    tol2_pass = basin_vol >= 0.80
    tol3_pass = result.tier.value in ("confirmation", "definitive")

    passes_tolerance = tol1_pass and tol2_pass and tol3_pass

    output = {
        "reproduction": "p25_canalization",
        "pattern": "P25",
        "anchor": "Waddington (1957) — canalized restoration / equifinality",
        "parameters": {
            "n_dims": params.n_dims,
            "basin_strength": params.basin_strength,
            "ic_spread": params.ic_spread,
            "n_ics": params.n_ics,
            "n_steps": params.n_steps,
            "dt": params.dt,
            "seed": params.seed,
        },
        "ic_variance": round(ic_var, 6),
        "final_variance": round(final_var, 10),
        "convergence_variance_ratio": round(conv_ratio, 8),
        "basin_volume": round(basin_vol, 4),
        "detector_tier": result.tier.value,
        "detector_confidence": round(result.confidence, 4),
        "detector_p_value": round(result.null_p_value, 6),
        "tolerance_checks": {
            "convergence_ratio_lt_0.10": {
                "value": round(conv_ratio, 8),
                "threshold": 0.10,
                "pass": tol1_pass,
            },
            "basin_volume_gte_0.80": {
                "value": round(basin_vol, 4),
                "threshold": 0.80,
                "pass": tol2_pass,
            },
            "detector_confirmation_plus": {
                "tier": result.tier.value,
                "pass": tol3_pass,
            },
        },
        "passes_tolerance": passes_tolerance,
    }

    return output


if __name__ == "__main__":
    result = run_reproduction()
    out_path = os.path.join(
        os.path.dirname(__file__), '..', 'outputs',
        'p25_canalization_reproduction.json',
    )
    os.makedirs(os.path.dirname(out_path), exist_ok=True)
    with open(out_path, 'w') as f:
        json.dump(result, f, indent=2)
    print(json.dumps(result, indent=2))
    print(f"\nWritten to {out_path}")
    print(f"passes_tolerance: {result['passes_tolerance']}")
