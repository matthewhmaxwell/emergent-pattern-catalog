"""Sprint 61 regression tests for P1 multi-cluster prerequisite (Schelling 1971).

Two tests:
1. Gradient substrate correctly rejected by the multi-cluster prerequisite.
2. Canonical Schelling positive still reaches CONFIRMATION (no native-domain regression).
"""

import numpy as np
import pytest

from epc.detectors.p1_aggregation import P1AggregationDetector
from epc.detector_result import DetectionTier
from epc.models.schelling import run_schelling
from epc.phase2a.synthetic import linear_gradient_field


def test_gradient_rejected_by_multi_cluster_prereq():
    """linear_gradient substrate must be rejected at SCREENING.

    A binarized gradient (left half type 0, right half type 1) has exactly
    1 connected component per type — a monotonic spatial partition, not
    multi-cluster Schelling aggregation. The multi-cluster prerequisite
    (Sprint 61, Schelling 1971) should catch this.
    """
    det = P1AggregationDetector(n_permutations=199)
    history = linear_gradient_field("grid", seed=42)
    result = det.detect(history)

    assert not result.detected, "gradient should NOT be detected as P1"
    assert result.tier == DetectionTier.SCREENING
    assert any("multi_cluster_failed" in w for w in result.warnings)


def test_canonical_schelling_positive_not_regressed():
    """Canonical Schelling positive (threshold=0.375, density=0.9) must
    still reach at least CONFIRMATION tier after multi-cluster prereq.

    Schelling segregation produces multiple disconnected same-type clusters
    (typically 6–28 per type on a 32×32 grid), so the prerequisite should
    pass trivially.
    """
    det = P1AggregationDetector(n_permutations=999)
    history = run_schelling(
        grid_size=32, density=0.9, threshold=0.375, n_steps=80, seed=0,
    )
    # run_schelling returns dicts with 'grid' as ndarray; add grid_dims
    adapted = []
    for s in history:
        grid = np.asarray(s["grid"])
        adapted.append({"grid": grid, "grid_dims": grid.shape})
    result = det.detect(adapted)

    assert result.detected, "Schelling positive must still be detected"
    assert result.tier.value >= DetectionTier.CONFIRMATION.value
    assert result.confidence >= 0.60
