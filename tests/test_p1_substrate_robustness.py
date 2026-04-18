"""Sprint 14 B.1: P1 substrate robustness.

Regression test for the carry-forward from Sprint 13: P1 previously
raised ``KeyError: "Need 'type_labels_at_pos' or 'grid' for 2D"`` when
invoked on a state_history lacking integer-grid observables (canonical
example: Gray-Scott exposes ``field`` but no ``grid`` / ``type_labels_at_pos``).

The Sprint 14 fix hardens P1's 2D branch to match the pattern used by
P11/P13/P22: emit a substrate-mismatch warning and return benign zeros
from ``_compute_primary`` so the detector short-circuits on screening
with ``detected=False``, ``tier=SCREENING``, and no traceback.
"""

from __future__ import annotations

import pytest

from epc.detectors.p1_aggregation import P1AggregationDetector
from epc.detector_result import DetectionTier


def test_p1_gs_substrate_rejected_cleanly():
    """Gray-Scott × P1 rejects cleanly without raising KeyError.

    Gray-Scott state dicts carry ``field`` (continuous float) but neither
    ``grid`` nor ``type_labels_at_pos``. The pre-Sprint-14 P1 raised
    KeyError from ``aggregation._extract_labels_and_adj``. Post-fix, P1
    returns a graceful detected=False result with a substrate warning.
    """
    from epc.models.gray_scott import GrayScott

    m = GrayScott(
        rows=64, cols=64,
        feed_rate=0.037, kill_rate=0.060,
        seed=42,
    )
    history = m.run(200)  # short run — we only exercise the substrate guard

    det = P1AggregationDetector(n_permutations=99)
    # Must not raise
    result = det.detect(history, model_metadata=m.get_metadata())

    assert not result.detected, \
        f"P1 should reject GS (no integer grid), got detected={result.detected}"
    assert result.tier == DetectionTier.SCREENING, \
        f"P1 should short-circuit at screening, got tier={result.tier}"
    assert any("P1 needs integer-labeled spatial data" in w for w in result.warnings), \
        f"P1 should warn about missing integer labels, got warnings={result.warnings}"
    # Primary metric should be benign zeros, not partially computed
    assert result.primary_metric["morans_i"] == 0.0
    assert result.primary_metric["n_unique_types"] == 0


def test_p1_empty_state_history_rejected_cleanly():
    """P1 on empty state_history must not raise."""
    det = P1AggregationDetector(n_permutations=9)
    result = det.detect([], model_metadata=None)
    assert not result.detected
    assert result.tier == DetectionTier.SCREENING


def test_p1_raw_field_only_dict_rejected_cleanly():
    """P1 on a minimal ``{'field': ..., 'grid_dims': ...}`` state.

    Direct unit test that doesn't depend on the Gray-Scott model plumbing;
    confirms the guard triggers on the observable-shape alone.
    """
    import numpy as np

    # Bare-minimum field-only state (no grid / type_labels_at_pos / cell_types)
    fake_history = [
        {"field": np.random.rand(16, 16).astype(float), "grid_dims": (16, 16)},
    ] * 10  # repeat to satisfy any run-length prereqs

    det = P1AggregationDetector(n_permutations=9)
    result = det.detect(fake_history, model_metadata={"name": "synthetic_field"})

    assert not result.detected
    assert result.tier == DetectionTier.SCREENING
    assert any("integer-labeled" in w for w in result.warnings)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
