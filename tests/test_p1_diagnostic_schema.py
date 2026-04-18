"""Sprint 14.5 D.3: P1 diagnostic-schema polish.

Ensures every P1 DetectorResult carries a non-trivial diagnostic story
even when screening rejects. Previously, SIR × P1 produced
``secondary_metrics = {}`` on early exit, with no in-result indication
of *why* screening rejected; and ``sustained_i_cv`` was reported as
``float('inf')`` when ``mean_i`` was near zero — correct numerically
but opaque about the underlying undefined-CV condition.

The D.3 polish:

* Adds ``screening_rejection_reason`` to ``primary_metric`` for every
  rejection path (uniform_state, below_expected, below_magnitude_floor,
  substrate_mismatch, empty_state_history). When the detector passes
  screening, the reason is ``"none"``.
* Adds ``sustained_i_cv_undefined`` to ``secondary_metrics``, which is
  ``True`` when ``|mean_i| < 0.01`` (the CV-undefined regime).
  ``sustained_i_cv`` is retained as ``float('inf')`` for backwards
  compatibility with the numeric ``cv > 0.3`` confirmation check.
"""

from __future__ import annotations

import math

import pytest

from epc.detectors.p1_aggregation import P1AggregationDetector
from epc.detector_result import DetectionTier


REJECTION_REASONS = {
    "none",
    "uniform_state",
    "below_expected",
    "below_magnitude_floor",
    "substrate_mismatch",
    "empty_state_history",
}


def test_p1_sir_rejects_with_uniform_state_reason():
    """SIR at t=T_end is uniform recovered → rejection_reason=uniform_state."""
    from epc.models.sir_epidemic import SIREpidemicModel

    m = SIREpidemicModel(
        rows=60, cols=60,
        infection_prob=0.2, recovery_prob=0.05,
        seed=42,
    )
    history = m.run(400)

    det = P1AggregationDetector(n_permutations=9)
    result = det.detect(history, model_metadata=m.get_metadata())

    assert not result.detected
    assert result.tier == DetectionTier.SCREENING
    assert "screening_rejection_reason" in result.primary_metric
    assert result.primary_metric["screening_rejection_reason"] == "uniform_state", (
        f"SIR post-wavefront should report uniform_state rejection, "
        f"got {result.primary_metric['screening_rejection_reason']!r}"
    )
    # Diagnostic: peak should be preserved even on screening reject
    assert result.primary_metric["morans_i_peak"] > 0.5, (
        "SIR wavefront peak Moran's I should be preserved as diagnostic"
    )


def test_p1_gs_rejects_with_substrate_mismatch_reason():
    """Gray-Scott (no integer grid) → rejection_reason=substrate_mismatch."""
    from epc.models.gray_scott import GrayScott

    m = GrayScott(rows=64, cols=64, feed_rate=0.037, kill_rate=0.060, seed=42)
    history = m.run(200)

    det = P1AggregationDetector(n_permutations=9)
    result = det.detect(history, model_metadata=m.get_metadata())

    assert result.primary_metric["screening_rejection_reason"] == "substrate_mismatch"


def test_p1_empty_history_rejects_with_empty_reason():
    """Empty state_history → rejection_reason=empty_state_history."""
    det = P1AggregationDetector(n_permutations=9)
    result = det.detect([], model_metadata=None)
    assert result.primary_metric["screening_rejection_reason"] == "empty_state_history"


def test_p1_rejection_reason_always_valid():
    """Every P1 result's screening_rejection_reason is in the allowed set."""
    from epc.models.sir_epidemic import SIREpidemicModel
    from epc.models.gray_scott import GrayScott

    det = P1AggregationDetector(n_permutations=9)

    # Three inputs exercising three different rejection paths
    cases = [
        # Empty
        [],
        # Substrate mismatch (GS)
        GrayScott(rows=32, cols=32, feed_rate=0.037, kill_rate=0.060, seed=42).run(100),
        # Uniform state (SIR)
        SIREpidemicModel(
            rows=40, cols=40, infection_prob=0.2, recovery_prob=0.05, seed=42
        ).run(300),
    ]

    for history in cases:
        result = det.detect(history, model_metadata=None)
        reason = result.primary_metric.get("screening_rejection_reason")
        assert reason in REJECTION_REASONS, (
            f"Invalid rejection reason {reason!r}; "
            f"must be one of {REJECTION_REASONS}"
        )


def test_p1_sustained_cv_undefined_flag_on_uniform_grid():
    """SIR post-collapse has |mean_i| < 0.01 → sustained_i_cv_undefined=True.

    Force-calls _compute_secondaries since SIR rejects at screening and
    secondaries are normally skipped. This is a direct test of the
    secondaries-dict schema on the undefined-CV regime.
    """
    from epc.models.sir_epidemic import SIREpidemicModel

    m = SIREpidemicModel(
        rows=60, cols=60,
        infection_prob=0.2, recovery_prob=0.05,
        seed=42,
    )
    history = m.run(400)

    det = P1AggregationDetector(n_permutations=9)
    sec = det._compute_secondaries(history, 10.0)

    assert "sustained_i_cv_undefined" in sec
    assert sec["sustained_i_cv_undefined"] is True, (
        f"SIR post-collapse should flag CV as undefined; "
        f"sustained_i_mean={sec['sustained_i_mean']:.6f}"
    )
    # The legacy numeric field stays at inf (for backwards-compat with
    # the `cv > 0.3` confirmation check)
    assert sec["sustained_i_cv"] == float("inf")


def test_p1_schelling_canonical_positive_has_none_reason():
    """Canonical positive Schelling × P1 must report rejection_reason=none."""
    from epc.models.schelling import run_schelling

    history = run_schelling(
        grid_size=50,
        density=0.9,
        threshold=0.375,
        n_steps=200,
        seed=42,
    )

    det = P1AggregationDetector(n_permutations=9)
    result = det.detect(history, model_metadata={"name": "schelling"})

    assert result.detected, "Schelling × P1 should pass screening"
    assert result.primary_metric["screening_rejection_reason"] == "none", (
        f"Canonical positive should report 'none' rejection; "
        f"got {result.primary_metric['screening_rejection_reason']!r}"
    )


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
