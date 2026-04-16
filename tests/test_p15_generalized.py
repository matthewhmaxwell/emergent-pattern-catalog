"""Tests for the generalized P15 persistent computation detector.

Tests verify that the substrate-independent P15 detector:
1. DEFINITIVE on GoL dense random IC (canonical positive)
2. Does not detect GH broken wave (deterministic but no diversity)
3. Does not detect SIR (stochastic — fails reproducibility)
4. Does not detect Schelling without step_fn
5. SCREENING on Nowak-May (deterministic, partial diversity)
6. Correctly catches stochastic models at intermediate checkpoints

Reference for GoL computation:
  Langton (1990), Rendell (2002), Lizier et al. (2012).
"""

import numpy as np
import pytest

from epc.detectors.p15_persistent_computation import (
    P15PersistentComputationDetector,
    make_step_fn_for_gol,
    make_step_fn_from_model,
    _classify_outcome,
)
from epc.detector_result import DetectionTier


# ============================================================
# Outcome classifier unit tests
# ============================================================


class TestOutcomeClassifier:
    """_classify_outcome must correctly identify pattern types."""

    def test_dead_grid(self):
        g = np.zeros((5, 5), dtype=int)
        assert _classify_outcome([g, g, g]) == "dead"

    def test_static_pattern(self):
        g = np.zeros((5, 5), dtype=int)
        g[1:3, 1:3] = 1
        assert _classify_outcome([g, g, g, g]) == "static"

    def test_period_2(self):
        g1 = np.zeros((5, 5), dtype=int)
        g1[2, 1:4] = 1
        g2 = np.zeros((5, 5), dtype=int)
        g2[1:4, 2] = 1
        assert _classify_outcome([g1, g2, g1, g2, g1]) == "period_2"

    def test_moving_object(self):
        frames = []
        for k in range(8):
            g = np.zeros((20, 20), dtype=int)
            g[3 + k, 3 + k] = 1
            frames.append(g)
        assert _classify_outcome(frames) == "moving"


# ============================================================
# GoL canonical positive
# ============================================================


class TestGoLP15Definitive:
    """GoL with dense random IC must achieve DEFINITIVE P15 detection.

    Dense random GoL is the canonical computing regime: deterministic,
    with diverse outcomes under input perturbation (still lifes,
    oscillators, moving structures from different random ICs).
    """

    def test_gol_dense_definitive(self):
        from epc.models.game_of_life import GameOfLife

        m = GameOfLife(
            rows=40, cols=40, init_mode="random",
            init_density=0.37, boundary="periodic", seed=42,
        )
        history = m.run(300, record_every=1)
        metadata = m.get_metadata()

        step_fn = make_step_fn_for_gol()
        det = P15PersistentComputationDetector(
            step_fn=step_fn, n_variations=12, seed=42,
        )
        result = det.detect(history, metadata)

        assert result.detected, "GoL dense random should trigger P15"
        assert result.tier >= DetectionTier.DEFINITIVE, \
            f"Expected DEFINITIVE, got {result.tier.name}"
        assert result.primary_metric["reproducibility"] == 1.0, \
            "GoL is deterministic — replay must match exactly"
        assert result.primary_metric["n_distinct_outcomes"] >= 3.0, \
            f"Expected ≥3 outcome classes, got {result.primary_metric['n_distinct_outcomes']}"

    def test_gol_reproducibility_is_exact(self):
        """GoL step function must reproduce the trajectory bit-exactly."""
        from epc.models.game_of_life import GameOfLife

        m = GameOfLife(
            rows=30, cols=30, init_mode="random",
            init_density=0.37, boundary="periodic", seed=99,
        )
        history = m.run(100, record_every=1)

        step_fn = make_step_fn_for_gol()
        det = P15PersistentComputationDetector(
            step_fn=step_fn, n_variations=4, seed=42,
        )
        result = det.detect(history, m.get_metadata())

        assert result.primary_metric["reproducibility"] == 1.0


# ============================================================
# Cross-model rejections
# ============================================================


class TestP15CrossDetection:
    """P15 must correctly reject non-computing models."""

    def test_gh_not_detected(self):
        """GH excitable waves are deterministic but show no outcome diversity.

        All perturbations of a broken wave produce the same spiral
        pattern → n_distinct = 1 → fails screening.
        """
        from epc.models.greenberg_hastings import GreenbergHastings

        m = GreenbergHastings(
            rows=60, cols=60, n_states=5, threshold=1,
            init_mode="broken_wave", boundary="periodic", seed=42,
        )
        history = m.run(200, record_every=1, vectorized=True)

        gh_step = make_step_fn_from_model(m)
        det = P15PersistentComputationDetector(
            step_fn=gh_step, n_variations=8, seed=42,
        )
        result = det.detect(history, m.get_metadata())

        assert not result.detected, "GH spirals should not be P15-detected"

    def test_sir_fails_reproducibility(self):
        """SIR is stochastic: replay at intermediate checkpoints must differ.

        The final state (all recovered) may match, but intermediate
        wavefront positions differ because the RNG path diverges.
        """
        from epc.models.sir_epidemic import SIREpidemicModel

        m = SIREpidemicModel(
            rows=40, cols=40, infection_prob=0.3, recovery_prob=0.1,
            init_mode="single_seed", seed=42,
        )
        history = m.run(150, record_every=1)

        sir_step = make_step_fn_from_model(m)
        det = P15PersistentComputationDetector(
            step_fn=sir_step, n_variations=4, seed=42,
        )
        result = det.detect(history, m.get_metadata())

        assert result.primary_metric["reproducibility"] < 0.9, \
            f"SIR replay should fail: repro={result.primary_metric['reproducibility']}"
        assert not result.detected

    def test_schelling_no_step_fn(self):
        """Schelling without step_fn has reproducibility=0 → not detected."""
        from epc.models.schelling import run_schelling

        history = run_schelling(grid_size=40, n_steps=200, seed=42)
        det = P15PersistentComputationDetector(
            step_fn=None, n_variations=4, seed=42,
        )
        result = det.detect(history, None)

        assert result.primary_metric["reproducibility"] == 0.0
        assert not result.detected

    def test_nowak_may_screening_only(self):
        """Nowak-May is deterministic but shows limited outcome diversity.

        Expected: SCREENING (repro=1.0, n_distinct≥2) but not CONFIRMATION
        (needs ≥3 distinct outcomes for that).
        """
        from epc.models.nowak_may import NowakMayModel

        m = NowakMayModel(rows=40, cols=40, b=1.8, init_mode="random", seed=42)
        m.setup()
        history = m.run(100)

        nm_step = make_step_fn_from_model(m)
        det = P15PersistentComputationDetector(
            step_fn=nm_step, n_variations=8, seed=42,
        )
        result = det.detect(history, m.get_metadata())

        assert result.primary_metric["reproducibility"] == 1.0, \
            "Nowak-May is deterministic — should reproduce exactly"
        assert result.tier <= DetectionTier.SCREENING, \
            f"Nowak-May should be at most SCREENING, got {result.tier.name}"
