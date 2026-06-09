"""End-to-end tests for P23 anti-coordination detection on the Minority Game.

Covers:
- Deterministic from seed
- DEFINITIVE on canonical efficient regime (α ≈ 0.63, m=6, N=101)
- Negative control: random-choice agents → not detected
- Detector reads via T1a adapter
"""

import numpy as np
import pytest

from epc.models.minority_game import (
    MinorityGame,
    MinorityGameParams,
    ElFarolBar,
    ElFarolParams,
    RandomChoiceBaseline,
)
from epc.detectors.p23_anticoordination import (
    P23AnticoordinationDetector,
    extract_observation_bundle,
)


class TestMinorityGameDeterminism:
    """Model produces identical results across runs with same seed."""

    def test_deterministic_from_seed(self):
        params = MinorityGameParams(n_agents=101, memory=6, n_rounds=200, seed=42)
        mg1 = MinorityGame(params)
        h1 = mg1.simulate()
        mg2 = MinorityGame(MinorityGameParams(n_agents=101, memory=6, n_rounds=200, seed=42))
        h2 = mg2.simulate()
        att1 = [r['attendance'] for r in h1]
        att2 = [r['attendance'] for r in h2]
        assert att1 == att2, "MinorityGame is not deterministic from seed"

    def test_different_seeds_differ(self):
        mg1 = MinorityGame(MinorityGameParams(n_agents=101, memory=3, n_rounds=100, seed=1))
        mg2 = MinorityGame(MinorityGameParams(n_agents=101, memory=3, n_rounds=100, seed=2))
        att1 = [r['attendance'] for r in mg1.simulate()]
        att2 = [r['attendance'] for r in mg2.simulate()]
        assert att1 != att2, "Different seeds should produce different results"


class TestP23OnMinorityGame:
    """P23 detector on canonical Minority Game regimes."""

    def test_definitive_on_efficient_regime(self):
        """α ≈ 0.63 (m=6, N=101): efficient phase → DEFINITIVE."""
        mg = MinorityGame(MinorityGameParams(
            n_agents=101, memory=6, n_rounds=1000, seed=42,
        ))
        history = mg.simulate()
        det = P23AnticoordinationDetector(n_permutations=199, seed=42)
        result = det.detect(history, mg.get_metadata())
        assert result.detected, f"P23 not detected: {result.notes}"
        assert result.tier.value == "definitive", (
            f"Expected definitive, got {result.tier.value}: {result.notes}"
        )

    def test_variance_below_random_baseline(self):
        """In efficient regime, σ²/N should be well below 0.25."""
        mg = MinorityGame(MinorityGameParams(
            n_agents=101, memory=6, n_rounds=1000, seed=42,
        ))
        history = mg.simulate()
        att = np.array([r['attendance'] for r in history[200:]])
        sv = np.var(att) / 101
        assert sv < 0.20, f"σ²/N = {sv:.4f}, expected < 0.20"

    def test_negative_autocorrelation(self):
        """Efficient regime should show negative lag-1 autocorrelation."""
        mg = MinorityGame(MinorityGameParams(
            n_agents=101, memory=6, n_rounds=1000, seed=42,
        ))
        history = mg.simulate()
        att = np.array([r['attendance'] for r in history[200:]], dtype=float)
        xm = att - np.mean(att)
        ac1 = float(np.sum(xm[:-1] * xm[1:]) / np.sum(xm ** 2))
        assert ac1 < 0, f"Expected negative ac1, got {ac1:.4f}"

    def test_detector_reads_via_adapter(self):
        """T1a: detector reads observation bundle, not model object."""
        mg = MinorityGame(MinorityGameParams(
            n_agents=101, memory=6, n_rounds=500, seed=42,
        ))
        history = mg.simulate()
        bundle = extract_observation_bundle(history)
        assert 'attendance' in bundle
        assert 'n_agents' in bundle
        assert len(bundle['attendance']) == 500

    def test_confidence_at_definitive(self):
        """Definitive should have confidence ≥ 0.75."""
        mg = MinorityGame(MinorityGameParams(
            n_agents=101, memory=6, n_rounds=1000, seed=42,
        ))
        det = P23AnticoordinationDetector(n_permutations=199, seed=42)
        result = det.detect(mg.simulate(), mg.get_metadata())
        assert result.confidence >= 0.75, (
            f"Expected confidence ≥ 0.75, got {result.confidence:.2f}"
        )


class TestRandomAgentNegativeControl:
    """Random-choice agents (no adaptation) → detector does NOT fire."""

    def test_random_agents_not_detected(self):
        """Random baseline: σ² = random baseline, no anti-correlation → no detection."""
        rb = RandomChoiceBaseline(n_agents=101, n_rounds=1000, seed=42)
        history = rb.simulate()
        det = P23AnticoordinationDetector(n_permutations=199, seed=42)
        result = det.detect(history, rb.get_metadata())
        assert not result.detected, (
            f"P23 false positive on random agents: {result.summary()}"
        )

    def test_random_variance_near_baseline(self):
        """Random agents: σ²/N ≈ 0.25 (Binomial variance)."""
        rb = RandomChoiceBaseline(n_agents=101, n_rounds=5000, seed=42)
        history = rb.simulate()
        att = np.array([r['attendance'] for r in history])
        sv = np.var(att) / 101
        assert abs(sv - 0.25) < 0.05, f"Random σ²/N = {sv:.4f}, expected ≈ 0.25"


class TestSymmetricPhaseNegativeControl:
    """In the symmetric (maladaptive) phase (very small α), σ²/N > 0.25.
    The detector should NOT fire at definitive."""

    def test_symmetric_phase_not_definitive(self):
        """m=1, α=0.02: symmetric phase → σ²/N > random baseline."""
        mg = MinorityGame(MinorityGameParams(
            n_agents=101, memory=1, n_rounds=1000, seed=42,
        ))
        history = mg.simulate()
        det = P23AnticoordinationDetector(n_permutations=199, seed=42)
        result = det.detect(history, mg.get_metadata())
        # Either not detected or at most screening
        if result.detected:
            assert result.tier.value != "definitive", (
                f"Symmetric phase should not reach definitive: {result.notes}"
            )


class TestElFarolCrossModel:
    """P23 must also fire on El Farol (T1b cross-model test)."""

    def test_el_farol_detected(self):
        """El Farol variant should be detected by P23."""
        ef = ElFarolBar(ElFarolParams(
            n_agents=100, capacity=60, n_rounds=1000, seed=42,
        ))
        history = ef.simulate()
        det = P23AnticoordinationDetector(n_permutations=199, seed=42)
        result = det.detect(history, ef.get_metadata())
        assert result.detected, (
            f"P23 not detected on El Farol: {result.notes}"
        )
        assert result.tier.value in ("confirmation", "definitive"), (
            f"Expected confirmation+, got {result.tier.value}: {result.notes}"
        )

    def test_el_farol_deterministic(self):
        """El Farol is deterministic from seed."""
        ef1 = ElFarolBar(ElFarolParams(n_agents=100, capacity=60, n_rounds=200, seed=42))
        ef2 = ElFarolBar(ElFarolParams(n_agents=100, capacity=60, n_rounds=200, seed=42))
        att1 = [r['attendance'] for r in ef1.simulate()]
        att2 = [r['attendance'] for r in ef2.simulate()]
        assert att1 == att2, "ElFarolBar is not deterministic from seed"
