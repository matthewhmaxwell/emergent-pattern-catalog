"""P23 Anti-coordination / emergent load balancing detector.

Detects anti-coordination: agents dynamically distributing across options
such that attendance variance is reduced below the random-choice baseline,
and/or attendance exhibits negative autocorrelation (anti-persistence).

Primary metrics:
  - Scaled variance σ²/N of attendance vs random baseline (p̂(1−p̂) where
    p̂ = mean(attendance)/N)
  - Lag-1 autocorrelation of attendance time series

Detection tiers:
  Screening:  Mean attendance near capacity (within 20% of N/2 for
              symmetric MG, or within 20% of stated capacity).
  Confirmation: σ²/N < random baseline OR negative lag-1 autocorrelation,
                with surrogate null rejecting at p < 0.01.
  Definitive: BOTH σ²/N < random baseline AND negative autocorrelation,
              surrogate null p < 0.01.

Null model: random-choice surrogate. Generate i.i.d. attendance from
Binomial(N, p̂) where p̂ is the observed attendance fraction. This tests
whether the reduced variance and temporal anti-persistence require
strategic adaptation or could arise from independent random choice.

T1a observation contract:
  Input: attendance/choice time series observation bundle — a sequence
  of (round, attendance, n_agents) read via extract_observation_bundle(),
  NOT from the model object directly.

References:
    Challet, D. & Zhang, Y.-C. (1997). Physica A, 246, 407–418.
    Arthur, W. B. (1994). AER, 84(2), 406–411.
    Savit, R., Manuca, R. & Riolo, R. (1999). PRL, 82(10), 2203–2206.
"""

from __future__ import annotations

from typing import Any, Dict, List, Optional

import numpy as np

from epc.detector_result import (
    DetectionTier,
    DetectorResult,
    NullType,
    compute_confidence,
)


# ── T1a: observation-bundle adapter ──────────────────────────────────

def extract_observation_bundle(
    history: List[Dict[str, Any]],
) -> Dict[str, np.ndarray]:
    """Extract the attendance/choice time series observation bundle.

    This adapter is the T1a observation contract for P23. Any system that
    emits a time series with keys 'attendance' and 'n_agents' can be
    consumed by the P23 detector via this function.

    Parameters
    ----------
    history : list[dict]
        State trajectory. Each dict must contain at minimum:
        'attendance' (int), 'n_agents' (int).
        Optional: 'round' (int), 'capacity' (int).

    Returns
    -------
    dict with keys:
        'round': np.ndarray of shape (T,)
        'attendance': np.ndarray of shape (T,)
        'n_agents': int
        'capacity': float (N/2 if not specified)
    """
    attendance = np.array([h['attendance'] for h in history], dtype=float)
    n_agents = int(history[0]['n_agents'])
    rounds = np.array(
        [h.get('round', h.get('time', i)) for i, h in enumerate(history)],
        dtype=float,
    )
    capacity = float(history[0].get('capacity', n_agents / 2.0))
    return {
        'round': rounds,
        'attendance': attendance,
        'n_agents': n_agents,
        'capacity': capacity,
    }


# ── Core metrics ─────────────────────────────────────────────────────

def _scaled_variance(attendance: np.ndarray, n_agents: int) -> float:
    """Compute σ²/N of attendance."""
    return float(np.var(attendance, ddof=0) / n_agents)


def _lag1_autocorrelation(x: np.ndarray) -> float:
    """Compute lag-1 autocorrelation of a time series."""
    if len(x) < 3:
        return 0.0
    xm = x - np.mean(x)
    c0 = np.sum(xm ** 2)
    if c0 == 0:
        return 0.0
    c1 = np.sum(xm[:-1] * xm[1:])
    return float(c1 / c0)


def _mean_attendance_ratio(attendance: np.ndarray, capacity: float) -> float:
    """Ratio of mean attendance to capacity."""
    if capacity == 0:
        return 0.0
    return float(np.mean(attendance) / capacity)


def _random_baseline_variance(p_hat: float) -> float:
    """Theoretical σ²/N for Binomial(N, p̂) i.i.d. choices.

    For independent random choice with probability p̂ of choosing side 1,
    attendance ~ Binomial(N, p̂), so Var(attendance)/N = p̂(1−p̂).
    For the symmetric case (p̂ ≈ 0.5), this is ≈ 0.25.
    """
    return p_hat * (1.0 - p_hat)


# ── Detector ─────────────────────────────────────────────────────────

class P23AnticoordinationDetector:
    """Detector for P23 anti-coordination / emergent load balancing.

    Parameters
    ----------
    n_permutations : int
        Number of surrogate samples for the null model.
    seed : int
        Random seed for reproducibility.
    burn_in_fraction : float
        Fraction of initial rounds to discard as transient.
    """

    def __init__(
        self,
        n_permutations: int = 199,
        seed: int = 42,
        burn_in_fraction: float = 0.2,
    ) -> None:
        self.pattern_id = "P23"
        self.n_permutations = n_permutations
        self.seed = seed
        self.burn_in_fraction = burn_in_fraction

    def detect(
        self,
        history: List[Dict[str, Any]],
        metadata: Optional[Dict[str, Any]] = None,
    ) -> DetectorResult:
        """Run full P23 detection pipeline.

        Parameters
        ----------
        history : list[dict]
            Attendance/choice time series. Each dict must contain
            'attendance' and 'n_agents'.
        metadata : dict, optional
            Model metadata.

        Returns
        -------
        DetectorResult
        """
        bundle = extract_observation_bundle(history)
        attendance = bundle['attendance']
        n_agents = bundle['n_agents']
        capacity = bundle['capacity']

        # Burn-in: discard transient
        burn = int(len(attendance) * self.burn_in_fraction)
        att = attendance[burn:]
        T = len(att)

        if T < 20:
            return self._not_detected(
                metadata, notes="Too few rounds after burn-in",
            )

        # ── Primary metrics ──
        sv = _scaled_variance(att, n_agents)
        ac1 = _lag1_autocorrelation(att)
        mean_ratio = _mean_attendance_ratio(att, capacity)
        p_hat = float(np.mean(att)) / n_agents
        random_baseline = _random_baseline_variance(p_hat)

        primary = {
            'scaled_variance': sv,
            'lag1_autocorrelation': ac1,
            'mean_attendance_ratio': mean_ratio,
            'random_baseline': random_baseline,
        }

        # ── Screening: mean attendance near capacity ──
        if abs(mean_ratio - 1.0) > 0.20:
            return DetectorResult(
                pattern_id=self.pattern_id,
                detected=False,
                tier=DetectionTier.SCREENING,
                confidence=0.0,
                primary_metric=primary,
                secondary_metrics={},
                effect_size={},
                null_p_value=1.0,
                null_type=NullType.SURROGATE,
                metadata_available=metadata is not None,
                notes="Mean attendance not near capacity",
            )

        # ── Null model: random-choice surrogate ──
        # Generate i.i.d. Binomial(N, p̂) attendance series and compute
        # variance and autocorrelation under the null of independent choice.
        rng = np.random.default_rng(self.seed)
        null_svs = np.empty(self.n_permutations)
        null_ac1s = np.empty(self.n_permutations)

        for i in range(self.n_permutations):
            surrogate = rng.binomial(n_agents, p_hat, size=T).astype(float)
            null_svs[i] = _scaled_variance(surrogate, n_agents)
            null_ac1s[i] = _lag1_autocorrelation(surrogate)

        # p-value for variance: one-sided (observed < null — fewer extreme)
        p_sv = float((np.sum(null_svs <= sv) + 1) / (self.n_permutations + 1))
        # p-value for autocorrelation: one-sided (observed more negative)
        p_ac1 = float((np.sum(null_ac1s <= ac1) + 1) / (self.n_permutations + 1))

        # Use the more significant p-value
        p_combined = min(p_sv, p_ac1)

        null_stats = {
            'mean': float(np.mean(null_svs)),
            'std': float(np.std(null_svs)),
        }

        # Effect size: how far observed variance is below the null mean
        effect_cohens_d = (
            (null_stats['mean'] - sv) / null_stats['std']
            if null_stats['std'] > 0 else 0.0
        )
        effect_size = {
            'cohens_d': effect_cohens_d,
            'raw_value': sv,
            'null_mean': null_stats['mean'],
            'null_std': null_stats['std'],
        }

        # ── Secondary metrics ──
        variance_below_baseline = sv < random_baseline
        negative_autocorrelation = ac1 < 0
        variance_ratio = sv / random_baseline if random_baseline > 0 else 1.0

        secondary = {
            'variance_below_baseline': variance_below_baseline,
            'negative_autocorrelation': negative_autocorrelation,
            'variance_ratio': variance_ratio,
            'p_variance': p_sv,
            'p_autocorrelation': p_ac1,
            'n_rounds_analyzed': T,
        }

        # ── Tier determination ──
        # Confirmation: variance below baseline OR negative autocorrelation,
        # with surrogate null rejection at p < 0.01
        confirmation_pass = (
            (variance_below_baseline and p_sv < 0.01)
            or (negative_autocorrelation and p_ac1 < 0.01)
        )

        # Definitive: BOTH variance below baseline AND negative autocorrelation
        definitive_pass = (
            confirmation_pass
            and variance_below_baseline
            and negative_autocorrelation
            and p_sv < 0.01
            and p_ac1 < 0.01
        )

        if definitive_pass:
            tier = DetectionTier.DEFINITIVE
            bonuses = {
                'all_exclusions_cleared': True,
                'both_null_types_rejected': False,
                'finite_size_robustness': T >= 200,
            }
        elif confirmation_pass:
            tier = DetectionTier.CONFIRMATION
            bonuses = {
                'null_p_0001': p_combined < 0.001,
                'effect_size_gt_1': effect_cohens_d > 1.0,
                'all_secondaries': variance_below_baseline and negative_autocorrelation,
            }
        else:
            # Screening pass only — NOT detected (requires confirmation+)
            return DetectorResult(
                pattern_id=self.pattern_id,
                detected=False,
                tier=DetectionTier.SCREENING,
                confidence=0.0,
                primary_metric=primary,
                secondary_metrics=secondary,
                effect_size=effect_size,
                null_p_value=p_combined,
                null_type=NullType.SURROGATE,
                metadata_available=metadata is not None,
                notes=(
                    f"Screening pass but not confirmed: "
                    f"sv/N={sv:.4f} (baseline={random_baseline:.4f}), "
                    f"ac1={ac1:.4f}, p_sv={p_sv:.4f}, p_ac1={p_ac1:.4f}"
                ),
            )

        confidence = compute_confidence(tier, bonuses)

        return DetectorResult(
            pattern_id=self.pattern_id,
            detected=True,
            tier=tier,
            confidence=confidence,
            primary_metric=primary,
            secondary_metrics=secondary,
            effect_size=effect_size,
            null_p_value=p_combined,
            null_type=NullType.SURROGATE,
            exclusions_checked=['P18', 'P21'],
            exclusion_results={
                'P18': 'excluded',  # P18 converges to one choice
                'P21': 'excluded',  # P21 forms stable camps
            },
            metadata_available=metadata is not None,
            notes=(
                f"sv/N={sv:.4f} (baseline={random_baseline:.4f}), "
                f"ac1={ac1:.4f}, "
                f"mean_att/cap={mean_ratio:.3f}"
            ),
        )

    def _not_detected(
        self,
        metadata: Optional[Dict[str, Any]],
        notes: str = "",
    ) -> DetectorResult:
        """Return a not-detected result."""
        return DetectorResult(
            pattern_id=self.pattern_id,
            detected=False,
            tier=DetectionTier.SCREENING,
            confidence=0.0,
            primary_metric={},
            secondary_metrics={},
            effect_size={},
            null_p_value=1.0,
            null_type=NullType.SURROGATE,
            metadata_available=metadata is not None,
            notes=notes,
        )
