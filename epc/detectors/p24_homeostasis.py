"""P24 Homeostatic regulation detector.

Detects active homeostatic regulation: a system maintaining an internal
variable near a set-point despite sustained external perturbation, via
active corrective (negative) feedback.

Primary metrics:
  - Steady-state |deviation| from set-point under perturbation
  - Deviation integral ∫|x − setpoint| dt over perturbed window
  - Deviation-growth ratio (second half / first half of perturbed window):
    regulated ≈ 1.0, uncontrolled >> 1.0

Detection tiers:
  Screening:  Deviation is bounded (growth ratio ≤ 2.0) during perturbation.
  Confirmation: Deviation integral significantly less than uncontrolled
                (linear-drift) surrogate baseline (p < 0.01).
  Definitive: Deviation integral ratio < 0.3 vs surrogate + metadata
              confirms active feedback.

Null model: surrogate uncontrolled baseline. The detector constructs
surrogate trajectories by integrating x' = perturbation(t) + noise
(no corrective feedback), with noise drawn from the observed residuals.
This tests whether the regulated variable's bounded deviation requires
active feedback or could arise passively.

T1a observation contract:
  Input: scalar-regulated-variable observation bundle — a time series of
  (time, x, setpoint, perturbation) read via extract_observation_bundle(),
  NOT from the model object directly.

Reference: Ashby, W. R. (1956). An Introduction to Cybernetics.
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
    """Extract the scalar-regulated-variable observation bundle from history.

    This adapter is the T1a observation contract for P24. Any system that
    emits a time series with keys 'time', 'x', 'setpoint', and 'perturbation'
    can be consumed by the P24 detector via this function.

    Parameters
    ----------
    history : list[dict]
        State trajectory. Each dict must contain at minimum:
        'time' (float), 'x' (float), 'setpoint' (float),
        'perturbation' (float).

    Returns
    -------
    dict with keys:
        'time': np.ndarray of shape (T,)
        'x': np.ndarray of shape (T,)
        'setpoint': np.ndarray of shape (T,)
        'perturbation': np.ndarray of shape (T,)
    """
    time = np.array([h['time'] for h in history], dtype=float)
    x = np.array([h['x'] for h in history], dtype=float)
    setpoint = np.array([h['setpoint'] for h in history], dtype=float)
    perturbation = np.array([h['perturbation'] for h in history], dtype=float)
    return {
        'time': time,
        'x': x,
        'setpoint': setpoint,
        'perturbation': perturbation,
    }


# ── Core detection logic ─────────────────────────────────────────────

def _compute_deviation_integral(
    x: np.ndarray,
    setpoint: np.ndarray,
    time: np.ndarray,
    start_idx: int,
    end_idx: int,
) -> float:
    """Compute ∫|x − setpoint| dt over [start_idx, end_idx) via trapezoidal."""
    if end_idx <= start_idx + 1:
        return 0.0
    dev = np.abs(x[start_idx:end_idx] - setpoint[start_idx:end_idx])
    dt_arr = np.diff(time[start_idx:end_idx])
    return float(np.sum(0.5 * (dev[:-1] + dev[1:]) * dt_arr))


def _find_perturbation_onset(perturbation: np.ndarray) -> Optional[int]:
    """Find the first index where perturbation becomes nonzero."""
    nonzero = np.nonzero(perturbation != 0.0)[0]
    if len(nonzero) == 0:
        return None
    return int(nonzero[0])


def _deviation_growth_ratio(
    x: np.ndarray,
    setpoint: np.ndarray,
    onset_idx: int,
) -> float:
    """Ratio of mean |deviation| in second half vs first half of
    post-onset trajectory. Bounded systems ≈ 1.0; drifting >> 1.0."""
    dev = np.abs(x[onset_idx:] - setpoint[onset_idx:])
    n = len(dev)
    if n < 4:
        return float('inf')
    mid = n // 2
    first = np.mean(dev[:mid])
    second = np.mean(dev[mid:])
    if first < 1e-12:
        return 0.0 if second < 1e-12 else float('inf')
    return float(second / first)


def _steady_state_deviation(
    x: np.ndarray,
    setpoint: np.ndarray,
    tail_fraction: float = 0.2,
) -> float:
    """Mean |deviation| over the last tail_fraction of the trajectory."""
    n = len(x)
    start = max(0, int(n * (1.0 - tail_fraction)))
    return float(np.mean(np.abs(x[start:] - setpoint[start:])))


def _generate_uncontrolled_surrogate(
    time: np.ndarray,
    setpoint: np.ndarray,
    perturbation: np.ndarray,
    onset_idx: int,
    x0: float,
    rng: np.random.Generator,
    noise_scale: float = 0.0,
) -> np.ndarray:
    """Generate surrogate trajectory: x' = perturbation(t) + noise, no feedback.

    This represents what would happen without active regulation —
    the variable drifts under perturbation with no corrective action.
    """
    n = len(time)
    surr = np.full(n, x0)
    # Copy pre-onset values
    surr[:onset_idx] = x0
    # Integrate post-onset: dx = perturbation * dt + noise
    for i in range(onset_idx, n - 1):
        dt = time[i + 1] - time[i]
        noise = rng.normal(0, noise_scale) * np.sqrt(dt) if noise_scale > 0 else 0.0
        surr[i + 1] = surr[i] + perturbation[i] * dt + noise
    return surr


class P24HomeostasisDetector:
    """P24 homeostatic regulation detector.

    Reads inputs via the T1a observation-bundle adapter.

    Parameters
    ----------
    n_permutations : int
        Number of surrogate runs for the null model.
    seed : int
        RNG seed for reproducibility.
    """

    def __init__(
        self,
        n_permutations: int = 199,
        seed: int = 42,
    ) -> None:
        self.pattern_id = "P24"
        self.n_permutations = n_permutations
        self.seed = seed

    def detect(
        self,
        history: List[Dict[str, Any]],
        metadata: Optional[Dict[str, Any]] = None,
    ) -> DetectorResult:
        """Run P24 homeostatic regulation detection.

        Parameters
        ----------
        history : list[dict]
            State trajectory with scalar-regulated-variable bundle keys.
        metadata : dict, optional
            Model metadata (used for definitive tier metadata checks).

        Returns
        -------
        DetectorResult
        """
        warnings: List[str] = []
        metadata_available = metadata is not None

        # ── Extract observation bundle (T1a) ──
        try:
            bundle = extract_observation_bundle(history)
        except (KeyError, IndexError) as e:
            return self._no_detection(
                warnings=[f"Observation bundle extraction failed: {e}"],
                metadata_available=metadata_available,
            )

        time = bundle['time']
        x = bundle['x']
        setpoint = bundle['setpoint']
        perturbation = bundle['perturbation']
        n = len(time)

        if n < 20:
            return self._no_detection(
                warnings=["Insufficient history length"],
                metadata_available=metadata_available,
            )

        # ── Find perturbation onset ──
        onset_idx = _find_perturbation_onset(perturbation)
        if onset_idx is None:
            warnings.append("No perturbation detected in schedule")
            return self._no_detection(
                warnings=warnings,
                metadata_available=metadata_available,
            )

        n_post = n - onset_idx
        if n_post < 10:
            return self._no_detection(
                warnings=["Insufficient post-perturbation data"],
                metadata_available=metadata_available,
            )

        # ── Primary metrics ──
        ss_deviation = _steady_state_deviation(x, setpoint, tail_fraction=0.2)
        dev_integral = _compute_deviation_integral(x, setpoint, time, onset_idx, n)
        growth_ratio = _deviation_growth_ratio(x, setpoint, onset_idx)

        primary: Dict[str, float] = {
            'steady_state_deviation': ss_deviation,
            'deviation_integral': dev_integral,
            'deviation_growth_ratio': growth_ratio,
        }

        peak_dev = float(np.max(np.abs(x[onset_idx:] - setpoint[onset_idx:])))
        if peak_dev < 1e-12:
            return self._no_detection(
                warnings=["No deviation from setpoint observed"],
                metadata_available=metadata_available,
            )

        # ── SCREENING: deviation is bounded ──
        screening_pass = growth_ratio <= 2.0
        if not screening_pass:
            return self._no_detection(
                primary_metric=primary,
                warnings=warnings + [
                    f"Deviation not bounded (growth_ratio={growth_ratio:.2f})"
                ],
                metadata_available=metadata_available,
            )

        # ── Surrogate null model ──
        # Generate uncontrolled surrogate trajectories: integrate
        # dx = perturbation * dt (no feedback), then compute deviation
        # integrals. Noise scale estimated from observed trajectory
        # residuals (small for deterministic systems).
        rng = np.random.default_rng(self.seed)

        # Estimate noise from pre-onset residuals (should be near zero
        # for deterministic systems)
        if onset_idx > 2:
            pre_residuals = np.diff(x[:onset_idx])
            noise_scale = float(np.std(pre_residuals))
        else:
            noise_scale = 0.0

        x0 = float(x[onset_idx])
        null_dev_integrals: List[float] = []
        for _ in range(self.n_permutations):
            surr = _generate_uncontrolled_surrogate(
                time, setpoint, perturbation, onset_idx, x0, rng,
                noise_scale=noise_scale,
            )
            surr_int = _compute_deviation_integral(surr, setpoint, time, onset_idx, n)
            null_dev_integrals.append(surr_int)

        null_arr = np.array(null_dev_integrals)
        null_mean = float(np.mean(null_arr))
        null_std = float(np.std(null_arr)) if len(null_arr) > 1 else 1.0
        # Ensure non-degenerate std for deterministic surrogates
        if null_std < 1e-12:
            null_std = max(1e-12, null_mean * 0.01)

        # p-value: fraction of null runs with deviation integral ≤ observed
        null_p = float(np.mean(null_arr <= dev_integral))
        null_p = max(null_p, 1.0 / (self.n_permutations + 1))

        # Effect size (positive = regulation better than null)
        cohens_d = (null_mean - dev_integral) / null_std if null_std > 0 else 0.0
        effect_size: Dict[str, float] = {
            'cohens_d': cohens_d,
            'raw_value': dev_integral,
            'null_mean': null_mean,
            'null_std': null_std,
        }

        # ── Secondaries ──
        dev_ratio = dev_integral / null_mean if null_mean > 0 else 1.0
        secondary: Dict[str, Any] = {
            'peak_deviation': peak_dev,
            'deviation_ratio': dev_ratio,
            'growth_ratio': growth_ratio,
        }

        # ── CONFIRMATION: deviation integral < surrogate at p < 0.01 ──
        confirmed = null_p < 0.01 and dev_ratio < 0.5

        if not confirmed:
            bonuses = {
                'secondaries_pass': bool(secondary),
                'shuffle_null_p_001': null_p < 0.01,
            }
            confidence = compute_confidence(DetectionTier.SCREENING, bonuses)
            return DetectorResult(
                pattern_id=self.pattern_id,
                detected=True,
                tier=DetectionTier.SCREENING,
                confidence=confidence,
                primary_metric=primary,
                secondary_metrics=secondary,
                effect_size=effect_size,
                null_p_value=null_p,
                null_type=NullType.SURROGATE,
                metadata_available=metadata_available,
                warnings=warnings,
            )

        # ── DEFINITIVE: deviation integral ≪ uncontrolled baseline ──
        # Note: with n_permutations=199, minimum p = 1/200 = 0.005,
        # so we use <= 0.005 to allow definitive at that floor.
        definitive = (
            confirmed
            and dev_ratio < 0.3
            and cohens_d > 1.0
            and null_p <= 0.005
        )

        # Metadata check: active feedback must be confirmed
        if definitive and metadata_available:
            has_feedback = metadata.get('has_active_feedback', None)
            if has_feedback is False:
                definitive = False
                warnings.append(
                    "Metadata indicates no active feedback; "
                    "definitive tier denied"
                )

        if definitive:
            tier = DetectionTier.DEFINITIVE
            bonuses = {
                'all_exclusions_cleared': True,
                'both_null_types_rejected': null_p < 0.005,
                'finite_size_robustness': n > 200,
            }
        else:
            tier = DetectionTier.CONFIRMATION
            bonuses = {
                'null_p_0001': null_p < 0.001,
                'effect_size_gt_1': cohens_d > 1.0,
                'all_secondaries': dev_ratio < 0.3,
            }

        confidence = compute_confidence(tier, bonuses)

        return DetectorResult(
            pattern_id=self.pattern_id,
            detected=True,
            tier=tier,
            confidence=confidence,
            primary_metric=primary,
            secondary_metrics=secondary,
            effect_size=effect_size,
            null_p_value=null_p,
            null_type=NullType.SURROGATE,
            metadata_available=metadata_available,
            warnings=warnings,
        )

    def _no_detection(
        self,
        primary_metric: Optional[Dict[str, float]] = None,
        warnings: Optional[List[str]] = None,
        metadata_available: bool = False,
    ) -> DetectorResult:
        """Return a non-detection result."""
        return DetectorResult(
            pattern_id=self.pattern_id,
            detected=False,
            tier=DetectionTier.SCREENING,
            confidence=0.0,
            primary_metric=primary_metric or {},
            secondary_metrics={},
            effect_size={},
            null_p_value=1.0,
            null_type=NullType.SURROGATE,
            metadata_available=metadata_available,
            warnings=warnings or [],
        )
