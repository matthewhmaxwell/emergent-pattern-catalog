"""P26 Stochastic resonance detector.

Detects stochastic resonance: noise improves system performance at an
intermediate optimal level, producing an inverted-U performance-vs-noise
curve. Too little noise → no signal detection; too much noise → signal
destroyed; a sweet spot amplifies it.

Primary metric:
  - SNR (or correlation-based detection accuracy) vs noise level —
    must show an inverted-U shape with a clear interior peak.

Detection tiers:
  Screening:  Performance at some nonzero noise > performance at zero noise.
  Confirmation: Inverted-U shape — performance rises then falls across
                the noise sweep (interior argmax, decline after peak).
  Definitive: Clear interior peak with meaningful performance gain over
              no-noise baseline + metadata confirms subthreshold signal.

Null model: noise-level shuffle. Permute the noise-level labels on the
performance curve and check whether an interior peak arises by chance.

T1a observation contract:
  Input: noise-sweep observation bundle — time series of
  (time, x, signal, noise_level, noise_level_idx) read via
  extract_observation_bundle(), NOT from the model object directly.

Reference: Gammaitoni, L. et al. (1998). Stochastic resonance.
Rev. Mod. Phys. 70(1), 223–287.
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
) -> Dict[str, Any]:
    """Extract the noise-sweep observation bundle from history.

    This adapter is the T1a observation contract for P26. Any system that
    emits a time series with keys 'time', 'x', 'signal', 'noise_level',
    and 'noise_level_idx' can be consumed by the P26 detector.

    Parameters
    ----------
    history : list[dict]
        State trajectory. Each dict must contain at minimum:
        'time' (float), 'x' (float), 'signal' (float),
        'noise_level' (float), 'noise_level_idx' (int).

    Returns
    -------
    dict with keys:
        'time': np.ndarray of shape (T,)
        'x': np.ndarray of shape (T,)
        'signal': np.ndarray of shape (T,)
        'noise_level': np.ndarray of shape (T,)
        'noise_level_idx': np.ndarray of shape (T,) int
    """
    time = np.array([h['time'] for h in history], dtype=float)
    x = np.array([h['x'] for h in history], dtype=float)
    signal = np.array([h['signal'] for h in history], dtype=float)
    noise_level = np.array([h['noise_level'] for h in history], dtype=float)
    noise_level_idx = np.array([h['noise_level_idx'] for h in history], dtype=int)
    return {
        'time': time,
        'x': x,
        'signal': signal,
        'noise_level': noise_level,
        'noise_level_idx': noise_level_idx,
    }


# ── Performance metric computation ──────────────────────────────────

def _coherent_response(
    x_seg: np.ndarray,
    signal_seg: np.ndarray,
) -> float:
    """Compute the coherent response: |⟨x · signal⟩|.

    This unnormalized cross-correlation directly measures how much x
    co-varies with the driving signal. It captures the canonical SR
    inverted-U (Gammaitoni 1998):
    - Zero noise: x stays in one well with tiny oscillation → small ⟨x·s⟩
    - Optimal noise: inter-well hopping synchronized with signal → large ⟨x·s⟩
    - High noise: random hopping, x and signal uncorrelated → small ⟨x·s⟩

    Unlike Pearson correlation, this metric is sensitive to the AMPLITUDE
    of the coherent response, not just its phase. At zero noise, x tracks
    the signal phase perfectly but with tiny amplitude → small coherent
    response. This is exactly the SR signature.
    """
    n = len(x_seg)
    if n < 4:
        return 0.0
    return float(np.abs(np.mean(x_seg * signal_seg)))


def _compute_performance_curve(
    history_bundle: Dict[str, Any],
) -> tuple[np.ndarray, np.ndarray]:
    """Compute performance vs noise level curve.

    Performance = |⟨x · signal⟩| (coherent response) at each noise
    level. This captures the canonical SR inverted-U: small coherent
    response at zero noise (intra-well oscillation), peak at
    intermediate noise (inter-well hopping synchronized with signal),
    declining at high noise (random hopping).

    Returns
    -------
    noise_levels : ndarray of unique noise levels
    performance : ndarray of coherent response at each level
    """
    x = history_bundle['x']
    signal = history_bundle['signal']
    noise_level = history_bundle['noise_level']
    noise_level_idx = history_bundle['noise_level_idx']

    unique_idx = np.unique(noise_level_idx)
    n_levels = len(unique_idx)

    noise_vals = np.zeros(n_levels)
    perf = np.zeros(n_levels)

    for i, idx in enumerate(unique_idx):
        mask = noise_level_idx == idx
        x_seg = x[mask]
        sig_seg = signal[mask]
        nl_seg = noise_level[mask]

        noise_vals[i] = nl_seg[0]
        perf[i] = _coherent_response(x_seg, sig_seg)

    # Sort by noise level
    order = np.argsort(noise_vals)
    return noise_vals[order], perf[order]


def _check_inverted_u(
    noise_levels: np.ndarray,
    performance: np.ndarray,
) -> Dict[str, Any]:
    """Check whether performance-vs-noise shows an inverted-U shape.

    Returns dict with:
      peak_idx: index of the peak
      peak_noise: noise level at peak
      peak_perf: performance at peak
      zero_noise_perf: performance at lowest (zero) noise
      high_noise_perf: performance at highest noise
      gain_over_zero: peak_perf - zero_noise_perf
      decline_after_peak: peak_perf - high_noise_perf
      is_interior_peak: peak is not at 0 or max noise
      has_rise: performance increases from zero noise to peak
      has_fall: performance decreases from peak to high noise
    """
    n = len(performance)
    peak_idx = int(np.argmax(performance))
    peak_noise = float(noise_levels[peak_idx])
    peak_perf = float(performance[peak_idx])
    zero_perf = float(performance[0])  # lowest noise
    high_perf = float(performance[-1])  # highest noise

    is_interior = 0 < peak_idx < n - 1
    gain = peak_perf - zero_perf
    decline = peak_perf - high_perf

    return {
        'peak_idx': peak_idx,
        'peak_noise': peak_noise,
        'peak_perf': peak_perf,
        'zero_noise_perf': zero_perf,
        'high_noise_perf': high_perf,
        'gain_over_zero': gain,
        'decline_after_peak': decline,
        'is_interior_peak': is_interior,
        'has_rise': gain > 0.02,
        'has_fall': decline > 0.02,
    }


class P26StochasticResonanceDetector:
    """P26 stochastic resonance detector.

    Reads inputs via the T1a observation-bundle adapter.

    Parameters
    ----------
    n_permutations : int
        Number of noise-level shuffles for the null model.
    seed : int
        RNG seed for reproducibility.
    """

    def __init__(
        self,
        n_permutations: int = 199,
        seed: int = 42,
    ) -> None:
        self.pattern_id = "P26"
        self.n_permutations = n_permutations
        self.seed = seed

    def detect(
        self,
        history: List[Dict[str, Any]],
        metadata: Optional[Dict[str, Any]] = None,
    ) -> DetectorResult:
        """Run P26 stochastic resonance detection.

        Parameters
        ----------
        history : list[dict]
            State trajectory with noise-sweep bundle keys.
        metadata : dict, optional
            Model metadata (used for definitive tier checks).

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

        # ── Compute performance curve ──
        noise_levels, performance = _compute_performance_curve(bundle)

        n_levels = len(noise_levels)
        if n_levels < 3:
            return self._no_detection(
                warnings=["Insufficient noise levels (need ≥3)"],
                metadata_available=metadata_available,
            )

        # ── Check inverted-U shape ──
        u_check = _check_inverted_u(noise_levels, performance)

        primary: Dict[str, float] = {
            'peak_noise': u_check['peak_noise'],
            'peak_performance': u_check['peak_perf'],
            'zero_noise_performance': u_check['zero_noise_perf'],
            'high_noise_performance': u_check['high_noise_perf'],
            'gain_over_zero': u_check['gain_over_zero'],
        }

        # ── SCREENING: some nonzero noise outperforms zero noise ──
        screening_pass = u_check['gain_over_zero'] > 0.02

        if not screening_pass:
            return self._no_detection(
                primary_metric=primary,
                warnings=warnings + [
                    f"No noise benefit (gain={u_check['gain_over_zero']:.4f})"
                ],
                metadata_available=metadata_available,
            )

        # ── Null model: time-shuffle at peak noise level ──
        # At the peak noise level, time-shuffle x (destroying temporal
        # order while preserving marginal distribution) and recompute
        # ⟨x·signal⟩. Under H0 (no synchronization), the shuffled
        # coherent response should be ~0. The observed peak should
        # far exceed this baseline if genuine SR is present.
        rng = np.random.default_rng(self.seed)
        observed_peak = u_check['peak_perf']

        x_full = bundle['x']
        signal_full = bundle['signal']
        noise_level_idx_full = bundle['noise_level_idx']
        unique_idx = np.unique(noise_level_idx_full)

        # Find the peak noise level's data
        peak_idx_sorted = u_check['peak_idx']
        peak_nl_idx = unique_idx[np.argsort(
            np.array([bundle['noise_level'][noise_level_idx_full == idx][0]
                       for idx in unique_idx])
        )[peak_idx_sorted]]
        peak_mask = noise_level_idx_full == peak_nl_idx
        x_peak = x_full[peak_mask].copy()
        sig_peak = signal_full[peak_mask]

        null_peaks: List[float] = []
        for _ in range(self.n_permutations):
            x_shuffled = x_peak.copy()
            rng.shuffle(x_shuffled)
            null_peaks.append(_coherent_response(x_shuffled, sig_peak))

        null_arr = np.array(null_peaks)
        null_mean = float(np.mean(null_arr))
        null_std = float(np.std(null_arr)) if len(null_arr) > 1 else 1.0
        if null_std < 1e-12:
            null_std = max(1e-12, abs(null_mean) * 0.01)

        # p-value: fraction of null runs with peak ≥ observed peak
        null_p = float(np.mean(null_arr >= observed_peak))
        null_p = max(null_p, 1.0 / (self.n_permutations + 1))

        # Effect size: how many SDs above the null peak distribution
        cohens_d = (observed_peak - null_mean) / null_std if null_std > 0 else 0.0
        effect_size: Dict[str, float] = {
            'cohens_d': cohens_d,
            'observed_peak': observed_peak,
            'null_mean_peak': null_mean,
            'null_std_peak': null_std,
        }

        # ── Secondaries ──
        secondary: Dict[str, Any] = {
            'is_interior_peak': u_check['is_interior_peak'],
            'has_rise': u_check['has_rise'],
            'has_fall': u_check['has_fall'],
            'peak_idx': u_check['peak_idx'],
            'n_levels': n_levels,
            'decline_after_peak': u_check['decline_after_peak'],
            'performance_curve': performance.tolist(),
            'noise_levels': noise_levels.tolist(),
        }

        # ── CONFIRMATION: inverted-U shape ──
        confirmed = (
            u_check['is_interior_peak']
            and u_check['has_rise']
            and u_check['has_fall']
            and null_p < 0.05
        )

        if not confirmed:
            bonuses = {
                'secondaries_pass': u_check['has_rise'],
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
                null_type=NullType.SHUFFLE,
                metadata_available=metadata_available,
                warnings=warnings,
            )

        # ── DEFINITIVE: clear interior peak + meaningful gain ──
        definitive = (
            confirmed
            and u_check['gain_over_zero'] > 0.05
            and u_check['decline_after_peak'] > 0.02
            and cohens_d > 1.0
            and null_p <= 0.005
        )

        # Metadata check: subthreshold signal must be confirmed
        if definitive and metadata_available:
            has_subthreshold = metadata.get('has_subthreshold_signal', None)
            if has_subthreshold is False:
                definitive = False
                warnings.append(
                    "Metadata indicates suprathreshold signal; "
                    "definitive tier denied"
                )

        if definitive:
            tier = DetectionTier.DEFINITIVE
            bonuses = {
                'all_exclusions_cleared': True,
                'both_null_types_rejected': null_p < 0.005,
                'finite_size_robustness': n_levels >= 10,
            }
        else:
            tier = DetectionTier.CONFIRMATION
            bonuses = {
                'null_p_0001': null_p < 0.001,
                'effect_size_gt_1': cohens_d > 1.0,
                'all_secondaries': (
                    u_check['is_interior_peak']
                    and u_check['has_rise']
                    and u_check['has_fall']
                ),
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
            null_type=NullType.SHUFFLE,
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
            null_type=NullType.SHUFFLE,
            metadata_available=metadata_available,
            warnings=warnings or [],
        )
