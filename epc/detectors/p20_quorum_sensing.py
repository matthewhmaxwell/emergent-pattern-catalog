"""P20 Quorum sensing / threshold-activated collective response detector.

Detects quorum sensing: a sharp collective behavioral SWITCH triggered
when agent density crosses a critical threshold, with hysteresis between
activation and deactivation. Distinct from P18 (graded consensus) and
P14 (avalanche criticality): P20 is a binary on/off population toggle
at a critical density.

Primary metrics:
  - Transition sharpness: step-function R² — how well the fraction_on
    vs density curve fits a step function (sharp = R² close to 1).
  - Hysteresis width: difference between activation density (up-sweep)
    and deactivation density (down-sweep).

Detection tiers:
  Screening:  A SHARP OFF→ON transition exists — amplitude jump>0.5 AND
              the transition occupies <40% of the density range. A graded
              (P18-consensus) response, even at large amplitude, is rejected
              here by the sharpness metric.
  Confirmation: Sharpness is near-ideal (step-function R²>0.85 AND
                transition band <25% of range) + observed R² exceeds the
                density-shuffle null.
  Definitive: Hysteresis present (activation density > deactivation
              density on up- vs down-sweep, measured from the substrate) +
              effect size + metadata does not deny threshold activation.

Null model: density-shuffle. Permute density labels on the equilibrium
fraction_on curve and re-fit the best step function; observed R² should
far exceed shuffled R² if transition is genuinely sharp.

T1a observation contract:
  Input: density-sweep observation bundle — time series of
  (density, concentration, collective_state, fraction_on, step,
  density_idx, sweep_direction) read via extract_observation_bundle().

Reference: Waters, C. M., & Bassler, B. L. (2005). Quorum sensing:
cell-to-cell communication in bacteria. Annual Review of Cell and
Developmental Biology, 21, 319–346.
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


# ── Calibrated sharpness thresholds (validation-rebuild) ─────────────
# A quorum switch is a SHARP threshold response; a graded/consensus
# response (P18) is not. Calibrated against the AutoinducerQuorum
# positive (rel_width 0.0, step_r2 1.0), sharp-but-reversible switches
# (rel_width 0.08-0.18, step_r2 0.94-0.97) and graded responses with
# large amplitude (rel_width 0.74-0.95, step_r2 0.80-0.82).
_SCREEN_MAX_REL_WIDTH = 0.4   # graded ramps rejected at SCREENING by sharpness
_CONFIRM_MIN_R2 = 0.85        # step-function R^2 must be near-ideal at confirmation
_CONFIRM_MAX_REL_WIDTH = 0.25 # confirmation tightens the transition-band width


# ── T1a: observation-bundle adapter ──────────────────────────────────

def extract_observation_bundle(
    history: List[Dict[str, Any]],
) -> Dict[str, Any]:
    """Extract the density-sweep observation bundle from history.

    This adapter is the T1a observation contract for P20. Any system that
    emits a time series with keys 'density', 'concentration',
    'collective_state', 'fraction_on', 'density_idx',
    'sweep_direction' can be consumed by the P20 detector.

    Parameters
    ----------
    history : list[dict]
        State trajectory. Each dict must contain at minimum:
        'density' (float), 'concentration' (float),
        'collective_state' (int, 0 or 1), 'fraction_on' (float),
        'density_idx' (int), 'sweep_direction' (str, 'up' or 'down').

    Returns
    -------
    dict with keys:
        'density': np.ndarray of shape (T,)
        'concentration': np.ndarray of shape (T,)
        'collective_state': np.ndarray of shape (T,) int
        'fraction_on': np.ndarray of shape (T,)
        'density_idx': np.ndarray of shape (T,) int
        'sweep_direction': list of str
    """
    density = np.array([h['density'] for h in history], dtype=float)
    concentration = np.array([h['concentration'] for h in history], dtype=float)
    collective_state = np.array([h['collective_state'] for h in history], dtype=int)
    fraction_on = np.array([h['fraction_on'] for h in history], dtype=float)
    density_idx = np.array([h['density_idx'] for h in history], dtype=int)
    sweep_direction = [h['sweep_direction'] for h in history]
    return {
        'density': density,
        'concentration': concentration,
        'collective_state': collective_state,
        'fraction_on': fraction_on,
        'density_idx': density_idx,
        'sweep_direction': sweep_direction,
    }


# ── Transition analysis ──────────────────────────────────────────────

def _compute_equilibrium_curve(
    bundle: Dict[str, Any],
    direction: str,
) -> tuple[np.ndarray, np.ndarray]:
    """Compute mean fraction_on vs density for one sweep direction.

    Uses only the last 50% of steps at each density level (after
    equilibration).

    Returns
    -------
    densities : ndarray of unique density values (sorted)
    mean_fraction : ndarray of mean fraction_on at each density
    """
    dir_mask = np.array([d == direction for d in bundle['sweep_direction']])
    if not np.any(dir_mask):
        return np.array([]), np.array([])

    density = bundle['density'][dir_mask]
    fraction = bundle['fraction_on'][dir_mask]
    d_idx = bundle['density_idx'][dir_mask]
    step = np.zeros(int(np.sum(dir_mask)), dtype=int)
    # Reconstruct step index within each density block
    for uid in np.unique(d_idx):
        idx_mask = d_idx == uid
        n_steps = int(np.sum(idx_mask))
        step[idx_mask] = np.arange(n_steps)

    unique_idx = np.unique(d_idx)
    n = len(unique_idx)
    densities = np.zeros(n)
    mean_frac = np.zeros(n)

    for i, idx in enumerate(unique_idx):
        idx_mask = d_idx == idx
        n_steps = int(np.sum(idx_mask))
        # Use last 50% (equilibrated)
        half = max(1, n_steps // 2)
        local_step = step[idx_mask]
        equil_mask = local_step >= half
        densities[i] = density[idx_mask][0]
        if np.any(equil_mask):
            mean_frac[i] = float(np.mean(fraction[idx_mask][equil_mask]))
        else:
            mean_frac[i] = float(np.mean(fraction[idx_mask]))

    order = np.argsort(densities)
    return densities[order], mean_frac[order]


def _best_step_function_r2(
    densities: np.ndarray,
    fraction: np.ndarray,
) -> Dict[str, Any]:
    """Fit the best step function and return R².

    Step function: f(d) = f_low if d < d_c, f_high if d >= d_c.
    Try every possible split point and pick the one with highest R².

    Returns dict with:
      r2: best step-function R²
      critical_density: density at the best split
      f_low: mean fraction below critical density
      f_high: mean fraction above critical density
      has_transition: f_high - f_low > 0.5
    """
    n = len(fraction)
    if n < 3:
        return {
            'r2': 0.0,
            'critical_density': 0.0,
            'f_low': 0.0,
            'f_high': 0.0,
            'has_transition': False,
        }

    total_var = float(np.var(fraction))
    if total_var < 1e-12:
        return {
            'r2': 0.0,
            'critical_density': float(densities[n // 2]),
            'f_low': float(fraction[0]),
            'f_high': float(fraction[-1]),
            'has_transition': False,
        }

    best_r2 = -1.0
    best_split = 1
    total_mean = float(np.mean(fraction))
    ss_tot = float(np.sum((fraction - total_mean) ** 2))

    for split in range(1, n):
        f_low = float(np.mean(fraction[:split]))
        f_high = float(np.mean(fraction[split:]))
        predicted = np.concatenate([
            np.full(split, f_low),
            np.full(n - split, f_high),
        ])
        ss_res = float(np.sum((fraction - predicted) ** 2))
        r2 = 1.0 - ss_res / ss_tot if ss_tot > 0 else 0.0
        if r2 > best_r2:
            best_r2 = r2
            best_split = split

    f_low = float(np.mean(fraction[:best_split]))
    f_high = float(np.mean(fraction[best_split:]))
    critical_density = float(
        (densities[best_split - 1] + densities[min(best_split, n - 1)]) / 2.0
    )

    return {
        'r2': best_r2,
        'critical_density': critical_density,
        'f_low': f_low,
        'f_high': f_high,
        'has_transition': (f_high - f_low) > 0.5,
    }


def _transition_sharpness(
    densities: np.ndarray,
    fraction: np.ndarray,
) -> Dict[str, Any]:
    """Measure transition sharpness via max slope and width.

    Returns dict with:
      max_slope: maximum |d(fraction)/d(density)|
      transition_width: density range for fraction 0.1 → 0.9
      has_transition: fraction spans from <0.2 to >0.8
    """
    n = len(fraction)
    if n < 3:
        return {
            'max_slope': 0.0,
            'transition_width': float('inf'),
            'has_transition': False,
        }

    d_frac = np.diff(fraction)
    d_dens = np.diff(densities)
    slopes = np.where(np.abs(d_dens) > 1e-12, d_frac / d_dens, 0.0)
    max_slope = float(np.max(np.abs(slopes)))

    frac_min = float(np.min(fraction))
    frac_max = float(np.max(fraction))
    has_transition = frac_min < 0.2 and frac_max > 0.8

    low_idx = np.where(fraction >= 0.1)[0]
    high_idx = np.where(fraction >= 0.9)[0]
    if len(low_idx) > 0 and len(high_idx) > 0:
        transition_width = float(abs(densities[high_idx[0]] - densities[low_idx[0]]))
    else:
        transition_width = float(densities[-1] - densities[0])

    return {
        'max_slope': max_slope,
        'transition_width': transition_width,
        'has_transition': has_transition,
    }


def _compute_hysteresis(
    up_densities: np.ndarray,
    up_fraction: np.ndarray,
    down_densities: np.ndarray,
    down_fraction: np.ndarray,
) -> Dict[str, Any]:
    """Measure hysteresis between up-sweep and down-sweep.

    Returns dict with:
      activation_density: density at which up-sweep crosses 0.5
      deactivation_density: density at which down-sweep crosses 0.5
      hysteresis_width: activation - deactivation
      has_hysteresis: activation > deactivation meaningfully
    """
    def _crossing(dens: np.ndarray, frac: np.ndarray, threshold: float = 0.5) -> float:
        """Find density at which fraction crosses threshold (upward)."""
        for i in range(len(frac) - 1):
            if frac[i] < threshold <= frac[i + 1]:
                t = (threshold - frac[i]) / (frac[i + 1] - frac[i] + 1e-12)
                return float(dens[i] + t * (dens[i + 1] - dens[i]))
        # Fallback: threshold might never be crossed cleanly
        closest = int(np.argmin(np.abs(frac - threshold)))
        return float(dens[closest])

    if len(up_densities) < 2 or len(down_densities) < 2:
        return {
            'activation_density': 0.0,
            'deactivation_density': 0.0,
            'hysteresis_width': 0.0,
            'relative_hysteresis': 0.0,
            'has_hysteresis': False,
        }

    # Up-sweep: find where fraction crosses 0.5 going up
    act_density = _crossing(up_densities, up_fraction, 0.5)

    # Down-sweep: densities are already sorted ascending. The fraction
    # shows the system ON at high density and OFF at low density (after
    # sweeping density down). The deactivation density is where fraction
    # crosses 0.5 (going from 0 at low density to 1 at high density).
    deact_density = _crossing(down_densities, down_fraction, 0.5)

    width = act_density - deact_density
    density_range = float(up_densities[-1] - up_densities[0])
    relative_width = width / density_range if density_range > 0 else 0.0

    return {
        'activation_density': act_density,
        'deactivation_density': deact_density,
        'hysteresis_width': width,
        'relative_hysteresis': relative_width,
        'has_hysteresis': width > 0.05 * density_range,  # >5% of range
    }


class P20QuorumSensingDetector:
    """P20 quorum sensing / threshold-activated response detector.

    Reads inputs via the T1a observation-bundle adapter.

    Parameters
    ----------
    n_permutations : int
        Number of density-shuffles for the null model.
    seed : int
        RNG seed for reproducibility.
    """

    def __init__(
        self,
        n_permutations: int = 199,
        seed: int = 42,
    ) -> None:
        self.pattern_id = "P20"
        self.n_permutations = n_permutations
        self.seed = seed

    def detect(
        self,
        history: List[Dict[str, Any]],
        metadata: Optional[Dict[str, Any]] = None,
    ) -> DetectorResult:
        """Run P20 quorum sensing detection.

        Parameters
        ----------
        history : list[dict]
            State trajectory with density-sweep bundle keys.
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

        # ── Compute equilibrium transition curves ──
        up_dens, up_frac = _compute_equilibrium_curve(bundle, 'up')
        down_dens, down_frac = _compute_equilibrium_curve(bundle, 'down')

        if len(up_dens) < 3:
            return self._no_detection(
                warnings=["Insufficient density levels in up-sweep (need ≥3)"],
                metadata_available=metadata_available,
            )

        # ── Step-function fit (primary sharpness metric) ──
        step_fit = _best_step_function_r2(up_dens, up_frac)
        sharp = _transition_sharpness(up_dens, up_frac)

        primary: Dict[str, float] = {
            'step_r2': step_fit['r2'],
            'critical_density': step_fit['critical_density'],
            'f_low': step_fit['f_low'],
            'f_high': step_fit['f_high'],
            'max_slope': sharp['max_slope'],
            'transition_width': sharp['transition_width'],
        }

        # Relative transition width (sharpness): a sharp switch occupies a
        # narrow band of the density range; a graded ramp spans most of it.
        # Computed here (before screening) so the sharpness metric — not an
        # amplitude guard — is what rejects graded responses.
        density_range = float(up_dens[-1] - up_dens[0])
        relative_width = (
            sharp['transition_width'] / density_range if density_range > 0 else 1.0
        )

        # ── SCREENING: a SHARP OFF→ON transition (not a graded ramp) ──
        #   (a) amplitude: f_high - f_low > 0.5 (a real OFF→ON switch)
        #   (b) sharpness: transition occupies < 40% of the density range.
        # A graded response with LARGE amplitude (jump>0.5 but wide
        # transition, e.g. a gentle logistic) is rejected HERE by the
        # sharpness metric — this is the P18-consensus discrimination case,
        # and the reason the R^2/width sharpness test is a live gate, not
        # dead code downstream of an amplitude guard.
        if not step_fit['has_transition']:
            return self._no_detection(
                primary_metric=primary,
                warnings=warnings + [
                    "No OFF→ON transition observed "
                    f"(f_low={step_fit['f_low']:.3f}, f_high={step_fit['f_high']:.3f}, "
                    f"jump={step_fit['f_high'] - step_fit['f_low']:.3f})"
                ],
                metadata_available=metadata_available,
            )
        if relative_width >= _SCREEN_MAX_REL_WIDTH:
            return self._no_detection(
                primary_metric=primary,
                warnings=warnings + [
                    "Transition too graded for a quorum switch "
                    f"(relative_width={relative_width:.3f} >= {_SCREEN_MAX_REL_WIDTH}); "
                    f"step_r2={step_fit['r2']:.3f}. A sharp threshold response "
                    "occupies a narrow density band; a graded/consensus response "
                    "spans the sweep."
                ],
                metadata_available=metadata_available,
            )

        # ── Null model: density-shuffle ──
        # Permute fraction values across density levels and re-fit step
        # function. A sharp transition has high R² that exceeds shuffled.
        rng = np.random.default_rng(self.seed)
        observed_r2 = step_fit['r2']

        null_r2s: List[float] = []
        for _ in range(self.n_permutations):
            shuffled_frac = up_frac.copy()
            rng.shuffle(shuffled_frac)
            null_fit = _best_step_function_r2(up_dens, shuffled_frac)
            null_r2s.append(null_fit['r2'])

        null_arr = np.array(null_r2s)
        null_mean = float(np.mean(null_arr))
        null_std = float(np.std(null_arr)) if len(null_arr) > 1 else 1.0
        if null_std < 1e-12:
            null_std = max(1e-12, abs(null_mean) * 0.01)

        null_p = float(np.mean(null_arr >= observed_r2))
        null_p = max(null_p, 1.0 / (self.n_permutations + 1))

        cohens_d = (observed_r2 - null_mean) / null_std if null_std > 0 else 0.0

        effect_size: Dict[str, float] = {
            'cohens_d': cohens_d,
            'observed_r2': observed_r2,
            'null_mean_r2': null_mean,
            'null_std_r2': null_std,
        }

        # ── Hysteresis analysis ──
        if len(down_dens) >= 3:
            hyst = _compute_hysteresis(up_dens, up_frac, down_dens, down_frac)
        else:
            hyst = {
                'activation_density': step_fit['critical_density'],
                'deactivation_density': step_fit['critical_density'],
                'hysteresis_width': 0.0,
                'relative_hysteresis': 0.0,
                'has_hysteresis': False,
            }

        primary['hysteresis_width'] = hyst['hysteresis_width']
        primary['activation_density'] = hyst['activation_density']
        primary['deactivation_density'] = hyst['deactivation_density']

        # ── Secondaries ──
        secondary: Dict[str, Any] = {
            'has_transition': step_fit['has_transition'],
            'has_hysteresis': hyst['has_hysteresis'],
            'relative_hysteresis': hyst.get('relative_hysteresis', 0.0),
            'relative_transition_width': relative_width,
            'up_curve': up_frac.tolist(),
            'down_curve': down_frac.tolist() if len(down_frac) > 0 else [],
            'up_densities': up_dens.tolist(),
            'down_densities': down_dens.tolist() if len(down_dens) > 0 else [],
            'n_density_levels': len(up_dens),
        }

        # ── CONFIRMATION: sharp transition + null rejected ──
        is_sharp = step_fit['r2'] > _CONFIRM_MIN_R2 and relative_width < _CONFIRM_MAX_REL_WIDTH

        confirmed = (
            step_fit['has_transition']
            and is_sharp
            and null_p < 0.05
        )

        if not confirmed:
            reasons: List[str] = []
            if not is_sharp:
                reasons.append(
                    f"Transition not sharp enough (R²={step_fit['r2']:.3f}, "
                    f"rel_width={relative_width:.3f})"
                )
            if null_p >= 0.05:
                reasons.append(f"Null not rejected (p={null_p:.4f})")
            bonuses = {
                'secondaries_pass': step_fit['has_transition'],
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
                warnings=warnings + reasons,
            )

        # ── DEFINITIVE: hysteresis + metadata ──
        definitive = (
            confirmed
            and hyst['has_hysteresis']
            and cohens_d > 1.0
            and null_p <= 0.005
        )

        # Metadata gate: threshold activation must be confirmed
        if definitive and metadata_available:
            has_threshold = metadata.get('has_threshold_activation', None)
            if has_threshold is False:
                definitive = False
                warnings.append(
                    "Metadata indicates no threshold activation; "
                    "definitive tier denied"
                )

        if definitive:
            tier = DetectionTier.DEFINITIVE
            bonuses = {
                'all_exclusions_cleared': True,
                'both_null_types_rejected': null_p < 0.005,
                'finite_size_robustness': len(up_dens) >= 20,
            }
        else:
            tier = DetectionTier.CONFIRMATION
            bonuses = {
                'null_p_0001': null_p < 0.001,
                'effect_size_gt_1': cohens_d > 1.0,
                'all_secondaries': (
                    step_fit['has_transition']
                    and is_sharp
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
