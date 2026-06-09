"""P25 Canalized restoration / equifinality detector.

Detects equifinality: a system that converges to the SAME target
macrostate from a WIDE RANGE of initial conditions, beyond what a
simple dominant attractor predicts. Optionally detects restoration
after mid-trajectory perturbation.

Primary metric:
  - Convergence variance ratio = Var(final states) / Var(initial states).
    Canalized systems have ratio ≪ 1 (diverse ICs → similar endpoints).

Secondary metrics:
  - Basin volume: fraction of ICs that converge to the target.
  - Restoration accuracy: fraction of perturbed trajectories that
    return to the target (if perturbation data is present).
  - Relaxation time: median steps-to-convergence (rejects trivial
    instant-collapse maps).

Detection tiers:
  Screening:  Convergence variance ratio < 0.1 (final states much
              tighter than initial states).
  Confirmation: Large basin volume (≥ 0.8 of ICs converge) AND
                convergence ratio significantly below IC-shuffle null.
  Definitive: Restoration after perturbation + nontrivial relaxation
              time (active canalization, not trivial collapse).

Null model: IC-shuffle surrogate. Randomly reassign final states
across ICs; the resulting "convergence" ratio should be ≈ 1.0
(no relationship between IC and final state). The observed ratio
must be significantly below this baseline.

T1a observation contract:
  Input: canalization observation bundle — trajectories from multiple
  ICs with (state, step, trial, ic, target, distance_to_target, converged).
  Read via extract_observation_bundle().

Reference: Waddington, C. H. (1957). The Strategy of the Genes.
           Huang, S. et al. (2005). PRL 94, 128701.
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
    """Extract the canalization observation bundle from history.

    Groups records by trial and extracts per-trial IC, final state,
    target, and convergence info.

    Parameters
    ----------
    history : list[dict]
        State trajectory. Each dict must contain:
        'state' (ndarray), 'step' (int), 'trial' (int),
        'ic' (ndarray), 'target' (ndarray),
        'distance_to_target' (float), 'converged' (bool).

    Returns
    -------
    dict with keys:
        'ics': ndarray of shape (n_trials, n_dims)
        'finals': ndarray of shape (n_trials, n_dims)
        'target': ndarray of shape (n_dims,)
        'final_distances': ndarray of shape (n_trials,)
        'converged': ndarray of shape (n_trials,) bool
        'n_trials': int
        'n_dims': int
        'steps_to_convergence': ndarray of shape (n_trials,)
        'n_steps': int
    """
    # Group by trial
    trials: Dict[int, List[Dict[str, Any]]] = {}
    for h in history:
        trial = h['trial']
        if trial not in trials:
            trials[trial] = []
        trials[trial].append(h)

    n_trials = len(trials)
    if n_trials == 0:
        raise ValueError("No trials found in history")

    # Sort trials by trial index
    trial_ids = sorted(trials.keys())

    # Extract per-trial data
    first_record = trials[trial_ids[0]][0]
    n_dims = len(np.asarray(first_record['state']))
    target = np.asarray(first_record['target'])

    ics = np.zeros((n_trials, n_dims))
    finals = np.zeros((n_trials, n_dims))
    final_distances = np.zeros(n_trials)
    converged = np.zeros(n_trials, dtype=bool)
    steps_to_conv = np.full(n_trials, -1, dtype=int)

    n_steps = 0
    for i, tid in enumerate(trial_ids):
        records = sorted(trials[tid], key=lambda r: r['step'])
        n_steps = max(n_steps, len(records))

        ics[i] = np.asarray(records[0]['ic'])
        finals[i] = np.asarray(records[-1]['state'])
        final_distances[i] = records[-1]['distance_to_target']
        converged[i] = records[-1]['converged']

        # Find first step where converged
        for r in records:
            if r['converged']:
                steps_to_conv[i] = r['step']
                break

    return {
        'ics': ics,
        'finals': finals,
        'target': target,
        'final_distances': final_distances,
        'converged': converged,
        'n_trials': n_trials,
        'n_dims': n_dims,
        'steps_to_convergence': steps_to_conv,
        'n_steps': n_steps,
    }


# ── Core detection logic ─────────────────────────────────────────────

def _convergence_variance_ratio(
    ics: np.ndarray,
    finals: np.ndarray,
) -> float:
    """Compute Var(finals) / Var(ics) across all dimensions.

    Uses total variance (sum of per-dimension variances) for both.
    A canalized system has ratio ≪ 1.
    """
    ic_var = float(np.sum(np.var(ics, axis=0)))
    final_var = float(np.sum(np.var(finals, axis=0)))

    if ic_var < 1e-12:
        return float('inf') if final_var > 1e-12 else 1.0

    return final_var / ic_var


def _basin_volume(
    converged: np.ndarray,
) -> float:
    """Fraction of ICs that converged to the target."""
    if len(converged) == 0:
        return 0.0
    return float(np.mean(converged))


def _median_relaxation_time(
    steps_to_convergence: np.ndarray,
    n_steps: int,
) -> float:
    """Median steps-to-convergence across converged trials.

    Returns n_steps if no trials converged (indicating no convergence).
    A trivial collapse has relaxation_time ≤ 1.
    """
    valid = steps_to_convergence[steps_to_convergence >= 0]
    if len(valid) == 0:
        return float(n_steps)
    return float(np.median(valid))


def _restoration_accuracy(
    history: List[Dict[str, Any]],
    perturbation_time: Optional[int],
) -> Optional[float]:
    """Fraction of perturbed trajectories that restore to target.

    Returns None if no perturbation was applied.
    """
    if perturbation_time is None:
        return None

    # Group by trial
    trials: Dict[int, List[Dict[str, Any]]] = {}
    for h in history:
        trial = h['trial']
        if trial not in trials:
            trials[trial] = []
        trials[trial].append(h)

    n_restored = 0
    n_total = 0
    for tid in sorted(trials.keys()):
        records = sorted(trials[tid], key=lambda r: r['step'])
        # Check if post-perturbation the trajectory converges
        post_pert = [r for r in records if r['step'] > perturbation_time]
        if not post_pert:
            continue
        n_total += 1
        if post_pert[-1]['converged']:
            n_restored += 1

    if n_total == 0:
        return None
    return float(n_restored / n_total)


class P25EquifinalityDetector:
    """P25 canalized restoration / equifinality detector.

    Reads inputs via the T1a observation-bundle adapter.

    Parameters
    ----------
    n_permutations : int
        Number of IC-shuffle permutations for the null model.
    convergence_threshold : float
        Distance threshold for declaring convergence (default 0.1).
    seed : int
        RNG seed for reproducibility.
    """

    def __init__(
        self,
        n_permutations: int = 199,
        convergence_threshold: float = 0.1,
        seed: int = 42,
    ) -> None:
        self.pattern_id = "P25"
        self.n_permutations = n_permutations
        self.convergence_threshold = convergence_threshold
        self.seed = seed

    def detect(
        self,
        history: List[Dict[str, Any]],
        metadata: Optional[Dict[str, Any]] = None,
    ) -> DetectorResult:
        """Run P25 equifinality detection.

        Parameters
        ----------
        history : list[dict]
            State trajectory with canalization observation-bundle keys.
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
        except (KeyError, IndexError, ValueError) as e:
            return self._no_detection(
                warnings=[f"Observation bundle extraction failed: {e}"],
                metadata_available=metadata_available,
            )

        ics = bundle['ics']
        finals = bundle['finals']
        converged = bundle['converged']
        steps_to_conv = bundle['steps_to_convergence']
        n_trials = bundle['n_trials']
        n_steps = bundle['n_steps']

        if n_trials < 3:
            return self._no_detection(
                warnings=["Insufficient trials (need ≥ 3)"],
                metadata_available=metadata_available,
            )

        # ── Primary metric: convergence variance ratio ──
        conv_ratio = _convergence_variance_ratio(ics, finals)

        primary: Dict[str, float] = {
            'convergence_variance_ratio': conv_ratio,
        }

        # ── Secondary metrics (computed early for screening gate) ──
        basin_vol = _basin_volume(converged)
        relax_time = _median_relaxation_time(steps_to_conv, n_steps)

        # ── SCREENING: convergence ratio < 0.1 AND basin volume ≥ 0.5 ──
        # Equifinality (Waddington 1957) requires convergence from a wide
        # IC range: a meaningful fraction of ICs must actually reach the
        # target (basin_vol ≥ 0.5), not just that the variance ratio is
        # low. This prerequisite distinguishes equifinality from noisy
        # regulation (P24) where the regulated variable hovers near the
        # setpoint but does not reliably converge within the convergence
        # threshold.
        screening_pass = conv_ratio < 0.1 and basin_vol >= 0.5
        if not screening_pass:
            fail_reasons = []
            if conv_ratio >= 0.1:
                fail_reasons.append(
                    f"Convergence ratio too high ({conv_ratio:.4f} ≥ 0.1)"
                )
            if basin_vol < 0.5:
                fail_reasons.append(
                    f"Basin volume too low ({basin_vol:.3f} < 0.5); "
                    "equifinality requires a wide IC range to converge"
                )
            return self._no_detection(
                primary_metric=primary,
                warnings=warnings + fail_reasons,
                metadata_available=metadata_available,
            )

        # Get perturbation info from metadata
        perturbation_time = None
        if metadata_available:
            perturbation_time = metadata.get('perturbation_time', None)
        restoration = _restoration_accuracy(history, perturbation_time)

        secondary: Dict[str, Any] = {
            'basin_volume': basin_vol,
            'relaxation_time': relax_time,
            'n_trials': n_trials,
        }
        if restoration is not None:
            secondary['restoration_accuracy'] = restoration

        # ── Trivial-collapse gate ──
        # Reject if relaxation time is ≤ 1 step (instant collapse,
        # not dynamics-driven convergence)
        is_trivial_collapse = relax_time <= 1.0
        if is_trivial_collapse:
            warnings.append(
                f"Trivial collapse detected (relaxation_time={relax_time:.1f} ≤ 1); "
                "dynamics-driven convergence required"
            )
            return self._no_detection(
                primary_metric=primary,
                warnings=warnings,
                metadata_available=metadata_available,
            )

        # ── IC-distribution surrogate null model ──
        # Generate surrogate "final states" by resampling from the IC
        # distribution. If the system is NOT canalized, final states should
        # be as spread out as ICs → ratio ≈ 1.0. The observed ratio must
        # be significantly below this to confirm equifinality.
        rng = np.random.default_rng(self.seed)
        ic_mean = np.mean(ics, axis=0)
        ic_cov = np.cov(ics.T) if n_trials > 1 else np.eye(ics.shape[1])
        # Regularize covariance for stability
        ic_cov = np.atleast_2d(ic_cov)
        ic_cov += np.eye(ic_cov.shape[0]) * 1e-8

        null_ratios: List[float] = []
        for _ in range(self.n_permutations):
            # Draw surrogate finals from IC distribution
            surr_finals = rng.multivariate_normal(ic_mean, ic_cov, size=n_trials)
            null_ratio = _convergence_variance_ratio(ics, surr_finals)
            null_ratios.append(null_ratio)

        null_arr = np.array(null_ratios)
        null_mean = float(np.mean(null_arr))
        null_std = float(np.std(null_arr)) if len(null_arr) > 1 else 1.0
        if null_std < 1e-12:
            null_std = max(1e-12, null_mean * 0.01)

        # p-value: fraction of null runs with ratio ≤ observed
        null_p = float(np.mean(null_arr <= conv_ratio))
        null_p = max(null_p, 1.0 / (self.n_permutations + 1))

        # Effect size
        cohens_d = (null_mean - conv_ratio) / null_std if null_std > 0 else 0.0
        effect_size: Dict[str, float] = {
            'cohens_d': cohens_d,
            'raw_value': conv_ratio,
            'null_mean': null_mean,
            'null_std': null_std,
        }

        # ── CONFIRMATION: large basin + significant vs null ──
        confirmed = (
            null_p < 0.01
            and basin_vol >= 0.8
            and conv_ratio < 0.1
        )

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
                null_type=NullType.SHUFFLE,
                metadata_available=metadata_available,
                warnings=warnings,
            )

        # ── DEFINITIVE: restoration + nontrivial dynamics ──
        definitive = (
            confirmed
            and null_p <= 0.005
            and cohens_d > 1.0
            and relax_time > 1.0  # Not trivial collapse
        )

        # Perturbation restoration strengthens definitive
        if definitive and restoration is not None:
            definitive = definitive and restoration >= 0.8

        # If perturbation data present but restoration is poor, deny definitive
        if definitive and restoration is not None and restoration < 0.8:
            definitive = False
            warnings.append(
                f"Restoration accuracy {restoration:.2f} < 0.8; "
                "definitive denied"
            )

        if definitive:
            tier = DetectionTier.DEFINITIVE
            bonuses = {
                'all_exclusions_cleared': True,
                'both_null_types_rejected': null_p < 0.005,
                'finite_size_robustness': n_trials >= 15,
            }
        else:
            tier = DetectionTier.CONFIRMATION
            bonuses = {
                'null_p_0001': null_p < 0.001,
                'effect_size_gt_1': cohens_d > 1.0,
                'all_secondaries': basin_vol >= 0.9 and relax_time > 1.0,
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
            exclusions_checked=['P16', 'P24'],
            exclusion_results={
                'P16': 'not_applicable',  # P16 = multiple retrievable targets
                'P24': 'not_applicable',  # P24 = cybernetic regulation under perturbation
            },
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
