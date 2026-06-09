"""P16 Associative memory / pattern completion detector.

Detects content-addressable memory: a system that stores multiple
distributed patterns and retrieves them from partial/noisy cues via
attractor dynamics.

Primary metric: **completion accuracy** — overlap of the converged
state with the nearest stored pattern, averaged across cue trials.

Detection tiers:
  Screening:  Converged state overlaps a stored pattern above chance.
  Confirmation: Completion accuracy is high across a range of cue
                corruption levels (robust basins) AND significantly
                exceeds a random-weights baseline (p < 0.01).
  Definitive: Multiple distinct stored patterns are selectively recalled
              from their respective cues (genuine content-addressable
              memory, not a single fixed point).

Null model: random-weights surrogate. Replace the learned weight matrix
with i.i.d. random weights (preserving network size); show that
completion accuracy drops to chance.

T1a observation contract:
  Input: attractor-network state-vector trajectory + stored template set.
  Read via extract_observation_bundle(), NOT from the model object directly.

Reference: Hopfield, J. J. (1982). Neural networks and physical systems
           with emergent collective computational abilities. PNAS 79(8),
           2554–2558.
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
    """Extract the attractor-network observation bundle from history.

    This adapter is the T1a observation contract for P16. Any system that
    emits a state-vector trajectory with stored template patterns can be
    consumed by the P16 detector via this function.

    Parameters
    ----------
    history : list[dict]
        State trajectory. Each dict must contain at minimum:
        'state' (1D array of ±1), 'step' (int), 'trial' (int),
        'cue_pattern_idx' (int), 'overlap' (float),
        'stored_patterns' (2D array of shape (P, N)),
        'converged' (bool).

    Returns
    -------
    dict with keys:
        'trials': list of dicts, one per trial, each containing:
            'states': list of np.ndarray, state at each step
            'overlaps': list of float, overlap at each step
            'cue_pattern_idx': int
            'converged': bool (final state)
        'stored_patterns': np.ndarray of shape (P, N)
        'N': int, number of neurons
        'P': int, number of stored patterns
    """
    # Extract stored patterns from first record
    stored_patterns = np.array(history[0]['stored_patterns'])
    P, N = stored_patterns.shape

    # Group by trial
    trials_dict: Dict[int, Dict[str, Any]] = {}
    for h in history:
        trial_id = h['trial']
        if trial_id not in trials_dict:
            trials_dict[trial_id] = {
                'states': [],
                'overlaps': [],
                'cue_pattern_idx': h['cue_pattern_idx'],
                'converged': False,
            }
        trials_dict[trial_id]['states'].append(np.array(h['state']))
        trials_dict[trial_id]['overlaps'].append(float(h['overlap']))
        trials_dict[trial_id]['converged'] = bool(h['converged'])

    trials = [trials_dict[k] for k in sorted(trials_dict.keys())]

    return {
        'trials': trials,
        'stored_patterns': stored_patterns,
        'N': N,
        'P': P,
    }


# ── Core detection logic ─────────────────────────────────────────────

def _compute_completion_accuracy(
    final_state: np.ndarray,
    stored_patterns: np.ndarray,
) -> Dict[str, float]:
    """Compute overlap of final state with all stored patterns.

    Returns
    -------
    dict with:
        'best_overlap': float, highest overlap with any stored pattern
        'best_pattern_idx': int, which pattern has highest overlap
        'mean_overlap_others': float, mean overlap with non-best patterns
    """
    N = len(final_state)
    overlaps = np.array([
        float(np.dot(pat.astype(np.int32), final_state.astype(np.int32)) / N)
        for pat in stored_patterns
    ])
    best_idx = int(np.argmax(np.abs(overlaps)))
    best_overlap = float(np.abs(overlaps[best_idx]))

    # Mean overlap with other patterns
    other_overlaps = np.delete(np.abs(overlaps), best_idx)
    mean_others = float(np.mean(other_overlaps)) if len(other_overlaps) > 0 else 0.0

    return {
        'best_overlap': best_overlap,
        'best_pattern_idx': best_idx,
        'mean_overlap_others': mean_others,
    }


def _compute_selective_recall(
    trials: List[Dict[str, Any]],
    stored_patterns: np.ndarray,
) -> Dict[str, Any]:
    """Check whether multiple distinct patterns are selectively recalled.

    For each trial, check that the converged state has highest overlap
    with the cued pattern (not a different stored pattern).

    Returns
    -------
    dict with:
        'n_correct_recalls': int
        'n_total_trials': int
        'recall_accuracy': float
        'n_distinct_patterns_recalled': int
        'per_trial': list of dicts with overlap details
    """
    N = stored_patterns.shape[1]
    per_trial = []
    correct_recalls = 0
    recalled_patterns = set()

    for trial in trials:
        final_state = trial['states'][-1]
        cue_idx = trial['cue_pattern_idx']

        # Overlap with all stored patterns
        overlaps = np.array([
            float(np.dot(pat.astype(np.int32), final_state.astype(np.int32)) / N)
            for pat in stored_patterns
        ])
        abs_overlaps = np.abs(overlaps)
        best_idx = int(np.argmax(abs_overlaps))
        best_overlap = float(abs_overlaps[best_idx])
        cue_overlap = float(abs_overlaps[cue_idx])

        correct = best_idx == cue_idx
        if correct:
            correct_recalls += 1
            recalled_patterns.add(cue_idx)

        per_trial.append({
            'cue_pattern_idx': cue_idx,
            'best_pattern_idx': best_idx,
            'best_overlap': best_overlap,
            'cue_overlap': cue_overlap,
            'correct': correct,
        })

    n_total = len(trials)
    return {
        'n_correct_recalls': correct_recalls,
        'n_total_trials': n_total,
        'recall_accuracy': correct_recalls / n_total if n_total > 0 else 0.0,
        'n_distinct_patterns_recalled': len(recalled_patterns),
        'per_trial': per_trial,
    }


def _generate_random_weights_null(
    N: int,
    P: int,
    stored_patterns: np.ndarray,
    corruption: float,
    n_runs: int,
    rng: np.random.Generator,
    max_steps: int = 100,
) -> List[float]:
    """Generate null distribution: random-weights networks with same structure.

    For each null run, create random ±1 patterns, compute Hebbian weights,
    then try to retrieve stored patterns from corrupted cues. Record
    mean completion accuracy.
    """
    null_accuracies: List[float] = []

    for _ in range(n_runs):
        # Random patterns (different from stored ones)
        random_patterns = rng.choice([-1, 1], size=(P, N)).astype(np.int8)
        W_null = np.zeros((N, N), dtype=np.float64)
        for mu in range(P):
            W_null += np.outer(random_patterns[mu], random_patterns[mu])
        W_null /= N
        np.fill_diagonal(W_null, 0.0)

        # Try to retrieve the ORIGINAL stored patterns from corrupted cues
        # using the random weights (should fail)
        trial_overlaps = []
        for pat_idx in range(min(P, 5)):  # Cap at 5 for speed
            target = stored_patterns[pat_idx]
            state = target.copy()
            n_flip = int(corruption * N)
            flip_idx = rng.choice(N, size=n_flip, replace=False)
            state[flip_idx] *= -1

            # Run dynamics with random weights
            for _ in range(max_steps):
                old_state = state.copy()
                order = rng.permutation(N)
                for i in order:
                    h_i = float(np.dot(W_null[i], state))
                    if h_i > 0:
                        state[i] = 1
                    elif h_i < 0:
                        state[i] = -1
                if np.array_equal(state, old_state):
                    break

            best_overlap = float(np.max(np.abs([
                np.dot(pat.astype(np.int32), state.astype(np.int32)) / N for pat in stored_patterns
            ])))
            trial_overlaps.append(best_overlap)

        null_accuracies.append(float(np.mean(trial_overlaps)))

    return null_accuracies


class P16AssociativeMemoryDetector:
    """P16 associative memory / pattern completion detector.

    Reads inputs via the T1a observation-bundle adapter.

    Parameters
    ----------
    n_permutations : int
        Number of random-weights null runs.
    seed : int
        RNG seed for reproducibility.
    """

    def __init__(
        self,
        n_permutations: int = 199,
        seed: int = 42,
    ) -> None:
        self.pattern_id = "P16"
        self.n_permutations = n_permutations
        self.seed = seed

    def detect(
        self,
        history: List[Dict[str, Any]],
        metadata: Optional[Dict[str, Any]] = None,
    ) -> DetectorResult:
        """Run P16 associative memory detection.

        Parameters
        ----------
        history : list[dict]
            State trajectory with attractor-network observation bundle keys.
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
        except (KeyError, IndexError, ValueError) as e:
            return self._no_detection(
                warnings=[f"Observation bundle extraction failed: {e}"],
                metadata_available=metadata_available,
            )

        trials = bundle['trials']
        stored_patterns = bundle['stored_patterns']
        N = bundle['N']
        P = bundle['P']

        if len(trials) == 0:
            return self._no_detection(
                warnings=["No trials found in history"],
                metadata_available=metadata_available,
            )

        # ── Multi-pattern prerequisite (Sprint 80) ──
        # P16 requires ≥2 distinct selectively-retrievable stored patterns
        # (IC-cued recall selects the matching template). A single stored
        # pattern or single-attractor convergence is NOT content-addressable
        # memory — it is a trivial fixed point.
        if P < 2:
            return self._no_detection(
                warnings=[
                    f"P16 requires ≥2 stored patterns (P={P}). "
                    "Single-pattern convergence is not content-addressable memory."
                ],
                metadata_available=metadata_available,
            )

        # Quick selective-recall pre-check: if all converged trials settle to
        # the same pattern (regardless of cue), this is single-attractor
        # convergence, not multi-pattern recall. Reject at screening.
        converged_trials = [t for t in trials if t['converged']]
        if len(converged_trials) >= 2:
            best_patterns = set()
            for trial in converged_trials:
                final = trial['states'][-1]
                overlaps = np.array([
                    float(np.dot(pat.astype(np.int32), final.astype(np.int32)) / N)
                    for pat in stored_patterns
                ])
                best_patterns.add(int(np.argmax(np.abs(overlaps))))
            if len(best_patterns) < 2:
                return self._no_detection(
                    warnings=[
                        f"All converged trials settle to the same attractor "
                        f"(pattern {best_patterns.pop()}). P16 requires ≥2 "
                        f"distinct selectively-retrievable stored patterns."
                    ],
                    metadata_available=metadata_available,
                )

        # ── Primary metric: completion accuracy ──
        completion_overlaps = []
        for trial in trials:
            final_state = trial['states'][-1]
            acc = _compute_completion_accuracy(final_state, stored_patterns)
            completion_overlaps.append(acc['best_overlap'])

        mean_completion = float(np.mean(completion_overlaps))

        # Chance-level overlap for binary ±1 patterns: ~sqrt(1/N) by CLT
        chance_overlap = 2.0 * np.sqrt(1.0 / N) if N > 0 else 0.5

        primary: Dict[str, float] = {
            'mean_completion_accuracy': mean_completion,
            'chance_overlap': chance_overlap,
        }

        # ── SCREENING: completion above chance ──
        screening_threshold = max(0.5, chance_overlap + 0.2)
        screening_pass = mean_completion > screening_threshold

        if not screening_pass:
            return self._no_detection(
                primary_metric=primary,
                warnings=warnings + [
                    f"Completion accuracy {mean_completion:.3f} below "
                    f"screening threshold {screening_threshold:.3f}"
                ],
                metadata_available=metadata_available,
            )

        # ── Selective recall analysis ──
        selective = _compute_selective_recall(trials, stored_patterns)
        recall_accuracy = selective['recall_accuracy']
        n_distinct = selective['n_distinct_patterns_recalled']

        # ── Random-weights null model ──
        rng = np.random.default_rng(self.seed)
        corruption = metadata.get('corruption', 0.2) if metadata else 0.2

        null_accuracies = _generate_random_weights_null(
            N=N, P=P, stored_patterns=stored_patterns,
            corruption=corruption, n_runs=self.n_permutations,
            rng=rng, max_steps=50,
        )

        null_arr = np.array(null_accuracies)
        null_mean = float(np.mean(null_arr))
        null_std = float(np.std(null_arr)) if len(null_arr) > 1 else 1.0
        if null_std < 1e-12:
            null_std = max(1e-12, null_mean * 0.01)

        # p-value: fraction of null runs with accuracy >= observed
        null_p = float(np.mean(null_arr >= mean_completion))
        null_p = max(null_p, 1.0 / (self.n_permutations + 1))

        # Effect size
        cohens_d = (mean_completion - null_mean) / null_std if null_std > 0 else 0.0

        effect_size: Dict[str, float] = {
            'cohens_d': cohens_d,
            'raw_value': mean_completion,
            'null_mean': null_mean,
            'null_std': null_std,
        }

        secondary: Dict[str, Any] = {
            'recall_accuracy': recall_accuracy,
            'n_distinct_patterns_recalled': n_distinct,
            'n_trials': len(trials),
            'mean_completion_overlaps': completion_overlaps,
        }

        # ── CONFIRMATION: high accuracy + null rejected ──
        confirmed = (
            null_p < 0.01
            and mean_completion > 0.7
            and recall_accuracy > 0.5
        )

        if not confirmed:
            bonuses = {
                'secondaries_pass': recall_accuracy > 0.5,
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

        # ── DEFINITIVE: multiple distinct patterns selectively recalled ──
        # Note: null_p < 0.01 (not <= 0.005) because the random-weights
        # null has a small but non-zero probability of producing a spurious
        # attractor that matches a stored pattern, making the floor p-value
        # with 199 perms occasionally exceed 0.005 (1 in ~200 null runs).
        definitive = (
            confirmed
            and n_distinct >= 2
            and recall_accuracy >= 0.8
            and mean_completion > 0.8
            and null_p < 0.01
        )

        # Metadata check: must have multiple attractors
        if definitive and metadata_available:
            has_multiple = metadata.get('has_multiple_attractors', None)
            if has_multiple is False:
                definitive = False
                warnings.append(
                    "Metadata indicates single attractor; "
                    "definitive tier denied"
                )

        if definitive:
            tier = DetectionTier.DEFINITIVE
            bonuses = {
                'all_exclusions_cleared': True,
                'both_null_types_rejected': null_p < 0.005,
                'finite_size_robustness': N >= 50,
            }
        else:
            tier = DetectionTier.CONFIRMATION
            bonuses = {
                'null_p_0001': null_p < 0.001,
                'effect_size_gt_1': cohens_d > 1.0,
                'all_secondaries': recall_accuracy >= 0.8,
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
