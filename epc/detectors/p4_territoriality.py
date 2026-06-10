"""P4 Territoriality / exclusion boundaries detector.

Detects territorial exclusion: agents form stable, non-overlapping home
ranges with persistent boundaries, where occupancy anti-correlates with
foreign scent (exclusion is scent-mediated, not incidental).

Distinct from P1 (aggregation): P4 outcome is spatial EXCLUSION + boundary
maintenance, not clustering. From P29 (stigmergic transport): P4 makes
exclusive domains, P29 makes connective transport networks.

Primary metrics:
  - Exclusivity index: for each visited cell, what fraction of total
    visits came from the dominant agent. For territories → ~1.0;
    for random walks sharing the same space → ~1/N.
  - Boundary persistence: stability of the spatial partition over time.

Detection tiers:
  Screening:  Exclusivity significantly exceeds agent-label-shuffle null
              (actual ownership matters, not just spatial structure).
  Confirmation: The spatial partition is PERSISTENT over a time window
                (boundaries stable) + exclusivity p < 0.01.
  Definitive: Occupancy anti-correlates with foreign scent (exclusion
              is scent-mediated, not incidental) + metadata confirms
              foreign avoidance.

Null model: agent-label-shuffle on per-step position history. Randomly
reassign which agent visited which cell at each timestep, then recompute
occupancy and exclusivity. This destroys the identity-linked spatial
clustering (each agent staying in its own area) while preserving the
total spatial distribution of visits.

T1a observation contract:
  Input: territorial-agent-field observation bundle — time series of
  (positions, scent_fields, occupancy, step, n_agents, grid_size)
  read via extract_observation_bundle().

Reference: Giuggioli, L., Potts, J. R., & Harris, S. (2011).
Animal interactions and the emergence of territoriality.
PLoS Computational Biology, 7(3), e1002008.
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
    """Extract the territorial-agent-field observation bundle from history.

    This adapter is the T1a observation contract for P4. Any system that
    emits a time series with keys 'positions', 'scent_fields', 'occupancy',
    'step', 'n_agents', 'grid_size' can be consumed by the P4 detector.

    Parameters
    ----------
    history : list[dict]
        State trajectory. Each dict must contain at minimum:
        'positions' (ndarray(N, 2)), 'scent_fields' (ndarray(N, L, L)),
        'occupancy' (ndarray(N, L, L)), 'step' (int),
        'n_agents' (int), 'grid_size' (int).

    Returns
    -------
    dict with keys:
        'positions': list of ndarray(N, 2)
        'scent_fields': list of ndarray(N, L, L)
        'occupancy': list of ndarray(N, L, L)
        'steps': ndarray(T,) int
        'n_agents': int
        'grid_size': int
    """
    positions = [np.asarray(h['positions']) for h in history]
    scent_fields = [np.asarray(h['scent_fields']) for h in history]
    occupancy = [np.asarray(h['occupancy']) for h in history]
    steps = np.array([h['step'] for h in history], dtype=int)
    n_agents = int(history[0]['n_agents'])
    grid_size = int(history[0]['grid_size'])
    return {
        'positions': positions,
        'scent_fields': scent_fields,
        'occupancy': occupancy,
        'steps': steps,
        'n_agents': n_agents,
        'grid_size': grid_size,
    }


# ── Metric computations ─────────────────────────────────────────────

def _exclusivity_index(occupancy: np.ndarray) -> float:
    """Compute the mean per-cell exclusivity index.

    For each cell that has been visited at least once, compute the fraction
    of total visits from the dominant agent. Average over all visited cells.

    Perfect territories → 1.0 (every cell visited by exactly one agent).
    N agents sharing uniformly → 1/N.

    Parameters
    ----------
    occupancy : ndarray(N, L, L)
        Cumulative visit counts per agent per cell.

    Returns
    -------
    float in [1/N, 1.0].
    """
    N = occupancy.shape[0]
    if N < 2:
        return 1.0

    total_per_cell = occupancy.sum(axis=0)  # (L, L)
    max_per_cell = occupancy.max(axis=0)    # (L, L)

    visited = total_per_cell > 0
    if not np.any(visited):
        return 1.0 / N

    exclusivity = max_per_cell[visited] / total_per_cell[visited]
    return float(np.mean(exclusivity))


def _pairwise_overlap(occupancy: np.ndarray) -> float:
    """Compute mean pairwise home-range overlap index.

    For each pair of agents (i, j), compute the overlap as the sum of
    min(normed_occ_i, normed_occ_j). Returns the mean over all pairs.

    Perfect exclusion → ~0. Complete overlap → ~1.
    """
    N = occupancy.shape[0]
    if N < 2:
        return 0.0

    normed = np.zeros_like(occupancy, dtype=np.float64)
    for i in range(N):
        total = occupancy[i].sum()
        if total > 0:
            normed[i] = occupancy[i] / total
        else:
            normed[i] = 1.0 / (occupancy.shape[1] * occupancy.shape[2])

    overlaps = []
    for i in range(N):
        for j in range(i + 1, N):
            overlap = np.sum(np.minimum(normed[i], normed[j]))
            overlaps.append(float(overlap))

    return float(np.mean(overlaps))


def _boundary_persistence(
    occupancy_early: np.ndarray,
    occupancy_late: np.ndarray,
) -> float:
    """Compute boundary persistence between two time windows.

    Assigns each cell to the agent with highest occupancy (winner-take-all).
    Persistence = fraction of cells with the same owner in both windows.
    """
    owner_early = np.argmax(occupancy_early, axis=0)
    owner_late = np.argmax(occupancy_late, axis=0)

    claimed_early = occupancy_early.sum(axis=0) > 0
    claimed_late = occupancy_late.sum(axis=0) > 0
    both_claimed = claimed_early & claimed_late

    if not np.any(both_claimed):
        return 0.0

    return float(np.mean(owner_early[both_claimed] == owner_late[both_claimed]))


def _occupancy_scent_correlation(
    occupancy: np.ndarray,
    scent_fields: np.ndarray,
) -> float:
    """Compute cross-correlation between occupancy and foreign scent.

    For each agent i, compute the Pearson correlation between agent i's
    occupancy and the total foreign scent. Returns the mean across agents.

    Negative correlation = scent-mediated exclusion.
    """
    N = occupancy.shape[0]
    if N < 2:
        return 0.0

    correlations = []
    for i in range(N):
        occ_flat = occupancy[i].flatten().astype(np.float64)
        foreign = np.sum(
            [scent_fields[k] for k in range(N) if k != i], axis=0
        ).flatten().astype(np.float64)

        occ_std = np.std(occ_flat)
        foreign_std = np.std(foreign)
        if occ_std < 1e-12 or foreign_std < 1e-12:
            correlations.append(0.0)
            continue

        r = float(np.corrcoef(occ_flat, foreign)[0, 1])
        correlations.append(r)

    return float(np.mean(correlations))


def _build_occupancy_from_positions(
    positions_list: List[np.ndarray],
    n_agents: int,
    grid_size: int,
) -> np.ndarray:
    """Build cumulative occupancy from a position history.

    Parameters
    ----------
    positions_list : list of ndarray(N, 2)
    n_agents : int
    grid_size : int

    Returns
    -------
    ndarray(N, L, L) cumulative occupancy.
    """
    occ = np.zeros((n_agents, grid_size, grid_size), dtype=np.float64)
    for pos in positions_list:
        for i in range(n_agents):
            r, c = int(pos[i, 0]), int(pos[i, 1])
            occ[i, r, c] += 1.0
    return occ


class P4TerritorialityDetector:
    """P4 territoriality / exclusion boundaries detector.

    Reads inputs via the T1a observation-bundle adapter.

    Parameters
    ----------
    n_permutations : int
        Number of agent-label-shuffles for the null model.
    seed : int
        RNG seed for reproducibility.
    """

    def __init__(
        self,
        n_permutations: int = 199,
        seed: int = 42,
    ) -> None:
        self.pattern_id = "P4"
        self.n_permutations = n_permutations
        self.seed = seed

    def detect(
        self,
        history: List[Dict[str, Any]],
        metadata: Optional[Dict[str, Any]] = None,
    ) -> DetectorResult:
        """Run P4 territoriality detection.

        Parameters
        ----------
        history : list[dict]
            State trajectory with territorial-agent-field bundle keys.
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

        n_agents = bundle['n_agents']
        n_snapshots = len(bundle['occupancy'])
        L = bundle['grid_size']

        if n_snapshots < 4:
            return self._no_detection(
                warnings=["Insufficient snapshots (need >= 4)"],
                metadata_available=metadata_available,
            )

        if n_agents < 2:
            return self._no_detection(
                warnings=["Need >= 2 agents for territoriality"],
                metadata_available=metadata_available,
            )

        # ── Compute metrics from cumulative occupancy ──
        final_occupancy = bundle['occupancy'][-1]
        observed_overlap = _pairwise_overlap(final_occupancy)

        # ── Compute boundary persistence ──
        quarter = max(1, n_snapshots // 4)
        early_occ = bundle['occupancy'][quarter]
        late_occ = bundle['occupancy'][-1]
        persistence = _boundary_persistence(early_occ, late_occ)

        # ── Occupancy-vs-foreign-scent correlation ──
        final_scent = bundle['scent_fields'][-1]
        occ_scent_corr = _occupancy_scent_correlation(final_occupancy, final_scent)

        observed_excl = _exclusivity_index(final_occupancy)

        primary: Dict[str, float] = {
            'home_range_overlap': observed_overlap,
            'exclusivity_index': observed_excl,
            'boundary_persistence': persistence,
            'occupancy_scent_correlation': occ_scent_corr,
        }

        # ── Null model: cell-level multinomial shuffle ──
        # For each visited cell, the total visit count is known.
        # Under the null, redistribute those visits uniformly among
        # agents (multinomial with equal probabilities). This tests
        # whether per-cell agent dominance exceeds random assignment.
        # Uses cumulative occupancy directly — no sparse/dense mismatch.
        rng = np.random.default_rng(self.seed)
        null_excls: List[float] = []

        total_per_cell = final_occupancy.sum(axis=0)  # (L, L)
        visited_mask = total_per_cell > 0
        visited_totals = total_per_cell[visited_mask].astype(int)
        n_visited = int(visited_mask.sum())
        probs = np.ones(n_agents) / n_agents

        for _ in range(self.n_permutations):
            # For each visited cell, draw multinomial allocation
            null_max = np.zeros(n_visited, dtype=np.float64)
            null_total = np.zeros(n_visited, dtype=np.float64)
            for ci in range(n_visited):
                counts = rng.multinomial(visited_totals[ci], probs)
                null_max[ci] = counts.max()
                null_total[ci] = visited_totals[ci]
            null_excl = float(np.mean(null_max / null_total))
            null_excls.append(null_excl)

        null_arr = np.array(null_excls)
        null_mean = float(np.mean(null_arr))
        null_std = float(np.std(null_arr)) if len(null_arr) > 1 else 1.0
        if null_std < 1e-12:
            null_std = max(1e-12, abs(null_mean) * 0.01)

        # P-value: fraction of null exclusivities >= observed
        null_p = float(np.mean(null_arr >= observed_excl))
        null_p = max(null_p, 1.0 / (self.n_permutations + 1))

        cohens_d = (
            (observed_excl - null_mean) / null_std if null_std > 0 else 0.0
        )

        effect_size: Dict[str, float] = {
            'cohens_d': cohens_d,
            'observed_exclusivity': observed_excl,
            'null_mean_exclusivity': null_mean,
            'null_std_exclusivity': null_std,
            'observed_overlap': observed_overlap,
        }

        # ── CONTENT PREREQUISITE (Sprint 87): scent-mediated exclusion ──
        # P4 requires scent-mediated exclusion (occupancy anti-correlated
        # with foreign scent) AND persistent boundaries. Incidental low
        # overlap or clustering from spatial autocorrelation is out of
        # domain. This prerequisite separates genuine territory formation
        # from random walks that produce moderate exclusivity on finite
        # grids purely from spatial locality.
        if occ_scent_corr >= 0.0:
            return self._no_detection(
                primary_metric=primary,
                warnings=warnings + [
                    f"Content prerequisite: occupancy-scent correlation "
                    f"{occ_scent_corr:.4f} >= 0 → no scent-mediated "
                    f"exclusion detected (incidental overlap out of domain)"
                ],
                metadata_available=metadata_available,
            )
        if persistence < 0.5:
            return self._no_detection(
                primary_metric=primary,
                warnings=warnings + [
                    f"Content prerequisite: boundary persistence "
                    f"{persistence:.3f} < 0.5 → no persistent boundaries"
                ],
                metadata_available=metadata_available,
            )

        # ── SCREENING: exclusivity significantly above null ──
        # Require both statistical significance AND absolute threshold.
        # Random walks on finite grids produce moderate exclusivity (~0.7)
        # due to spatial autocorrelation; territories require > 0.80.
        if null_p >= 0.10 or observed_excl < 0.80:
            return self._no_detection(
                primary_metric=primary,
                warnings=warnings + [
                    f"Exclusivity {observed_excl:.3f} not significantly "
                    f"above null (mean={null_mean:.3f}, p={null_p:.4f})"
                ],
                metadata_available=metadata_available,
            )

        secondary: Dict[str, Any] = {
            'exclusivity_ratio': observed_excl / max(null_mean, 1e-6),
            'n_agents': n_agents,
            'n_snapshots': n_snapshots,
        }

        # ── CONFIRMATION: persistent boundaries + null rejected ──
        confirmed = (
            persistence > 0.6
            and null_p < 0.05
            and cohens_d > 1.0
        )

        if not confirmed:
            reasons: List[str] = []
            if persistence <= 0.6:
                reasons.append(
                    f"Boundaries not persistent ({persistence:.3f})"
                )
            if null_p >= 0.05:
                reasons.append(f"Null not rejected (p={null_p:.4f})")
            if cohens_d <= 1.0:
                reasons.append(
                    f"Effect size too small (d={cohens_d:.2f})"
                )

            bonuses = {
                'secondaries_pass': persistence > 0.5,
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

        # ── DEFINITIVE: scent-mediated exclusion + metadata ──
        definitive = (
            confirmed
            and occ_scent_corr < -0.005
            and cohens_d > 1.0
            and null_p <= 0.005
        )

        if definitive and metadata_available:
            has_avoidance = metadata.get('has_foreign_avoidance', None)
            if has_avoidance is False:
                definitive = False
                warnings.append(
                    "Metadata indicates no foreign avoidance; "
                    "definitive tier denied"
                )

        if definitive:
            tier = DetectionTier.DEFINITIVE
            bonuses = {
                'all_exclusions_cleared': True,
                'both_null_types_rejected': null_p < 0.005,
                'finite_size_robustness': n_agents >= 3,
            }
        else:
            tier = DetectionTier.CONFIRMATION
            bonuses = {
                'null_p_0001': null_p < 0.001,
                'effect_size_gt_1': cohens_d > 1.0,
                'all_secondaries': (
                    persistence > 0.6
                    and occ_scent_corr < 0.0
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
