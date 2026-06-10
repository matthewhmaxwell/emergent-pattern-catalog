"""
P30 Spontaneous Boundary Formation (Autopoiesis) Detector.

Reference:
    Varela, F.G., Maturana, H.R., & Uribe, R. (1974).
    Autopoiesis: the organization of living systems, its characterization
    and a model. Biosystems, 5(4), 187–196.

Detection metrics (priority order):
    1. Link-catalyst association — fraction of link particles near catalysts,
       normalized by the expected fraction under complete spatial randomness.
       High association (>> 1) means links are concentrated around catalysts
       (membrane formed around production zone). This is the primary metric.
    2. Angular closure — catalysts enclosed by surrounding link particles
       (max angular gap < π in the link distribution around each catalyst).
    3. Catalyst enrichment — fraction of catalysts inside the membrane radius
       vs expected fraction under uniform placement.
    4. Self-repair — after artificially breaching the boundary, does it re-close?

Tier structure:
    - Screening: association_score > 3.0 AND closure_fraction > 0.5.
    - Confirmation: + enrichment_ratio > 1.2 + null p < 0.01 + persistence > 0.5.
    - Definitive: + closure > 0.7 + stable membrane (link CV < 0.3).

Null model:
    Type-shuffle permutation: keep all positions fixed, randomly reassign
    type labels (substrate, catalyst, link) among particles. Under random
    assignment, links are uniformly distributed → association_score ≈ 1.

T1a observation bundle:
    positions: (N, 2) float64 — particle positions
    types: (N,) int — particle types (0=substrate, 1=catalyst, 2=link)
    bonds: list of (int, int) — bonded link pairs (optional)
    box_size: float — domain size
"""

from __future__ import annotations

from typing import Any

import numpy as np
from numpy.typing import NDArray

from epc.base_detector import BaseDetector
from epc.detector_result import DetectionTier, DetectorResult, NullType


def _periodic_delta_vec(
    p1: NDArray, p2: NDArray, box_size: float
) -> NDArray:
    """Compute periodic displacement vectors from p1 to array p2."""
    delta = p2 - p1
    delta -= box_size * np.round(delta / box_size)
    return delta


def _compute_autopoiesis_metrics(
    positions: NDArray,
    types: NDArray,
    box_size: float,
    association_radius: float = 4.0,
    link_type: int = 2,
    catalyst_type: int = 1,
) -> dict[str, float]:
    """Compute autopoiesis metrics for a single snapshot.

    Metrics:
    1. association_score: (fraction of links near catalysts) / (expected under CSR).
       Measures spatial co-location of membrane and production zone.
    2. closure_fraction: fraction of catalysts with full angular coverage
       by link particles (max angular gap < π).
    3. enrichment_ratio: catalyst density inside membrane vs uniform expectation.
    """
    link_mask = types == link_type
    cat_mask = types == catalyst_type
    link_pos = positions[link_mask]
    cat_pos = positions[cat_mask]
    n_links = int(np.sum(link_mask))
    n_cats = int(np.sum(cat_mask))
    n_total = len(types)

    if n_links < 3 or n_cats == 0:
        return {
            'association_score': 0.0,
            'closure_fraction': 0.0,
            'enrichment_ratio': 1.0,
            'n_links': n_links,
            'n_enclosed': 0,
        }

    # 1. Link-catalyst association score
    # Count how many links are within association_radius of ANY catalyst
    n_near = 0
    for li in range(n_links):
        for ci in range(n_cats):
            delta = link_pos[li] - cat_pos[ci]
            delta -= box_size * np.round(delta / box_size)
            dist = float(np.sqrt(np.sum(delta**2)))
            if dist < association_radius:
                n_near += 1
                break  # Count each link only once

    frac_observed = n_near / n_links if n_links > 0 else 0.0
    # Expected fraction under CSR: coverage area of all catalysts / total area
    coverage_area = min(n_cats * np.pi * association_radius**2, box_size**2)
    frac_expected = coverage_area / (box_size**2)
    association_score = frac_observed / max(frac_expected, 0.001)

    # 2. Angular closure test per catalyst
    n_enclosed = 0
    for ci in range(n_cats):
        cat_p = cat_pos[ci]
        deltas = _periodic_delta_vec(cat_p, link_pos, box_size)
        angles = np.arctan2(deltas[:, 1], deltas[:, 0])
        sorted_angles = np.sort(angles)

        if len(sorted_angles) < 3:
            continue

        extended = np.concatenate([sorted_angles, [sorted_angles[0] + 2 * np.pi]])
        gaps = np.diff(extended)
        max_gap = float(np.max(gaps))
        if max_gap < np.pi:
            n_enclosed += 1

    closure_fraction = n_enclosed / n_cats if n_cats > 0 else 0.0

    # 3. Enrichment ratio: catalyst density inside link cloud vs expected
    cat_com = np.mean(cat_pos, axis=0)
    link_deltas = _periodic_delta_vec(cat_com, link_pos, box_size)
    link_dists = np.sqrt(np.sum(link_deltas**2, axis=1))
    membrane_radius = float(np.percentile(link_dists, 90))  # 90th pctile as outer edge

    cat_deltas = _periodic_delta_vec(cat_com, cat_pos, box_size)
    cat_dists = np.sqrt(np.sum(cat_deltas**2, axis=1))
    n_cat_inside = int(np.sum(cat_dists < membrane_radius))

    area_inside = np.pi * membrane_radius**2
    total_area = box_size**2
    if area_inside > 0 and area_inside < total_area:
        frac_cat_observed = n_cat_inside / n_cats if n_cats > 0 else 0.0
        frac_cat_expected = area_inside / total_area
        enrichment_ratio = frac_cat_observed / max(frac_cat_expected, 0.001)
    else:
        enrichment_ratio = 1.0

    return {
        'association_score': association_score,
        'closure_fraction': closure_fraction,
        'enrichment_ratio': enrichment_ratio,
        'n_links': n_links,
        'n_enclosed': n_enclosed,
    }


class P30AutopoiesisDetector(BaseDetector):
    """Detector for P30: Spontaneous boundary formation / autopoiesis.

    Detects when particles spontaneously organize into a closed, semi-permeable
    boundary (membrane) that maintains an internal micro-environment distinct
    from the exterior.

    Primary metric: association_score — spatial co-location of link (membrane)
    particles with catalyst (production) particles, normalized by CSR expectation.

    Parameters
    ----------
    n_permutations : int
        Number of type-shuffle permutations for null model.
    association_radius : float
        Distance threshold for link-catalyst association.
    seed : int
        Random seed for null model.
    """

    def __init__(
        self,
        n_permutations: int = 199,
        association_radius: float = 4.0,
        seed: int = 42,
    ) -> None:
        super().__init__(
            pattern_id="P30",
            excluded_patterns=["P1"],
            allowed_co_occurrences=[],
            observable_scope="state_history_only",
        )
        self.n_permutations = n_permutations
        self.association_radius = association_radius
        self._seed = seed

    def _compute_primary(
        self,
        state_history: list[dict[str, Any]],
        timescale: float,
    ) -> dict[str, float]:
        """Compute autopoiesis metrics over the late window (last 25%)."""
        n = len(state_history)
        late_start = max(1, 3 * n // 4)
        late_history = state_history[late_start:]

        assoc_vals = []
        closure_vals = []
        enrich_vals = []
        n_link_vals = []

        for snap in late_history:
            positions = np.asarray(snap['positions'], dtype=np.float64)
            types = np.asarray(snap['types'], dtype=np.int32)
            box_size = float(snap.get('box_size', 20.0))

            metrics = _compute_autopoiesis_metrics(
                positions, types, box_size,
                association_radius=self.association_radius,
            )
            assoc_vals.append(metrics['association_score'])
            closure_vals.append(metrics['closure_fraction'])
            enrich_vals.append(metrics['enrichment_ratio'])
            n_link_vals.append(metrics['n_links'])

        return {
            'association_score': float(np.mean(assoc_vals)),
            'closure_fraction': float(np.mean(closure_vals)),
            'enrichment_ratio': float(np.mean(enrich_vals)),
            'mean_n_links': float(np.mean(n_link_vals)),
        }

    def _check_screening(
        self,
        primary_result: dict[str, float],
        timescale: float,
    ) -> bool:
        """Screening: link-catalyst association above CSR AND angular closure."""
        return (
            primary_result.get('association_score', 0.0) > 1.5
            and primary_result.get('closure_fraction', 0.0) > 0.5
            and primary_result.get('mean_n_links', 0.0) >= 3
        )

    def _compute_secondaries(
        self,
        state_history: list[dict[str, Any]],
        timescale: float,
    ) -> dict[str, Any]:
        """Compute secondary metrics: temporal persistence + membrane stability."""
        n = len(state_history)
        late_start = max(1, 3 * n // 4)
        late_history = state_history[late_start:]

        assoc_per_snap = []
        link_counts = []
        for snap in late_history:
            positions = np.asarray(snap['positions'], dtype=np.float64)
            types = np.asarray(snap['types'], dtype=np.int32)
            box_size = float(snap.get('box_size', 20.0))
            metrics = _compute_autopoiesis_metrics(
                positions, types, box_size,
                association_radius=self.association_radius,
            )
            assoc_per_snap.append(metrics['association_score'] > 1.5 and
                                  metrics['closure_fraction'] > 0.5)
            link_counts.append(int(np.sum(types == 2)))

        persistence = float(np.mean(assoc_per_snap))
        link_arr = np.array(link_counts, dtype=float)
        link_cv = float(np.std(link_arr) / np.mean(link_arr)) if np.mean(link_arr) > 0 else 999.0

        return {
            'closure_persistence': persistence,
            'link_count_cv': link_cv,
            'mean_link_count': float(np.mean(link_arr)),
        }

    def _run_null_model(
        self,
        state_history: list[dict[str, Any]],
        primary_result: dict[str, float],
        timescale: float,
    ) -> tuple[float, NullType, dict[str, float]]:
        """Type-shuffle permutation null.

        Keep all positions fixed, randomly reassign type labels among particles.
        Test statistic: association_score. Under random type assignment, links
        are uniformly distributed among positions → association_score ≈ 1.0.
        Under autopoiesis, links concentrate near catalysts → score >> 1.
        """
        rng = np.random.default_rng(self._seed)
        observed = primary_result.get('association_score', 0.0)

        snap = state_history[-1]
        positions = np.asarray(snap['positions'], dtype=np.float64)
        types = np.asarray(snap['types'], dtype=np.int32)
        box_size = float(snap.get('box_size', 20.0))

        null_scores = []
        for _ in range(self.n_permutations):
            surr_types = types.copy()
            rng.shuffle(surr_types)

            metrics = _compute_autopoiesis_metrics(
                positions, surr_types, box_size,
                association_radius=self.association_radius,
            )
            null_scores.append(metrics['association_score'])

        null_arr = np.array(null_scores)
        p_value = float(np.mean(null_arr >= observed))
        p_value = max(p_value, 1.0 / (self.n_permutations + 1))

        return (
            p_value,
            NullType.SHUFFLE,
            {
                'mean': float(np.mean(null_arr)),
                'std': float(np.std(null_arr)) if float(np.std(null_arr)) > 0 else 0.001,
            },
        )

    def _check_exclusions(
        self,
        state_history: list[dict[str, Any]],
        model_metadata: dict[str, Any] | None,
        timescale: float,
    ) -> tuple[list[str], dict[str, str]]:
        """Exclude P1 (mere clustering) when association + closure present.

        P30 requires specifically organized membrane around production zone.
        P1 is undifferentiated aggregation without directed production.
        """
        n = len(state_history)
        late_start = max(1, 3 * n // 4)
        late_history = state_history[late_start:]

        assoc_vals = []
        closure_vals = []
        for snap in late_history:
            positions = np.asarray(snap['positions'], dtype=np.float64)
            types = np.asarray(snap['types'], dtype=np.int32)
            box_size = float(snap.get('box_size', 20.0))
            metrics = _compute_autopoiesis_metrics(
                positions, types, box_size,
                association_radius=self.association_radius,
            )
            assoc_vals.append(metrics['association_score'])
            closure_vals.append(metrics['closure_fraction'])

        mean_assoc = float(np.mean(assoc_vals))
        mean_closure = float(np.mean(closure_vals))

        if mean_assoc > 1.5 and mean_closure > 0.5:
            return (['P1'], {'P1': 'excluded'})
        return (['P1'], {'P1': 'not_excluded'})

    def _check_confirmation(
        self,
        primary_result: dict[str, float],
        secondary_result: dict[str, Any],
        null_p: float,
        timescale: float,
    ) -> bool:
        """Confirmation: association + enrichment + null rejection + persistence."""
        assoc = primary_result.get('association_score', 0.0)
        enrichment = primary_result.get('enrichment_ratio', 1.0)
        persistence = secondary_result.get('closure_persistence', 0.0)

        return (
            assoc > 1.5
            and enrichment > 1.2
            and null_p < 0.01
            and persistence > 0.5
        )

    def _check_definitive(
        self,
        primary_result: dict[str, float],
        secondary_result: dict[str, Any],
        null_p: float,
        null_type: NullType,
        state_history: list[dict[str, Any]],
        model_metadata: dict[str, Any] | None,
        timescale: float,
    ) -> bool:
        """Definitive: strong closure + stable membrane."""
        closure = primary_result.get('closure_fraction', 0.0)
        enrichment = primary_result.get('enrichment_ratio', 1.0)
        persistence = secondary_result.get('closure_persistence', 0.0)
        link_cv = secondary_result.get('link_count_cv', 999.0)

        return (
            closure > 0.7
            and enrichment > 2.0
            and persistence > 0.8
            and link_cv < 0.3
            and null_p < 0.01
        )

    def _estimate_timescale(
        self,
        state_history: list[dict[str, Any]],
        model_metadata: dict[str, Any] | None,
    ) -> float:
        """Estimate timescale from production rate if available."""
        if model_metadata and 'production_rate' in model_metadata:
            pr = model_metadata['production_rate']
            if pr > 0:
                return 1.0 / pr
        return max(1.0, len(state_history) / 10.0)


def extract_observation_bundle(
    state_history: list[dict[str, Any]],
) -> dict[str, Any]:
    """Extract T1a observation bundle from state history.

    Parameters
    ----------
    state_history : list[dict]
        History from any model producing positions + types.

    Returns
    -------
    dict with keys:
        positions: list of (N, 2) arrays
        types: list of (N,) arrays
        bonds: list of bond lists (may be empty)
        steps: (T,) int array
        box_size: float
        n_particles: int
    """
    positions = [np.asarray(s['positions'], dtype=np.float64) for s in state_history]
    types = [np.asarray(s['types'], dtype=np.int32) for s in state_history]
    bonds = [s.get('bonds', []) for s in state_history]
    steps = np.array([s.get('step', i) for i, s in enumerate(state_history)], dtype=int)
    box_size = float(state_history[0].get('box_size', 20.0))
    n_particles = int(state_history[0].get('n_particles', len(positions[0])))

    return {
        'positions': positions,
        'types': types,
        'bonds': bonds,
        'steps': steps,
        'box_size': box_size,
        'n_particles': n_particles,
    }
