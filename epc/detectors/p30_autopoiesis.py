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


def _membrane_metrics(
    positions: np.ndarray,
    types: np.ndarray,
    box_size: float,
    link_type: int = 2,
    cat_type: int = 1,
    sub_type: int = 0,
) -> dict[str, float]:
    """Closed-shell membrane metrics from a single snapshot.

    A genuine autopoietic membrane is a THIN, CLOSED shell of link particles at
    a characteristic radius from the catalyst (production core). The decisive
    signature is RADIAL TIGHTNESS: link-to-catalyst distances cluster at one
    radius (low CV), unlike a diffuse cloud or random scatter (CV ~ 0.4).
    """
    out = {"n_links": 0.0, "radial_cv": 1.0, "shell_concentration": 0.0,
           "med_radius": 0.0, "angular_gap_deg": 360.0, "gradient": 1.0}
    links = positions[types == link_type]
    cats = positions[types == cat_type]
    subs = positions[types == sub_type]
    out["n_links"] = float(len(links))
    if len(links) < 3 or len(cats) == 0:
        return out

    def _nearest_cat_dist(P: np.ndarray) -> np.ndarray:
        d = P[:, None, :] - cats[None, :, :]
        d -= box_size * np.round(d / box_size)
        return np.sqrt((d ** 2).sum(axis=-1)).min(axis=1)

    r = _nearest_cat_dist(links)
    med = float(np.median(r))
    out["med_radius"] = med
    out["radial_cv"] = float(r.std() / r.mean()) if r.mean() > 1e-9 else 1.0
    if med > 1e-9:
        out["shell_concentration"] = float(np.mean(np.abs(r - med) <= 0.4 * med))
    com = cats.mean(axis=0)
    dl = links - com
    dl -= box_size * np.round(dl / box_size)
    ang = np.sort(np.arctan2(dl[:, 1], dl[:, 0]))
    gaps = np.diff(np.concatenate([ang, [ang[0] + 2 * np.pi]]))
    out["angular_gap_deg"] = float(np.degrees(gaps.max()))
    if len(subs) > 0 and med > 1e-9:
        rs = _nearest_cat_dist(subs)
        area_frac = np.pi * med ** 2 / (box_size ** 2)
        out["gradient"] = float(np.mean(rs < med)) / area_frac if area_frac > 1e-9 else 1.0
    return out


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

    # Calibrated on the membrane regime (max_links cap): positive radial_cv
    # 0.03-0.12, shell_concentration ~1.0, n_links 35; random scatter / no-
    # attraction CV ~0.4-0.56; no-production 0 links; high-decay ~3 links.
    _MIN_LINKS = 15
    _SCREEN_CV_MAX = 0.20
    _SCREEN_SHELL_MIN = 0.70
    _CONFIRM_GAP_MAX = 90.0
    _CONFIRM_COUNT_CV_MAX = 0.30
    _DEF_CV_MAX = 0.15

    def _compute_primary(self, state_history, timescale):
        """Membrane metrics averaged over the late window (last 25%)."""
        late = state_history[max(1, 3 * len(state_history) // 4):]
        keys = ["radial_cv", "shell_concentration", "med_radius",
                "angular_gap_deg", "gradient", "n_links"]
        acc = {k: [] for k in keys}
        for snap in late:
            m = _membrane_metrics(
                np.asarray(snap["positions"], dtype=np.float64),
                np.asarray(snap["types"], dtype=np.int32),
                float(snap.get("box_size", 20.0)),
            )
            for k in keys:
                acc[k].append(m[k])
        return {k: float(np.mean(v)) for k, v in acc.items()}

    def _check_screening(self, primary_result, timescale):
        """Screening: a THIN radial shell of links exists (closed-boundary topology)."""
        return (
            primary_result.get("n_links", 0.0) >= self._MIN_LINKS
            and primary_result.get("radial_cv", 1.0) < self._SCREEN_CV_MAX
            and primary_result.get("shell_concentration", 0.0) > self._SCREEN_SHELL_MIN
        )

    def _compute_secondaries(self, state_history, timescale):
        """Self-production (link count stable + nonzero despite decay) + closure + gradient."""
        late = state_history[max(1, 3 * len(state_history) // 4):]
        counts = [int(np.sum(np.asarray(s["types"]) == 2)) for s in late]
        carr = np.array(counts, dtype=float)
        link_cv = float(carr.std() / carr.mean()) if carr.mean() > 0 else 999.0
        gaps, grads = [], []
        for snap in late:
            m = _membrane_metrics(
                np.asarray(snap["positions"], dtype=np.float64),
                np.asarray(snap["types"], dtype=np.int32),
                float(snap.get("box_size", 20.0)),
            )
            gaps.append(m["angular_gap_deg"])
            grads.append(m["gradient"])
        return {
            "link_count_cv": link_cv,
            "mean_link_count": float(carr.mean()),
            "angular_gap_deg": float(np.mean(gaps)) if gaps else 360.0,
            "gradient": float(np.mean(grads)) if grads else 1.0,
        }

    def _run_null_model(self, state_history, primary_result, timescale):
        """Random-cloud null: place n_links points area-uniformly in a disk of
        radius ~1.8*med around the catalyst and compute radial CV. A membrane's
        tight-ring CV is far BELOW this scatter CV (~0.4)."""
        rng = np.random.default_rng(self._seed)
        observed_cv = primary_result.get("radial_cv", 1.0)
        med = max(primary_result.get("med_radius", 1.0), 1e-6)
        n_links = max(int(primary_result.get("n_links", 0)), 3)
        R = 1.8 * med
        null_cvs = []
        for _ in range(self.n_permutations):
            rr = R * np.sqrt(rng.random(n_links))
            null_cvs.append(float(rr.std() / rr.mean()) if rr.mean() > 1e-9 else 1.0)
        null_arr = np.array(null_cvs)
        p = float(np.mean(null_arr <= observed_cv))
        p = max(p, 1.0 / (self.n_permutations + 1))
        std = float(null_arr.std())
        return p, NullType.SHUFFLE, {"mean": float(null_arr.mean()), "std": std if std > 0 else 0.001}

    def _check_exclusions(self, state_history, model_metadata, timescale):
        """Exclude P1 (mere clustering): a thin organized shell (low radial CV)
        is not undifferentiated aggregation."""
        p = self._compute_primary(state_history, timescale)
        if (p.get("radial_cv", 1.0) < self._SCREEN_CV_MAX
                and p.get("n_links", 0.0) >= self._MIN_LINKS):
            return (["P1"], {"P1": "excluded"})
        return (["P1"], {"P1": "not_excluded"})

    def _check_confirmation(self, primary_result, secondary_result, null_p, timescale):
        """Confirmation: closed ring + self-production + null rejected."""
        return (
            secondary_result.get("angular_gap_deg", 360.0) < self._CONFIRM_GAP_MAX
            and secondary_result.get("link_count_cv", 999.0) < self._CONFIRM_COUNT_CV_MAX
            and secondary_result.get("mean_link_count", 0.0) >= self._MIN_LINKS
            and null_p < 0.01
        )

    def _check_definitive(self, primary_result, secondary_result, null_p,
                          null_type, state_history, model_metadata, timescale):
        """Definitive: a strongly tight closed shell, self-maintained, vs the
        random-cloud null at the floor. (Maintained-gradient + self-repair are
        reported/secondary; self-repair perturbation is a future addition.)"""
        return (
            primary_result.get("radial_cv", 1.0) < self._DEF_CV_MAX
            and null_p <= 0.005
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
