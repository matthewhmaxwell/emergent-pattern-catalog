"""Per-detector primary-metric invariance flags (Phase-2a spec v1.2).

These flags are panel-harness metadata: when ``primary_metric_permutation_invariant``
is True, the harness skips Class A's ``permutation_shuffled_positive`` substrate
for that detector (it is degenerate-by-construction); when
``primary_metric_time_shuffle_invariant`` is True, ``time_shuffled_positive`` is
skipped.

Centralized declaration (Sprint 34 Approach 2) keeps panel metadata out of the
detector modules — the flags are consumed only by ``epc/phase2a/panel.py`` and
have nothing to do with the detector's own contract. Greppable from one place
when audits or v1.x revisions retune the table.

Authoritative initial values from v1.2 spec §"Change 2". The
``primary_metric`` column is informational; the flags are what the harness
reads.
"""

from __future__ import annotations

from dataclasses import dataclass
from typing import Dict


@dataclass(frozen=True)
class InvarianceFlags:
    permutation_invariant: bool
    time_shuffle_invariant: bool
    primary_metric: str           # short label, informational
    rationale: str = ""           # optional one-line rationale


DETECTOR_INVARIANCE_FLAGS: Dict[str, InvarianceFlags] = {
    # Sprint 49 batch update — flags derived from Sprints 39-48 empirical Class A FP evidence.
    "P1":  InvarianceFlags(False, True,  "Moran's I + same-type neighbor fraction",
                            rationale="Spatial autocorrelation depends on adjacency (perm_inv=False). "
                                      "P1's Moran's I is computed per-frame; each segregated frame has "
                                      "high I regardless of temporal order → time_shuffled is "
                                      "degenerate-by-construction. C-p1-time-shuffle-fp (Sprint 43)."),
    "P2":  InvarianceFlags(True,  False, "two_phase_score (spatial density ratio)",
                            rationale="P2's two_phase_score is a spatial-distribution statistic "
                                      "invariant under cell-index permutation (it measures cluster area "
                                      "fraction, not cell identity) → permutation_shuffled is "
                                      "degenerate-by-construction. C-p2-perm-shuffled-fp (Sprint 46)."),
    "P3":  InvarianceFlags(False, True,  "radial FFT peak on final field snapshot",
                            rationale="P3 computes spatial FFT per frame; each Gray-Scott frame "
                                      "contains the full Turing pattern regardless of temporal "
                                      "ordering — time_shuffled is degenerate-by-construction."),
    "P5":  InvarianceFlags(True,  True,  "heading order parameter |⟨e^iθ⟩|",
                            rationale="perm_inv=True: mean of unit vectors is invariant under particle "
                                      "relabelling (aggregate over headings). time_shuffle_inv=True: "
                                      "each Vicsek frame has high φ independent of temporal order — "
                                      "time_shuffled is degenerate-by-construction. "
                                      "C-p5-time-shuffle-fp (Sprint 46)."),
    "P6":  InvarianceFlags(True,  True,  "group angular momentum |L|",
                            rationale="perm_inv=True: |L| is a sum over particles, invariant under "
                                      "particle relabelling. time_shuffle_inv=True: each milled frame "
                                      "has high |L| regardless of temporal order → time_shuffled is "
                                      "degenerate-by-construction. C-p6-time-shuffle-fp (Sprint 46)."),
    "P8":  InvarianceFlags(True,  True,  "stopped_fraction (time-averaged jamming density)",
                            rationale="P8's stopped_fraction is a time-average statistic invariant "
                                      "under both cell-index permutation (cells anonymous: NaSch is "
                                      "translation-invariant) and temporal shuffle (the mean is "
                                      "preserved by any reordering) → both Class A substrates are "
                                      "degenerate-by-construction. "
                                      "C-p8-perm-shuffled-fp + C-p8-time-shuffle-fp (Sprint 47)."),
    "P9":  InvarianceFlags(True,  True,  "Kuramoto order parameter r",
                            rationale="Aggregate over phases, final-state."),
    "P10": InvarianceFlags(False, False, "local-coherence partitioning",
                            rationale="Coherent vs incoherent regions are spatial. "
                                      "C-p10-perm-shuffled-fp (Sprint 47): the FP arises from "
                                      "catalog-adapter binarization preserving bimodal phase "
                                      "structure — this is an adapter artifact, NOT mathematical "
                                      "permutation invariance of local_r. Flags intentionally "
                                      "unchanged (Sprint 49 decision: wrong fix to auto-flag)."),
    "P11": InvarianceFlags(False, True,  "population oscillation period / amplitude",
                            rationale="Spatial well-mixed; trajectory order matters."),
    "P12": InvarianceFlags(False, False, "spiral morphology / species lag",
                            rationale="Spatial + temporal structure."),
    "P13": InvarianceFlags(False, False, "wavefront propagation speed",
                            rationale="Spatial + temporal structure."),
    "P14": InvarianceFlags(True,  True,  "avalanche-size power-law exponent",
                            rationale="Aggregate over event list, order-free."),
    "P15": InvarianceFlags(False, False, "TE across collisions + functional reproducibility",
                            rationale="Information transfer depends on adjacency + ordering."),
    "P17": InvarianceFlags(False, False, "collective chemotactic index",
                            rationale="Direction + trajectory matter."),
    "P18": InvarianceFlags(True,  True,  "convergence to consensus (max f_k)",
                            rationale="Aggregate fraction, final-state."),
    "P19": InvarianceFlags(False, False, "influence-asymmetry TE ratio",
                            rationale="Spatial + temporal information flow."),
    "P21": InvarianceFlags(True,  False, "dip test on opinion distribution",
                            rationale="perm_inv=True: Hartigan dip test operates on the sorted opinion "
                                      "value histogram — permuting the N=100 opinion values from a "
                                      "bimodal HK positive preserves the bimodal shape exactly → "
                                      "permutation_shuffled is degenerate-by-construction. "
                                      "time_shuffle_inv=False: shuffling temporal order mixes early "
                                      "unimodal pre-convergence steps with late bimodal steps; the HK "
                                      "trajectory is NOT time-order-independent, so time_shuffled "
                                      "remains a meaningful test (not skipped). "
                                      "C-p21-perm-shuffled-fp (Sprint 48); C-p21-time-shuffled-fp "
                                      "remains open (HK convergence-timing issue, not invariance)."),
    "P22": InvarianceFlags(False, False, "cascade size / propagation speed",
                            rationale="Network-temporal structure."),
    # P27 time-shuffle flag is PROVISIONAL — see C-p27-time-shuffle-invariance carry-forward.
    # If Sprint 35+ P27 panel run reveals the assumption is wrong, change the flag and re-run.
    "P27": InvarianceFlags(False, True,  "cooperation fraction time-series",
                            rationale="Provisional per v1.2 spec Change 2; validate in Sprint 35+. C-p27-time-shuffle-invariance."),
    "P28": InvarianceFlags(True,  True,  "wealth-distribution Gini / cluster index",
                            rationale="Distributional + time-aggregated."),
    "P31": InvarianceFlags(False, False, "DG monotonicity / avg_wandering_range",
                            rationale="Sequence ordering is the signal."),
}


# --- Lookup helpers ---------------------------------------------------------

_DEFAULT = InvarianceFlags(False, False, "(unknown)", rationale="default fallback")


def get_flags(pattern_id: str) -> InvarianceFlags:
    """Return the invariance flags for ``pattern_id``. Defaults to (False, False)
    for unknown patterns — the conservative choice that runs every substrate.
    """
    return DETECTOR_INVARIANCE_FLAGS.get(pattern_id, _DEFAULT)


def is_permutation_invariant(pattern_id: str) -> bool:
    return get_flags(pattern_id).permutation_invariant


def is_time_shuffle_invariant(pattern_id: str) -> bool:
    return get_flags(pattern_id).time_shuffle_invariant
