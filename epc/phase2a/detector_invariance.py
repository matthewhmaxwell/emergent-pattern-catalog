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
    "P1":  InvarianceFlags(False, False, "Moran's I + same-type neighbor fraction",
                            rationale="Spatial autocorrelation depends on adjacency."),
    "P3":  InvarianceFlags(False, True,  "radial FFT peak on final field snapshot",
                            rationale="P3 computes spatial FFT per frame; each Gray-Scott frame "
                                      "contains the full Turing pattern regardless of temporal "
                                      "ordering — time_shuffled is degenerate-by-construction."),
    "P5":  InvarianceFlags(True,  False, "heading order parameter |⟨e^iθ⟩|",
                            rationale="Aggregate over headings; final-state metric."),
    "P6":  InvarianceFlags(False, False, "group rotational dynamics",
                            rationale="Trajectory-shape detector."),
    "P9":  InvarianceFlags(True,  True,  "Kuramoto order parameter r",
                            rationale="Aggregate over phases, final-state."),
    "P10": InvarianceFlags(False, False, "local-coherence partitioning",
                            rationale="Coherent vs incoherent regions are spatial."),
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
    "P21": InvarianceFlags(False, False, "dip test on opinion distribution",
                            rationale="Dip test is distributional but permuting opinion values "
                                      "from the canonical HK positive preserves the bimodal "
                                      "cluster structure → permutation_shuffled and time_shuffled "
                                      "are expected FPs (carry-forwards C-p21-perm-shuffled-fp, "
                                      "C-p21-time-shuffled-fp). Sprint 48: flags corrected from "
                                      "(True,True) per brief Notes."),
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
