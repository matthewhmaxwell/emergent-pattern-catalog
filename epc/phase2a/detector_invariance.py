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
    "P7":  InvarianceFlags(True,  False, "lane order parameter (lateral directional segregation)",
                            rationale="perm_inv=True: lane order parameter φ_lane is a spatial-strip "
                                      "statistic (Nowak & Schadschneider 2012) — it bins agents by "
                                      "y-position and measures directional purity per strip; permuting "
                                      "agent indices does not change strip membership or label counts → "
                                      "permutation_shuffled is degenerate-by-construction. "
                                      "time_shuffle_inv=False: P7 confirmation requires temporal "
                                      "stability (CV < 0.3 over measurement window) and definitive "
                                      "requires throughput gain vs early transient; shuffling temporal "
                                      "order destroys both metrics → time_shuffled is a meaningful test."),
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
    "P17": InvarianceFlags(True,  False, "collective chemotactic index",
                            rationale="perm_inv=True: group CI is the distance from the group "
                                      "center-of-mass to the field peak — permuting agent indices "
                                      "does not change the CoM → permutation_shuffled is "
                                      "degenerate-by-construction. time_shuffle_inv=False: "
                                      "gradient climbing is a temporal trajectory; shuffling "
                                      "temporal order destroys the directional approach signal."),
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
    "P24": InvarianceFlags(True,  True,  "deviation integral / recovery trajectory",
                            rationale="perm_inv=True: scalar time series has a single regulated "
                                      "variable — permuting 'agents' is N/A (degenerate-by-construction "
                                      "for a single scalar). time_shuffle_inv=True: the primary metric "
                                      "(deviation integral ∫|x−sp|dt) is a sum of absolute deviations "
                                      "weighted by constant dt — invariant under temporal reordering. "
                                      "The growth_ratio gates only screening; post-screening metrics "
                                      "are order-invariant → time_shuffled is degenerate-by-construction. "
                                      "Sprint 73 panel run confirmed: time_shuffled FP at definitive."),
    # P27 time-shuffle flag is PROVISIONAL — see C-p27-time-shuffle-invariance carry-forward.
    # If Sprint 35+ P27 panel run reveals the assumption is wrong, change the flag and re-run.
    "P27": InvarianceFlags(False, True,  "cooperation fraction time-series",
                            rationale="Provisional per v1.2 spec Change 2; validate in Sprint 35+. C-p27-time-shuffle-invariance."),
    "P28": InvarianceFlags(True,  True,  "wealth-distribution Gini / cluster index",
                            rationale="Distributional + time-aggregated."),
    "P31": InvarianceFlags(False, False, "DG monotonicity / avg_wandering_range",
                            rationale="Sequence ordering is the signal."),
    "P16": InvarianceFlags(True,  True,  "completion accuracy (overlap with stored patterns)",
                            rationale="perm_inv=True: completion accuracy = (1/N)|Σ ξ_i s_i| "
                                      "is invariant under a consistent permutation of neuron "
                                      "indices applied to both state and stored patterns "
                                      "(Hopfield 1982, §II) → permutation_shuffled is "
                                      "degenerate-by-construction. time_shuffle_inv=True: "
                                      "P16 operates on fixed-point (converged) states; the "
                                      "converged state appears in multiple history entries per "
                                      "trial (all post-convergence steps identical), so "
                                      "time_shuffled preserves the signal → degenerate-by-"
                                      "construction."),
    "P25": InvarianceFlags(True,  True,  "convergence variance ratio Var(finals)/Var(ICs)",
                            rationale="perm_inv=True: convergence variance ratio is computed "
                                      "over the IC ensemble — permuting trial indices does not "
                                      "change per-dimension variances → permutation_shuffled is "
                                      "degenerate-by-construction. time_shuffle_inv=True: "
                                      "the observation bundle sorts records by (trial, step); "
                                      "ICs and finals are the first/last records per trial. "
                                      "Shuffling temporal order of list entries preserves step "
                                      "labels, so sorted extraction recovers original ICs/finals → "
                                      "time_shuffled is degenerate-by-construction."),
    "P23": InvarianceFlags(True,  True,  "σ²/N of attendance + lag-1 autocorrelation",
                            rationale="perm_inv=True: σ²/N is computed over attendance counts, "
                                      "which are scalar sums over agents — permuting agent indices "
                                      "does not change attendance → permutation_shuffled is "
                                      "degenerate-by-construction. time_shuffle_inv=True: "
                                      "σ²/N is the primary confirmation-level signal and is "
                                      "preserved by time shuffling (it is a distribution statistic). "
                                      "Lag-1 AC is a secondary signal too weak to reliably gate "
                                      "confirmation at the efficient phase (AC ≈ −0.01 to −0.05, "
                                      "p > 0.01 for most seeds). Sprint 77 panel run confirmed: "
                                      "time_shuffled FP at confirmation via σ²/N alone. Same "
                                      "pattern as P5 (Sprint 46), P1 (Sprint 43)."),
    "P29": InvarianceFlags(False, False, "weight-distance correlation + network efficiency",
                            rationale="perm_inv=False: weight-distance correlation depends on "
                                      "node-edge structure — permuting node indices changes edge-to-"
                                      "node correspondence, altering the correlation metric. "
                                      "time_shuffle_inv=False: P29 reads the FINAL snapshot's "
                                      "edge weights, which are the cumulative result of reinforcement "
                                      "dynamics — shuffling temporal order changes which snapshot is "
                                      "'final', potentially presenting an early-stage unreinforced "
                                      "network → correlation drops. Sprint 89 panel."),
    "P4":  InvarianceFlags(True,  True,  "exclusivity index + boundary persistence",
                            rationale="perm_inv=True: exclusivity index = mean per-cell "
                                      "max(visits_i)/total_visits is invariant under consistent "
                                      "agent-index relabelling (same max, same total) → "
                                      "permutation_shuffled is degenerate-by-construction. "
                                      "time_shuffle_inv=True: detector reads CUMULATIVE occupancy "
                                      "from each snapshot (monotonically increasing); boundary "
                                      "persistence compares winner-take-all ownership between "
                                      "early and late cumulative snapshots. After time shuffling, "
                                      "a late cumulative snapshot still contains the same ownership "
                                      "structure → persistence is preserved. Sprint 87 panel "
                                      "confirmed: time_shuffled FP at confirmation (d=1.2) when "
                                      "flag was False. Changed to True per cumulative-occupancy "
                                      "implementation semantics."),
    "P20": InvarianceFlags(False, True,  "step-function R² on density-response curve",
                            rationale="perm_inv=False: P20's sharpness metric is computed over "
                                      "a density sweep — permuting fraction_on values across "
                                      "density levels destroys the step-function fit → "
                                      "permutation_shuffled is a meaningful test. "
                                      "time_shuffle_inv=True: P20's equilibrium curves are "
                                      "reconstructed from density_idx and sweep_direction tags "
                                      "attached to each history dict — shuffling the list order "
                                      "does not destroy the signal because grouping is tag-based, "
                                      "not position-based. Sprint 84 panel run confirmed: "
                                      "time_shuffled FP at definitive."),
    "P26": InvarianceFlags(True,  True,  "coherent response |⟨x·signal⟩| over noise sweep",
                            rationale="perm_inv=True: noise-sweep timeseries has a single scalar "
                                      "variable x — permuting 'agents' is N/A (degenerate-by-construction "
                                      "for a single scalar, same as P24). time_shuffle_inv=True: "
                                      "time_shuffled preserves per-noise-level (x, signal) pairs "
                                      "(noise_level_idx travels with each dict), and coherent response "
                                      "|mean(x·signal)| is order-invariant within each group → "
                                      "the inverted-U curve is identical after shuffling. "
                                      "Sprint 75 panel run confirmed: time_shuffled FP at definitive."),
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
