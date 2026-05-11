# Phase-2a Standard Negative Panel — Specification v1.2

This is a revision of v1.1 (`docs/phase2a_panel_spec.md`, dated Sprint 31). v1.2 supersedes v1.1. The v1.1 file is archived to `docs/archive/phase2a_panel_spec_v1_1.md`.

The Sprint 32 P9 Kuramoto result and the Sprint 33 P14 BTW result both surfaced the same failure mode: two substrates in Class A (`permutation_shuffled_positive` and `time_shuffled_positive`) preserve aggregate distributions over the canonical positive, so detectors whose primary metric is order-invariant correctly fire on them — making the substrates degenerate tests for those detectors.

v1.2 closes carry-forward C-class-a-permutation-degenerate by allowing detectors to declare invariance properties of their primary metric. The harness skips degenerate-by-construction substrates accordingly.

---

## What changed and why

### Change 1 — Detectors declare primary-metric invariance flags

**v1.1 problem:** Class A contained `permutation_shuffled_positive` (shuffles cells of the canonical positive's final state) and `time_shuffled_positive` (shuffles timesteps of the canonical positive's trajectory). Both substrates preserve aggregate statistics. Detectors whose primary metric is an aggregate statistic (Kuramoto r, sandpile power-law exponent, voter consensus fraction, opinion-distribution dip test) fire on them — *correctly*. The substrate is degenerate-by-construction for those detectors. Forcing a TNR computation on degenerate substrates produced misleading PARTIAL verdicts in Sprint 32 (P9) and Sprint 33 (P14).

**v1.2 fix:** Each detector declares two invariance flags in its metadata:

- `primary_metric_permutation_invariant: bool` — True when the detector's primary metric is invariant to spatial permutations of cell positions. Examples: Kuramoto r, voter consensus fraction max f_k, sandpile avalanche-size power-law exponent, opinion-distribution dip test.
- `primary_metric_time_shuffle_invariant: bool` — True when the detector's primary metric is invariant to temporal shuffling of timesteps. Examples: any detector that reduces a trajectory to a final-state statistic or an unordered set of events (sandpile avalanche-size distribution).

When the harness runs the panel against a detector:

- If `primary_metric_permutation_invariant` is True, the harness skips `permutation_shuffled_positive` and excludes it from TNR computation. The substrate appears in the panel JSON output with `verdict: "SKIPPED-degenerate-by-construction"` and a `skip_reason` field.
- If `primary_metric_time_shuffle_invariant` is True, the harness skips `time_shuffled_positive` similarly.
- Per-class TNR is computed over substrates that actually ran, not over the original Class A size. The panel JSON records both `class_a_size_total` and `class_a_size_evaluated`.

Defaults are `False` for both flags. If a detector doesn't declare them, both substrates run, and if the detector turns out to be invariant under one or both, a degenerate-by-construction FAIL will surface — which is the correct signal to add the flag.

### Change 2 — Documented detector-flag assignments (initial)

The catalog enumerates the initial flag assignments for currently-implemented detectors. This is the authoritative source for harness configuration; detector cards should reference this section.

| Pattern | Detector primary metric | `permutation_invariant` | `time_shuffle_invariant` | Rationale |
|---|---|---|---|---|
| P1 (Schelling) | Moran's I + same-type neighbor fraction | False | False | Spatial autocorrelation depends on adjacency. |
| P3 (Gray-Scott) | Spot/stripe morphology metrics | False | False | Pattern formation requires spatial coherence. |
| P5 (Vicsek) | Heading order parameter φ = |⟨e^iθ⟩| over headings | True | False | Aggregate over headings; final-state metric, not trajectory. |
| P6 (D'Orsogna mill) | Group rotational dynamics | False | False | Trajectory-shape detector. |
| P9 (Kuramoto) | Order parameter r = |⟨e^iφ⟩| | True | True | Aggregate over phases, final-state. |
| P10 (chimera) | Local coherence partitioning | False | False | Spatial structure required (coherent vs incoherent regions). |
| P11 (LV) | Population oscillation period / amplitude | False | True | Time-series shape matters (spatial doesn't for well-mixed); but trajectory order matters. |
| P12 (RPS) | Spiral morphology / species lag | False | False | Spatial + temporal structure. |
| P13 (GH) | Wavefront propagation speed | False | False | Spatial + temporal structure. |
| P14 (BTW) | Avalanche-size power-law exponent | True | True | Aggregate over event list, order-free. |
| P15 (GoL) | TE across collisions + functional reproducibility | False | False | Information transfer depends on adjacency + ordering. |
| P17 (Berdahl) | Collective chemotactic index | False | False | Direction + trajectory matter. |
| P18 (voter) | Convergence to consensus (max f_k) | True | True | Aggregate fraction, final-state. |
| P19 (Couzin) | Influence asymmetry TE ratio | False | False | Spatial + temporal information flow. |
| P21 (HK) | Dip test on opinion distribution | True | True | Distributional, final-state. |
| P22 (SIR) | Cascade size / propagation speed | False | False | Network-temporal structure. |
| P27 (Nowak-May) | Cooperation fraction time-series | False | True | Spatial structure matters; temporal order also matters for stability — but flag as time_shuffle_invariant pending detector audit. |
| P28 (Yard-Sale) | Wealth-distribution Gini / cluster index | True | True | Distributional + time-aggregated. |
| P31 (Zhang) | DG monotonicity / avg_wandering_range | False | False | Sequence ordering is the signal. |

This table is illustrative for the implementation; the *authoritative* source is the per-detector code declaration. Harness reads flags from detector metadata, not from this table.

**P27 caveat:** the flag for P27 is provisional. If Sprint 35's P27 panel run reveals the time-shuffle invariance assumption is wrong, change the flag and re-run — same `do not modify the detector to make it pass` rule still applies.

### Change 3 — Verdict labels handle skipped substrates

The verdict labels (PASS / PASS-with-weakness / PARTIAL / FAIL) are unchanged in v1.2. But the TNR denominators may now be smaller when substrates are skipped. The PASS criterion (≥95% TNR overall, soft expectation ≥90% per class with ≥5 evaluated substrates) operates on *evaluated* substrates only.

If skipping reduces Class A's evaluated size below 5, per-class TNR for Class A is reported as `class_a_tnr_advisory` rather than gating PASS — same rule as v1.1's small-class handling.

### Change 4 — Ratify P15 canonical positive

Sprint 33 discovered mid-sprint that P15's canonical positive needed to switch from R-pentomino to dense-random GoL (init_density=0.37) for the detector to fire. R-pentomino is too sparse to produce the structural-diversity signal P15's screening relies on.

v1.2 ratifies this change. The P15 canonical positive of record is now **dense-random GoL at init_density=0.37**. R-pentomino remains a valid positive for *qualitative* P15 demonstration (it produces gliders, blocks, and blinkers) but is not the panel's canonical positive.

This change is documented in this spec and should be reflected in `REPLICATION_NOTES.md`'s P15 section with rationale: R-pentomino's small initial activation produces a high-variance state trajectory where the structural-diversity metric is dominated by noise during the early activation transient. Dense-random IC produces a stable high-activity GoL with diverse structures consistently across seeds.

### Change 5 — Carry-forward changes

- **C-class-a-permutation-degenerate**: CLOSED by Change 1.
- **C-p14-class-c-borderline** (Sprint 33): remains OPEN, low priority, separate from v1.2.
- **C-p27-time-shuffle-invariance** (NEW): the P27 flag in Change 2 is provisional and must be validated in Sprint 35+.

---

## Composition (v1.2)

Same as v1.1 for Classes B and C. Class A composition is unchanged at 10 substrates, but per-detector evaluation may skip 1–2 of them based on the invariance flags.

---

## PASS criterion (v1.2)

Unchanged from v1.1. PASS / PASS-with-weakness / PARTIAL / FAIL with ≥95% TNR overall and ≥1.0 Cohen's d.

---

## Harness output (v1.2)

Schema additions:

- `panel_version: "1.2"`.
- `class_a_size_total: 10` and `class_a_size_evaluated: N` (N ≤ 10).
- Substrates that were skipped appear with `verdict: "SKIPPED-degenerate-by-construction"` and `skip_reason: "primary_metric_permutation_invariant" | "primary_metric_time_shuffle_invariant"`.

---

## Migration path

1. Sprint 34 (this spec, code-led application): land the spec, archive v1.1, update harness to read detector invariance flags, add the flags to existing detectors per the table above. NO new panel runs.
2. Sprint 35 (code-led): re-run P9 and P14 under v1.2. Both expected to PASS. Then begin lattice_2d batch (3–4 patterns per sprint).
3. Subsequent sprints: continue the batch across all PARTIAL patterns until all dim 4 cells are PASS or PASS-with-weakness.
