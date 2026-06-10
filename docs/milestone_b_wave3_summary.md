# Milestone B — Wave 3 Summary (P16 / P25 / P20)

**Date:** 2026-06-10 (Sprint 85 checkpoint)
**Wave 3 sprints:** 79–84
**Outcome:** 3/3 patterns to AT-DEPTH

---

## Counts

| Metric | Value |
|--------|-------|
| Implemented patterns | 28 |
| AT-DEPTH | 27 |
| GAP | 1 (P12 — documented finite-size measurement limitation) |
| Wave 3 patterns attempted | 3 |
| Wave 3 patterns to AT-DEPTH | 3 |

---

## Wave 3 Outcome Table

| Column | P16 Associative memory | P25 Canalized restoration | P20 Quorum sensing |
|--------|------------------------|---------------------------|-------------------|
| **Model file** | `epc/models/hopfield_network.py` | `epc/models/canalization.py` | `epc/models/quorum_sensing.py` |
| **Detector file** | `epc/detectors/p16_associative_memory.py` | `epc/detectors/p25_canalization.py` | `epc/detectors/p20_quorum_sensing.py` |
| **dim1 reproduction** | Amit-Gutfreund-Sompolinsky (1985) storage capacity: α sweep at N=500, 10 seeds/point. Perfect retrieval at α ≤ 0.10. Transition midpoint α ≈ 0.173 (tolerance [0.10, 0.20] PASS). High-load overlap 0.373 < 0.50 (PASS). All 3 tolerance checks PASS | Waddington (1957) canalized landscape: convergence variance ratio ≈ 0.0 (tolerance < 0.10, PASS), basin volume = 1.0 (tolerance ≥ 0.80, PASS). DEFINITIVE | Waters & Bassler (2005) quorum sensing: step-function R² = 1.000 (tolerance > 0.9, PASS), hysteresis width = 1.190 (tolerance > 0.1, PASS). DEFINITIVE. All 3 tolerance checks PASS |
| **dim1 tolerance verdict** | PASS | PASS | PASS |
| **dim1 anchor strength** | **Published quantitative.** AGS (1985) storage capacity curve with specific α_c prediction (~0.138 at N→∞). Finite-size shifted midpoint (0.173 at N=500) is well within expected finite-N range. Tolerances reference published predictions. | **Internal threshold.** Waddington (1957) describes qualitative canalization. The convergence-variance-ratio tolerance (< 0.10) and basin-volume threshold (≥ 0.80) are internally defined. Anchor validates the phenomenon, not a specific published measurement. | **Published qualitative + internal thresholds.** Waters & Bassler (2005) describe sharp threshold activation + hysteresis as canonical quorum-sensing signatures. Specific R² and hysteresis-width thresholds are internal. |
| **dim4 panel TNR** | 1.000 | 1.000 | 1.000 |
| **dim4 Cohen's d** | +inf | +inf | +inf |
| **AT-DEPTH** | Yes | Yes | Yes |
| **T1a adapter** | Yes — `extract_observation_bundle()`, `state_vector` format | Yes — `extract_observation_bundle()`, `canalization_bundle` format | Yes — `extract_observation_bundle()`, `density_sweep_timeseries` format |
| **T1b cross-model test** | Yes — `BooleanGRN` | Yes — `MultiBasinGRN` | Yes — `FractionThreshold` |

---

## Content Prerequisites Added During Panels

| Pattern | Sprint | Prerequisite | Literature grounding |
|---------|--------|-------------|---------------------|
| P16 | 80 | ≥2 distinct selectively-retrievable stored patterns required at confirmation | Hopfield (1982): content-addressable memory is defined by multi-pattern retrieval; single-attractor collapse is not associative memory |
| P25 | 82 | basin_volume ≥ 0.5 at screening (wide-IC convergence) | Waddington (1957): equifinality requires convergence from diverse initial conditions |
| P20 | 84 | Invariance flags: perm_inv=False, time_shuffle_inv=True (tag-based grouping preserves signal) | Waters & Bassler (2005): density-label permutation destroys step-function fit; equilibrium curves are reconstructed from density/direction tags |

---

## dim1 Anchor-Strength Assessment

Wave 3 anchor quality is mixed:

- **P16 (strong):** Published quantitative curve (AGS 1985) with specific α_c prediction. The storage-capacity transition is one of the most precisely characterized results in neural network theory. Finite-size shift is well-understood and within expected range.
- **P20 (moderate):** Sharp threshold activation and hysteresis are published qualitative signatures of quorum sensing (Waters & Bassler 2005), but the specific R² and hysteresis-width thresholds are internal. Validates the phenomenon class.
- **P25 (weak):** Waddington (1957) is the conceptual source but provides no quantitative predictions for computational reproduction. All tolerances are internal — functionally a "does the model work as designed" check.

---

## Open Carry-Forwards

| ID | Pattern | Description | Priority | Since |
|----|---------|-------------|----------|-------|
| C-p7-time-shuffled-fp | P7 | `time_shuffled` FP at screening; each frame preserves lane structure independently of temporal order | Low | Sprint 66 |
| C-p19-bias-zero-chance-alignment | P19 | 1/5 bias_zero Class C FP at confirmation by chance alignment (4% rate) | Low | Sprint 70 |
| P12 dim1 | P12 | Documented finite-size measurement limitation (λ ∝ √M scaling); accepted after 4 attempts (Sprints 54/58/59/63) | Accepted | Sprint 54 |

---

## Wave 3 Sprint History

| Sprint | Pattern | Work | Outcome |
|--------|---------|------|---------|
| 79 | P16 | dim1–dim3: AGS (1985) storage-capacity reproduction, 20-seed campaign, methods note, T1a/T1b OOD-readiness | PASS (all 3 dims) |
| 80 | P16 | dim4: Phase-2a panel v1.2 + multi-pattern prerequisite | PASS (TNR=1.000, d=+inf) → **AT-DEPTH** |
| 81 | P25 | dim1–dim3: Waddington canalized landscape, 20-seed campaign, methods note, T1a/T1b OOD-readiness | PASS (all 3 dims) |
| 82 | P25 | dim4: Phase-2a panel v1.2 + basin-volume prerequisite | PASS (TNR=1.000, d=+inf) → **AT-DEPTH** |
| 83 | P20 | dim1–dim3: Waters-Bassler quorum sensing, 20-seed campaign, methods note, T1a/T1b OOD-readiness | PASS (all 3 dims) |
| 84 | P20 | dim4: Phase-2a panel v1.2 + invariance flags | PASS (TNR=1.000, d=+inf) → **AT-DEPTH** |
| 85 | — | Wave 3 summary (this document) | Non-blocking checkpoint |

---

## Cumulative Milestone B Progress (Waves 1–3)

| Wave | Patterns | Sprints | Outcome |
|------|----------|---------|---------|
| Wave 1 | P7, P17, P19 | 65–70 | 3/3 AT-DEPTH |
| Wave 2 | P24, P26, P23 | 72–77 | 3/3 AT-DEPTH |
| Wave 3 | P16, P25, P20 | 79–84 | 3/3 AT-DEPTH |
| **Total (Waves 1–3)** | **9 patterns** | **65–84** | **9/9 AT-DEPTH** |
| Wave 4 (remaining) | P4, P29, P32, P30 | 86+ | Pending — final 4 patterns → 32/32 |

Milestone B target: 13 new patterns (P4, P7, P16, P17, P19, P20, P23, P24, P25, P26, P29, P30, P32).
Waves 1–3 complete: 9 of 13 done. Wave 4 = final 4 patterns remaining for 32/32 catalog coverage.
