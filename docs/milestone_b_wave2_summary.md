# Milestone B — Wave 2 Summary (P24 / P26 / P23)

**Date:** 2026-06-09 (Sprint 78 checkpoint)
**Wave 2 sprints:** 72–77
**Outcome:** 3/3 patterns to AT-DEPTH

---

## Counts

| Metric | Value |
|--------|-------|
| Implemented patterns | 25 |
| AT-DEPTH | 24 |
| GAP | 1 (P12 — documented finite-size measurement limitation) |
| Wave 2 patterns attempted | 3 |
| Wave 2 patterns to AT-DEPTH | 3 |

---

## Wave 2 Outcome Table

| Column | P24 Homeostatic regulation | P26 Stochastic resonance | P23 Anti-coordination |
|--------|---------------------------|--------------------------|----------------------|
| **Model file** | `epc/models/homeostasis.py` | `epc/models/stochastic_resonance.py` | `epc/models/minority_game.py` |
| **Detector file** | `epc/detectors/p24_homeostasis.py` | `epc/detectors/p26_stochastic_resonance.py` | `epc/detectors/p23_anticoordination.py` |
| **dim1 reproduction** | Ashby (1956) proportional homeostat: deviation ratio = 0.0027 (tolerance < 0.30) | Gammaitoni (1998) / Collins (1995): peak coherent response = 0.918, gain = 0.855 (> 0.05), decline = 0.811 (> 0.02), interior peak confirmed | Savit, Manuca & Riolo (1999) σ²/N vs α curve: interior minimum at α ≈ 0.32, σ²/N = 0.077 (below random baseline 0.25); symmetric phase σ²/N up to 1.45. All 3 tolerance checks PASS |
| **dim1 tolerance verdict** | PASS | PASS | PASS |
| **dim1 anchor strength** | **Internal threshold.** Ashby (1956) describes qualitative homeostatic regulation. The deviation-ratio tolerance (< 0.30) is an internally defined quantitative threshold, not a reproduction of a published figure or table. | **Published qualitative + internal thresholds.** The inverted-U (non-monotone peak in performance vs noise) is the canonical published SR signature (Gammaitoni 1998, Rev. Mod. Phys.). Specific numerical tolerances (gain > 0.05, decline > 0.02) are internal. | **Published quantitative.** The σ²/N vs α curve with interior minimum below random baseline is a well-characterized published quantitative result from Savit et al. (1999, Phys. Rev. Lett.). Tolerances directly reference the paper's predictions (below-baseline minimum, symmetric-phase excess). |
| **dim4 panel TNR** | 1.000 | 1.000 | 1.000 |
| **dim4 Cohen's d** | +inf | +inf | 14.504 |
| **AT-DEPTH** | Yes | Yes | Yes |
| **T1a adapter** | Yes — `extract_observation_bundle()`, `scalar_timeseries` format | Yes — `extract_observation_bundle()`, `noise_sweep` format | Yes — `extract_observation_bundle()`, `choice_timeseries` format |
| **T1b cross-model test** | Yes — `IntegralHomeostat` in `tests/test_homeostasis_p24_e2e.py` | Yes — `ThresholdUnit` in `tests/test_stochastic_resonance_p26_e2e.py` | Yes — `ElFarolBar` in `tests/test_minority_game_p23_e2e.py` |

---

## Content Prerequisites Added During Panels

| Pattern | Sprint | Prerequisite | Literature grounding |
|---------|--------|-------------|---------------------|
| P24 | 73 | Invariance flags only (perm=True, time=True); deviation integral is inherently order-invariant for constant dt | Ashby (1956): homeostatic deviation is a cumulative measure |
| P26 | 75 | **Inverted-U shape** required at screening: gain > 0.02, interior peak, rise-then-fall | Gammaitoni (1998): SR is defined by non-monotone noise-performance relationship |
| P23 | 77 | **Non-degenerate variance** (σ² > 0) AND **variance below random baseline** (σ²/N < p̂(1−p̂)) at confirmation | Savit, Manuca & Riolo (1999): efficient-phase anti-coordination is characterized by sub-random fluctuations |

---

## dim1 Anchor-Strength Assessment

Wave 2 anchor quality is mixed:

- **P23 (strong):** Published quantitative curve (Savit et al. 1999) with specific numeric predictions. All three tolerance checks reference paper-derived quantities. This is the strongest dim1 anchor of the Wave 2 set.
- **P26 (moderate):** The inverted-U shape is a published qualitative signature (Gammaitoni 1998), and any SR system should exhibit it. However, the specific coherent-response thresholds are internal. The anchor validates the phenomenon but not a specific published measurement.
- **P24 (weak):** Ashby (1956) is the conceptual source but does not provide quantitative predictions for computational reproduction. The deviation-ratio tolerance is entirely internal. This is the weakest dim1 anchor in the catalog — functionally a "does the model work as designed" check rather than a literature reproduction.

---

## Open Carry-Forwards

| ID | Pattern | Description | Priority | Since |
|----|---------|-------------|----------|-------|
| C-p7-time-shuffled-fp | P7 | `time_shuffled` FP at screening; each frame preserves lane structure independently of temporal order | Low | Sprint 66 |
| C-p19-bias-zero-chance-alignment | P19 | 1/5 bias_zero Class C FP at confirmation by chance alignment (4% rate) | Low | Sprint 70 |
| P12 dim1 | P12 | Documented finite-size measurement limitation (λ ∝ √M scaling); accepted after 4 attempts (Sprints 54/58/59/63) | Accepted | Sprint 54 |

---

## Wave 2 Sprint History

| Sprint | Pattern | Work | Outcome |
|--------|---------|------|---------|
| 72 | P24 | dim1–dim3: Ashby homeostat reproduction, 20-seed campaign, methods note | PASS (all 3 dims) |
| 73 | P24 | dim4: Phase-2a panel v1.2 | PASS (TNR=1.000, d=+inf) → **AT-DEPTH** |
| 74 | P26 | dim1–dim3: Gammaitoni/Collins reproduction, 20-seed campaign, methods note | PASS (all 3 dims) |
| 75 | P26 | dim4: Phase-2a panel v1.2 + inverted-U prerequisite | PASS (TNR=1.000, d=+inf) → **AT-DEPTH** |
| 76 | P23 | dim1–dim3: Savit curve reproduction, 25-seed campaign, methods note, T1a/T1b OOD-readiness | PASS (all 3 dims) |
| 77 | P23 | dim4: Phase-2a panel v1.2 + below-baseline prerequisite | PASS (TNR=1.000, d=14.504) → **AT-DEPTH** |
| 78 | — | Wave 2 summary (this document) | Non-blocking checkpoint |
