# T2c — OOD Validation Suite (instrument generalization)

**Seeds per system:** 3. Battery pointed at held-out systems via `analysis/battery_profile.profile_observation`. Recognition systems are independent re-implementations the detectors were never tuned on; novelty systems are emergent but out-of-catalog; null systems are non-emergent. Probes (reported last) expose known emergence-indicator coverage gaps and are excluded from the headline rates.

## Headline rates — raw (MATCH at any tier, current keystone)

| metric | value | reading |
|---|---|---|
| recognition recall | 0.833 | alt-model fires the correct pattern (MATCH) |
| recognition top-1 | 1.0 | correct pattern ranked #1 (incl. below-threshold) |
| novelty abstention | 0.583 | out-of-catalog emergence → EMERGENT-UNCLASSIFIED |
| null specificity | 1.0 | non-emergent → NO-EMERGENCE |
| **false-MATCH rate** | **0.208** | novelty/null wrongly claimed as a pattern (lower=better) |
| false-novel rate | 0.0 | null wrongly flagged emergent (lower=better) |

Counts: {"recognition": {"recognized": 15, "top1_only": 3}, "novelty": {"false_match": 5, "abstained": 7}, "null": {"correct": 12}}

## Headline rates — strict (MATCH requires ≥ confirmation tier)

Screening-tier detections are demoted to the emergence-driven verdict. This is a pure instrument-layer gate (no detector is changed) and shows how much raw false-MATCH is low-rigor screening noise.

| metric | value | reading |
|---|---|---|
| recognition recall | 0.833 | alt-model fires the correct pattern (MATCH) |
| recognition top-1 | 1.0 | correct pattern ranked #1 (incl. below-threshold) |
| novelty abstention | 1.0 | out-of-catalog emergence → EMERGENT-UNCLASSIFIED |
| null specificity | 1.0 | non-emergent → NO-EMERGENCE |
| **false-MATCH rate** | **0.0** | novelty/null wrongly claimed as a pattern (lower=better) |
| false-novel rate | 0.0 | null wrongly flagged emergent (lower=better) |

Counts: {"recognition": {"recognized": 15, "top1_only": 3}, "novelty": {"abstained": 12}, "null": {"correct": 12}}


## Recognition arm (held-out re-implementations)

| system | expect | verdicts (across seeds) | outcome | mean em |
|---|---|---|---|---|
| recog_p16_boolean_grn | P16 | MATCH×3  [→ P16] | recognized | 0.089 |
| recog_p20_fraction_threshold | P20 | MATCH×3  [→ P20] | recognized | 0.089 |
| recog_p23_el_farol | P23 | MATCH×3  [→ P23] | recognized | 0.089 |
| recog_p24_integral | P24 | NO-EMERGENCE×3 | top1_only | 0.089 |
| recog_p25_multibasin | P25 | MATCH×3  [→ P25] | recognized | 0.949 |
| recog_p26_threshold_unit | P26 | MATCH×3  [→ P26] | recognized | 0.089 |

## Novelty arm (out-of-catalog emergent)

| system | expect | verdicts (across seeds) | outcome | mean em |
|---|---|---|---|---|
| nov_dla | — | MATCH×2, EMERGENT-UNCLASSIFIED  [→ P18] | false_match | 0.822 |
| nov_keller_segel | — | EMERGENT-UNCLASSIFIED×3 | abstained | 1.0 |
| nov_active_nematic | — | EMERGENT-UNCLASSIFIED×3 | abstained | 0.895 |
| nov_eden | — | MATCH×3  [→ P18] | false_match | 0.887 |

## Null arm (non-emergent)

| system | expect | verdicts (across seeds) | outcome | mean em |
|---|---|---|---|---|
| null_spatial_noise | — | NO-EMERGENCE×3 | correct | 0.065 |
| null_random_walk | — | NO-EMERGENCE×3 | correct | 0.077 |
| null_uncoupled_phases | — | NO-EMERGENCE×3 | correct | 0.131 |
| null_frozen_noise | — | NO-EMERGENCE×3 | correct | 0.124 |

## Probe arm (coverage-limit demonstrations — not scored)

| system | expect | verdicts (across seeds) | outcome | mean em |
|---|---|---|---|---|
| probe_static_blob | — | EMERGENT-UNCLASSIFIED×3 | EMERGENT-UNCLASSIFIED | 0.623 |
| probe_percolation | — | EMERGENT-UNCLASSIFIED×3 | EMERGENT-UNCLASSIFIED | 0.686 |
| probe_phase_ordering | — | MATCH×3  [→ P3] | MATCH | 0.807 |
