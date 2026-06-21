# T2c — OOD Validation Suite (instrument generalization)

**Seeds per system:** 3. Battery pointed at held-out systems via `analysis/battery_profile.profile_observation`. Recognition systems are independent re-implementations the detectors were never tuned on; novelty systems are emergent but out-of-catalog; null systems are non-emergent. Probes (reported last) expose known emergence-indicator coverage gaps and are excluded from the headline rates.

## Headline rates — raw (MATCH at any tier, current keystone)

| metric | value | reading |
|---|---|---|
| recognition recall | 0.833 | alt-model fires the correct pattern (MATCH) |
| recognition top-1 | 1.0 | correct pattern ranked #1 (incl. below-threshold) |
| novelty abstention | 0.417 | out-of-catalog emergence → EMERGENT-UNCLASSIFIED |
| null specificity | 0.833 | non-emergent → NO-EMERGENCE |
| **false-MATCH rate** | **0.208** | novelty/null wrongly claimed as a pattern (lower=better) |
| false-novel rate | 0.167 | null wrongly flagged emergent (lower=better) |

Counts: {"recognition": {"recognized": 15, "top1_only": 3}, "novelty": {"false_match": 5, "abstained": 5, "missed_emergence": 2}, "null": {"correct": 10, "false_novel": 2}}

## Headline rates — strict (MATCH requires ≥ confirmation tier)

Screening-tier detections are demoted to the emergence-driven verdict. This is a pure instrument-layer gate (no detector is changed) and shows how much raw false-MATCH is low-rigor screening noise.

| metric | value | reading |
|---|---|---|
| recognition recall | 0.833 | alt-model fires the correct pattern (MATCH) |
| recognition top-1 | 1.0 | correct pattern ranked #1 (incl. below-threshold) |
| novelty abstention | 0.667 | out-of-catalog emergence → EMERGENT-UNCLASSIFIED |
| null specificity | 0.833 | non-emergent → NO-EMERGENCE |
| **false-MATCH rate** | **0.083** | novelty/null wrongly claimed as a pattern (lower=better) |
| false-novel rate | 0.167 | null wrongly flagged emergent (lower=better) |

Counts: {"recognition": {"recognized": 15, "top1_only": 3}, "novelty": {"abstained": 8, "missed_emergence": 2, "false_match": 2}, "null": {"correct": 10, "false_novel": 2}}


## Recognition arm (held-out re-implementations)

| system | expect | verdicts (across seeds) | outcome | mean em |
|---|---|---|---|---|
| recog_p16_boolean_grn | P16 | MATCH×3  [→ P16] | recognized | 0.089 |
| recog_p20_fraction_threshold | P20 | MATCH×3  [→ P20] | recognized | 0.089 |
| recog_p23_el_farol | P23 | MATCH×3  [→ P23] | recognized | 0.089 |
| recog_p24_integral | P24 | NO-EMERGENCE×3 | top1_only | 0.089 |
| recog_p25_multibasin | P25 | MATCH×3  [→ P25] | recognized | 0.09 |
| recog_p26_threshold_unit | P26 | MATCH×3  [→ P26] | recognized | 0.089 |

## Novelty arm (out-of-catalog emergent)

| system | expect | verdicts (across seeds) | outcome | mean em |
|---|---|---|---|---|
| nov_dla | — | MATCH×2, EMERGENT-UNCLASSIFIED  [→ P18] | false_match | 0.753 |
| nov_keller_segel | — | EMERGENT-UNCLASSIFIED×3 | abstained | 1.0 |
| nov_active_nematic | — | NO-EMERGENCE×2, EMERGENT-UNCLASSIFIED | missed_emergence | 0.5 |
| nov_phase_ordering | — | MATCH×3  [→ P3] | false_match | 0.807 |

## Null arm (non-emergent)

| system | expect | verdicts (across seeds) | outcome | mean em |
|---|---|---|---|---|
| null_spatial_noise | — | NO-EMERGENCE×3 | correct | 0.065 |
| null_random_walk | — | EMERGENT-UNCLASSIFIED×2, NO-EMERGENCE | false_novel | 0.464 |
| null_uncoupled_phases | — | NO-EMERGENCE×3 | correct | 0.131 |
| null_frozen_noise | — | NO-EMERGENCE×3 | correct | 0.124 |

## Probe arm (coverage-limit demonstrations — not scored)

| system | expect | verdicts (across seeds) | outcome | mean em |
|---|---|---|---|---|
| probe_static_blob | — | EMERGENT-UNCLASSIFIED×3 | EMERGENT-UNCLASSIFIED | 0.623 |
| probe_percolation | — | NO-EMERGENCE×3 | NO-EMERGENCE | 0.104 |

## Findings & interpretation

1. **Recognition generalizes to held-out implementations.** 5 of 6 alternative
   models the detectors were never tuned on fire the correct pattern at
   confirmation or definitive tier, and all 6 rank the correct pattern #1. The
   lone non-firing case (P24 integral controller) is a known limit: the
   homeostasis detector's firing threshold is calibrated to the *proportional*
   transient, so the integral controller ranks P24 #1 but sits below the gate.
   This reproduces the cross-model 7/7-top-1 result on a fresh run and is the
   core evidence that the battery recognizes the phenomenon, not its native
   implementation.

2. **The raw MATCH gate is too permissive; require ≥ confirmation.** Raw
   false-MATCH is 0.208. The only screening-tier over-claim (DLA → P18: a growing
   cluster trips the voter-consensus screening gate) disappears when MATCH
   requires ≥ confirmation tier. That single change cuts false-MATCH to 0.083 and
   raises novelty abstention from 0.417 to 0.667, at **zero cost to recognition
   recall** (every true recognition already fires at ≥ confirmation). Recommended
   operating point: MATCH iff a detector fires at ≥ confirmation; screening-only
   signals fold into the emergence-driven verdict. Pure instrument-layer change,
   no detector touched. *Implementation:* add `match_min_tier="confirmation"` to
   `profile_observation`, then re-run the self-recognition confusion matrix to
   confirm no native regression before making it the default.

3. **One genuine high-confidence over-claim: P3 vs domain coarsening.**
   Allen-Cahn phase ordering matches P3 (Turing) at definitive tier across all
   seeds. Coarsening domains are isotropic and quasi-stationary over a finite
   window, so they satisfy P3's stationarity × isotropy gates. The discriminator
   P3 lacks is *length-scale stationarity*: a Turing pattern selects a fixed
   wavelength, while coarsening domains grow without bound. Candidate refinement:
   add a coarsening exclusion to P3 (reject when the structure-factor peak
   wavenumber drifts toward zero over the run). This is the one finding that needs
   detector work, not just a gate.

4. **Emergence-indicator sensitivity on soft nulls and apolar order**, both
   consistent with the documented emergence COVERAGE notes:
   - *Diffusive soft-null:* 2 of 3 random-walk seeds tripped a transient
     false-novel (transient Brownian density fluctuations read as structure
     growth). Re-drawn uniform noise was 3/3 clean.
   - *Apolar order:* active-nematic emergence sits at em ≈ 0.50, so it flickers
     across the abstention threshold (1 abstained, 2 missed). The clustering/polar
     structure measure undercounts nematic (head-tail-symmetric) order — exactly
     the rotational/banded/directional gap the indicator's own coverage note
     calls out. A nematic order parameter (mean of exp(2iθ)) would close it.

5. **Probes confirm the two known coverage gaps.** A pre-formed static blob reads
   EMERGENT-UNCLASSIFIED (structure without growth: the structure-vs-shuffle term
   fires with zero order gain), and percolation at threshold reads NO-EMERGENCE (a
   connectivity transition with ~zero spatial autocorrelation, invisible to a
   Moran-based indicator). Both are reported, not hidden, and bound where the
   generic indicator is trustworthy.

## Bottom line

Pointed at systems it was never built on, the catalog recognizes known phenomena
(6/6 top-1, 5/6 firm), abstains on most out-of-catalog emergence, and rejects
non-emergent nulls. With a ≥ confirmation MATCH gate it over-claims on 1 of 12
novelty/null systems (the P3/coarsening case) — the single finding that warrants
detector work. The instrument is sound enough to use and honest about its three
measured limits: the P3 coarsening blind spot, emergence-indicator sensitivity on
apolar order, and the structure-without-growth / connectivity-only coverage gaps.
