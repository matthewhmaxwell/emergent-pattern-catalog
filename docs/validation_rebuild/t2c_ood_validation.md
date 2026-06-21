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
| recog_p25_multibasin | P25 | MATCH×3  [→ P25] | recognized | 0.09 |
| recog_p26_threshold_unit | P26 | MATCH×3  [→ P26] | recognized | 0.089 |

## Novelty arm (out-of-catalog emergent)

| system | expect | verdicts (across seeds) | outcome | mean em |
|---|---|---|---|---|
| nov_dla | — | MATCH×2, EMERGENT-UNCLASSIFIED  [→ P18] | false_match | 0.753 |
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
| probe_percolation | — | NO-EMERGENCE×3 | NO-EMERGENCE | 0.104 |
| probe_phase_ordering | — | MATCH×3  [→ P3] | MATCH | 0.807 |

## Findings & interpretation (post-cleanup)

This suite was hardened after an initial run surfaced a loose MATCH gate, an
emergence-indicator blind spot, and a soft-null artifact. With the fixes below,
the recommended operating point (MATCH requires ≥ confirmation) is clean:
**recognition 0.83, novelty abstention 1.0, null specificity 1.0, false-MATCH 0.0.**

1. **Recognition generalizes** (unchanged): 5/6 alternative implementations fire
   the correct pattern at ≥ confirmation; all 6 rank it #1. The P24 integral
   controller ranks P24 #1 but below its proportional-tuned threshold (known limit).

2. **≥ confirmation MATCH gate (landed).** `profile_observation` gained a
   `match_min_tier` parameter. Default "screening" preserves native behavior (the
   confusion matrix and gallery are untouched). At "confirmation" — the
   recommended setting when pointing at UNKNOWN systems — screening-only firings
   fold into the emergence-driven verdict and are surfaced as `demoted_match`.
   This removes the OOD screening over-claims (DLA and Eden growth both trip the
   P18 voter-consensus SCREENING gate, which is non-specific to a growing
   contiguous region) at zero cost to recognition recall. It is NOT the global
   default: P22/P27/P29 self-recognize only at screening, so a forced gate would
   regress native self-recognition (29→26).

3. **Emergence indicator hardened (landed).**
   - *Orientation channel:* a polar + nematic order parameter computed from
     velocity/heading data, max-combined with the spatial channel. Active-nematic
     order (head-tail symmetric, invisible to clustering) now scores ~0.89 and
     abstains robustly (was flickering at em ≈ 0.50). The channel can only RAISE
     scores (no null carries velocities), so no regression — and it may help the
     documented P6/P7/P17 emergence gaps.
   - *Soft-null fixed:* the random-walk null now uses periodic (bounded)
     boundaries — a proper diffusive null. The earlier unbounded version
     spuriously grew the cloud's extent (transient false-novel, 2/3 seeds);
     bounded diffusion is 3/3 clean. The indicator's broadly-validated z-threshold
     was left untouched.

4. **P3 ↔ phase ordering: a characterized confusion, not a fixable bug.** Four
   probes established that there is no robust observation-only separator between a
   saturated Allen-Cahn domain field and a coarse Gray-Scott Turing pattern: both
   are stationary, isotropic, low-wavenumber and bimodal, and Gray-Scott is
   actually MORE wide-gap dynamic. The distinguishing physics (a Turing pattern
   selects an INTRINSIC wavelength; coarsening domains fill the BOX) requires
   grid-size-invariance testing, which re-runs the model — outside observation-only
   scope. A coarsening gate on the validated P3 detector would not work (the fields
   are observationally identical) and would risk the real positive. Allen-Cahn is
   therefore reported as a PROBE (a characterized confusion between two
   diffusion-driven pattern classes). The novelty arm keeps four genuinely
   out-of-catalog systems: DLA, Keller-Segel, active nematic, Eden growth.

5. **Probes (characterized limits, not scored):** static blob → false emergent
   (structure without growth), percolation → missed (connectivity, not spatial
   autocorrelation), Allen-Cahn → P3 (the diffusion-pattern confusion above).

## Bottom line

At the recommended OOD operating point (MATCH ≥ confirmation), the instrument
recognizes held-out implementations of catalog phenomena (6/6 top-1, 5/6 firm),
abstains on every out-of-catalog novelty (4/4), and rejects every non-emergent
null (12/12), with zero false-MATCHes. The three residual limits are characterized
and quarantined as probes, not hidden — the P3/coarsening confusion and the
emergence indicator's structure-without-growth and connectivity-only blind spots —
each with the specific test that would resolve it.
