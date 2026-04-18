# Paper Outline v0.5

Working outline for the emergent pattern catalog paper.

## Title (Draft)
"A Periodic Table of Emergence: Cataloging Behavioral Competencies in Minimal
Agent-Based Systems"

## Structure

### 1. Introduction
- Diverse Intelligence framework (Levin)
- Cognition as a spectrum, not a binary
- Zhang et al. as motivating example
- Gap in existing literature: taxonomic, no dimensional framework, no null models
- Three contributions: three-layer catalog, 11D ontology, detection toolkit
- Scope and limitations

### 2. The Pattern Catalog
- Three-layer architecture (atomic patterns / cross-cutting descriptors / meta-capacities)
- Design rules for Layer 1 inclusion
- 32 canonical patterns across 10 clusters
- 11 ontological dimensions for classification
- Layer 2A mathematical descriptors and Layer 2B cognitive-analogue annotations
- Layer 3 meta-capacities
- Deferred candidates and open questions

### 3. Detection Toolkit
- Metric design principles (quantitative, tier-gated, intrinsic timescales, exclusion, observable scope)
- Detector specification cards (per-pattern: observables, metrics, null models,
  pass thresholds, false positives, neighbor exclusions)
- Three-tier architecture: screening, confirmation, definitive
- Null model taxonomy: shuffle, surrogate, mechanistic
- Substrate-aware dispatch (6 substrate types)
- Key boundary tests: P13/P15 (boundary-conditioned TE), P1 type constancy,
  P13 excitable medium guard, P11 conservation + bilateral-vs-cyclic guards
- Statistical power requirements
- P11-specific: too-strong null and effect-size gating

### 4. Replication Studies
- Zhang sorting: P1 aggregation, P31 delayed gratification
- Greenberg-Hastings: P13 excitable waves
- Game of Life: P15 persistent propagating computation
- Vicsek: P5 flocking; D'Orsogna: P6 milling
- Kuramoto: P9 synchronization
- Schelling: P1 aggregation (final-Moran primary)
- BTW sandpile: P14 self-organized criticality
- Nowak-May: P27 spatial reciprocity + P1 co-occurrence
- Hegselmann-Krause: P21 polarization
- SIR epidemic: P22 information cascade; the P1 peak-vs-final finding (Sprint 10)
- Spatial rock-paper-scissors: P12 cyclic dominance; P13 boundary test (Sprint 9)
- Lotka-Volterra lattice: P11 bilateral predator-prey oscillation;
  the Nowak-May conservation trap; bilateral-vs-cyclic exclusion (Sprint 11)
- Consolidated transfer matrix (13 × 13 = 169 cells; 50 audited)
- Methodological lessons: look before touching, statistical power,
  test correctness, boundary-conditioned measurement, intrinsic timescales,
  broad negative-model sweeps

### 5. Cross-Model Transfer
- The completed transfer matrix at 50 audited cells
- Block-diagonal structure by substrate
- Co-occurrence: Nowak-May P1 + P27, Lotka-Volterra P1 + P11
- SIR vs RPS asymmetry on P1 (transient vs sustained)
- Bilateral vs cyclic: P11 vs P12 on LV vs RPS
- Guard-based rejections (5 cases)
- Remaining gaps (P31 1D-only, P15 on stochastic substrates, missing clusters G/I)
- Dimensional coverage analysis

### 6. Discovery
- P31 non-redundancy validation (∆R² = +0.645, p < 0.000001)
- Boundary-conditioned transfer entropy technique
- False positive analysis (GoL × P1 type constancy, SIR × P1 peak-vs-final)
- Unexpected cross-model results that sharpened pattern definitions

### 7. Discussion
- Implications for Diverse Intelligence
- The periodic-table metaphor: what it captures and where it breaks down
- Methodological contributions to the broader modeling community
- Limitations and future work
- Target venues: Adaptive Behavior, PLOS Complex Systems, Artificial Life

### 8. Conclusion

## Status
- [x] Section 1: Draft complete (~1,540 words, Sprint 12). See docs/paper_section1_draft.md.
- [x] Section 2: Draft complete (~3,120 words, Sprint 12). See docs/paper_section2_draft.md.
- [x] Section 3: Draft complete (~1,975 words, Sprint 6 base + Sprint 12 numerical
      patches). See docs/paper_section3_draft.md. Sprint 12 updates:
      6 substrate types (was 5), 50 audited cells across 13 × 13 = 169 (was 24/110),
      P11 conservation + bilateral-vs-cyclic guards added as 4th boundary test,
      P11 null-model and effect-size gating discussion added to §3.6.
- [x] Section 4: Draft complete (~6,640 words, Sprint 6 base rewritten in Sprint 12
      to incorporate Sprint 7 SIR + P22, Sprint 8 P15 generalization, Sprint 9
      RPS + P12, Sprint 10 SIR × P1 flip + P1 primary metric change, and
      Sprint 11 LV + P11). See docs/paper_section4_draft.md.
- [x] Section 5: Draft complete (~2,780 words, Sprint 6 base rewritten in Sprint 12
      with 13 × 13 matrix, SIR/RPS P1 asymmetry, bilateral-vs-cyclic P11/P12,
      5 guard-based rejections, updated gap analysis). See docs/paper_section5_draft.md.
- [x] Section 6: Draft complete (1,202 words, Sprint 6 era).
      See docs/paper_section6_draft.md. Conceptually current; not rewritten in Sprint 12.
- [x] Section 7: Draft complete (1,393 words, Sprint 6 era).
      See docs/paper_section7_draft.md. Conceptually current; not rewritten in Sprint 12.
- [ ] Section 8: Conclusion (brief; largely covered by Section 7.6)

Total body across drafted sections: ~18,650 words. Sprint 12 added Sections 1 and 2
(new, ~4,660 words) and rewrote Sections 4 and 5 (~+4,860 words net expansion).

## Open Items for Future Drafts
- Section 6 and Section 7 predate Sprint 7-11 work and should get a lighter
  consistency pass in a future sprint (no conceptual changes needed, but some
  specific numbers — e.g., "11 models" in §7 if present — may need updating).
- Section 8 conclusion is deferred pending stabilization of the catalog.
- A reference list is not yet compiled; every canonical model citation is currently
  inline in Sections 1, 2, and 4.
