# Cross-Model Generalization (T2c) — recognizes the phenomenon, not the implementation

**Date:** 2026-06-16. **HEAD:** 96d63f7+. For each pattern with an independent
SECOND implementation of the same phenomenon, a positive was built from the
ALTERNATIVE model (one the detector was NOT tuned on) and run through the full
battery (`analysis/battery_profile`). A detector that recognises the phenomenon
(not its native implementation) should still rank its own pattern #1.

**Result: 7/7 rank the correct pattern #1; 6/7 are firm MATCHes (detector fired).**

| pattern | tuned on → tested on | top-1 | tier | verdict |
|---|---|---|---|---|
| P16 associative memory | Hopfield → **BooleanGRN** | P16 | confirmation | MATCH |
| P26 stochastic resonance | DoubleWell → **ThresholdUnit** | P26 | definitive | MATCH |
| P29 trail networks | Ant → **Physarum** | P29 | definitive | MATCH |
| P20 quorum sensing | Autoinducer → **FractionThreshold** | P20 | definitive | MATCH |
| P25 equifinality | Canalized landscape → **MultiBasinGRN** | P25 | definitive | MATCH |
| P23 anti-coordination | Minority Game → **El Farol** | P23 | confirmation | MATCH |
| P24 homeostasis | Proportional → **Integral** controller | P24 | none | (ranked top, below threshold) |

**Honest note (P24):** the integral controller regulates with a different
transient than the proportional controller the detector learned; P24 is still the
top-ranked pattern on it, but the detector stays below its firing threshold. A
genuine generalization limit, surfaced rather than hidden — candidate for a
controller-agnostic restoring-feedback signature.

Together with the self-recognition confusion matrix (29/31 firm, ~nil
cross-fire), this is the core evidence that the EPC battery is a measurement
instrument, not a lookup table over its native models.
