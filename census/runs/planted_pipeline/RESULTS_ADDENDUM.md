# Planted pipeline: addendum to RESULTS.md (2026-11-06)

RESULTS.md is written by the frozen `planted_pipeline.report`. This addendum adds what the plan lists as descriptive,
the verification-path losses (verification_path.json), and notes needed to read the table.

## Verification path, where each behavior is lost
Screen 3 (N08, N13, N14); namer 1 (N01); verification 7 of the 9 with cards: compute cap 3 (N03, N04, N11), decomposes 1
(N02), artifact 1 (N17), single mechanism without a flagged verdict, so no Stage 3, 2 (N12, N15); Stage 3 (b) 2 (N07,
N19: all four cards CLOSELY RELATED). Undocumented candidates: 0 of 20.

## Stage 3
Stage 3 (a) (frozen strict prompt with the simulation evidence, opus, target blocked): N07 KNOWN-VARIANT, N19 F0-F2
KNOWN. Stage 3 (b) (one Claude research agent per card, census brief plus "count only sources dated before" the target's
v1 date, not told why): CLOSELY RELATED for all four, overturning the three KNOWN verdicts. Two of the four agents found
the target paper itself (N07; N19 F0), dated on the cutoff day, and excluded it under the date rule; the N07 agent wrote
that the label would otherwise likely be DOCUMENTED. The agents found misattributed or mischaracterized citations in the
Stage 3 replies (N19: Painter 2009 described as a review, a 2025 PRL cited without authors and for a claim it does not
support, the nutrient-as-repellent equivalence presented as documented). These are the labels the census's own Stage 3
phenomena mostly received (61 of 84 CLOSELY RELATED).

## Descriptive (plan "Scoring")
- Presence at census scale (presence.json): present 4 (N07, N11, N13, N15), absent 4 (N03, N08, N14, N19), not
  observable in one run at one setting 4 (N01, N02, N12, N17), undetermined 1 (N04).
- Card features against the paper's measure (feature_check/output.json, one Claude agent, blind to verdicts): YES 0,
  PARTLY 11, NO 1 (N04); four papers' measures need a sweep, a size series or a comparison system.
- Hand-written planted cards for the same nine models (S11, primary arm, by card majority under the plan's reading
  rule): correct for N04 and N17; the pipeline cards: 0 of 9.
- Novelty dispute recomputation (N14, N17, N19 disputed): 0 of 10 in scope on every path.
- Fidelity: RESULTS.md shows the module's `pass` field, which for N11 is the plan's original stripes check (False);
  under deviation P2 the deciding check is agreement with the authors' anisotropy, which passed (FIDELITY.md).

## Notes
- N19 at 4x area: f_moran and f_peak_sharp fail persistence because the knock-out field, uniform to about 1e-5,
  crosses the fingerprint's zero-variance cutoff at 4x area (Moran's I of the residual about 1.0, against 0 at census
  size); the real runs have f_moran 0.98 in both. Frozen classes unchanged; f_std and f_change are the robust features.
- Compute: the 38 verification jobs ran 2026-11-04 11:38 to 2026-11-05 13:56 UTC under P6; N04 census and N11 3x length
  were killed at 24 CPU-h as projected.
