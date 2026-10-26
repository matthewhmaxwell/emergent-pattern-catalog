# Literature-gate controls: scoring rules (fixed 2026-10-26, before any control card was run)

Post-registration, exploratory analysis added after peer review. The frozen census verdicts are not changed.

## Runs
1. **Primed (frozen protocol):** all 40 cards through `census.litcheck_runner.run_one` with the frozen instruction
   text and the census first-pass model (`claude-sonnet-5-5`, web search, at most 10 searches).
2. **Second pass, as in the census:** every first-pass KNOWN-VARIANT verdict goes to the frozen strict second-pass
   prompt (`claude-opus-5-5`).
3. **Neutral:** the 20 fabricated cards again, with one sentence removed from the instructions: "Most such behaviours
   are known; "known" is the expected answer and a perfectly good one." Nothing else changes.

Cards carry no hint that they are controls; `truth.json` is never shown to a model.

## Scoring
- **Known control, identified:** the verdict is KNOWN or KNOWN-VARIANT, and the standard name or citations name the
  intended phenomenon or its defining model. **Partly identified:** the right model family but a different behaviour.
  **Missed:** LITERATURE-SILENT, or a wrong model.
- **Fabricated control, correct:** the verdict is LITERATURE-SILENT, or KNOWN-VARIANT where the reply says that the
  observed behaviour differs from or contradicts the documented behaviour. **False known:** KNOWN (the reply asserts
  the reported behaviour is documented).
- One adjudicator (Claude) scores each verdict against these rules, quoting the reply. Every verdict file is kept.

## What is reported
- Known recall (identified / 20), with an exact 95% binomial interval.
- Fabricated false-known rate (KNOWN / 20), with an exact 95% interval, primed versus neutral.
- The effect of the second pass on both sets.

## Addendum (2026-10-26, written after seeing the first 31 primed verdicts, before this arm was run)
The primed first pass labelled every fabricated card seen so far KNOWN, naming the system's documented behaviour.
In the census protocol only KNOWN-VARIANT verdicts reach the strict, behaviour-focused second pass, so these cards
never would. A diagnostic arm, **strictall**, therefore runs the frozen strict second-pass prompt (claude-opus-5-5) on
all 40 primed verdicts. It is scored with the same rules: known controls should stay KNOWN or KNOWN-VARIANT with the
right phenomenon, and fabricated controls should not be KNOWN. This arm was added after seeing results and is
reported as such.
