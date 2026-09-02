# Skill Benchmark: de-slop

**Model**: claude-fable-5-1 (subagents inherit the session model)
**Date**: 2026-09-02T12:10:13Z
**Evals**: 1, 2, 3, 4, 5 (1 runs each per configuration)

## Summary

| Metric | New Skill | Old Skill | Delta |
|--------|------------|---------------|-------|
| Pass Rate | 99% ± 2% | 93% ± 7% | +0.06 |
| Time | 489.4s ± 150.4s | 326.1s ± 136.3s | +163.4s |
| Tokens | 123934 ± 15205 | 89114 ± 11594 | +34819 |

## Notes

- One run per configuration per eval (n=1), so the deltas are directional, not statistical. Both configurations ran on the same model; the comparison isolates the skill text and the scanner, not model capability.
- The old skill was already strong on this model: 58 of 63 assertions passed without the revision. The revision's gains concentrate where the new instructions add a mechanical step (scanner, compare mode) or a stated judgment (the 'Kept on purpose' block, the protect-list collision report).
- Assertions that pass only with the new skill (4): eval 1: Process: user_notes.md reports running scripts/slop_scan.py on the draft or rewr; eval 1: [judged] The rewrite reads as institutional prose a person wrote: no new formula; eval 3: Process: user_notes.md reports running scripts/slop_scan.py; eval 5: Process: user_notes.md reports running scripts/slop_scan.py on the draft.
- Assertions that pass only with the old skill (0): none.
- The three 'Process: ran scripts/slop_scan.py' assertions can only pass with the new skill (the script did not exist before). They measure instruction-following, not output quality; discount them and the new skill still leads by one output assertion (eval 1, no new formulas introduced) with no regressions.
- Eval 1, vendor name: both configurations kept 'Microsoft Teams'. The new run kept it deliberately, listed it under 'Kept on purpose' with the reason (memo section rather than deck body; the user asked that nothing factual be lost) and handed the call to the user. The assertion is stricter than the profile rule, which names deck body content; treat this as a judgment call rather than a miss.
- Cost: the new skill used more tokens (mean 123934 vs 89114) and more wall time (mean 489s vs 326s), almost all of it in running the scanner before and after the rewrite and reading current-register.md. The eval-2 and eval-5 runs also spent turns working around scanner false positives that have since been fixed (emoji bullets on plain surfaces; salutation and sign-off lines counted as sentences).
- Eval 3 (light touch on a bylined paragraph) is the clearest behavioral difference: the new run changed 1.8% of the words, kept every protected move, and surfaced the ', not by the SOP' flag as a protect-list collision for the author to decide, instead of editing it. The old run also passed the light-touch bar but without the collision report.
- Eval 4 (ingest) passed every assertion in both configurations; the new run additionally deduped against the provisional 2026-09 batch and proposed promoting the 'aphoristic closer paragraph' entry on its first real catch, which is the maintenance behavior the revision specifies.
- Assertions that pass in both configurations (facts intact, no dashes, no contrasts, domain terms kept) do not discriminate between versions on this model; they guard against regressions rather than measure improvement.
