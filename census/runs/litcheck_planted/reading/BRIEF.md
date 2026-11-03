# Reading brief: planted positive controls, KNOWN-VARIANT replies (2026-11-03)

One reader (a fresh Claude agent) reads every KNOWN-VARIANT reply that the scoring rule requires to be read. Reply IDs
are shuffled; the reader is not told which arm, run or prompt produced a reply.

## Background
A literature gate (a language model with web search) was shown a card describing a simulated behavior and asked
whether it is documented. The behaviors come from papers posted in September and October 2026 whose authors present
them as new. Each input file has, for one card: the behavior the authors report as new, the card the gate saw, and
several gate replies. Every reply here has the verdict KNOWN-VARIANT.

## The rule (frozen before any reply existed)
"Correct: LITERATURE-SILENT, or KNOWN-VARIANT whose reply states that the reported behavior itself (not only rule
details) is not described. Miss: KNOWN, or KNOWN-VARIANT that says the behavior is described. Scoring is done by
reading each reply against these rules; every decision is recorded with the quoted sentence that decided it."

## How to apply it
- CORRECT: the reply states that the reported behavior (its central emergent feature, as given under "Reported
  behavior" and on the card) was not found, or is not described, in the literature the reply found. This holds even if
  the reply says the model family, the ingredients or the mechanism are documented, and even if it adds that the
  behavior is plausible or expected.
- MISS: the reply states or concludes that the reported behavior itself is documented, for this or an equivalent
  system, and what it says it could not find is only numbers, parameter values, rule details, or features that are not
  central to the reported behavior.
- The JSON `verdict` and `behaviour_documented` fields are not decisive (one prompt's template shows `true` as the
  example value for that field). Decide from what the reply says about the reported behavior.
- If a reply points both ways (for example "the behaviour is documented" and "the specific final-state structure was
  not found"), decide by which statement concerns the central feature of the reported behavior, and mark close: true.
- Judge only what the reply states. Do not search the literature and do not judge whether the reply is right.
- Be consistent: the same kind of statement gets the same decision on every card.

## Output
A JSON list, one object per reply: {"reply": "R07", "card": "planted:N01", "correct": true|false, "close": true|false,
"quote": "<deciding sentence, verbatim, at most 40 words>", "reason": "<one sentence>"}.
