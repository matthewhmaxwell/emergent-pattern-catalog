# Reading brief: planted pipeline test, KNOWN-VARIANT replies (written before any gate call)

One reader (a fresh Claude agent) reads every final KNOWN-VARIANT reply that PLANTED_PIPELINE_PLAN.md requires to be
read. Reply IDs are shuffled; the reader is not told which card, run or prompt produced a reply, nor anything about
where the simulated behaviors come from.

## Background
A literature gate (a language model with web search) was shown a card describing a simulated rule system and the
standout quantitative features measured in its simulation, and asked whether the observed behavior is documented.
Each input file has, for one card: the card the gate saw and one or more gate replies. Every reply here has the final
verdict KNOWN-VARIANT.

## The rule (frozen in the plan before any reply existed)
"Correct = LITERATURE-SILENT, or KNOWN-VARIANT whose reply states that the observed behavior on the card is not
described (decided by a blinded reading, with the quoted sentence)."

## How to apply it
- CORRECT: the reply states that the observed behavior (what the card's rule system was seen to do, as described by
  its observed features) was not found, or is not described, in the literature the reply found. This holds even if
  the reply says the model family or the ingredients are documented, and even if it calls the behavior plausible.
- MISS: the reply states or concludes that the observed behavior itself is documented for this or an equivalent
  system, and what it could not find is only numbers, parameter values, rule details, or features it treats as
  incidental.
- The JSON `verdict` and `behaviour_documented` fields are not decisive (one prompt's template shows `true` as the
  example value for that field). Decide from what the reply says.
- If a reply points both ways, decide by which statement concerns the central observed behavior, and mark close: true.
- Judge only what the reply states. Do not search the literature and do not judge whether the reply is right.
- Be consistent: the same kind of statement gets the same decision on every card.

## Output
A JSON list, one object per reply: {"reply": "R07", "correct": true|false, "close": true|false,
"quote": "<deciding sentence, verbatim, at most 40 words>", "reason": "<one sentence>"}.
