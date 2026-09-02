---
name: de-slop
description: >
  Detect, rewrite, and memorialize AI writing patterns (slop) in prose. Use it
  whenever text will be read by a human: client deliverables, decks, memos,
  emails, LinkedIn posts, essays, research papers, documentation, and Claude's
  own drafts before they go out. Trigger on "de-slop", "humanize", "does this
  read like AI", "sounds like ChatGPT", "clean this up", "tighten this", "edit
  my draft", "polish", "proofread", "remove AI tells", or any request to make
  writing sound like a person wrote it. Also run it silently before delivering
  any client-facing, published, or bylined prose, even when nobody asks. When
  the user pastes text marked "slop:" or says "ingest this", run the ingest
  flow so the pattern becomes a permanent rule. Ships a mechanical scanner
  (scripts/slop_scan.py) for dashes, colon rhythm, cadence, tricolons, binary
  contrasts, tiered vocabulary, fact drift, and the 2026 model register.
license: MIT
metadata:
  version: "1.1.0"
  revised: "2026-09-02"
  upstream: kjmagnan1s/anti-slop v0.1.0
  merged: [ehmo/slopbeth, adenaufal/anti-slop-writing, jalaalrd/anti-ai-slop-writing]
  pairs-with: matt-voice
---

# de-slop

One maintained skill for removing AI writing patterns and for capturing new
ones as they appear. The rule lists are commodity. The asset is the living
corpus (`references/living-corpus.md`): dated tells caught in real work, tagged
with the mechanism that produces them, including the patterns Matt has flagged
more than once.

The 2026-09 revision changed two things. A scanner (`scripts/slop_scan.py`)
now does the counting that a model does badly by eye, before and after every
rewrite. And `references/current-register.md` records what the current model
generation does once the older tells are trained out: the vocabulary it moved
into, the punctuation the dash rhythm displaced into, and the ways a model
editing its own draft goes wrong.

## Before you start

Decide three things and say them to yourself in one line.

- **Mode.** rewrite (default), detect (flag only), or ingest (capture a new
  tell). Details under Modes.
- **Profile.** The strictness matrix is in `references/patterns.md`.
  Auto-detect from content when none is stated. Anything client-facing and in
  doubt is client-deliverable, the strictest profile. Research papers, emails,
  and plain-text messages have their own profiles in the same file.
- **Whose draft.** A human's draft gets a light hand; strong text gets a light
  edit or none. Your own draft gets the self-edit section below, because your
  tells read as natural to you.

## The pass

Run these in order on every rewrite. Skipping the scan or the fact lock is how
a clean-looking rewrite ships a wrong number.

1. **Scan.** Write the draft to a file and run the scanner from this skill's
   base directory (the path announced when the skill loads):

   ```
   python3 scripts/slop_scan.py draft.md --profile client-deliverable
   python3 scripts/slop_scan.py draft.txt --surface plain      # email, chat
   ```

   It reports P0, P1, and P2 findings with line numbers, plus the rhythm
   numbers: sentence-length spread, the longest run in the 17-23 word band,
   opener repetition, colons per 100 sentences, paragraph shape. Every flag is
   a candidate for judgment. A tricolon can be a real list of three; a colon
   can be the right mark. When Python is unavailable, do the two 30-second
   tests by hand (first word of each sentence; three consecutive sentences of
   17-23 words) and count the dashes and colons yourself.

2. **Lock the facts.** Inventory what the text claims before touching style:
   entities, numbers, dates, units, qualifiers, causal claims, scope limits,
   quotations, the author's stance. These survive intact. If the source is
   vague and no supporting material exists, enter evidence-bound mode: tighten
   and de-hype, never invent the missing specifics. Bracket what the author
   must supply (`[metric]`, `[owner]`). The contract is in
   `references/patterns.md` under "Fact preservation".

3. **Protect list on a byline.** Load `references/protect-list.md` when the
   text is Matt's. Never strip a protected signature at its stated dose. The
   scanner will flag some protected moves (a rhetorical question, a "Y, not X"
   verdict, a deliberate fragment); check each flag against the list before
   acting. On a collision between a floor rule and a protected item, surface
   it and do not auto-edit.

4. **Structure pass.** Structure is the first detection signal, above
   vocabulary. Cadence uniformity, opener repetition, the 17-23 band, the
   fragment/long seesaw, over-fragmented paragraphs, dash rhythm displaced into
   colons, the binary-contrast family, templated lists, the aphoristic closer.
   See "2026 structural tells" in `references/patterns.md` and sections 3 and
   4 of `references/current-register.md`.

5. **Vocabulary pass.** The vocabulary pass is tiered. Tier 1 always replace, Tier 2
   flag in clusters, Tier 3 flag by density. Then the current register
   (`references/current-register.md`, section 1): the words the model moved
   into after delve and robust were trained out. Each is fine alone. Two in a
   paragraph is a tell. Replace with the plain word or, better, with the
   specific thing the word was standing in for.

6. **Seam pass.** At each paragraph boundary generate two or three candidate
   openers, always including "no transition, start cold". Choose by fit. Then
   a monotony read: if two boundaries lean on the same move, break one. Do not
   add transitions to make the text flow; that is the seam problem in reverse.

7. **Adversarial read.** Read the rewrite as the recipient, not the author,
   and answer one question: what would make an editor say a model wrote this?
   Name the three biggest remaining tells. Fix those and stop. An editor who
   keeps going past the third produces the over-polished uniform profile that
   reads as AI.

8. **Post-check.** Run the scanner on the rewrite, and the compare mode
   against the original:

   ```
   python3 scripts/slop_scan.py --compare original.md rewrite.md --profile client-deliverable
   ```

   Compare reports missing or changed numbers, invented numbers, dropped
   qualifiers, missing entities, edit size, register migration (idioms the
   rewrite introduced), flattened rhythm, and a new closer. Any P0 in compare
   blocks delivery. Any P1 gets fixed or gets a stated reason in the change
   summary.

9. **Quality gate.** The rewrite obeys the same rules as the original. No new
   formulas: clipped aphorisms, tidy triads, forced contrast, consultant voice.
   Noun-swap test on every sentence that feels smooth: if it survives
   transplant to another topic with only noun swaps, it needs a fact, scene,
   number, decision, or consequence. If the original was already strong, say
   so and make a light edit or none. Over-editing is a failure.

## Modes

**rewrite** (default). Run the pass. Return the clean text, then a change
summary in this shape:

```
Changed (n)
- what changed, in plain words, and the rule or corpus entry behind it
Kept on purpose
- anything the scanner flagged that stays, with the reason (protected move,
  real list of three, the right colon)
Facts
- n locked, all present   (or: what drifted and how it was restored)
```

Keep the summary plain. A summary that reads "tightened cadence, removed
tells, preserved facts" is a tricolon and has slop in it. If nothing needed
changing, say so in one line and return the original.

**detect**. Flag only, no rewriting. For someone else's text, published text,
or a quick scan. Group by severity with exact spans quoted; a critique that
cannot point to a span is too vague to act on.

```
P0 credibility   "exact span" (line n): label. Why it reads as AI.
P1 obvious       ...
P2 polish        ...
Scan             dashes n, contrasts n, colons n per 100 sentences,
                 longest 17-23 run n, opener paragraphs n
Verdict          one line: post as is / fix P0 and P1 first / leave it alone
```

Apply the profile before flagging. End-of-line emoji and hashtags on a
LinkedIn post are not tells; a bolded lead-in on a deck slide is convention.
Flagging convention as slop costs credibility with the person who asked.

**ingest**. The capture flow in `references/ingestion.md`. Triggered by text
marked `slop:`, by "ingest this", or by catching your own slop after Matt
flags it. Run the scanner on the pasted span first; it names the mechanical
pattern. Produce a dated corpus entry with mechanism tag, tier, before/after,
and the generalized rule, filed into `references/living-corpus.md` after a
protect-list cross-check. A deliberate flag gets shown and filed on a one-word
ok; self-caught slop gets auto-filed with a revertible diff. This is how a
one-off correction becomes a permanent rule.

## Editing your own draft

Most invocations are this case: the skill runs silently before a deliverable,
an email, or a post goes out, on text you just wrote. It is the hardest case,
for one reason. Your tells do not look like tells to you. They look like
fluent prose, because they are what your fluency produces. That has three
consequences.

- **Trust the counts over your ear.** You will read a paragraph with three
  20-word sentences and hear variety. The scanner will not. When the scan and
  your ear disagree on anything countable, the scan wins.
- **Every rule you obey moves the impulse somewhere.** Ban the em dash and
  the aside becomes a colon, then a label sentence ("The catch: ..."), then a
  "which"-tail. Ban the tricolon and it becomes a pair ("clear and
  consistent"), then a two-word abstract balance. Ban "let me know if" and it
  becomes "Happy to draft that if useful." The displacement table in
  `references/current-register.md` (section 4) lists where each rule pushes
  the next draft. Read it before the post-check, then look for those
  specifically.
- **The rewrite can be worse than the draft.** Four ways, in order of
  frequency: the rewrite swaps banned words for the current register
  ("load-bearing", "non-trivial", "the delta"), so it sounds like this year's
  model instead of last year's; the rewrite flattens into short declaratives
  with no subordination, which reads as a telegram; the rewrite rounds a
  number or drops a qualifier while smoothing a sentence; the rewrite gains a
  one-line closer the draft never had. The compare mode checks all four.
  Section 5 of the current-register file has the full list.

When the draft was already fine, the right output is the draft. Say so.

## Quick checks the scanner cannot do

- Is any sentence transplantable to another topic with noun swaps only? Add
  the specific.
- Is an inanimate thing doing a human's verb ("the data tells us", "the
  framework ensures")? Name the human.
- Does each tricolon have three real members? Drop test: if one could go
  unnoticed, cut to the true item or find the real count.
- Is a colon the right mark here, or the dash rhythm in a new coat? Most
  become periods. Not all.
- Is a defined term being synonym-cycled for variety? Repeat the defined term.
- Is a protected signature at its dose, or past it into a tic?
- Did any qualifier get stronger or weaker ("may" to "will", "some" to
  "many")? Revert it.
- Did a specific detail appear that no source supports? Remove or bracket it.
- Is this a plain-text surface (email, message) carrying markdown? Convert to
  prose.
- Would "leave this alone" be the honest answer?

## Scoring

Rate 1-10 each and revise below 35/50: Directness (statements, not
announcements), Rhythm (varied, not metronomic), Trust (respects the reader),
Authenticity (a person wrote it), Density (nothing cuttable). Fidelity is a
pass-or-fail gate. Any claim that differs from the source blocks delivery
regardless of total.

## Maintenance

The update trigger is "Matt spotted a tell", not a schedule. When that
happens, run ingest, grow the corpus, and re-save the skill so the installed
copy carries the entry. Tells age. Re-tier or retire entries as models change,
and date everything so the rise and fall of a tell stays visible. Entries
sourced from a model's self-report (the 2026-09 batch) are provisional until
one is caught in Matt's real work; promote or retire them on first contact.

## References

- `scripts/slop_scan.py`: the scanner. `--json` for machine-readable output,
  `--compare` for original versus rewrite, `--profile` and `--surface` to set
  strictness, `--info` to list single Tier 2 and register hits.
- `references/patterns.md`: the full rule library: tiers, content patterns,
  structure patterns, 2026 structural tells, fact preservation, rewrite
  guards, context profiles (including client-deliverable, research-paper, and
  email), severity tiers. Table of contents at the top.
- `references/current-register.md`: dated 2026-09. The post-delve vocabulary,
  sincerity intensifiers, structural tells that survived the earlier lists,
  the displacement table, self-edit failure modes.
- `references/living-corpus.md`: dated tells with mechanism tags, seeded with
  Matt's repeatedly flagged patterns plus the 2026-09 self-report batch.
- `references/protect-list.md`: what must never be stripped from Matt's
  writing, plus stricter-than-floor settings.
- `references/ingestion.md`: the capture flow for memorializing new slop, with
  the mechanism tags.

## Credits

Fork of anti-slop by kjmagnan1s (MIT), itself a consolidation of
avoid-ai-writing (Conor Bronsdon, MIT), humanizer (MIT, from Wikipedia's
"Signs of AI writing", CC BY-SA 4.0), and stop-slop (Hardik Pandya, MIT).
Merged material from slopbeth (ehmo), anti-slop-writing (adenaufal), and
anti-ai-slop-writing (jalaalrd). The 2026-09 revision adds the scanner and the
current-register file. See CREDITS.md.
