# Current model register

Dated 2026-09-02. Re-read the date before trusting a section.

What the current model generation writes like once the older tells are gone.
Provenance: self-report by Claude (Fable 5.1) during the 2026-09 revision of
this skill, checked against that model's own drafts in the same session. None
of these had been flagged by Matt at the time of writing. Each is provisional
until one is caught in real work; then it is promoted into the corpus at a
confirmed tier, or retired. The vocabulary in section 1 will age fastest. The
displacement table and the self-edit failures will age slowest, because they
describe the mechanism rather than the surface.

Why this file exists. The 2024-25 tells (delve, tapestry, robust, the em dash)
were trained out of newer models by exactly the kind of list this skill
carries. The impulse behind them was not trained out. It refilled from a
different corpus: technical Twitter, Hacker News, consulting decks, product
management writing. The result is prose that passes every old checklist and
still reads, to a careful reader, as a model. Section 1 catalogs that
register. Sections 3 and 4 do the same for structure.

Sections: 1 vocabulary, 2 sincerity intensifiers, 3 structural tells,
4 displacement table, 5 self-edit failure modes, 6 scanner coverage,
7 retirement.

## 1. The post-delve vocabulary (Tier 2: flag in clusters)

Each of these is a fine word in a human's mouth. The tell is density: two or
more distinct items from this list in one paragraph, or the same one twice in
a short piece. The fix is not a synonym. It is the plain word, or the specific
thing the word was standing in for.

**Praise for prose, design, or thinking.** crisp, clean, tight, sharp,
opinionated, principled, thoughtful, intentional, deliberate, elegant,
quietly (as in "quietly one of the best"), legible.
Plain: say what is good about it. "The memo is crisp" says nothing; "the memo
makes one ask in the first line" does.

**Analysis-flavored weight words.** load-bearing, non-trivial, orthogonal,
first-order, second-order, at the margin, the delta, high-signal, low-signal,
high-leverage, meaningfully, materially, substantively, directionally,
durable.
Plain: important, separate, biggest, the difference, a lot, some. Better: the
number.

**Verbs of handling.** surface (a risk), flag, pressure-test, stress-test,
sanity-check, gut-check, tease apart, disentangle, reach for, lean on, lean
into, land (a point lands), bake in, map onto, reframe, zoom out, step back,
paper over, unpack.
Plain: raise, check, separate, use, prefer, include, hide, explain.

**Nouns from the operator toolkit.** the ask, the unlock, the lift, the
through-line, the north star, table stakes, the long pole, a forcing function,
a footgun, the happy path, the failure mode, the mental model, the lens, the
framing, the shape of, the contours of, the texture of, surface area, the
flywheel, scaffolding, primitives, the plumbing, hygiene, the bottleneck (when
not literal).
Plain: the request, what it makes possible, the effort, the theme, the goal,
the minimum, the slowest part, a deadline, a trap, the normal case, how it
breaks, how I think about it.

**Idioms of weight and balance.** does real work, earns its place, pulls its
weight, carries weight, holds the tension, sits in tension with, there's a
tension between, cuts both ways, double-edged, there's a version of this
where, if anything, it turns out, there's a reason X, it's no accident that,
this is where X comes in, some of this is X and some of it is Y, part of the
answer is.
Plain: matters, is needed, conflicts with. Most of these are stalling before
a claim. Delete the idiom and the claim is usually already there.

**Speaker frames.** I'd argue, I'd push back on, I'd gently push back, I'd
note, I'd add, I'd offer, fair enough, that's fair, reasonable people can
disagree, the honest version, the short version, the long version, the boring
version, the uncomfortable version.
Plain: state the claim. "I'd argue X" is X with a hedge glued on. "The short
version is X" is X.

Collisions. Some of these are terms of art on some surfaces. "failure mode",
"happy path", "primitive", and "orthogonal" are exact in engineering writing.
"operating posture" is a matt-voice preference and is not listed here. The
scanner exempts the engineering terms on the technical-blog and docs
profiles. On a byline, the protect list wins.

## 2. Sincerity intensifiers (Tier 3: flag by density)

actual, actually, real, really, genuine, genuinely, truly, literally, and the
pairs specific / specifically and concrete / concretely when they praise the
writer's own claim instead of being specific.

Mechanism: reward-tuning. Training rewards prose that sounds sincere and
grounded, and the cheapest way to sound grounded is "actual" in front of a
noun. "Name the actual risk" and "name the risk" make the same request; the
first implies a fake risk was on the table. A cousin of performed honesty:
marking sincerity implies the alternative.

Rule: above roughly six per thousand words, cut them all and re-read. If a
sentence feels weaker without its "actually", the sentence needed evidence,
not an adverb.

## 3. Structural tells that survived the 2026-06 list

The 2026-06 list in patterns.md (cadence uniformity, the bimodal seesaw, colon
displacement, over-fragmented paragraphs) still holds. These are additional
shapes, each with the sentence that shows it.

- **Label-colon sentence.** "The catch: nobody owns the intake queue."
  "Translation: the committee will not meet." "The result: three workflows."
  The dash aside, having lost the dash and then the colon-as-aside, reappears
  as a label. Tier 1. Fix: a subject and a verb. "Nobody owns the intake
  queue."
- **Which-tail verdict.** "..., which is exactly the problem." "..., which is
  the point." "..., which matters." A verdict smuggled in as a relative clause
  so it does not have to be argued. Tier 1. Fix: make it a sentence with a
  subject, or cut it; the reader reached the verdict already.
- **One-word opener.** "Yes." "No." "Correct." "Short answer: no."
  Directness training produced a directness tic. Tier 2. Fix: begin with the
  first substantive sentence; the yes or no is inside it.
- **Enumerator opener.** "Two things." "A few notes." "Three observations."
  Announces a count instead of starting. Tier 2. Fix: start with the first
  thing. If the count matters, it will be visible.
- **Imperative staging.** "Consider the intake form." "Take activation."
  "Look at the Q2 numbers." A stage direction before the point. Tier 2. Fix:
  the point.
- **Restatement chain.** The claim, then "In other words," then "Put
  simply,". Helpfulness training rewards saying it again more simply. Tier 1
  when two or more appear. The second restatement is always the cut, and
  usually the first version was the clearer one.
- **Aphoristic closer paragraph.** A final one-sentence paragraph delivering
  the verdict: "That's the whole game." "That is the job." "Nothing more to
  it." The banned "In conclusion" came back as an epigram. Tier 1 on
  deliverables, Tier 2 elsewhere. Fix: end on the last real point, even if it
  lands mid-thought.
- **Parenthetical wink.** "(yes, really)", "(this is the part people skip)",
  "(it isn't)", "(spoiler: no)". Chat register leaking into prose. Tier 1.
  Fix: cut the parenthesis; if it carried a claim, make the claim.
- **Tail-qualifier stack.** "at least for now", "for what it's worth", "as
  far as I can tell", "in my experience", hung off the ends of sentences,
  often two in a paragraph. Uncertainty conditioning that landed at the tail
  instead of in the claim. Tier 2. Fix: put the hedge inside the claim once
  ("Two of the three sites have not reported") or commit.
- **"The [adjective] version."** "The honest version is...", "the short
  version", "the uncomfortable version". A restatement announcer wearing
  candor. Tier 2; Tier 1 when the adjective is honest, real, or candid, since
  that is performed honesty.
- **Question as heading.** "Why does this matter?" "What changed?" as
  section headings, answered in the section. Instruction-tuned document
  scaffolding. Tier 2 on deliverables. Fix: a heading that states the answer.
- **Offer closer.** "Happy to draft the summary if useful." "If helpful, I
  can pull the numbers." "Want me to expand on any of these?" The banned "let
  me know if you need anything" came back as an offer. Tier 1 in deliverables
  and email. Fix: cut. If an action is needed, ask for the decision.
- **Asker-evaluation opener.** "Good instinct." "Right call." "Fair push."
  "Worth asking." The banned "great question" came back as a compliment on
  the asker's judgment. Tier 1. Fix: cut; answer.
- **"Because" fragment answer.** A question, then "Because that's where the
  risk is." as its own sentence. Tier 2. Fix: join it to the claim it
  answers, or let the claim stand alone.
- **Some / some parallel.** "Some of this is process. Some of it is trust." A
  tricolon with one member removed. Tier 2. Fix: say which part is which and
  how much.
- **Sentence-initial And, But, So as a rhythm crutch.** Fine at a low rate.
  When a third of the sentences in a section begin with one of them, the
  writer is faking spoken cadence. Tier 3 by density.

## 4. Displacement table: where the impulse goes when a rule bites

Naming the mechanism is what lets a tell be predicted before it is common.
Each rule this skill enforces removes a surface form and leaves the impulse in
place. This is where each one goes next, and what catches it.

| Rule obeyed | Where the same impulse reappears | Catch it by |
|---|---|---|
| No em or en dashes | colons; comma appositives; parentheses; "which" tails; one-sentence paragraphs; sentence-initial And, But, So | counting interruptions per sentence, not glyphs |
| Fewer colons | label-colon sentences ("The catch: ..."); "Meaning," and "Translation," openers; period-fragments ("The fix. Stop."); semicolons | label-colon and fragment checks; a period that splits one thought |
| No binary contrast | "X, but Y" concessions; "X. Y is also true."; "less about X, more about Y"; "Y, not X."; "rather than X, Y"; "X matters less than Y" | delete the foil: if nothing is lost, it was a contrast |
| No tricolons | pairs ("clear and consistent"); two-word abstract balances ("real but bounded"); four-item lists with one filler; "X, Y, Z, and more" | count the true members; a pair of abstractions is still a pattern |
| No delve, leverage, robust, ensure | dig into, lean on, solid, make sure; then the whole post-delve register (section 1) | register cluster check; the register is the tell, not the word |
| No hedging stacks | false precision ("roughly 40%" with no source); confident vagueness ("materially", "meaningfully"); tail qualifiers | evidence-bound mode: every number has a source or a bracket |
| No bullet lists | paragraphs that are bullets with periods: one sentence each, same template, no connectives | subordinators per sentence; templated-paragraph check |
| No bold | colon-led inline labels ("Timing: ..."); the lead sentence turned into a heading | same test as inline-header lists |
| No "In conclusion" | the aphoristic closer paragraph; "That's the whole game." | last-paragraph length and novelty |
| No "Great question" | "Good instinct.", "Right call.", "Fair push.", opening with "Yes." | any opener that evaluates the asker or the ask |
| No "Let me know if" | "Happy to X if useful.", "If helpful, I can...", "Want me to...?" | any closer that offers |
| No passive voice | "you" as a universal actor when a specific person acted; "we" where it was "I" | name the actor |
| No "It is important to note" | "Worth noting", "One thing to flag", "A caveat:", "Note that" | announcer list |
| No uniform sentence length | the bimodal seesaw: fragment, long sentence, fragment, long sentence | alternation count |
| No long paragraphs | one-sentence paragraphs everywhere | median sentences per paragraph |
| No chatbot artifacts | the parenthetical wink; "(yes, really)" | parenthesis scan |

The scanner implements the right-hand column where it is countable. The rest
is the adversarial read.

## 5. Self-edit failure modes

When the model edits its own draft under these rules, the rewrite can be
worse than the draft. In rough order of frequency:

1. **Register migration.** Flagged words are swapped for section 1 idioms.
   The rewrite now sounds like this year's model. Compare mode reports the
   idioms the rewrite introduced.
2. **Flattened rhythm.** Every sentence short, no subordination, telegraphic.
   Compare mode flags a mean sentence length that fell by a third or more, or
   subordinators per sentence that halved. Fix: reconnect related thoughts
   with because, which, although, when.
3. **Fact drift while smoothing.** A number rounded, a qualifier dropped,
   "may" hardened to "will", "some" inflated to "many". Compare mode reports
   missing or changed numbers and dropped qualifiers. Revert each.
4. **A new closer.** The rewrite ends on a one-line verdict the draft never
   had. Compare mode flags it. Cut it.
5. **The summary block slops.** "Tightened cadence, removed tells, preserved
   facts." The change summary is prose too. Write it as specific edits with
   the rule behind each.
6. **Reorganization into the arc.** The model restructures the author's
   paragraphs into claim, evidence, implication, every paragraph. Keep the
   author's order and paragraph count unless asked.
7. **Seam transitions added for flow.** "That said," "In practice," "To that
   end," inserted at paragraph boundaries. The seam pass generates candidates
   that include "start cold"; usually that one wins.
8. **Jargon stripping.** A defined term or a term of art replaced with a
   plain-language synonym, which changes the meaning. The protect list covers
   Matt's domains; the general rule is that a term defined once stays that
   term.
9. **Zero-tolerance colons.** The colon rule is about rhythm density, not the
   mark. A rewrite that replaces every colon with a period produces fragment
   chains, which is failure mode 2.
10. **Over-editing strong text.** Changing more than about a fifth of a good
    draft without a stated cause. "Leave it alone" is a valid output and
    sometimes the correct one.

## 6. What the scanner sees and what it cannot

It sees: dashes, colon and semicolon rate, sentence-length spread and runs,
opener repetition, paragraph shape, fragment triads, tricolon candidates, the
binary-contrast shapes, tiered vocabulary and the register in clusters,
announcers, honesty markers, padding, chatbot and tool artifacts, markdown on
plain surfaces, heading case, and in compare mode the drift, invention, edit
size, register migration, flattening, and new-closer checks.

It cannot see: whether a tricolon is a real list of three, whether a colon is
the right mark, whether a claim is true, whether a protected move is at dose,
whether a sentence says anything (the noun-swap test), or whether the text was
already fine. Those are the pass.

## 7. Retirement

Re-read this file when the model generation changes. For each section 1
group, ask whether the words still show up at density in fresh model output.
Move confirmed items into the corpus with a mechanism tag. Delete the ones
that faded. Date the edit at the top.
