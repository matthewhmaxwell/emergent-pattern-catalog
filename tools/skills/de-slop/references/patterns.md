# Patterns: the rule library (the floor)

The deduped union of avoid-ai-writing, humanizer, and stop-slop. Where all three
overlapped (~70%), there is one entry. Examples are rewritten in our own words
rather than copied (see CREDITS.md on the Wikipedia lineage). Patterns that
collide with a writer's voice are cross-referenced to `protect-list.md`; on a
byline with a protect list, the protect list wins.

Order of operations: structure first (the #1 detection signal), then vocabulary,
then formatting. Fixing words while leaving robotic rhythm still reads as AI.

Contents: tier system and word tables; template phrases; content patterns;
structure patterns (binary contrasts, false agency, other tells);
communication and filler patterns; style and formatting; context-profile
matrix; severity tiers; over-polishing warning; 2026 structural tells; fact
preservation and evidence-bound mode; rewrite-quality guards; client-
deliverable profile; research-paper profile; email and plain-message profile.
The 2026-09 additions (post-delve register, displacement table, self-edit
failures) live in `current-register.md`.

## Tier system (vocabulary)

- **Tier 1: always replace.** Appears 5-20x more in AI text than human text.
- **Tier 2: flag in clusters.** Fine alone; 2+ in one paragraph is a strong tell.
- **Tier 3: flag by density.** Normal words AI overuses. Flag only when they hit
  roughly 3%+ of the text.

### Tier 1: always replace

| Replace | With |
|---|---|
| delve / delve into | explore, dig into, look at |
| deep dive / dive into | look at, examine, explore |
| unpack / unpacking | explain, break down, walk through |
| landscape (metaphor) | field, space, industry, world |
| realm | area, field, domain |
| tapestry | (describe the actual complexity) |
| paradigm | model, approach, framework |
| embark | start, begin |
| beacon | (rewrite entirely) |
| nestled | is in, sits in, is located in |
| testament to | shows, proves, demonstrates |
| robust | strong, reliable, solid |
| comprehensive | thorough, complete, full |
| cutting-edge | latest, newest, advanced |
| leverage (verb) | use |
| utilize | use |
| pivotal | important, key, critical |
| underscores | highlights, shows |
| meticulous / meticulously | careful, detailed, precise |
| seamless / seamlessly | smooth, easy, without friction |
| game-changer / game-changing | say what specifically changed and why |
| hit differently / hits different | say what changed, or cut |
| watershed moment | turning point, shift |
| marking a pivotal moment | (state what happened) |
| the future looks bright | (cut, or say something specific) |
| only time will tell | (cut, or say something specific) |
| vibrant | (describe what makes it active, or cut) |
| thriving / bustling | growing, busy, active (or cite a number) |
| despite challenges... continues to thrive | name the challenge and response, or cut |
| showcasing | showing, demonstrating (or cut the clause) |
| intricate / intricacies / complexities | complex, detailed (name the specific complexity) |
| ever-evolving | changing, growing (or describe how) |
| enduring | lasting, long-running (or cite how long) |
| daunting | hard, difficult, challenging |
| holistic / holistically | complete, full, whole |
| actionable | practical, useful, concrete |
| impactful | effective, significant (or describe the impact) |
| learnings | lessons, findings, takeaways |
| thought leader / thought leadership | expert, authority (or describe the contribution) |
| best practices | what works, proven methods, standard approach |
| at its core | (cut, just state the thing) |
| synergy / synergies | (describe the actual combined effect) |
| interplay | relationship, connection, interaction |
| in order to | to |
| due to the fact that | because |
| serves as | is |
| features (verb) | has, includes |
| boasts | has |
| presents (inflated) | is, shows, gives |
| commence | start, begin |
| ascertain | find out, determine, learn |
| endeavor | effort, attempt, try |
| keen (as intensifier) | interested, eager (or cut) |
| symphony (metaphor) | (describe the actual coordination) |
| embrace (metaphor) | adopt, accept, use, switch to |
| garner | get, earn, attract |
| adept | skilled, good at |
| commendable | (say what was done well, or cut) |
| noteworthy | (state the fact; if it is notable, the fact shows it) |

The last four rows come from Liang et al. 2024 ("Mapping the increase of LLM
usage in scientific papers"), the words with the sharpest post-2022 rise in
arXiv abstracts. The rest of that list (realm, intricate, showcasing, pivotal,
meticulously, delve, underscores, boasts) was already in the table.

### Tier 2: flag when 2+ appear in one paragraph

| Replace | With |
|---|---|
| harness | use, take advantage of |
| navigate / navigating | work through, handle, deal with |
| foster | encourage, support, build |
| elevate | improve, raise, strengthen |
| unleash | release, enable, unlock |
| streamline | simplify, speed up |
| empower | enable, let, allow |
| bolster | support, strengthen, back up |
| spearhead | lead, drive, run |
| resonate / resonates with | connect with, appeal to, matter to |
| revolutionize | change, transform, reshape |
| facilitate / facilitates | enable, help, allow, run |
| underpin / underpinnings | support, basis, foundation |
| nuanced | specific, subtle (or name the actual nuance) |
| crucial | important, key, necessary |
| multifaceted | (describe the actual facets, or cut) |
| ecosystem (metaphor) | system, community, network, market |
| myriad / plethora | many, numerous (or give a number) |
| encompass | include, cover, span |
| catalyze | start, trigger, accelerate |
| reimagine | rethink, redesign, rebuild |
| galvanize | motivate, rally, push |
| augment | add to, expand, supplement |
| cultivate | build, develop, grow |
| illuminate | clarify, explain, show |
| elucidate | explain, clarify, spell out |
| juxtapose | compare, contrast, set side by side |
| transformative / transformation | (describe what changed and how) |
| paradigm-shifting | (describe what actually shifted) |
| cornerstone | foundation, basis, key part |
| paramount | most important, top priority |
| poised (to) | ready, set, about to |
| burgeoning | growing, emerging (or cite a number) |
| nascent | new, early-stage, emerging |
| quintessential | typical, classic, defining |
| overarching | main, central, broad |

### Tier 3: flag only at high density

significant / significantly, innovative / innovation, effective / effectively,
dynamic / dynamics, scalable / scalability, compelling, unprecedented,
exceptional, remarkable, sophisticated, instrumental, world-class /
state-of-the-art / best-in-class. These are normal words. Flag only when the
text leans on them in place of specifics: numbers, comparisons, examples.

### Template / slot-fill phrases

If a phrase has a blank where any noun or adjective would fit and still sound the
same, it was generated, not written.

- "a [adjective] step toward [adjective] X" → name the specific capability or outcome.
- "Whether you're [X] or [Y]" → false breadth. Pick the real audience, or cut.
- "I recently had the pleasure of [verb]-ing" → just say what happened.

## Content patterns

Structural tells, not single words. Examples are illustrative.

1. **Significance / legacy inflation.** Routine facts dressed as history.
   Before: "Founded in 1989, it marked a pivotal moment in the evolution of
   regional governance." After: "Founded in 1989 to publish regional statistics."
2. **Notability name-dropping.** Stacked prestige citations.
   Before: "Her work has appeared in the Times, the BBC, and the FT." After: "In
   a 2024 Times interview she argued that regulation should target outcomes."
3. **Superficial -ing analyses.** Present-participle tails faking depth.
   Before: "...the palette resonates with the region, symbolizing renewal,
   reflecting community pride." After: state the actual reason the colors were
   chosen, or cut.
4. **Promotional / tourism-brochure prose.** Before: "Nestled in the breathtaking
   foothills, a vibrant hub of culture." After: "A town in the Gonder region with
   a weekly market."
5. **Vague attributions / weasel words.** "Experts believe", "studies show",
   "industry leaders agree" with nobody named. Cite the specific source or drop
   the attribution and state the claim.
6. **Formulaic challenges sections.** "Despite challenges, it continues to
   thrive." Name the actual challenge and the actual response, or cut.
7. **Copula avoidance.** Before: "Gallery 825 serves as the exhibition space and
   boasts 3,000 square feet." After: "Gallery 825 is the exhibition space. It has
   3,000 square feet." Default to is / has.
8. **Negative parallelism and tailing negations.** "Not just a song, it's a
   statement." "The options come from the item, no guessing." State it directly.
   FATAL on a byline whose voice spec bans it (zero, not "max one"; see
   protect-list.md).
9. **Rule-of-three overuse.** Forced triads to sound comprehensive. Vary the
   grouping: two items, four, or a full sentence.
10. **Synonym cycling (elegant variation).** "developers... engineers...
    practitioners... builders" in one paragraph. Repeat the clearest word.
11. **False ranges.** "from the Big Bang to dark matter", "from ancient
    civilizations to modern startups." List the actual topics, or pick the one
    that matters.
12. **Novelty inflation.** Treating an established concept as a discovery ("she
    coined the term", "a failure mode nobody names"). Describe what the person
    did with the concept, not that they invented it.
13. **Emotional flatline.** Claiming the feeling instead of earning it ("what
    surprised me most", "I was fascinated to discover"). If it is surprising, the
    content should show it. Otherwise cut the claim.

## Structure patterns (the #1 signal)

Weight these above vocabulary. Detection tooling scores structural regularity
higher than word choice.

### Binary contrasts (the fatal family)

Telegraphed reversals. State the positive claim directly. Variants to catch:

| Pattern | |
|---|---|
| "Not because X. Because Y." | "It's not just X, it's Y." |
| "X isn't the problem. Y is." | "The answer isn't X. It's Y." |
| "It feels like X. It's actually Y." | "The question isn't X. It's Y." |
| "stops being X and starts being Y" | "doesn't mean X, but actually Y" |
| "is about X but not Y" | "not just X but also Y" |
| "more than just X" / "goes beyond X" | "Less X, more Y." |

FATAL on a byline whose voice spec bans it. Drop the negation; assert Y.

### False agency (inanimate things doing human verbs)

AI uses this to avoid naming the actor. Name the human, or use "you".

| Before | Why it's wrong |
|---|---|
| "the complaint becomes a fix" | someone fixed it |
| "a bet lives or dies in days" | someone ships or kills it |
| "the decision emerges" | someone decides |
| "the culture shifts" | people change behavior |
| "the data tells us" | someone read it and drew a conclusion |
| "the market rewards" | buyers pay for things |

### Other structure tells

- **Negative listing.** "Not a X... Not a Y... A Z." State Z; skip the runway.
- **Dramatic fragmentation.** "Noun. That's it. That's the thing." Complete
  sentences instead. (Protected in deliberate doses on bylines that list fragment
  beats and short-short-long stacking. See protect-list.md.)
- **Rhetorical setups.** "What if...?", "Here's what I mean:", "Think about it:",
  "And that's okay." Make the point; let the reader conclude.
- **Narrator-from-a-distance.** "Nobody designed this." "People tend to..." Put
  the reader in the room: "You don't sit down one day and decide to..."
- **Passive voice / subjectless fragments.** "Mistakes were made." "No config
  needed." Find the actor; move them to the front.
- **Wh- openers as a crutch.** "What makes this hard is..." becomes "The
  constraint is..." (Protected as hand-the-ball-back closers on bylines that list
  them; see protect-list.md.)
- **Sentence- and paragraph-length uniformity.** The metronome. Mix short (3-8
  words) with long (20+). Some one-sentence paragraphs. Read-aloud test: if a TTS
  engine could read it without sounding odd, it is too uniform.

## Communication / filler patterns

- **Chatbot artifacts.** "I hope this helps", "Certainly!", "Great question",
  "Feel free to reach out", "Let me know if you need anything." Strip entirely.
- **Cutoff disclaimers.** "As of my last update", "While specific details are
  limited based on available information." Find the fact or remove the hedge.
- **Sycophancy.** "You're absolutely right!", "That's a really insightful
  observation." Remove.
- **Acknowledgment loops.** "To answer your question", "The question of whether",
  restating the prompt before answering. Just answer.
- **Confidence-calibration adverbs.** "Notably", "Interestingly", "Importantly",
  "Surprisingly", "Undoubtedly." Let the fact carry its own weight. Flag by
  density (one in 2,000 words is fine; three in 500 is a tell).
- **Filler phrases.** "It is important to note that" → state it. "In terms of" →
  rewrite. "The reality is that" → cut.
- **Excessive hedging.** "It could potentially possibly be argued that it might."
  → "It may." (A documented hedge-then-commit signature is the exception; see
  protect-list.md.)
- **Generic positive conclusions.** "The future looks bright", "exciting times
  ahead", "a step in the right direction." Cut, or make it specific.
- **Reasoning-chain artifacts.** "Let me think step by step", "Breaking this
  down", "Step 1:", "Here's my thought process." Scaffolding the reader does not
  need. State the conclusion, then the evidence.
- **Signposting and "let's" openers.** "Let's dive in", "Let's explore", "Here's
  what you need to know", "without further ado." Announcing the move instead of
  making it. Start with the point.

## Style / formatting patterns

- **Em dashes.** Target zero. Catch both the Unicode glyph and the double-hyphen
  substitute, in headings too. Rewrite with commas, periods, colons, or
  parentheses.
- **Boldface overuse.** One bolded phrase per major section at most. If it
  matters enough to bold, lead the sentence with it instead.
- **Inline-header vertical lists.** Bullets that open with a bold header
  repeating itself ("Performance: performance improved by..."). Strip the header,
  write the point, or make it a paragraph.
- **Title case in headings.** Use sentence case for subheadings. Title case only
  for the top title, if at all.
- **Emoji in headers.** Remove. Social posts may use one or two at end-of-line.
- **Quotes.** Straight quotes on plain-text and markdown surfaces (chat,
  email, code, README). Curly quotes are correct typography inside a Word or
  PowerPoint deliverable; do not straighten them there.
- **Excessive bullet lists.** Convert bullet-heavy prose to paragraphs. Bullets
  only for genuinely list-like content.
- **Fragmented headers.** A heading followed by a one-line restatement of itself
  before the real content. Cut the warm-up line.

## Context-profile matrix

Adjusts strictness per surface. Rules not listed apply at full strength
everywhere.

| Rule | linkedin | blog | technical-blog | investor-email | docs | casual |
|------|----------|------|----------------|----------------|------|--------|
| Em dashes | 2/post OK | strict | strict | strict | relaxed | skip |
| Bold overuse | hooks OK | strict | strict | strict | relaxed | skip |
| Emoji in headers | 1-2 end-of-line | strict | strict | strict | skip | skip |
| Excessive bullets | skip | strict | relaxed | strict | skip | skip |
| Hedging | strict | strict | relaxed | strict | relaxed | skip |
| Word tables | strict | strict | partial | strict | relaxed | P0 only |
| Promotional | relaxed | strict | strict | extra strict | strict | skip |
| Copula avoidance | skip | strict | relaxed | strict | skip | skip |
| Generic conclusions | skip | strict | strict | extra strict | skip | skip |

Technical-blog word exceptions (legit technical meaning): robust, comprehensive,
seamless, ecosystem, leverage (platform/API), facilitate, underpin, streamline.
Still flag: delve, tapestry, beacon, embark, testament to, game-changer, harness.

Auto-detect when no profile is passed: <300 words + hashtags = linkedin; code
blocks = technical-blog; salutation + fundraising = investor-email; step-by-step
+ params = docs; steering-committee / governance / pre-read vocabulary =
client-deliverable; "we show" / "et al." = research-paper; salutation without
fundraising = email; else blog. The scanner's `--profile` flag takes the same
names.

Matt's stricter-than-floor settings override this matrix wherever they
conflict: em dashes are zero on every surface, including the LinkedIn cell
above; tricolons on deliverables are Tier 1. See protect-list.md.

## Severity tiers (for triage)

- **P0, credibility killers:** cutoff disclaimers, chatbot artifacts, vague
  attributions, significance inflation on routine events.
- **P1, obvious AI smell:** Tier 1 word hits, template phrases, "let's" openers,
  synonym cycling, formulaic openings, bold overuse, em dash frequency, the
  binary-contrast family.
- **P2, stylistic polish:** generic conclusions, rule of three, uniform paragraph
  length, copula avoidance, transition phrases.

Quick passes do P0+P1. Full audit covers all three.

## Over-polishing warning

Applying every rule at maximum strictness sands writing back toward the uniform
statistical profile that reads as AI. Deliberate fragments, the occasional "And"
opener, idiosyncratic word choice, and uneven pacing are what keep text human.
The goal is to sound like a person, not like clean prose. When in doubt on a
byline, defer to the protect list.

## 2026 structural tells (merged from adenaufal/anti-slop-writing)

The legacy vocabulary tells (delve, tapestry, vibrant) are being trained out of
newer models. Their absence proves nothing. What survives model updates is
structural. Weight this section at least as heavily as the vocabulary tiers.

- **Cadence uniformity is the current #1 tell.** Sentences landing at 18-24
  words, one after another. Two 30-second tests an editor applies by eye:
  1. First word of each sentence in a paragraph. More than half starting with
     "The", "This", "It", or "In" reads as LLM-assisted.
  2. Three or more consecutive sentences in the 17-23 word band. Same verdict.
- **The bimodal trap.** Newer models fake burstiness by mechanically alternating
  a punchy fragment with a very long sentence. Repeated, that seesaw is itself a
  fingerprint. Human writing clusters around medium length with occasional
  irregular swings in both directions. Vary irregularly: two mediums, a
  fragment, a long one, a medium. Not a metronome, not a seesaw.
- **The parataxis / fragment tension, resolved.** Chains of blunt declaratives
  ("Short sentence. Then another. Then another.") read as AI. So does zero
  fragmentation. The rule: fragments are a spice, not a syntax. Most related
  thoughts should connect through subordination, conjunction, or punctuation
  that shows how the ideas relate (causation, contrast, qualification).
- **Over-fragmentation at paragraph level.** Newer models over-correct into many
  tiny 1-2 sentence paragraphs plus bullet spam. A real writer runs a paragraph
  to 7-8 sentences when the argument needs it. Uniform brevity is as
  machine-like as uniform length.
- **Punctuation displacement into colons.** With em dashes suppressed, the same
  interruptive rhythm moves into colons (measured at ~4x human rate in recent
  Claude output) and semicolons (~3x). A colon every few sentences is a tell
  even with zero dashes. Replace most with periods. En dashes: same ban as em
  dashes; plain hyphens only for ranges and compound adjectives.
- **High-signal current vocabulary:** "ensuring / ensures" as padding (the
  strongest single AI word in current corpus data), "plays a crucial role in
  shaping", "conversely" (heavily overrepresented; use "but" or restructure),
  "in essence / essentially / fundamentally" (cut), hedging-adverb stacks
  ("typically... often... potentially" in one passage; commit or quantify).
- **Break the clean arc.** Claim, evidence, implication, every paragraph, is an
  AI signature. Start some paragraphs mid-thought or with a specific detail.
  End some before the expected "so what". Mix in a genuine question or an
  imperative; AI writes almost exclusively in declaratives.

Continued in `current-register.md` (dated 2026-09): the post-delve vocabulary,
sincerity intensifiers, label-colon sentences, which-tail verdicts, aphoristic
closers, the displacement table, and the self-edit failure modes.

## Fact preservation and evidence-bound mode (merged from ehmo/slopbeth)

A de-slop pass that changes what the text claims has failed, whatever it did
for style. Run this contract before any rewrite.

**Preserve (lock before editing):** named entities, numbers, thresholds, dates,
units, prices, versions; URLs, citations, quotations, identifiers; causal
claims and scope limits; explicit uncertainty and qualifiers that carry
meaning; the author's stance. Critical failures that block delivery: wrong
number, changed entity, removed necessary qualifier, false causality, invented
evidence, dropped constraint.

**Evidence-bound mode.** When the source text is vague and no supporting
material is present, the rewrite may tighten claims, remove hype, improve
rhythm, and name what proof is missing. It may NOT invent the missing story:
no invented metrics, owners, deadlines, features, customer facts, workflow
steps, or outcomes. "Faster decisions" does not become "cut review time 40%"
unless a source says so. Prefer bracketed slots (`[metric]`, `[owner]`) or a
direct request for the missing fact. Density cannot come from fictional
specificity; fabricated specificity is worse than honest vagueness.

**Label added detail** when concreteness is introduced: `source` (present in
the draft or materials), `inference` (follows narrowly without new facts),
`placeholder` (user must verify). Anything that would read as fact but is not
supported stays out.

## Rewrite-quality guards (merged from slopbeth + jalaalrd)

- **The noun-swap test.** If a sentence survives transplantation to another
  topic with only noun swaps, it is filler regardless of clean vocabulary. It
  needs a fact, scene, number, decision, contradiction, or consequence.
- **Do not replace slop with a new formula.** Clipped aphorisms, tidy triads,
  forced contrast, dramatic fragments, and generic consultant voice are the
  standard failure mode of naive de-slopping. Check the rewrite against the
  same rules as the original.
- **Over-editing is a failure.** Already-strong human text gets a light edit or
  none. "Leave this alone" is a valid and sometimes correct output. Changing
  more than about a fifth of a good draft without clear cause is a warning sign.
- **Tool artifacts.** Scan for leaked machinery: `oaicite`, `turn0search`,
  `contentReference`, `utm_source=`, placeholder URLs, stray markdown asterisks
  rendering literally. Strip on sight.
- **No markdown in plain-text surfaces.** Headers, bold, and bullet syntax in
  emails, texts, and chat messages are an instant tell. Prose only.
- **Span-level critique.** In detect mode on long or high-stakes text, point to
  exact spans with a label and reason. A critique that cannot point to a span
  is too vague to act on.

## Client-deliverable context profile (Matt)

Applies on top of the matrix above for consulting deliverables (decks, memos,
reports, proposals). Strictest profile in the skill.

- All floor rules at full strength. Binary contrast: FATAL (zero, not max one).
- Em dashes inside sentences: zero. (Standing rule predating this skill.)
- Sentence case for all headings and slide titles.
- No vendor, product, or consultant names in client-facing deck body content.
- Tricolons: treat as Tier 1. Filler triads are a known personal irritant.
- Parallel-construction slogans and meta-slogans about how rules or frameworks
  "should function": cut. State what the thing does.
- Abstract nouns where a verb works ("the implementation of" vs "implementing",
  "alignment" vs "agree on"): rewrite with the verb.
- Repeated structural templates across list items (every bullet shaped
  "X: which does Y, enabling Z"): vary or collapse to prose.
- Anything that reads like a LinkedIn post or strategy-memo aphorism: cut.
- Engagement-specific style rules (naming order conventions, terminology
  fences) live in the engagement's own notes; load them alongside this profile
  when working those documents. Example from current work: RWJBH before
  Rutgers, and never "partnership/JV" for the affiliation itself.

## Research-paper context profile

For manuscripts, preprints, methods notes, and technical reports. The register
has conventions that look like tells elsewhere and are not tells here, and one
class of tell that is specific to it.

- Hedges are required here. "suggests", "is consistent with", "we
  observe", "may indicate" carry meaning about evidence strength. Do not
  commit them to "shows" or "proves". The hedge-stack rule applies at a higher
  threshold; a hedge that repeats inside one sentence is still cut.
- "significant" is reserved for statistical significance with a test behind
  it. Anywhere else, replace it with the size of the effect.
- Passive voice is acceptable in methods where the actor is the procedure.
  Elsewhere, "we" is the actor.
- Terms of art stay: robust (robustness), comprehensive (of a search or
  survey), dynamics, paradigm in its Kuhnian sense, novel when a claim of
  novelty is being made and defended, orthogonal, non-trivial in its
  mathematical sense. The scanner exempts these on this profile.
- Conventions that are fine: "In this paper, we", "To the best of our
  knowledge" once, numbered contribution lists when the count is real, colons
  before enumerations and definitions, heading case per the venue.
- The tells specific to this register (Liang et al. 2024 measured these rising
  in LLM-assisted abstracts): delve, intricate, meticulously, commendable,
  pivotal, underscores, showcasing, realm, garner, adept, notably, crucial.
  Also "Interestingly," and "Importantly," at sentence starts, "It is worth
  noting that", significance inflation in the introduction ("a paradigm
  shift", "transformative"), and a contributions list whose third item
  restates the first two.
- Numbers, units, thresholds, seeds, sample sizes, and citation keys are
  locked facts. Every number in prose must match its table. Rewrite nothing
  in a results sentence without re-checking the table it cites.
- Em dashes: zero, per the standing rule, even where the venue tolerates
  them. Colons before enumerations are fine; colons as asides are not.

## Email and plain-message context profile

For anything sent as an email, chat message, or text. The surface is plain
text, so the formatting rules become P0.

- No markdown: no headers, bold, bullet syntax, or backticks. Prose only. A
  list of three items goes in a sentence, or on three short lines without
  bullets.
- No dashes. Straight quotes.
- No chatbot openers or closers: "I hope this finds you well", "I wanted to
  reach out", "Please don't hesitate", "Happy to X if useful", "Let me know
  if you have any questions". Open with the point. Close with the ask or the
  next step and a plain sign-off.
- One ask per email where possible, stated as a question the reader can
  answer in one line.
- Length follows the reader. A busy executive gets the slip, the reason, and
  the question in under 120 words.
- The register is the sender's, not the assistant's. No "Certainly", no
  "Great question", no "I'd be happy to".
