# Living corpus

Dated AI tells caught in the wild, each tagged with the mechanism that produces
it. This is the part of the skill that compounds. Grow it with the ingest flow
(`ingestion.md`). Entry format is defined there.

Tells age. Re-tier or retire entries as the models change. Note the date so we
can see how fast a tell rises and falls.

---

### anti-em-dash displacement

- Added: 2026-06-21  |  Tier: context-dependent
- Mechanism: displacement
- Context: A second-order tell. Now that every slop skill flags em dashes, models
  shunt the same "punchy aside" rhythm into colons and comma-spliced
  appositives. The punctuation changed; the metronome did not. Flag when the
  colon-or-appositive rhythm repeats across several sentences, not on a single
  clean use.
- Before: "The fix is simple: stop. It works, a clean little loop, every time."
- After: "The fix is simple. Stop. It runs as a clean loop every time."
- Rule: Removing em dashes is not enough. Check whether the same interruptive
  rhythm just moved into colons or paired commas. Vary the sentence shape, not
  only the punctuation mark.
- Source: upstream (kjmagnan1s/anti-slop)

---

### concession reflex

- Added: 2026-06-21  |  Tier: 2
- Mechanism: reward-tuning
- Context: "To be fair, X. That said, Y." dropped in to simulate balance when
  there is no real counterpoint. Fine once when there is a genuine tradeoff; a
  tell when it is reflexive.
- Before: "To be fair, the tool has limitations. That said, it is still useful."
- After: "The tool is useful for X and weak at Y." (or just state the limitation)
- Rule: Do not manufacture balance. If there is a real counterpoint, name it
  specifically. If there is not, drop the concession and make the claim.
- Source: upstream (kjmagnan1s/anti-slop)

---

### em-dash aphorism

- Added: 2026-07-19  |  Tier: 1 (FATAL on client deliverables)
- Mechanism: reward-tuning
- Context: The "X is Y, not Z" verdict sentence, classically set off with an em
  dash, now appearing with colons or commas after dash suppression. A subtype
  of the binary-contrast family with a specific rhythm: short setup, punchy
  reversal, delivered as if profound.
- Before: "Governance is a routing problem, not a committee problem."
- After: State how routing solves it, or name the specific committee failure.
  Assert the positive claim with its evidence; skip the reversal.
- Rule: Any aphoristic contrast delivered as a verdict gets rewritten as a
  direct claim with a specific. Repeated corrections on this since early 2026.
- Source: my own output, flagged repeatedly by Matt

---

### parallel-construction slogan

- Added: 2026-07-19  |  Tier: 1 on deliverables, 2 elsewhere
- Mechanism: instruction-tuning
- Context: Matched-clause slogans ("One enterprise. Two institutions. Shared
  research.") and mirrored headline pairs. Reads as advertising copy inside a
  working document.
- Before: "Aligned governance. Shared infrastructure. Unified science."
- After: Name the governance body, the infrastructure, and the operating
  change, in a sentence.
- Rule: If clauses share a template and could be reordered without loss, they
  carry no information order. Rewrite as prose with specifics.
- Source: my own output, flagged repeatedly by Matt

---

### filler tricolon

- Added: 2026-07-19  |  Tier: 1 on deliverables
- Mechanism: instruction-tuning
- Context: Triads generated to sound complete rather than because three things
  exist. Distinct from a genuine three-item list, which is fine.
- Before: "clear, consistent, and repeatable"
- After: Pick the one adjective that is true and necessary, or cite the
  behavior that makes it so.
- Rule: For every triad, ask whether the content genuinely has three members.
  If any member could be dropped without anyone noticing, it is filler.
- Source: my own output, flagged repeatedly by Matt

---

### abstract-noun stack

- Added: 2026-07-19  |  Tier: 2
- Mechanism: pretraining-register
- Context: Nominalizations where a verb works: "the implementation of",
  "the establishment of", "alignment", "enablement", "operationalization".
  Consulting register amplifies this; client deliverables are the highest-risk
  surface.
- Before: "the operationalization of the routing framework requires alignment"
- After: "to run the routing framework, the co-chairs need to agree on"
- Rule: When a sentence's weight sits in -tion/-ment/-ance nouns, find the verb
  and the actor and rewrite around them.
- Source: my own output, flagged repeatedly by Matt

---

### meta-slogan about rules

- Added: 2026-07-19  |  Tier: 1
- Mechanism: reward-tuning
- Context: Sentences that describe how a rule or framework should function
  instead of stating the rule. Common in governance and process writing.
- Before: "The rulebook should serve as a living document that evolves with
  the enterprise."
- After: State the update trigger and who owns it: "The co-chairs review the
  rulebook each quarter and after any routing dispute."
- Rule: If a sentence is about the document rather than the operation, replace
  it with the operation.
- Source: my own output, flagged repeatedly by Matt

---

### repeated list-item template

- Added: 2026-07-19  |  Tier: 1 on deliverables, 2 elsewhere
- Mechanism: instruction-tuning
- Context: Every bullet in a list shaped identically ("**X:** which does Y,
  enabling Z"), often with a bolded lead-in. The template is visible before
  the content is read.
- Before: Four bullets each reading "[Noun]: [capability], enabling [outcome]."
- After: Vary the shape per item, or collapse the list into a paragraph that
  states what actually differs between the items.
- Rule: If the list's syntax is more regular than its content, the syntax is
  lying. Break the template or drop the list.
- Source: my own output, flagged repeatedly by Matt

---

### two-word abstract balance

- Added: 2026-07-19  |  Tier: 1
- Mechanism: reward-tuning
- Context: Paired abstractions posing as insight: "real but bounded",
  "necessary but insufficient", "ambitious but achievable", "material yet
  manageable". Corporate-MBA register. The two words sound balanced and feel
  analytical; neither says anything specific.
- Before: "The integration risk is real but bounded."
- After: Name the risk and what bounds it: "Two systems double-book the same
  coordinators today; the shared calendar closes that in Q3."
- Rule: Any "[abstract adjective] but/yet [abstract adjective]" pair gets
  replaced with the specific thing on each side of the tension.
- Source: matt-voice tone spec ("What the voice never does")

---

### performed honesty

- Added: 2026-07-19  |  Tier: 1
- Mechanism: assistant-persona
- Context: Phrases that mark honesty: "to be honest", "the honest answer is",
  "let me be clear", "the real truth", "candidly". Marking honesty implies
  the possibility of dishonesty elsewhere. Related to but distinct from the
  fake-authenticity openers ("Here's the thing:"); this family specifically
  performs candor.
- Before: "To be honest, the timeline is at risk."
- After: "The timeline is at risk."
- Rule: Delete the honesty marker and keep the claim. If the claim feels
  weaker without it, the claim needs evidence, not a marker.
- Source: matt-voice tone spec ("Never perform honesty")

---

### speaker-centering contrast

- Added: 2026-07-19  |  Tier: 2
- Mechanism: reward-tuning
- Context: Framings that put the speaker's reactions at the center of the
  sentence: "I have no problem with X. I have a problem with Y." "I am not
  quite sure what to do with that." Doubles as a binary-contrast variant when
  paired. Distinct from the protected self-aware admission, which discloses
  the speaker's state when the state itself is the relevant fact.
- Before: "I have no problem with the concept. I have a problem with the
  rollout."
- After: "The concept is sound; the rollout is the issue worth raising."
- Rule: Restructure to put the subject at the center unless the speaker's
  state is genuinely the point (see protect-list, self-aware admission).
- Source: matt-voice tone spec ("Never elevate the speaker's reactions")

---

## 2026-09 batch: model self-report

The entries below were added during the 2026-09 revision from the current
model's account of its own tells, not from a catch in Matt's work. They are
provisional. Promote one to a confirmed tier the first time it is caught in a
real draft; retire it if a year passes without a catch. Overview and the
displacement table are in `current-register.md`.

---

### post-delve register cluster

- Added: 2026-09-02  |  Tier: 2 (flag in clusters)
- Mechanism: register-migration
- Context: The vocabulary newer models moved into after the Tier 1 list was
  trained out: load-bearing, non-trivial, orthogonal, the delta, meaningfully,
  crisp, surface (verb), reach for, lean on, the ask, through-line, forcing
  function, earns its place, does real work. Each is fine alone and some are
  exact in engineering writing. Two in a paragraph reads as this year's model.
- Before: "The load-bearing change is non-trivial, but it's where the
  leverage is."
- After: "The change that matters is moving coverage analysis ahead of
  budget. It is hard because two offices have to swap order."
- Rule: Replace with the plain word or, better, with the specific the idiom
  stood in for. Full list in current-register.md, section 1.
- Source: model self-report (Fable 5.1, 2026-09-02 revision); provisional until caught in Matt's work

---

### sincerity intensifier

- Added: 2026-09-02  |  Tier: 3 (by density)
- Mechanism: reward-tuning
- Context: actual, actually, real, really, genuine, genuinely, truly,
  literally in front of nouns and claims: "the actual risk", "what really
  matters", "genuinely useful". Marking sincerity implies the alternative was
  on the table. A cousin of performed honesty.
- Before: "The actual problem is that nobody genuinely owns the intake queue."
- After: "Nobody owns the intake queue."
- Rule: Above roughly six per thousand words, cut them all. If a sentence
  weakens without its "actually", it needed evidence, not an adverb.
- Source: model self-report (Fable 5.1, 2026-09-02 revision); provisional until caught in Matt's work

---

### label-colon sentence

- Added: 2026-09-02  |  Tier: 1
- Mechanism: displacement
- Context: "The catch: ...", "Translation: ...", "The result: ...",
  "Concretely: ...". Third-generation displacement of the em-dash aside: dash
  to colon to label.
- Before: "The catch: nobody owns the intake queue."
- After: "Nobody owns the intake queue."
- Rule: A label is not a subject. Give the sentence a subject and a verb.
- Source: model self-report (Fable 5.1, 2026-09-02 revision); provisional until caught in Matt's work

---

### which-tail verdict

- Added: 2026-09-02  |  Tier: 1
- Mechanism: displacement
- Context: A verdict appended as a relative clause: ", which is exactly the
  problem." ", which is the point." ", which matters." It avoids having to
  argue the verdict.
- Before: "Both offices report to different chairs, which is exactly the
  problem."
- After: "Both offices report to different chairs. No one can order them to
  sequence the work."
- Rule: If the tail asserts a verdict, argue it in its own sentence or cut it.
- Source: model self-report (Fable 5.1, 2026-09-02 revision); provisional until caught in Matt's work

---

### aphoristic closer paragraph

- Added: 2026-09-02  |  Tier: 1 on deliverables, 2 elsewhere
- Mechanism: reward-tuning
- Context: A final one-line paragraph delivering the verdict: "That's the
  whole game." "That is the job." It replaced "In conclusion" after that was
  banned; same function, more swagger.
- Before: "...the co-chairs meet monthly. That is the whole game."
- After: "...the co-chairs meet monthly and settle disputes in that meeting."
- Rule: End on the last real point. If the closer could be pasted onto
  another document, cut it.
- Source: model self-report (Fable 5.1, 2026-09-02 revision); provisional until caught in Matt's work

---

### restatement chain

- Added: 2026-09-02  |  Tier: 1 when two or more, 2 for one
- Mechanism: reward-tuning
- Context: "In other words," "Put simply," "Which is to say," following a
  claim that was already clear. Helpfulness training rewards saying it again.
- Before: "Activation depends on sequencing. In other words, the order of
  steps matters more than the steps. Put simply: sequence first."
- After: "Activation depends on the order of the steps more than on the
  steps."
- Rule: Keep the clearer of the two versions and delete the other; the
  announcer goes with it.
- Source: model self-report (Fable 5.1, 2026-09-02 revision); provisional until caught in Matt's work

---

### one-word opener

- Added: 2026-09-02  |  Tier: 2
- Mechanism: instruction-tuning
- Context: "Yes." "No." "Correct." "Short answer: no." as the first sentence
  of a reply or section. Directness training that became a tic; reads as a
  chatbot even inside a memo.
- Before: "No. The committee cannot approve the budget without the coverage
  analysis."
- After: "The committee cannot approve the budget without the coverage
  analysis."
- Rule: Begin with the substantive sentence. The yes or no is inside it.
- Source: model self-report (Fable 5.1, 2026-09-02 revision); provisional until caught in Matt's work

---

### enumerator opener

- Added: 2026-09-02  |  Tier: 2
- Mechanism: instruction-tuning
- Context: "Two things." "A few notes." "Three observations." as an opener.
  Announces the count, then the list.
- Before: "Two things. The budget is late, and the coverage analysis is not
  started."
- After: "The budget is late, and the coverage analysis is not started."
- Rule: Start with the first thing. A visible count needs no announcement.
- Source: model self-report (Fable 5.1, 2026-09-02 revision); provisional until caught in Matt's work

---

### offer closer

- Added: 2026-09-02  |  Tier: 1 in deliverables and email, 2 elsewhere
- Mechanism: assistant-persona (displacement of "let me know if")
- Context: "Happy to draft the summary if useful." "If helpful, I can pull
  the numbers." "Want me to expand on this?" Same function as the banned
  closer, reworded as an offer.
- Before: "Happy to pull the site-level numbers if that would help."
- After: cut, or ask for the decision: "Do you want the site-level numbers
  in the pre-read?"
- Rule: A closer that offers is cut. If an action is needed, ask for the
  decision that gates it.
- Source: model self-report (Fable 5.1, 2026-09-02 revision); provisional until caught in Matt's work

---

### asker-evaluation opener

- Added: 2026-09-02  |  Tier: 1
- Mechanism: assistant-persona (displacement of "great question")
- Context: "Good instinct." "Right call." "Fair push." "Worth asking."
  Compliments the asker's judgment before answering.
- Before: "Good instinct. The Q2 numbers are the weak point."
- After: "The Q2 numbers are the weak point."
- Rule: Any opener that evaluates the asker or the ask is cut.
- Source: model self-report (Fable 5.1, 2026-09-02 revision); provisional until caught in Matt's work

---

### flattened rhythm after de-slopping

- Added: 2026-09-02  |  Tier: context-dependent
- Mechanism: self-edit
- Context: The rewrite obeys every rule and reads like a telegram: every
  sentence short, no because or although or which, mean length down by a
  third. Produced by treating the colon and dash rules as zero-tolerance and
  by the parataxis reflex.
- Before: "The budget is late. The coverage analysis is not started. The
  sites are waiting. The committee meets Thursday."
- After: "The budget is late and the coverage analysis has not started, so
  the sites are waiting on both before the committee meets Thursday."
- Rule: Compare mode flags a mean sentence length that fell by a third or
  subordinators that halved. Reconnect related thoughts; fragments are a
  spice.
- Source: model self-report (Fable 5.1, 2026-09-02 revision); provisional until caught in Matt's work

---

### tail-qualifier stack

- Added: 2026-09-02  |  Tier: 2
- Mechanism: uncertainty-conditioning
- Context: "at least for now", "for what it's worth", "as far as I can tell",
  "in my experience" hung off the ends of sentences, two or more in a
  paragraph.
- Before: "Two sites have not reported, at least for now. The numbers look
  fine, for what it's worth."
- After: "Two sites have not reported. The other twelve look fine."
- Rule: One hedge, inside the claim, or commit.
- Source: model self-report (Fable 5.1, 2026-09-02 revision); provisional until caught in Matt's work

---

### parenthetical wink

- Added: 2026-09-02  |  Tier: 1
- Mechanism: assistant-persona
- Context: "(yes, really)", "(this is the part people skip)", "(it isn't)",
  "(spoiler: no)". Chat register inside prose.
- Before: "The form is not the problem (it never was)."
- After: "The form is not the problem."
- Rule: Cut the parenthesis. If it carried a claim, make the claim in a
  sentence.
- Source: model self-report (Fable 5.1, 2026-09-02 revision); provisional until caught in Matt's work

---

### imperative staging

- Added: 2026-09-02  |  Tier: 2
- Mechanism: instruction-tuning
- Context: "Consider the intake form." "Take activation." "Look at the Q2
  numbers." A stage direction before the point.
- Before: "Consider the intake form. Three offices use three versions of it."
- After: "Three offices use three versions of the intake form."
- Rule: Delete the direction; the point is the next sentence.
- Source: model self-report (Fable 5.1, 2026-09-02 revision); provisional until caught in Matt's work
