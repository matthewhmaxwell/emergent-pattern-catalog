# Ring 3 — competency-search comprehensiveness audit (and search compass)

**Purpose.** Make the claim "our competency search is comprehensive" *auditable* (not "we tried a lot"),
AND use the coverage map as a **compass**: the principled-empty cells are the next targets. Comprehensive
is scoped — *comprehensive within this taxonomy, at toy substrate scale*. Exhaustive is impossible.

## The generating taxonomy (axes that partition the competency space)
- **A. Level** — object (world-facing) / meta (rule-facing) / multi-agent (other-facing) / composition.
- **B. State** — reactive / state-requiring.
- **C. What the state holds** — commitment (held intention) · storage (held value, incl progress) ·
  accumulation (running tally→threshold) · belief (distribution over a latent rule) · externalized
  (state written to the WORLD) · none (reactive).
- **D. Multi-agent coordination MEDIUM** — channel · observation · environment(stigmergy) · [reputation?].
- **E. Motive** — cooperative · competitive/mixed-motive.

## Coverage map (16 verified + 2 negatives → cells)
| competency | level | state-content | medium/motive | status |
|---|---|---|---|---|
| navigation | object | commitment | — | ✓ |
| memory | object | storage | — | ✓ |
| counting | object | accumulation | — | ✓ |
| sequencing | object | storage(progress) | — | ✓ |
| delayed-grat | object | none(reactive) | — | ✓ |
| regulation | object | none(reactive) | — | ✓ |
| rule-inference | meta | belief(1-shot) | — | ✓ |
| rule-tracking | meta | belief(changing) | — | ✓ |
| conditional-selection | meta | storage(cue)→switch | — | ✓ |
| communication | multi | — | channel/coop | ✓ |
| role-division | multi | — | observation/coop | ✓ |
| contest-resolution | multi | — | observation/competitive | ✓ |
| fair-alternation | multi | none(reactive) | (none)/competitive | ✓ |
| count-communication | composition | accumulation×channel | — | ✓ |
| stigmergy | multi | externalized | environment/coop | ✓ |
| tool-use | object | externalized(construct) | self-instrumental | ✓ |
| anticipation | object | (prediction) | — | ✗ NEGATIVE (reactive shortcut dominates) |
| deception | multi | — | channel/competitive | ✗ NEGATIVE (self-defeating in minimal signaling) |

## Saturation curve (new competency-KIND per probe, in discovery order)
commitment→**NEW**, storage→**NEW**, accumulation→**NEW**, reactive→**NEW**, belief/meta→**NEW**,
multi-channel→**NEW**, multi-observation→**NEW**, composition→**NEW**, environment/externalized→**NEW**
... then: contest=observation-variant, fair-alternation=reactive-variant, tool-use=externalized-variant,
anticipation=negative, deception=negative → **the last 5 diverse probes added ~0 new KINDS.** The curve is
flattening (like a species-accumulation curve), but stigmergy (#15) DID add a kind recently → not yet
provably flat.

## Taxonomy red-team (against the cross-disciplinary competency vocabulary)
Science's competency kinds vs our coverage. COVERED: perception-via-marks, memory, counting/numerosity,
sequencing, decision/stopping, control, meta-learning, communication, cooperation, competition, tool-use,
composition. **NAMED KINDS WE HAVE NOT PROBED (the genuine gaps):**
1. **Abstraction / systematic generalization** (meta) — apply a learned RELATION to HELD-OUT instances
   (relational reasoning), not memorize instances. Cleanly demandable via a train/test split. **UNTOUCHED.**
2. **Reputation / indirect reciprocity** (multi) — coordinate via a tracked HISTORY of others' past
   behavior = a **4th coordination MEDIUM (social memory)**, completing axis D. Also the mechanism that
   *stabilizes* honest signaling — directly addresses the deception negative. Needs iteration. **UNTOUCHED.**
3. **Imitation / social learning** (multi) — acquire behavior by observing a demonstrator. Demandable
   (demonstrator-present vs absent). **UNTOUCHED.**
4. Harder/richer (likely need scale): theory-of-mind, causal reasoning, multi-step planning, teaching,
   coalition formation (3+ agents). Flagged as principled-empty-at-toy-scale.

## Comprehensiveness scorecard (the 5 signals)
1. Taxonomy + coverage map: **built** (this doc) — and it exposes 3 untouched demandable cells.
2. Saturation curve: **flattening, not flat** — last 5 probes ~0 new kinds, but a recent probe (stigmergy)
   added one, so saturation is not yet established.
3. Substrate-expressivity bound: **argued** (composition only under forced pressure; precise recurrent
   computation needs gradient; deception self-defeats) — but not yet the formal write-up.
4. Red-team: **partial** — 2 forced-novelty probes + 2 negatives done; the *taxonomy* red-team (this doc)
   now names 3 untested cells, so red-team coverage is incomplete until those are probed.
5. Unbiased generator: **yes** for the object grammar (openworld samples rules without baking the answer);
   not yet demonstrated for the meta/multi-agent probes (hand-built).

**Verdict: NOT YET comprehensive — three named, demandable cells are genuinely untouched.** This is the
compass working: it converts "diminishing returns" into a *specific* finish line. Probing cells 1–3 will
either (a) add new kinds → the search was NOT saturated (keep going), or (b) recur to known / negative →
saturation CONFIRMED and we can declare comprehensive-at-toy-scale with all 5 signals satisfied.

## Guided next search (in priority order)
1. **Abstraction / systematic generalization** — train a relational rule on a colour subset, test on
   held-out colours; the train→test generalization gap is the competency. Distinct, clean, untouched.
2. **Reputation / indirect reciprocity** — iterated game where an agent tracks a partner's past honesty
   and conditions cooperation on it; completes coordination-medium axis D AND tests whether reputation
   stabilizes the deception that minimal signaling could not.
3. **Imitation / social learning** — learner reaches a hidden goal only by copying a demonstrator;
   demonstrator-ablation is the diagnostic.
Each is logged to the registry. After these three, the saturation curve + red-team + expressivity bound +
unbiased-generator signals would jointly license a *comprehensive-at-toy-scale* declaration.
