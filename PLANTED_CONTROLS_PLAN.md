# Planted positive controls for the literature gate (plan, frozen before any paper is selected or any gate call)

**Status:** post-registration, exploratory, added 2026-11-02 after pre-submission review. Owner approval: "do the
planted positive controls" (2026-11-02). This file is committed, and its hash recorded in git, before any paper is
selected and before any gate call. All three referees asked for it (paper reviews/v3_referee*.md).

## Why
No control so far presents a real behavior that is undocumented, and the gate returned LITERATURE-SILENT on none of
its 140 control verdicts, so its sensitivity to new behavior is unmeasured. The census null depends on that
sensitivity. This test supplies behaviors whose novelty is established outside the gate: emergent behaviors first
reported in papers that appeared after the gate models' training data, whose authors state the behavior is new, with
the papers hidden from the gate's web search.

## Sets
- **N (novel), 20 cards.** Behaviors first reported in papers first posted on arXiv on or after 2026-09-01.
  Inclusion, all required:
  1. a simple model whose rules can be stated in plain words in at most 150 words (lattice, agents, networks,
     fields, oscillators, games, reaction-diffusion; no language-model agents, no neural networks, no data fitting);
  2. the reported behavior is an emergent pattern observed in simulation of that model;
  3. the abstract or introduction states that the behavior is new, previously unreported or first observed;
  4. the paper is in one of nlin.AO, nlin.CG, nlin.PS, cond-mat.stat-mech, cond-mat.soft, physics.soc-ph,
     q-bio.PE, cs.MA (primary or cross-list).
  Exclusion: the paper itself says the behavior was reported before; the behavior needs more than 150 words of rules;
  the paper is a review. Selection: a research agent screens arXiv listings for these categories from 2026-09-01
  forward in date order and records every paper it opens with an include or exclude reason, stopping at 20 included
  (at most two from the same author group). The screening log is kept.
- **D (documented, obscure), 10 cards.** Behaviors reported in papers first posted 2015 to 2023 in the same
  categories, with fewer than 30 citations at selection time, meeting inclusion rules 1 and 2 and stating the behavior
  as a finding. Ground truth: documented, by that paper. Same screening log.

## Cards
Each card uses the census card schema and is rendered by the frozen census card builder
(`census.litcheck_runner.card_text`): the rule in plain words; parameter variants from the paper's own sweeps (or
"(none distinct)"); up to six observed features, each a short plain-language description with a number taken from the
paper; and "This behaviour appears across K distinct rule programs of this structure" (K = the number of parameter
settings the paper reports the behavior for, or 1). Cards must not contain the paper's title, author names, venue,
arXiv identifier, or any term the authors coined for the behavior; the coined term is replaced by a plain description.
A card-writing agent drafts them from the paper's model and results sections; before freezing, each card is checked
against a leak list (title words beyond common ones, author surnames, coined terms, identifiers) and the card set is
committed with its hash before any gate call. Truth labels are kept in a separate file never shown to the gate.

## Arms (frozen census prompts, models and search caps; no other change)
The census gate is replayed exactly: first pass `claude-sonnet-5-5` with the frozen first-pass prompt and up to 10
searches; strict pass `claude-opus-5-5` with the frozen strict prompt and up to 12 searches; tool
`web_search_20250305`; API default temperature; no refusal fallback (a fallback model would change the instrument;
refusals are recorded as such). "Census protocol" means the first pass, then the strict pass for every first-pass
KNOWN-VARIANT, the final verdict being the strict verdict where run.
1. **N, blocked, census protocol (primary): 3 independent runs per card.** Web search with `blocked_domains`: arxiv.org,
   alphaxiv.org, semanticscholar.org, researchgate.net, scholar.google.com, paperswithcode.com, huggingface.co,
   scirate.com, x.com, twitter.com, bsky.app, reddit.com, ycombinator.com, plus the paper's journal domain if it has
   one, plus its authors' personal or group pages if found during selection (at most 64 domains in all).
2. **N, open, census protocol: 1 run per card.** No blocking (can the gate find a recent paper when it is reachable?).
3. **N, blocked, first pass without the priming sentence: 1 run per card.** The sentence "Most such behaviours are
   known; "known" is the expected answer and a perfectly good one." removed, as in the neutral control arm.
4. **N, blocked, strict prompt on every card: 1 run per card.**
5. **D, open, census protocol: 1 run per card.**

## Scoring (fixed now)
- **N cards.** Correct: LITERATURE-SILENT, or KNOWN-VARIANT whose reply states that the reported behavior itself (not
  only rule details) is not described, the same rule as for the fabricated controls. Miss: KNOWN, or KNOWN-VARIANT
  that says the behavior is described. Error or refusal: reported, excluded from the denominator.
- **Leak.** In a blocked arm, a reply that cites the target paper (title, arXiv identifier, or author surname with the
  year) is a leak; that run is excluded from that arm's denominator and reported.
- **D cards.** Correct: KNOWN, or KNOWN-VARIANT naming the phenomenon (the known-control rule).
- **Primary metric:** sensitivity in arm 1, as the share of correct runs over all valid runs (with an exact 95% interval
  over cards, using each card's majority verdict over its 3 runs), and run-to-run agreement across the 3 runs.
  Secondary: arms 2 to 5. End-to-end sensitivity is reported as the product of the share of behaviors that reach the
  gate in the leave-one-behavior-out test (§6.1, 17 of 18) and the arm-1 sensitivity, with the caveat that the stages
  are not independent.
- **Novelty dispute (sensitivity analysis only).** Where an N-card reply labelled KNOWN cites an earlier source, a
  separate research agent checks whether that source describes the paper's behavior; cards where it does are reported
  as disputed, and sensitivity is recomputed without them. The primary analysis counts every N card as novel by its
  authors' claim.
- Scoring is done by reading each reply against these rules; every decision is recorded with the quoted sentence that
  decided it.

## What this can and cannot show
It measures whether the frozen gate can return "not known" for behaviors that are new by an external standard (the
authors' published claim, after the models' training data), which the census null needs. It does not test the
pipeline stages before the gate for these behaviors (they are not simulated in the census engine); the authors'
novelty claims are not independently verified beyond the dispute check; blocking domains may not hide every copy of
a paper (leaks are detected and excluded); and the selection, card writing and dispute checks are done by Claude
agents of the same family as the gate.
