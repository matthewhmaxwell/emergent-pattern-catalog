# Credits

anti-slop consolidates and extends three prior open-source writing skills. Each
contributed something distinct, and each is credited here. The point of this
file is honesty about lineage: the rule lists are not original to anti-slop.

## Sources

### avoid-ai-writing
- Author: Conor Bronsdon
- License: MIT
- Contributed: the tiered vocabulary system (Tier 1 always-replace, Tier 2
  flag-in-clusters, Tier 3 flag-by-density), the context-profile matrix
  (linkedin / blog / technical-blog / investor-email / docs / casual), and the
  severity tiers (P0/P1/P2).
- Note: its vocabulary tiering was itself adapted from `brandonwise/humanizer`
  (github.com/brandonwise/humanizer), which is credited transitively here.

### humanizer
- License: MIT
- Basis: Wikipedia, "Signs of AI writing", maintained by WikiProject AI Cleanup.
  Wikipedia content is licensed CC BY-SA 4.0.
- Contributed: the content-pattern catalog (significance inflation, superficial
  -ing analyses, promotional language, false ranges, copula avoidance,
  negative parallelisms) and the adversarial self-audit step ("what makes this
  obviously AI? now fix it").

### stop-slop
- Author: Hardik Pandya (https://hvpandya.com)
- License: MIT
- Contributed: the false-agency rule (inanimate things doing human verbs), the
  expanded binary-contrast variant table, and the 5-dimension scoring rubric
  (directness / rhythm / trust / authenticity / density).

## What anti-slop adds

- One deduplicated rule library instead of three overlapping ones.
- The ingestion flow: a defined process for memorializing new tells, tagged to
  the generative mechanism that produces them.
- The living corpus: dated, mechanism-tagged tells caught in the wild.
- The protect-list seam: a per-byline overlay so the floor never strips a
  writer's genuine signatures.

## Licensing for distribution

anti-slop's own prose is MIT (see `LICENSE`). Example text in the humanizer /
Wikipedia lineage has been rewritten rather than copied verbatim, but the
underlying material is CC BY-SA 4.0. The Wikipedia-derived portions stay clearly
attributed here; if you adapt them further, keep the attribution and the
share-alike terms. This file serves as the NOTICE for that lineage. It is not
legal advice.

---

## de-slop fork additions (2026-07-19)

This fork (de-slop, by Matt Maxwell) merges into the kjmagnan1s/anti-slop base:

- **slopbeth** by ehmo (MIT): fact-preservation contract, evidence-bound mode,
  noun-swap test, over-editing warning, no-new-formula guard, tool-artifact
  scan, span-level critique.
- **anti-slop-writing** by adenaufal (MIT): the 2026 structural-tell material
  (cadence uniformity, 30-second tests, bimodal burstiness trap, paragraph
  over-fragmentation, colon/semicolon displacement data, current
  high-signal vocabulary).
- **anti-ai-slop-writing** by jalaalrd (MIT): anti-parataxis rule, accuracy and
  honesty rules, plain-text formatting rules.

Personal corpus entries and the client-deliverable profile are original.

---

## 2026-09 revision (v1.1.0)

- `scripts/slop_scan.py`: a standard-library scanner for the countable tells
  (dashes, colon rhythm, cadence bands, openers, paragraph shape, tricolon
  candidates, binary-contrast shapes, tiered vocabulary, announcers, artifacts,
  markdown on plain surfaces) and a compare mode for fact drift, invented
  numbers, dropped qualifiers, edit size, register migration, flattened
  rhythm, and new closers. Original to this fork.
- `references/current-register.md`: the post-delve vocabulary, sincerity
  intensifiers, structural tells that survived the 2026-06 list, the
  displacement table, and the self-edit failure modes. Written from the
  current model's self-report and marked provisional.
- Fourteen corpus entries from the same self-report, filed under a dated
  batch header with the provisional rule.
- Research-paper and email profiles in patterns.md. The scientific-paper
  vocabulary (garner, adept, commendable, noteworthy, plus the words already
  listed) follows Liang et al. 2024, "Mapping the Increase of LLM Usage in
  Scientific Papers".
- Fixes: the scoring arithmetic (five scored dimensions, Fidelity as a gate),
  the LinkedIn em-dash cell that Matt's zero rule already overrode, the
  curly-quote rule for typeset deliverables, and two register idioms in the
  skill's own prose.
