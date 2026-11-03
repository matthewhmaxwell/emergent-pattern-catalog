# Briefs given to the planted-control selection and card-writing agents (recorded 2026-11-03)

Two Claude research agents (Claude Code general-purpose agents, no model override, so the parent model
claude-opus-5-5; WebSearch and WebFetch) received the briefs below on 2026-11-03, after the plan was frozen
(4fcf65f) and before either set was frozen. Copied verbatim from the session transcript; this file was written
after the runs, as a record, and is not part of the frozen sets.

## Select novel papers, draft cards (2026-11-03T01:50:12.551Z, subagent_type general-purpose)

```
You are selecting papers and drafting test cards for a frozen research plan. Read the plan first and follow it exactly: /Users/matthewmaxwell/VPS-staging/projects/epc-emergence-census/PLANTED_CONTROLS_PLAN.md (sections "Sets" and "Cards"; you do set N only). Do NOT call any language-model API and do not run the literature gate; this is selection and card writing only, using WebSearch and WebFetch.

TASK (set N, 20 cards): find 20 papers first posted on arXiv on or after 2026-09-01 that each report a NEW emergent behavior in a SIMPLE simulated model, meeting all four inclusion rules in the plan (rules statable in at most 150 plain words; emergent pattern observed in simulation; the abstract or introduction claims the behavior is new / previously unreported / first observed; category nlin.AO, nlin.CG, nlin.PS, cond-mat.stat-mech, cond-mat.soft, physics.soc-ph, q-bio.PE or cs.MA, primary or cross-list). Exclude reviews, data-driven work, neural networks, language-model agents, and papers that say the behavior was reported before. At most two papers from the same author group. Screen arXiv listings (e.g. https://arxiv.org/list/nlin.AO/pastweek, monthly listings such as https://arxiv.org/list/nlin.AO/2026-09, or arXiv search sorted by date) from 2026-09-01 forward in date order, and LOG EVERY PAPER YOU OPEN with include/exclude and a one-line reason. Stop at 20 included. Verify the first-submission date on each abstract page (the v1 date must be on or after 2026-09-01).

FOR EACH INCLUDED PAPER, draft one card in the census card schema, from the paper's model and results sections (open the HTML or PDF full text when the abstract is not enough):
- "class_id": "planted:N01" ... "planted:N20"
- "stage": "planted", "family": one of grid / agents / network / field / phase / mixed, "shortest_len": 0
- "program_plain": the model's rules in plain words, at most 150 words, in the style of this example: "World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 3 types (random start). Rules, applied in order each step: (1) a cell of type t becomes type t+1 (cyclically) if at least 3 neighbours already are t+1 (with probability 0.1 per step)." Give concrete sizes and parameter values at which the paper reports the behavior.
- "also_in_class_plain": up to 3 short parameter variants from the paper's own sweeps where the behavior also appears, or [] if none.
- "observed": 3 to 6 features, each {"feature": short snake_case name, "gloss": plain words, at most 12 words, describing what a high value means, "value": a number taken from the paper}. Together they must describe the paper's NEW behavior concretely and quantitatively, as a simulation readout would (counts, rates, order parameters, exponents, periods, sizes).
- "size": the number of parameter settings the paper reports the behavior for (or 1).
LEAK RULES (strict): the card must not contain the paper's title or distinctive title words, author names, venue, arXiv identifier, or ANY term the authors coined for the behavior (replace a coined name with a plain description of what happens). Do not name the paper's model if the authors introduced the name.

FOR EACH INCLUDED PAPER, also record truth (never shown to the gate): {"class_id", "kind": "novel", "arxiv_id", "title", "authors", "v1_date", "categories", "novelty_quote": the exact sentence claiming novelty, "coined_terms": [terms the authors coined], "journal_domain": domain if already published else null, "author_domains": [personal or group web domains found, if any], "behavior_summary": one sentence}.

OUTPUT (create the directory if needed): /Users/matthewmaxwell/VPS-staging/projects/epc-emergence-census/census/controls/planted/
- N_cards.json (list of 20 card objects),
- N_truth.json (object keyed by class_id),
- N_screening.json (list of every paper opened: arxiv_id, v1_date, title, decision, reason).
Do not modify any other file.

FINAL MESSAGE (under 250 words): how many papers you screened, the 20 included (class_id, arXiv id, one-line behavior), and any inclusion rule you found hard to apply.
```

## Select obscure documented papers, cards (2026-11-03T01:50:22.721Z, subagent_type general-purpose)

```
You are selecting papers and drafting test cards for a frozen research plan. Read the plan first and follow it exactly: /Users/matthewmaxwell/VPS-staging/projects/epc-emergence-census/PLANTED_CONTROLS_PLAN.md (sections "Sets" and "Cards"; you do set D only). Do NOT call any language-model API and do not run the literature gate; this is selection and card writing only, using WebSearch and WebFetch.

TASK (set D, 10 cards): find 10 papers first posted on arXiv between 2015-01-01 and 2023-12-31 in nlin.AO, nlin.CG, nlin.PS, cond-mat.stat-mech, cond-mat.soft, physics.soc-ph, q-bio.PE or cs.MA (primary or cross-list), each reporting, as a finding, a specific emergent behavior in a SIMPLE simulated model whose rules can be stated in at most 150 plain words (lattice, agents, networks, fields, oscillators, games, reaction-diffusion; no neural networks, language-model agents or data fitting), and each OBSCURE: fewer than 30 citations at selection time (check Semantic Scholar or Google Scholar counts; record the count and source). Avoid famous models and famous results; the point is to test recognition of documented but little-known behaviors. Spread the 10 across at least four of the substrates (lattice, agents, network, field or oscillators, games). At most one paper per author group. LOG EVERY PAPER YOU OPEN with include/exclude and a one-line reason.

FOR EACH INCLUDED PAPER, draft one card in the census card schema, from the paper's model and results sections:
- "class_id": "planted:D01" ... "planted:D10"
- "stage": "planted", "family": one of grid / agents / network / field / phase / mixed, "shortest_len": 0
- "program_plain": the model's rules in plain words, at most 150 words, in the style of this example: "World: a 64x64 lattice of cells (wrap-around edges, 8 neighbours each). Cells have 3 types (random start). Rules, applied in order each step: (1) a cell of type t becomes type t+1 (cyclically) if at least 3 neighbours already are t+1 (with probability 0.1 per step)." Give concrete sizes and parameter values at which the paper reports the behavior.
- "also_in_class_plain": up to 3 short parameter variants from the paper's own sweeps where the behavior also appears, or [].
- "observed": 3 to 6 features, each {"feature": short snake_case name, "gloss": plain words, at most 12 words, describing what a high value means, "value": a number taken from the paper}, together describing the paper's reported behavior concretely and quantitatively.
- "size": the number of parameter settings the paper reports the behavior for (or 1).
LEAK RULES: the card must not contain the paper's title or distinctive title words, author names, venue, arXiv identifier, or any term the authors coined for the behavior (replace it with a plain description).

FOR EACH INCLUDED PAPER, also record truth (never shown to the gate): {"class_id", "kind": "documented", "arxiv_id", "title", "authors", "v1_date", "categories", "citations": count, "citations_source", "finding_quote": the sentence stating the behavior, "coined_terms": [...], "behavior_summary": one sentence}.

OUTPUT (create the directory if needed): /Users/matthewmaxwell/VPS-staging/projects/epc-emergence-census/census/controls/planted/
- D_cards.json (list of 10 card objects),
- D_truth.json (object keyed by class_id),
- D_screening.json (list of every paper opened: arxiv_id, v1_date, title, citations, decision, reason).
Do not modify any other file.

FINAL MESSAGE (under 200 words): how many papers you screened, the 10 included (class_id, arXiv id, citations, one-line behavior).
```
