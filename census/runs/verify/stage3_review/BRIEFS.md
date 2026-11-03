# Briefs given to the Stage 3 (b) review agents (deviation V4), 2026-11-01

Five Claude research agents (claude-opus-5-5, Claude Code, web search and page fetching, read access to the
repository) received the following. Groups A to D got the common brief with their own input file and a one-line
description of their model families; group E got the deep-review brief. The Tier-3 literature pass brief for qd:060
is at the end.

## Common brief (groups A to D)

You are doing the "Stage 3 (b)" source review for a scientific verification study. The standard is "must survive
hostile peer review". Work carefully and record evidence.

CONTEXT. A census of simple simulated rule-programs (lattice cells, agents, networks, fields) flagged behaviours. A
simulation check confirmed, for each phenomenon below, that certain measured features are real: caused by interaction
(they differ from an interaction-free knock-out by at least 3 robust spreads) and persistent at 4x area and 3x run
length; for multi-rule programs, no smaller sub-program reproduces all of them ("combination effect"). Then a
language model (claude-opus-5-5 with web search, "Stage 3 (a)") judged whether each behaviour is documented: KNOWN,
KNOWN-VARIANT or LITERATURE-SILENT, with citations. Your job is Stage 3 (b): review that judgement against the
sources and assign the final label.

Known failure mode you must guard against: in controls, the same kind of checker recognized the model family and
labelled behaviours KNOWN even when the reported behaviour contradicted the documented system (18 of 20 fabricated
behaviours labelled KNOWN). So check that a source describes THIS behaviour (the robust features on the card, at the
measured values, and for a combination effect the way the full program differs from its sub-programs), not merely
the model family or each rule on its own.

FINAL LABELS (the study's frozen definitions):
- DOCUMENTED: a source describes the behaviour in this or an equivalent system. For a combination effect, a source
  must describe the combined behaviour (the coupling of these mechanisms producing this effect), or the effect must
  follow directly and explicitly from a documented result about that coupling. Sources that each document one rule
  alone are not enough.
- CLOSELY RELATED: a source describes a near relative (same mechanisms or coupling in a nearby setting, or the same
  general effect in a related model), and you can state the difference.
- UNDOCUMENTED CANDIDATE: after at least 3 targeted searches for the behaviour itself, no source describes it or a
  near relative. Record the searches. (This label is a candidate for further checks, never a discovery.)
A trivial outcome that follows from textbook results (e.g. an irreversible spreading rule filling the lattice,
fixation in an absorbing state) counts as DOCUMENTED if a source states that result for this or an equivalent system.

PROCEDURE, for each phenomenon in your input file:
1. Read the card (rule in plain words + simulation evidence) and the Stage 3 reply (verdict, standard name, notes,
   citations, full text).
2. Pick the 1 to 3 load-bearing citations. Verify each with WebSearch / WebFetch: does it exist (arXiv abstract, DOI
   or publisher page, Google Scholar-style listing; Wikipedia is acceptable only for existence and basic facts of
   classic models), is the attribution correct (authors, year, venue), and does it support the specific claim it is
   cited for? Reuse sources you have already verified for other phenomena in your group; do not re-fetch them.
3. If the cited sources do not establish the label, search yourself for a source that does.
4. Assign the final label and say why in one or two sentences.

INPUT: inputs/<group>.md. OUTPUT: outputs/<group>.json, one object per phenomenon with keys pid, stage3_verdict,
final_label, agrees_with_stage3, behaviour_reviewed, basis, difference, citations_checked (citation, url, exists,
attribution_correct, supports_claim, note), citation_problems, searches_run, confidence.

Group descriptions: A, spatial evolutionary games (Nowak-May prisoner's dilemma, snowdrift) with voter noise,
conformity, mobility or rewiring (25); B, cyclic cellular automata and rock-paper-scissors lattice models (12);
C, Schelling segregation, voter and majority-vote or Potts coarsening, adaptive voter networks, bootstrap percolation
(26); D, self-propelled particles, chemotaxis and autochemotaxis, Kuramoto oscillators (16).

## Deep-review brief (group E, the five KNOWN-VARIANT replies)

Same context. Labels as above, with UNDOCUMENTED CANDIDATE requiring at least 5 targeted searches for the specific
twist plus the citation trail of the closest papers. "Be adversarial in both directions: the Stage 3 model may have
under-searched (so the twist may in fact be documented), and it is known to downplay or explain away unfamiliar
behaviour (so a twist it called 'minor' may be a genuinely undocumented phenomenology). Judge the behaviour (the
measured features and how the full program differs from its sub-programs), not the exact rule combination."

## Tier-3 second literature pass (qd:060)

A fresh agent received only `census/runs/verify/tier3/qd060_litpass_input.md` (the plan-conforming card, the control
results and short notes on the code; it did not see the Stage 3 replies or the review) and was asked whether the
behaviour, or a near relative, is documented, with at least 6 targeted searches, every relied-on source opened and
verified, and one label (DOCUMENTED / CLOSELY RELATED / UNDOCUMENTED CANDIDATE), plus a judgement of whether it is a
genuine emergent pattern or a trivial superposition of known effects.
