# Simulation verification of the census literature verdicts (plan, frozen before any verification run)

**Status:** post-registration, added 2026-10-27 after peer review. It carries out, in a simulation-based form, the
registered but deferred step of prereg §4.3 (a KNOWN verdict should be checked against the documented model before it
counts). Owner approval: 2026-10-27 ("do it now"; scope: all phenomena). This file is committed and its SHA-256 is
uploaded to the private OSF project 67u4e before the first verification run.

## Why
Controls (census/controls/SCORING.md) showed the language-model literature gate recognizes model families but labels
behaviors that contradict a documented system as KNOWN, explaining the contradiction as a measurement artifact
(fabricated behaviors labelled KNOWN in 18 of 20). About 130 census verdicts note a mismatch or artifact and still say
KNOWN. Whether those observations are artifacts or real behaviors can be decided by simulation, not by a reader.

## Scope
Every phenomenon that reached a literature verdict: the 60 distinct unnamed census phenomena (both grammars, including
the 6 grammar-2 twins, each verified on its own representative), the 73 atypical-look phenomena of the post-run review,
and the 81 phenomena of the corrected quality-diversity arm. 214 phenomena in all. Each is tested on its representative
program (the shortest member of its group), decoded from its stored bit string with its own grammar, and measured on the
view in which it was flagged (a grid program has a whole-run and a transient view with the same feature names). The
representative and its view come from the census class list for the census batch. For the other two batches they come
from rebuilding each batch's grouping step from its stored labels, in the original code's order and without any model
call; every rebuilt representative matched the program and behaviours stored with its verdict (154 of 154).

## Stage 1: what is the discrepancy (no language model)
- **Card features:** the standout fingerprint features on the phenomenon's literature card. The post-run review and
  quality-diversity cards were not stored with their verdicts; for those two batches the card features are recovered
  with the card builder's own rule (the six largest-magnitude features of the representative's stored fingerprint).
- **Flagged:** the verdict text (first pass, plus strict pass where run) matches the mismatch pattern used in the
  controls. Verdict text is the full text and notes; census first-pass verdicts were stored without full text, so for
  them it is the stored phenomenon description, model family and standard name. Citations and search strings are left
  out, so that a cited title cannot trip the pattern. The pattern is the controls' pattern with the British spelling
  artefact added: contradict, inconsistent, does not match, doesn't match, not expected, unexpected, implausib(le),
  impossible, artifact, artefact, anomal(y/ous), cannot produce, should not, would not produce, at odds, mismatch,
  opposite (case-insensitive).
- **Discrepant features:** card features named verbatim in the verdict text of a flagged phenomenon; if none is named,
  all card features. For unflagged phenomena, all card features.
- **Mechanics check before freezing:** three phenomena (one-rule, one-rule grammar 2, two-rule) were simulated once to
  confirm that every condition runs and returns the view's features on 10 seeds; only timings, errors and array shapes
  were read, and those outputs are not part of the results.

## Stage 2: simulation tests (frozen census code; only system size and run length are changed)
For the representative program R and its interaction knock-out K (the frozen census knock-out), 10 seeds each, under
three conditions: census settings (64 x 64 lattice and fields, 400 agents, 200 nodes, 1,000 steps); 4x area (128 x 128,
agents x4 with the agent box doubled in side, nodes x4 at the same mean degree); 3x length (3,000 steps). For a program
with two or more rules, every proper non-empty sub-program (rules removed, parameters kept) is also run at census
settings on 10 seeds. Features are the frozen census fingerprint features of the relevant view, averaged over the late
window as in the census.

Definitions, per feature f (medians over seeds; s = 1.4826 x pooled median absolute deviation of the two seed sets
compared, floored at max(0.02 x (|m1| + |m2|), 0.001)):
- **Interaction effect** e = median(R) - median(K). **Interaction-driven** if |e| >= 3 s at census settings.
- **Persistent** if, at 4x area and at 3x length (each against its own knock-out), e keeps its sign and |e| >= 0.5 x
  the census-setting |e|.
- **Robust** = interaction-driven and persistent.
- **Replicates the card** if |median(R) - card value| <= max(3 s, 0.25 x |card value|) (reported, not used to classify).
- A sub-program P **reproduces** f if |median(R) - median(P)| <= 3 s (s pooled over R and P).

## Classification (frozen)
1. **ARTIFACT**: no discrepant feature is robust.
2. **DECOMPOSES**: R has two or more rules, at least one discrepant feature is robust, and some proper sub-program
   reproduces every robust discrepant feature. The behavior is that of a smaller set of documented mechanisms.
3. **COMBINATION EFFECT**: R has two or more rules, at least one discrepant feature is robust, and no proper
   sub-program reproduces all of them. -> Stage 3.
4. **SINGLE MECHANISM**: R has one rule (one textbook mechanism) and at least one discrepant feature is robust. ->
   Stage 3 if the phenomenon was flagged (the gate dismissed a robust, interaction-driven feature as an artifact);
   otherwise reported as a documented single mechanism.
5. **NOT RUNNABLE**: the representative cannot be decoded or simulated; reported with the reason.

## Stage 3: candidates
For each phenomenon sent to Stage 3: (a) the strict literature prompt (claude-opus-5-5, web search) is given the
program, the robust features with their values, and the simulation evidence (10 seeds, 4x area, 3x length, knock-out,
sub-programs), and asked whether this behavior, in this or an equivalent system, is documented; (b) the author reviews
the evidence and the cited sources. Final labels: **DOCUMENTED** (a source describes the behavior), **CLOSELY RELATED**
(a source describes a near relative; difference stated), **UNDOCUMENTED CANDIDATE** (no source found after (a) and
(b); reported with all evidence as a candidate, never as a discovery). Any undocumented candidate then gets the
registered Tier-3 checks: null and provably symmetric controls and a second independent literature pass.

## Reporting
Counts per class for each batch, with every phenomenon's numbers in a table; the plan's thresholds are not changed
after any result is seen. If a rule above turns out to be unworkable when first run, the change is logged as a
deviation from this plan with its reason before any classification is read.
