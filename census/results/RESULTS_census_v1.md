# EPC Emergence Census — results report (version 4)

**Answers pre-registration OSF `e5hkf`** (project `67u4e`), registered 2026-10-11 before any census program was run.
Pre-registration frozen file `PREREG_census_v1.md`, SHA-256 `7a398029068d04e14d4be3253269799aee906017460d36b6d74900db72ab8427`.
Frozen code: EPC repository branch `emergence-census`, commit `52cee0f38d4672d3eb30b187299adc80a3ae8bd9`
(plus `cca391d` for Deviation 1 and `56a13eb` for Deviation 2). Author: Matt Maxwell.

> **Correction notice.** Version 1 of this report (2026-10-24) stated that every item the pre-registration promised
> had been delivered. That was wrong. It omitted registered analyses (prereg §3 avalanche confirmation, §4.5 review
> of atypical-look views, §5 robustness, minimal-program knock-outs and the false-flag audit, §10 ten-seed recheck),
> it presented the quality-diversity arm as a completed search when its archive was too small to search (§3 below),
> and it overstated three points (grammar independence, the independence of the validation seed sets, and which test
> caught which detector fault). Version 2 corrects all of these. The missing analyses were run afterwards on the
> frozen code; they are reported in §5 with that timing stated. Version 3 adds the corrected quality-diversity arm
> (Deviation 2, post-registration, exploratory; §3). Version 4 (2026-10-26) corrects the first-pass
> literature-check tally in §2(d): one phenomenon was checked twice, so the first pass gave KNOWN 39, not 40, over
> the 60 distinct phenomena. No verdict changed. It also names the second deviation's commit in the header. Earlier
> versions remain in the file history.

---

## 1. Summary

The census enumerated every simple rule-program up to a fixed size in two grammars — **289,798** programs in G1
(≤ 22 bits) and **98,082** in G2 (≤ 21 bits; a finer-grained re-encoding of G1's mechanism vocabulary). Each program
was simulated, screened for emergence, confirmed to depend on interaction by a knock-out test, named against a
reference library, and — where it could not be named — sent to a literature check.

Two results were promised regardless of outcome:

1. **An emergence-complexity map:** for each of 10 known behaviours, the simplest program that produces it in each
   grammar, and every distinct route to it. The simplicity ranking is robust to encoding granularity (Spearman
   ρ = 0.768 between G1 and G2).
2. **A result on novelty, bounded to the exhaustive census:** of the 60 distinct unnamed phenomena (177 fingerprint
   groups), a literature check and a stricter adversarial second pass found **0 literature-silent**; the registered
   review of the 37,860 atypical-look views (73 distinct phenomena) also found **0 literature-silent** (§5.4). Within
   the census bound, **no new-to-science behaviour was found.** A corrected quality-diversity arm, run after the
   census results were seen (Deviation 2, exploratory), explored **19,096** distinct coupled programs — 64 % of them
   longer than the census bound, up to 76 bits — across all six coupling classes, and also found **0
   literature-silent** (§3). That arm is a heuristic search, not an exhaustive one.

---

## 2. Answers to the pre-registration's promised reports (prereg §8)

### (a) The map — emergence-complexity table

Simplest program (length in bits) producing each behaviour, per grammar, with the number of distinct G1 routes.
Every one of the 20 minima (10 behaviours × 2 grammars) was re-checked on 10 seeds and holds at its census length
(§5.2); none had to be replaced.

| Behaviour | Catalog | Simplest G1 | Simplest G2 | G1 routes | G1 programs |
|---|---|---|---|---|---|
| Network synchronization | P9 | 9 | 11 | 10 | 4,002 |
| Lattice phase locking | P9 | 9 | 11 | 24 | 9,356 |
| Network consensus | P18 | 10 | 14 | 32 | 6,198 |
| Network fragmentation | P34 | 11 | 14 | 9 | 2,751 |
| Domain coarsening | P18 | 11 | 14 | 68 | 23,542 |
| Segregation with vacancies | P1 | 11 | 18 | 18 | 5,608 |
| Cyclic waves and spirals | — | 13 | 18 | 12 | 1,309 |
| Flocking | P5 | 14 | 17 | 21 | 19,080 |
| Power-law avalanches | P14 | 15 | 18 | 8 | 151 |
| Spatial game chaos | P27 | 15 | 16 | 17 | 688 |

### (b) Grammar dependence

All 10 behaviours appear in both grammars; the shortest-length rankings correlate **Spearman ρ = 0.768** (n = 10).
G2 is built from the same mechanism vocabulary as G1 at a finer granularity, so this shows the ranking survives a
change of encoding granularity. It does not test independence from the mechanism vocabulary, which would need a
grammar built from different primitives.

### (c) Census statistics

All arms completed with **0 errors and 0 fingerprint errors** (after Deviation 1).

| Arm | Programs / evals | Flagged | Named | Atypical-look (review) | Unnamed (groups) |
|---|---|---|---|---|---|
| G1 (≤ 22 bits) | 289,798 | 141,585 | 39,114 | 31,126 | 75,231 (70) |
| G2 (≤ 21 bits) | 98,082 | 21,172 | 4,122 | 4,024 | 13,131 (61) |
| Quality-diversity arm | 20,200 | 18,030 | 7,457 | 2,710 | 10,357 (46) |

### (d) Every unnamed group with its literature verdict

Each of the 60 distinct phenomena went to a language-model research agent with web search under the frozen
`LITCHECK_PROTOCOL`. First pass (60 distinct phenomena): **KNOWN 39, KNOWN-VARIANT 15, COVERED-BY-TWIN 6, LITERATURE-SILENT 0.** (The
attraction-only agents class was checked twice, in the pilot wave under the label `g1:agents-attract` and again as
its group's representative `g1:agents-1`; both checks returned KNOWN. Versions 1 to 3 counted it twice, giving 40.) The 14
known-variant verdicts that noted the exact rule combination was not found verbatim were re-examined by an
adversarial second pass (stronger model, prompt aimed at behaviour-novelty rather than rule-combination-novelty):
**KNOWN 13, KNOWN-VARIANT 1, LITERATURE-SILENT 0**; 10 of the 14 were judged trivial or artifactual; 9 overturned the
first pass, all toward more-known. Full verdicts with citations accompany this report (`litcheck/`). Limits of this
gate are stated in §4.

### (e) Every candidate with its vetting outcome

**None** — not from the unnamed groups (§2d), the atypical-look review (§5.4), or the corrected quality-diversity
arm (§3). The prereg's novelty-vetting path (§4.4) was therefore not triggered.

### (f) Error rates (prereg §6) and the post-run false-flag audit

Before registration, on the frozen code: the gate passed with negatives drawn fresh on four separate seed sets —
interacting-but-trivial negatives flagged 0/1,243 across the four; real library examples flagged 474/476 (one
library build); named behaviours surviving the knock-out 0/474; fingerprint errors 0. Placebo test
(interaction-free vs. interaction-free; any flag false by construction): 0 false flags / 6,642 views. Three detector
faults were found and fixed before freezing — one by the gate (a screen channel firing on pure noise), two by the
placebo (a complexity detector that was a coin-flip on oscillator phases, and its leak into the knock-out
comparison). **Post-run false-flag audit (§5.1): 0 false flags in 5,972 views over 3,000 fresh census programs**
(rate ≤ 0.05 %, 95 % confidence).

### (g) Deviations

**Deviation 1** (2026-10-21; commit `cca391d`). All 1,206 G2 sandpile (`ADD_GRAIN`) programs crashed because the G2
operator recorded avalanche sizes but not durations. The fix — recording durations as the G1 operator does — is
confined to that G2-only branch, was verified result-identical on 16 already-run programs (0 changed), and the
1,206 were re-run on their original seeds with 0 errors. Made after seeing the affected data.

**Deviation 2** (2026-10-24; commit `56a13eb`). The quality-diversity arm re-run with its archive construction
fixed, approved by the owner after the census results were seen; reported as exploratory in §3.

**Post-run analyses run late.** The analyses in §5 were registered but were run only after the first version of
this report; they use the frozen code and the registered definitions, and the one choice the registration left
open (the 10-seed pass criterion) is stated in §5.2.

---

## 3. The quality-diversity arm was uninformative

The registration included a CVT-MAP-Elites search over longer, coupled programs (budget 150,000 evaluations, or stop
after 20,000 evaluations with no new niche). As implemented, its niche archive was sized from the 200-program
initialization (k = ⌊flagged/2⌋ cells per substrate combination), giving **25 cells in total**; all 25 were filled
almost immediately, so the stopping rule fired mechanically at **20,200** evaluations. Two coupling classes
(cells + agents, cells + fields) had fewer than four flagged programs in the initialization and received no cells,
so programs of those classes could never enter the archive. The arm therefore did not search the region it was
meant to search. It ran as registered; the defect is in the registered design. Its evaluations are counted in the
census statistics, but **no novelty conclusion is drawn from it about programs beyond the exhaustive bound.**

**Corrected re-run (Deviation 2; post-registration, exploratory).** With the owner's approval (2026-10-24), the arm
was re-run with only the archive construction changed (`census/novelty_v2.py`, commit `56a13eb`): a fixed 64 cells
per substrate combination, built from a targeted initialization of up to 2,500 random programs per combination;
search space, evaluation pipeline, quality, descriptors, mutation operators, budget and stop rule as registered.
Because it was decided after the census results were seen, it is reported separately and as exploratory.

| | Registered arm | Corrected arm (Deviation 2) |
|---|---|---|
| Archive cells | 25 | **331** (64 each for CA, AF, CAF, C, N; 11 for CF) |
| Coupling classes with cells | 4 of 6 | **6 of 6** |
| Evaluations | 20,200 | **25,880** (19,096 distinct programs) |
| Programs beyond the 22-bit census bound | — | **64 %** (lengths 13–76 bits, median 26) |
| Errors / fingerprint errors | 0 / 0 | 0 / 0 |

Cells + fields programs rarely produce interaction-driven emergence: 22 of 2,520 random ones were flagged during the
initialization (cells + agents: 128 of 2,000), so that class received 11 cells. The arm stopped under the registered
rule (20,000 evaluations without a new niche) at 25,880 evaluations: the CVT centroids are derived from the
initialization's own programs, which filled 322 of the 331 cells at once, and the search phase filled no new cell
(it improved elites within cells). It therefore explored broadly around its initialization rather than growing its
archive.

**Result.** 20,595 flagged views (5,215 named, 3,268 atypical-look, 12,112 unnamed). The behaviours reached are 8 of
the census's 10 (none new; power-law avalanches and network fragmentation were not reached in this coupled space).
The unnamed and atypical-look views form 81 distinct phenomena; each was literature-checked under the same protocol,
41 with the strict second pass: **KNOWN 78, KNOWN-VARIANT 3, LITERATURE-SILENT 0.** The three known-variants are a
Snowdrift-payoff form of Nowak–May spatial chaos, fixation in the spatial prisoner's dilemma, and a kinetically
slowed Schelling segregation with vacancies — each a documented family with a noted twist.

---

## 4. Limits of the novelty result

- **Language-model literature gate.** Both passes were run by language-model agents with web search. The first
  pass was told, per the frozen protocol, that "known" is the expected answer; the second pass was shown the first
  verdict. Both may lean toward matching a known model.
- **Coarse descriptions.** Each phenomenon reached the gate as summary fingerprint statistics, so a genuinely new
  behaviour sharing a known behaviour's summary statistics would likely be matched to it.
- **One member per group.** Only the shortest member of each fingerprint group was checked; some groups are large
  (up to ~9,000 programs), so within-group heterogeneity is not examined.
- **Size bound.** The exhaustive census covers programs up to 22 bits (G1) / 21 bits (G2). The corrected
  quality-diversity arm (§3) reaches longer coupled programs, but it is a heuristic, post-registration search whose
  archive never grew beyond its initialization; it extends the null only in that weaker, exploratory sense.

---

## 5. Post-run registered analyses

Run on 2026-10-24/25, after version 1 of this report, on the frozen code with the registered definitions
(`census/postrun.py`). Full output: `postrun/POSTRUN.md` and its JSON files.

### 5.1 False-flag audit (prereg §5)

The placebo test repeated on a fresh random sample of **3,000** census programs of all lengths (both grammars):
each program's interaction-free knock-out is run twice on independent seeds, one run treated as "real" and one as its
knock-out, so any flag is false by construction. **0 false flags in 5,972 views; 0 errors.** False-flag rate
≤ 0.05 % (95 %, rule of three).

### 5.2 Minimal programs re-checked on 10 seeds (prereg §10) and rule knock-outs (prereg §5)

Pass criterion (the one choice the registration left open, fixed before the recheck ran): flagged under the census
rule with 10 seeds (≥ 2 seeds emergent and interaction-driven by the paired knock-out), the behaviour named or shown
by the namer, and emergent in at least half of the 10 seeds. **All 20 minimal programs pass at their census length**
(G1 and G2, all 10 behaviours); no fallback was needed. Rule knock-outs: every G1 minimal program consists of a
single rule, so that rule is necessary by construction (removing it leaves an interaction-free program, which the
census never flags).

### 5.3 Robustness (prereg §5) and avalanche confirmation (prereg §3)

Robustness = the fraction of all census programs with the minimal program's rule structure (substrates, type
counts, field count, rule templates; every parameter free) whose run is flagged and shows the behaviour:

| Behaviour | G1 fraction (n) | G2 fraction (n) |
|---|---|---|
| Segregation with vacancies | 1.00 (6/6) | 0.25 (12/48) |
| Flocking | 0.86 (186/216) | 0.51 (756/1,488) |
| Domain coarsening | 0.57 (8/14) | 0.57 (8/14) |
| Network fragmentation | 0.57 (8/14) | 0.57 (8/14) |
| Network consensus | 0.36 (5/14) | 0.29 (4/14) |
| Lattice phase locking | 0.29 (4/14) | 0.15 (30/196) |
| Network synchronization | 0.29 (4/14) | 0.12 (24/196) |
| Spatial game chaos | 0.11 (4/36) | 0.23 (19/84) |
| Power-law avalanches | 0.07 (1/14) | 0.07 (1/14) |
| Cyclic waves and spirals | 0.06 (7/112) | 0.13 (6/48) |

Some behaviours are generic in their rule structure (segregation, flocking); others need specific parameter values
(avalanches need a very slow drive; cyclic spirals need particular thresholds; game chaos needs particular payoffs).
G1 denominators are small because its parameter grid is coarse within 22 bits.

Avalanche confirmation (three shortest power-law-avalanche programs per grammar, 5,000 steps, confirmation tier,
≥ 2 of 3 seeds): **5 of 6 confirmed**, including both grammars' minimal programs (3/3 seeds each). One non-minimal G2
program (`C kC=4 | C[3] ALWAYS() -> ADD_GRAIN() p=0.003`, 18 bits) did not reach the confirmation tier in any seed;
under prereg §3 it does not carry the "power-law avalanches" name. The map's minima are unaffected.

### 5.4 Review of atypical-look views (prereg §4.5)

The 37,860 views that passed a known behaviour's textbook measure but did not look like its library examples were
grouped like the unnamed views (fixed-radius fingerprint grouping per family and shown behaviour) into **155 groups /
73 distinct phenomena**, each literature-checked under the same protocol, with the strict second pass on any
literature-silent or exact-combination-not-found verdict. **KNOWN 69, KNOWN-VARIANT 4, LITERATURE-SILENT 0.** (One
reply was truncated twice by the length cap; the cap was raised and the check re-run, giving KNOWN after the second
pass.)

### 5.5 Catalog expansion (prereg §4.3)

Deferred. Under §4.3 a KNOWN verdict leads to a new library class with its own gate and placebo re-run, after which
stored results are re-named. No reported result depends on it (it would move phenomena from "unnamed" to "named",
not change any verdict), so it is reported here as not yet done.

---

## Reproducibility

Code commit `52cee0f` (+ `cca391d`; post-run analysis scripts `census/postrun.py`, commits `627548a` and later; corrected arm `census/novelty_v2.py`, commit `56a13eb`, using the frozen
modules unchanged). SHA-256 of every module, the reference library (`census_v3`) and both program lists are in the
pre-registration §12. Any run regenerates exactly from the hash of its program's bits.
