# EPC Emergence Census — pre-registration v1

**Owner:** Matt Maxwell · **Status:** DRAFT for owner approval — nothing is registered or run until approved
**Code:** branch `emergence-census` of the EPC repository; frozen commit and file hashes in §12

## 1. Plain summary
We list every simple rule-system up to a fixed size, run each one, and ask three questions in a fixed order:
1. **Did anything emerge?** A pattern counts only if it disappears when the rule-system's interactions are removed.
2. **Is it a behaviour we already know?** It is named only if it looks like our reference examples of that behaviour
   AND passes that behaviour's textbook measure. Otherwise it stays unnamed.
3. **If unnamed, is it known to science?** An independent literature check decides. Only a behaviour the literature
   is silent about becomes a novelty candidate, and it must then survive replication and knock-out tests.

Two results are promised whatever happens: (a) a **map** — for each known behaviour, the simplest rule-system that
produces it and every other route to it; (b) an honest count of what was flagged, named, unnamed, and why.
A discovery is possible but not promised.

## 2. What is being searched (frozen)
- **Grammar G1** (`census/grammar_g1.py`): 31 mechanism-level rule templates on four substrates — lattice cells,
  moving agents, network nodes, diffusing fields — and their couplings. Every program has an exact length in bits
  (prefix-free code; coarse parameter values are cheaper than finely tuned ones). Equivalent programs (type
  relabelings) are counted once.
- **Census (exhaustive):** every G1 program of length **≤ 22 bits = 289,798 programs**, shortest first.
- **Second grammar G2** (`census/grammar_g2.py`): the same capabilities built from finer parts (subject, condition,
  action, rate). Exhaustive to **≤ 21 bits = 98,082 programs**. Purpose: test whether "simplest" depends on the grammar.
- **Novelty arm:** quality-diversity search (CVT-MAP-Elites, `census/novelty.py`) over G1 programs that couple at
  least two ingredients, mostly longer than 22 bits; fixed budget of **150,000 evaluations** (or earlier if no new
  niche is filled in 20,000 consecutive evaluations).
- **Simulation** (`census/sim.py`): 64×64 lattice and fields; 400 agents (20×20 box when agents are alone, the 64×64
  lattice when coupled); 200-node random network (mean degree 4); 1,000 steps; 3 seeds per program; fixed random
  initial conditions; seed = hash of the program's bit string (any run can be reproduced exactly).

## 3. Decision rules (frozen)
**FLAG** (`census/filter.py`, `census/knockout.py`). A view of a program (lattice, agents, network, phases, a field,
avalanches) is flagged when (i) the screen says emergent in ≥ 2 of 3 seeds, and (ii) the interaction knock-out says
the pattern depends on interaction: re-running with interaction rules removed (same seeds) either drops the screen's
evidence by ≥ 0.25 and ≥ 2× its seed-to-seed spread, or a textbook order measure is clearly higher with interaction
(≥ 0.1, ≥ 4× spread, ≥ 25% relative). Programs with no interaction rule are never flagged.
- *Screen:* generic emergence score ≥ 0.5, or consensus gain ≥ 0.2, or the model-free complexity detector (not used
  on oscillator-phase views, where it fires at random; network node states use their own three channels).
- *Evidence in the knock-out comparison:* the continuous channels only (generic score, consensus gain). The yes/no
  complexity detector is never evidence: the placebo test (§6) showed it fires at in-between rates on several kinds
  of interaction-free data, and one random "yes" among three seeds is enough to fake a drop.

**NAME** (`census/namer.py`, reference library §12). Each seed of a flagged view is compared with the library: it
gets class *c* only if all 5 nearest library examples are *c*, the nearest is within *c*'s own typical distance (90th
percentile), every other class is ≥ 1.5× farther, *c* has ≥ 3 parameter variants, and the run passes *c*'s textbook
measure. A view's primary name needs ≥ 2 seeds to agree. Names are reported as **behaviours**: classes that share a
textbook measure are one behaviour (voter and majority coarsening → "domain coarsening"); the mechanism is the
program's own rule set (its **route**). A view also lists every other behaviour whose textbook measure it passes.

**Status of a flagged view:** NAMED · KNOWN-BEHAVIOUR-ATYPICAL-LOOK (passes a known behaviour's measure but does not
look like the library — goes to review) · UNNAMED (goes to the literature check). Avalanche views have no
fingerprint; they are flagged by the catalog's SOC detector and reported as "power-law avalanches"; the three
shortest are confirmed with a 5,000-step run at the detector's confirmation tier before the name is used.

## 4. Unnamed results → literature check → candidates (frozen)
1. Unnamed views are grouped per substrate family on robust-z fingerprints by fixed-radius grouping, radius 7
   (`census/triage.py`): in census order, a view founds a group unless a group founder is within the radius; every
   view then joins its nearest founder, so each group is represented by its shortest program. (A fixed radius keeps
   its meaning at any run size; the pilot's Ward clustering did not — the same data copied 40× split into 5× more
   groups — so it was replaced before freezing. At rehearsal size the two give the same number of groups.)
2. Each group, shortest first, gets the literature-check protocol (`census/LITCHECK_PROTOCOL.md`): an independent
   research agent with web search, told that "known" is the expected answer, must return KNOWN / KNOWN-VARIANT /
   LITERATURE-SILENT with verified citations.
3. **KNOWN → catalog expansion:** the behaviour becomes a new library class (≥ 3 parameter variants, its own
   textbook measure). The validation gate and the placebo test (§6) are re-run with the enlarged library on fresh
   seeds and must pass before the new name is used; all stored results are then re-named (no re-simulation needed).
4. **LITERATURE-SILENT → novelty candidate:** replication (10 seeds, two system sizes), rule knock-outs (which
   ingredients are necessary), the EPC interventional debunking battery (null and provably symmetric controls), a
   second independent literature pass, and owner review. Only then any claim.
5. KNOWN-BEHAVIOUR-ATYPICAL-LOOK views are reviewed in the same way as unnamed groups (a new route or a new twist
   on a known behaviour can itself be a finding).

## 5. Analysis plan (frozen)
- **Emergence-complexity table:** for each behaviour, the shortest program in G1 (≤ 22 bits) and in G2 (≤ 21 bits);
  behaviours not reached are reported as "> L".
- **Grammar dependence:** Spearman rank correlation between G1 and G2 shortest lengths over behaviours found in
  both; reported with its value whatever it is.
- **Robustness:** for each behaviour's minimal program, the fraction of all parameter-grid settings of the same
  rule structure that still produce the behaviour.
- **Routes:** every distinct mechanism (substrates + rule templates) that produces each behaviour; ingredients
  common to all routes; rule knock-outs on each minimal program.
- **Census statistics:** flagged, named, atypical and unnamed rates by program length and substrate; number of
  unnamed groups; literature-check outcomes.
- **False-flag audit:** after the run, the placebo test (§6) is repeated on a fresh random sample of 3,000 census
  programs of all lengths; its false-flag count is reported with the results, whatever it is.
- **Website:** each catalog behaviour's page gains its simplest census program and its routes (deploy only with
  owner approval).

## 6. Validation evidence (done before this registration)
The pipeline had to pass a gate whose criteria were fixed before each run; after every rule change the whole gate was
re-run on fresh random seeds. Twice a version passed and then failed the next fresh-seed run by a single case. Both
times the cause was the same kind of fault: a detector that fires at random on some kind of unordered data, which can
make a negligible interaction look decisive. A single pass was therefore not accepted as evidence, and a placebo
test was added that looks for that fault directly and at scale.

Final evidence on the frozen code (§12), every run on seeds never used before:
- **Gate, four separate fresh seed sets** (`census/validation/gate2_1..4`), all four PASS:
  - real examples (textbook-verified library examples) flagged: **474 / 476 (99.6%)**;
  - interacting-but-trivial programs flagged: **0 / 1,243** (rate below about 0.25% with 95% confidence);
  - names the run does not show: **0**; named behaviours that survive the knock-out: **0 / 474**;
  - exact primary name: **414 / 474 (87%)**; the other 60 are sent to review, none gets a wrong name;
  - fingerprint errors: **0** (the gate fails on any; an earlier silent failure of the agent fingerprint is why).
- **Placebo test, two independent runs** (`census/placebo.py`, `census/validation/placebo2_a, placebo2_b`): 2,074
  programs (all of ≤ 14 bits in both grammars plus a random sample of 600 of ≤ 22 bits) have their interaction rules
  removed and are run twice on independent seeds; one run is treated as the "real" program and the other as its
  knock-out. Nothing interacts in either, so any flag is false by construction. **False flags: 0 / 6,642 views**
  (rate below about 0.05% with 95% confidence). Before the two fixes below the same test gave 23 / 2,444.
- Second-grammar check: 8 known behaviours written in G2 flagged and named correctly; 3 G2 negatives not flagged.
- Dress rehearsal, full pipeline end to end on every program of ≤ 14 bits: G1 914 programs — 407 flagged views =
  191 named + 57 review + 159 unnamed in 18 groups; G2 564 programs — 73 flagged views = 23 named + 14 review + 36
  unnamed in 6 groups; 0 errors and 0 fingerprint errors in both. Eight behaviours already appear in the G1 map
  (shortest: synchronization at 9 bits).

Failures on the way to this version (all fixed before freezing, listed so the record is complete):
- Network node states: a generic screen channel fired on pure noise (gate failed 1 / 313). Replaced by three
  network-specific channels; knock-out test made paired over all seeds.
- Oscillator phases: the model-free complexity detector fired on about 40% of uncoupled oscillators (second
  fresh-seed gate failed 1 / 311; placebo 23 / 2,444, all on phase views). Phase views now use the generic score only.
- Knock-out evidence: the same yes/no detector fires at in-between rates on a third of interaction-free field
  programs and on non-moving agents (placebo 1 / 5,690 after the phase fix alone). The knock-out comparison now uses
  continuous evidence only. No real example was lost (474 / 476 flagged, up from 465 / 474).
- Grouping of unnamed views: Ward clustering with a fixed threshold splits the same data into more groups as the
  run grows and needs memory quadratic in the number of views. Replaced by fixed-radius grouping (§4.1).

## 7. Stopping rules
- Census: stops when every program in §2 has run. No early stopping, no extension after seeing results.
- Novelty arm: 150,000 evaluations, or 20,000 consecutive evaluations with no new niche.
- A crashed program is re-run with the same seed; persistent errors are reported, not dropped silently.

## 8. What will be reported regardless of outcome
The map, the statistics of §5, every unnamed group with its literature verdict, every candidate with its vetting
outcome (including those that fail), all deviations, and the error rates of §6.

## 9. Deviations
Any change to code, thresholds, library or procedure after registration is logged with date, reason, and whether it
was made before or after seeing the affected data (`docs/census/DEVIATIONS.md`). Catalog expansion (§4.3) is a
planned procedure, not a deviation, but every expansion is logged with its gate result.

## 10. Known limits (stated up front)
- "Simplest" is relative to a grammar; that is why G2 exists, and instability across grammars is itself a result.
- The reference library covers 9 behaviours plus avalanches. Most unnamed groups are expected to be known science
  that is simply not in the library yet; the literature check and catalog expansion exist for that.
- The namer abstains on roughly one flagged result in eight rather than risk a wrong name.
- Error rates are measured, not zero: with 0 false flags in 6,642 placebo views and 0 in 1,243 trivial programs, up
  to a few hundred false flags among several hundred thousand census views cannot be excluded. A false flag cannot
  become a named behaviour without also passing that behaviour's textbook measure, and cannot become a candidate
  without the 10-seed replication and knock-outs of §4.4.
- Three seeds per program make borderline cases (an effect right at a threshold) unstable from seed to seed; the
  map reports the shortest program for which a behaviour is flagged, and each minimal program is re-checked on 10
  seeds before it is reported as the minimum.
- Runs last 1,000 steps on small worlds: slow or large-scale phenomena can be missed. Gray-Scott-type patterns and
  chemotaxis need longer programs than 22 bits and are reachable only through the novelty arm.
- Forty years of work on small rule systems make a genuinely new behaviour unlikely; the design does not depend on one.

## 11. Compute and visibility
- VPS, 5 worker processes, niced, resumable, as user matthewhmaxwell in tmux. Measured cost: about 12 s per census
  program and 15 s per novelty-arm evaluation (3 seeds each, knock-out included). Estimated 8 days (G1) + 3 days (G2)
  + 5 days (novelty arm) = about 16 days. Storage about 1 GB (summaries only; any run can be regenerated from its seed).
- A daily digest (`census/digest.py`): progress, speed, errors, fingerprint errors, the behaviour map so far, new
  unnamed groups, and — at the top — any candidate.

## 12. Frozen artifacts
- **Code:** EPC repository, branch `emergence-census`, commit `52cee0f38d4672d3eb30b187299adc80a3ae8bd9` (the census
  modules below, plus the catalog's screen and SOC detector they call, are all at this commit).
- **Reference library:** `census/reflib/census_v3/library.json` — 10 classes, 93 parameter variants, seed sets 30–31,
  476 textbook-verified examples.
- **Program lists** (the census is exactly these programs, in this order): `census/programs/programs_g1_22.txt`
  (289,798 programs) and `census/programs/programs_g2_21.txt` (98,082 programs).
- **Environment:** Python 3.12.3, numpy 2.5.0, scipy 1.18.0, networkx 3.6.1, scikit-learn 1.9.0, powerlaw 2.0.0.
- **Validation reports:** `census/validation/gate2_1..4/VALIDATION.md`, `census/validation/placebo2_a, placebo2_b/
  PLACEBO.md`; the failed runs are kept too (`final`, `final3`, `placebo_before`, `placebo_b`).
- **SHA-256 of the frozen files:**

```
513a955f5fa064d6cc617e23a0a7e251a4ce41772e6ed8a2bef3bcbb0c8bbeac  census/coding.py
a2959bcb33bb61440a10cfb95af23c88518ddfc467f209ca3a3bf0f60c44f3b3  census/grammar_g1.py
4d8549cbb16e9740426f6ae2b692fa09050d207c6f32923efd0a5a4d42e99a21  census/grammar_g2.py
729156d9a247c12a56b2678ac312c0c42a25fb1e821f839b82e9b6ccc751bf83  census/sim.py
3a5d885de474f326fe9c847afc99a1dafcb8c94779c863a7bb30be0a7ddc8d80  census/filter.py
50abcacd40202740c68f623baaed8dd57247b3091c6f3c91bed3f6bfa42df6e9  census/fingerprint.py
f53d977794d99530ace7e02bf5fc6db146e6df4c5d1bb152c9d77e88e88ae5ac  census/knockout.py
0989dbb828cdfeefef5f34ed39c21cf6ab9e9f6f4233485c29edc3e68c340b96  census/reflib.py
edd13715ab6b3c384aea4d3a33f2bde7114fa35c32337e2900d5ad7353e00869  census/namer.py
a36b853772173845a0af23193ca8b29b51c80d02cac81ccc9d6158fe184f4a8a  census/runner.py
d44ef0d56313dba67279c6f6924c4a07605bdec327ba5c3bfb2ae471455cf2a7  census/triage.py
a15ef49a6a283c37c5e89638174d662fb3d7b832baad219120df2136ad9060db  census/novelty.py
b487cf6c452ba68876d8a239912b1fa7d3654179ddeb0381d50ffdbb51091ef3  census/digest.py
10853d943e457f0e77e6dc0c837d2eb5575f81550d987967b29e07be095ad2e0  census/placebo.py
cdf2f4742d6f69e6c245a59e2fb7ebec443b40f44dd08de9e50d37dc9272f7ca  census/validate4.py
f924a81c648138f0156ebe21d4cb2fd2f9c3566ac11204ce4dc8c50b3a5a61f9  census/LITCHECK_PROTOCOL.md
f45f9e1bb3eeff24b74fde40d95a8181fc4ac451b1466a3b1385141a2164e029  census/reflib/census_v3/library.json
5a3c3f844d7411f6318bbb45e56edab8e69bd6debce70557edb02c5b89da2962  census/programs/programs_g1_22.txt
3e3694841288d18c0fc6080716a9d758c30ffe301376a01ec9a665d2f51018e1  census/programs/programs_g2_21.txt
```
