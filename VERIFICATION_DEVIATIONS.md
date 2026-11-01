# Deviations from VERIFICATION_PLAN.md

Plan SHA-256 8c91c3c8792c9167f491ece52a2aecaa0a3bffe34b9cba7bbc70e52ceed56beb (frozen at commit c223562; uploaded to
the private OSF project 67u4e on 2026-10-28 08:40 UTC, before the first run). The first run finished 2026-10-28 12:20 UTC
(214 of 214, no errors). Both deviations below were found, approved by the owner and logged on 2026-10-28, **before any
classification was computed or read**. Only the following were inspected before logging: whether each tested feature
could be measured in each condition and arm, which views each record contains, and seed counts. No effect, scale or
class was computed.

| # | change | reason | approval |
|---|---|---|---|
| V1 | `census/verify.py` `_set_size` now also sets `census.fingerprint.W` (and `census.filter.W`, `NA`, `NN`, which feed metadata only) to the condition's size. The 4x-area condition is re-run for all 214 phenomena with the fix (`python -m census.verify rerun-area4`), on the same fresh seeds. The census-setting and 3x-length conditions and all sub-programs run at the fingerprint's own size of 64, so they are unaffected and are not re-run. The first-run records are kept unchanged in `census/runs/verify/sims_v1/`. | A bug in the verification harness, not in the census code. `fingerprint.py` holds its own grid size (`W = 64`), which the harness did not update for the 128 x 128 condition. Effects at 4x area: every field feature failed (the radial power spectrum is built for a 64-wide grid, so the fingerprint raised and returned only an error marker; 11 field-view phenomena had no 4x-area values), and the correlation length `g_corrlen` was capped at 32 instead of 64 (a tested feature in 104 phenomena). Agent, network and phase features do not use this constant (the agent box size is passed with the data). Found by checking measurability before classifying: 91 of 1,035 tested features were unmeasurable somewhere. | Owner, 2026-10-28 ("Fix + rerun all") |
| V2 | When the knock-out never produces the flagged view (no seed has it) while the real program does, the effect is undefined under the plan's formula. The census's own frozen knock-out convention applies (`census/knockout.py`: "interaction-driven <=> knock-out view absent OR ..."): such a feature is **interaction-driven**, and **persistent** if at 4x area and at 3x length the real program still has the view and the knock-out still lacks it. Sub-program reproduction is tested as in the plan. | Affects 3 phenomena, all on the transient grid view (`C.transient`, which exists only when a run collapses to a single type after at least 20 recorded frames): census g1 grid-6, census novelty grid-4, quality-diversity 004. Without interaction these runs never collapse, so the knock-out has no transient view. Under the frozen rule they would count as not robust because a quantity could not be computed, not because of any measured result. | Owner, 2026-10-28 ("Use census convention") |

## Logged after Stage 3 ran (found by an independent fidelity audit, 2026-10-28)

| # | change | reason | status |
|---|---|---|---|
| V3 | The Stage 3 script (`census/verify_stage3.py`, committed `85ceedb` before any class was read) added to each card a line "FIRST PASS concluded: <stored verdict>" and the first 1,500 characters of the stored verdict text. The plan's Stage 3 (a) lists only the program, the robust features with their values and the simulation evidence. The stored verdict is the final one: for 57 of the 84 cards it is KNOWN from an earlier strict pass by the same model (claude-opus-5-5) with the same instructions, while that phenomenon's actual first pass was KNOWN-VARIANT (61 of the 84 first passes were KNOWN-VARIANT). The first Stage 3 run therefore largely re-asked the same model with its own earlier verdict in view. | An implementation choice that went beyond the plan, not approved as a deviation and not logged before the run. | Disclosed. Stage 3 is re-run as the plan specifies (`--blind`: no earlier verdict or text on the card; the strict instructions are unchanged, so their sentence about a first pass judging a variant remains), output in `census/runs/verify/stage3_blind/`. Both runs are reported. |

| V4 | Stage 3 (b), the review of the evidence and the cited sources, is carried out by Claude (claude-opus-5-5, through Claude Code research agents with web search and page fetching) instead of the author. Each of the 84 plan-conforming Stage 3 replies is checked: do the load-bearing citations exist, are they attributed correctly, and do they describe this behaviour in this or an equivalent system. Each phenomenon gets a final label under the plan's definitions (documented / closely related / undocumented candidate), and the five known-variants get a deeper search for their stated twist. | The owner cannot do the review (2026-11-01: "i cant review anything - you review"). Limitation: the reviewer is the same model family that wrote the Stage 3 replies and the census verdicts, so the review checks sources and reasoning but is not independent of the model. | Owner, 2026-11-01. Logged before any review result was read. |

**Flag computation (a flaw in the frozen inventory, not changed).** The mismatch pattern searches the stored full text
of each verdict. For strict-pass replies that text ends with a JSON block containing the key `trivial_or_artifact`,
which the pattern `artifact` matches. Ten verdicts are flagged only because of the key (census g1 grid-25, census
novelty grid-4, review 035, 050, 058, quality-diversity 033, 040, 049, 057, 077). Thirteen more are flagged only by
replies that repeat the prompt's own question ("genuine or artifact?"), and the pattern also matches negated uses. The
frozen flags and classes stand. A post-hoc sensitivity analysis (`census/verify_sensitivity.py`, `SENSITIVITY.md`)
removes the key name: 130 flagged (census batch 27), 55 naming card features, of which 25 are artifacts; re-classifying
with all card features the phenomena whose tested set the false flag narrowed changes one class (review 035, artifact
on 3 of 6 features, becomes a single mechanism on all 6, unflagged, so not a Stage 3 case) and no Stage 3 membership.

## Notes (not deviations: the plan's rules applied as written)

- **Flagged view absent in the real run.** For 4 phenomena (census g1 grid-11, grid-17, grid-29 and quality-diversity
  022) the flagged transient view did not appear on any of the 10 fresh seeds at census settings: the representative
  did not collapse to a single type. No feature can then be robust, so the plan classifies them ARTIFACT. They are
  reported with the reason "flagged view did not recur on fresh seeds".
- **Features the representative does not have.** For 5 agent phenomena (census g1 agents-3, agents-8, novelty
  agents-6, agents-8, agents-17) the literature card listed type features (`a_type_entropy`, `a_type_segregation`)
  computed over the class, while the representative has a single agent type. These features are undefined for it and
  are dropped, as the plan states for undefined features; the other card features are tested.
- **Partial seeds.** In one arm (census g1 grid-17, real program, 3x length) the flagged view appeared on 7 of 10
  seeds; medians use those 7, as the plan's drop rule implies.
- **Sub-programs without the view.** 30 sub-program runs in 16 phenomena never produce the flagged view: the phase
  view exists only when the oscillator rule is kept, and the transient view only when the run collapses to one type.
  Such a sub-program does not reproduce the feature (the plan's reproduction test cannot be met without values). No
  sub-program run raised an error. Checked after the 4x-area re-run, before classification: all 214 records carry the
  fix, with no 4x-area errors, and 37 of 1,035 tested features remain unmeasurable somewhere, all in the cases listed
  above.

## Found during the Stage 3 (b) review (2026-11-01; corrections and limitations, not deviations)

- **Feature windows were misdescribed.** The plan and the Stage 3 cards say every feature is "averaged over the late
  part of the run". In the frozen census fingerprint, grid activity, type entropy and interface are late-window means;
  Moran's I, correlation length, largest domain and number of domains come from the last frame; reversals, cyclic
  changes and the Moran trend use the whole run; the activity trend compares early and late activity. The
  computation is the census's own and was not changed; the description is corrected in the paper.
- **Frame spacing depends on run length.** The fingerprint thins each view to about 101 frames, so reversals and
  cyclic changes compare grid frames 10 steps apart in a 1,000-step run and 30 steps apart at three times the length.
  Real program and knock-out are thinned alike, so each comparison is like with like, but for these two features the
  persistence test compares slightly different quantities. Grid activity is computed before thinning, between frames
  two steps apart (the card's "per step" gloss is wrong).
- **Numerically diverged runs.** The census records an `unstable` flag (a field above 1e6 or non-finite stops the
  run) but does not exclude such runs. A post-hoc check (`census/verify_stability.py`, `stability.json`) found 2 of
  the 214 representatives diverge at census settings: census g1 field-2 (an artifact) and quality-diversity 071,
  which stops at step 6; its robust features and its "combination effect" describe that truncation.
- **Other measurement notes from the reviewers.** quality-diversity 065 is measured after the field its agents climb
  has been used up (late-run values are a frozen state); reversals are saturated at 1.0 whenever cells change with
  two types, so that feature only separates a program from a frozen knock-out.
- **Tier-3 for the only review-stage candidate.** quality-diversity 060 was vetted as fixed in VERIFICATION_TIER3.md
  (commits a3d4ad0, 312a973): the decomposition hypothesis was supported (10 pinned flippers reproduce the ordering;
  rules 1+2 the activity decline), the uniform-field control removed the ordering, and an independent literature pass
  found it closely related, a superposition of documented effects. Closed as closely related. The plan expected the
  symmetric control to reproduce no robust feature; it reproduced three weak or saturated ones (number of domains,
  correlation length, reversals) but not the two that carry the effect.
