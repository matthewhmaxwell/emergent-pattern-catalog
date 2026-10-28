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
