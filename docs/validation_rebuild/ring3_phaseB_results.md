# Ring 3 Phase B — open-ended hunt with a capable learner: results

**Tip:** origin/main (analysis/ring3_competency/openworld_ppo.py, openworld_ppo_probe.py).
**Bound unchanged:** 0 competencies new to science. The contribution is a **map extension** (known
competencies, new to *our* governing-dynamics map) produced by the capable learner where the
evolutionary learner could not — plus the validated open-ended pipeline.

## The catch first (7th over-claim, caught before any claim)
The first Phase-B run reported "7/8 environments demand competency" — but the memoryless baselines
failed *uniformly* (0.01–0.09), including on trivially reactive tasks. That was a **broken baseline**,
not a result: the memoryless agent did memoryless *rollouts* but the PPO *update* re-ran the GRU with
memory (a train/rollout mismatch), so it never learned, inflating the MCC. Fixed it to a true
feedforward MLP control and re-ran. (Also ran durably under `vps-job` and survived two session
teardowns + is resumable — infra working as intended.)

## Corrected sweep (8 sampled environments, true MLP baseline)
| env | rule (sampled, not designed) | PPO | memoryless | gap | verdict |
|---|---|---|---|---|---|
| 0 | count(t2)≥2 | 0.96 | 0.94 | +0.02 | reactive ✓ |
| 1 | t0-before-t1 AND count(t1)≥3 | 0.72 | 0.10 | +0.62 | **DEMANDING** |
| 2 | count(t1)≥2 AND collect(all) | 0.79 | 0.57 | +0.22 | borderline |
| 3 | t0-before-t1 | 0.86 | 0.49 | +0.38 | **DEMANDING** |
| 4 | collect(t0,t2) | 0.92 | 0.92 | +0.00 | reactive ✓ |
| 5 | count(t1)==3 AND collect(all) | 0.41 | 0.05 | — | PPO under-solved |
| 6 | count(t1)==3 | 0.79 | 0.14 | +0.65 | **DEMANDING** |
| 7 | collect(t1,t2) | 0.80 | 0.92 | −0.12 | reactive ✓ |

The corrected baseline is honest and discriminating: **collection and count≥k are reactively solvable**
(memoryless 0.92–0.94, gap ~0), while **ordering and exact-count genuinely demand internal state.**
3/8 demanding: envs 1, 3, 6.

## Interventional verification (agent-observer protocol — hypothesize, intervene, debunk)
- **env 3 — sequencing / ordering (`t0-before-t1`): GENUINE.** Across seeds the agent reliably touches
  type-0 *before* type-1 (even while wandering through type-2), then goes to goal; memoryless caps at
  0.49 (can't track "prerequisite done") vs PPO 0.86. State-requiring.
- **env 6 — accumulation / counting (`count==3`): GENUINE (imperfect, 0.79).** Decisive supply-invariance
  test: when object supply DOUBLES (5→10), mean type-1 collected stays **2.83 → 3.17** (a greedy
  collect-all would scale ~5→~10). It tracks a count toward 3 and stops, regardless of supply — the
  exact test that debunked sweep-1's greedy mirage, now passed. (On overshoot to 4 it correctly never
  goes to goal, since `==3` is then impossible — it "knows" it failed.)
- **env 1 — combination (sequencing + count≥3): GENUINE.** Touches type-0 before the type-1s and
  collects ≥3 type-1.

## What this adds to the governing-dynamics map
The map was: state-requiring forms {navigation = *commitment*, memory = *storage*} vs reactive forms
{delayed-gratification, regulation}. Phase B's capable learner realizes two faces evolution could not:
- **Accumulation** (increment-a-tally-and-stop-at-threshold) — a **candidate third state-mode**,
  functionally distinct from storage's hold-and-recall. This is the mode we hypothesized and the
  evolutionary learner repeatedly failed to produce; the gradient learner produces it (env 6,
  supply-invariant).
- **Sequencing / ordering** (track a prerequisite / order of events) — a state-requiring competency;
  mechanistically a storage-of-progress (hold "t0 done") applied to the agent's own history.

Both are **known to science** (counting, sequencing) → **0 new-to-science**, the honest bound holds.
But they are **new to our map**, and per the Platonic framing a known competency realized in a
structurally-distinct substrate is a genuine data point that refines the governing law — exactly the
Tier-2 payoff the escalation was for.

## Bottom line
The Phase-A diagnosis is confirmed at the system level: with a **capable gradient learner**, the
open-ended generator (rules sampled, nothing baked in) + MCC filter + interventional verification
**produces genuinely state-demanding competencies the weak evolutionary learner could not** —
sequencing and accumulation — and the self-debunking discipline held (caught the broken baseline; the
supply-invariance test separated counting from greedy). Next: enrich the grammar toward competencies
that are *not* obviously known (conditional/rule-switching, in-episode rule inference) — the genuine
remaining route to off-map, now on a learner strong enough to reach it.
