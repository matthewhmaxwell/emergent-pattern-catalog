# Ring 3 — the commonality hunt: an unknown REGULARITY across emergent competencies

**Reframe (user directive, 2026-06-30):** keep searching for a genuinely-new competency (Track 1: scale),
AND — since the toy regime is exhausted for *individual* competencies — treat the verified competencies as
DATA and hunt for a **cross-cutting regularity / law** that the literature does not already name (Track 2).
The novelty would be a *principle of emergence*, not a new behavior. Same discipline: any candidate law
must survive the literature gate (most regularities are already named).

## The corpus (verified competencies, A–D)
| # | competency | level | state? | what state holds | diagnostic perturbation (collapses it) | learner that reached it | lit name |
|---|---|---|---|---|---|---|---|
| 1 | navigation | object | yes | commitment (held detour) | away-from-goal barrier | directed search | path-planning |
| 2 | memory | object | yes | storage (a held value) | cue/feedback ablation → chance | ES (robust) | working memory |
| 3 | accumulation/counting | object | yes | a running tally + threshold | supply-doubling (counter is supply-invariant) | **PPO only** (ES/SNES failed) | counting |
| 4 | sequencing | object | yes | storage-of-progress | order violation | PPO | ordered execution |
| 5 | delayed-gratification | object | **no** | — (reactive) | — | reactive | optimal stopping |
| 6 | regulation | object | **no** | — (reactive) | — | reactive | homeostasis/control |
| 7 | rule inference (1-shot) | meta | yes | belief over a latent rule | feedback ablation → chance | PPO | meta-learning / RL² |
| 8 | rule tracking (continual) | meta | yes | belief over a *changing* latent | feedback ablation | PPO | contextual/meta-RL |
| 9 | conditional selection (cued) | meta | yes | held cue → policy switch | cue removal | PPO | contextual policy |
| 10 | communication | multi | n/a | a shared symbol convention | channel scramble → chance | joint PPO (≥2 agents) | Lewis signaling |
| 11 | role-division | multi | n/a | partner's state (observed) | blind-partner → collapse | joint PPO | role allocation |
| 12 | contest-resolution | multi | n/a | relative-position asymmetry | blind-partner | joint PPO | Bourgeois / uncorrelated asymmetry |
| 13 | fair-alternation | multi | **no** | — (reactive on last action) | NEITHER channel nor memory | joint PPO (iterated) | turn-taking equilibrium |

## Candidate cross-cutting regularities (each literature-gated)
1. **Minimal-mechanism law.** The learner converges to the *minimal* sufficient mechanism (e.g.,
   fair-alternation solved by a reactive last-action flip, not the channel/memory it was "supposed" to
   need). → **KNOWN** (simplicity / Occam bias of optimization; NN simplicity bias, MDL).
2. **State-required ⇔ hidden-latent.** A competency needs internal state IFF the goal-relevant variable is
   NOT in the instantaneous observation (counting's tally, memory's cue, meta's latent rule are hidden;
   alternation's "whose turn" is recoverable from last-action → reactive). → **KNOWN** (POMDP belief-state
   theory: state needed iff partially observable).
3. **Diagnostic-perturbation uniqueness.** Every genuine competency has a UNIQUE perturbation that collapses
   it to baseline, and that perturbation *is* the operational definition of what the competency is FOR
   (its goal-content). → partly Levin/causal-intervention (KNOWN method), but the systematic
   *form ↔ diagnostic-perturbation* map as a classification is **under-named** — worth gating carefully.
4. **What-state-holds ladder.** State content climbs in abstraction: commitment (intention) → storage
   (value) → accumulation (count) → belief (distribution over a latent rule). → memory taxonomies exist;
   a *unified cross-level* ladder spanning object→meta may be a fresh organizing axis (**candidate**).
5. **★ Substrate-reachability ordering (the most promising).** The competency *form* (functional signature)
   is substrate-invariant — counting "stop at k regardless of supply" is the same target however reached —
   but *which learner can reach it* is structured: counting was UNREACHABLE by FSM-mutation, OpenAI-ES, AND
   SNES, yet REACHABLE by gradient PPO. Across our corpus the state-requiring competencies show a
   **reachability hierarchy over optimizers** (gradient ≻ natural-ES ≻ vanilla-ES ≻ FSM-mutation), and the
   *threshold* at which a competency becomes reachable appears to track its state-content abstraction
   (ladder #4). A quantified, empirical "competency-reachability spectrum" — *the same Platonic form, a
   lawful ordering of which optimizers can ingress it* — is the candidate least obviously covered by an
   existing named law. (Learning-theory has hardness results, but an empirical reachability map of
   *emergent competencies across optimizers*, tied to a state-content ladder, is not a thing I can name.)

## Verdict + plan
Most candidate regularities are KNOWN (simplicity bias, POMDP, causal intervention). The two that may be
genuinely under-named are **#4 (what-state-holds ladder)** and especially **#5 (substrate-reachability
ordering of a substrate-invariant form)** — and they may be the same law viewed two ways: *the more
abstract the state-content a competency requires, the stronger the optimizer needed to ingress its (fixed,
Platonic) form.*

**Next concrete step (Track 2):** rigorously test #5. Take 3–4 competencies spanning the state-content
ladder (memory → counting → sequencing → rule-inference) and run EACH through the SAME battery of learners
{FSM-mutation, OpenAI-ES, SNES, gradient PPO} on matched tasks, measuring the reachability threshold per
competency. If the reachability order is lawful and tracks the state-content ladder, we have a candidate
regularity to literature-gate hard. This is cheap (reuses existing harnesses) and directly attacks the
unknown-commonality target. Track 1 (scale for a new competency) proceeds in parallel if resourced.
