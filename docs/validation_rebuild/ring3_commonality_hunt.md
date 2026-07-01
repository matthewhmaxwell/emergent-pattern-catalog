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

6. **★ Sufficiency-pressure / minimal→compositional transition (NEW, from the composition probe #14).**
   The "minimal-mechanism" law (#1) is sharper than "agents find the simplest rule": **agents find the
   minimal SUFFICIENT mechanism, and compositional structure emerges ONLY when no single minimal mechanism
   is sufficient for reward.** Demonstrated *causally + reversibly* with one dial: in the count→communicate
   probe, peaked (Binomial) counts make a 1-bit code sufficient → the agent uses a single non-compositional
   symbol (0.63, 1/3 slots informative); uniform counts make 1 bit insufficient → a fully **compositional**
   3-slot code emerges (0.84, 3/3 slots). So compositionality is not a property of the agent or the task
   alone but a function of **sufficiency-pressure**, and it is *controllable*. Open, measurable question:
   **is the minimal→compositional transition a sharp PHASE TRANSITION in sufficiency-pressure, with a
   lawful threshold, measurable uniformly across the corpus (not just language)?** Compositionality-
   emergence *conditions* are studied in emergent-comms, but a phase-transition law parameterized by a
   sufficiency-pressure dial, as a general property of emergent competency, is a candidate I cannot name.

## Verdict + plan
Most candidate regularities are KNOWN (simplicity bias, POMDP, causal intervention). The candidates that
may be genuinely under-named: **#4 (what-state-holds ladder)**, **#5 (substrate-reachability ordering)**,
and **#6 (sufficiency-pressure / compositionality phase-transition)**. #4 and #5 may be one law (*more
abstract state-content ⇒ stronger optimizer needed to ingress the fixed Platonic form*); #6 is the newest
and the only one demonstrated *causally and reversibly* with a controllable dial — making it the most
directly testable. Two concrete Track-2 experiments now queued: (a) the **reachability battery** (#5), and
(b) the **compositionality-transition sweep** (#6: vary count distribution peaked→uniform, locate the
threshold, test sharpness/lawfulness). Both literature-gated hard before any claim.

7. **Where the state lives — internalized vs EXTERNALIZED (from stigmergy #15).** The what-state-holds
   ladder (#4) implicitly assumed state is *inside the agent*. Stigmergy is a clean counter-case: the
   coordination "memory" lives in the **environment** (the shared trail), not the agents — verified
   inter-agent (own-trail-only 5.7 < combined 8.1; no-trail 1.4). So the ladder needs a prior axis: *is
   the task-relevant state held in the agent or in the world?* Externalized state (stigmergy, niche
   construction) is a recognized idea (extended mind / extended phenotype), but as a *coordinate* on the
   competency map — orthogonal to the state-content ladder — it sharpens the taxonomy.
8. **A taxonomy of coordination MEDIA (multi-agent cluster).** The three genuine multi-agent coordination
   competencies use three distinct media: a symbol **channel** (#10 communication), direct **observation**
   of the partner (#11 role-division, #12 contest), and the **environment** (#15 stigmergy). Each is
   isolated by a *different* ablation (channel-scramble / blind-partner / no-trail). "Coordination via
   channel vs observation vs stigmergy" is a known distinction in the literature, but the clean
   ablation-separable 3-way decomposition in one framework is a tidy structural result.

**Next concrete step (Track 2):** rigorously test #5. Take 3–4 competencies spanning the state-content
ladder (memory → counting → sequencing → rule-inference) and run EACH through the SAME battery of learners
{FSM-mutation, OpenAI-ES, SNES, gradient PPO} on matched tasks, measuring the reachability threshold per
competency. If the reachability order is lawful and tracks the state-content ladder, we have a candidate
regularity to literature-gate hard. This is cheap (reuses existing harnesses) and directly attacks the
unknown-commonality target. Track 1 (scale for a new competency) proceeds in parallel if resourced.

## Corpus note — competency #18 (abstraction) adds a GRADED datapoint
#18 (relational abstraction) is the FIRST **partial/graded** competency in the corpus: train 0.88 vs
held-out 0.59 (chance 0.33). All prior competencies were binary (ablation collapses it, or it doesn't). This
introduces a possible commonality-hunt angle: **competency as a graded generalization property**, where the
train->test GAP measures HOW FAR a learned mechanism transfers beyond its training support. Candidate law
(under-named, needs the hard literature gate): for a generic (non-relational-architecture) learner, emergent
relations transfer to held-out instances **partially and with a characteristic gap** rather than fully or
not-at-all. Testable via a capacity / relational-inductive-bias sweep (does the gap shrink with a comparison-
structured net?). Filed alongside #5 (substrate-reachability ordering) and #6 (sufficiency-pressure). Track 2
stays deferred per directive ('keep searching first'); this is logged, not pursued.

## Corpus note — audit sweep complete (19 competencies); the minimal-mechanism pattern holds
With #17-#19 the audit sweep is done: 19 catalogued competencies, ALL known to science. The central
regularity survives the extended search: agents find the MINIMAL SUFFICIENT mechanism for the demanded task,
and the minimal mechanism is always something science already named. #19 (imitation) is a fresh datapoint --
the learner solves a hidden-goal task by the minimal available means (copy the demonstrator's choice / goal-
emulation) rather than any richer social inference (no theory-of-mind needed). This strengthens candidate law
#6 (sufficiency-pressure): richer competencies (ToM, causal, planning) do not emerge because the toy tasks do
not make them the MINIMAL sufficient solution -- emulation suffices, so emulation is what appears. Prediction:
those richer competencies require worlds where the minimal sufficient solution IS the richer mechanism
(e.g. partners whose goals can only be predicted, not observed) -- i.e. the SCALE fork. Track 2 stays deferred
per directive; logged, not pursued.
