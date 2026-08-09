# Ring-3 catalog reconciliation — source-of-truth handback for Claude Code

One canonical numbered catalog, ready to apply to the repo. Everything gates KNOWN; the 0-new-to-science result holds.
The machine-readable catalog is `canonical_catalog.json` (array of {id, name, description, mechanism, diagnostic,
metric, gate, lit_name, state_required, state_mode}); this doc carries the numbering decision, the §5C text, and the
demo list.

## 1. Numbering decision (Task 1 + F)
**Clean-renumbered the six substrate faces from old #30–35 to #25–30.** #25–29 are NOT left vacant. Reason: at
publication a contiguous 1–30 is far cleaner than a 1–24 + gap + 30–35 catalog, and nothing external pins the old
numbers except the research-env findings `.md` files, which the map below reconciles.

**old→new map** (apply when ingesting the findings files):

| old | new | name |
|----|----|------|
| 30 | 25 | metabolic-commons foraging |
| 31 | 26 | niche construction / terraforming |
| 32 | 27 | non-stationary co-player adaptation |
| 33 | 28 | continuous-physics momentum control |
| 34 | 29 | morphology-conditioned specialization |
| 35 | 30 | compositional attention (joint conditioning) |

#1–24 are unchanged from your repo registry (your repo has #1–21; #22 division-of-labor, #23 group-size-regulation,
#24 stigmergic-coordination are the three you're missing and are in the canonical JSON).

## 2. Canonical catalog (Task A) — 30 competencies, all gate=known
Full fields in `canonical_catalog.json`. Summary (# — name — diagnostic ablation that collapses it):
1 navigation — away-from-goal barrier · 2 memory — cue ablation→chance · 3 accumulation/counting — supply-doubling
· 4 sequencing — order scramble · 5 delayed-gratification — reactive (state buys little) · 6 regulation — reactive
· 7 rule-inference(one-shot) · 8 rule-tracking(continual) · 9 conditional-selection(cued) · 10 communication —
mute channel · 11 role-division — block mutual obs · 12 contest-resolution · 13 fair-alternation · 14
count-communication(compositional) — mute · 15 stigmergy — block mark trail · 16 tool use/construction · 17
reputation/indirect-reciprocity — memwipe · 18 abstraction/systematic-gen(partial) · 19 imitation/social-learning ·
20 predictive-intention-reading(ToM) — freeze mover motion · 21 cumulative-culture-ratchet · 22 division-of-labor —
noid · 23 group-size-regulation — noid · 24 stigmergic-coordination — nomark · **25 metabolic-commons — zero
own-energy** · **26 niche-construction — revert edits** · **27 non-stationary-adaptation — zero opp-history** ·
**28 momentum-control — zero own-velocity** · **29 morphology-specialization — zero own-body (noid does NOT
collapse)** · **30 compositional-attention — nobody AND freezeblock both necessary**.

## 3. Learner boundaries (Task B) — LB1–LB3, NOT competencies
Full fields in `canonical_catalog.json["learner_boundaries"]`.
- **LB1 mixed-motive cooperation** — independent PPO defects (V 0.203); LOLA opponent-shaping → tit-for-tat (V 0.994).
- **LB2 periodic synchronization** — independent learner (K=0) fails; coupled-oscillator dynamics required.
- **LB3 one-channel necessity budget (recurrent bottleneck)** — a GRU holds ~one necessary own-state channel and
  drops the other; a Transformer holds both. This is the boundary complement of competency #30: #30 is the reached
  behavior, LB3 is the mapped limit and its resolution. Gated KNOWN (transformer-vs-recurrent POMDP memory lit).

## 4. known_reference_set.json (Task 3) — stays a fingerprint LIBRARY
35 entries, incl. negative controls (`relay-null`, `relay-multiseed-null`) and declined probes
(`deception-declined`, `anticipation-declined`, `tool-declined`). These are detector reference vectors, NEVER
promoted to competencies. Keep as-is; do not renumber to match the catalog.

## 5. §5C.4 map update (Task C) — state-requirement × state-mode
Each face placed on the two map axes. The budget-artifact meta-finding and #30's joint-conditioning mode are the
two structural additions.

**Reactive class (state_required = no):** delayed-gratification (#5), regulation (#6), niche-construction (#26).
**State-requiring, by state-mode:**
- *commitment*: navigation (#1)
- *storage*: memory (#2), reputation/social-memory (#17)
- *accumulation*: counting (#3)
- *belief (other-model)*: predictive-intention-reading (#20), non-stationary-adaptation (#27)
- *internal-state (own body/energy/velocity)*: metabolic-commons (#25), momentum-control (#28),
  morphology-specialization (#29)
- *externalized (world-state)*: stigmergy (#15), stigmergic-coordination (#24)
- **joint-conditioning (NEW mode)**: compositional-attention (#30) — the first face requiring TWO simultaneously-
  necessary channels; reachable only by changing the learner architecture (see LB3).

**Budget-artifact meta-finding (add to §5C.4 honest-status):** two faces (#28 momentum, #29 morphology) were
initially logged "untrainable" at modest budget (1500it / 800it) and both train cleanly at 8000it. The
reachable-competent set is a BUDGET-DEPENDENT lower bound, not a fixed frontier — the map's boundaries move with
training budget, and short-run negatives are not trustworthy negatives.

## 6. §5C.5 update text (Task D) — MCC co-evolution + bottleneck refinement
§5C.5 names environment/agent co-evolution as the most likely route to a novel competency. That experiment is now
run. We built Minimal Criterion Coevolution (Brant & Stanley 2017 + resource-limitation): two populations —
environments as evolvable genomes and GRU agents with a PPO inner loop — each reproducing only on a binary minimal
criterion, with no reward or objective. Over 400 generations the populations turned over continuously and produced
many competent behaviors, but a fair shared-channel re-screen against the catalog found ZERO novel-to-catalog (an
in-loop screen first flagged 1356 "candidates," all a support-mismatch artifact of comparing sparse MCC fingerprints
to the dense catalog; the honest re-screen resolves every one to an existing behavior). Removing the human from
objective-specification did not escape the learner-bounded reachable set.

This refines §5C's central claim. The limiting factor on discovering new competency is not the detector, and not the
environment alone: it is the **environment AND the learner architecture jointly**. LB3 is the clean demonstration.
Holding the world fixed and swapping GRU for a Transformer changes which behaviors are reachable, so the instrument
(the learner) is part of what determines the reachable phenomena. The un-crossed lever toward novelty is a richer
learner class, not a richer objective or world.

## 7. Track-2 literature-gate framing (Task 4) — incorporate into §5C.6
The genuine contribution is the interventional, self-debunking competency instrument, gated as genuinely UNDER-named.
Frame it as the interventional successor to Ronald, Sipper & Capcarrère 1999 "Design, Observation, Surprise!", which
replaces observer-surprise with survives-a-collapse-ablation as the emergence test. The killer demonstration is that
the instrument caught its own operator's ~9 over-claims (the auditor corrections logged across Ring-3). Differentiate
from: Levin (competency-by-perturbation, the philosophical basis), Schaeffer 2023 (emergent-abilities mirage),
Lowe/Foerster 2019 (pitfalls of measuring emergent communication), and Hamilton 2022 (eval-time ablation overstates
dependence, so we use fair separately-trained baselines).

**Honesty corrections (apply to §5C.6):**
- DROP "closed set / completeness" language — unfalsifiable at n≈14. Replace with: minimal-sufficient mechanisms are
  densely catalogued in the literature, so a toy-scale search is structurally biased toward rediscovering them.
- The minimal-mechanism "law" is ESTABLISHED (shortcut learning, Geirhos 2020; simplicity bias, Shah 2020;
  specification gaming, Krakovna; compression–expressivity, Kirby). Present the Ring-3 result as CONSOLIDATION of
  these, not discovery. Bound it explicitly to toy scale (tension with novelty-search/OEE, Lehman & Stanley 2011).
- Secondary candidates — present as PROPOSED, not validated: (a) coordination-media cost hierarchy +
  communication-niche (partially stated; needs a dialable-information-structure experiment); (b) shortcut-cascade-
  depth (under-named measure; needs a well-definedness proof); (c) temporal-scale axis (the axis is textbook
  Baldwin/Kirby/Henrich; only the cross-scale-invariance conjunction is under-named).

## 8. Demo-worthy list (Task E) — best candidates for a Ring-3 animated learned-behavior track
Ranked by visual legibility of the emergent behavior (separate track from the P-series physics gallery):
1. **#22 division-of-labor** — agents visibly split to distinct targets from a symmetric start (unstable attractor; show the seed spread honestly).
2. **#29 morphology-specialization** — red/blue bodies sorting onto matching food; occupancy 90–97% is a clean visual.
3. **#30 compositional-attention** — agent tracking BOTH its body-type and the opponent's switch; the joint-conditioning is the headline.
4. **#24 stigmergic-coordination** — blind agents covering a torus via mark trails only; the trail is inherently animatable.
5. **#28 momentum-control** — puck driven under inertia into a goal; continuous-physics motion reads well.
6. **#26 niche-construction** — agents building/locking shelters and ramps; before/after world-state is vivid.
(#25 metabolic-commons and #27 non-stationary-adaptation are scientifically clean but visually subtle — lower priority.)
