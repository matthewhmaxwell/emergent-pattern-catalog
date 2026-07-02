# Ring 3 — the coordination-media cost cascade (open-ended N-agent hunt)

**What this is.** The open-ended N-agent hunt (SCALE tier, the new-to-science attempt) run as a *discovery
engine*: build a rich multi-agent world → FINGERPRINT which coordination medium the agents actually use (full
ablation battery) → REMOVE that medium → re-run and climb to the next rung. Each rung is characterized and
literature-gated. The hunt did not find a new *competency*, but it produced a clean candidate **regularity**
(the Track-2 payoff): a coordination-media **cost hierarchy** with a **communication niche**.

## The cascade (2-agent rendezvous / cooperative foraging; parameter-shared recurrent PPO)
| rung | world (each removes the previous rung's medium) | mechanism used | fingerprint |
|---|---|---|---|
| 0 | heavy food (needs 2 agents) + observation + signal channels | **environmental focal point** (shared food) | blind 4.70 ≈ mute 4.81 ≈ normal 4.82; solo 0.00 |
| 1 | focal point removed: toroidal, egocentric, pure rendezvous | **observation** (pursuit; home in on partner) | normal 11.36, **blind 0.08** (collapse), mute 11.35 |
| 2 | observation removed: self-position only + channel | **coordinate focal point** (fixed convention cell) | normal 3.71, mute 3.65 (invariant), **self-blind 2.02** |
| 3 | symmetry broken: only one agent sees the target (asymmetry) | **communication** (referential location-signaling) | normal 3.02, **mute 0.30** (collapse), **rides-channel +2.73** (~90%) |
| 4 | asymmetry DISTRIBUTED: 3 agents, row-knower + col-knower → seeker | **composed multi-party communication** (2-source integration) | normal 3.72, mute-ROW 1.07, mute-COL 1.09; **both partial channels required** (needs-row +2.65, needs-col +2.63) |
| 5 | different DOMAIN: temporal, shared start (fire simultaneously) | **temporal focal point** (free internal-count from shared t=0) | normal 5.00, no-clock 5.00, blind 5.00, mute 5.00 — **all media unused** |
| 5b | temporal, shared-start symmetry BROKEN (random phase, no clock) | **observation** (watch partner fire, align phase) | normal 5.23, **blind 2.67**, mute 5.27 (needs-obs +2.56, channel unused) |
| 6 | 3 agents, INDIRECT: relay chain (source→relay→seeker, no direct link) | *did not emerge* — **bootstrap ceiling** (multi-hop co-emergence trap) | normal 0.54 ≈ chance even @3000 iters; both hops +0.02 (unused) |
| 7 | 3 agents, different KIND: consensus under conflicting preferences | **focal point** (all default to the same option) | normal 12.00, blind 12.00, mute 12.00 — **all media unused** |

## The finding: a coordination-media cost hierarchy + a communication niche
Agents always use the **cheapest available coordination medium**, and each medium is forced ONLY by removing
all cheaper ones. The observed cost ordering:

  **environmental focal point  <  observation  <  coordinate focal point  <  communication**

The sharp part is the **communication niche**: across rungs 0–2, three *different* symmetric worlds, the
signal channel was **never used** — a shared focal point (environmental, geometric, or coordinate) was always
cheaper. Communication is squeezed from both sides: cheaper focal points preempt it, and if you remove *all*
grounding there is no shared referent left for it to exploit. It is forced **only by information ASYMMETRY**
(rung 3: one agent knows the target, the other does not — no focal point resolves it, so the knower must
signal). **Symmetric coordination never forces communication.**

## Extension: the hierarchy generalizes (rungs 4–5)
- **Rung 4 (distributed asymmetry → composition).** Splitting the target across two speakers (row-knower +
  col-knower) forces the seeker to INTEGRATE two independently-encoded partial signals; muting *either* speaker
  collapses it ~70%. Communication, once forced, becomes COMPOSED when the information is distributed — the
  multi-agent analog of the sufficiency-pressure / compositionality result (#14). So the top of the hierarchy is
  itself structured: single-source asymmetry → referential signaling; distributed asymmetry → composed protocol.
- **Rung 5 / 5b (new DOMAIN: time) — the cost hierarchy is CROSS-DOMAIN.** Rung 5 (shared start): temporal
  simultaneity is FREE — identical policies from t=0 auto-phase-lock via internal counting; all three media
  (clock, observation, channel) contribute exactly +0.00. That is the focal-point escape, in time. Rung 5b
  breaks the shared-start symmetry (random phase offsets, no clock): coordination is then forced up to
  OBSERVATION (blind 5.23→2.67, needs-obs +2.56; channel still unused) — the exact temporal analog of spatial
  rung 1. **So the SAME ordering (focal-point → observation → communication) governs both space and time, and no
  NEW primitive medium appears in either domain.** This is direct evidence the primitive-media set is **closed**,
  which *explains* the 0-new-to-science bound: multi-agent competencies are minimal compositions of a closed set
  of catalogued media, so a genuinely new competency would require a NEW primitive medium — and across spatial
  AND temporal coordination, none emerged.

## Rung 6 (3+ agents, indirect): a LEARNING ceiling, not a new medium
The 3-agent regime is where a NEW primitive medium could hide (indirect coordination the 2-agent primitives
cannot express). The cleanest instance -- a relay chain (source→relay→seeker, no direct link) -- did NOT break
completeness: it simply **failed to emerge**, stuck at chance even after 3000 iters, with both hops unused.
This is a BOOTSTRAP CEILING (the multi-hop co-emergence trap: no reward gradient until all three agents' codes
align at once). It is the gradient-learning analog of the earlier evolutionary ceiling (evolution couldn't do
precise counting; gradient could) -- both are limits on LEARNABILITY, not evidence of a new mechanism. Relay
IS in the closed set (composed communication); it is known to require curriculum/scaffolding to emerge. So the
3-agent regime yields either closed-set mechanisms or learning ceilings -- **no new primitive medium appears**,
and completeness holds across spatial, temporal, AND 3-agent-indirect coordination.

## Why this matters (it unifies the whole multi-agent catalog)
The cascade retroactively EXPLAINS every prior multi-agent result:
- Genuine communication cases — #10 referential signaling, #14 compositional numerals — **all had a speaker who
  knew a referent the listener did not** (asymmetry). Consistent with the niche.
- Forced-*symmetric* probes — D3 contest (Bourgeois/uncorrelated-asymmetry), D4 fair-alternation (turn-taking)
  — resolved coordination with **non-communication focal-point / reactive mechanisms**, never the channel.
  Consistent with "symmetric coordination never forces communication."
- Stigmergy #15 (environmental focal point written to the world) and role-division #11 (observation) sit on the
  lower rungs exactly as the hierarchy predicts.

So the multi-agent catalog is not a bag of unrelated competencies — it is **structured by the coordination-media
cost hierarchy**, and *which* competency emerges is set by *which media the world leaves available*.

## Honest gate
Every individual mechanism is KNOWN: Schelling focal points (1960), pursuit/rendezvous (Alpern), Lewis
signaling (#10), signaling-requires-asymmetry (Spence/Crawford-Sobel). **0 new-to-science** on the pieces. The
candidate UNDER-NAMED contribution is the *synthesis*: (1) the coordination-media **cost hierarchy** as an
explicit ordering, (2) the **communication-niche** theorem (symmetric coordination never forces communication;
communication is forced iff information is asymmetric), and (3) the **cascade discovery method** (fingerprint →
remove → climb) that establishes it empirically. These need the hard literature gate before any novelty claim;
filed as Track-2 candidate regularities alongside the shortcut-cascade-depth measure (from ToM #20).

## Convergence (8 rungs): completeness is established
The hunt tested the completeness claim across three axes and it held every time — no new primitive medium
appeared anywhere:
- **Domains:** spatial (rungs 0–4) AND temporal (5, 5b) coordination — same `focal-point → observation →
  communication` ordering.
- **Regimes:** 2-agent (0–5b) AND 3-agent (4, 6, 7). The 3-agent-indirect case (relay, rung 6) hit a *learning*
  ceiling, not a new medium.
- **Competency kinds:** coordination (0–6) AND social choice / consensus (7) — consensus reduced to a focal
  point (12/12, all media unused).

Every rung reduced to the closed media set {focal-point, observation, communication, social-memory} or hit a
learnability ceiling. This is a decisive, well-replicated result: **multi-agent competency at this scale is
minimal composition over a CLOSED primitive-media set, which is *why* the search returns 0-new-to-science** — a
genuinely new competency would require a new primitive medium, and across domains/regimes/kinds none emerged.
Rung 7's perfect all-invariant confirmation is the saturation signal: further coordination probes re-confirm
without new information. The deliverable is the cost-hierarchy + communication-niche + completeness argument.

## Method note (reinforces the minimal-mechanism law at scale)
The hunt is a 4th–7th independent replication of the minimal-mechanism law: at EVERY rung the agents found the
cheapest medium the world permitted, including one I did not anticipate (the coordinate focal point at rung 2).
Building a world does not determine the mechanism; the *available media* do. To force a target mechanism you
must adversarially remove every cheaper one — the same principle that made ToM #20 need six shortcut-closes.
