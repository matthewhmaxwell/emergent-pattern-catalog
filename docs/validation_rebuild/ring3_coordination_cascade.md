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

## Method note (reinforces the minimal-mechanism law at scale)
The hunt is a 4th–7th independent replication of the minimal-mechanism law: at EVERY rung the agents found the
cheapest medium the world permitted, including one I did not anticipate (the coordinate focal point at rung 2).
Building a world does not determine the mechanism; the *available media* do. To force a target mechanism you
must adversarially remove every cheaper one — the same principle that made ToM #20 need six shortcut-closes.
