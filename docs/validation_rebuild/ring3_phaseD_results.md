# Ring 3 Phase D — the richer world (multi-agent): results

Validated rationale (after retracting the "low odds" over-claim): the single-agent toy ceiling is
*local* and does not transfer; multi-agent worlds add competency CLASSES that cannot exist with one
agent, so this is the highest-expected-value direction, not a long shot. Same loop as before
(open-ended generator → capable learner → discriminating baseline/ablation → agent-observer → literature
gate). Honest bound unchanged: literature gate before any novelty claim.

## D1 — substrate validation via emergent communication (`commworld_ppo.py`)
A referential signaling game in a gridworld: 3 colored goals; one is the hidden TARGET. The SPEAKER
sees the target color and emits a symbol (K=3); the LISTENER sees the symbol + directions to each
colored goal and navigates; shared reward for reaching the target-color goal. No built-in symbol→color
map — speaker and listener must CO-ADAPT a convention (Lewis signaling). Joint PPO trains both MLP
policies. A single agent cannot do this.

**Result (800 iters):** success 0.05 → **0.63** (plateau ~0.6–0.66). **Channel-scramble debunk** (feed
the listener a random symbol): success → **0.31 = chance (1/3)**. The collapse under scramble proves the
symbol genuinely carried the target information — **genuine emergent communication**, not a shortcut.
Imperfect (0.63, i.e. the co-adapted convention is partial/noisy), but real.

**Classification:** emergent communication / Lewis signaling — **known to science** (large emergent-comms
literature) → **0 new-to-science**. The significance is that it is the FIRST genuinely multi-agent
competency the pipeline has reached: it empirically confirms the richer world produces qualitatively new
competency classes (communication) the single-agent regime structurally could not, validating the
substrate + the channel-ablation discriminator for multi-agent work.

## Next — D2: the open-ended multi-agent hunt
Generalize from the hand-built referential game to a GENERATOR over multi-agent task structures
(coordination on a shared goal; anti-coordination / role-division; conditional/division-of-labor;
adversarial/deception variants; with/without a channel), train joint PPO, and have the agent-observer
NAME whatever multi-agent competency emerges (debunk via channel-ablation / role-swap / partner-swap),
classify vs the map or flag OFF-MAP, literature-gate any survivor. This is where the qualitatively
larger competency space actually gets searched for something not nameable in the existing literature.

## D2 — open-ended multi-agent sweep (`commhunt_ppo.py`)
Generator over 2-mobile-agent task structures; shared dual-head (move+symbol) policy; joint PPO;
discriminators = channel-scramble (collapse ⇒ uses COMMUNICATION) and blind-partner (collapse ⇒ uses
the PARTNER). 800 iters/task:
| task | success | scramble | blind | reading |
|---|---|---|---|---|
| referential | 0.34 | 0.33 | 0.32 | **under-solved** (chance) |
| coordination | 1.00 | 1.00 | 1.00 | shared **Schelling convention** — no real multi-agent demand |
| role_div | 0.99 | 1.00 | **0.44** | **GENUINE partner-dependent coordination** |
| independent | 0.34 | 0.33 | 0.32 | **under-solved** (chance) |

**The genuine finding — `role_div`:** scramble-invariant (not communication) but blind-partner collapses
it 0.99→0.44, so the agents genuinely **observe each other's positions and split roles** (anti-coordinate
by mutual adjustment). A real multi-agent competency, verified by the discriminating ablation, and
*distinct* from D1's communication (this one rides on observation, not the channel).

**Honest caveats:** `coordination` confirmed the Schelling-convention prediction (both head to the same
goal-index — no interaction needed; both ablations leave it at 1.00). `referential` and `independent`
stalled at chance: both require **private-target-conditioned navigation** (select the goal matching a
private color), which the shared 2-mobile-agent policy did not learn in budget — a tooling/difficulty
gap, not a comms-impossible result (D1 already established genuine communication in the clean
speaker–listener setup).

**Classification + bound:** role-division / mutual-adjustment is KNOWN (MARL role allocation) → **0
new-to-science**, a Tier-2 multi-agent face. Across D1+D2 the richer world has now yielded two genuine,
verified multi-agent competencies (communication; partner-dependent role-division), cleanly separated by
the scramble-vs-blind discriminators — empirically refuting the retracted "low odds" claim that the
richer world wouldn't produce qualitatively new competency classes.

**Next (toward off-map) — defeat the easy shortcuts:** both ways results stayed *known* are escapable by
**shared conventions** (coordination) or **simple mutual observation** (role_div). A probe designed to
defeat both — **anonymous agents** (remove the agent-id that lets them pre-split) with a **contested
resource** (both want one high-value goal; collision penalised) — forces genuine symmetry-breaking /
negotiation that neither a shared deterministic rule nor passive observation can solve. That is the next
place a multi-agent competency might resist every name in the literature. (It also requires fixing the
private-target-navigation tooling gap, or using D1-style asymmetric roles.)

## D3 — symmetry-breaking under contest (`contest_ppo.py`)
The probe built to defeat both shortcuts: ANONYMOUS agents (no id), one HIGH goal (+2) vs two LOW (+0.5),
both-high collides (−1), shared team reward. To win, who-takes-high must be decided WITHOUT id and WITHOUT
a trivially-shared rule. 1000 iters:
| metric | normal | chan-scramble | blind-partner |
|---|---|---|---|
| capture (one-on-high) | **0.63** | 0.65 | **0.14** |
| collision (both-high) | 0.02 | 0.01 | 0.16 |

**Result — genuine but KNOWN.** The agents resolve the contest by **breaking symmetry on relative
position**: collision is driven to 0.02, and the resolution is partner-geometry-dependent (blind collapses
capture 0.63→0.14 and raises collision) and NOT communication (scramble-invariant). This is exactly
**contest resolution via an uncorrelated asymmetry** — Maynard Smith's "Bourgeois" convention (evolutionary
game theory, 1973). KNOWN → **0 new-to-science**. It is also only *partial* (capture 0.63: they are
over-cautious, often both yielding so nobody claims high) — a sub-optimal Bourgeois convention, not a new
competency.

## D4 — forced negotiation / emergent fairness (`fairness_ppo.py`)
Removed the geometric asymmetry entirely: an abstract ITERATED contest between two FULLY SYMMETRIC agents
(shared recurrent policy, identical obs), selfish per-agent payoff (claim → +2 / +0.5; both-claim → −1;
iterated Hawk–Dove over 14 rounds). With no id and no geometry, the only symmetry-breakers are the channel
or memory.
- Greedy (deterministic) play: efficiency 0.00 — confirming the **symmetry-breaking obstruction** (a
  theorem: symmetric agents cannot coordinate deterministically without an asymmetry/correlation device).
- Stochastic play: **EMERGENT FAIR ALTERNATION** — efficiency 0.93, fairness 0.98, flip-rate 0.96 (each
  round exactly one claims; the prize split ~50/50). A genuine self-organized **fairness norm** from
  symmetric selfishness.
- Mechanism (ablations): channel-scramble 0.93 (NOT communication) AND memory-wipe 0.88 (NOT deep memory).
  The agents found an even simpler rule than designed for — **reactive anti-coordination**: flip my own
  last action, and the partner's last action (in my observation) keeps us out of phase, bootstrapped by a
  stochastic first round. A known repeated-game alternation/turn-taking equilibrium → **0 new-to-science**.

D4's twist sharpens the whole story: a probe built to *force* negotiation through the channel was solved by
the agents with a minimal reactive last-action rule — not the channel, not deep memory. They keep finding
the simplest mechanism, and the simplest mechanism is always named.

## Honest pattern across A–D (now including TWO probes built to force novelty)
Every probe — object-level (counting, memory, sequencing), meta-level (rule inference/tracking), and
multi-agent (communication, role-division, Bourgeois contest-resolution, emergent fair-alternation) — has
produced a GENUINE, verified competency that the literature already names. The two probes designed
specifically to *force* something unnamed — D3 (anonymous contested resource) and D4 (symmetric forced
negotiation) — BOTH still landed on classical results (uncorrelated-asymmetry resolution; repeated-game
turn-taking), because the agents find the MINIMAL mechanism that solves the task, and minimal-competency
space is densely mapped by science (ethology, control theory, game theory, ML). This is an empirical
observation across ~9 probes, not an extrapolation. The instrument works (genuine discovery + correct classification + 8+ self-caught
over-claims); the new-to-science target is genuinely hard for a reason we can now state: **simple,
minimal competencies are exactly the ones science has already catalogued.**

## Genuine remaining levers (a real fork, honestly stated)
1. **Force the channel** — remove the geometric asymmetry too (give both agents identical relative info),
   so the ONLY symmetry-breaker is the channel → forces *negotiated* turn-taking between truly symmetric
   agents (harder; may need stochastic protocols). Iterated contests could surface a *fairness/alternation
   norm* (richer than Bourgeois).
2. **Scale / complexity** — 3+ agents, iterated interaction, open-ended artifacts (tool construction),
   where *combinations* might exceed any single named competency. The genuine but expensive lever.
3. **Consolidate** — the validated instrument + a competency map spanning object/meta/multi-agent levels +
   a rigorous, well-evidenced account of why minimal competencies are all known, is a strong, honest,
   publishable contribution on its own.

