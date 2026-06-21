# T2d — External Demonstration: the instrument pointed at a real multi-agent LLM swarm

**Date:** 2026-06-21. The headline external application (instrument roadmap T2d):
point the full catalog battery, **blind**, at a genuinely external multi-agent
system — a swarm of independent LLM agents — and report its emergent-pattern
profile. One well-executed external application is the proof that turns a
validated taxonomy into a working instrument.

## Setup

`analysis/t2d_llm_swarm.py`: N = 24 independent LLM agents (Claude Haiku 4.5) move
on a 2-D torus. Each agent, every round, receives **only its local egocentric
neighborhood** as text — the bearings, distances and relative headings of the
agents within its sensing radius — and outputs a single turn. No agent ever sees
the global state, and the only coordination is whatever the agents do in response
to the local rule in their prompt. Any global order is therefore genuinely
**emergent from local LLM decisions**, not imposed.

Three conditions differ only in the one-line local rule. Each resulting trajectory
was profiled by `analysis/battery_profile.profile_observation` **without telling the
battery which condition produced it**, at the T2c-recommended OOD operating point
(`match_min_tier="confirmation"`) and at the native gate.

## Result — three swarms, three verdicts

| swarm | local rule | what emerged | instrument verdict |
|---|---|---|---|
| **align** | "match your neighbors' average heading" | polar order 0.31 → **0.97** (flocked) | **MATCH P5 (flocking)** — fires at screening; top candidate under ≥confirmation |
| **segregate** | "steer toward your own type; stay when surrounded by your kind" | same-type fraction 0.58 → **1.0** (segregated) | **EMERGENT-UNCLASSIFIED**, emergence **1.0** — strong self-organization, no substrate-matching detector |
| **random** | "ignore neighbors, turn randomly" | polar order ~0.15 (disordered) | **NO-EMERGENCE**, emergence 0.17 |

This exercises **all three instrument verdicts on a real LLM multi-agent system**:

- **MATCH** — a specific catalogued pattern: the alignment swarm is identified as
  flocking (P5). The signal is unambiguous: polarization 0.96, Cohen's d 32.2,
  shuffle-null p 0.005, calibrated confidence 0.91.
- **EMERGENT-UNCLASSIFIED** — the "none-of-the-above" / discovery channel: the
  segregation swarm self-organizes unmistakably (same-type fraction → 1.0,
  emergence 1.0) but matches no specific detector's substrate, so the instrument
  flags real emergence without forcing a label.
- **NO-EMERGENCE** — the random control is rejected; the instrument does not
  hallucinate a pattern where none exists.

The instrument also **discriminates by emergence kind**: orientational order
(flocking) for align vs positional clustering (aggregation) for segregate vs
neither for random — read off the same trajectories blind.

## Honest notes

1. **Flocking detection is real but reaches screening tier on this short trace.**
   P5's confirmation gate requires ≥ 10 crossing-times of data; 91 LLM rounds give
   ~6–14 depending on speed, enough to screen but borderline for confirmation. So
   at the conservative ≥confirmation OOD gate (the setting that gave T2c a zero
   false-MATCH rate), align is reported as strong emergence with P5 as the top
   (demoted) candidate rather than a firm MATCH. This is the same conservatism by
   design, applied honestly to an unfamiliar system; a longer run would lift P5 to
   confirmation.
2. **The aggregation swarm is EMERGENT-UNCLASSIFIED by substrate, not by failure.**
   The catalog's aggregation detector (P1, Schelling) requires a **lattice**
   substrate ("integer-labeled spatial data"); this swarm is in **continuous
   space**, so P1 rejects at prerequisite (`substrate_mismatch`). A spatial-binning
   observation adapter (T1a) was tried, but a 24-agent grid is too sparse to
   trigger P1. The honest reading: the swarm clearly self-organizes (emergence 1.0)
   and the instrument says so, while correctly declining to name a pattern whose
   substrate does not match — exactly the discovery signal the catalog is meant to
   emit. A denser swarm or a continuous-space aggregation detector would convert
   this to a firm MATCH P1.
3. **The control is clean.** Random turns give no order and the instrument returns
   NO-EMERGENCE (0.17) — no false-novel, no false-MATCH.

## Reproduce

- Simulator: `analysis/t2d_llm_swarm.py --condition {align,segregate,random} ...`
  (pure stdlib + anthropic; reads ANTHROPIC_API_KEY; Haiku by default).
- Profiler: `analysis/t2d_profile.py <trace.json> ...` (runs the battery at both
  gates). Traces preserved in the T2d savepoint.

## Bottom line

Pointed blind at real LLM agent swarms, the EPC battery reads their emergent
behavior — identifying flocking, flagging unclassified self-organization, and
rejecting noise. This is the external demonstration that turns the validated
32-pattern catalog from a taxonomy into a measurement instrument for multi-agent
AI systems.
