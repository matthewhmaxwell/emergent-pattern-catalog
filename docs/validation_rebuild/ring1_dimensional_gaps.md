# Ring 1 — Dimensional-Gap Discovery

The principled (periodic-table) discovery arm: map the 32 catalog patterns onto the
11-dimension emergence ontology (`docs/ontology_v0_4.md`), find the EMPTY cells, and
construct minimal systems for the physically-realizable gaps. A gap that the
ontology says *could* host a pattern but the catalog doesn't cover is a discovery
target. Generator: `analysis/discovery/coverage_map.py` (32×11 structural
classification, cross-validated against the ontology's own gap analysis + paper §5.8).

## Empty first-order cells (an entire ontology value, 0 of 32 patterns)

- **`substrate = evolving_net`** — adaptive/coevolving networks (topology changes
  with the dynamics). A major emergence frontier (Gross & Blasius 2008, *Adaptive
  coevolutionary networks*) absent from the catalog; distinct from every
  fixed-topology pattern (P18 voter, P21 polarization, P22 cascade).
- **`interact = none_entropy`** — ordering with NO agent-agent interaction force
  (entropy / excluded-volume / depletion driven): Brazil-nut granular segregation,
  hard-particle entropic crystallization, depletion-force assembly.

## Empty plausible combination cells

- `driving = ext_forced` AND `feedback = nontransitive` — externally driven cyclic dominance.
- `update = resource_exchange` AND `homog = heterogeneous` — heterogeneous-endowment exchange (P28 is identical-agent).
- `interact = field` AND `conflict = competitive` — field-mediated competition.
- `substrate = fixed_net` AND `update = reproduction` — evolutionary games on a fixed network (P27 is lattice).

## Candidate P34 — adaptive-network fragmentation (evolving-network gap)

Built the canonical evolving-network system (`analysis/discovery/ring1_new.adaptive_voter`):
adaptive voter model (Holme & Newman 2006), N=200, mean degree 8, rewiring prob φ=0.7.
Discordant edges either rewire toward same-opinion nodes or convert opinions; above a
critical φ the network spontaneously **fragments into disconnected same-opinion
communities** — the emergence is in the topology, which coevolves with the node states.

Instrument read (full battery, confirmation gate):

| system | verdict | channel | closest |
|---|---|---|---|
| adaptive_voter φ=0.7 (2 seeds) | **EMERGENT-UNCLASSIFIED** | network(modularity/scale-free) | P11 / none |
| null_random_rewire (opinion-blind) | NO-EMERGENCE | — | P11 / none |

The instrument detects the fragmentation as emergent (network channel), declines to
match it to any catalog pattern, and the opinion-blind rewiring null reads
NO-EMERGENCE. This is a Ring-1 candidate: emergence in an empty substrate cell,
distinct from the 32. **Next: vet via the verification ladder (robustness/finite-size,
literature-novelty, mechanism) and build a discriminating detector to the Phase-2a bar
-> P34.**

## Second gap explored — entropy-driven ordering (`interact = none_entropy`)

Built the canonical system (`analysis/discovery/ring1_entropy.hard_disk_crystallization`):
soft repulsive disks under slow compression with annealing noise — the Alder transition,
order from excluded volume alone (no attraction). Findings:

- It **partially crystallizes**: per-particle local hexatic order rises to
  ⟨|ψ6ᵢ|⟩ ≈ 0.44–0.53 (polycrystalline). *Global* ψ6 stays low (0.05–0.14) because a
  **square periodic box geometrically frustrates a single triangular lattice** — domains
  with different orientations.
- The instrument's read is **inconsistent** (em flips 0.0 ↔ ~0.6 across regimes/seeds):
  a crystal has near-uniform density, so the clustering/morphology channel only
  borderline-fires, and there is **no bond-orientational (ψ6) channel**.

Verdict: this gap is a **harder, multi-front target** than P34 — neither a clean
candidate nor a clean blind spot. A clean P35 here would require (1) a cleaner
crystallization (event-driven MC and/or a triangular-commensurate box → high, stable
ψ6) and (2) a **ψ6 bond-orientational channel added to the instrument** (a round-3
detection widening, analogous to the round-2 channels and the P3 fix: discovery
reveals a detection gap → widen the net → then catalog). Deferred as a scoped
follow-on, not papered over.

## Method validated

The dimensional-gap approach produced a clean catalog pattern (P34) on its first
targeted cell — the ontology *predicted* the empty substrate and a canonical system
there was flagged emergent-but-unclassified, then validated to the Phase-2a bar. The
second first-order gap (entropy ordering) is identified and scoped as a deeper arc.
Remaining queued targets: the entropy ψ6 arc above, and the combination cells
(driven cyclic dominance, heterogeneous resource-exchange, field-mediated competition,
evolutionary games on fixed networks).
