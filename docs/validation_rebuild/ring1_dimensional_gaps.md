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

## Method validated

The dimensional-gap approach produced a candidate on its first targeted cell: the
ontology *predicted* where to look (empty substrate), and a canonical system there was
flagged emergent-but-unclassified. Remaining gap targets queued: entropy-driven
ordering (`interact=none_entropy`) and the combination cells above.
