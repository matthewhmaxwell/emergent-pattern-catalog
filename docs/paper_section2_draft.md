# Section 2: The Pattern Catalog

This section describes the pattern catalog: its three-layer architecture
(2.1), design rules for Layer 1 inclusion (2.2), the 11-dimensional
ontology used for classification (2.3), and a survey of the 32 atomic
patterns organized by cluster (2.4). Complete per-pattern entries,
including canonical references, detection metrics, and dimensional
profiles, are provided in the accompanying catalog document.

## 2.1 Three-Layer Architecture

The catalog is organized into three layers, each serving a distinct
function.

**Layer 1 — Atomic patterns.** The core inventory. Each entry is a
single, distinct, detectable macro-level behavior with at least one
canonical minimal model and a quantitative detection signature. This is
the periodic-table proper. Layer 1 currently contains 32 patterns
distributed across 10 clusters. Every detection, every transfer-matrix
cell, and every replication study in this paper operates on Layer 1
entries.

**Layer 2 — Cross-cutting descriptors.** Annotations applied across
Layer 1 rows, not row-defining criteria. Layer 2 splits into two bands.

Layer 2A contains *mathematical descriptors* — structural phenomena
such as phase transitions, symmetry breaking, critical regimes,
attractor dynamics, and nontransitive interaction — that describe *how*
atomic patterns appear, disappear, or relate to dynamical-systems
theory. These apply across multiple Layer 1 rows (e.g., phase
transitions annotate P5, P8, P9, P20, P21) and are used for
classification rather than as standalone catalog entries. A
demonstration that a model exhibits a phase transition is not a Layer 1
detection; the phase transition annotates whichever Layer 1 behavior is
undergoing the transition.

Layer 2B contains *cognitive-analogue descriptors* — interpretive
lenses from the Diverse Intelligence framework (memory, goal-seeking,
decision-making, perception, gratification deferral, collective
intelligence, self/other distinction) that annotate atomic patterns
with cognitive interpretations. These are *observer-dependent*
attributions applied across rows, not row-defining criteria. The
decision to treat cognitive-analogue labels as Layer 2B annotations
rather than Layer 1 entries is deliberate: behaviorally identical
patterns should occupy the same catalog row regardless of whether a
particular observer finds a cognitive interpretation useful. Layer 2B
bridges the behavioral ontology to Levin's TAME framework without
embedding cognitive claims in the detector logic.

**Layer 3 — Meta-capacities.** Higher-order organizational principles
that describe how systems *generate* atomic patterns rather than being
single detectable behaviors. Current Layer 3 entries include stigmergic
coordination (the interaction modality enabling several Layer 1
patterns that use environmental traces), multiscale competency (the
central organizing concept of Levin's TAME framework, explaining why
atomic patterns compose into collective intelligence at higher scales),
emergent goal shift (system-level objectives diverging from individual
agent goals), and productive forgetting (information loss that improves
system performance). Layer 3 is the theoretical context of the catalog,
not part of the detectable inventory.

**Boundary rule.** Layer 2 terms are descriptors that can be assigned
row-by-row across many atomic patterns. Layer 3 terms are explanatory
principles about generation, composition, or interpretation and are not
assigned as row attributes. This boundary rule disambiguates edge
cases: "critical regime" is Layer 2A (it describes patterns P14 and
P15), while "multiscale competency" is Layer 3 (it explains how Layer 1
patterns compose, but is not a property of any single pattern).

## 2.2 Design Rules for Layer 1 Inclusion

The catalog admits new entries only under five design rules, applied
jointly.

*A behavior, not a mechanism.* Entries must describe a macro-level
*behavior* — what the system does at the collective scale — not a
mechanism, interaction type, system property, or operating regime.
"Phase transition" fails this rule: it describes how a behavior
appears, not the behavior itself. "Polarization" passes: the behavior
is a durable split of the population into distinct camps.

*Grounded in at least one minimal model.* The pattern must be
observable in at least one minimal computational model (cellular
automaton, particle system, grid agents, sorting network, Boolean
network, or similar). This excludes candidates that appear only in
elaborate domain-specific simulations whose behaviors depend on many
interacting mechanisms.

*Quantitatively detectable.* The pattern must have at least one
plausible quantitative detection metric or computational signature.
Metrics may be direct (e.g., polarization order parameter for
flocking), indirect (e.g., likelihood ratio against a power-law null
for self-organized criticality), or structural (e.g., transfer entropy
across collisions for computational dynamics). A candidate whose
signature is only qualitative — observable by eye but not by metric —
is deferred until a detection protocol is proposed.

*Non-redundant.* Each entry must be distinct from every other entry by
at least one of: different observable signature, different necessary
structural conditions, or different region of ontological space.
Candidates that collapse to existing patterns under the current
detector definitions are documented as mechanistic variants of the
parent pattern (e.g., entropic segregation is a variant of P1
similarity-driven aggregation). Explicit non-redundancy tests are
required for borderline cases; P31 (delayed gratification) was
initially flagged provisional and only admitted after a 600-run
cross-validation test established that P31 features capture sorting
dynamics information orthogonal to P1 features (Section 6).

*Documented in peer-reviewed literature.* Each entry cites at least one
canonical demonstration. This is an empirical check against theoretical
patterns that sound plausible but have never been operationalized in a
working model. It does not guarantee correctness of the canonical
demonstration itself — several patterns have published demonstrations
we could not replicate exactly, as discussed in Section 4 — but it
ensures the pattern is not invented de novo for this catalog.

## 2.3 The Ontology: 11 Dimensions

Each Layer 1 pattern is classified along 11 ontological dimensions.
Dimensional profiles in the catalog list the key dimensions relevant to
each entry rather than exhaustively specifying all 11 axes for every
entry; most entries are characterized by four to six salient
dimensions.

1. **Spatial scale.** Local (single-agent behavior), mesoscale
   (clusters, neighborhoods), or global (system-wide).

2. **Temporal character.** Steady-state attractor, transient,
   oscillatory, or propagating. The choice is behaviorally important:
   P13 excitable waves are propagating; P16 associative memory is
   steady-state attractor; P11 predator-prey oscillation is oscillatory.

3. **Interaction type.** Direct (agent-to-agent), indirect-stigmergic
   (through environmental modification), global field (response to
   shared external signal), or none (entropy-driven, no interaction).

4. **Interaction substrate.** Continuous space, discrete lattice, fixed
   network, evolving network, or well-mixed compartment. This
   dimension is the most operationally important for the detection
   toolkit, since it determines detector dispatch (Section 3.4).

5. **Agent homogeneity.** Identical or heterogeneous. Heterogeneity
   can be in type (predator vs prey in P11), in parameters (informed
   vs naive agents in P19), or in information (different internal
   states).

6. **Goal structure.** Goalless (no objective function), implicit
   (behavior consistent with goal-seeking without explicit objective),
   or explicit (programmed objectives).

7. **Feedback structure.** Positive (amplifying), negative
   (dampening/corrective), mixed, or nontransitive (cyclic dominance,
   A > B > C > A).

8. **Memory.** Memoryless, local state (agent carries internal state),
   or environmental trace (information persists outside agents).

9. **Conflict structure.** Cooperative, competitive, mixed-motive, or
   nontransitive.

10. **External driving.** Autonomous, externally forced, periodically
    driven, or field-coupled.

11. **Update mode.** Motion-transport (agents relocate),
    state-transition (agents change state on fixed substrate),
    reproduction-selection (strategies replicate differentially), or
    resource exchange (conserved quantity transfers between agents).

**Dimensional coverage analysis.** The dimensional profile of the
current 32-entry Layer 1 is not uniformly distributed. Several regions
of the dimensional space are underpopulated or empty:

- The *evolving network* substrate has no current Layer 1 entries.
  Network co-evolution with behavior is a known research frontier with
  no agreed canonical minimal model.

- The combination of *externally forced* plus *nontransitive feedback*
  has no entries. Whether systems with periodic driving and cyclic
  dominance produce qualitatively distinct patterns is an open question.

- *Resource exchange* under *heterogeneous agents* has no entries —
  P28 (wealth condensation) uses identical agents, but
  heterogeneous-endowment variants may produce distinct signatures.

- The combination of *heterogeneous agents*, *global field*, and *no
  agent-agent interaction* contains only the deferred Brazil Nut Effect
  candidate, suggesting either a genuine gap or that patterns in this
  region are sensitive to substrate details that make them difficult to
  operationalize.

These gaps inform future discovery efforts: a new candidate pattern
that fills an empty region is more likely to be genuinely novel than
one that lands on a densely populated dimensional profile.

## 2.4 Layer 1: The 32 Atomic Patterns

The 32 atomic patterns are organized into 10 clusters by dominant
behavioral character. Cluster assignment is a convenience for
presentation, not a structural claim; a pattern's true identity is its
dimensional profile and detection signature, not its cluster. The
complete per-pattern entries with canonical references, detection
metrics, key dimensions, and distinctness statements from neighboring
patterns are provided in the accompanying catalog document.

### Cluster A: Spatial Organization

- **P1 — Similarity-driven aggregation.** Agents with similar
  properties cluster spatially without explicit grouping rules, driven
  by local preference. Canonical models: Schelling segregation; Zhang
  cell-view sorting. Detection: spatial autocorrelation (Moran's I),
  cluster size distribution, segregation index.

- **P2 — Activity-induced phase separation (MIPS).** Dense and dilute
  phases emerge among self-propelled agents without attractive forces,
  because agents accumulate where motion slows. Canonical model: active
  Brownian particles with density-dependent speed reduction.

- **P3 — Turing pattern formation.** Spatial periodic patterns
  (stripes, spots, labyrinths) arise from reaction-diffusion dynamics.
  Canonical models: Gierer-Meinhardt activator-inhibitor; Gray-Scott.
  Detection: 2D FFT dominant spatial frequency, wavelength stability.

- **P4 — Territoriality / exclusion boundaries.** Mobile agents using
  transient environmental marks create stable mutually exclusive
  domains. Canonical model: territorial random walks with scent
  marking and decay (Giuggioli et al., 2011).

### Cluster B: Collective Motion

- **P5 — Translational alignment (flocking).** Agents achieve global
  directional alignment through local velocity-matching interactions
  only. Canonical models: Reynolds Boids; Vicsek. Detection:
  polarization order parameter φ, group velocity relative to
  individual speed.

- **P6 — Milling / vortex formation.** A moving group spontaneously
  organizes into a rotating torus with high angular order and near-zero
  net translation. Canonical model: D'Orsogna-Chuang-Bertozzi-Chayes
  self-propelled particles with Morse-potential interaction. Distinct
  from P5 in that global order is circulatory, not translational.

- **P7 — Lane formation in counterflow.** Oppositely moving agents
  self-segregate into same-direction lanes. Canonical model: social-force
  pedestrian models; counterflow cellular automata.

- **P8 — Self-organized jamming / stop-go waves.** Local fluctuations
  amplify into backward-traveling congestion waves without any
  bottleneck. Canonical model: Nagel-Schreckenberg traffic CA.

### Cluster C: Temporal Dynamics and Synchronization

- **P9 — Temporal synchronization.** Agents with individual oscillatory
  dynamics spontaneously entrain rhythms without a central clock.
  Canonical models: Kuramoto coupled oscillators; Peskin pulse-coupled
  oscillators. Detection: Kuramoto order parameter r(t), critical
  coupling strength K_c.

- **P10 — Chimera states.** Identical coupled oscillators split into
  coexisting synchronized and desynchronized subpopulations despite
  homogeneous rules. Canonical model: nonlocally coupled Kuramoto
  oscillators on a ring. Distinct from P9 in coexistence of coherence
  and incoherence rather than uniform entrainment.

- **P11 — Predator-prey oscillations.** Coupled populations exhibit
  sustained oscillatory dynamics without external forcing. Canonical
  model: stochastic lattice Lotka-Volterra (Mobilia-Georgiev-Täuber
  2007). Detection: anti-phase cross-correlation of species fractions
  at nonzero lag, FFT peak-to-mean of population time series.

- **P12 — Cyclic dominance / coexistence waves.** Nontransitive local
  interactions (A beats B, B beats C, C beats A) generate rotating
  dominance and long-lived coexistence that would collapse under
  well-mixed conditions. Canonical model: spatial rock-paper-scissors
  (Reichenbach-Mobilia-Frey 2007). Distinct from P11 in that the
  interaction graph is nontransitive (≥3 species, cyclic), not
  bilateral.

### Cluster D: Wave Propagation and Excitable Dynamics

- **P13 — Excitable spiral and target waves.** Propagating trigger
  waves and rigidly rotating spirals self-sustain in excitable media
  with refractory-recovery dynamics. Canonical models: Greenberg-Hastings
  excitable CA; Belousov-Zhabotinsky. Distinct from P15 in that the
  core object is propagation geometry and refractory dynamics, not
  information processing.

- **P14 — Self-organized criticality.** System self-tunes to a critical
  state with power-law-distributed avalanches without parameter
  tuning. Canonical models: Bak-Tang-Wiesenfeld sandpile;
  Drossel-Schwabl forest fire. Detection: power-law exponent of
  avalanche size distribution, 1/f temporal noise.

### Cluster E: Information Processing and Collective Cognition

- **P15 — Persistent propagating computation.** System sustains
  information-bearing structures whose collisions transform, route, or
  compute. Canonical models: Conway's Game of Life; random Boolean
  networks at critical connectivity. Detection: measurable information
  transfer between spatial regions (transfer entropy > 0 across
  collisions), sustained non-trivial transient length. Distinct from
  P13 in that propagating structures must carry and transform
  information, not merely propagate.

- **P16 — Associative memory / pattern completion.** System stores and
  retrieves distributed patterns through attractor dynamics. Canonical
  models: Hopfield network; Boolean gene regulatory networks with
  multiple fixed-point attractors.

- **P17 — Distributed sensing / collective gradient detection.** A
  group extracts directional information that no individual agent can
  detect alone. Canonical models: Berdahl et al. mobile-group gradient
  sensing; Camley et al. collective chemotaxis.

### Cluster F: Decision-Making and Social Dynamics

- **P18 — Collective consensus / decision-making.** Group converges to
  a shared choice without central authority. Canonical models:
  honeybee nest-site selection; voter model with competing options.

- **P19 — Emergent leadership / minority guidance.** A small informed
  subset steers a larger group's collective behavior without explicit
  signaling. Canonical model: Couzin et al. (2005) informed-minority
  flocking.

- **P20 — Quorum sensing / threshold-activated collective response.**
  System-level behavioral switch triggered at critical density — a
  binary toggle, not graded consensus. Canonical models: bacterial
  quorum sensing; agent threshold-activation. Distinct from P18 in
  being a density-dependent on/off toggle rather than resolution among
  alternatives.

- **P21 — Polarization / fragmentation.** Population splits into
  persistent internally coherent camps rather than converging or
  contagion-spreading. Canonical model: Hegselmann-Krause
  bounded-confidence opinion dynamics. Detection: persistent
  multimodality, between-cluster distance exceeding confidence bound.

- **P22 — Information cascade / social contagion.** Behavior or
  information propagates in self-amplifying waves that overshoot the
  original trigger then burn out. Canonical models: SIR/SIS epidemic
  models on networks; Watts threshold contagion. Distinct from P21 in
  being one-directional spreading that consumes available fuel rather
  than a self-stabilized partition.

- **P23 — Anti-coordination / emergent load balancing.** Agents
  repeatedly self-distribute across options to avoid overcrowding,
  producing utilization near capacity without central assignment.
  Canonical models: Minority Game (Challet-Zhang 1997); El Farol Bar
  (Arthur 1994).

### Cluster G: Resilience and Regulation

- **P24 — Homeostatic regulation.** System maintains an internal
  variable near a set-point despite sustained external perturbation
  through active corrective feedback.

- **P25 — Canalized restoration (equifinality).** System converges to
  the same target macrostate from a wide range of initial conditions
  or after diverse perturbations. Canonical models: planarian
  regeneration; Waddington landscape; developmental GRN attractor
  basins. Distinct from P24 in that P24 is cybernetic (continuous
  regulation) while P25 is geometric (convergence from arbitrary
  starting points).

- **P26 — Stochastic resonance.** Noise improves system performance at
  an intermediate optimal level; too little noise yields no signal
  detection and too much destroys it, with a sweet spot that amplifies.

### Cluster H: Competition and Cooperation

- **P27 — Spatial reciprocity / emergent cooperation.** Locally
  interacting selfish agents generate persistent islands of cooperation
  that survive exploitation, even when well-mixed payoff structure
  predicts cooperation's extinction. Canonical model: Nowak-May (1992)
  spatial Prisoner's Dilemma on a lattice. Distinct from P1 in that
  clusters emerge through differential survival on a static lattice,
  not physical motion.

- **P28 — Wealth condensation / spontaneous inequality.** Symmetrical
  fair exchanges of a conserved resource between identical agents
  produce highly skewed winner-take-all distributions. Canonical models:
  Yard-Sale (Boghosian 2017); Bouchaud-Mézard.

### Cluster I: Emergent Structure Formation

- **P29 — Trail / network formation.** Efficient transport networks
  emerge from collective movement and environmental modification
  without global planning. Canonical models: ant trail optimization;
  Physarum polycephalum network models; active walker models.

- **P30 — Spontaneous boundary formation (autopoiesis).** Agents
  organize into closed, semi-permeable topology maintaining an internal
  micro-environment distinct from the exterior. Canonical models:
  Varela-Maturana-Uribe SCL; computational lipid bilayer simulations.

### Cluster J: Emergent Agent-Level Competencies

- **P31 — Delayed gratification.** A statistically significant fraction
  of agents accept short-term cost (moving away from local goal) to
  achieve long-term benefit (global optimum), without explicit planning.
  Canonical model: Zhang cell-view sorting (2024). Detection: the DG
  index (fraction of steps moving away from eventual final position),
  cross-validated against baseline P1 features via the non-redundancy
  test (Section 6).

- **P32 — Emergent specialization (division of labor).** Identical
  agents spontaneously differentiate into distinct functional roles
  that increase collective efficiency. Canonical models: response
  threshold models (Bonabeau 1996); foraging task allocation in ant
  colonies.

## 2.5 Layer 2 Descriptors and Layer 3 Meta-Capacities

**Layer 2A — Mathematical descriptors.** Five annotations currently
apply across Layer 1 rows: phase transition / symmetry breaking (P5,
P8, P9, P20, P21); critical regime (P14, P15); positive feedback
amplification (P1, P2, P8, P22, P28); attractor dynamics (P16, P25);
nontransitive interaction (P12). These annotations are used for
classification and are not standalone detectable patterns.

**Layer 2B — Cognitive-analogue descriptors.** Seven annotations apply
across Layer 1 rows: memory (P4, P16); goal-directedness (P24, P25);
decision-making (P18, P23); perception / inference (P17, P26);
gratification deferral (P31); collective intelligence (P17, P18, P32);
self/other distinction (P30). Future work will formalize each
annotation as a predicate sheet specifying plain-language definition,
minimum evidence standard, disqualifying cases, exemplar rows, and
confidence level.

**Layer 3 — Meta-capacities.** Four entries: stigmergic coordination
(interaction modality enabling P4, P20, P29, P32); multiscale
competency (the central organizing concept of Levin's TAME framework);
emergent goal shift (system-level objectives diverging from individual
agent goals; pending operationalization); productive forgetting
(information loss improving system performance; manifests across P29,
P32). Layer 3 names the theoretical context of the catalog and is not
a detectable inventory.

## 2.6 Deferred Candidates and Open Questions

Six candidate patterns are currently deferred. Competitive coarsening /
Ostwald ripening exhibits t^(1/3) scaling of aggregate-level dynamics;
the literature treats it as a kinetic regime of clustering, and the
current catalog treats it as a mechanistic variant of P1. Swarm
morphogenesis / programmable self-assembly (Rubenstein et al. 2014
Kilobots) is the strongest deferred candidate, but its canonical
demonstrations require target specifications, making it unclear whether
the pattern is target-free self-assembly or goal-directed assembly.
Active nematics / topological defects — where dense elongated agents
form ±1/2 defects acting as mobile super-agents — is a major focus of
the active-matter literature but may be too substrate-specific.
Selfish-routing inefficiency / Braess paradox is at the edge of our
scope; canonical agents are often one-shot path selectors rather than
ongoing agents. Entropic segregation produces a macro-signature similar
to P1 and is treated as a mechanistic variant under behavior-first
rules. The Brazil Nut Effect / kinematic sorting occupies a rare
dimensional region (heterogeneous agents, global field, no agent-agent
interaction) but may be too domain-specific.

Four open architectural questions remain. The P31 non-redundancy
protocol established that P31 captures information orthogonal to P1
features (Section 6); the same protocol should be applied to other
borderline cases before they are admitted or declined. The P13/P15
boundary is operationalized via the boundary-conditioned transfer
entropy test (Section 6) and should be validated on additional
excitable and computational substrates as they are added. The Layer 2B
cognitive-analogue annotations require formalization as predicate
sheets. The detector specification cards — each pattern's required
observables, metrics, null model, pass threshold, false positives, and
neighbor exclusions — are the next major catalog deliverable; several
are drafted in `docs/detector_cards.md` and are refined in Section 3.

## 2.7 Summary

Layer 1 contains 32 atomic patterns across 10 clusters, each grounded
in at least one canonical minimal model and a quantitative detection
signature. Layer 2 provides 5 mathematical descriptors and 7
cognitive-analogue annotations applied across Layer 1 rows. Layer 3
names 4 meta-capacities that describe how atomic patterns are generated
and composed. The 11-dimensional ontology supports dimensional coverage
analysis and underpins the substrate-aware dispatch system of the
detection toolkit described in Section 3.
