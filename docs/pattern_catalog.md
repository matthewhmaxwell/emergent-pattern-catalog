# Emergent Pattern Catalog v0.4

A systematic catalog of emergent behavioral competencies in minimal agent-based
systems. Each entry describes a macro-level behavior that arises without being
explicitly programmed, is observable in at least one minimal computational model,
and is detectable via quantitative metrics.

This version incorporates four rounds of external review across multiple models.

---

## Catalog Architecture

The catalog is organized into three layers:

**Layer 1 — Atomic patterns.** The core entries. Each is a single, distinct,
detectable macro-behavior with a canonical minimal model and at least one
computational signature. This is the "periodic table" proper. Current count: 32.

**Layer 2 — Cross-cutting descriptors.** Two sub-bands:

- **2A: Mathematical descriptors.** Structural phenomena (phase transitions,
  symmetry breaking, attractor dynamics) that describe *how* atomic patterns
  appear, disappear, or relate to dynamical systems theory. They apply across
  multiple entries and are used for classification, not as standalone catalog rows.
- **2B: Cognitive-analogue descriptors.** Interpretive lenses from the Diverse
  Intelligence framework (memory, goal-seeking, gratification, inference,
  decision-making) that annotate atomic patterns with cognitive interpretations.
  These are observer-dependent attributions applied *across* rows, not
  row-defining criteria. This band bridges the catalog to Levin's TAME framework
  without embedding cognitive claims into the behavioral ontology.

**Layer 3 — Meta-capacities.** Higher-order organizational principles that
describe *how systems generate* atomic patterns rather than being single
detectable behaviors. They frame the theoretical context of the catalog,
especially its connection to Levin's Diverse Intelligence program.

**Boundary rule:** Layer 2 terms are descriptors that can be assigned row-by-row
across many atomic patterns. Layer 3 terms are explanatory principles about
generation, composition, or interpretation and are not assigned as row attributes.

### Design Rules for Layer 1 Inclusion

1. Must be a macro-level *behavior*, not a mechanism, interaction type, system
   property, or operating regime
2. Must be observable in at least one minimal computational model (cellular
   automaton, particle system, grid agents, sorting network, Boolean network, etc.)
3. Must have a plausible quantitative detection metric or computational signature
4. Must be non-redundant — distinct from every other entry by at least one of:
   different observable signature, different necessary structural conditions, or
   different region of ontological space
5. Must be documented in peer-reviewed literature with at least one canonical
   demonstration

### Architecture Note: Cognitive Analogues as Descriptors

The Levin/Diverse Intelligence framework motivates interpreting certain patterns
through cognitive lenses (memory, goal-seeking, gratification, inference). In
v0.4, cognitive analogues are formalized as Layer 2B cross-cutting descriptors
applied across atomic patterns. This means: Layer 1 rows are behaviorally
defined; the cognitive interpretation is an annotation applied across rows. This
preserves the DI framing without muddying the behavioral ontology.

Future versions should formalize 2B annotations as predicate sheets with:
plain-language definition, minimum evidence standard, disqualifying cases,
exemplar rows, and confidence level.

---

## Ontological Dimensions

Each pattern is situated along the following axes. Dimensional profiles in the
catalog list key dimensions relevant to each entry rather than exhaustively
specifying all axes.

1. **Spatial scale:** local / mesoscale / global
2. **Temporal character:** steady-state attractor / transient / oscillatory / propagating
3. **Interaction type:** direct / indirect-stigmergic / global field / none (entropy-driven)
4. **Interaction substrate:** continuous space / discrete lattice / fixed network / evolving network / well-mixed compartment
5. **Agent homogeneity:** identical / heterogeneous
6. **Goal structure:** goalless / implicit / explicit
7. **Feedback structure:** positive / negative / mixed / nontransitive
8. **Memory:** memoryless / local state / environmental trace
9. **Conflict structure:** cooperative / competitive / mixed-motive / nontransitive
10. **External driving:** autonomous (fully endogenous) / externally forced / periodically driven / field-coupled
11. **Update mode:** motion-transport / state transition / reproduction-selection / resource exchange

---

## Layer 1 — Atomic Patterns

### Cluster A: Spatial Organization

#### P1 — Similarity-driven aggregation
Agents with similar properties cluster spatially without explicit grouping instructions, driven by local preference rules.
- **Canonical model:** Schelling segregation model; Zhang et al. cell-view sorting
- **Detection metric:** Spatial autocorrelation (Moran's I), cluster size distribution, segregation index
- **Key dimensions:** mesoscale / steady-state / direct / lattice / motion-transport / autonomous
- **Note:** Entropic segregation (size-based sorting with zero agent preference) produces a similar macro-signature and is treated as a mechanistic variant of P1.
- **Cognitive-analogue annotation:** Preference, collective intelligence

#### P2 — Activity-induced phase separation (MIPS)
Dense and dilute phases emerge among self-propelled agents even without attractive forces, because agents accumulate where motion slows.
- **Canonical model:** Active Brownian particles with density-dependent speed reduction
- **Detection metric:** Bimodal density distribution, phase coexistence with density-dependent velocity correlation, absence of attractive interaction in the rule set
- **Key dimensions:** mesoscale / steady-state / no direct attraction / continuous space / motion-transport / autonomous
- **Distinctness from P1:** No similarity preference or attraction rule. Clustering is purely kinetic.

#### P3 — Turing pattern formation
Spatial periodic patterns (stripes, spots, labyrinthine structures) arise from reaction-diffusion dynamics among agents or chemical species.
- **Canonical model:** Turing reaction-diffusion; Gierer-Meinhardt activator-inhibitor; Gray-Scott
- **Detection metric:** Dominant spatial frequency via 2D FFT, wavenumber selection, pattern wavelength stability
- **Key dimensions:** mesoscale / steady-state / local interaction / continuous space or lattice / state transition / autonomous
- **Distinctness from P1/P2:** Produces periodic spatial order (regular spacing), not mere clustering.

#### P4 — Territoriality / exclusion boundaries
Mobile agents using transient environmental marks create stable, mutually exclusive spatial domains with persistent boundaries.
- **Canonical model:** Territorial random walks with scent marking and decay (Giuggioli et al., 2011)
- **Detection metric:** Home range overlap index, boundary persistence over time, spatial cross-correlation between occupancy and foreign scent
- **Key dimensions:** mesoscale / steady-state / indirect-stigmergic / continuous space / motion-transport / environmental trace memory / autonomous
- **Distinctness from P1:** Outcome is spatial exclusion and boundary maintenance, not clustering.
- **Cognitive-analogue annotation:** Memory, goal-seeking

### Cluster B: Collective Motion

#### P5 — Translational alignment (flocking)
Agents achieve global directional alignment through local velocity-matching interactions only.
- **Canonical model:** Reynolds Boids; Vicsek model
- **Detection metric:** Polarization order parameter (mean alignment), group velocity magnitude relative to individual speed
- **Key dimensions:** global / steady-state / direct / continuous space / motion-transport / memoryless / autonomous

#### P6 — Milling / vortex formation
A moving group spontaneously organizes into a rotating torus or mill with strong angular order but near-zero net translation.
- **Canonical model:** D'Orsogna-Chuang-Bertozzi-Chayes self-propelled particles with attraction-repulsion
- **Detection metric:** Net angular momentum around group center, ring-like radial density profile, low center-of-mass velocity despite high individual speed
- **Key dimensions:** global / steady-state / direct / continuous space / motion-transport / autonomous
- **Distinctness from P5:** Global order is circulatory, not translational.

#### P7 — Lane formation in counterflow
Oppositely moving agents spontaneously segregate into same-direction lanes, increasing throughput and reducing collisions.
- **Canonical model:** Social-force pedestrian model; counterflow cellular automata
- **Detection metric:** Lane order parameter (local directional segregation), head-on collision frequency reduction, throughput gain vs mixed baseline
- **Key dimensions:** mesoscale / transient-dynamic / direct / continuous space or lattice / heterogeneous (opposing goals) / motion-transport / autonomous
- **Distinctness from P5:** Population does not converge on one heading — it self-sorts into oppositely directed streams.
- **Distinctness from P29:** Lanes are dynamic flow structures, not persistent built paths.

#### P8 — Self-organized jamming / stop-go waves
Small local fluctuations amplify into backward-traveling congestion waves or dynamic arrest without any bottleneck or central obstruction.
- **Canonical model:** Nagel-Schreckenberg traffic cellular automaton
- **Detection metric:** Backward wave propagation speed, bimodal speed distribution, oscillatory density bands in spacetime diagrams, breakdown in flow-density fundamental diagram
- **Key dimensions:** mesoscale / propagating / direct / lattice / motion-transport / positive feedback / autonomous
- **Distinctness from P14:** Signature is a coherent traveling congestion wave, not a scale-free avalanche distribution.

### Cluster C: Temporal Dynamics and Synchronization

#### P9 — Temporal synchronization
Agents with individual oscillatory dynamics spontaneously entrain rhythms without a central clock or pacemaker.
- **Canonical model:** Kuramoto coupled oscillators; Peskin pulse-coupled oscillators (firefly model)
- **Detection metric:** Kuramoto order parameter r(t), time to synchronization, critical coupling strength
- **Key dimensions:** global / steady-state / direct / network or continuous space / state transition / local state memory / autonomous

#### P10 — Chimera states
Identical coupled oscillators spontaneously split into coexisting synchronized and desynchronized subpopulations, despite homogeneous rules and coupling.
- **Canonical model:** Nonlocally coupled Kuramoto oscillators on a ring (Kuramoto & Battogtokh, 2002; Abrams & Strogatz, 2004)
- **Detection metric:** Coexistence of high and low local Kuramoto order parameters, spatial entropy of phase coherence, stable partial synchrony
- **Key dimensions:** mesoscale / steady-state / direct (nonlocal) / lattice or network / state transition / autonomous
- **Distinctness from P9:** Mixed order — coexistence of coherence and incoherence — not uniform entrainment.

#### P11 — Predator-prey oscillations
Coupled populations exhibit sustained oscillatory dynamics in abundance without external forcing.
- **Canonical model:** Stochastic Lotka-Volterra on a lattice; individual-based predator-prey models
- **Detection metric:** Period regularity of population cycles, amplitude stability, phase relationship between species, persistence time
- **Key dimensions:** global / oscillatory / direct / lattice or continuous space / heterogeneous (two types) / reproduction-selection / autonomous

#### P12 — Cyclic dominance / coexistence waves
Nontransitive local interactions (A beats B, B beats C, C beats A) generate rotating dominance, spatial spiral waves, and long-lived species coexistence that would collapse under well-mixed conditions.
- **Canonical model:** Spatial rock-paper-scissors (Reichenbach, Mobilia & Frey, 2007)
- **Detection metric:** Phase-lagged species abundance oscillations, spiral/domain morphology, coexistence time vs mobility, extinction probability
- **Key dimensions:** mesoscale / oscillatory-propagating / direct / lattice / heterogeneous (3+ types) / reproduction-selection / nontransitive feedback / autonomous
- **Distinctness from P11:** Interaction graph is nontransitive (cyclic, ≥3 species), not bilateral. Spatial structure is essential for coexistence.

### Cluster D: Wave Propagation and Excitable Dynamics

#### P13 — Excitable spiral and target waves
Propagating trigger waves and rigidly rotating spiral waves self-sustain in excitable media where agents have discrete refractory states. The core object is propagation geometry and refractory-recovery dynamics.
- **Canonical model:** Greenberg-Hastings excitable cellular automaton; Belousov-Zhabotinsky reaction-diffusion
- **Detection metric:** Wavefront propagation speed, spiral tip trajectory and topological charge, refractory-tail structure, persistence of self-sustaining wave sources
- **Key dimensions:** mesoscale / propagating / local interaction / lattice or continuous space / state transition / local state memory / autonomous
- **Distinctness from P9:** No requirement for phase locking across all agents.
- **Distinctness from P15:** P13 is about propagation geometry and refractory dynamics. P15 requires information-bearing signal interactions that transform, route, or compute. Discriminating metric: Transfer Entropy across structure collisions is approximately 0 for P13 (waves annihilate), > 0 for P15 (outputs depend on inputs).

#### P14 — Self-organized criticality
System self-tunes to a critical state characterized by power-law distributed avalanches/events, without external parameter tuning.
- **Canonical model:** Bak-Tang-Wiesenfeld sandpile; Drossel-Schwabl forest fire model
- **Detection metric:** Power-law exponent of avalanche size distribution, 1/f noise in temporal signal, finite-size scaling collapse
- **Key dimensions:** global / steady-state / local interaction / lattice / state transition / mixed feedback (slow drive + fast relaxation) / autonomous

### Cluster E: Information Processing and Collective Cognition

#### P15 — Persistent propagating computation
System sustains information-bearing structures (gliders, signals, logical gates) whose collisions and interactions transform, route, or compute — enabling open-ended information processing through ongoing structural dynamics. The key requirement is not merely persistence of moving structures, but that interactions between structures carry and transform information.
- **Canonical model:** Conway's Game of Life (glider collisions implementing logic gates); Langton's lambda-tuned cellular automata; random Boolean networks at critical connectivity
- **Detection metric:** Collision-based logic gate construction, measurable information transfer between spatial regions (transfer entropy > 0 across collisions), sustained non-trivial transient length relative to system size, demonstrated signal routing or transformation
- **Key dimensions:** global / propagating / local interaction / lattice / state transition / autonomous
- **Distinctness from P13:** P13 requires only propagation geometry and refractory dynamics. P15 requires that propagating structures interact in ways that carry, transform, or route information — computational function, not just wave geometry.
- **Note:** The critical regime that enables this behavior is treated as a Layer 2A cross-cutting descriptor, not as the behavior itself.

#### P16 — Associative memory / pattern completion
System stores and retrieves distributed patterns through attractor dynamics, completing partial inputs to stored templates.
- **Canonical model:** Hopfield network; Boolean gene regulatory networks with multiple fixed-point attractors
- **Detection metric:** Pattern completion accuracy from partial/noisy cues, basin of attraction size, storage capacity
- **Key dimensions:** global / steady-state attractor / direct / network / state transition / local state memory / autonomous
- **Cognitive-analogue annotation:** Memory, recall

#### P17 — Distributed sensing / collective gradient detection
A group extracts directional or environmental information that no individual agent can detect alone.
- **Canonical model:** Berdahl et al. mobile-group gradient sensing; Camley et al. collective chemotaxis
- **Detection metric:** Group-level chemotactic index or navigation accuracy that rises with group size while isolated individuals remain at chance
- **Key dimensions:** global / transient / direct / continuous space / motion-transport / field-coupled
- **Distinctness from P5:** Motion alone is insufficient — the key behavior is group-level inference about an external field.
- **Cognitive-analogue annotation:** Perception, collective inference

### Cluster F: Decision-Making and Social Dynamics

#### P18 — Collective consensus / decision-making
Group converges to a shared choice among alternatives without central authority, often exceeding individual accuracy.
- **Canonical model:** Honeybee nest-site selection (Seeley et al.); voter model with competing options
- **Detection metric:** Convergence time, decision accuracy vs individual baseline, sensitivity to quality differences between options
- **Key dimensions:** global / steady-state / direct / network or well-mixed / state transition / local state memory / autonomous
- **Cognitive-analogue annotation:** Decision-making, collective intelligence

#### P19 — Emergent leadership / minority guidance
A small informed subset steers a much larger group's collective behavior without explicit signaling, centralized control, or stable role assignment.
- **Canonical model:** Couzin et al. (2005) informed-minority flocking model
- **Detection metric:** Influence asymmetry: directional cross-correlation or transfer entropy between minority and majority; guidance efficacy of informed fraction relative to its population share
- **Key dimensions:** global / transient or steady-state / direct / continuous space / heterogeneous (informed vs naive) / motion-transport / autonomous
- **Distinctness from P18:** Disproportionate directional influence by a subset, not symmetric pooling.
- **Distinctness from P32:** Leadership role is transient and context-dependent, not a stable task identity.

#### P20 — Quorum sensing / threshold-activated collective response
System-level behavioral switch triggered when a critical density or agreement level is reached — a binary state toggle, not a graded consensus.
- **Canonical model:** Bacterial quorum sensing models; agent threshold-activation models
- **Detection metric:** Sharp transition in collective behavior at critical agent density, hysteresis in activation/deactivation thresholds
- **Key dimensions:** global / steady-state / indirect-stigmergic / well-mixed or continuous space / state transition / environmental trace memory / autonomous
- **Distinctness from P18:** Binary density-dependent toggle (on/off), not resolution among competing alternatives.

#### P21 — Polarization / fragmentation
A population splits into persistent, internally coherent camps under local interaction rules, rather than converging to consensus or spreading via contagion.
- **Canonical model:** Hegselmann-Krause bounded-confidence opinion dynamics
- **Detection metric:** Persistent multimodality in opinion distribution, between-cluster distance exceeding confidence bound, nonzero final variance after dynamics converge
- **Key dimensions:** global / steady-state / direct / network or well-mixed / state transition / autonomous
- **Distinctness from P18:** Macrostate is durable disagreement, not shared consensus.
- **Distinctness from P22:** Self-stabilized partition, not a one-way spreading event.

#### P22 — Information cascade / social contagion
Behavior or information propagates through a network in self-amplifying waves that overshoot the original trigger, then burn out.
- **Canonical model:** SIR/SIS epidemic models on networks; threshold models of social contagion (Watts, 2002)
- **Detection metric:** Cascade size distribution, propagation speed, final adoption fraction, critical threshold for global spread
- **Key dimensions:** global / transient / direct / network / state transition / memoryless or binary state / autonomous
- **Distinctness from P21:** One-directional spreading that consumes available fuel, not a self-stabilized partition.

#### P23 — Anti-coordination / emergent load balancing
Agents repeatedly self-distribute across options to avoid overcrowding, producing collective utilization near capacity without central assignment. The group neither converges to one shared choice nor fragments into stable camps — it dynamically balances.
- **Canonical model:** Minority Game (Challet & Zhang, 1997); El Farol Bar problem (Arthur, 1994)
- **Detection metric:** Mean attendance near capacity, reduced variance of attendance versus random-choice baseline, negative autocorrelation in option attendance, convergence to Nash equilibrium utilization
- **Key dimensions:** global / oscillatory or steady-state / direct or indirect / well-mixed or network / state transition / local state memory / autonomous
- **Distinctness from P18:** Group does not converge on one shared choice — it distributes across options.
- **Distinctness from P21:** Split is not stable camp formation — allocation is dynamically adaptive.
- **Distinctness from P32:** Allocation is continuously re-negotiated, not role-fixed.
- **Cognitive-analogue annotation:** Decision-making, adaptive resource allocation

### Cluster G: Resilience and Regulation

#### P24 — Homeostatic regulation
System maintains internal variable near a set-point despite sustained external perturbation, through active corrective feedback.
- **Canonical model:** Various adversarial perturbation protocols; thermostat-like agent models; bioelectric homeostasis models
- **Detection metric:** Recovery time after perturbation, steady-state deviation from set-point, deviation integral over time
- **Key dimensions:** global / steady-state / mixed interaction / state transition / negative feedback / local state memory / externally forced
- **Cognitive-analogue annotation:** Homeostasis, goal-maintenance

#### P25 — Canalized restoration (equifinality)
System converges to the same target macrostate from a wide range of initial conditions or after diverse perturbations, beyond what a simple dominant attractor would predict.
- **Canonical model:** Planarian regeneration models; Waddington landscape models; developmental GRN attractor basins
- **Detection metric:** Convergence variance across initial conditions, basin of attraction volume, restoration accuracy after disruption at different points
- **Key dimensions:** global / steady-state attractor / mixed interaction / state transition / local state memory / autonomous
- **Distinctness from P24:** P24 is cybernetic (continuous regulation). P25 is geometric (convergence from arbitrary starting points). Active vs passive resilience.
- **Cognitive-analogue annotation:** Goal-directedness, robustness

#### P26 — Stochastic resonance (computational-performance behavior)
Noise improves system performance at an intermediate optimal level — too little noise yields no signal detection, too much noise destroys the signal, but a sweet spot amplifies it.
- **Canonical model:** Bistable systems with subthreshold periodic signal plus noise; Collins et al. sensory neuron models
- **Detection metric:** Performance vs noise curve showing inverted-U shape, signal-to-noise ratio peak at nonzero noise
- **Key dimensions:** global / steady-state / local interaction / various substrates / state transition / periodically driven
- **Note:** Computational-performance behavior rather than a spatial, temporal, or social pattern. Retained due to crisp detection signature and large literature.

### Cluster H: Competition and Cooperation

#### P27 — Spatial reciprocity / emergent cooperation
Locally interacting selfish agents generate persistent islands of cooperation that survive exploitation, even when the payoff structure predicts cooperation should go extinct in a well-mixed population.
- **Canonical model:** Nowak & May (1992) spatial Prisoner's Dilemma on a lattice
- **Detection metric:** Cooperator survival fraction above well-mixed baseline, positive spatial autocorrelation of cooperative strategies, cluster survival time, interface dynamics
- **Key dimensions:** mesoscale / steady-state / direct / lattice / reproduction-selection / mixed-motive / autonomous
- **Distinctness from P1:** No physical movement — clusters emerge through differential survival on a static lattice.
- **Distinctness from P18:** Stable social regime under conflict, not convergence to shared choice.

#### P28 — Wealth condensation / spontaneous inequality
Symmetrical, fair exchanges of a conserved resource between identical agents inevitably lead to highly skewed, winner-take-all distributions.
- **Canonical model:** Yard-Sale model (Boghosian, 2017); Bouchaud-Mézard model
- **Detection metric:** Gini coefficient trajectory toward 1.0, emergence of Pareto power-law tail, oligarchic fraction holding majority of resources
- **Key dimensions:** global / steady-state (absorbing) / direct / well-mixed or network / resource exchange / positive feedback / local state memory / autonomous
- **Distinctness from P32:** Resource accumulation, not functional role differentiation.
- **Distinctness from P14:** Not a critical state with scale-free avalanches but convergence toward an absorbing boundary of extreme inequality.

### Cluster I: Emergent Structure Formation

#### P29 — Trail / network formation
Efficient transport networks emerge from collective movement and environmental modification without global planning.
- **Canonical model:** Ant trail optimization; Physarum polycephalum network models; active walker models
- **Detection metric:** Network efficiency relative to minimum spanning tree, Steiner ratio, fault tolerance, comparison to optimal network
- **Key dimensions:** mesoscale / steady-state / indirect-stigmergic / continuous space / motion-transport / environmental trace memory / positive feedback / autonomous

#### P30 — Spontaneous boundary formation (autopoiesis)
Agents spontaneously organize into a closed, semi-permeable topology that maintains an internal micro-environment distinct from the exterior.
- **Canonical model:** Varela, Maturana & Uribe SCL model; computational lipid bilayer simulations
- **Detection metric:** Topological closure (Euler characteristic), persistent concentration gradient between interior and exterior, self-repair after breach
- **Key dimensions:** mesoscale / steady-state / local interaction / continuous space / heterogeneous / state transition / autonomous
- **Distinctness from P1:** Not mere clustering — requires topological closure and functional inside/outside distinction.
- **Cognitive-analogue annotation:** Self/other distinction, individuality

### Cluster J: Emergent Agent-Level Competencies

#### P31 — Delayed gratification (provisional — pending non-redundancy test)
A systematic population-level pattern in which a statistically significant fraction of agents accept short-term cost (moving away from local goal) to achieve long-term benefit (global optimum), without any explicit planning mechanism.
- **Canonical model:** Zhang et al. cell-view sorting (2024)
- **Detection metric:** DG index (fraction of steps where agent moves away from eventual final position), population distribution of DG index values, correlation between population-level DG index and global sorting efficiency, comparison to null model (random walk) DG rates
- **Key dimensions:** local-to-mesoscale / transient / direct / lattice / motion-transport / memoryless / autonomous
- **Provisional status:** Retained pending empirical non-redundancy test. Protocol: (1) baseline model using only aggregation features (Moran's I, cluster count, cluster size distribution, path length, final segregation index, endpoint quality); (2) extended model adding DG features (mean DG index, DG variance, quantiles, temporal concentration of DG events, DG-event clustering across agents); (3) ablation model where DG temporal ordering is shuffled but marginal movement statistics are preserved. P31 survives only if DG features add out-of-sample predictive power after controlling for baseline, and if that gain disappears when DG timing structure is destroyed.
- **Cognitive-analogue annotation:** Gratification deferral, implicit planning

#### P32 — Emergent specialization (division of labor)
Identical agents spontaneously differentiate into distinct functional roles that increase collective efficiency, without pre-assigned task allocation.
- **Canonical model:** Response threshold models (Bonabeau et al., 1996); foraging task allocation in ant colonies
- **Detection metric:** Behavioral entropy reduction per agent over time, efficiency gain vs non-specialized baseline, task-switching frequency
- **Key dimensions:** mesoscale / steady-state / direct or indirect / various substrates / state transition / local state memory / autonomous
- **Merges former entries:** G1 (role differentiation) + N9 (division of labor)
- **Cognitive-analogue annotation:** Specialization, adaptive task allocation

---

## Layer 2 — Cross-Cutting Descriptors

### 2A: Mathematical Descriptors

- **Phase transition / symmetry breaking** — Applies to: P5, P8, P9, P20, P21
- **Critical regime** — Applies to: P14, P15
- **Positive feedback amplification** — Applies to: P1, P2, P8, P22, P28
- **Attractor dynamics** — Applies to: P16, P25
- **Nontransitive interaction** — Applies to: P12

### 2B: Cognitive-Analogue Descriptors

Observer-dependent interpretive lenses from the Diverse Intelligence framework.
Future versions: formalize as predicate sheets with plain-language definition,
minimum evidence standard, disqualifying cases, exemplar rows, confidence level.

- **Memory** — Annotates: P4, P16
- **Goal-directedness** — Annotates: P24, P25
- **Decision-making** — Annotates: P18, P23
- **Perception / inference** — Annotates: P17, P26
- **Gratification deferral** — Annotates: P31
- **Collective intelligence** — Annotates: P17, P18, P32
- **Self/other distinction** — Annotates: P30

---

## Layer 3 — Meta-Capacities

- **Stigmergic coordination** — Interaction modality enabling P4, P20, P29, P32.
- **Multiscale competency** — Central concept in Levin's TAME framework. Organizing principle for why atomic patterns compose into collective intelligence at higher scales.
- **Emergent goal shift** — System-level objectives diverge from individual agent goals. Pending operationalization.
- **Productive forgetting** — Information loss improving system performance. Manifests across P29, P32.

---

## Summary Statistics

- **Atomic patterns (Layer 1):** 32
- **Clusters:** 10
- **Mathematical descriptors (Layer 2A):** 5
- **Cognitive-analogue descriptors (Layer 2B):** 7
- **Meta-capacities (Layer 3):** 4
- **Ontological dimensions:** 11

---

## Deferred Candidates

- **Competitive coarsening / Ostwald ripening** — t^(1/3) scaling of aggregate-level dynamics. Literature treats as kinetic regime of clustering. One reviewer recommended promotion; two recommended deferral.
- **Swarm morphogenesis / programmable self-assembly** — Agents forming non-periodic bounded shapes (Rubenstein et al. 2014 Kilobots). Strong but canonical demonstrations require target specifications. Top candidate if target-free shape formation gets cleaner canonical model.
- **Active nematics / topological defects** — Dense elongated agents forming +/-1/2 defects acting as mobile super-agents. Major active matter focus but may be too substrate-specific.
- **Selfish-routing inefficiency / Braess paradox** — Locally optimal path selection producing globally inferior equilibria. Canonical agents often one-shot path selectors, at edge of scope.
- **Entropic segregation** — Macro-signature similar to P1; mechanistic variant under behavior-first rules.
- **Brazil Nut Effect / kinematic sorting** — Rare dimensional region but potentially too domain-specific.

---

## Open Questions

1. **P31 non-redundancy test.** See three-stage protocol in P31 entry. Primary v0.5 validation target.
2. **P13/P15 boundary.** Discriminating metric: Transfer Entropy across collisions. Validate empirically.
3. **Cognitive-analogue formalization.** Next step: predicate sheets for each 2B annotation.
4. **Detector specification cards.** Primary v0.5 deliverable. Each entry needs: required observables, metrics, null model, pass threshold, false positives, neighbor exclusions.

---

## Changelog

### v0.4 (current)
- Added P23 (anti-coordination / emergent load balancing)
- Formalized cognitive analogues as Layer 2B annotation band
- Added Layer 2/Layer 3 boundary rule
- Added update mode and external driving ontological dimensions
- Tightened P15 (requires information-bearing signal interactions)
- Tightened P13/P15 distinctness via Transfer Entropy criterion
- Added provisional status and three-stage non-redundancy test for P31
- Added swarm morphogenesis and active nematics to deferred candidates
- Fixed count mismatches and cross-references
- Standardized dimensional profiles as key dimensions
- Added cognitive-analogue annotations to relevant entries

### v0.3
- Added P13 (excitable spiral/target waves) and P19 (emergent leadership)
- Reframed edge-of-chaos as P15 (persistent propagating computation)
- Added interaction substrate ontological dimension
- Flagged cognitive analogues as open architectural question

### v0.2
- Expanded from 9 to 28 patterns via systematic literature review
- Introduced three-layer architecture
- Split P1 into similarity-driven and MIPS
- Merged role differentiation and division of labor
- Demoted phase transitions to Layer 2, stigmergy to Layer 3

### v0.1
- Initial catalog of 9 confirmed patterns (P1-P9)
- 5 gap analysis predictions
- Basic ontology with 8 dimensions
