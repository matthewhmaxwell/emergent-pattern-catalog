# Section 1: Introduction

## 1.1 Motivation

Minimal agent-based models produce macro-level behaviors that no
individual agent is programmed to perform. Flocks align without a
conductor, traffic jams propagate backward against traffic flow,
cooperators survive in payoff landscapes that predict their extinction,
opinions polarize into durable camps rather than converging on consensus.
These emergent patterns are the empirical substance of the complex-systems
literature, and they have been studied for decades across ecology,
physics, economics, sociology, and biology.

Two recent developments renew the scientific importance of cataloging
these patterns systematically. First, the Diverse Intelligence research
program (Levin, 2019, 2022) argues that cognition is not a binary property
of neural systems but a continuum of problem-solving competencies
exhibited by biological systems at every scale — from morphogenetic
cells navigating bioelectric landscapes to bacterial populations
coordinating antibiotic resistance to social insect colonies solving
optimization problems. A concrete step toward operationalizing this claim
is to build an inventory of collective competencies whose signatures can
be measured across substrates. Second, a growing class of experiments in
synthetic morphology, swarm robotics, and computational embryology
reports behaviors that look superficially like planning, memory, or
preference in systems whose individual components have none of these
capabilities. Interpreting such findings requires shared tools for
distinguishing a genuinely novel competency from a previously documented
pattern with a different name.

The motivating example for this work is Zhang, Goldstein & Levin (2024),
who reimplemented the classical Bubble, Insertion, and Selection sort
algorithms as cell-autonomous agents on a shared one-dimensional array.
The striking finding was not that decentralized sorting succeeds — that
was the expected result — but that the trajectories exhibit *delayed
gratification*: a statistically significant fraction of steps move agents
away from their eventual final positions, in exchange for global progress
toward a sorted state. The authors further demonstrated that in chimeric
arrays mixing algorithm types, agents spontaneously segregate by
algorithm, producing spatial aggregation without any similarity-preference
rule. Both observations pose a diagnostic challenge: what pattern, in
what catalog, has been observed here? The aggregation signature resembles
Schelling-style segregation; the delayed-gratification signature resembles
credit-assignment in reinforcement learning; the overall trajectory
resembles quenched annealing. Without a shared catalog with quantitative
detectors, the temptation is to invent new terminology for each
observation, or to assimilate findings into whichever nearest-neighbor
pattern the observer happens to know best.

## 1.2 The Gap

Reviews of the emergent-phenomena literature have enumerated patterns
for over three decades. Camazine et al. (2001) catalogued self-organized
biological systems; Ball (2009) surveyed pattern formation across
physics, chemistry and biology; Vicsek & Zafeiris (2012) reviewed
collective motion; Castellano, Fortunato & Loreto (2009) reviewed
statistical physics of social dynamics. These surveys are valuable
reference works, but they share three limitations that constrain their
use as diagnostic tools.

*Taxonomic rather than operational.* Most existing surveys describe
patterns at the level of mechanisms and canonical examples rather than
quantitative detectors. A reader can learn that Schelling segregation
produces clustering and that Nowak-May produces spatial reciprocity, but
no shared metric — applicable across both — is specified that would
return a structured result on an arbitrary lattice model. Without such a
metric, cross-substrate claims ("this biological system exhibits the same
competency as Schelling") remain qualitative.

*No cross-cutting dimensional framework.* Reviews typically organize
patterns by disciplinary origin (physics, biology, sociology) or
mechanism family (diffusion, threshold, game-theoretic), which makes it
difficult to identify which ontological dimensions distinguish patterns
from each other or to predict where an underspecified pattern should
live. The result is a collection of well-described instances without a
predictive scaffolding that can flag gaps or catch redundancies.

*No null-model discipline.* Published demonstrations of a pattern
typically report positive findings only. Whether the same model produces
the same detection signature under negative controls — randomized
interactions, shuffled topology, ablated mechanisms — is frequently not
tested, and across-pattern exclusion ("this system exhibits P5 but not
P6") is almost never reported. Without null models and exclusion tests,
confident cross-pattern claims are unwarranted.

These three gaps are independent, and closing them requires different
work: taxonomy requires catalog construction, dimensional analysis
requires an ontology, and null-model discipline requires a detection
toolkit with explicit statistical protocol. This paper addresses all
three.

## 1.3 Contribution

We present three interlocking artifacts.

*A three-layer pattern catalog of 32 atomic emergent behaviors.* Layer 1
contains 32 atomic patterns across 10 clusters, each defined by an
observable macro-behavior, a canonical minimal model, and at least one
quantitative detection metric. Layer 2 adds cross-cutting descriptors —
mathematical (phase transitions, critical regimes, attractor dynamics)
and cognitive-analogue (memory, goal-seeking, preference) — that apply
across rows as annotations rather than row-defining criteria. Layer 3
names meta-capacities (stigmergic coordination, multiscale competency)
that describe *how* atomic patterns are generated and composed, bridging
the catalog to the Diverse Intelligence framework without embedding
interpretive claims in the behavioral ontology. The design rules for
Layer 1 inclusion and the boundary between Layers 2 and 3 are stated
explicitly in Section 2.

*An ontology of 11 dimensions for dimensional classification.* Each
pattern is situated along axes covering spatial scale, temporal
character, interaction type, substrate, agent homogeneity, goal
structure, feedback, memory, conflict, external driving, and update
mode. Dimensional coverage analysis identifies underexplored regions of
the space — candidates for future discovery — and the ontology underpins
the cross-exclusion graph used by the detection toolkit.

*A detection toolkit of 13 quantitative detectors with explicit null
models, tier-gated significance criteria, and exclusion logic.* Every
detector returns a structured `DetectorResult` containing primary and
secondary metrics, a null p-value under a specified null model, a
detection tier (screening, confirmation, or definitive), confidence
capped by tier, and the outcomes of nearest-neighbor exclusion tests.
Detectors declare their substrate compatibility and their observable
scope (whether they operate on state history alone or require model
metadata), enabling systematic dispatch across the model inventory.

We validate the toolkit by implementing 20 registered models across
19 distinct canonical model families (Zhang cell-view sorting has
sequential and threaded variants of the same family that share
substrate, observables, and primary patterns and are folded as a
single display row in the consolidated transfer matrix) spanning
nine pattern clusters and seven substrates (lattice_1d, lattice_2d,
lattice_2d_continuous, continuous_2d, oscillator, opinion_space,
and scalar_wealth). Each model replicates published quantitative
results and each of 19 implemented detectors produces the correct
tier assignment on both positive and negative controls. The
resulting transfer matrix contains 79 substrate-and-observable
compatible cells out of 380 registry cells (19 displayed rows × 19
detector columns yields 361 cells in the folded display); the
remaining 301 registry cells (284 in the displayed table) are
correctly eliminated without empirical testing — 274 by substrate
mismatch and 27 by detector–observable incompatibility. Of the 79
audited cells, 19 produce canonical DEFINITIVE detections (one per
primary model family), pinned in dedicated end-to-end test files;
the other 195 audited model–detector outcomes — covering screening
admits, screening rejections, and confirmation-tier co-occurrences
— are pinned in the cross-detection regression suite. Every
canonical positive reaches at least confirmation tier on its
primary detector and every cross-pattern negative control is
correctly rejected by substrate match, observable match,
prerequisite guard, or primary-metric screening. The transfer
matrix figures cited here are programmatically derived from the
model and detector registries by `scripts/count_transfer_matrix.py`
and pinned by `tests/test_transfer_matrix_counts.py` so that any
future registry change that drifts from the cited numbers fails the
test suite before propagating silently into the manuscript.

## 1.4 Scope and Limitations

Three scope boundaries bound the claims of this paper.

*Atomic patterns, not composite ones.* Layer 1 entries are single
observable macro-behaviors with canonical minimal models. Composite
phenomena that combine multiple atomic patterns — morphogenetic programs
building target shapes, neural-network-like recurrent computation —
are not Layer 1 entries and are not detected directly. The catalog
supports analysis of such phenomena through co-occurrence of atomic
signatures, but it does not provide composite detectors.

*Detection, not explanation.* The toolkit answers "does pattern X appear
in system Y?" quantitatively; it does not answer "why does pattern X
appear in system Y?" Mechanistic explanations are returned to the
modeler. This is a deliberate scoping choice: detection that does not
depend on knowing the generative mechanism is precisely what is needed
when the long-term goal is to apply these methods to systems whose rules
are unknown.

*Twenty model families, not an exhaustive census.* The transfer matrix
at the time of writing covers 20 registered models (19 distinct model
families, with two Zhang sorting variants folded as one) across the
seven substrate types lattice_1d (Zhang sorting; Nagel-Schreckenberg
traffic), lattice_2d (nine cellular automata: Greenberg-Hastings, Game
of Life, Schelling, Nowak-May, BTW sandpile, SIR epidemic, spatial
rock-paper-scissors, lattice Lotka-Volterra, voter model),
lattice_2d_continuous (Gray-Scott reaction-diffusion), continuous_2d
(Vicsek, D'Orsogna, active Brownian particles), oscillator (all-to-all
Kuramoto, non-local Kuramoto), opinion_space (Hegselmann-Krause), and
scalar_wealth (Yard-Sale wealth exchange). Several documented patterns —
lane formation (P7), territoriality (P4), associative memory (P16),
distributed sensing (P17), leadership (P19), quorum sensing (P20),
anti-coordination (P23), homeostatic regulation (P24), canalized
restoration (P25), stochastic resonance (P26), trail formation (P29),
autopoiesis (P30), and emergent specialization (P32) — do not yet have
canonical-model implementations in our toolkit. We treat the current
matrix as a working core, not a closed table, and identify
dimensional-coverage gaps in Section 5.

## 1.5 Organization

Section 2 presents the pattern catalog: its three-layer architecture,
the design rules for Layer 1 inclusion, the 11-dimensional ontology,
and a survey of the 32 atomic patterns organized by cluster. Section 3
specifies the detection toolkit: design principles, detector
architecture, null-model taxonomy, substrate-aware dispatch, key
boundary tests, and statistical power requirements. Section 4 reports
replication studies for the 19 canonical model families, verifying
each against published results and reporting detector tier assignments
on positive and negative controls. Section 5 presents the consolidated
transfer matrix over 214 audited cells (195 in the cross-detection
regression suite plus 19 canonical positives pinned in dedicated
end-to-end files), analyzes its block-diagonal structure by substrate,
and examines four cross-model findings — P1 co-occurrence with P27 on
Nowak-May, the asymmetric P1 signature on SIR versus RPS, the
bilateral-versus-cyclic exclusion between P11 and P12, and the
same-substrate content-level discrimination of voter coarsening from
the four other lattice_2d-with-grid patterns — that sharpen pattern
definitions beyond their initial specifications. Section 6 reports
findings that were not anticipated at the outset: the P31
non-redundancy result, the boundary-conditioned transfer entropy
technique, and resolved false positives that refined operational
definitions. Section 7 discusses implications for the Diverse
Intelligence framework, the periodic-table metaphor, methodological
contributions, limitations, and future work. Section 8 concludes.

All source code, parameter files, detector specifications, and test
suites are available at [repository URL]. The complete inventory
includes 20 model implementations (19 distinct model families),
19 detectors, 214 audited model × detector cells (195 in the
cross-detection regression suite and 19 canonical positives in
dedicated e2e files), and a fast-test suite of approximately 478
canonical regression tests (plus several slow-marked tests for
higher-power replications and finite-size studies).
