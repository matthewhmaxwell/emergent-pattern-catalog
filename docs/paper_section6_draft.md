# Section 6: Emergent Findings

The systematic application of quantitative detectors across the model
inventory produced several findings that were not anticipated at the
outset. These fall into four broad categories: validation of a
contested pattern (P31's non-redundancy with P1); a novel
methodological technique (boundary-conditioned transfer entropy);
unexpected cross-model results that sharpened existing pattern
definitions (GoL × P1, GoL × P13, Nowak-May × P1, SIR × P1 vs RPS ×
P1, LV × P11 vs RPS × P12, BTW sandpile's 1/f signature); and a
recurring meta-finding that reshaped how new detectors are designed —
the pattern that the catalog's initial "obvious" statistical recipes
frequently fail empirically on the actual substrate, forcing a pivot
to mechanism-derived metrics and a corresponding architectural
investment in substrate-content and metadata-mechanism discrimination
classes. This section documents the most load-bearing findings in
each category.

## 6.1 P31 Non-Redundancy Validation

P31 (delayed gratification) was initially flagged as "provisional" in the
pattern catalog because it was unclear whether the DG signal in Zhang's
sorting model captured information not already present in P1 (aggregation).
If DG were merely a temporal signature of the same spatial aggregation
process, it would not merit a separate pattern entry.

The non-redundancy test addresses this directly. A baseline regression model
uses only P1-type features (Moran's I, cluster count, cluster size
distribution, path length, final segregation index, endpoint quality) to
predict sorting efficiency across 600 chimeric sorting runs. An extended
model adds DG features (mean DG index, DG variance, quantiles, temporal
concentration of DG events, DG-event clustering across agents). If DG
features add no information beyond aggregation, the extended model should
show no improvement.

The result is decisive: the baseline model achieves R² = −0.02 (P1 features
alone cannot distinguish algorithms — all three produce identically sorted
endpoints), while the extended model achieves R² = +0.63. The improvement
ΔR² = +0.645 (p < 0.000001, 10-fold cross-validation) establishes that DG
captures temporal structure in the sorting trajectory that is orthogonal to
the spatial aggregation endpoint. Ablation — shuffling DG features while
preserving marginal statistics — destroys the signal (R² = −0.03),
confirming that temporal ordering is the information-bearing component.

P31 is therefore elevated from provisional to confirmed as a distinct
pattern. The non-redundancy protocol itself is generalizable: any candidate
pattern suspected of redundancy with a neighboring pattern can be tested by
comparing predictive power of baseline (neighbor features only) vs. extended
(baseline + candidate features) models on a discriminative task.

## 6.2 Boundary-Conditioned Transfer Entropy

The most significant methodological finding is the boundary-conditioned
transfer entropy technique developed for the P13/P15 discriminator. This
emerged from a failed attempt to use raw average TE to distinguish excitable
wave propagation from computational dynamics on lattices.

**The problem.** Both Greenberg-Hastings (excitable waves) and Game of Life
(computation) produce spatial dynamics on 2D grids. A natural hypothesis is
that computation involves greater information transfer between cells. Raw
average TE — computed over all cells — gives the opposite result: GH
produces higher average TE than GoL. This occurs because deterministic wave
propagation creates trivially high TE at interior cells (a cell's next state
is perfectly predicted by its excited neighbor's current state), dominating
the average.

**The solution.** Restricting TE measurement to boundary cells — cells
adjacent to the transition between active and inactive regions — isolates
the signal that matters. At boundaries, excitable waves produce minimal TE
(the wave simply passes through) while computational dynamics produce high
TE (glider collisions, structure interactions). The boundary-conditioned
TE ratio cleanly separates the two systems: GoL produces ratios of 15–16×
above a permutation null, while GH produces ratios of 1–2×.

**Significance.** This technique is applicable to any system where a bulk
average obscures boundary phenomena. The general principle — that the
information-theoretic signature of computation is concentrated at structural
interaction boundaries rather than in propagation interiors — may extend to
biological systems where the distinction between stereotyped signal relay
and adaptive information processing is diagnostically important.

## 6.3 False Positive Analysis

Two resolved false positives produced insights that sharpened pattern
definitions beyond their original specifications.

**GoL × P1: what counts as aggregation.** The P1 detector initially
produced a confirmation-tier detection on Game of Life because B3/S23
dynamics generate genuine spatial autocorrelation — alive cells cluster
due to the survival rule's neighbor requirements. The resolution required
defining what P1 actually detects more precisely: not spatial
autocorrelation from any source, but the spatial aggregation of *persistent
agent type labels* through *relocation dynamics*. This led to two guards:
a type constancy check (values at each position must be stable over the
measurement window, which GoL's flickering states fail) and a temporal
convergence check (Moran's I must show monotonic increase toward a plateau,
the signature of agents actively relocating to preferred neighborhoods).
The false positive forced a conceptual refinement that makes P1 more
precisely defined and less prone to spurious detections on systems with
incidental spatial structure.

**GoL × P13: what counts as excitable media.** Game of Life passed P13's
screening threshold because synchronous binary updates produce wavefront
speed CV = 0.0 — every cell has identical inter-excitation intervals.
The resolution was straightforward (requiring n_states ≥ 3 for excitable
dynamics), but the lesson is general: metrics designed for one class of
system can produce artifactually clean statistics on structurally different
systems. Binary automata are not merely "low-state excitable media" — they
lack the refractory mechanism entirely. The guard enforces the structural
prerequisite, not just a threshold on the metric.

## 6.4 Cross-Pattern Co-Occurrence

The detection of P1 (aggregation) at confirmation tier on Nowak-May
was the most interesting unexpected cross-model result. Cooperator clusters
in the spatial Prisoner's Dilemma (Moran's I = 0.898, segregation index =
0.753, p = 0.005) exhibit the same spatial autocorrelation signature as
Schelling segregation (Moran's I = 0.423, segregation index = 0.652,
p = 0.001), despite arising from entirely different mechanisms.

In Schelling's model, agents physically relocate to achieve neighborhood
homogeneity — a preference-driven movement process. In Nowak-May's model,
no agent moves; instead, agents change strategy by imitating their
highest-payoff neighbor. Cooperator clusters form because mutual cooperation
produces higher payoffs than exploitation at boundaries, causing cooperator
strategy to spread inward from cluster cores while defector strategy invades
at edges.

The co-occurrence of P1 and P27 on Nowak-May suggests that aggregation is
a broader phenomenon than originally conceived. The pattern catalog defines
P1 in terms of spatial autocorrelation and type segregation, not in terms
of the mechanism (movement vs. imitation) that produces it. This is
deliberate: the catalog's patterns are defined by observable signatures,
not generative mechanisms, precisely because the long-term goal is to detect
patterns in systems where the mechanism may be unknown.

The difference in tier (Nowak-May reaches confirmation, not definitive)
correctly captures the mechanistic distinction: imitation-based clustering
lacks the convergence dynamics (monotonic I increase toward plateau) that
characterize movement-based aggregation. The tier system thus provides
a natural vocabulary for expressing "same pattern, different mechanism"
without either rejecting the detection or over-claiming equivalence.

## 6.5 The 1/f Noise Correction

The BTW sandpile's 1/f noise signature was initially unmeasurable
(spectral exponent β = −0.17), casting doubt on a key published claim.
Investigation revealed that the spectral analysis was applied to the
avalanche size sequence — which is approximately IID, with each event
being largely independent of the previous one. The temporal correlations
in SOC systems appear in the *total system energy* E(t) = Σh(x,t), which
integrates over the spatial state rather than recording individual events.
Measuring PSD of E(t) yielded β = 1.41, correctly in the 1/f range.

This error would have been difficult to catch without the quantitative
replication standard: not merely observing power-law avalanches (which were
correct), but verifying all published claims including the 1/f temporal
signature. The correction also illustrates a general principle: the
observable that produces a clean metric may not be the one that first comes
to mind. Careful matching between the theoretical prediction and the
measured quantity is essential.

## 6.6 The SIR versus RPS Asymmetry on P1

Both SIR (spatial epidemic) and spatial rock-paper-scissors produce
strong peak Moran's I during their dynamics — SIR at 0.89 during
wavefront propagation, RPS at 0.55+ continuously during spiral
evolution. An initial version of the P1 detector used peak Moran's I
as its primary metric and fired at screening tier on both. This
conflated two qualitatively distinct patterns: SIR's clustering is
*transient* — a wavefront that vanishes once the epidemic recovers —
while RPS's clustering is *sustained* — spiral domains rotate
indefinitely but never dissipate.

The Sprint 10 characterization across six models revealed the
distinguishing observable: final-state Moran's I. SIR drops to I_final
≈ 0.02 while RPS maintains I_final ≈ 0.55. Changing the primary
metric from peak to final correctly flips SIR × P1 from screening to
rejected while leaving Schelling, Nowak-May, and RPS detections
intact.

The finding illustrates a general phenomenon relevant to any
signature-based catalog: the same peak value can correspond to
qualitatively different patterns, and distinguishing them requires
observables that capture temporal persistence. A peak-sensitive metric
answers "has the system ever clustered?" — which SIR does during its
wavefront. A final-state metric answers "does the system remain
clustered?" — which SIR does not. The two questions have different
answers across transient and sustained dynamics, and only the second
is the question P1 is designed to ask. This is a broader instance of
the same conceptual issue that §6.5's BTW investigation surfaced: the
*observable* rather than the *metric* needs to be matched to the
theoretical prediction.

## 6.7 Pattern-Catalog-Obvious Recipes That Failed Empirically

The most consequential meta-finding across Sprints 10, 13, 15, 16, 17,
and 18 was a recurring structural failure in the catalog's initial
detector specifications. For multiple patterns, the first "obvious"
statistical recipe — the one that appeared in the literature's
textbook characterization of the pattern — failed empirically when
applied to the canonical model. The failures share a common structure:
the obvious metric produced numerically indistinguishable values
between the target pattern and at least one well-defined negative
control, so no threshold on that metric could discriminate between
them.

Five instances, by chronological sprint order:

*Sprint 13, P3 Turing-wavelength, peak-to-mean radial-FFT metric.*
Spatial RPS at low mobility produces raw-grid peak-to-mean ≈ 23.10,
exceeding Gray-Scott labyrinths at 18.75. No empirical threshold can
separate these systems. Resolution (Decision 37): substrate-content
gate — P3 requires a `field` observable and n_unique_values ≥ 50.

*Sprint 15, P8 Traffic Jamming, stopped-car fraction.* At saturation
density ρ = 0.80 with p = 0, pigeonhole geometry forces stopped-
fraction = 0.750 — far above the screening threshold — with no
emergent stop-go dynamics. Resolution (Decision 41): secondary metric
— 95th-percentile consecutive v = 0 run length, which is 13 for
emergent jamming and 4 for density-saturation.

*Sprint 16, P2 MIPS, Hartigan dip test.* Particle-level local densities
are integer counts divided by a constant area, producing discrete
distributions that are trivially non-uniform by Hartigan's test
regardless of underlying physics. Dip p-values floored at the
bootstrap minimum across every tested regime, including known-uniform
dilute and known-one-phase stuck. Resolution (Decision 44): primary
metric — two-phase coexistence score min(f_gas, f_liquid).

*Sprint 17, P28 Wealth Condensation, Pareto α Hill estimator.* The
canonical Pareto range 1 < α < 2 is reached only in a narrow transient
window that shifts with stake fraction f; at long time α drops below
1 (degenerate Pareto) and eventually approaches 0 (δ-on-winner). No
fixed-α gate discriminates condensation from non-condensation.
Resolution (Decision 47): primary metric — Gini coefficient at the
final frame of the measurement window.

*Sprint 18, P10 Chimera States, per-window coexistence.* Ordinary
Kuramoto near K_c produces per-window local-r statistics numerically
indistinguishable from a chimera — both systems exhibit persistently
coherent and persistently incoherent windows, because the mean-field
model's internal ω-sort creates frequency-structured window
heterogeneity that mimics the chimera's position-structured
heterogeneity. Resolution (Decision 50): primary metric — spatial
autocorrelation of per-oscillator phase velocity at ring lag 4, which
separates position-organized (chimera) from frequency-organized
(Kuramoto) window structure.

Each resolution replaced a pattern-descriptive metric (one that asks
"does the system produce the pattern's characteristic *output*?") with
a mechanism-derived metric (one that asks "does the system produce
the pattern via its characteristic *mechanism*?"). This is not merely
a collection of ad-hoc detector fixes: it is the recurring empirical
finding that *specification-by-output* and *specification-by-mechanism*
are genuinely different detector architectures, and that mechanistic
specifications are more robust against cross-model confounders. The
"look before touching" discipline — always running broad empirical
characterization before locking detector thresholds — surfaced each
of these cases before they entered the production threshold set. A
catalog built on pattern-descriptive metrics alone would have
accumulated false positives silently.

## 6.8 The Three-Class Discrimination Framework

The recurring failures in §6.7 and the increasingly dense transfer
matrix together drove the emergence of a three-class discrimination
framework that became the catalog's default architecture for detector
specification (Decisions 37, 41, 43, 49, 52):

*Class 1: substrate-type discrimination (registry).* The orchestration
layer's MODEL_REGISTRY / DETECTOR_REGISTRY pair defines a block-
diagonal structure on the transfer matrix. Detectors register a
required substrate (`lattice_1d`, `lattice_2d`, `lattice_2d_continuous`,
`continuous_2d`, `oscillator`, `opinion_space`, `scalar_wealth`) and
models register theirs; cross-substrate detector application is
prevented automatically. No false positives reach the detector
pipeline across substrate boundaries.

*Class 2: substrate-content discrimination (observable values).* Within
a shared substrate, detectors declare required observable keys and,
where needed, content-level prerequisites on those values — P3's
`n_unique_values ≥ 50` gate, P8's 1D-integer-velocities gate, P11's
two-species prerequisite. These content-level gates separate models
that share substrate but differ in the kind of data they produce.
Decisions 37 (P3 continuous-field gate), 41 (P8 integer-velocity
gate), and 39 (P11 two-species prerequisite) are the canonical
instances.

*Class 3: metadata-mechanism discrimination (rule flags).* Within a
shared substrate and matching observable type, detectors query
`model_metadata` flags that assert the presence or absence of the
specific mechanism the detector is built to identify. Decision 43
(ABP × P2 with `has_density_dependent_speed = True`,
`has_alignment_rule = False`, `has_attraction_rule = False`),
Decision 49 (Yard-Sale × P28 with the four-flag wealth-dynamics gate),
and Decision 52 (Kuramoto-nonlocal × P10 with `has_nonlocal_coupling
= True`, `has_frequency_heterogeneity = False`) all employ this
pattern. The metadata gate blocks DEFINITIVE tier when the mechanism
flag is false even if empirical signatures superficially match.

The three classes are complementary rather than alternative. Most
detector registrations use only Class 1 (substrate is the primary
filter); detectors dealing with multi-observable substrates add
Class 2; detectors dealing with same-substrate-same-observable model
pairs add Class 3. P2, P28, and P10 each employ all three. The
framework is substrate-independent: the same metadata-mechanism
architecture works on continuous_2d (P2), scalar_wealth (P28), and
oscillator (P10), suggesting it generalizes to any future substrate
where multiple dynamically distinct models coexist.

## 6.9 Bilateral versus Cyclic: the LV/RPS Boundary

The sharpest within-substrate cross-pattern discrimination in the
catalog is between P11 (bilateral predator-prey oscillation,
canonical on Lotka-Volterra) and P12 (cyclic dominance, canonical on
spatial rock-paper-scissors). LV and RPS share the `lattice_2d`
substrate, share the `grid` observable type, and both produce
time-oscillating population fractions with strong spatial
autocorrelation. The pattern-catalog-obvious discriminator — "P11 uses
cross-correlation between two species; P12 uses intransitivity score
across three species" — fails to capture the structural boundary,
because each detector's numerical output on the wrong model is not
always well-separated from its output on the right model.

The architectural resolution is to move the discrimination to a
prerequisite on species identity, not on primary-metric magnitude.
P11 requires `n_unique_species_observed == 2` (exactly two non-zero
species in the grid value histogram); P12 requires `n_species ≥ 3`.
Models with a strict 2-species conservation law (Nowak-May's
cooperator / defector) trigger P11's prerequisite but fail the
`total_std` content-level guard (Decision 35: a nontrivial empty
reservoir is required — in a strictly conservation-locked 2-species
system, one species' fraction is the complement of the other's, so
anti-phase cross-correlation is mechanically guaranteed and carries
no information about predator-prey coupling). Models with 3+ species
(RPS) trigger P12 cleanly and fail P11's `n_species == 2` check.

The LV/RPS boundary is therefore *not* a boundary in detector-metric
space. It is a boundary in *species-count space*, enforced at
prerequisite-evaluation time before any primary metric is computed.
This is the Class 2 (substrate-content) discrimination framework of
§6.8 applied at its sharpest: the two patterns are structurally
mutually exclusive on their defining content-level property
(2-species vs 3+-species), and neither detector needs to worry about
the other's numerical behavior on the wrong model. The same
architectural strategy could be applied to any future pattern pair
that shares substrate and shares observable type but differs in
a structural property discoverable from state alone.
