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

## 6.10 Within-Substrate Discrimination Without Metadata: the Voter Case

Sprint 20 added the voter model and P18 coarsening-to-consensus to a
lattice_2d-with-grid block that already contained eight models (GH,
GoL, BTW, Schelling, Nowak-May, SIR, RPS, LV) and six detectors (P1,
P11, P12, P13, P15, P22). Voter shares substrate, observable type,
and even surface phenomenology (cluster formation, coarsening) with
several existing models. Yet P18 discriminates voter from each of the
nearby patterns without using a single metadata flag — the
discrimination is entirely metric-based at content level, in
contrast to Decisions 43 (P2), 49 (P28), and 52 (P10) which all
required a metadata-mechanism gate (Class 3 of §6.8) to separate
same-substrate same-observable models.

The voter discriminator architecture relies on three observations
about how the canonical voter trajectory differs from the four
nearest lattice_2d-with-grid patterns:

*Versus P13 (excitable waves) on GH spirals.* GH with a broken-wave
initial condition produces a Moran's I trajectory that is
*stationary* at ≈ 0.87 from the moment the wave stabilizes. Spearman
ρ over time is exactly zero for a constant signal. Voter, by
contrast, produces a Moran's I trajectory that *grows monotonically*
from ≈ 0 to ≈ 0.5 over the first τ ~ L / 2 sweeps. The screening
gate `moran_spearman_early > 0.70` rejects GH spiral immediately
because its early-window Spearman is 0; voter passes.

*Versus P13 on GH random transient.* GH with a random initial
condition is the harder case. Sparse surviving excited cells form
local clusters as excitation propagates and dies out, producing a
genuine early-time Moran growth (Spearman ≈ +0.93) and a genuine
early-time wall decay. GH random reaches the confirmation tier of
P18 — the screening and confirmation gates do not separate it from
voter. The discrimination instead happens at the definitive tier:
voter's wall density plateaus at ≈ 0.21 (significantly above zero
because it never reaches consensus on a finite torus in the
characterization window), while GH random's wall density collapses
to ≈ 0.011 once excitation extinguishes. The definitive gate
`wall_final_qtr_mean > 0.05` excludes GH random from DEFINITIVE
while still allowing it to be flagged as having "some coarsening
signal" at confirmation. This is a deliberate use of the
three-tier framework: the lower tiers correctly indicate that GH
random does undergo cluster formation in some sense, while the
definitive tier preserves the specific identity of "voter-like
coarsening to a balanced two-opinion plateau".

*Versus P15 (persistent computation) on GoL.* Game of Life with a
random initial condition decays to a sparse landscape of still-
lifes and oscillators. The Moran's I of the surviving alive cells
does grow in the early window (Spearman ≈ +0.87) but the plateau
saturates at ≈ 0.27, below the 0.30 screening floor of P18. GoL with
an r-pentomino seed starts at high Moran (≈ 0.35 due to the
compactness of the initial pattern) and shows essentially no
early-time growth. Both fail the screening gates without requiring
any metadata-mechanism gate.

*Versus P1 (similarity aggregation) on Schelling.* Schelling
segregation does coarsen, in a superficial sense — agents move
toward like-typed clusters. But the agents do not *flip* their
type labels; the coarsening is movement-driven, not copying-driven.
The voter dynamics' identifying mechanism is opinion *copying* with
local imitation, and this is exposed via the model metadata key
`update = 'asynchronous_copy_neighbor'`. The P18 detector's P1
exclusion uses this metadata key as a check, not as a hard gate
(the metric-level gates already exclude Schelling at screening
because Schelling's wall-density decay is geometry-bounded rather
than copying-driven, and its `moran_growth` over the early window
is below P18's 0.20 floor). The metadata key would only become
load-bearing if a future model implemented a copying-driven
mechanism that also satisfied the geometry-bounded metric profile
— i.e., it is a *defense in depth* rather than the primary
discriminator.

The voter case demonstrates that Class 3 (metadata-mechanism)
discrimination is not strictly required even for same-substrate
same-observable model pairs, provided the metric-level signatures
are sharply enough separated. Decisions 54-56 record the design
choices that achieved this metric-level sharpness:

  - **Decision 54 (Sprint 20).** P18 uses a full random permutation
    null on the Moran's I trajectory rather than a circular-shift
    null. The circular-shift null was the natural starting choice
    because it preserves time-series autocorrelation, which under H0
    of "no monotonic trend" should be preserved. But Moran's I in
    the voter dynamics is so strongly autocorrelated (consecutive
    values differ by less than 0.05 typically) that circular shifts
    fail to spread the null distribution adequately: the null
    Spearman ρ retains substantial mass at large positive values,
    inflating p-values above the 0.01 confirmation gate even on
    canonical positive runs. Replacing the circular-shift null with
    a full random permutation destroys both the trend and the
    autocorrelation, and the null Spearman distribution becomes
    centered near zero with std ≈ 0.16. Under this null all 10
    voter seeds tested at L = 64 produce p < 0.01, and all four
    discriminator scenarios produce p ≈ 1. This decision parallels
    Sprint 11 ADR 36 (circular-shift autocorrelation preservation
    in the P11 cross-correlation null) in a different detector
    context. The lesson is that null-model design must consider
    not just what the null preserves but what it must *destroy* to
    achieve adequate spread under the test statistic of interest.

  - **Decision 55 (Sprint 20).** P18's secondary metric is the
    wall-density Spearman ρ computed over t ≤ τ, not over the
    full trajectory. The voter wall-density trajectory has two
    qualitatively distinct regimes: a sharp decay over t ∈ [0, τ]
    where the wall density drops from ≈ 0.5 to ≈ 0.27, and a slow
    random-walk drift afterward at the plateau level. A
    full-window Spearman is dominated by the larger number of
    samples in the late regime, where the random-walk noise can
    push the empirical Spearman positive on individual seeds —
    causing false-negative confirmation-tier rejections. The
    early-window Spearman is robustly ≤ −0.83 across all seeds at
    every L tested (64, 128, 256). The lesson is that
    double-regime trajectories require careful window-restriction
    in the secondary metric, with the window aligned to the
    coarsening-active regime rather than to the full run.

  - **Decision 56 (Sprint 20).** All Sprint 20 voter
    characterization uses the canonical asynchronous Glauber-like
    voter dynamics: one Monte Carlo sweep is N elementary site-
    updates, each picking a uniformly random site and copying a
    uniformly random neighbor. A checkerboard parallelization of
    the dynamics was prototyped during development for ≈ 4×
    speedup at L = 256; it preserves the late-time coarsening
    exponent within statistical noise but produces early-time
    wall-density trajectories that differ from canonical async by
    more than 3σ at t = 10 sweeps. Because all P18 detector gates
    are calibrated against the early-time canonical trajectories
    (via Decisions 54 and 55), the speedup did not justify the
    quantitative drift, and the canonical async dynamics is the
    only dynamics used in characterization, in detector
    calibration, and in the slow-test pinning. The lesson is that
    parallel-update approximations to inherently asynchronous
    dynamics must be validated for trajectory-level equivalence,
    not only for asymptotic-exponent equivalence, before being
    adopted in a detector pipeline.

These three decisions together complete the voter discrimination
without requiring a metadata gate. The voter case therefore
represents a fourth class of within-substrate discrimination,
complementing Class 1 (substrate registry), Class 2 (observable-
content prerequisite), and Class 3 (metadata-mechanism flag): pure
metric-level discrimination via thresholds calibrated against a
dense same-substrate discriminator ensemble. This pure-metric
class is preferable when achievable because it avoids the
maintenance burden of metadata flags that must be kept in sync
between models and detectors. It is achievable when the
target pattern's metric trajectory is sharply enough separated
from each near-neighbor pattern's trajectory; this is more often
true for early-time signatures than for late-time plateaus, which
is why P18 is built primarily on early-window metrics.
