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
on the tier itself.

The Sprint 21 5-seed audit at L = 64, density = 0.9, threshold = 0.375
(Schelling 1971's canonical 3/8) confirmed that no seed reaches
CONFIRMATION: four of five seeds fail screening because their
moran_final_qtr saturates at ≈ 0.24–0.29 (below the 0.30 floor),
and the remaining seed scrapes through screening with moran_final_qtr
≈ 0.301 but is rejected at confirmation because its wall_final_qtr
≈ 0.36 exceeds the 0.30 ceiling. Note that the rejection mechanism
is *not* a low wall plateau — Schelling's three-state grid {0, 1, 2}
counts every empty-cell boundary as a wall under the Moore-neighborhood
metric, yielding a wall plateau of ≈ 0.36 across all five seeds, far
above the DEFINITIVE wall_final > 0.05 floor. The discriminator is
instead the *combination* of moderate Moran plateau and elevated wall
plateau: voter coarsens both Moran *up* (≈ 0.55) and wall *down*
(≈ 0.21), while Schelling at threshold = 0.375 coarsens Moran weakly
*up* and wall weakly *down* but stops at a partition with too much
remaining boundary.

The metadata key `update` containing `copy`/`imitation`/`voter` is
therefore *defense in depth* rather than the primary discriminator at
canonical Schelling parameters — it would only become load-bearing
if a future model implemented a copying-driven mechanism that also
satisfied the geometry-bounded metric profile.

*A caveat at non-canonical Schelling parameters.* The Sprint 21 audit
also tested Schelling at threshold = 0.5 (the strong-segregation
parameter sometimes used in textbook Schelling expositions). At this
parameter, all five characterized seeds reach P18 DEFINITIVE with P1
exclusion marked "inconclusive" because Schelling's registry metadata
lacks any copy/imitation/voter `update` key. This is a metric-level
false positive: at threshold = 0.5, Schelling's wall plateau drops
to ≈ 0.27 (just below the 0.30 confirmation ceiling) and its Moran
plateau rises to ≈ 0.39 (within the [0.30, 0.75] definitive window).
The Class 4 pure-metric discrimination claim is therefore parameter-
contingent for the Schelling × P18 pair: it holds at canonical
threshold = 0.375 but breaks at threshold = 0.5. The recovery path
is either a stricter wall-density confirmation ceiling, a definitive-
tier downgrade keyed on the inconclusive P1 exclusion, or both —
recorded as Sprint 21 carry-forward #20b for a follow-up science
sprint.

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
without requiring a metadata gate at canonical Schelling parameters
(threshold = 0.375). The voter case therefore represents a fourth
class of within-substrate discrimination, complementing Class 1
(substrate registry), Class 2 (observable-content prerequisite),
and Class 3 (metadata-mechanism flag): pure metric-level
discrimination via thresholds calibrated against a dense
same-substrate discriminator ensemble. This pure-metric class is
preferable when achievable because it avoids the maintenance burden
of metadata flags that must be kept in sync between models and
detectors. It is achievable when the target pattern's metric
trajectory is sharply enough separated from each near-neighbor
pattern's trajectory; this is more often true for early-time
signatures than for late-time plateaus, which is why P18 is built
primarily on early-window metrics.

A consequence of the Sprint 21 audit, however, is that Class 4
purity should be understood as parameter-contingent rather than
unconditional. The claim "metric gates alone exclude Schelling
from P18 confirmation" is empirically true only inside the
discriminator-ensemble parameter range used at calibration time.
Outside that range — at Schelling threshold = 0.5, the
strong-segregation parameter — the same metric gates admit a
DEFINITIVE false positive. The honest characterization of Class 4
is therefore that pure-metric discrimination is *valid against the
parameter ensemble it was calibrated against*, and that any
extension of the catalog to a wider parameter neighborhood requires
either re-calibrating the gates or supplementing them with a
metadata-mechanism (Class 3) check. This is the cost of choosing
purity over defense-in-depth, and a future iteration of the P18
detector may revisit the trade-off.

## 6.11 Aggregate Grading Status

Across the 20-pattern inventory (Sprints 1–66), **19 / 20** patterns have
reached AT-DEPTH grade and **1 / 20** remains at GAP. The nineteen AT-DEPTH
patterns are P1 (aggregation), P2 (MIPS), P3 (Turing pattern formation),
P5 (flocking), P6 (milling), P7 (lane formation), P8 (self-organized jamming),
P9 (synchronization), P10 (chimera states),
P11 (predator-prey oscillations), P13 (excitable waves),
P14 (self-organized criticality), P15 (persistent computation),
P18 (coarsening-to-consensus), P21 (polarization), P22 (information cascade),
P27 (spatial reciprocity), P28 (wealth condensation), and
P31 (delayed gratification).

Five of the six AT-DEPTH patterns reached that grade via the Phase-2a
standard negative panel (§3.7): P9 returned PASS-with-weakness, P15
returned PASS, P18 returned PASS, and P27 returned PASS (Sprint 40),
each providing an independent confirmation of the AT-DEPTH grade assigned
via content-level discriminators; P28 and P31 reached AT-DEPTH prior to
the panel via
mechanistic null gates and the P31 non-redundancy test respectively.
Sprint 39 ran panels for P22 (SIR cascade) and P27 (Nowak-May spatial
reciprocity), returning PARTIAL and FAIL respectively; neither pattern
advanced to AT-DEPTH. Sprint 40 re-ran both panels after targeted fixes.
P27 advanced to PASS (TNR = 1.000, Cohen's d = 2.95) via an observable-
prerequisite guard (see §3.5, §4.8); P22 improved from overall TNR 0.519
→ 0.889 (Class C fixed: failed-regime TNR 0.000 → 1.000) but remained
PARTIAL due to catalog false positives (Lotka-Volterra and RPS). P27 now
has PASS/PASS/PASS/PASS across all four dimensions and advances to
AT-DEPTH. Sprint 41 added a literature-anchored irreversibility
prerequisite to P22 (see §3.5, §4.10), closing the Class B false positives
and advancing P22 dim4 to PASS (TNR = 1.000, Cohen's d = +∞); however,
dims 1–3 remain PARTIAL so P22 does not yet reach AT-DEPTH. The AT-DEPTH
count remains 6 / 19.

The transfer-matrix aggregate figures stand at: 20 models / 19 detectors
/ 79 compatible cells / 274 substrate-mismatch rejections / 27
observable-mismatch rejections / 19 display rows / 361 displayed cells /
77 displayed compatible / 284 displayed rejections. These figures are
unchanged from Sprint 26 and are locked by `scripts/count_transfer_matrix.py`
whose output is pinned in the pre-flight and post-flight checks at each
sprint boundary.

Sprint 42 ran the Phase-2a panel (v1.2) against P1 (Schelling segregation)
and attempted P3 (Gray-Scott). P1 returned PARTIAL (TNR=0.593, Cohen's
d=1.298): Moran's I fires on `linear_gradient` and `time_shuffled` synthetic
substrates (Class A) and on P11/P15/P12 lattice_2d catalog mates (Class B),
reflecting the broad-spectrum nature of spatial autocorrelation as a primary
metric. P3 panel was paused due to the lattice_2d_continuous substrate-undercount
condition (0 Class B mates; threshold is 3). AT-DEPTH count unchanged: **6 / 19**.

Sprint 43 addressed both patterns. For P3: two `lattice_2d_continuous`
supplements were added (`smooth_random_field`, `sinusoidal_traveling_wave`),
unblocking the Class B panel; the P3 detector was parameterised with
`stability_stride=5` to ensure ≥5 stability frames from 100-snapshot
histories; and P3 was reclassified as `time_shuffle_invariant=True` (each
Gray-Scott frame contains the complete stationary Turing pattern regardless of
temporal order). P3 panel v1.2 returned PASS (TNR = 1.000, Cohen's d = +∞);
P3 advances to AT-DEPTH. For P1: the type-constancy guard (Schelling 1971) was
extended to gate CONFIRMATION tier (CV threshold tightened to 0.01), closing
C-p1-class-b-lattice2d-fp; cat TNR advances to 1.000. Residual syn and fai
FPs remain unresolved; P1 panel returns PARTIAL (TNR = 0.704). AT-DEPTH count:
**7 / 19**.

Sprint 44 ran Phase-2a panels for P12 (spatial RPS, cyclic dominance) and P13
(Greenberg-Hastings, excitable waves) — both lattice_2d dim4 batch 3. P12 panel
v1.2 returned PASS (TNR = 1.000, Cohen's d = +∞): all five RPSSpatialModel
positives at mobility = 10⁻⁴ reached CONFIRMATION (confidence = 0.700); 10
high-mobility extinction regimes at mobility ∈ [5×10⁻³, 5×10⁻²] (11–111× M_c)
and 7 catalog lattice_2d mates were all correctly rejected. P12 dim4 advances
from PARTIAL to PASS; grade remains GAP (dim1: λ ∝ √M wavelength scaling not
replicated, carry-forward C3 open; dim2: single-seed characterisation). P13 panel
v1.2 returned PASS (TNR = 1.000, Cohen's d = +∞): all five GreenbergHastings
positives (n_states = 8, threshold = 1, Moore neighbourhood, density = 0.3)
reached SCREENING (confidence = 0.500); 10 low-density init regimes at density ∈
[0.01, 0.10] and 7 catalog mates were all correctly rejected. Positives at
SCREENING rather than CONFIRMATION reflects that 300-step trajectories on 50×50
grids accumulate wavefront statistics sufficient for the CV gate but not yet 50+
spiral-tip rotations; the class separation is sharp (0.500 vs. 0.000), yielding
Cohen's d = +∞. P13 dim4 advances PARTIAL→PASS; all four dimensions now PASS →
P13 advances to AT-DEPTH. AT-DEPTH count: **8 / 19**.

**Sprint 45 (P11 panel, AT-DEPTH +0).**

The P11 predator-prey oscillation detector (Sprint 11) was subjected to the Phase-2a
panel v1.2 with five seeds of LotkaVolterraLattice at canonical coexistence parameters
(predation_rate = 4.0, μ = σ = 1.0, L = 100, 1200 steps). Panel v1.2 returned PASS
(TNR = 1.000, Cohen's d = +∞): all five positives reached DEFINITIVE (confidence =
0.900); 10 predator-extinction regimes at predator_death_rate ∈ linspace(2.0, 5.0, 10)
(at/above the mean-field extinction boundary μ ≥ λ = 2.0) and 7 catalog lattice_2d
mates were all correctly rejected. P11 dim4 advances from PARTIAL to PASS. Grade
remains GAP (dim1 PARTIAL: no specific Fig/table from Mobilia-Georgiev-Täuber 2007
reproduced with stated tolerance). AT-DEPTH count unchanged: **8 / 19**.

**Sprint 46 (P5/P2/P6 panels, AT-DEPTH +1).**

Sprint 46 completes the continuous_2d dim4 batch by running the Phase-2a v1.2 panel
against P5 (Vicsek flocking), P2 (ABP / MIPS), and P6 (D'Orsogna milling). P5 returned
PASS-with-weakness (TNR = 0.957, Cohen's d = 4.987): all five positives reached DEFINITIVE;
the `time_shuffled` Class A substrate fires at DEFINITIVE (φ is per-frame invariant to temporal
order — carry-forward C-p5-time-shuffle-fp). P2 returned PASS (TNR = 0.958, Cohen's d = 3.401):
three of five positives reached DEFINITIVE; the `permutation_shuffled` substrate fires at
SCREENING (two_phase_score is spatial-distribution invariant to particle-index permutation;
carry-forward C-p2-perm-shuffled-fp per brief instruction not to auto-flip flag). P6 returned
PASS (TNR = 0.958, Cohen's d = 5.087): all five positives DEFINITIVE; `time_shuffled` fires at
DEFINITIVE (|L| per frame is temporal-order invariant — carry-forward C-p6-time-shuffle-fp).
All Class C failed-regime TNRs = 1.000. P5 dim4 advances PARTIAL→PASS with all other dims
already PASS → P5 advances to **AT-DEPTH**. P2 and P6 dim4 advance PARTIAL→PASS; other dims
remain PARTIAL. AT-DEPTH count: **9 / 19**.

**Sprint 47 (P8/P10 panels, AT-DEPTH +1).**

Sprint 47 closes the lattice_1d (P8 Nagel-Schreckenberg) and oscillator (P10 chimera)
dim4 depth-gaps by running the Phase-2a v1.2 panel. P8 returned PARTIAL (TNR = 0.652,
Cohen's d = 1.751): Class A yields two FPs (`permutation_shuffled` and `time_shuffled`
both fire at SCREENING — stopped_fraction is time-average-invariant to both operations;
P8 absent from the invariance registry; carry-forwards C-p8-perm-shuffled-fp and
C-p8-time-shuffle-fp). Class C yielded 6 FPs for regimes at rho ≥ 0.1167, all at or
above the NS jamming onset at p_slow=0.3; carry-forward C-p8-class-c-near-onset (restrict
next sweep to rho ∈ linspace(0.01, 0.09, 10)). P8 dim4 remains PARTIAL; escalate.

P10 returned PASS (TNR = 0.957, Cohen's d = 9.679): all five positives reached DEFINITIVE;
one Class A FP (`permutation_shuffled` at SCREENING — permuting oscillator indices preserves
the velocity distribution; carry-forward C-p10-perm-shuffled-fp). Class B catalog mate
(P9_kuramoto) and Class C (10 ordinary all-to-all Kuramoto above K_c) all correctly rejected.
P10 dim4 advances PARTIAL→PASS; all four dimensions now PASS → P10 advances to **AT-DEPTH**.
AT-DEPTH count: **10 / 19**.

Sprint 48 ran the Phase-2a panel v1.2 against P21 (Hegselmann-Krause opinion dynamics / polarization). The panel required new `opinions` detector_format plumbing added to the harness (synthetic.py, catalog.py, panel.py). All five canonical positives reached DEFINITIVE tier (confidence = 0.950). Class B and C TNR = 1.000. Class A TNR = 0.800 due to two expected FPs (`permutation_shuffled` at CONFIRMATION, `time_shuffled` at SCREENING; both carry-forwards C-p21-perm-shuffled-fp and C-p21-time-shuffled-fp document the invariance issue). Overall TNR = 0.913, Cohen's d = 4.543 → PARTIAL. P21 dim4 remains PARTIAL; AT-DEPTH count unchanged at **10 / 19**.

**Sprint 49 (invariance-flag batch, AT-DEPTH +0).**

Sprint 49 applies a batch of invariance-flag updates to `epc/phase2a/detector_invariance.py` for P1, P2, P5, P6, P8, and P21, derived from the Class A false-positive patterns documented in Sprints 39–48. The sprint re-runs all six affected panels under v1.2. Flags are grounded in the mathematical structure of each primary metric (see §3.5): permutation invariance for aggregate metrics over agents/cells; time-shuffle invariance for per-frame metrics whose mean is preserved under temporal reordering. The P10 carry-forward (C-p10-perm-shuffled-fp) is documented as an adapter artifact and intentionally not addressed by invariance flagging.

Re-run outcomes: P2 advances to TNR=1.000 (clean PASS; `permutation_shuffled` SKIPPED); P5 advances to TNR=1.000 (clean PASS; both substrates SKIPPED); P6 advances to TNR=1.000 (clean PASS; both substrates SKIPPED); P21 advances from PARTIAL to PASS-with-weakness (TNR 0.913→0.955, Cohen's d 4.543→5.487; `permutation_shuffled` SKIPPED; `time_shuffled` FP at CONFIRMATION (0.850) persists — convergence-timing issue, not invariance; C-p21-time-shuffled-fp remains open). P1 and P8 remain PARTIAL: their residual FPs are Class C calibration issues (C-p1-linear-gradient-fp, C-p1-class-c-subthreshold-fp, C-p8-class-c-near-onset) outside the invariance-flag mechanism's scope.

P21 dim4 advances from PARTIAL to PASS. No pattern reaches AT-DEPTH in Sprint 49 (P21 dims 1–3 remain PARTIAL). **AT-DEPTH count: 10 / 19** (unchanged). Seven carry-forwards closed (C-p1-time-shuffle-fp, C-p2-perm-shuffled-fp, C-p5-time-shuffle-fp, C-p6-time-shuffle-fp, C-p8-perm-shuffled-fp, C-p8-time-shuffle-fp, C-p21-perm-shuffled-fp).

**Sprint 50 (P11 dim1 closure, AT-DEPTH +1).**

Sprint 50 closes the P11 dim1 depth gap via numerical reproduction of
Mobilia-Georgiev-Täuber (2007) Sec. III (Fig. 3). The O(1/L)
oscillation-amplitude scaling law is confirmed: power-law fit across
L ∈ {30, 50, 100, 150} (3 seeds each) yields exponent −0.967
(R²=0.990; published −1.0; relative error 3.3%; tolerance ±0.20). The
coexistence + oscillatory focus regime is verified at L=100 (5 seeds;
mean predator density 0.081, FFT peak-to-mean 48.9). Large deviations
from the mean-field fixed point (ρ_prey measured 0.589 vs MF 0.250) are
a confirmed published finding of the paper (spatial correlations
dominant; Sec. III), not an implementation error. dim1 PARTIAL→PASS;
all four dimensions now PASS → P11 advances to **AT-DEPTH**.
**AT-DEPTH count: 11 / 19** (P3, P5, P9, P10, **P11**, P13, P15, P18,
P27, P28, P31).

**Sprint 51 (P22 dim1 closure, AT-DEPTH +0).**

Sprint 51 closes the P22 (information cascade / SIR epidemic CA) dim1 depth
gap via numerical reproduction of Datta & Acharyya (2021) §3.1.1 / Fig. 11
(arXiv:2104.10456). The paper's most quantitative claim is the wavefront
speed (linear-fit slope of R(t)): 0.4405 ± 0.0008 cells/step. The paper
uses a fixed-duration SIR CA (t_τ=4 deterministic steps per infection, Von
Neumann neighbourhood, p0=0.25, p1=0.97, p2=0.10 re-infection) rather than
stochastic geometric recovery. The paper's exact CA rules were implemented
inline in `analysis/reproductions/p22_dattaacharyya2005.py`; running 20 seeds
on a 200×200 lattice yields measured speed 0.4612 ± 0.0164 (relative error
4.7%, tolerance <15%). All 20 seeds achieve R² > 0.995 confirming linear
(superdiffusive) propagation. dim1 PARTIAL→PASS. P22 dims 2–3 remain PARTIAL
→ P22 grade remains GAP; AT-DEPTH count unchanged at **11 / 19**.

**Sprint 52 (P2 dim1 closure, AT-DEPTH +0).**

Sprint 52 closes the P2 (MIPS / active Brownian particles) dim1 depth gap via
numerical reproduction of Fily & Marchetti (2012) Fig. 2 (PRL 108, 235702). At
(φ=0.5, Pe=100), N=800 particles, rho_star=4.0, 2500 steps (500 burn-in), 5
seeds: (1) seed-mean two_phase_score = 0.1237 ± 0.077 (tolerance ≥ 0.10 from
Fig. 2 dense-cluster fractions f_gas ≈ 0.20–0.30, f_liquid ≈ 0.70–0.80);
(2) density-speed Pearson r = −0.958 ± 0.020 (|r| ≥ 0.70 confirming v(ρ) =
v₀(1 − ρ/ρ*) coupling); (3) thermal Pe=5 score = 0.052 < 0.08 (single-phase
homogeneous regime). All three tolerance checks PASS. Nucleation stochasticity
note: 2/5 seeds show low scores due to cluster-formation lag within the 2000-step
measurement window; the Pearson r is strong (|r| ≥ 0.944) in all five seeds,
confirming the mechanistic coupling is active throughout. dim1 PARTIAL→PASS.
P2 dims 2–3 remain PARTIAL → P2 grade remains GAP. dim1 PARTIAL count 3→2
(P12, P21 remain). AT-DEPTH count unchanged at **11 / 19**.

**Sprint 54 (P12 dim1 attempt, AT-DEPTH +0).**

Sprint 54 attempts to close the sole remaining dim1 gap (P12, spatial RPS)
by reproducing the spiral-wavelength scaling law λ ~ M^(1/2) from
Reichenbach, Mobilia & Frey (2007) Fig. 2c. L=100 lattice, σ=μ=1,
M ∈ {3×10⁻⁴, 4×10⁻⁴, 5×10⁻⁴}, radial ACF wavelength estimator,
10 seeds per M value. Measured log-log slope = 0.37 (target 0.5, tolerance
[0.4, 0.6]); outside tolerance. Wavelengths qualitatively match the formula
(within 15%, rank order correct) but the narrow M range (1.67×) with ~10%
per-point variance cannot confirm the exponent. P12 dim1 remains PARTIAL.
AT-DEPTH count unchanged at **11 / 19**.

**Sprint 55 (P14 dim2 closure, AT-DEPTH +1).**

Sprint 55 closes the P14 (BTW sandpile, self-organized criticality) dim2 depth gap
via a 20-seed multi-seed extension (L=32, n_drive=30,000, n_burn=3,000). The
power-law exponent τ = 1.2914 ± 0.0012 (CV = 0.09%) across all 20 seeds,
confirming that the SOC signature is robust to stochastic variation. The
narrow spread (τ ∈ [1.2895, 1.2932]) reflects the deterministic self-organization
of BTW dynamics: once the critical state is reached, the avalanche power-law is
seed-insensitive. dim2 PARTIAL→PASS; all four dimensions now PASS → P14 advances
to **AT-DEPTH**. **AT-DEPTH count: 12 / 19** (P3, P5, P9, P10, P11, P13, **P14**,
P15, P18, P27, P28, P31).

**Sprint 56 (P6 + P12 dim2 closures, AT-DEPTH +1).**

Sprint 56 closes dim2 on two patterns simultaneously. For P6 (D'Orsogna milling):
20 seeds with random initialisation (uniform random positions in [−R, R]², random
headings; canonical parameters N=100, α=1.0, β=0.5, Morse potential; warmup=2500
steps) yield mean |L| = 0.9818 ± 0.0301 (CV = 3.1%). All 20 seeds converge to
stable milling (min |L| = 0.884), confirming that the milling attractor is globally
reachable. dim2 PARTIAL→PASS; all four dimensions now PASS → P6 advances to
**AT-DEPTH**. For P12 (spatial RPS, cyclic dominance): 20 seeds at fixed M=10⁻⁴
yield mean spiral wavelength λ = 52.1 ± 10.4 (CV = 20.0%; n_valid=20/20). All
seeds produce a measurable spiral wavelength, confirming coexisting spiral state
robustness. dim2 PARTIAL→PASS; dim1 remains PARTIAL (λ ∝ √M scaling unresolved)
→ P12 grade stays GAP. **AT-DEPTH count: 13 / 19** (P3, P5, **P6**, P9, P10, P11,
P13, P14, P15, P18, P27, P28, P31).

**Sprint 53 (P21 dim1 closure, AT-DEPTH +0).**

Sprint 53 closes the P21 (polarization / fragmentation, Hegselmann-Krause)
dim1 depth gap via numerical reproduction of the canonical cluster-count
vs. ε curve from Hegselmann & Krause (2002) Fig. 2. N=100 agents, uniform
U[0,1] initial opinions, synchronous averaging, 20 seeds per ε. All 8 ε
points pass: fragmentation at ε ≤ 0.15 (median ≥ 3 clusters), polarization
at ε = 0.20 (median 2 clusters, 19/20 seeds), consensus at ε ≥ 0.30 (median
1 cluster, all seeds). At ε = 0.25 (the 2→1 transition boundary), 14/20 seeds
reach consensus and 6/20 reach two clusters; the published range is widened
to [1, 3] consistent with the paper's finite-N transition discussion. dim1
PARTIAL→PASS. P21 dims 2–3 remain PARTIAL → grade remains GAP. dim1 PARTIAL
count 2→**1** (P12 only). AT-DEPTH count unchanged at **11 / 19**.

**Sprint 57 (P2 + P21 + P22 dim3 closure, AT-DEPTH +0).**

Sprint 57 is a pure documentation sprint closing the dim3 (methods note)
depth gap for three patterns simultaneously: P2 (activity-induced phase
separation / MIPS), P21 (polarization / fragmentation, Hegselmann-Krause),
and P22 (information cascade / SIR). Full methods notes are authored in
`docs/methods_notes/` (p2_methods.md, p21_methods.md, p22_methods.md).
P2's note documents the Fily-Marchetti ABP overdamped Langevin equations,
the two_phase_coexistence_score primary metric, the Hartigan-dip failure
(ADR 44), and the FFT structure-factor relationship. P21's note documents
the synchronous HK update rule, L∞ convergence criterion, and sorted-gap
cluster counting (gap = ε/2). P22's note documents S/I/R state encoding
(0/1/2), independent-neighbours infection, the Sprint 41 irreversibility
prerequisite (S→I→R monotone, eliminates LV/RPS false positives), and
the Moore-neighbourhood percolation threshold context (p_c ≈ 0.038 at
q = 0.1). C5 carry-forward (methods notes thin for HK and SIR) is
**CLOSED**. All three patterns remain GAP (dim2 still PARTIAL for P2, P21,
P22). **AT-DEPTH count: 13 / 19** (unchanged).

**Sprint 58 (P12 dim1 wide-sweep attempt, AT-DEPTH +0).**

Sprint 58 extends the P12 Reichenbach-Mobilia-Frey (2007) Fig. 2c λ ∝ M^½
reproduction to a wide log-spaced M sweep [10^{−5}, 5×10^{−4}] (7 points,
15 seeds each, L=100, T_eq=1000 gen). The sweep obtains log-log slope = 0.107
(R²=0.769), outside the [0.40, 0.60] acceptance band. The failure reveals that
the analytical formula λ = 0.8·L·√(M/M_c) breaks down far below M_c: at
M ≤ 5×10^{−5} (M/M_c ≤ 0.11), measured wavelengths are flat (~42–44 lattice
units) versus the formula's prediction of 12–27, relative errors 59–269%. The
formula is derived via linearized theory near the extinction threshold and is
valid only in the near-M_c region (M/M_c ≳ 0.4). C2/C3 carry-forwards are
updated: the path forward is a dense near-M_c sweep (M ∈ [2×10^{−4}, 5×10^{−4}])
with ≥30 seeds per point and ≥5 M values, or use of a structure-factor estimator
more robust than ACF first-zero in the low-M regime. P12 dim1 remains PARTIAL.
**AT-DEPTH count: 13 / 19** (unchanged).

**Sprint 59 (P12 dim1 near-M_c dense sweep, AT-DEPTH +0).**

Sprint 59 narrows the P12 Reichenbach (2007) Fig. 2c reproduction to the formula-valid
near-M_c regime: M ∈ [2×10^{−4}, 5×10^{−4}] (7 points, 30 seeds each, L=100), with the
log-log fit restricted to M ∈ [2×10^{−4}, 4.5×10^{−4}]. The fit yields slope = 0.244
(R² = 0.918, n_fit_points = 6), outside the [0.40, 0.60] acceptance band. The R² is now
above 0.90 (confirming a clean linear relationship in log-log space), but the measured
slope is roughly half the published exponent (0.5). The systematic bias is attributed
to finite-size ACF compression: the radial ACF first-zero estimator underestimates λ
when λ/L > 0.5, and this bias grows nonlinearly with M near M_c. Reproducing the full
exponent requires L ≥ 200 (the original paper used L = 256) or a structure-factor
estimator less sensitive to finite-size boundary effects. P12 dim1 remains PARTIAL;
C2/C3 carry-forwards are updated with the finite-size diagnosis.
**AT-DEPTH count: 13 / 19** (unchanged).

**Sprint 60 (P2 + P21 + P22 dim2 closure, AT-DEPTH +3).**

Sprint 60 closes dimension 2 (multi-seed dispersion) for three patterns via 20-seed
campaigns. P2 (MIPS / ABP) at the canonical state (φ=0.5, Pe=100, N=800): two_phase_score
= 0.1134 ± 0.0790 (CV = 69.7%, reflecting nucleation stochasticity) with density-speed
Pearson r = −0.9585 ± 0.0196 (CV = 2.1%). P21 (Hegselmann-Krause) at ε=0.20 (N=100):
cluster count = 1.90 ± 0.31 (CV = 16.2%, median = 2). P22 (SIR wavefront speed) using
the paper-exact Datta-Acharyya (2021) CA (L=200): speed = 0.4606 ± 0.0163 cells/step
(CV = 3.5%, 19/20 valid seeds). All three patterns advance dim2 PARTIAL → PASS and,
with all four dimensions now PASS, advance to AT-DEPTH.
**AT-DEPTH count: 16 / 19** (+3: P2, P21, P22).

**Sprint 61 (P1 dim4 closure, AT-DEPTH +1).**

Sprint 61 resolves both P1 dim4 carry-forwards. (1) The `linear_gradient`
Class A false positive was a monotonic spatial partition (one connected
component per type) with high Moran's I from structure alone, not from
aggregation. A multi-cluster prerequisite grounded in Schelling (1971) —
genuine segregation produces multiple disconnected same-type clusters from
random initial conditions — cleanly excludes gradients while the canonical
positive (10–20 components per type) passes trivially. (2) The Class C regime
parameters were brief-author errors: threshold ∈ [0.161, 0.250] at
density=0.9 were above the empirical critical threshold (≈0.13), making them
true positives mislabeled as negatives. Corrected: threshold ∈ linspace(0.01,
0.10, 10), grid_size 32→50. Panel re-run: overall TNR = 1.000, syn = 1.000,
cat = 1.000, fai = 1.000, Cohen's d = +∞. P1 dim4 PARTIAL→PASS; all four
dimensions now PASS → P1 advances to AT-DEPTH.
**AT-DEPTH count: 17 / 19** (+1: P1). Remaining gaps: P8 (dim4), P12 (dim1).

**Sprint 62 (P8 dim4 closure, AT-DEPTH +1).**

Sprint 62 resolves the P8 dim4 carry-forward C-p8-class-c-near-onset. The
Class C failed-regime densities linspace(0.05, 0.20, 10) placed 6/10 regimes
at ρ ≥ 0.1167, above the empirical jamming onset ρ_c ≈ 0.10 for v_max=5,
p_slow=0.3 at L=1000 (stopped_fraction jumps from ~0 at ρ=0.09 to ~0.045 at
ρ=0.11, multi-seed validated). These regimes genuinely jam — mislabeled
negatives (brief-author error, same class as Sprint 40 P22 and Sprint 61 P1
corrections). Corrected: density ∈ linspace(0.02, 0.07, 10), all well below
onset in the free-flow phase (stopped_fraction = 0.0 across all densities and
seeds). Panel re-run: overall TNR = 1.000, syn = 1.000, cat = 1.000,
fai = 1.000, Cohen's d = +∞. P8 dim4 PARTIAL→PASS; all four dimensions now
PASS → P8 advances to AT-DEPTH.
**AT-DEPTH count: 18 / 19** (+1: P8). Remaining gap: P12 (dim1).

**Sprint 63 (P12 dim1 final attempt — accepted limitation, AT-DEPTH +0).**

Sprint 63 is the fourth and final attempt at the P12 Reichenbach (2007) Fig. 2c
λ ∝ √M scaling-law reproduction. Using L=200 (4× cells vs L=100) with a
zero-padded FFT structure-factor ring-peak estimator (replacing the ACF first-zero
method), the sweep yields slope = 0.161 (R² = 0.792), outside the [0.40, 0.60]
acceptance band. The per-M data reveal a fundamental constraint: the linearized
√M formula is valid only for M/M_c ∈ [0.44, 1.00], but in this narrow range
spiral wavelengths approach L (even at L=200), and near-M_c extinction effects
flatten the curve. Only 2 of 5 M values fall in the usable regime — insufficient
for slope determination. The scaling-law reproduction is accepted as a
**documented finite-size measurement limitation**: P12 is validated via Phase-2a
panel PASS (Sprint 44), dim2 multi-seed PASS (Sprint 56), and qualitative spiral
presence across all tested M values. Carry-forwards C2/C3 reclassified as
closed-as-documented-limitation (accepted). P12 dim1 remains PARTIAL; grade
remains GAP. This is a deliberate stopping point — the scaling law requires
L ≥ 500 with many points in the narrow formula-valid window, beyond this
project's compute budget.
**AT-DEPTH count: 18 / 19** (unchanged). P12 dim1 is the sole remaining gap,
accepted as a documented limitation.

**Sprint 65 (P7 implementation, AT-DEPTH +0).**
New pattern P7 (lane formation in counterflow) implemented: Helbing & Molnár
(1995) social-force model + Nowak & Schadschneider (2012) lane order parameter
detector. dim1 PASS (φ_final = 0.92, gain = 0.42, throughput = 99.8%), dim2 PASS
(20-seed: φ = 0.897 ± 0.091, CV = 10.2%), dim3 PASS (methods note). dim4 pending.
Pattern inventory: 20 models × 20 detectors. P7 grade: GAP (dim4 pending).
**AT-DEPTH count: 18 / 20** (inventory grew 19→20; AT-DEPTH unchanged).

**Sprint 66 (P7 dim4 closure, AT-DEPTH +1).**
Phase-2a panel v1.2 for P7 (lane formation): PASS-with-weakness.
Overall TNR = 0.955, Cohen's d = 6.932. All 5 positives reach CONFIRMATION
(0.700). Catalog mates (P2, P5, P6) correctly rejected — none carry population
labels. Class C: 8 weak-repulsion regimes (A ∈ [0.1, 0.8]) + 2 single-population
regimes, all rejected (fai TNR = 1.000). Content prerequisite added: counterflow
requires ≥2 opposing populations with minority ≥10% (Helbing & Molnár 1995).
Sole weakness: `time_shuffled` synthetic FP at screening — each frame preserves
per-frame φ_lane independently of temporal ordering (carry-forward
C-p7-time-shuffled-fp; same class as C-p21-time-shuffled-fp).
dim4 pending→PASS; all four dimensions PASS → P7 advances to **AT-DEPTH**.
**AT-DEPTH count: 19 / 20** (+1: P7). Remaining gap: P12 (dim1).

**Sprint 67 (P17 implementation, AT-DEPTH +0).**
New pattern P17 (distributed sensing / collective gradient detection)
implemented: Berdahl et al. (2013) speed-modulation model + collective
chemotactic index detector. dim1 PASS (CI slope vs log(N) = 0.133, Spearman
ρ = 0.90), dim2 PASS (20-seed CI = 0.394 ± 0.130, CV = 33.0%), dim3 PASS
(methods note: speed-modulation mechanism, SNR ∝ √N derivation, distinctness
from P5). dim4 pending. Pattern inventory: 21 models × 21 detectors.
P17 grade: GAP (dim4 pending).
**AT-DEPTH count: 19 / 21** (inventory grew 20→21; AT-DEPTH unchanged).

**Sprint 68 (P17 dim4 closure, AT-DEPTH +1).**
Phase-2a panel v1.2 for P17 (collective sensing): PASS.
Overall TNR = 1.000, Cohen's d = 11.117. Three literature-grounded
prerequisites: (1) field_samples present in history (rejects non-field
substrates), (2) individual SNR ≤ 3.0 (group advantage required — Berdahl
2013), (3) social cohesion ratio ≤ 0.20 (group must be cohesive). Class A:
9/9 rejected (1 skipped — permutation_invariant=True). Class B: P2, P5,
P6, P7 all lack field_samples → prereq rejection. Class C: 5 social_off +
5 field_too_strong, all rejected.
dim4 pending→PASS; all four dimensions PASS → P17 advances to **AT-DEPTH**.
**AT-DEPTH count: 20 / 21** (+1: P17). Remaining gap: P12 (dim1).

**Sprint 69 (P19 implementation, AT-DEPTH +0).**
New pattern P19 (emergent leadership / minority guidance) implemented: Couzin
et al. (2005) informed-minority Vicsek model + label-shuffle directional pull
detector. dim1 PASS (accuracy vs ρ: Spearman ρ = 1.0; all 5 tolerance checks
PASS), dim2 PASS (20-seed accuracy = 1.000 ± 0.000), dim3 PASS (methods note:
Vicsek+bias dynamics, pull metric, P5/P17/P18 distinctness). dim4 pending.
Pattern inventory: 23 models × 22 detectors, 104 compatible pairs.
P19 grade: GAP (dim4 pending).
**AT-DEPTH count: 20 / 22** (inventory grew 21→22; AT-DEPTH unchanged).

**Sprint 70 (P19 dim4 closure, AT-DEPTH +1; Wave 1 complete).**
Phase-2a panel v1.2 for P19 (emergent leadership): PASS.
Overall TNR = 0.960, Cohen's d = 5.418. Content prerequisite (Sprint 70):
early-window informed→naive leadership gap — in the convergence phase
(10–40% of trajectory), informed agents' alignment with θ_pref must exceed
naive agents' alignment, verifying that the minority leads the process
(Couzin 2005). Class A: 10/10 rejected (time_shuffled now blocked by
early-leadership prereq — scrambled temporal order destroys the convergence
signature). Class B: 5/5 rejected (P2, P5, P6, P7, P17 all lack
informed_mask → prereq rejection). Class C: 9/10 rejected (5 rho_zero +
4/5 bias_zero; 1 bias_zero FP at confirmation by chance alignment —
carry-forward C-p19-bias-zero-chance-alignment).
dim4 pending→PASS; all four dimensions PASS → P19 advances to **AT-DEPTH**.
**AT-DEPTH count: 21 / 22** (+1: P19). Remaining gap: P12 (dim1).
Completes Milestone B Wave 1 (P7 + P17 + P19 all AT-DEPTH).

**Sprint 72 (P24 implementation, AT-DEPTH +0).**
P24 (homeostatic regulation) — first Wave-2 pattern. ProportionalHomeostat +
IntegralHomeostat models, P24 detector with surrogate uncontrolled null, T1a
observation-bundle adapter, T1b cross-model generalization. New substrate type:
scalar_timeseries. dim1: deviation ratio 0.0027 (PASS). dim2: 20-seed CV=0.8%
(PASS). dim3: methods note (PASS). dim4 pending (Sprint 73).
**AT-DEPTH count: 21 / 23** (inventory grew 22→23; AT-DEPTH unchanged).

**Sprint 73 (P24 dim4 closure, AT-DEPTH +1; Wave 2 first pattern).**
P24 Phase-2a panel v1.2: PASS (TNR=1.000, Cohen's d=+inf). New
`scalar_timeseries` format added to all 10 Class A generators. Invariance
flags: permutation_invariant=True (single scalar variable — permutation
is degenerate), time_shuffle_invariant=True (deviation integral is
order-invariant for constant dt; confirmed by Sprint 73 first-run FP).
Class A: 8/8 evaluated substrates rejected (uncontrolled drift → screening
fail). Class B: 0 catalog mates + 2 scalar_timeseries supplements
(passive_ou_decay, uncontrolled_random_walk_scalar) both rejected. Class C:
2/2 (gain_zero_drift, no_perturbation) both rejected. All 5 positives
DEFINITIVE (0.900).
dim4 pending→PASS; all four dimensions PASS → P24 advances to **AT-DEPTH**.
**AT-DEPTH count: 22 / 23** (+1: P24). Remaining gap: P12 (dim1).
Completes Wave 2 first pattern.

**Sprint 74 (P26 implementation, AT-DEPTH +0).**
P26 (stochastic resonance) — second Wave-2 pattern. BistableDoubleWell +
ThresholdUnit models, P26 detector with coherent-response metric and
time-shuffle null, T1a noise-sweep observation-bundle adapter, T1b
cross-model generalization. New substrate type: noise_sweep_timeseries.
dim1: gain = 0.855 over zero noise (PASS). dim2: 20-seed campaign: peak
coherent response = 0.897 ± 0.049 (CV=5.4%), gain = 0.833 ± 0.049 (CV=5.8%),
all 20 seeds DEFINITIVE (PASS). dim3: methods note (PASS). dim4 pending
(future sprint).
**AT-DEPTH count: 22 / 24** (inventory grew 23→24; AT-DEPTH unchanged).

**Sprint 75 (P26 dim4 closure, AT-DEPTH +1; Wave 2 second pattern).**
P26 (stochastic resonance) Phase-2a panel v1.2 PASS (TNR=1.000, d=+inf).
New `noise_sweep` detector format. Invariance flags: permutation_invariant=True
(scalar noise-sweep, single variable), time_shuffle_invariant=True (coherent
response order-invariant within noise-level groups; confirmed by Sprint 75
first-run FP). Content prerequisite: inverted-U shape at screening (Gammaitoni
1998: SR requires interior peak with rise-then-fall; monotone or flat
performance is not SR). Class A: 8/8 TN. Class B: 2/2 TN (advisory).
Class C: 2/2 TN (suprathreshold_signal, extreme_noise_only). All 5 positives
DEFINITIVE (0.900).
**AT-DEPTH count: 23 / 24** (+1: P26). Remaining gap: P12 (dim1).

**Sprint 76 (P23 anti-coordination / emergent load balancing — Milestone B Wave 2).**
Implemented P23 end-to-end: Minority Game (Challet & Zhang 1997) + El Farol Bar
(Arthur 1994) models, P23AnticoordinationDetector with random-choice surrogate null.
New substrate type: choice_timeseries. T1a observation-bundle adapter for
attendance/choice time series. T1b cross-model: El Farol detected at confirmation
(independent implementation validates phenomenon recognition). dim1: Savit curve
reproduced (σ²/N vs α, interior minimum at α≈0.32, σ²/N=0.077 << baseline 0.25).
dim2: 25-seed campaign (σ²/N = 0.075 ± 0.006, CV=8.7%, all detected). dim3:
methods note. dim4 pending (Sprint 77).
**AT-DEPTH count: 23 / 25** (inventory grew 24→25; P23 GAP pending dim4).

**Sprint 77 (P23 dim4 Phase-2a panel — Milestone B Wave 2 closure).**
Phase-2a v1.2 panel for P23 anti-coordination: **PASS** (TNR=1.000, Cohen's
d=14.504). Invariance flags: perm_inv=True, time_shuffle_inv=True (σ²/N is
the primary confirmation signal and is a distribution statistic preserved by
time shuffling). Content prerequisite: nondegenerate variance + variance below
random baseline at confirmation (Savit et al. 1999). Class A: 8/8 evaluated
synthetics rejected. Class B: 2/2 supplements rejected (consensus_herding,
random_choice). Class C: 2/2 failed regimes rejected (random_agents, herding).
P23 dim4 pending→PASS; all 4 dims PASS → **AT-DEPTH**. Milestone B Wave 2
(P24, P26, P23) complete.
**AT-DEPTH count: 24 / 25** (P12 sole remaining GAP).

**Sprint 80 — P16 associative memory dim4 Phase-2a panel.** Panel v1.2 PASS
(TNR=1.000, d=+inf). New `state_vector` detector format for attractor-network
histories. Invariance: completion accuracy invariant under consistent neuron-index
permutation (Hopfield 1982); converged fixed-point state replicated across post-
convergence steps → both shuffled substrates skipped. Content prerequisite: ≥2
distinct selectively-retrievable stored patterns; single-pattern (P=1) convergence
and single-attractor collapse rejected. Over-capacity (α=1.0) rejected at screening
(best-overlap 0.464 < 0.500). 5 positives DEFINITIVE. P16 dim4 pending→PASS;
all 4 dims PASS → **AT-DEPTH**. Milestone B Wave 3 first dim4 closure.
**AT-DEPTH count: 25 / 26** (P12 sole remaining GAP).

**Sprint 81 — P25 canalized restoration / equifinality implementation.** New
P25 detector + CanalizedLandscape + MultiBasinGRN models. New substrate:
`canalization_landscape` (12th). Convergence variance ratio as primary metric
with IC-distribution surrogate null. Trivial-collapse gate prevents false
positives from instant constant maps. Two negative controls: DiffusiveDynamics
(ICs diverge → screening fails) and TrivialCollapse (instant map → trivial-
collapse gate rejects). dim1: ratio ≈ 0 PASS, basin 1.0 PASS, DEFINITIVE.
dim2: 20 seeds, all DEFINITIVE. dim3: methods note. T1b: MultiBasinGRN
DEFINITIVE. Registry: 31 models × 27 detectors = 837 cells, 112 compatible
pairs. dim4 pending (Sprint 82).
**AT-DEPTH count: 25 / 27** (P12, P25 GAP — P25 dim4 pending).

**Sprint 82 — P25 equifinality dim4 Phase-2a panel.** Panel v1.2 PASS
(TNR=1.000, d=+inf). New `canalization_bundle` detector format for multi-IC
observation bundles. Invariance: convergence variance ratio permutation-invariant
over IC ensemble (aggregate statistic) + time-shuffle degenerate (step labels
preserved during extraction). Content prerequisite: basin_volume ≥ 0.5 at
screening (Waddington 1957: equifinality requires wide-IC convergence),
discriminating from noisy homeostatic regulation (P24). Class A: 8/8 rejected.
Class B: 2/2 supplements rejected (diffusive_multi_ic, homeostatic_regulation_bundle).
Class C: 2/2 failed regimes rejected (narrow_basin, divergent_dynamics). 5
positives DEFINITIVE. P25 dim4 pending→PASS; all 4 dims PASS → **AT-DEPTH**.
**AT-DEPTH count: 26 / 27** (P12 sole remaining GAP).

**Sprint 83 — P20 quorum sensing / threshold-activated response.** New pattern
P20 implemented end-to-end: AutoinducerQuorum model (mean-field bacterial QS
ODE with positive-feedback hysteresis), P20QuorumSensingDetector (step-function
R² primary metric, density-shuffle null), FractionThresholdModel (T1b cross-
model), GradedResponseModel (negative control). New substrate type
`density_sweep_timeseries` (13th). T1a observation bundle: density-sweep
with up/down sweeps for hysteresis detection. dim1: step R²=1.000, hysteresis
width=1.190, DEFINITIVE. dim2: 20 seeds all DEFINITIVE. dim3: methods note.
T1b: FractionThresholdModel DEFINITIVE. Negative control (graded response)
rejected at screening. Registry: 33 models × 28 detectors = 924 cells, 114
compatible pairs. dim4 pending (Sprint 84).
**AT-DEPTH count: 26 / 28** (P12, P20 GAP — P20 dim4 pending).

**Sprint 84 — P20 dim4 Phase-2a panel (Milestone B Wave 3 completion).** Phase-2a
panel v1.2 **PASS** (TNR=1.000, Cohen's d=+inf). All 5 positives DEFINITIVE
(0.900); all negatives correctly rejected. Class A: 9/9 TN (time_shuffled
SKIPPED — tag-based grouping makes primary metric time-shuffle-invariant).
Class B: 0 catalog mates + 2 supplements (smooth_sigmoid_density_sweep,
always_off_density_sweep — both rejected at screening). Class C: 2/2 TN
(sub_threshold, graded_response). Key P18 discrimination: graded consensus
response lacks sharp threshold and hysteresis. Invariance: perm_inv=False,
time_shuffle_inv=True. P20 advances GAP→AT-DEPTH. Completes Wave 3.
**AT-DEPTH count: 27 / 28** (P12 sole GAP).

**Sprint 85 — Milestone B Wave 3 summary (non-blocking checkpoint).** No new
code; summary sprint documenting Wave 3 completion status. Wave 3 delivered
P16, P25, P20 (all AT-DEPTH). AT-DEPTH count unchanged.
**AT-DEPTH count: 27 / 28**.

**Sprint 86 — P4 territoriality / exclusion boundaries implementation.** New
pattern P4 (Giuggioli et al. 2011 scent-mediated territorial exclusion).
ScentMarkingModel (canonical, N=4 agents on L=48 torus, 20000 steps,
two-stage movement: foreign-scent avoidance + own-scent-weighted softmax),
PheromoneRepulsionModel (T1b, hard-threshold avoidance), PlainRandomWalkModel
(negative control). P4 detector: exclusivity index + pairwise overlap +
boundary persistence + occupancy-scent correlation; cell-level multinomial
shuffle null. T1a observation bundle: territorial_agent_field format.
dim1: Giuggioli 2011 reproduction (excl=0.902, overlap=0.034,
persistence=0.865, d=157.5). dim2: 20-seed campaign. dim3: methods note.
dim4 pending (Sprint 87). Registry: 35 models × 29 detectors, 116 compatible
pairs.
**AT-DEPTH count: 28 / 29** (+1: P4, provisionally AT-DEPTH pending dim4).

**Sprint 87 — P4 dim4 Phase-2a panel.** Phase-2a panel v1.2 **PASS**
(TNR=1.000, Cohen's d=4.153). New `territorial_agent_field` detector format
wired into panel harness. Invariance flags: permutation_invariant=True
(exclusivity invariant under agent-index relabelling),
time_shuffle_invariant=True (cumulative occupancy preserves ownership
structure under temporal reordering — corrected from brief's prescribed
False). Content prerequisite (Sprint 87): scent-mediated exclusion required
at screening — occ-scent correlation < 0 AND persistence ≥ 0.5
(Giuggioli 2011: genuine territorial exclusion requires negative
occupancy-scent coupling). First-run FPs fixed: (1) mixing-adequate null
histories (3000 internal steps for O(L²) random-walk coverage on L=32 grid),
(2) content prerequisite gates incidental spatial autocorrelation,
(3) Class C home_attraction=0 for high-tolerance regimes (prevents
own-scent self-organization without foreign avoidance). Class A: 8/8 TN
(permutation_shuffled + time_shuffled SKIPPED). Class B: 0 catalog mates +
2 supplements (random_walk_territory, clustering_agents_territory — both
rejected at screening). Class C: 10/10 TN (5 high-tolerance + 5 fast-decay).
4/5 canonical positives detected (seed 2 fails occ-scent prereq;
mean_score=0.60). P4 dim4 pending→PASS; all 4 dims PASS → **P4 AT-DEPTH**
confirmed.
**AT-DEPTH count: 28 / 29** (unchanged; P4 provisionally AT-DEPTH since
Sprint 86). Remaining gap: P12 (dim1).

**Sprint 89 (P29 dim4 closure, AT-DEPTH +1; Wave 4 first dim4 closure).**
P29 trail / network formation (Tero 2010) Phase-2a panel v1.2 PASS
(TNR=1.000, Cohen's d=10.550). trail_network detector format. Invariance:
permutation_invariant=False, time_shuffle_invariant=False (10/10 Class A
evaluated, no skips). Two content prerequisites added (Sprint 89):
(1) temporal reinforcement dynamics — edge weights must change over time,
blocking static/precomputed graphs (Tero 2010: conductance evolves from
uniform to near-MST); (2) weight accumulation — late total edge weight must
exceed early by ≥50%, blocking i.i.d. random-weight-per-step nulls and
time-shuffled positives (Deneubourg 1990: pheromone reinforcement produces
net accumulation on preferred edges). Class A: 10/10 TN. Class B: 0 catalog
mates + 2 supplements (`static_mst_graph` rejected at screening via
temporal-dynamics prereq; `uniform_traffic_graph` rejected at screening via
weight-accumulation prereq). Class C: 5 high-evaporation
(evaporation_rate ∈ [0.80, 0.99]) + 5 no-reinforcement (NoReinforcementModel),
all 10 rejected. 4/5 canonical positives detected (1 DEFINITIVE, 2 SCREENING
at 0.600, 1 SCREENING at 0.500; seed 0 not detected due to stochastic layout).
P29 dim4 pending→PASS; all 4 dims PASS → **P29 AT-DEPTH** confirmed.
**AT-DEPTH count: 29 / 30** (+1: P29). Remaining gap: P12 (dim1).
