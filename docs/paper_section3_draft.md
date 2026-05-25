# Section 3: Detection Toolkit

The core methodological contribution of this work is not the pattern catalog
itself — taxonomies are inexpensive to propose — but the accompanying
detection toolkit: a library of quantitative detectors, each with explicit
null models, tier-gated significance criteria, and nearest-neighbor exclusion
logic. This section describes the design principles, detection architecture,
and key methodological choices that make the toolkit reliable.

## 3.1 Design Principles

Five principles govern the toolkit's design.

*Quantitative, not qualitative.* Every detector returns a structured result
with numerical metrics, p-values, effect sizes, and confidence scores —
never a bare boolean. This enables cross-model comparison and meta-analysis
across the catalog.

*Tier-gated significance.* Detection is organized into three tiers —
screening, confirmation, and definitive — with progressively stricter
requirements. Screening catches candidates at the cost of some false
positives. Confirmation requires secondary metrics and a null-model test
(p < 0.01). Definitive requires all confirmation criteria plus
nearest-neighbor exclusion tests and (where available) a mechanistic null
model. Confidence scores are capped by tier: a screening result cannot
exceed 0.60 regardless of signal strength, preventing overconfidence in
preliminary assessments.

*Intrinsic timescales.* Fixed timestep thresholds are fragile across
implementations and parameter choices. All persistence and duration
requirements use system-intrinsic timescales: domain crossing time T_cross
for particle systems, grid propagation time T_prop = L/v for cellular
automata, oscillation period T_osc for coupled oscillators, and convergence
time T_sort for sorting processes. This ensures that a persistence
requirement of "≥ 10 T_cross" adapts automatically to the system's natural
dynamics.

*Nearest-neighbor exclusion.* Many patterns share observational signatures.
Flocking (P5) and milling (P6) both involve coordinated motion; excitable
waves (P13) and persistent computation (P15) both produce spatial dynamics
on lattices. Every detector specifies which neighboring patterns must be
excluded before a detection is considered definitive. The exclusion graph
is asymmetric — P5 must exclude P6/P7, but P6 need only exclude P5 — and
is derived from the pattern catalog's ontological dimensions.

*Observable scope transparency.* Detectors declare whether they operate on
state history alone, benefit from model metadata, or require it. A detector
marked "state_history_only" can be applied to any system that produces
the right observables, even without knowing the model's rules. A detector
marked "model_metadata_required" needs access to the rule structure
(e.g., verifying that a game has Prisoner's Dilemma payoff ordering).
This distinction matters because the long-term goal is to apply these
detectors to systems whose rules are unknown.

## 3.2 Detector Architecture

Each detector follows a common pipeline implemented in a base class:

1. **Timescale estimation.** The detector estimates or receives the
   system-intrinsic timescale from model metadata or state history
   analysis.

2. **Prerequisite validation.** Guards check that the state history meets
   minimum requirements: sufficient run length (typically ≥ 5τ for screening,
   ≥ 10τ for confirmation), adequate system size (N ≥ 50 for meaningful
   statistics), and required observable keys. Failures produce warnings
   that propagate to the final result. Certain prerequisites — such as the
   requirement that excitable media have at least three states — are hard
   gates that prevent false detection.

3. **Primary metric computation.** The defining observable for the pattern:
   polarization φ for flocking, Moran's I for aggregation, angular momentum
   |L| for milling, wavefront speed CV for excitable waves. The primary
   metric determines whether the screening threshold is met.

4. **Screening check.** A relaxed threshold on the primary metric. Designed
   to have high sensitivity (few false negatives) at the cost of moderate
   specificity (some false positives).

5. **Secondary metrics.** Additional measurements that distinguish the
   target pattern from confounders: group speed ratio R distinguishes
   flocking from milling, segregation index distinguishes aggregation from
   random clustering, ring density profile distinguishes milling from
   random rotation.

6. **Null model test.** A permutation, surrogate, or mechanistic null model
   generates a reference distribution for the primary metric under the
   hypothesis that the pattern is absent. Three null types are distinguished:
   shuffle nulls permute observed data (e.g., random heading assignment
   preserving positions); surrogate nulls generate synthetic data matching
   marginal statistics; mechanistic nulls modify the model's rules (e.g.,
   zero coupling). Mechanistic nulls provide the strongest evidence but
   require access to the model, so they are optional for confirmation and
   required only for definitive tier.

7. **Confirmation and definitive checks.** Stricter thresholds on the primary
   metric, combined with secondary metric requirements, null p-value
   thresholds, and nearest-neighbor exclusion results.

The pipeline produces a DetectorResult dataclass containing the pattern ID,
detection tier, confidence score, all computed metrics, null model results,
exclusion outcomes, and warnings. This structured output enables downstream
analysis of detection patterns across the model inventory.

## 3.3 Null Model Taxonomy

The choice of null model determines the strength of evidence a detection
provides. We distinguish three categories with increasing inferential power.

**Shuffle nulls** permute observed data while preserving spatial or temporal
structure. For P5 (flocking), headings are replaced with independent uniform
random draws while positions and speeds are preserved. For P1 (aggregation),
agent type labels are randomly permuted across spatial positions. Shuffle
nulls test whether the observed metric value could arise from random
arrangement of the observed elements, without any coupling.

**Surrogate nulls** generate synthetic trajectories matching marginal
statistics. For TE estimation, phase-randomized surrogates preserve the
power spectrum while destroying temporal coupling. Surrogate nulls test
whether the observed temporal structure is consistent with uncoupled dynamics
having the same spectral properties.

**Mechanistic nulls** modify the model's rules to remove the hypothesized
mechanism. For P9 (synchronization), coupling is set to zero while preserving
natural frequencies. For P13 (excitable waves), the refractory period is
eliminated, producing immediate re-excitation rather than wave propagation.
Mechanistic nulls provide the strongest evidence because they directly test
the causal mechanism, but they require access to the model's rule structure.

The tier system maps onto this taxonomy: screening requires only a threshold
check, confirmation requires a shuffle or surrogate null at p < 0.01, and
definitive requires all confirmation criteria plus (where the observable
scope allows) a mechanistic null.

## 3.4 Substrate-Aware Dispatch

Not every detector should be applied to every model. The P5 flocking
detector expects continuous-space positions and headings; applying it to a
lattice cellular automaton is meaningless. The orchestration layer prevents
these cross-substrate mismatches through a two-level compatibility check.

The first level is **substrate type matching.** Seven substrate types are
defined: lattice_1d (Zhang sorting; Nagel-Schreckenberg traffic),
lattice_2d (Schelling, Nowak-May, SIR, RPS, lattice Lotka-Volterra,
voter, plus the cellular-automata models Greenberg-Hastings, Game of
Life, BTW sandpile), lattice_2d_continuous (Gray-Scott reaction-
diffusion, distinguished from lattice_2d by the `n_unique_values ≥ 50`
content prerequisite that separates continuous-field models from
integer-state cellular automata), continuous_2d (Vicsek, D'Orsogna,
active Brownian particles), oscillator (all-to-all Kuramoto, non-local
Kuramoto), opinion_space (Hegselmann-Krause), and scalar_wealth
(Yard-Sale wealth exchange). LV's prey/predator/empty type information
distinguishes it from cyclic-competition systems via the
`n_unique_species_observed == 2` content prerequisite of P11 rather
than via a dedicated substrate. Each detector declares which substrate
types it is compatible with.

The second level is **observable matching.** Within a compatible
substrate, the detector requires specific state keys. P27 (spatial
reciprocity) requires the `coop_fraction` observable that only
Prisoner's Dilemma models produce, preventing it from firing on other
lattice_2d systems like Greenberg-Hastings or Game of Life. P14
(self-organized criticality) requires `avalanche_sizes`, restricting
it to sandpile-type models. P11 (predator-prey oscillation) requires
prey and predator time series plus a nontrivial-empty-reservoir
prerequisite (`std(species_A + species_B) > 0.005`), which both
distinguishes it from conservation-locked two-state systems like
Nowak-May and — combined with the `n_unique_species_observed == 2`
gate — separates it from three-species cyclic systems like RPS.

Of 380 possible model × detector pairs in the current inventory
(20 registered models × 19 detectors), 79 are substrate-compatible
and observable-compatible. The remaining 301 cells are correctly
eliminated by substrate mismatch (274) or by detector-observable
incompatibility (27, primarily P31 which requires `lattice_1d` with
`cell_types`, P14 which requires `avalanche_sizes`, and P27 which
requires `coop_fraction`). The cross-detection regression suite
(`tests/test_cross_detection_matrix.py::EXPECTED_OUTCOMES`) pins 195
audited model–detector outcomes — covering compatible-and-tested
cells, expected substrate-mismatch rejections, and observable-
prerequisite rejections — and the 19 canonical DEFINITIVE positives
are pinned separately in dedicated end-to-end test files. This
block-diagonal structure ensures that detectors operate only in
domains where their metrics are meaningful.

## 3.5 Key Boundary Tests

Three boundary tests illustrate the toolkit's discrimination capabilities
and the care required to avoid false positives.

**P13/P15 Transfer Entropy discriminator.** Excitable waves (P13) and
persistent computation (P15) both produce spatial dynamics on lattice_2d
substrates, but differ in information-processing character. We developed a
boundary-conditioned transfer entropy test that cleanly separates them:
Greenberg-Hastings spiral waves produce TE ratios of 1–2× above a
permutation null, while Game of Life produces 15–16×. Critically, raw
average TE gives the *wrong* ordering (GH > GoL) because deterministic wave
propagation creates trivially high TE at interior cells. Restricting
measurement to boundary cells — where structures interact — isolates the
computational signal. This methodological contribution emerged from the
replication process itself.

**P1 type constancy guard.** Game of Life generates significant spatial
autocorrelation from B3/S23 dynamics — alive cells cluster — initially
triggering P1 (aggregation). This was resolved by recognizing that P1
detects aggregation of *persistent agent types*, not merely spatial
autocorrelation from any source. GoL alive/dead states change every
timestep and are therefore not persistent type labels. A type constancy
check verifies that the distribution of values at each spatial position is
stable over the measurement window, combined with a temporal convergence
guard that verifies monotonic increase toward a plateau (the signature of
genuine aggregation dynamics).

**P13 excitable medium guard.** Models with fewer than three states (binary
alive/dead for GoL, cooperator/defector for Nowak-May) can produce spurious
wavefront speed statistics — GoL's synchronous update gives every cell
identical inter-excitation intervals, producing CV = 0 and passing the
screening threshold. A hard prerequisite guard requires n_states ≥ 3
(corresponding to the resting, excited, and refractory states necessary
for excitable dynamics), preventing detection on structurally incompatible
models.

**P11 conservation and bilateral-vs-cyclic guards.** The P11
predator-prey oscillation detector combines two prerequisites that
emerged from broad negative-model sweeps. The anti-phase
cross-correlation primary metric (ρ_anti) false-positives on strictly
conserved two-state systems — Nowak-May's cooperator and defector
fractions sum to 1 exactly and are therefore anti-correlated at
ρ = −1 at short lag by algebra alone, independent of any dynamical
coupling. A total-std prerequisite (`std(species_A + species_B) >
0.005`, requiring a nontrivial empty reservoir) catches this: Nowak-May
has total_std = 0.000; Lotka-Volterra has total_std ≈ 0.034. The
primary metric also triggers strongly on three-species cyclic systems:
spatial RPS produces ρ_anti ≈ −0.94 on any pair of species, *stronger*
than Lotka-Volterra's ρ_anti ≈ −0.86 on the predator-prey pair. The
n_unique_species_observed prerequisite (`== 2`) keeps P11 specific to
bilateral coupling rather than firing on cyclic dominance. Both
prerequisites were identified by testing the proposed primary metric
against every existing two-species lattice model before locking the
detector, not against the planned negatives alone.

**P27 coop_fraction observable prerequisite (Sprint 40).** P27's primary
metric — the late-stage cooperation fraction (`f_C`) and its spatial
autocorrelation — is definitionally absent on substrates that do not
track a cooperator/defector distinction. However, the v1.2 panel
(Sprint 39) found that P27 fired at SCREENING on 8/9 generic lattice_2d
substrates: because the panel runner pre-computed `coop_fraction =
(grid == 0).mean()` to avoid a KeyError, any substrate with ≥2%
zero-valued cells (GoL, GH, RPS, voter, etc.) satisfied the screening
criterion `f_C > 0.02 AND n_gen > 100`. Sprint 40 added a hard
prerequisite guard at the top of `detect_p27`: if the state history
lacks the `coop_fraction` key (which `NowakMayModel` always provides
natively), the detector short-circuits immediately with
`detected=False, confidence=0.0`. Concurrently the panel runner's
synthetic augmentation of `coop_fraction` was removed. This guard is
the P27 analog of P11's `total_std` conservation prerequisite: a
content-level domain restriction that prevents out-of-domain misfires
without altering the detector's behavior on its native substrate
(Nowak-May with the cooperator/defector observable). After the fix, the
Phase-2a panel v1.2 returned PASS for P27 (TNR = 1.000, Cohen's d =
2.95, Sprint 40).

**P22 irreversibility prerequisite (Sprint 41).** Following Datta &
Acharyya (2021), P22 requires that observed state transitions exhibit
the irreversible S→I→R flow: once a cell enters the Recovered state
(integer value 2), it never returns to Susceptible (0) or Infected (1).
Substrates without this property — e.g., the Mobilia-Georgiev-Täuber
(2007) LV lattice, where predator death produces the backward transition
PREDATOR(2)→EMPTY(0), or the Reichenbach (2007) spatial RPS, where cyclic
dominance sends species through EMPTY(0) repeatedly — trigger a hard
short-circuit before primary-metric evaluation, returning
`detected=False, confidence=0.0` at SCREENING tier. The Phase-2a panel
v1.2 re-run (Sprint 41) confirmed that this guard eliminates the two
Class B false positives (P11 LV and P12 RPS) that persisted through
Sprint 40, bringing P22 to TNR = 1.000 (Cohen's d = +∞). This prereq
is the P22 analog of P11's `total_std` conservation guard and P27's
`coop_fraction` observable guard: a content-level domain restriction
grounded in published characterization of the canonical model, not an
ad-hoc threshold adjustment.

**P1 type-constancy guard extension to CONFIRMATION (Sprint 43).** Following
Schelling (1971) "Dynamic Models of Segregation", P1 specifically detects
aggregation of *intrinsic type labels* — agent identities that never change
across the simulation; only positions change via tolerance-threshold moves.
The Sprint 42 Phase-2a panel found three Class B false positives at
CONFIRMATION tier: P11 Lotka-Volterra, P15 Game of Life, and P12 spatial RPS.
All three are lattice_2d models whose cell identities transition dynamically
(predator/prey cycling, alive/dead flipping, RPS dominance) — they are
fundamentally outside P1's domain per Schelling (1971). Sprint 43 extended the
existing type-constancy guard (previously applied only at DEFINITIVE tier) to
also gate CONFIRMATION: if the coefficient of variation of non-background type
counts exceeds 0.01 across the trajectory, the detector short-circuits at
SCREENING with `detected=False`. Schelling segregation has CV = 0.000
(agents are perfectly conserved); LV, GoL, and binarized-RPS (catalog adapter
converts 4-state RPS grid to occupied/empty) have CV ≥ 0.014 > 0.01. The
Phase-2a panel re-run (Sprint 43) confirmed cat TNR advances from 0.571 to
1.000, with C-p1-class-b-lattice2d-fp CLOSED.

## 3.6 Statistical Power Requirements

Permutation-based significance testing imposes minimum sample sizes that we
enforce across the toolkit. With n permutations, the smallest achievable
p-value is 1/(n + 1). For a significance threshold of p < 0.01, at least
99 permutations are required; for p < 0.005, at least 199. Our standard
configuration uses 999 permutations for P1 (yielding p floor of 0.001),
199 for P5/P6/P9/P22/P27 (p floor 0.005), 199 for P12 CONFIRMATION and
≥499 for DEFINITIVE, 99 for P11 (where the null p-value is not a tier
gate — see below), and 99 for TE-based tests (p floor 0.01). The P31
non-redundancy test requires at least 500 model runs for stable 10-fold
cross-validation with 8+ features. For P11 specifically, the canonical
positive requires ≥1,200 generations at L = 100 for DEFINITIVE tier
(effect-size Cohen's d ≤ −1.5 is the discriminator); shorter runs drop
to CONFIRMATION.

P11's null-model behavior is an instructive special case. Its
circular-shift null preserves each series' autocorrelation and FFT
magnitude spectrum, so on Lotka-Volterra the null frequently produces
extreme anti-correlations (null 5th percentile ≈ observed value; p in
the range 0.05–0.15 even at Cohen's d ≈ −2). This is not a bug — the
null correctly reflects that oscillating series are autocorrelated, and
autocorrelation alone is insufficient evidence for predator-prey
coupling. Tier gating therefore uses Cohen's d against the null
distribution rather than the one-sided p-value. The p-value is
reported as a diagnostic only. This pattern — a null model that is
intentionally too strong for clean p-values, with effect-size gating as
the compensating discriminator — may recur for other detectors whose
primary metric is a phase-relationship rather than a magnitude. P18
(Sprint 20) encountered the same circular-shift autocorrelation
preservation issue on its Moran's I trajectory; rather than switching
to effect-size gating, P18 resolved it by replacing the circular-shift
null with a full random permutation that destroys both the trend and
the autocorrelation, achieving p < 0.01 reliably across all canonical
voter seeds. The two resolutions — Cohen's-d gating with the
autocorrelation-preserving null, or full random permutation that
destroys autocorrelation — are both legitimate; the choice depends on
whether the autocorrelation structure is itself part of the pattern's
substrate (P11) or merely a property of the test statistic (P18).
See §6.10 (Decision 54) for the P18 reasoning.

These are not arbitrary choices — in multiple cases, underpowered initial
tests produced incorrect results that were only resolved by increasing
statistical power (Section 4.14). We now treat minimum permutation counts
as hard requirements, not suggestions.

## 3.7 Phase-2a Standard Negative Panel

### Motivation

The cross-detection matrix (§3.4) provides broad coverage — every audited
cell is tested against every other model in the catalog. However, the
matrix has a structural limitation: it tests each detector only against
models already in the catalog, which were selected for their diversity of
emergent behaviors rather than for their diversity as non-positives. A
detector may pass the full cross-detection matrix while failing badly
against substrate types not represented in the catalog. Sprint 28's audit
documented that dim4 (broad negative sweep — testing specificity against
substrate-diverse non-positives) was the dominant coverage gap across 15
of 19 patterns. The cross-detection matrix is necessary but not sufficient
for specificity: it tests against known-pattern positives rather than
known-pattern negatives drawn from a broad substrate space, and a detector
that produces no false positives against its catalog neighbors can still
produce false positives at scale against novel substrates.

### Panel Composition

The Phase-2a standard negative panel consists of 30 substrates organized
into three classes of 10. **Class A** substrates are synthetic null
substrates, held fixed across all patterns: they include maximally
randomized grids, uniform-state configurations, temporally shuffled
trajectories, periodic synthetic signals, and other substrates constructed
to exercise the detector's primary and secondary metrics in the absence of
any emergent pattern. Class A is substrate-type-agnostic and provides a
universal baseline across all 19 detectors. **Class B** substrates are
catalog-derived non-positives, selected per substrate type under v1.1 of
the panel spec. For each pattern, Class B draws from models in the catalog
that share the same substrate type as the canonical positive but do not
exhibit the target pattern — for example, a lattice_2d detector is tested
against other lattice_2d models (Schelling, BTW, GoL, SIR, etc.) running
in parameter regimes where the target pattern is absent. The
substrate-typed selection ensures Class B exercises the detector against
the same observable types it was designed for, catching specificity
failures that cross-substrate tests miss. **Class C** substrates are
pattern-specific failed-regime substrates: parameter configurations of the
canonical positive model in which the target pattern is absent or
suppressed — for example, the Kuramoto model below K_c for P9
synchronization, or the BTW sandpile with dissipation for P14 SOC. Class
C has an N/A escape hatch introduced in v1.1: if the canonical model has
no parameter regime in which the target pattern is cleanly suppressed
(because the pattern is topologically enforced or structurally inherent),
the Class C slot is marked N/A with a logged explanation.

### PASS Criterion and Invariance Flags

The panel's primary criterion is overall true-negative rate (TNR) ≥ 0.95
across all 30 non-positive substrates. Per-class TNRs are additionally
reported, with an advisory weak-class threshold of 0.90 that triggers a
PASS-with-weakness verdict rather than a hard fail; a per-class TNR below
0.90 triggers FAIL regardless of overall. A secondary criterion, Cohen's
d ≥ 1.0 between the canonical-positive score distribution and the
pooled-negative score distribution, ensures that passing verdicts reflect
genuine distributional separation rather than fortuitously low
false-positive rates on the specific 30 substrates sampled. The
5-substrate gating floor makes per-class TNRs advisory below 5 substrates
(Class C N/A reduces the effective Class C size).

Version 1.2, introduced in Sprint 34, adds two invariance flags per
detector: `primary_metric_permutation_invariant` and
`primary_metric_time_shuffle_invariant`. For detectors whose primary
metric is invariant under spatial permutation or temporal shuffling of the
input, Class A substrates constructed by permutation or time-shuffling
will trivially produce the same primary metric value as the canonical
positive — making false-positive classification a certainty rather than a
failure of specificity. The v1.2 invariance flags allow the panel to skip
these degenerate-by-construction Class A substrates for the affected
detector, preventing spurious FAIL verdicts where the panel design itself
is the source of the false positive. The flags do not skip Class B or
Class C substrates, nor do they apply to detectors with
permutation-sensitive primary metrics.

### Spec Evolution as Methodology Contribution

The v1.0 → v1.1 → v1.2 trajectory of the Phase-2a panel spec is itself a
methodological story. The v1.0 panel was prototype-run against seven
patterns in Sprint 28 and surfaced two systemic issues: (1) cross-format
adapter contamination — Class B substrates were constructed with a
different history-format adapter than the canonical positive, causing
metric artifacts that appeared as false positives attributable to the
detector rather than to substrate construction; (2) substrate-type
conflation — Class B drew from the full catalog regardless of substrate
type, so a lattice_1d detector was tested against continuous_2d history
objects that its observable prerequisites would correctly reject. Sprint
30 resolved both issues in v1.1: Class B became substrate-typed and the
cross-format adapter was standardized. Sprint 32's P9 and P18 re-runs
under v1.1 then revealed a second generation of issues: both detectors
encountered Class A substrates (`constant_field` for P9, `all_same_state`
for P18) that their primary metrics flag as trivially maximal — not
because a pattern is present, but because the substrate is constructed to
be maximally uniform. These are degenerate-by-construction false positives
where the detector is not wrong to flag them but the panel spec is wrong
to present them as evidence of poor specificity. Sprint 34's v1.2
invariance flags resolve this by allowing affected detectors to annotate
which primary metrics are susceptible to degenerate Class A substrates, so
the panel can skip those slots rather than score them as failures. The
v1.0 → v1.1 → v1.2 evolution was not designed up-front; it emerged from
running the panel against real detectors and encountering failure modes
that the spec did not anticipate. This is an honest account of how a
methodology under active development refines itself through use.

### Reproducibility Note

The Phase-2a standard negative panel spec is versioned in
`docs/phase2a_panel_spec.md`. The
panel harness — substrate construction, detector dispatch, TNR computation,
and JSON output — lives in `epc/phase2a/`. Results from v1.0 and v1.1
runs are archived under `analysis/outputs/archive/`. Per-pattern panel
outputs for all completed runs are at
`analysis/outputs/p<i>_phase2a_panel.json`, where `<i>` is the pattern
number. Panel design decisions are recorded as Architecture Decision
Records in `REPLICATION_NOTES.md` alongside the per-sprint run sections.

### §3.5 Sprint 49 — Invariance-flag batch update

The Sprint 49 batch update applies the invariance-flag mechanism introduced in §3.4 (v1.2) to six additional detectors — P1, P2, P5, P6, P8, and P21 — whose Class A false-positive patterns were documented across Sprints 39–48. Each flag assignment is grounded in the mathematical structure of the primary metric rather than in ad-hoc per-detector tuning.

**Literature-first grounding.** Each flag pair (perm_inv, time_shuffle_inv) follows from the metric's definition: (i) *Permutation invariance* holds when the primary metric is a symmetric aggregate over agents or cells — e.g., P21's Hartigan dip test operates on sorted opinion values, so permuting agent indices leaves the test statistic unchanged; P2's `two_phase_score` measures cluster area fractions on a coarse density grid, which are also invariant to particle relabelling; P5's mean polar order parameter φ = |⟨e^iθ⟩| and P6's angular momentum |L| = |Σ r_i × v_i|/N are both sums over particles. (ii) *Time-shuffle invariance* holds when the primary metric is evaluated per-frame rather than over temporal trajectories — e.g., P1's Moran's I is computed on each lattice snapshot independently, P5's φ is evaluated per frame, P6's |L| per frame, and P8's `stopped_fraction` is a time-average statistic whose mean is preserved under any temporal reordering. In each case, the degenerate Class A substrate is *correctly* identified by the primary metric as anomalous, but the anomaly is structural (the substrate was constructed to preserve the metric) rather than reflecting a false positive in any scientifically meaningful sense. The flag annotates the latter class.

**P10 exception.** The P10 carry-forward C-p10-perm-shuffled-fp arises from a catalog-adapter artifact: the adapter binarizes continuous phases in a way that preserves bimodal structure for P10's chimera-state input, not from mathematical permutation invariance of `local_r`. The correct fix is an adapter refinement; auto-flagging permutation_invariant=True would be conceptually wrong. P10 flags remain (False, False) and the artifact is documented.

**Empirical outcomes (Sprint 49 re-runs).** Applying the six flag updates and re-running the affected panels produces the following advances: P2 and P6 reach TNR=1.000 (clean PASS); P5 reaches TNR=1.000 (clean PASS, up from PASS-with-weakness); P21 advances from PARTIAL to PASS-with-weakness (TNR 0.913→0.955; `time_shuffled` FP at CONFIRMATION persists — a convergence-timing issue, not a mathematical invariance, so `time_shuffle_invariant` is correctly left False for P21); P1 and P8 remain PARTIAL as expected (their residual FPs — C-p1-linear-gradient-fp and C-p8-class-c-near-onset — are Class C calibration issues outside the invariance-flag mechanism's scope). The AT-DEPTH count is unchanged at **10 / 19**.

### §3.6 Sprint 50 — P11 dim1 closure (Mobilia-Georgiev-Täuber 2007 reproduction)

Sprint 50 closes the sole remaining depth gap for P11 (predator-prey
oscillation / Lotka-Volterra) by numerically reproducing the O(1/L)
oscillation-amplitude scaling law from Mobilia, Georgiev & Täuber
(2007) *J. Stat. Phys.* 128, 447–483 (arXiv: q-bio/0512039). The
scaling law — that the amplitude of predator-prey density oscillations
decays as L^{−1} in finite stochastic lattice systems (Sec. III / Fig.
3 of the paper) — is the paper's primary quantitative result for the
coexistence (focus) phase. It arises from the van Kampen system-size
expansion of the master equation and is independent of the specific
reaction-rate ratios, making it reproducible at lattice sizes accessible
to our pure-Python implementation.

The reproduction runs at λ=4.0, σ=μ=1.0 (our canonical coexistence
parameters) across L ∈ {30, 50, 100, 150} with 3 seeds each and yields
scaling exponent −0.967 (R²=0.990), within the ±0.20 tolerance of the
published value −1.0 (relative error 3.3%). Coexistence is confirmed at
L=100 (5 seeds; mean predator density 0.081, FFT peak-to-mean ratio
48.9). A secondary finding — measured mean densities deviate
substantially from mean-field predictions (ρ_prey measured 0.589 vs MF
0.250) — is itself a published result of the paper (spatial correlations
dominate in the single-occupation lattice; Sec. III therein), providing
an additional qualitative confirmation of correctness.

**Dim1 reproduction table (cumulative through Sprint 50):**

| Pattern | Paper | Reproduced observable | Relative error | Sprint |
|---------|-------|-----------------------|---------------|--------|
| P11 LV | Mobilia-Georgiev-Täuber (2007) | Amplitude scaling exponent (−0.967 vs −1.0) | 3.3% | 50 |
| P28 YS | Chakraborti-Boghosian (2002/2014) | Gini convergence | <4% | 17 |
| P31 Zhang | Zhang et al. (2024) | Swap counts + insertion DG | <4% | 16 |
| P3 GS | Pearson (1993) | Turing wavelength (T=8000 regime) | N/A (detector-level) | 13 |

Patterns P2, P12, P21, P22 retain dim1 PARTIAL (no named figure reproduced
with stated tolerance). P11 dim1 PARTIAL→PASS advances P11 to **AT-DEPTH**
(11/19 patterns).

### §3.6 Sprint 51 — P22 dim1 closure (Datta-Acharyya 2021 reproduction)

Sprint 51 closes the P22 (information cascade / SIR) dim1 gap by reproducing
the wavefront speed measurement from Datta & Acharyya (2021) §3.1.1 / Fig. 11
(arXiv:2104.10456; *Int. J. Mod. Phys. C* 33, 2250094). The paper uses a
fixed-duration SIR CA (t_τ=4 steps per infection, per-neighbour infection
probability p0=0.25, Von Neumann neighbourhood, re-infection probability
p2=0.10) on a 500×500 lattice with a single seed at the centre. The most
specific numerical claim is the linear-fit slope of the wavefront radius
R(t): **0.4405 ± 0.0008 cells/step**.

Because `epc.models.sir_epidemic` uses stochastic geometric recovery rather
than the paper's fixed infection duration, the exact paper CA was implemented
inline in `analysis/reproductions/p22_dattaacharyya2005.py`. Running 20 seeds
on a 200×200 lattice (L=200 adequate; wavefront traverses only ~44 cells in
the 100-step fit window) yields:

**Dim1 reproduction table (cumulative through Sprint 51):**

| Pattern | Paper | Reproduced observable | Relative error | Sprint |
|---------|-------|-----------------------|---------------|--------|
| P22 SIR | Datta-Acharyya (2021) §3.1.1/Fig.11 | Wavefront speed (0.4612 vs 0.4405 cells/step) | 4.7% | **51** |
| P11 LV | Mobilia-Georgiev-Täuber (2007) | Amplitude scaling exponent (−0.967 vs −1.0) | 3.3% | 50 |
| P28 YS | Chakraborti-Boghosian (2002/2014) | Gini convergence | <4% | 17 |
| P31 Zhang | Zhang et al. (2024) | Swap counts + insertion DG | <4% | 16 |
| P3 GS | Pearson (1993) | Turing wavelength (T=8000 regime) | N/A (detector-level) | 13 |

Patterns P2, P12, P21 retain dim1 PARTIAL. P22 dim1 PARTIAL→PASS; dims 2–3
remain PARTIAL → P22 grade stays GAP. AT-DEPTH count unchanged at **11/19**.
