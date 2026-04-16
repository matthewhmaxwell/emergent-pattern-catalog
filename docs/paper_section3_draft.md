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

The first level is **substrate type matching.** Five substrate types are
defined: lattice_1d (chimeric sorting), lattice_2d (cellular automata,
Schelling, Nowak-May), continuous_2d (Vicsek, D'Orsogna), oscillator
(Kuramoto), and opinion_space (Hegselmann-Krause). Each detector declares
which substrate types it is compatible with.

The second level is **observable matching.** Within a compatible substrate,
the detector requires specific state keys. P27 (spatial reciprocity)
requires the `coop_fraction` observable that only Prisoner's Dilemma models
produce, preventing it from firing on other lattice_2d systems like
Greenberg-Hastings or Game of Life. P14 (self-organized criticality)
requires `avalanche_sizes`, restricting it to sandpile-type models.

Of 110 possible model × detector pairs in the current inventory
(11 models × 10 detectors), 24 pass both checks. The remaining 86 are
correctly rejected — 82 by substrate mismatch and 4 by missing observables.
This block-diagonal structure ensures that detectors operate only in domains
where their metrics are meaningful.

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

## 3.6 Statistical Power Requirements

Permutation-based significance testing imposes minimum sample sizes that we
enforce across the toolkit. With n permutations, the smallest achievable
p-value is 1/(n + 1). For a significance threshold of p < 0.01, at least
99 permutations are required; for p < 0.005, at least 199. Our standard
configuration uses 999 permutations for P1 (yielding p floor of 0.001),
199 for P5/P6/P9/P27 (p floor 0.005), and 99 for TE-based tests (p floor
0.01). The P31 non-redundancy test requires at least 500 model runs for
stable 10-fold cross-validation with 8+ features.

These are not arbitrary choices — in multiple cases, underpowered initial
tests produced incorrect results that were only resolved by increasing
statistical power (Section 4.11). We now treat minimum permutation counts
as hard requirements, not suggestions.
