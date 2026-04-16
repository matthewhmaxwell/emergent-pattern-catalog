# Section 5: Cross-Model Transfer

The detection toolkit's value lies not only in confirming that canonical
models produce their expected patterns, but in revealing which patterns
appear — and which do not — when detectors are applied systematically
across the full model inventory. This section reports the completed
transfer matrix and analyzes the cross-model structure it reveals.

## 5.1 The Completed Transfer Matrix

We tested all 24 compatible model × detector pairs identified by the
substrate-aware dispatch system (Section 3.4). Of 110 total cells in the
11 × 10 matrix, 86 are correctly eliminated by substrate or observable
mismatch, leaving 24 cells requiring empirical evaluation. Two additional
cells (Nowak-May × P15 and Schelling × P15) are substrate-compatible but
untestable due to an implementation limitation: the P15 detector currently
uses an internal Game of Life stepper rather than accepting external model
histories.

The completed matrix (Table 1) contains 10 definitive detections
(one per canonical model-pattern pair), 3 confirmations, 2 screening-level
detections, 1 rejection by guard, and 6 correct non-detections on cross-
substrate negative controls. Every canonical positive model reaches at
least confirmation tier on its primary detector, and every cross-pattern
negative control produces correct rejection.

**Table 1: Completed Transfer Matrix**

|                  | P1   | P5   | P6   | P9   | P13  | P14  | P15  | P21  | P27  | P31  |
|------------------|------|------|------|------|------|------|------|------|------|------|
| Zhang sorting    | S    | ×    | ×    | ×    | ×    | ×    | ×    | ×    | ×    | C    |
| Schelling        | C    | ×    | ×    | ×    | rej  | ×    | †    | ×    | ×    | —    |
| Vicsek (ordered) | ×    | D    | rej  | ×    | ×    | ×    | ×    | ×    | ×    | ×    |
| D'Orsogna (mill) | ×    | rej  | D    | ×    | ×    | ×    | ×    | ×    | ×    | ×    |
| Kuramoto (sync)  | ×    | ×    | ×    | D    | ×    | ×    | ×    | ×    | ×    | ×    |
| GH spiral        | S    | ×    | ×    | ×    | C    | ×    | P13  | ×    | ×    | —    |
| GoL R-pentomino  | rej  | ×    | ×    | ×    | rej  | ×    | D    | ×    | ×    | —    |
| BTW sandpile     | ×    | ×    | ×    | ×    | ×    | D    | †    | ×    | ×    | ×    |
| Nowak-May (1.8)  | C    | ×    | ×    | ×    | rej  | ×    | †    | ×    | D    | ×    |
| HK (ε=0.2)      | ×    | ×    | ×    | ×    | ×    | ×    | ×    | D    | ×    | ×    |

D = definitive, C = confirmation, S = screening, rej = rejected by detector
guard, × = substrate/observable mismatch, † = P15 implementation-limited,
P13 = classified as excitable (not computational) by TE discriminator,
— = untested (P31 requires 1D substrate).

## 5.2 Block-Diagonal Structure

The matrix is approximately block-diagonal by substrate type, with each
detector firing only within its target substrate. This structure is
enforced by the dispatch system but validated empirically: within compatible
substrate blocks, detectors produce correct positive detections on canonical
models and correct rejections on non-canonical models.

The lattice_2d block is the densest, containing 5 models (Greenberg-Hastings,
Game of Life, BTW sandpile, Schelling, Nowak-May) and 5 potentially
applicable detectors (P1, P13, P14, P15, P27). Observable-level filtering
further restricts this block: P14 fires only on BTW (requires avalanche
data), P27 fires only on Nowak-May (requires cooperation fraction), and P15
is implementation-restricted to GoL.

The continuous_2d block (Vicsek, D'Orsogna) demonstrates clean cross-pattern
discrimination: P5 fires definitively on ordered Vicsek (φ = 0.988) and
correctly rejects milling D'Orsogna (φ = 0.046), while P6 fires definitively
on milling D'Orsogna (|L| = 0.996) and correctly rejects ordered Vicsek
(|L| = 0.031). These cross-exclusion results validate the detector design
at the metric level — the rejections are not merely substrate mismatches
but genuine false negative verifications where the wrong pattern's primary
metric falls far below its screening threshold.

## 5.3 Co-Occurrence: Nowak-May Aggregation

The most scientifically interesting transfer result is the detection of P1
(aggregation) on Nowak-May at confirmation tier (Moran's I = 0.898,
segregation index = 0.753, p = 0.005). Cooperator and defector clusters
in the spatial Prisoner's Dilemma exhibit the same spatial autocorrelation
signature as Schelling segregation, even though the underlying mechanism
is payoff imitation rather than preference-driven relocation.

This co-occurrence is expected: P1's design allows co-occurrence with P27
(spatial reciprocity), recognizing that cooperator clustering in spatial
games is simultaneously an instance of aggregation (like types collect)
and spatial reciprocity (clustering sustains cooperation). The detection
confirms that these are genuinely overlapping patterns, not competing
descriptions of the same phenomenon. A system can exhibit both, and the
detectors correctly identify both.

Importantly, the P1 detector reaches confirmation rather than definitive
tier on Nowak-May. This reflects a genuine conceptual nuance: Nowak-May
agents do not "move" to aggregate — they change strategy in place through
imitation. The aggregation arises from strategy dynamics on a fixed
substrate rather than from spatial relocation. The confirmation tier
correctly signals that the aggregation signature is real and statistically
significant while leaving room for the mechanistic distinction between
Schelling-style movement-based aggregation and Nowak-May-style
imitation-based clustering.

## 5.4 Guard-Based Rejections

Three cross-model tests produced rejections via detector guards rather than
simple metric failures. These are instructive because they reveal cases
where a naive metric application would produce misleading results.

**GoL × P1 (type constancy guard).** Game of Life generates significant
spatial autocorrelation from B3/S23 dynamics. Without the type constancy
guard, P1 would detect aggregation — correctly, at the statistical level —
on a system that has no persistent agent types to aggregate. The guard
catches this by verifying that the values at each spatial position are
stable over the measurement window. GoL alive/dead states flip every step,
failing this check.

**GoL × P13 (excitable medium guard).** GoL's synchronous binary update
produces wavefront speed statistics that pass P13's screening threshold
(speed CV = 0.0, since every cell has identical inter-excitation intervals).
The n_states ≥ 3 guard prevents this false positive by enforcing that the
model has the minimum state complexity for excitable dynamics
(resting/excited/refractory).

**Nowak-May × P13 (speed instability).** Even after passing the n_states
guard (which would correctly reject Nowak-May's 2-state system), the
wavefront speed CV of 0.47 — far above the 0.2 screening threshold —
independently confirms that Nowak-May's strategy dynamics do not produce
coherent wavefront propagation. Multiple independent rejection mechanisms
increase confidence in the result.

## 5.5 Remaining Gaps

Three categories of gaps remain in the matrix.

**P31 cross-tests.** P31 (delayed gratification) requires lattice_1d
substrate, limiting it to Zhang sorting variants. Testing P31 on other
systems would require either 1D variants of existing models (Schelling on a
ring, GH on a 1D chain) or generalizing the DG metric to 2D substrates.

**P15 generalization.** The P15 detector's internal GoL stepper prevents
cross-testing against other lattice_2d models. Generalizing P15 to accept
external model histories — and defining what "persistent propagating
computation" means for Schelling or Nowak-May — would require either a
general-purpose functional test framework or model-specific collision
libraries.

**Missing pattern clusters.** The current inventory covers 7 of 10 pattern
clusters. Clusters G (Resilience: P24–P26), I (Structure Formation:
P29–P30), and parts of Cluster F (P22 information cascades, P23
anti-coordination) have no implemented models or detectors. Expanding
coverage to these clusters would require new models (SIR epidemics, ant
trail formation, Physarum networks) and corresponding detection metrics.

## 5.6 Dimensional Coverage

The 11 implemented models span 8 of 11 ontological dimensions with at least
two distinct values each. Coverage gaps concentrate in three dimensions:
interaction type (no indirect-stigmergic models, such as ant trail
formation), memory (no environmental trace models), and external driving
(all models are autonomous). These gaps correspond precisely to the missing
pattern clusters (G and I), suggesting that filling in the model inventory
will naturally expand dimensional coverage.

Within the covered dimensions, the transfer matrix demonstrates that
detection is robust across dimensional variation. P1 fires correctly on
both lattice (Schelling) and game-theoretic (Nowak-May) substrates. P5/P6
discriminate correctly despite identical substrate, interaction type, and
spatial scale — the distinguishing dimension is the pattern's geometric
character (translational vs. rotational alignment). P13 and P15 share
substrate and observable type, requiring a cross-dimensional test (transfer
entropy) to separate them.
