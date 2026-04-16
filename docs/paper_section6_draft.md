# Section 6: Emergent Findings

The systematic application of quantitative detectors across the model
inventory produced several findings that were not anticipated at the outset.
These fall into three categories: validation of a contested pattern,
a novel methodological technique, and unexpected cross-model results that
sharpen the definitions of existing patterns.

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
