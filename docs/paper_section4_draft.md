# Section 4: Replication Studies

To validate the detection toolkit, we implemented seven canonical models
spanning five pattern clusters and verified that (a) each model reproduces
published quantitative results and (b) the corresponding detectors produce
correct tier assignments on both positive and negative controls. This section
reports results for each model family, followed by a consolidated cross-model
transfer matrix.

All models, metrics, and detectors are implemented in Python/NumPy without
heavy-framework dependencies. Source code, parameter files, and test suites
are available at [repository URL]. Every replication target is accompanied
by explicit pass/fail criteria defined before running the experiment. Where
results deviate from published values, we investigate and document the root
cause.

## 4.1 Zhang Cell-View Sorting (Clusters A, J)

**Reference:** Zhang, T., Goldstein, A. & Levin, M. (2024). Classical sorting
algorithms as a model of morphogenesis. *Adaptive Behavior*, 33, 25–54.

Zhang et al. reimplemented Bubble, Insertion, and Selection sort as
cell-autonomous agents operating on a shared 1D array (N=100), demonstrating
that decentralized sorting exhibits delayed gratification (P31) — temporary
decreases in local order that enable global progress — and that chimeric
arrays (mixed algotypes) spontaneously segregate by algorithm (P1).

**Model validation.** We provide two implementations: a deterministic
sequential-activation model for reproducible metric development, and a
faithful Python-threading model matching Zhang's architecture. Swap counts
match within 4% across all three algorithms (Bubble: 2,521 vs 2,449;
Insertion: 2,536 vs 2,483; Selection: 977 vs 1,096). Insertion DG matches
within 4% (1.06 vs 1.1). Bubble DG is consistently ~35% higher (0.33 vs
0.24), traced to an irreducible environmental difference: Linux/Python 3.12
thread scheduling produces different lock-contention patterns than Zhang's
macOS/Python 3.10 environment, directly affecting the monotonicity trajectory
shape that the DG metric measures.

**P31 detection.** The delayed gratification detector identifies DG in all
three algorithms. The critical validation is the non-redundancy test: does P31
capture information not already present in P1 (aggregation)? With properly
powered statistics (600 runs, 10-fold cross-validation), P31 survives
non-redundancy decisively (ΔR² = +0.645, p < 0.000001). P1 aggregation
features alone cannot distinguish algorithms (R² = −0.02 — all produce the
same sorted endpoint), while adding DG features explains 63% of sorting
efficiency variance. Ablation (shuffling DG features) destroys the signal
(R² = −0.03), confirming that temporal structure is the information-bearing
component.

**Methodological lesson.** Earlier underpowered tests (48–120 runs) produced
false negatives for the non-redundancy test. Statistical power requirements
for cross-validation with 8+ features demand ≥500 runs.

## 4.2 Greenberg-Hastings Excitable CA (Cluster D)

**References:** Greenberg, J.M. & Hastings, S.P. (1978). Spatial patterns for
discrete models of diffusion in excitable media. *SIAM Journal on Applied
Mathematics*, 34(3), 515–523.

The Greenberg-Hastings (GH) model is a κ-state excitable cellular automaton
where resting cells become excited if enough neighbors are excited, then
pass through refractory states before returning to rest. It produces spiral
waves and target patterns characteristic of P13 (excitable waves).

**Model validation.** Six canonical results replicated: winding number
conservation, wavefront propagation speed of exactly 1.0 cell/step, absence
of dispersion (speed independent of frequency), topological charge
conservation in periodic boundary conditions, spontaneous spiral
self-organization from random initial conditions, and threshold dependence
(waves extinguish above critical threshold).

**P13 detection.** The excitable wave detector reaches CONFIRMATION tier
(p = 0.005, 199 null trajectories) on both spiral and random-initialization
runs. Wavefront speed CV < 0.15 and persistent spiral tip trajectories
satisfy the confirmation criteria.

## 4.3 Conway's Game of Life (Cluster E)

**Reference:** Gardner, M. (1970). The fantastic combinations of John Conway's
new solitaire game "life." *Scientific American*, 223(4), 120–123.

The Game of Life (GoL) is the canonical model for P15 (persistent propagating
computation). Gliders and other structures propagate, collide, and produce
input-dependent outputs — the hallmark of computational dynamics.

**Model validation.** Cell-by-cell match against independent reference
implementations for all standard patterns: still lifes (block, beehive,
loaf), oscillators (blinker, toad, beacon), glider displacement c/4,
LWSS displacement c/2, and R-pentomino evolution to generation 1103 with
final population 116 — exact match with LifeWiki reference.

**P15 detection.** The P15 discriminator uses a two-stage approach to
distinguish computational dynamics (GoL) from merely excitable dynamics (GH):

*Stage 1 — Transfer Entropy.* Boundary-conditioned TE cleanly separates
the two systems: GoL produces TE ratios of 15–16× above permutation null,
while GH produces ratios of 1–2×. Critically, raw (non-boundary-conditioned)
average TE gives the wrong ordering (GH > GoL), because GH's deterministic
wave propagation creates trivial TE that dominates the average. Restricting
to boundary cells — where structures interact — isolates the computational
signal. This is a novel methodological contribution.

*Stage 2 — Functional test.* Deterministic replay of glider collisions at
different relative phases produces reproducibility = 1.000 (identical initial
conditions always produce identical outcomes) with 2 distinct coarse outcome
types (annihilation and still life). The combination of perfect
reproducibility with input-dependent outcome variation demonstrates
functional computation, not merely stochastic dynamics.

**P1 false positive — resolved.** GoL produces genuine spatial
autocorrelation from B3/S23 dynamics (Moran's I is significant vs
label-shuffle null), initially triggering P1. This was resolved through a
two-layer defense: (1) a type constancy check — GoL alive/dead states change
every step and therefore are not persistent agent type labels, making P1
structurally inapplicable; (2) substrate-aware dispatch that verifies the
model provides a 'types' observable before running P1.

## 4.4 Vicsek Flocking and D'Orsogna Milling (Cluster B)

**References:** Vicsek, T. et al. (1995). Novel type of phase transition in
a system of self-driven particles. *Physical Review Letters*, 75(6),
1226–1229. D'Orsogna, M.R. et al. (2006). Self-propelled particles with
soft-core interactions. *Physical Review Letters*, 96, 104302.

**Vicsek model validation (P5).** The standard Vicsek model — N point
particles with constant speed, angular noise, and local alignment
interaction — replicates all six published claims: kinetic phase transition
from disorder (φ = 0.051 ≈ 1/√300) to order (φ = 0.9995 at η = 0.1),
continuous transition controlled by noise η, density dependence (higher
density promotes order), and absence of rotational order (|L| = 0.04,
negative control for P6).

**D'Orsogna model validation (P6).** The D'Orsogna self-propelled particle
model with Morse attraction-repulsion potential produces stable milling
at published parameters: mean speed 1.396 ≈ v_eq = √(α/β) = 1.414
(99% match), angular momentum |L| = 0.996 (coherent rotation),
polarization φ = 0.008 (no net translation), and ring density profile
with empty core (hollowness = 0.000).

**Cross-detection transfer matrix.**

|               | Vicsek (η=0.5)                   | D'Orsogna (milling)                |
|---------------|----------------------------------|-------------------------------------|
| P5 (flocking) | **DEFINITIVE** (φ=0.988, p=0.005)| not detected (φ=0.008)             |
| P6 (milling)  | not detected (|L|=0.017)         | **DEFINITIVE** (|L|=0.996, p=0.005)|

Perfect discrimination: no false positives and no missed detections. The
heading-shuffle null uses random uniform headings rather than permutation of
existing headings — permuting near-identical headings in an ordered flock
preserves φ ≈ 1 and gives uninformative p-values (architecture decision #20).

## 4.5 Kuramoto Coupled Oscillators (Cluster C)

**Reference:** Kuramoto, Y. (1975). Self-entrainment of a population of
coupled non-linear oscillators. *International Symposium on Mathematical
Problems in Theoretical Physics*, Lecture Notes in Physics 39, 420–422.

**Model validation.** The all-to-all Kuramoto model with Lorentzian
frequency distribution (γ = 0.5) replicates all six published claims:
disordered baseline r = 0.071 = 1/√200 (ratio 1.00 exact), phase
transition onset at K_c = 1.0 for N ≥ 300 matching the analytical
K_c = 2γ, ordered-state r(K) matching the analytical r = √(1 − K_c/K)
within 0.04 at K = 1.2, 2.0, and 4.0, and frequency entrainment
with 97% of oscillators locked at K = 6K_c.

An initial report of a "27% K_c shift" was traced to a threshold artifact
on a coarse parameter grid, not a genuine discrepancy. Fine-resolution
scanning at N = 500 confirmed r(K = 1.2) matches analytical predictions
within 0.039.

**P9 detection.** The synchronization detector reaches DEFINITIVE tier
on Kuramoto at K = 8K_c (r = 0.963, 119 oscillation periods, p = 0.005,
199 random-phase null comparisons). Below K_c, the detector correctly
reports no detection (r = 0.087, p = 0.185). P10 (chimera states) is
excluded via uniform local order parameter check (CV = 0.037, indicating
globally uniform synchronization rather than coexisting synchronized
and desynchronized domains).

**KSG Transfer Entropy.** The KSG estimator for continuous-variable
TE — required because plug-in discretization destroys phase information —
validates on analytical ground truths (Gaussian mutual information error
0.013 nats) and correctly shows TE increasing with Kuramoto coupling
strength (p = 0.02 at K = 2K_c and K = 6K_c).

## 4.6 Schelling Segregation (Cluster A)

**Reference:** Schelling, T.C. (1971). Dynamic models of segregation.
*Journal of Mathematical Sociology*, 1(2), 143–186.

**Model validation.** The minimal Schelling model (50×50 grid, 90% density,
threshold 3/8) converges from random initial placement (Moran's I = 0.005)
to a segregated steady state (I = 0.415) within approximately 20 steps.

**P1 detection.** The full P1 aggregation detector reaches CONFIRMATION tier
(p = 0.001, 999 label-shuffle permutations, confidence 0.700). Moran's
I = 0.423 against a null mean of −0.0002 produces Cohen's d = 49.87.
Segregation index = 0.652 (well above the 0.4 confirmation threshold) with
sustained_i CV = 0.000 (perfect plateau stability). The temporal convergence
guard confirms genuine aggregation: I_initial = 0.005 → I_late = 0.415,
with the guard passing via has_gain AND is_plateaued (Schelling converges
fast then plateaus rather than showing monotonic increase throughout).

CONFIRMATION rather than DEFINITIVE is the correct tier: with 999
shuffle permutations, the floor p-value is 1/1000 = 0.001, which does not
satisfy the strict p < 0.001 definitive criterion. Reaching DEFINITIVE
would require either ≥1999 permutations or a mechanistic null model.

**Negative controls.** Random grids (same composition, no dynamics) produce
p = 0.452 (not confirmed). GoL is correctly rejected by the type constancy
check (alive/dead states are not persistent agent types).

## 4.7 BTW Sandpile (Cluster D)

**Reference:** Bak, P., Tang, C. & Wiesenfeld, K. (1987). Self-organized
criticality: An explanation of the 1/f noise. *Physical Review Letters*,
59(4), 381–384.

**Model validation.** The 2D BTW sandpile (L = 64, 100,000 driving events
after 10,000 burn-in) reaches the critical state with all heights below z_c
(max z = 3 = z_c − 1, mean height 2.098). Avalanche sizes span 4.3 decades
(range [1, 20,972]) with heavy-tailed distribution (mean/median = 12.1).

**P14 detection.** The SOC detector reaches DEFINITIVE tier (confidence 0.850).
The power-law exponent τ = 1.247 (MLE with x_min = 1) matches the published
2D BTW value of τ ≈ 1.20 within 0.05, cross-validated by log-binned PDF
slope (τ = 1.241). Power-law is strongly preferred over exponential
(likelihood ratio R = +80.6, p < 0.001). Duration scaling T ~ s^γ with
γ = 0.642 provides a passing secondary metric.

The dissipative null model (p_diss = 0.2, bulk grain loss per toppling)
produces exponentially distributed avalanches (LR R = −6.0, max size 68
vs 20,972) and is correctly not detected as SOC.

**Resolved caveat.** Log-normal is preferred over simple power-law
(R = −76.2) — a known property of the 2D BTW universality class, which
exhibits multifractal scaling with logarithmic corrections rather than a
clean simple power law. The 1/f noise signature was initially unmeasurable
(β = −0.17) because the spectral analysis used avalanche sizes (approximately
IID). Measuring the power spectral density of total energy E(t) = Σh(x,t)
instead yields β = 1.41, correctly in the 1/f range [0.5, 1.5]. This
confirms the expected temporal correlations in the SOC state.

## 4.8 Nowak-May Spatial Prisoner's Dilemma (Cluster H)

The Nowak-May spatial PD (1992) places cooperators and defectors on a lattice
where each cell plays a one-shot PD with its Moore neighbors and imitates the
highest-payoff neighbor synchronously. With benefit parameter b=1.8 (the
chaotic coexistence regime), the model produces fractal-like C/D boundaries
with cooperation sustained at f_C ≈ 0.41 through spatial clustering (Moran's
I = 0.497). This confirms that cooperation persistence arises from spatial
reciprocity rather than individual strategy.

The P27 detector achieves DEFINITIVE tier at b=1.8, with the Moran's I
permutation test (199 permutations, p=0.005) confirming that the spatial
clustering of cooperators is significantly stronger than random. PD payoff
structure is verified (T>R>P≥S). Negative controls: b=2.0 (C extinct, no
cooperation to cluster) and b=1.0 (all C, no dilemma) both correctly produce
no detection. The use of `coop_fraction` rather than `grid` as the required
observable prevents false compatibility with other lattice_2d models (GH, GoL,
Schelling) that produce spatial structure for unrelated reasons.

## 4.9 Hegselmann-Krause Bounded-Confidence Opinion Dynamics (Cluster F)

The Hegselmann-Krause model (2002) places N agents with continuous opinions in
[0,1] and updates each to the mean opinion of all agents within confidence
bound ε. At ε=0.2, the initially uniform distribution polarizes into 2 stable
clusters — the classic polarization outcome. At ε=0.1, fragmentation produces
4 clusters; at ε=0.05, 7 clusters; at ε=0.5, consensus (1 cluster).

The P21 detector achieves DEFINITIVE tier at ε=0.2 using Hartigan's dip test
(bootstrap with 1000 samples, p=0.001) to confirm multimodality in the final
opinion distribution, combined with verification that the initial condition was
unimodal (ruling out pre-existing structure). The persistence check confirms
clusters are stable for ≥10 steps. P18 (consensus) exclusion ensures
detection only fires when the population genuinely splits rather than converging.

## 4.10 Consolidated Transfer Matrix

The following matrix summarizes detection outcomes across all model-detector
pairs tested through Sprint 5. Entries show the highest achieved tier
(D = definitive, C = confirmation, S = screening, rej = rejected by guard,
— = untested compatible pair, × = substrate mismatch).

|                  | P1   | P5   | P6   | P9   | P13  | P14  | P15     | P21  | P27  | P31  |
|------------------|------|------|------|------|------|------|---------|------|------|------|
| Zhang sorting    | S    | ×    | ×    | ×    | ×    | ×    | ×       | ×    | ×    | C    |
| Schelling        | C    | ×    | ×    | ×    | ×    | ×    | ×       | ×    | ×    | —    |
| Vicsek (ordered) | ×    | **D**| —    | ×    | ×    | ×    | ×       | ×    | ×    | ×    |
| D'Orsogna (mill) | ×    | —    | **D**| ×    | ×    | ×    | ×       | ×    | ×    | ×    |
| Kuramoto (sync)  | ×    | ×    | ×    | **D**| ×    | ×    | ×       | ×    | ×    | ×    |
| GH spiral        | S    | ×    | ×    | ×    | C    | ×    | P13     | ×    | ×    | —    |
| GoL R-pentomino  | rej  | ×    | ×    | ×    | —    | ×    | **D**   | ×    | ×    | —    |
| BTW sandpile     | ×    | ×    | ×    | ×    | ×    | **D**| —       | ×    | ×    | ×    |
| Nowak-May (1.8)  | ×    | ×    | ×    | ×    | ×    | ×    | ×       | ×    | **D**| ×    |
| HK (ε=0.2)      | ×    | ×    | ×    | ×    | ×    | ×    | ×       | **D**| ×    | ×    |

Key observations:

1. **No false positives across substrate boundaries.** The substrate-aware
   dispatch system (5 substrate types: lattice_1d, lattice_2d, continuous_2d,
   oscillator, opinion_space) correctly prevents cross-substrate detector
   application. Of 110 possible model × detector pairs, 24 are compatible
   and 86 are correctly rejected.

2. **Clean within-cluster discrimination.** P5/P6 show perfect discrimination
   between flocking and milling. P13/P15 show clean separation between
   excitable waves and computation via the two-stage TE + functional
   discriminator. P1 correctly fires on aggregation models (Schelling,
   chimeric sorting) and is correctly rejected on GoL (type constancy guard).

3. **Observable-level filtering beyond substrate.** P27 requires `coop_fraction`
   (not just `grid`), preventing false matches with GH, GoL, and Schelling
   despite shared lattice_2d substrate. P14 requires `avalanche_sizes`,
   restricting to BTW. This demonstrates that substrate type alone is
   insufficient — fine-grained observable matching is essential.

4. **Tier ceilings are meaningful.** Several detections are capped at
   CONFIRMATION (Schelling P1, GH P13) rather than DEFINITIVE. In each case,
   the ceiling arises from specific methodological constraints (permutation
   count floor for Schelling, absence of mechanistic null for GH) rather
   than weak signals. This demonstrates that the tier system correctly
   distinguishes evidence strength from effect magnitude.

5. **Negative results are informative.** GoL P1 "rejection" required a
   two-layer defense (type constancy + substrate check) to resolve a genuine
   conceptual false positive. The resolution deepened our understanding of
   what P1 actually detects: persistent-type-label spatial sorting, not
   merely spatial autocorrelation from any source.

## 4.11 Methodological Lessons

Several methodological insights emerged from the replication process:

**Statistical power.** Underpowered tests produced incorrect results in
multiple cases: P31 non-redundancy with 48 runs gave false negatives that
were only resolved at 600 runs; P5 on Vicsek required 5000+ recorded steps
for definitive detection; KSG TE at T = 2000 produced false positive
(p = 0.03) on independent data, requiring T ≥ 5000 for correct non-significance.
We now enforce minimum power requirements for all permutation-based tests
(≥99 permutations for p = 0.01 floor, ≥199 for p = 0.005).

**Test correctness.** In multiple cases, initial test failures were caused by
errors in the test, not the code: P15 fidelity used overly fine-grained
outcome signatures; P15 collision test placed gliders too far apart to
collide within the measurement window; P1 guard sampled from step 7 instead
of step 0, missing fast Schelling convergence. The lesson: when a test fails,
the first response should be to verify the test is correct, not to weaken
the assertion.

**Boundary-conditioned measurement.** The P13/P15 TE discriminator required
restricting measurement to boundary cells rather than averaging over the
full lattice. Raw average TE gives the wrong ordering (GH > GoL) because
deterministic wave propagation creates trivially high TE at interior cells.
Boundary conditioning isolates the interaction signal that distinguishes
computational from merely excitable dynamics.

**Intrinsic timescales.** D'Orsogna milling detection failed when using the
initial particle spread to estimate the crossing time T_cross — the mill
compacts significantly during formation (from radius 5.0 to diameter 3.02).
Using the measured group diameter after equilibration gave the correct
timescale. Similarly, Kuramoto entrainment must measure frequency dispersion
among locked oscillators only; Lorentzian tail outliers that never lock
dominate the all-oscillator statistic.
