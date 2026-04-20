# Section 4: Replication Studies

To validate the detection toolkit, we implemented 14 canonical model
files covering 13 distinct model families (Zhang cell-view sorting has
sequential and threaded variants of the same model) spanning eight
pattern clusters and verified that (a) each model
reproduces published quantitative results and (b) the corresponding
detectors produce correct tier assignments on both positive and
negative controls. This section reports results for each model family,
followed by a consolidated cross-model transfer matrix with 37 audited
model × detector pairs.

All models, metrics, and detectors are implemented in Python/NumPy
without heavy-framework dependencies. Source code, parameter files, and
test suites are available at [repository URL]. Every replication target
is accompanied by explicit pass/fail criteria defined before running
the experiment. Where results deviate from published values, we
investigate and document the root cause. Where replication required
reshaping the detector around an unexpected empirical finding — as
occurred in Sprints 9, 10, and 11 — the finding is reported explicitly
rather than quietly absorbed.

## 4.1 Zhang Cell-View Sorting (Clusters A, J)

**Reference:** Zhang, T., Goldstein, A. & Levin, M. (2024). Classical
sorting algorithms as a model of morphogenesis. *Adaptive Behavior*,
33, 25–54.

Zhang et al. reimplemented Bubble, Insertion, and Selection sort as
cell-autonomous agents operating on a shared 1D array (N=100),
demonstrating that decentralized sorting exhibits delayed gratification
(P31) — temporary decreases in local order that enable global progress
— and that chimeric arrays (mixed algotypes) spontaneously segregate by
algorithm (P1).

**Model validation.** We provide two implementations: a deterministic
sequential-activation model for reproducible metric development, and a
faithful Python-threading model matching Zhang's architecture. Swap
counts match within 4% across all three algorithms (Bubble: 2,521 vs
2,449; Insertion: 2,536 vs 2,483; Selection: 977 vs 1,096). Insertion
DG matches within 4% (1.06 vs 1.1). Bubble DG is consistently ~35%
higher (0.33 vs 0.24), traced to an irreducible environmental
difference: Linux/Python 3.12 thread scheduling produces different
lock-contention patterns than Zhang's macOS/Python 3.10 environment,
directly affecting the monotonicity trajectory shape that the DG metric
measures.

**P31 detection.** The delayed gratification detector identifies DG in
all three algorithms. The critical validation is the non-redundancy
test: does P31 capture information not already present in P1
(aggregation)? With properly powered statistics (600 runs, 10-fold
cross-validation), P31 survives non-redundancy decisively (ΔR² =
+0.645, p < 0.000001). P1 aggregation features alone cannot distinguish
algorithms (R² = −0.02 — all produce the same sorted endpoint), while
adding DG features explains 63% of sorting efficiency variance.
Ablation (shuffling DG features while preserving marginal movement
statistics) destroys the signal (R² = −0.03), confirming that temporal
structure is the information-bearing component.

**Methodological lesson.** Earlier underpowered tests (48–120 runs)
produced false negatives for the non-redundancy test. Statistical power
requirements for cross-validation with 8+ features demand ≥500 runs.

## 4.2 Greenberg-Hastings Excitable CA (Cluster D)

**References:** Greenberg, J.M. & Hastings, S.P. (1978). Spatial
patterns for discrete models of diffusion in excitable media. *SIAM
Journal on Applied Mathematics*, 34(3), 515–523.

The Greenberg-Hastings (GH) model is a κ-state excitable cellular
automaton where resting cells become excited if enough neighbors are
excited, then pass through refractory states before returning to rest.
It produces spiral waves and target patterns characteristic of P13
(excitable waves).

**Model validation.** Six canonical results replicated: winding number
conservation, wavefront propagation speed of exactly 1.0 cell/step,
absence of dispersion (speed independent of frequency), topological
charge conservation in periodic boundary conditions, spontaneous spiral
self-organization from random initial conditions, and threshold
dependence (waves extinguish above critical threshold).

**P13 detection.** The excitable wave detector reaches CONFIRMATION
tier (p = 0.005, 199 null trajectories) on both spiral and
random-initialization runs. Wavefront speed CV < 0.15 and persistent
spiral tip trajectories satisfy the confirmation criteria.

## 4.3 Conway's Game of Life (Cluster E)

**Reference:** Gardner, M. (1970). The fantastic combinations of John
Conway's new solitaire game "life." *Scientific American*, 223(4),
120–123.

The Game of Life (GoL) is the canonical model for P15 (persistent
propagating computation). Gliders and other structures propagate,
collide, and produce input-dependent outputs — the hallmark of
computational dynamics.

**Model validation.** Cell-by-cell match against independent reference
implementations for all standard patterns: still lifes (block, beehive,
loaf), oscillators (blinker, toad, beacon), glider displacement c/4,
LWSS displacement c/2, and R-pentomino evolution to generation 1103
with final population 116 — exact match with LifeWiki reference.

**P15 detection.** The P15 discriminator uses a two-stage approach to
distinguish computational dynamics (GoL) from merely excitable dynamics
(GH):

*Stage 1 — Transfer Entropy.* Boundary-conditioned TE cleanly
separates the two systems: GoL produces TE ratios of 15–16× above
permutation null, while GH produces ratios of 1–2×. Critically, raw
(non-boundary-conditioned) average TE gives the wrong ordering (GH >
GoL), because GH's deterministic wave propagation creates trivial TE
that dominates the average. Restricting to boundary cells — where
structures interact — isolates the computational signal. This is a
novel methodological contribution, discussed further in Section 6.

*Stage 2 — Functional test.* Deterministic replay of glider collisions
at different relative phases produces reproducibility = 1.000
(identical initial conditions always produce identical outcomes) with
2 distinct coarse outcome types (annihilation and still life). The
combination of perfect reproducibility with input-dependent outcome
variation demonstrates functional computation, not merely stochastic
dynamics.

**P1 false positive — resolved.** GoL produces genuine spatial
autocorrelation from B3/S23 dynamics (Moran's I is significant vs
label-shuffle null), initially triggering P1. This was resolved through
a two-layer defense: (1) a type constancy check — GoL alive/dead states
change every step and therefore are not persistent agent type labels,
making P1 structurally inapplicable; (2) substrate-aware dispatch that
verifies the model provides a `types` observable before running P1.

## 4.4 Vicsek Flocking and D'Orsogna Milling (Cluster B)

**References:** Vicsek, T. et al. (1995). Novel type of phase
transition in a system of self-driven particles. *Physical Review
Letters*, 75(6), 1226–1229. D'Orsogna, M.R. et al. (2006).
Self-propelled particles with soft-core interactions. *Physical Review
Letters*, 96, 104302.

**Vicsek model validation (P5).** The standard Vicsek model — N point
particles with constant speed, angular noise, and local alignment
interaction — replicates all six published claims: kinetic phase
transition from disorder (φ = 0.051 ≈ 1/√300) to order (φ = 0.9995 at
η = 0.1), continuous transition controlled by noise η, density
dependence (higher density promotes order), and absence of rotational
order (|L| = 0.04, negative control for P6).

**D'Orsogna model validation (P6).** The D'Orsogna self-propelled
particle model with Morse attraction-repulsion potential produces
stable milling at published parameters: mean speed 1.396 ≈ v_eq =
√(α/β) = 1.414 (99% match), angular momentum |L| = 0.996 (coherent
rotation), polarization φ = 0.008 (no net translation), and ring
density profile with empty core (hollowness = 0.000).

**Cross-detection transfer matrix.**

|               | Vicsek (η=0.5)                   | D'Orsogna (milling)                |
|---------------|----------------------------------|-------------------------------------|
| P5 (flocking) | **DEFINITIVE** (φ=0.988, p=0.005)| not detected (φ=0.008)             |
| P6 (milling)  | not detected (|L|=0.017)         | **DEFINITIVE** (|L|=0.996, p=0.005)|

Perfect discrimination: no false positives and no missed detections.
The heading-shuffle null uses random uniform headings rather than
permutation of existing headings — permuting near-identical headings in
an ordered flock preserves φ ≈ 1 and gives uninformative p-values.

## 4.5 Kuramoto Coupled Oscillators (Cluster C)

**Reference:** Kuramoto, Y. (1975). Self-entrainment of a population of
coupled non-linear oscillators. *International Symposium on
Mathematical Problems in Theoretical Physics*, Lecture Notes in Physics
39, 420–422.

**Model validation.** The all-to-all Kuramoto model with Lorentzian
frequency distribution (γ = 0.5) replicates all six published claims:
disordered baseline r = 0.071 = 1/√200 (ratio 1.00 exact), phase
transition onset at K_c = 1.0 for N ≥ 300 matching the analytical K_c
= 2γ, ordered-state r(K) matching the analytical r = √(1 − K_c/K)
within 0.04 at K = 1.2, 2.0, and 4.0, and frequency entrainment with
97% of oscillators locked at K = 6K_c.

An initial report of a "27% K_c shift" was traced to a threshold
artifact on a coarse parameter grid, not a genuine discrepancy.
Fine-resolution scanning at N = 500 confirmed r(K = 1.2) matches
analytical predictions within 0.039.

**P9 detection.** The synchronization detector reaches DEFINITIVE tier
on Kuramoto at K = 8K_c (r = 0.963, 119 oscillation periods, p =
0.005, 199 random-phase null comparisons). Below K_c, the detector
correctly reports no detection (r = 0.087, p = 0.185). P10 (chimera
states) is excluded via uniform local order parameter check (CV =
0.037, indicating globally uniform synchronization rather than
coexisting synchronized and desynchronized domains).

**KSG Transfer Entropy.** The KSG estimator for continuous-variable
TE — required because plug-in discretization destroys phase information
— validates on analytical ground truths (Gaussian mutual information
error 0.013 nats) and correctly shows TE increasing with Kuramoto
coupling strength (p = 0.02 at K = 2K_c and K = 6K_c).

## 4.6 Schelling Segregation (Cluster A)

**Reference:** Schelling, T.C. (1971). Dynamic models of segregation.
*Journal of Mathematical Sociology*, 1(2), 143–186.

**Model validation.** The minimal Schelling model (50×50 grid, 90%
density, threshold 3/8) converges from random initial placement
(Moran's I = 0.005) to a segregated steady state (I = 0.415) within
approximately 20 steps.

**P1 detection.** The full P1 aggregation detector reaches CONFIRMATION
tier (p = 0.001, 999 label-shuffle permutations, confidence 0.700).
Moran's I = 0.423 against a null mean of −0.0002 produces Cohen's d =
49.87. Segregation index = 0.652 (well above the 0.4 confirmation
threshold) with sustained_i CV = 0.000 (perfect plateau stability).
The temporal convergence guard confirms genuine aggregation: I_initial
= 0.005 → I_late = 0.415, with the guard passing via has_gain AND
is_plateaued (Schelling converges fast then plateaus rather than
showing monotonic increase throughout).

The P1 primary metric is `morans_i_final` (the Moran's I at the end of
the trajectory), not peak or sustained Moran's I. This choice was
established empirically in Sprint 10 after a six-model characterization
revealed that peak Moran's I cannot distinguish transient wavefronts
from sustained spiral clustering (Section 4.10). On Schelling, peak and
final agree (the system converges quickly to a plateau), so the metric
choice is immaterial; the distinction matters elsewhere.

CONFIRMATION rather than DEFINITIVE is the correct tier: with 999
shuffle permutations, the floor p-value is 1/1000 = 0.001, which does
not satisfy the strict p < 0.001 definitive criterion. Reaching
DEFINITIVE would require either ≥1999 permutations or a mechanistic
null model.

**Negative controls.** Random grids (same composition, no dynamics)
produce p = 0.452 (not confirmed). GoL is correctly rejected by the
type constancy check (alive/dead states are not persistent agent
types).

## 4.7 BTW Sandpile (Cluster D)

**Reference:** Bak, P., Tang, C. & Wiesenfeld, K. (1987).
Self-organized criticality: An explanation of the 1/f noise. *Physical
Review Letters*, 59(4), 381–384.

**Model validation.** The 2D BTW sandpile (L = 64, 100,000 driving
events after 10,000 burn-in) reaches the critical state with all
heights below z_c (max z = 3 = z_c − 1, mean height 2.098). Avalanche
sizes span 4.3 decades (range [1, 20,972]) with heavy-tailed
distribution (mean/median = 12.1).

**P14 detection.** The SOC detector reaches DEFINITIVE tier (confidence
0.850). The power-law exponent τ = 1.247 (MLE with x_min = 1) matches
the published 2D BTW value of τ ≈ 1.20 within 0.05, cross-validated
by log-binned PDF slope (τ = 1.241). Power-law is strongly preferred
over exponential (likelihood ratio R = +80.6, p < 0.001). Duration
scaling T ~ s^γ with γ = 0.642 provides a passing secondary metric.

The dissipative null model (p_diss = 0.2, bulk grain loss per
toppling) produces exponentially distributed avalanches (LR R = −6.0,
max size 68 vs 20,972) and is correctly not detected as SOC.

**Resolved caveat.** Log-normal is preferred over simple power-law (R
= −76.2) — a known property of the 2D BTW universality class, which
exhibits multifractal scaling with logarithmic corrections rather than
a clean simple power law. The 1/f noise signature was initially
unmeasurable (β = −0.17) because the spectral analysis used avalanche
sizes (approximately IID). Measuring the power spectral density of
total energy E(t) = Σh(x,t) instead yields β = 1.41, correctly in the
1/f range [0.5, 1.5]. This confirms the expected temporal correlations
in the SOC state.

## 4.8 Nowak-May Spatial Prisoner's Dilemma (Cluster H)

The Nowak-May spatial PD (1992) places cooperators and defectors on a
lattice where each cell plays a one-shot PD with its Moore neighbors
and imitates the highest-payoff neighbor synchronously. With benefit
parameter b=1.8 (the chaotic coexistence regime), the model produces
fractal-like C/D boundaries with cooperation sustained at f_C ≈ 0.41
through spatial clustering. This confirms that cooperation persistence
arises from spatial reciprocity rather than individual strategy.

The P27 detector achieves DEFINITIVE tier at b=1.8, with the Moran's I
permutation test (199 permutations, p=0.005) confirming that the
spatial clustering of cooperators is significantly stronger than
random. PD payoff structure is verified (T>R>P≥S). Negative controls:
b=2.0 (C extinct, no cooperation to cluster) and b=1.0 (all C, no
dilemma) both correctly produce no detection. The use of
`coop_fraction` rather than `grid` as the required observable prevents
false compatibility with other lattice_2d models (GH, GoL, Schelling)
that produce spatial structure for unrelated reasons.

**P1 co-occurrence.** Nowak-May also fires P1 at CONFIRMATION tier
(Moran's I_final = 0.49 at 100×100, segregation_index = 0.70). This
co-occurrence is expected — cooperator clusters exhibit the same
spatial autocorrelation signature as Schelling segregation, even
though the underlying mechanism is payoff imitation rather than
preference-driven relocation. P1 and P27 mark each other as
co-occurrence candidates rather than mutually excluding.

## 4.9 Hegselmann-Krause Bounded-Confidence Opinion Dynamics (Cluster F)

The Hegselmann-Krause model (2002) places N agents with continuous
opinions in [0,1] and updates each to the mean opinion of all agents
within confidence bound ε. At ε=0.2, the initially uniform distribution
polarizes into 2 stable clusters — the classic polarization outcome.
At ε=0.1, fragmentation produces 4 clusters; at ε=0.05, 7 clusters; at
ε=0.5, consensus (1 cluster).

The P21 detector achieves DEFINITIVE tier at ε=0.2 using Hartigan's
dip test (bootstrap with 1000 samples, p=0.001) to confirm
multimodality in the final opinion distribution, combined with
verification that the initial condition was unimodal (ruling out
pre-existing structure). The persistence check confirms clusters are
stable for ≥10 steps. P18 (consensus) exclusion ensures detection only
fires when the population genuinely splits rather than converging.

## 4.10 SIR Epidemic CA and P22 Information Cascade (Cluster F)

**Primary reference:** Datta, A. & Acharyya, M. (2022). Modelling the
spread of an epidemic in presence of vaccination using cellular
automata. *International Journal of Modern Physics C* 33, 2250094
(arXiv:2104.10456).

An earlier version of this work cited Fuks & Lawniczak (2002) as the
primary reference for the SIR CA; we corrected this in Sprint 8 after
discovering that their paper describes a lattice gas cellular
automaton in which individuals physically move between sites — a
fundamentally different model from our fixed-site SIR. The correct
primary reference is Datta & Acharyya, who use the same fixed-site
three-state automaton with independent-neighbor infection probability
that our implementation realizes.

**Model validation.** Four replication targets were tested on 100×100
to 201×201 lattices. *Linear wavefront growth:* across six parameter
combinations (Moore and Von Neumann neighborhoods, p ∈ [0.2, 0.5], q ∈
[0.05, 0.1]), the wavefront radius R(t) grows linearly with R² > 0.98
at all configurations, confirming the characteristic spatial spreading
behavior. *Percolation transition:* sharp critical thresholds at p_c ≈
0.10 (Von Neumann) and p_c ≈ 0.038 (Moore), with 50% percolation
observed at the critical point. *Kermack-McKendrick curve shape:*
unimodal I(t) with exactly one sign change in dI/dt, peak infection
fraction 23% of grid interior, final I = 0. *Exact conservation:* S(t)
+ I(t) + R(t) = N holds to integer precision at every timestep across
300-step runs.

The mean-field R₀ approximation R₀ ≈ (1 − (1−p)ⁿ)/q substantially
overestimates the true lattice reproductive number, giving R₀_approx ≈
3.4 (Von Neumann) and 2.7 (Moore) at the measured critical thresholds
rather than the mean-field prediction R₀_c = 1. This discrepancy is
the expected consequence of spatial correlation: each infected cell
depletes its local susceptible pool before the wavefront propagates,
a well-known feature of spatial epidemic models (Grassberger, 1983).

**P22 detection.** SIR achieves DEFINITIVE tier on P22 (information
cascade) at 80×80, p=0.30, q=0.20, single-seed initialization:
Moran's I of the infection time-map I(t_peak) = 0.987 (cells that
became infected near-simultaneously are spatially clustered, the
signature of a coherent wavefront), Cohen's d = 109.5 versus
spatial-shuffle null, p = 0.005 (n_permutations = 199). Negative
controls on GoL (R-pentomino reach < 5% of grid), Nowak-May (no
wavefront; Moran's I of infection time ≈ 0), Schelling (no wavefront;
Moran's I of time ≈ 0.06), and GH spiral (reentrant, not cascade;
rejected at screening) all produce correct non-detections.

**P13 exclusion.** SIR and GH share the same substrate (lattice_2d
with n_states = 3), but P13 correctly rejects SIR. The hard
prerequisite n_states ≥ 3 passes for SIR, but the screening metric —
wavefront speed stability — fails because SIR's single-pass dynamics
produce no re-excitation: the wavefront propagates once, cells become
immune (R), and activity dies before the 5×T_prop persistence
threshold is reached. P13's discrimination between GH and SIR is the
core scientific finding of the Sprint 7–8 replication: even with three
discrete states and visually similar wavefronts, re-entrant excitable
dynamics and single-pass epidemic dynamics produce quantitatively
different detector outputs.

**P1 asymmetry — Sprint 10 finding.** An initial version of the P1
detector used the maximum of peak and final Moran's I as the primary
metric. Under that definition, SIR passed P1 at screening tier because
its infection wavefront produces peak Moran's I = 0.89 — a strong
spatial clustering signal, equal to or exceeding Schelling's plateau
value. A six-model empirical characterization revealed that this was
classifying a transient signal as a sustained one:

| Model              | I_peak | I_final | I_sustained | sustained CV |
|--------------------|--------|---------|-------------|--------------|
| Schelling τ=0.375  | +0.414 | +0.414  | +0.414      | 0.00         |
| Nowak-May b=1.8    | +0.794 | +0.530  | +0.500      | 0.07         |
| SIR β=0.20 γ=0.3   | +0.892 | **+0.019** | +0.175   | **0.99**     |
| RPS M=1e-4         | +0.582 | +0.550  | +0.562      | 0.016        |
| GH n=8 random      | +0.412 | +0.204  | +0.204      | 0.00         |
| Random grid (noise)| +0.028 | +0.015  | −0.007      | inf          |

SIR's I_peak of 0.89 occurs because the infection wavefront
transiently clusters infected cells; by the end of the run, the
population has converged to near-uniform recovered (R) state and the
clustering has vanished (I_final = 0.02, sustained CV = 0.99
indicating wild variation in the sustained window). Changing the
primary metric to `morans_i_final` (Decision 32) flips SIR × P1 from
`screening` to `rejected` while leaving Schelling, Nowak-May, and GH
unchanged. The same change preserves the RPS × P1 detection (Section
4.11) because RPS has I_peak ≈ I_final ≈ 0.55 (sustained spiral
clustering). A concurrent bugfix added a 0.05 magnitude floor on
I_final (Decision 33), fixing a pre-existing false positive on pure
noise (random grid I = 0.015 was passing the trivial `I > expected_I`
check via sampling variance).

This is a worked example of the "look before touching" principle: the
sprint plan initially proposed a peak→sustained swap, which the
characterization showed would not have distinguished SIR from RPS
(sustained = 0.175 vs 0.562 — both nonzero, both ambiguous). The
peak→final swap emerged from the empirical data, not from the plan.

## 4.11 Spatial Rock-Paper-Scissors and P12 Cyclic Dominance (Cluster C)

**Primary reference:** Reichenbach, T., Mobilia, M. & Frey, E. (2007).
Mobility promotes and jeopardizes biodiversity in rock-paper-scissors
games. *Nature* 448, 1046–1049 (arXiv:q-bio/0702032).

The spatial RPS model places three species {A, B, C} on a 2D lattice
with three reactions: selection (σ, A beats B, B beats C, C beats A,
loser dies leaving empty site), reproduction (µ, species expands into
adjacent empty site), and exchange (ε, two adjacent individuals swap).
Reichenbach et al. show that biodiversity — the ability of all three
species to coexist — depends critically on mobility M = 2εa²/N: below a
critical mobility M_c ≈ 4.5 × 10⁻⁴, spiral patterns stabilize
coexistence; above M_c, mobility destroys the spatial structure and
drives two species to extinction.

**Model validation.** Three replication targets were tested. *Biodiversity
regimes:* At L = 50, M = 10⁻⁵ produces deep coexistence (all species >
10% density after 200 generations, mean fractions near 1/3 each); M =
10⁻² ≈ 20 M_c produces deep extinction (minimum species fraction drops
below 5% within 500 generations, maximum exceeds 70%); minimum
species fraction decreases monotonically with mobility across M ∈
{10⁻⁵, 10⁻³, 5 × 10⁻³}. This matches the qualitative phase diagram
structure of Reichenbach Fig. 4, though we do not pin M_c precisely
(that would require a fine mobility sweep with ≥20 seeds per point,
roughly an hour of compute per data point on our grid). *Conservation:*
total site count L² is preserved at every snapshot. *Dominance-only
selection:* in a striped A/B lattice with ε = 0 and µ ≈ 0, A cells
(which nothing dominates in this two-species subset) are never killed;
B cells drop >50% in 20 generations.

We did not replicate the spiral-wavelength scaling law λ ∝ √M. A full
replication would require Fourier analysis or spiral-tip tracking at
multiple mobilities with long runs, and the resulting fit constant
depends on finite-size effects not controlled in our small L ≤ 60 test
configurations. The P12 detector does not key off wavelength; the more
diagnostic observable is the neighbor-conditional transition ratio
described below.

**P12 detection.** The P12 cyclic-dominance detector's primary metric is
the `intransitivity_score`, defined as log₁₀ of the minimum
neighbor-conditional replacement ratio over cyclic triples of species:

ρ(X,Y) = P(cell → Y | had Y-neighbor) / P(cell → Y | no Y-neighbor),
intransitivity_score = log₁₀ max_{triples} min_edge ρ(X,Y).

A large intransitivity_score means that replacements of species X by
Y occur strongly preferentially when Y-neighbors are present, across
the full cyclic triple — the spatial fingerprint of cyclic dominance.
At the canonical positive (L=30, M=10⁻⁴, 80 generations, n_perm=199),
RPS achieves CONFIRMATION (intransitivity_score = 1.83, which
corresponds to ρ_min ≈ 67.6; coexistence_fraction = 1.0;
direction_stable = True; identified_triple = [1,3,2] matching the
model's dominance map; null p-value = 0.005 at the permutation floor;
Cohen's d > 200 against a spatial-shuffle null with mean ≈ 0.01, std
≈ 0.01). At a stronger positive (L=40, M=10⁻⁵, 150 generations,
n_perm=499), RPS achieves DEFINITIVE (score = 2.06, p = 0.002).

**Negative controls.** GH (n_states = 3 and n_states = 5, threshold =
1): intransitivity_score = 0.0 exactly (ρ = 1 for all edges —
clock-driven transitions have no species preference). SIR
(infection_prob = 0.5, recovery_prob = 0.1): score = 0.0. GoL and
Nowak-May: rejected by prerequisite (n_candidate_species = 2 < 3).

**P13 boundary test — the Sprint 9 scientific headline.** The Sprint 9
prompt predicted that P13 might false-positive on RPS, since RPS
satisfies P13's obvious prerequisites: n_states ≥ 3 (passes hard
guard), persistent wavefronts (not died_out), cells are re-excited
(unlike SIR). We measured: P13 cleanly rejects RPS at screening across
three mobilities and three seeds. The rejection is on
`wavefront_speed_cv`, which lands in the range [0.59, 0.68] on RPS
versus ≈ 0.05–0.15 on GH — far exceeding the 0.2 screening threshold.

The mechanism is diagnostic. Excitable media have clock-driven
transitions: once excited, a cell ticks through refractory → rest on a
deterministic schedule, producing uniform wavefront speed. RPS
transitions are entirely neighbor-driven and stochastic: a cell
changes state only when a specific-species neighbor is selected by the
asynchronous scheduler, producing wavefronts that look spiral-like on
visual inspection but have much greater per-cell speed variability. The
CV metric captures this precisely. P12 and P13 therefore agree on every
tested case through independent mechanisms: P12 detects RPS and finds
P13 excluded (ρ_min > 10 rules out the clock-driven mechanism), while
P13 independently rejects RPS via the CV test.

## 4.12 Lotka-Volterra Lattice and P11 Predator-Prey Oscillation (Cluster C)

**Primary reference:** Mobilia, M., Georgiev, I.T. & Täuber, U.C.
(2007). Phase transitions and spatio-temporal fluctuations in
stochastic lattice Lotka-Volterra models. *Journal of Statistical
Physics* 128, 447–483 (arXiv:q-bio/0512039).

The stochastic lattice Lotka-Volterra model in the single-occupation
variant places prey (B) and predator (A) on a 2D lattice with three
reactions: predator death (A → ∅, rate µ), prey reproduction into an
empty neighbor (B + ∅ → B + B, rate σ), and predation (A + B → A + A,
rate λ). The single-occupation constraint means each lattice site
holds at most one individual, so the empty state ∅ is a third site
label and prey + predator densities do not sum to unity — a feature
that becomes crucial in the P11 detector design below.

**Parameter regime finding.** The Mobilia paper reports canonical
parameters λ = 0.2, σ = µ = 0.1 at L = 512. Our first attempt used
λ = 2.0, σ = µ = 1.0 at L = 100 — rates with identical ratios to the
published values, differing only in overall time-unit scaling.
Predators went extinct within 100–200 generations at these parameters,
a finite-size effect at L = 100 that shrinks with L. A λ sweep at
σ = µ = 1 and L = 100 revealed three regimes:

| λ (σ=µ=1, L=100) | Regime              | Predator dynamics                            |
|------------------|---------------------|----------------------------------------------|
| 2.0, 2.5         | Extinction          | All runs extinct by t ~ 100–200 generations  |
| 3.0 – 4.5        | Coexistence (focus) | Sustained population with erratic oscillations |
| 5.0              | Stationary nodes    | Clusters persist but no oscillation          |

Our canonical coexistence choice is λ = 4.0, σ = µ = 1.0, L = 100,
stable across seeds 42, 7, and 123 for ≥1200 generations with no
extinctions observed.

**Oscillation characterization.** Short runs (200 generations) show a
single deterministic-like initial swing (predator crash + prey
recovery) followed by noisy quasi-stationary fluctuations. This is
misleading: it looks like no persistent oscillation. Long runs (2000
generations) reveal the actual signal:

| L   | λ | Prey std (SS) | Pred std (SS) | FFT peak-to-mean | Dominant period |
|-----|---|---------------|---------------|------------------|-----------------|
| 30  | 4 | 0.113         | 0.029         | 22.3             | ~100 gens       |
| 100 | 4 | 0.034         | 0.009         | 22.6             | ~80 gens        |
| 128 | 4 | 0.019         | 0.005         | 22.4             | ~350 gens       |

Two facts emerge. First, oscillation *amplitudes* shrink as L
increases at the expected rate O(1/√N) of resonant demographic
fluctuations, matching Mobilia Fig. 3. Second, despite shrinking
amplitudes, the FFT peak-to-mean ratio stays robustly above 10 across
the L range tested — the *detection signal* is scale-robust even though
the amplitude is not. This distinction mattered for the detector
design: an amplitude-based metric would have been fragile, but a
spectral-concentration metric is robust.

**P11 detector — two empirical course corrections.** The Sprint 11
prompt proposed "max cross-correlation at positive lag τ > 0, with
predator lagging prey by quarter period" as the P11 primary metric,
following the textbook classical-Lotka-Volterra phase relationship. Two
measurements on the spatial stochastic model overrode this.

*Correction 1.* At L = 100, λ = 4, seed = 42, the cross-correlation at
τ = +40 (approximately the expected T/4 quarter-period lag) is only
+0.16 — a weak positive-lag signal. But the cross-correlation at τ =
−15 is **−0.85** — a strong anti-phase coupling at short lag. The
quarter-period phase shift combined with strong anti-phase coupling
means that the minimum of the cross-correlogram occurs at |τ| ~ 10–20
rather than at τ = T/4 ~ 40; the anti-phase coupling dominates over the
phase-lag signature at finite amplitudes. We redesigned the primary
metric as

ρ_anti = min_{|τ| ≥ 5} Pearson(prey(t), pred(t + τ))

with measured range −0.72 to −0.88 across seeds 42, 7, and 123 at the
canonical positive. The |τ| ≥ 5 floor was added as a safeguard
(Correction 2 below).

*Correction 2.* Before locking the detector, we ran a broad
negative-model sweep to test for false positives. The Nowak-May
cooperator/defector time series produced ρ_anti = −0.979 at lag +3 — a
stronger anti-correlation than any LV configuration. The cause is
algebraic, not dynamical: Nowak-May has cooperator_fraction +
defector_fraction = 1 exactly (every cell is one or the other, no
empty state), so the two fractions are identically anti-correlated
with correlation −1 at *every* lag by conservation. The detector was
correctly identifying strict anti-correlation but was doing so in a
model where that anti-correlation reflects conservation, not
predator-prey coupling. We added a prerequisite: `std(species_A +
species_B) > 0.005` — a nontrivial empty reservoir is required. LV has
std(prey + predator) ≈ 0.034; Nowak-May has exactly 0.000. The
|τ| ≥ 5 floor in the primary metric is a defense-in-depth against the
same class of artifact (conservation-driven correlations appear at
very short lags).

**Null model and tier gating.** The P11 null is a circular shift of
one species' time series, which preserves each series' autocorrelation
and FFT magnitude spectrum while destroying the cross-series phase
relationship. This null turned out to be intentionally too strong to
produce a clean p-value: at the canonical LV positive, the null
ρ_anti distribution has mean −0.46 and 5th percentile −0.869
(essentially equal to the observed −0.863), producing a one-sided p of
about 0.07. The reason is that LV's slow-mode autocorrelation survives
circular shifting and occasionally produces deep anti-correlations in
the null distribution. This is not a bug — it correctly reflects that
"oscillating series are autocorrelated" is insufficient evidence for a
predator-prey coupling. But it means the p-value cannot be a clean
tier gate. Our resolution: P11 does not gate tier assignment on
null_p (Decision 36). Discrimination relies on ρ_anti magnitude (|LV|
~ 0.8, |noise| ~ 0.1, null mean ~ 0.4) and Cohen's d of the observed
ρ_anti against the null distribution (|LV| between −1.75 and −2.21
across seeds). The p-value is reported as a diagnostic field but is
not used in the tier decision.

**Tier assignment outcomes.** At L = 100, λ = 4, seeds {42, 7, 123},
1200–1500 generations with burn-in 100:

| Seed | n_steps | ρ_anti | |τ_anti| | FFT p2m | Cohen's d | Tier        |
|------|---------|--------|----------|---------|-----------|-------------|
| 42   | 1500    | −0.863 | 16       | 25.1    | −2.21     | DEFINITIVE  |
| 42   | 1200    | −0.863 | 16       | 21.6    | −1.75     | DEFINITIVE  |
| 7    | 1200    | −0.716 | 16       | 13.6    | −1.77     | DEFINITIVE  |
| 123  | 1200    | −0.780 | 14       | 17.6    | −1.95     | DEFINITIVE  |

**Negative controls — bilateral-vs-cyclic separation.** The most
important finding of the Sprint 11 negative sweep is that spatial RPS
produces a *stronger* ρ_anti than Lotka-Volterra on any pair of species:

| Model             | ρ_anti | |τ| | FFT p2m | species std       | total std | Tier                            |
|-------------------|--------|-----|---------|-------------------|-----------|---------------------------------|
| LV seed=42        | −0.86  | 16  | 25.1    | 0.037, 0.009      | 0.034     | **DEFINITIVE**                  |
| LV seed=7         | −0.72  | 16  | 13.6    | 0.024, —          | 0.022     | **DEFINITIVE**                  |
| LV seed=123       | −0.78  | 14  | 17.6    | 0.027, —          | 0.024     | **DEFINITIVE**                  |
| RPS (A vs B)      | −0.94  | 28  | 77      | large             | ≠ 0       | rejected (n_species=3)          |
| RPS (B vs C)      | −0.94  | 35  | 89      | large             | ≠ 0       | rejected (n_species=3)          |
| Nowak-May         | −0.98  | 3   | 12.6    | 0.26, 0.26        | **0.000** | rejected (total_std prereq)     |
| Schelling τ=0.375 | 0.00   | —   | 0.0     | **0.000**, 0.000  | 0.000     | rejected (species_std prereq)   |
| SIR (post-BI)     | 0.00   | —   | 1.9     | **0.0002**, 0.0002| 0.000     | rejected (species_std prereq)   |
| White noise       | −0.08  | —   | 2.7     | —                 | —         | rejected (ρ_anti screen)        |

Three species pairs from RPS produce ρ_anti ≈ −0.94, stronger than LV
on the primary metric alone. What keeps P11 specific to *bilateral*
predator-prey coupling rather than firing on cyclic-dominance systems
is not a signal-level discriminator but a prerequisite gate:
`n_unique_species_observed == 2`. This is the single most important
design-level finding of Sprint 11. A detector keyed off ρ_anti alone
would false-positive on RPS, SIR (after the empty-reservoir
prerequisite), and any other multi-species anti-correlated system.
The n_species prerequisite converts the detector from a signal
classifier to a substrate classifier on the n_species axis.

## 4.13 Gray-Scott Reaction-Diffusion and P3 Turing Wavelength (Cluster D)

**Primary reference:** Pearson, J.E. (1993). Complex patterns in a
simple system. *Science* 261, 189–192.

The Gray-Scott model is a two-species reaction-diffusion system on a
continuous lattice: an activator u and inhibitor v satisfy
∂_t u = Dᵤ∇²u − uv² + F(1 − u) and ∂_t v = Dᵥ∇²v + uv² − (F + k)v, with
feed rate F and kill rate k. Depending on (F, k), solutions form
stationary spots, stripes, labyrinths, or self-replicating blobs. In
the detection framework, Gray-Scott is the first model outside
integer-grid cellular automata: it exposes a continuous `field`
observable (float-valued v in [0, 1]) rather than a categorical `grid`.
This introduces a new substrate type, `lattice_2d_continuous`, and
requires a detector whose primary signal is a spectral property of the
continuous field rather than a categorical adjacency statistic.

**Canonical positive characterization.** Pearson's Fig. 1 panel "spots"
is frequently cited with (F, k) = (0.062, 0.0609). At N = 128, dt = 1,
Dᵤ = 0.16, Dᵥ = 0.08 — grid-scale coefficients chosen for numerical
stability (dt · max(Dᵤ, Dᵥ) ≤ 0.25) — this parameter pair produces a
pattern with characteristic wavelength ≈ 64 pixels, half the domain.
At N = 256 the same parameters resolve into proper spots with
wavelength ≈ 12 pixels, matching the published figure. The N = 128
behavior is therefore a domain-size artifact, not a genuine spot regime,
and was unsuitable for routine testing at practical grid sizes. The
canonical positive was instead chosen to be the labyrinth regime (F =
0.037, k = 0.060), which produces wavelength ≈ 12 pixels already at
N = 128 and remains invariant to N = 64 (peak wavenumber k = 5, 7, 10,
16 at N = 64, 96, 128, 192; wavelength in pixels ≈ 12 throughout). The
pattern selection transient is long — roughly 3500–4000 timesteps at
N = 128 — so canonical runs are T = 4000.

**P3 detector — primary metric.** P3 measures the concentration of
spectral power at a non-trivial wavenumber via the radial-averaged 2D
FFT power spectrum. For a v-field snapshot V of shape (N, N), let
S(k) = ⟨|F[V](k_x, k_y)|²⟩ be the power averaged over annuli in k-space
at radius k. The primary metric is peak-to-mean = max_k S(k) /
mean_k S(k) over k in [k_min = 2, N/2]. At the canonical labyrinth
positive, peak-to-mean = 18.75 across seeds {42, 7, 123}, with Cohen's
d against the spatial-shuffle null ≈ 103. The shuffle null randomizes
pixel positions, preserving the one-point distribution of v-values
while destroying all spatial correlation; under the null, peak-to-mean
collapses to 1.0 ± 0.04.

**The RPS false-positive trap.** Before locking P3's primary
threshold, we ran the radial-FFT on every existing integer-grid model
with the same pipeline. Every model produced peak-to-mean near 1 —
except RPS at low mobility, which produced peak-to-mean ≈ 23.10,
numerically exceeding Gray-Scott's 18.75. The cause is mechanistic:
RPS spiral domains have a characteristic size set by the cyclic
reaction-diffusion wavelength λ ∝ √M, and the raw integer grid, when
treated as a float array, has enough spectral concentration at the
spiral wavenumber to spoof Turing-pattern-like peak-to-mean. This was
the load-bearing Sprint 13 finding: any P3 variant gated only on
peak-to-mean would false-positive on RPS, producing a silent failure
mode where a cyclic-dominance system masquerades as Turing pattern
formation.

**Substrate-level discrimination.** The resolution was to move
discrimination from empirical thresholding to substrate-level content
gates. P3 requires, as prerequisites, both (a) the presence of a
`field` observable — which RPS exposes only `grid` and so fails
immediately, and (b) n_unique_values ≥ 50 on the field snapshot. At
N = 128, Gray-Scott's continuous field has ≥ 16,000 distinct float
values; every integer-grid model has at most the number of discrete
states (typically 2–10). This pair of prerequisites separates the two
classes adversarially: even when RPS is deliberately re-labeled as a
`field` (the `TestAdversarialDiscreteFieldRejected` case), it rejects
at n_unique = 4 < 50 regardless of its peak-to-mean = 23.10. The
discriminator is therefore content-level, not signal-level, and this
is the architecturally cleaner pattern: where two canonical models
produce overlapping signal on the intended primary metric, look for a
content-level property that separates them before reaching for fragile
thresholds.

**Null model and tier gating.** The P3 null permutes pixel positions
within each v-field snapshot, computing peak-to-mean under the null
for 199 permutations. The lower bound on the null-p achievable at
n = 199 is 1/(199 + 1) = 0.005, so the p < 0.01 CONFIRMATION gate is
achievable but the p < 0.001 DEFINITIVE gate is not. Since P3's effect
sizes are enormous (Cohen's d ≈ 100), the extra permutations above 99
are effectively free at 1–2 seconds per detector call. Tier thresholds:
peak-to-mean > 5.0 for SCREENING; peak-to-mean > 10.0 with d > 10 and
peak-k coefficient-of-variation < 0.15 and null-p < 0.01 for
CONFIRMATION; peak-to-mean > 15.0 with d > 30 and peak-k CV < 0.05 for
DEFINITIVE.

**Canonical-positive outcomes.** At F = 0.037, k = 0.060, N = 128,
T = 4000, seeds {42, 7, 123}:

| Seed | peak-to-mean | peak-k | wavelength (px) | Cohen's d | peak-k CV | Tier        |
|------|--------------|--------|-----------------|-----------|-----------|-------------|
| 42   | 18.75        | 10     | 12.8            | 103       | 0.00      | DEFINITIVE  |
| 7    | 18.50        | 10     | 12.8            | 102       | 0.00      | DEFINITIVE  |
| 123  | 19.01        | 10     | 12.8            | 104       | 0.00      | DEFINITIVE  |

Grid-scaling: DEFINITIVE is reached at N = 64 (T = 3000), N = 96
(T = 3000), and N = 128 (T = 4000). The short-wavelength-spots regime
(F = 0.030, k = 0.062) reaches CONFIRMATION at peak-to-mean = 13.35,
below the DEFINITIVE threshold of 15.0 but above CONFIRMATION's 10.0;
whether to lower the DEFINITIVE threshold to admit both regimes as
DEFINITIVE remains a future-work question, since doing so would narrow
the discrimination margin against spurious signals.

**Transfer-matrix additions.** Gray-Scott × P3 is the Sprint 13
positive (DEFINITIVE). Every integer-grid model × P3 rejects at the
`field` prerequisite (rej; seven cells). Every 2D-grid-consuming
detector × Gray-Scott rejects at its own `grid` prerequisite (rej; P1,
P11, P12, P13, P22 — five cells). P15 × Gray-Scott is substrate-
compatible but not detected (nd; Gray-Scott is deterministic with no
stochastic step_fn for functional replay). The Gray-Scott × P1 cell
required a small Sprint 14 hardening: P1's 2D branch previously
raised `KeyError` on missing `grid` / `type_labels_at_pos` observables,
and was updated to return a graceful substrate-warning rejection
matching the pattern used by P11, P13, and P22. Thirteen new cells
(one D + seven rej + five rej + zero nd) join the audited transfer
matrix, bringing the total to 63 audited cells.

## 4.14 Nagel-Schreckenberg Traffic CA and P8 Traffic Jamming (Cluster D)

**Primary reference:** Nagel, K. & Schreckenberg, M. (1992). A cellular
automaton model for freeway traffic. *J. Phys. I France* 2, 2221–2229.

**Secondary reference:** Bette, H.M., Habel, L., Emig, T. &
Schreckenberg, M. (2017). Mechanisms of jamming in the Nagel-
Schreckenberg model for traffic flow. *Phys. Rev. E* 95, 012311.

The Nagel-Schreckenberg model is a minimal one-dimensional cellular
automaton for single-lane freeway traffic. L cells are arranged on a
ring (periodic boundaries); each cell is either empty or occupied by
one car; each car carries an integer velocity v ∈ {0, 1, …, v_max}.
Each time step is a parallel update applying four rules in order:
(1) accelerate v ← min(v + 1, v_max); (2) slow down to the gap
v ← min(v, d) where d is the number of empty cells to the next car
ahead; (3) randomize v ← max(v − 1, 0) with probability p; (4) move
x ← (x + v) mod L. Despite its extreme simplicity — no continuous
dynamics, no forces, one parameter each for maximum velocity and
driver imperfection — the model reproduces the textbook freeway flow-
density fundamental diagram and the spontaneous formation of traffic
jams from homogeneous initial conditions. Nagel-Schreckenberg occupies
the `lattice_1d` substrate, previously held only by Zhang cell-view
sorting, and is the first `lattice_1d` model in the catalog with an
integer velocity observable.

**Fundamental-diagram replication.** At L = 1000, v_max = 5, p = 0.3 —
Nagel-Schreckenberg's original illustrative parameter choice — our
implementation reproduces the published flow-density diagram
quantitatively: peak flow ≈ 0.46 at ρ ∈ [0.10, 0.12], free-flow mean
velocity ⟨v⟩ = 4.70 in the dilute limit (matching the analytic
prediction v_max − p = 4.7), and the characteristic asymmetric shape
in which flow rises linearly from ρ = 0 to the peak and decays roughly
linearly back to zero at ρ = 1. At p = 0 the transition sharpens to a
discontinuity at ρ_c = 1/(v_max + 1) = 1/6 ≈ 0.167, consistent with
the analytic result. At p = 0.5 the peak shifts to lower density
(ρ ≈ 0.08) and lower magnitude (flow ≈ 0.32), also matching the
published trend. All three regimes (p = 0, p = 0.3, p = 0.5)
replicate within < 3% of the published values at three seeds.

**P8 detector — primary metric.** P8's primary metric is the stopped-
car fraction P(v = 0) = ⟨1[v_i(t) = 0]⟩, averaged over post-burn-in
timesteps and all cars. This is the Bette-Habel-Emig-Schreckenberg
(2017) order parameter for the NS jamming transition; the BHES paper
also decomposes the jammed-car fraction P(v ∈ {0, 1}) into three
physically meaningful factors (jamming rate × jam lifetime × jam size)
and derives random-walk scaling exponents near the critical density.
Our characterization at L = 1000, 1000 burn-in + 2000 measurement,
three seeds, gives ρ = 0.05 → stopped = 0.000, ρ = 0.12 → stopped =
0.082 (onset), ρ = 0.15 → stopped = 0.181, ρ = 0.30 → stopped = 0.431
(deep jam), with seed-to-seed standard deviation < 0.005 at every
density. The stopped-fraction order parameter has tiny finite-size
fluctuation and a sharp transition around ρ ≈ 0.10–0.12 — the same
region where flow peaks.

**The density-saturation false-positive.** Before locking P8's
thresholds, we ran the detector across every NS parameter regime
including deterministic p = 0 at saturation density ρ = 0.80. Under
pigeonhole constraint, a ring holding 0.80 L cars-worth of traffic
cannot have every car moving: at minimum, a fraction (v_max ρ − 1)/
v_max of cars must be standing still at each timestep simply because
there isn't room for every car to advance. At ρ = 0.80, p = 0 this
gives stopped-fraction = 0.750 — high above the screening threshold
(0.05) and the definitive threshold (0.15). A P8 variant gated only on
stopped-fraction would false-positive here, reporting DEFINITIVE jamming
in a regime with no stochastic dynamics and no emergent stop-go waves —
only saturation geometry. The sharp analogue of Sprint 13's RPS-vs-
GrayScott false-positive trap.

**The jam-lifetime p95 discriminator.** The resolution was to introduce
a secondary metric — the 95th-percentile of per-car consecutive-v = 0
run lengths — and make it a hard gate at the CONFIRMATION tier.
Emergent NS jamming at canonical parameters produces heavy-tailed
lifetime distributions: p95 = 13 at ρ = 0.15, with a maximum lifetime
exceeding 50 timesteps. Pigeonhole density-saturation at ρ = 0.80,
p = 0 produces uniform short stops — every car is blocked by its front
neighbor for at most a few timesteps before geometry releases it —
giving p95 = 4 and a maximum lifetime of 6 with no heavy tail. The
three-fold separation between the two regimes (p95 = 13 vs p95 = 4)
motivates the CONFIRMATION threshold of 5; at no parameter choice does
genuine NS jamming give p95 ≤ 5, and at no deterministic saturation
choice does density-only stopping give p95 > 5. This is the load-
bearing Sprint 15 finding: a confirmation-tier gate cannot be the
primary's effect-size alone when the primary has a trivial geometric
inflation mode; a secondary that measures the distribution's tail
structure is the correct discriminator.

**Substrate-level discrimination.** P8 is gated on two content-level
prerequisites in addition to lattice_1d substrate registration:
presence of a 1D integer `velocities` observable, and values in
[0, 64]. These gates cleanly reject all fourteen other catalog models
without empirical thresholding. Zhang cell-view sorting shares the
lattice_1d substrate but exposes `array` and `cell_types` rather than
`velocities`; it rejects at observable-prereq. Vicsek and D'Orsogna
expose `velocities` but as continuous 2D vectors; they reject at the
non-1D / non-integer check. All lattice_2d and lattice_2d_continuous
models reject at substrate_mismatch. This matches the Sprint 13
Decision 37 architectural pattern: discrimination by substrate content
beats discrimination by empirical threshold.

**Null model and tier gating.** P8's null permutes each car's velocity
time-series independently, preserving per-car marginal v-distribution
(so stopped-fraction is invariant under the null) while destroying the
temporal persistence that produces heavy jam-lifetime tails. Under the
null, p95 collapses to the geometric-distribution p95 at the observed
stopped marginal, typically 2 timesteps at canonical parameters. With
199 permutations the null-p floor is 1/200 = 0.005; Cohen's d is
effectively infinite because the null distribution is concentrated
(null-std ≈ 0 across 199 draws from nearly identical shuffles). For
cost control the null subsamples up to 80 cars per permutation; the
per-car locality of the shuffle makes p95 extensively insensitive to
car count above ~50. Tier thresholds: stopped-fraction > 0.05 for
SCREENING; jam-lifetime-p95 > 5 with null-p < 0.01 for CONFIRMATION;
stopped-fraction > 0.15 with jam-lifetime-max > 20 for DEFINITIVE.

**Canonical-positive outcomes.** At L = 1000, v_max = 5, p = 0.3,
ρ = 0.15, 1000 burn-in + 2000 measurement, seeds {42, 123, 2024}:

| Seed | stopped | lt_p95 | lt_max | null-p | Cohen's d | Tier       |
|------|---------|--------|--------|--------|-----------|------------|
| 42   | 0.178   | 13     | 57     | 0.005  | > 1e6     | DEFINITIVE |
| 123  | 0.182   | 13     | 63     | 0.005  | > 1e6     | DEFINITIVE |
| 2024 | 0.185   | 14     | 68     | 0.005  | > 1e6     | DEFINITIVE |

**Transfer-matrix additions.** Nagel-Schreckenberg × P8 is the
Sprint 15 positive (DEFINITIVE). Every other model × P8 rejects: Zhang
at observable-prereq (no `velocities`); all lattice_2d, continuous_2d,
oscillator, opinion_space, and lattice_2d_continuous models at
substrate_mismatch (nine cells). Every 2D-substrate detector × NS
rejects at its own observable prereq: P1 needs `grid`, P3 needs
`field`, P11/P12/P13/P15/P22 need `grid`. Seventeen new cells (one D +
sixteen rej) join the audited transfer matrix, bringing the total to
80 audited cells across 15 models × 14 detectors.

## 4.15 Consolidated Transfer Matrix

The following matrix summarizes detection outcomes across all 80
audited model × detector pairs. Entries show the highest achieved tier
(D = definitive, C = confirmation, S = screening), reported non-
detections (rej = rejected by prerequisite or screening guard, nd =
substrate-compatible but not detected), or substrate/observable
incompatibility (× = substrate mismatch, — = detector requires
different substrate). D* denotes Game of Life with dense random IC (the
canonical configuration for the generalized P15 detector; R-pentomino
gives a lower P15 tier because it lacks the diversity of outcome
classes the generalized detector requires).

|                   | P1  | P3  | P5 | P6 | P8  | P9 | P11 | P12 | P13 | P14 | P15 | P21 | P22 | P27 | P31 |
|-------------------|-----|-----|----|----|-----|----|-----|-----|-----|-----|-----|-----|-----|-----|-----|
| Zhang sorting     | S   | ×   | ×  | ×  | rej | ×  | ×   | ×   | ×   | ×   | ×   | ×   | ×   | ×   | C   |
| Schelling         | C   | rej | ×  | ×  | ×   | ×  | rej | —   | rej | ×   | nd  | ×   | rej | ×   | —   |
| Vicsek (ordered)  | ×   | ×   | D  | rej| ×   | ×  | ×   | ×   | ×   | ×   | ×   | ×   | ×   | ×   | ×   |
| D'Orsogna (mill)  | ×   | ×   | rej| D  | ×   | ×  | ×   | ×   | ×   | ×   | ×   | ×   | ×   | ×   | ×   |
| Kuramoto (sync)   | ×   | ×   | ×  | ×  | ×   | D  | ×   | ×   | ×   | ×   | ×   | ×   | ×   | ×   | ×   |
| GH spiral         | S   | ×   | ×  | ×  | ×   | ×  | ×   | rej | C   | ×   | rej | ×   | rej | ×   | —   |
| GoL (R-pent/rand) | rej | ×   | ×  | ×  | ×   | ×  | ×   | rej | rej | ×   | D*  | ×   | rej | ×   | —   |
| BTW sandpile      | ×   | ×   | ×  | ×  | ×   | ×  | ×   | ×   | ×   | D   | nd  | ×   | ×   | ×   | ×   |
| Nowak-May (b=1.8) | C   | rej | ×  | ×  | ×   | ×  | rej | rej | rej | ×   | S   | ×   | rej | D   | ×   |
| HK (ε=0.2)        | ×   | ×   | ×  | ×  | ×   | ×  | ×   | ×   | ×   | ×   | ×   | D   | ×   | ×   | ×   |
| SIR epidemic      | rej | rej | ×  | ×  | ×   | ×  | rej | rej | rej | ×   | nd  | ×   | D   | ×   | —   |
| RPS spatial       | S   | rej | ×  | ×  | ×   | ×  | rej | C   | rej | ×   | nd  | ×   | S   | ×   | —   |
| Lotka-Volterra    | C   | rej | ×  | ×  | ×   | ×  | D   | rej | rej | ×   | nd  | ×   | S   | ×   | —   |
| Gray-Scott        | rej | D   | ×  | ×  | ×   | ×  | rej | rej | rej | ×   | nd  | ×   | rej | ×   | ×   |
| Nagel-Schreck.    | rej | rej | ×  | ×  | D   | ×  | rej | rej | rej | ×   | rej | ×   | rej | ×   | ×   |

Row count: 15 distinct model families. Column count: 15 displayed detector
slots (14 formally registered in the orchestration layer plus P11 which
was implemented in Sprint 11 but not retrofitted into the detector
registry — a pre-existing documentation-vs-registry gap carried forward
from Sprints 11–14). Out of 225 displayed cells, 147 are substrate
mismatches (×) and 7 are detector-substrate incompatibilities (—,
primarily P31 which requires lattice_1d). The remaining 71 cells are
substrate-compatible and empirically audited: 23 produce detections
(10 at DEFINITIVE — adding Nagel-Schreckenberg × P8 to the Sprint 13
tally of 9, 1 at DEFINITIVE with dense random IC (D*, Game of Life ×
P15), 6 at CONFIRMATION, 6 at SCREENING), 42 are rejected by
prerequisite or screening guard (rej; Sprint 15 added sixteen such
cells — nine existing-model × P8 rejections and seven NS × 2D-
substrate-detector rejections), and 6 run but do not fire (nd —
typically P15 on stochastic lattice models where the functional replay
test fails due to irreproducibility, or Gray-Scott's deterministic PDE
which has no stochastic step_fn). Of these 71 displayed audited cells,
68 appear in the cross-detection-matrix `EXPECTED_OUTCOMES` regression
table (`tests/test_cross_detection_matrix.py`); the remaining canonical
positives (e.g., Vicsek × P5, Kuramoto × P9, Gray-Scott × P3,
Nagel-Schreckenberg × P8) have their DEFINITIVE tier pinned by
dedicated e2e test files rather than the cross-detection matrix. Nine
additional audited cells in the LV × P11 row and P11 column (the
Sprint 11 work) live outside the orchestration registry but inside the
test suite.

Eight observations about the matrix.

*No false positives across substrate boundaries.* The substrate-aware
dispatch system (seven substrate types: lattice_1d, lattice_2d,
lattice_2d_continuous, continuous_2d, oscillator, opinion_space,
predator_prey) correctly prevents cross-substrate detector application.
Substrate mismatches account for the bulk of the blank cells and never
fire a detector. The Sprint 13 addition of lattice_2d_continuous for
Gray-Scott and the Sprint 14 hardening of P1's graceful-reject path
for missing-grid substrates extended this guarantee to the
continuous-field case.

*Clean within-cluster discrimination.* P5 / P6 show perfect
cross-exclusion (D'Orsogna milling rejected by P5 at φ = 0.046;
Vicsek flocking rejected by P6 at |L| = 0.031). P11 / P12 show the
bilateral-versus-cyclic boundary: LV triggers P11 cleanly and is
rejected by P12 (intransitivity_score = 0.24 << 1.0, because the
dominance graph is not cyclic); RPS triggers P12 at CONFIRMATION and
is rejected by P11 (n_unique_species_observed = 3 ≠ 2). P13 / P15
separate via the two-stage TE + functional discriminator. P3 /
aggregation detectors separate via substrate-level content gates:
continuous fields vs integer-grid observables, and
n_unique_values ≥ 50 as the secondary prerequisite.

*Co-occurrence rather than mutual exclusion.* Three models fire
multiple detectors: Nowak-May (P1 confirmation + P27 definitive + P15
screening); LV (P1 confirmation + P11 definitive); Schelling (P1
confirmation alone, with GH-style rejections on P13 and P22). In each
case the co-occurring patterns have compatible mechanisms (Nowak-May's
cooperator clusters exhibit both aggregation signatures and
reciprocity; LV's spatial domains exhibit both aggregation and
predator-prey oscillation). The exclusion graph captures structural
incompatibilities, not the combinatorial space of co-occurrence.

*Observable-level filtering beyond substrate.* P27 requires
`coop_fraction`, preventing false matches with GH, GoL, Schelling,
SIR, RPS, and LV despite shared lattice_2d substrate (extended to
include the new lattice_2d-like predator_prey class). P14 requires
`avalanche_sizes`, restricting to BTW. P11 requires `prey_count` and
`predator_count` (or an explicit species hint), restricting meaningful
detection to two-species systems with nontrivial empty reservoir. P3
requires a `field` observable and n_unique_values ≥ 50, restricting
detection to genuinely continuous spatial fields rather than
integer grids re-interpreted as floats (Sprint 13, Decision 37).

*The SIR × P1 rejection is informative.* Section 4.10 describes the
characterization that flipped this cell from `S` (under a peak-based
primary metric) to `rej` (under the final-state primary). The flip
preserves detection on Schelling, Nowak-May, and RPS while correctly
rejecting SIR — and this non-obvious combination was only achievable
because the characterization tested the alternative hypothesis
(peak → sustained) on actual data rather than adopting it on prior
grounds. The lesson generalizes: a detector's primary metric is a
hypothesis about what signal distinguishes the target pattern, and
hypotheses need empirical testing.

*Guard-based rejections are informative.* Six rejections required
detector guards beyond simple threshold checks: GoL × P1 (type
constancy), GoL × P13 (excitable medium guard, n_states < 3),
Nowak-May × P11 (total_std prerequisite, conservation trap), RPS ×
P11 (n_species prerequisite, bilateral-vs-cyclic separation),
SIR × P1 (final-Moran primary, transient-vs-sustained separation),
and Gray-Scott × P1 (substrate-warning graceful-reject path added in
Sprint 14 B.1). Each guard sharpened the operational definition of
the target pattern.

*The RPS × P3 rejection required substrate-level discrimination.* The
Sprint 13 characterization found that RPS at low mobility produces
raw-grid radial-FFT peak-to-mean ≈ 23.10, numerically exceeding the
Gray-Scott labyrinth value of 18.75. No empirical threshold on
peak-to-mean can separate these two systems. Discrimination therefore
operates at the content level: P3 requires a `field` observable (RPS
exposes only `grid`) and n_unique_values ≥ 50 (RPS grid values are
limited to the number of species, typically 3–4). When RPS is
adversarially re-labeled as a `field`, it rejects at the n_unique
prereq regardless of its p/m = 23.10. This is the architecturally
cleaner pattern: where two canonical models produce overlapping signal
on the intended primary metric, look for a content-level property
that separates them before reaching for fragile thresholds.

*Tier ceilings are meaningful.* Several detections are capped at
CONFIRMATION rather than DEFINITIVE (Schelling × P1, GH × P13, Nowak-May
× P1, LV × P1, RPS × P12). In each case, the ceiling arises from
specific methodological constraints — permutation-count floors (999
for P1, hitting p = 0.001 but not p < 0.001), absence of mechanistic
null (GH), intrinsic signal asymmetry (Nowak-May's imitation-based
clustering is weaker than Schelling's preference-based clustering on
segregation index) — rather than weak signals.

## 4.16 Methodological Lessons

Six methodological insights emerged from the replication work across
Sprints 1–11.

**Look before touching.** The most important pattern from Sprints 9,
10, and 11 is that empirical characterization on the target and
adjacent models must precede detector design. Sprint 10's SIR × P1
resolution emerged from a six-model characterization that would have
been wasted as a post-hoc verification; Sprint 11's P11 detector
reshaped its primary metric (positive-lag → anti-phase) and acquired a
critical prerequisite (conservation check, from Nowak-May's
ρ_anti = −0.98 trap) only because both were observed before the
detector code was written. A detector that is designed from a prior
and then validated on data may pass its own tests and still be wrong.

**Statistical power.** Underpowered tests produced incorrect results
in multiple cases: P31 non-redundancy with 48 runs gave false negatives
that were only resolved at 600 runs; P5 on Vicsek required ≥5,000
recorded steps for definitive detection; KSG TE at T = 2,000 produced
a false positive (p = 0.03) on independent data, requiring T ≥ 5,000
for correct non-significance. We now enforce minimum power
requirements for all permutation-based tests (≥99 permutations for p
= 0.01 floor, ≥199 for p = 0.005, ≥999 for p = 0.001) and minimum
run lengths tied to the detector's sensitivity (P11 requires ≥1,200
generations at L = 100 for DEFINITIVE; shorter runs give CONFIRMATION).

**Test correctness.** In multiple cases, initial test failures were
caused by errors in the test, not the code: P15 fidelity used overly
fine-grained outcome signatures; P15 collision test placed gliders too
far apart to collide within the measurement window; P1 guard sampled
from step 7 instead of step 0, missing fast Schelling convergence; the
initial LV parameter choice (λ = 2) was in the extinction regime. The
lesson: when a test fails, the first response should be to verify the
test is correct, not to weaken the assertion.

**Boundary-conditioned measurement.** The P13/P15 TE discriminator
required restricting measurement to boundary cells rather than
averaging over the full lattice. Raw average TE gives the wrong
ordering (GH > GoL) because deterministic wave propagation creates
trivially high TE at interior cells. Boundary conditioning isolates the
interaction signal. Section 6 develops the general principle.

**Intrinsic timescales.** D'Orsogna milling detection failed when
using the initial particle spread to estimate the crossing time T_cross
— the mill compacts significantly during formation (from radius 5.0 to
diameter 3.02). Using the measured group diameter after equilibration
gave the correct timescale. Similarly, Kuramoto entrainment must
measure frequency dispersion among locked oscillators only; Lorentzian
tail outliers that never lock dominate the all-oscillator statistic.
LV exhibits the same pattern: the initial 100-generation transient is
deterministic-like and must be excluded (burn-in) before the
quasi-stationary oscillation signal can be measured.

**Broad negative-model sweeps before locking thresholds.** Sprint 11's
most valuable routine check was running the proposed P11 primary
metric against every existing two-species lattice model before
finalizing the detector. The ρ_anti = −0.98 result on Nowak-May
(caused by strict conservation, not predator-prey coupling) would not
have been caught by testing only the planned negatives (Schelling,
SIR). A detector passing its intended positive and the specific
negatives the designer anticipated is insufficient: the adversary is
not only the absent pattern but also the confounders one hasn't
thought of. The catalog's internal model inventory makes this sweep
cheap; we now treat it as a standard step between detector design and
detector acceptance.
