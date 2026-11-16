# P33–P37 on the full Phase-2a panel (v12)

**Date:** 2026-11-15. **Branch:** `panels-p33-p37` (from `validation-rebuild` @ b00271e).
**Code commit the panels ran on:** 8802759 (recorded as `head_commit` in every v12 JSON).
**Outputs:** `analysis/outputs/p3{3..7}_*_phase2a_panel_v12.json` (new files; the
discovery-phase JSONs `p3{3..7}_*_phase2a_panel.json` are untouched).

The five discovery-phase detectors (P33 active nematic, P34 adaptive-network
fission, P35 entropy crystallization, P36 Pareto exchange, P37 resource
competition) had only standalone validations: 3–6 hand-picked negative types × 5
seeds, no Class A, and for P37 no Class A or B. This pass runs each one through
the same harness the other 31 detectors went through (`epc/phase2a/panel.py::
run_panel`, spec v1.2): 5 canonical positives, Class A structureless surrogates,
Class B catalog mates / look-alikes, Class C failed regimes of the target's own
model. A negative counts as rejected only if the detector returns *not detected*
at every tier.

**The detectors are frozen.** No file under `epc/detectors/`, no threshold, gate
or factory default was changed (`git diff b00271e -- epc/detectors` is empty).
Each panel calls the battery factory `make_p3X_detector_fn()` with its defaults.
Everything below is what those frozen detectors did.

## Results at a glance

| | P33 active nematic | P34 adaptive network | P35 entropy crystal | P36 Pareto exchange | P37 resource competition |
|---|---|---|---|---|---|
| detector format | orientation_field (new) | adjacency_timeseries (new) | positions (new) | scalar_wealth (existing, P28's) | abundance_resource_timeseries (new) |
| positives detected | 5/5 definitive | 5/5 definitive | 5/5 definitive | 5/5 definitive | 5/5 definitive |
| Class A TNR | 9/9 = 1.000 (1 skipped) | 10/10 = 1.000 | 9/9 = 1.000 (1 skipped) | 8/8 = 1.000 (2 skipped) | **1/8 = 0.125** (2 skipped) |
| Class B TNR | **6/8 = 0.750** | 5/5 = 1.000 | 9/9 = 1.000 | 5/5 = 1.000 | **1/4 = 0.250** (advisory, n<5) |
| Class C TNR | 8/8 = 1.000 | 8/8 = 1.000 | **7/8 = 0.875** | 10/10 = 1.000 | 10/10 = 1.000 |
| overall TNR | 23/25 = **0.920** | 23/23 = 1.000 | 25/26 = **0.962** | 23/23 = 1.000 | 12/22 = **0.545** |
| continuous d (harness) | 7.86 | 26.99 | 6.66 | 22.37 | **−0.09** |
| discrete d (tier confidence) | 7.98 | NaN (perfect separation) | 8.03 | NaN | 1.36 |
| harness verdict | PARTIAL | TNR-PASS-EFFECT-UNDEFINED | PASS-with-weakness | TNR-PASS-EFFECT-UNDEFINED | PARTIAL |
| false positives | 2 (Class B, screening) | 0 | 1 (Class C, confirmation) | 0 | 10 (7 A definitive, 3 B) |
| standalone validation said | TNR 20/20, d 4.81 | TNR 25/25, d 16.74 | TNR 30/30, d 8.31 | TNR 25/25, d 13.87 | TNR 15/15, d 15.09 |

"continuous d" is the harness's `cohens_d_continuous`: Cohen's d of the detector's
own continuous discriminating metric (the `CONTINUOUS_METRIC` entry, already
registered for P33–P37), 5 positives vs every evaluated negative, sign-oriented so
positive = correct separation, sample-variance pooling. It is not the same
statistic as the standalone scripts' d (population variances, different negative
set), so the two rows are not directly comparable. "TNR-PASS-EFFECT-UNDEFINED" is
the harness label for TNR ≥ 0.95 with a constant-vs-constant discrete score (the
discrete d is NaN by design); the continuous d is the effect size to read.

Headline: P34 and P36 hold up fully on the larger panel. P33 and P35 each let
through a small number of negatives (P33: two catalog mates read through its
particle-binning path; P35: one sub-freezing compression seed). P37 does not hold
up: its detector fires on most structureless temporal noise, on a monotone drift,
on a period-2 alternation, and on oscillating look-alikes that have no
field-mediated competition, and its continuous metric no longer separates
positives from the pooled negatives (d ≈ −0.1).

## How the panels were built

### Positives (all five)
Exactly the standalone validations' positives, same model call and parameters,
seeds 0–4, each starting from a disordered state (`build_p3X_positives` in
`analysis/run_phase2a_panel.py`):

| | model call | start |
|---|---|---|
| P33 | `active_nematic_field(seed, eta=0.18)` (G=64, 400 steps, 60 frames) | uniformly random director |
| P34 | `adaptive_voter(seed, N=200, mean_deg=4, phi=0.99, max_steps=400000, n_frames=80)` | connected ER graph, random opinions |
| P35 | `hard_disk_crystallization(seed, N=196, L=14, eta_end=0.74, n_steps=2500, n_frames=80, noise0=0.16)` | uniform random positions |
| P36 | `kinetic_exchange(seed)` (N=3000, λ_i ~ U[0,1), 5000 sweeps, 40 frames) | equal wealth |
| P37 | `resource_competition(seed)` (Huisman–Weissing 5×5, t=3000, 80 frames) | low seed-varied abundances |

The harness passes the seed-0 metadata to every substrate (standard). For P35 see
decision 3.

### Class A — formats and what the surrogates are
Four detectors read inputs no existing format produced, so each got a new
`detector_format` branch (`epc/phase2a/panel.py`) and generators in
`epc/phase2a/synthetic_discovery.py`, dispatched from the ten generators in
`synthetic.py`. Each surrogate is written in exactly the frame schema of the
pattern's positive. P36 reads conserved scalar wealth, which is the existing
`scalar_wealth` format (P28's); its Class A branches are reused unchanged.

| generator | orientation_field (P33; 64×64 director lattice) | adjacency_timeseries (P34; N=200) | positions (P35; N=196, L=14) | abundance_resource_timeseries (P37; 5 species, 5 resources) |
|---|---|---|---|---|
| random_uniform | director i.i.d. U[0,π) per site, re-drawn every frame; heading = director + random head/tail flip | fresh ER graph (mean degree 4) + fresh random opinions every frame | ideal gas re-drawn every frame | abundances i.i.d. U(0,1), resources i.i.d. U(0,1)·S |
| random_gaussian | director i.i.d. (π/2)·N(0,1) mod π (std = half the range, the phases-format analog) | fresh Chung–Lu graph, expected degrees \|N(4,2)\| | Gaussian blob (std L/4) re-drawn every frame | abundances N(1,0.3) clipped at 0 |
| random_binary | director Bernoulli(0.5) on {0, π/2} | fresh dense Bernoulli(0.5) graph | random half-occupation of a 20×20 square lattice, re-drawn | abundances ∈ {0.5, 1.5} |
| spatial_white_noise | as random_uniform (different draw) | as random_uniform | as random_uniform | abundances i.i.d. Exp(1) |
| temporal_white_noise | each site's director an independent random walk (step N(0,0.3)) | edge-flip random walk at stationary density + independent opinion flips (q=0.05/frame) | independent Brownian walkers (step 0.1) | independent log-random walks (step 0.05) |
| permutation_shuffled | positive's final frame, sites permuted (director, heading together), held static | positive's final graph with its edges re-placed at random (same edge count), held static — trail_network convention | consistent particle relabelling (degenerate; skipped) | consistent species relabelling (degenerate; skipped) |
| time_shuffled | frame order of the positive shuffled (generic harness branch, all formats) | same | same | same |
| constant | all directors 0, heading = director | complete graph, one opinion | one random configuration frozen | all abundances 1, resources = supply |
| linear_gradient | director rotating by π across the box (defect-free texture), static | static ring lattice (degree 4), opinion step along the ring | frozen configuration with density rising linearly in x | abundances rising linearly in time, resources falling linearly (monotone drift, the timeseries-format convention) |
| checkerboard | directors 0 / π/2 on a checkerboard, static | static bipartite ring (offsets 1, 3) with alternating opinions | frozen perfect 14×14 square lattice | period-2 alternation 1 ± 0.5, neighbouring species in antiphase |

All are structureless: i.i.d. draws, uncoupled per-unit walks, or fixed
deterministic textures with no dynamics. None contains any catalog mechanism.

### Invariance flags (decided before the panels ran)
Set in `epc/phase2a/detector_invariance.py` from each detector's metric
definition, with the rationale stored next to the flag. They were not flipped
after seeing results. Every skipped substrate was still scored, outside the TNR,
by `analysis/validation_rebuild/p33_p37_panel_diagnostics.py`:

| | permutation_invariant | time_shuffle_invariant | skipped substrate scored out-of-TNR (5 seeds) |
|---|---|---|---|
| P33 | False (site adjacency matters) | True (detection metric = per-frame scores averaged over the second half) | time_shuffled: fires 5/5 (4 screening, 1 definitive) |
| P34 | False (edge permutation changes Q) | False (Q_early vs Q_late is order) | none skipped |
| P35 | True (psi6 invariant to particle labels) | False (psi6 early vs late is order) | permutation_shuffled: fires 5/5 definitive (it is the positive, relabelled) |
| P36 | True (sorted-wealth statistics) | True (late-quarter set statistics, no growth check) | both fire 5/5 definitive |
| P37 | True (per-species CV median) | True (late-third set statistics) | both fire 5/5 definitive |

For P33 the time_shuffled result is informative beyond degeneracy: P33's
detection tier has no temporal-order requirement, so a frame-shuffled active
nematic is still "detected"; only its confirmation tier (emergence from the first
frame) reads order. For P36 and P37 the time_shuffled firing reflects that neither
detector checks how the state arose.

### Class B — typing
P33–P37 are not in `MODEL_REGISTRY`, so `catalog.DISCOVERY_CLASS_B` declares each
pattern's substrate type, any extra cross-type catalog mates, and the Class B'
look-alikes (`epc/phase2a/discovery_lookalikes.py`). The patterns are deliberately
not added to `PATTERN_TO_SUBSTRATE_ID`, because that would make them mates in the
already-committed panels of P2/P5/P6/P7/P17/P19/P28. Checked: `class_b_for_pattern`
returns identical output for every pre-existing pattern before and after this change.

New catalog plumbing: a generator for `P28_yard_sale` (the P28 canonical-positive
parameters, seed 0, which is not a P28 positive seed) and five format adapters.
The new adapters never tile a short trajectory (tiling would replay the transient
mid-run). Longer natives are subsampled uniformly, shorter ones kept at their own
length, and each frame carries the substrate's own box.

### Class C — ground truth, not detector verdicts
Each Class C regime is the target's own model at a setting where the pattern
does not form, or with its mechanism disabled. "Does not form" was established
from each model's own observables by
`analysis/validation_rebuild/p33_p37_class_c_groundtruth.py` (output
`analysis/outputs/p33_p37_class_c_groundtruth.json`), before any panel ran. Two
seeds per regime, seeds 200–209, distinct from the positive seeds 0–4.

## Per pattern

### P33 — active nematic (+-1/2 defect gas)
**Panel: 5 positives; A 9 (+1 skipped); B 8; C 8; 25 negatives.**

| class | negative | what it is / why it belongs | result |
|---|---|---|---|
| A | 9 surrogates above (time_shuffled skipped) | structureless director fields | all none |
| B mate | P2_abp | active Brownian MIPS (200 particles, box 16) | none |
| B mate | P5_vicsek | polar flock (the P5 look-alike, real model) | none |
| B mate | P6_dorsogna | milling (+1 vortex, real model); translated into [0, box) | none |
| B mate | **P7_lane_formation** | counterflow lanes, 20×4 corridor, box 20 | **screening (FP)** |
| B mate | **P17_collective_sensing** | gradient-climbing group, 50 agents, box 20 | **screening (FP)** |
| B mate | P19_informed_minority | led flock | none |
| B' | p33_polar_flock_field | standalone neg_polar_flock (lattice) | none |
| B' | p33_milling_vortex_field | standalone neg_milling (lattice) | none |
| C | uniform_nematic η=0.05 aligned ×2 | standalone neg_uniform_nematic: ordered, zero defects at 400 and 4000 steps | none, none |
| C | ordered_nematic η=0.10 aligned ×2 | S_loc 0.96, zero defects at 400 and 4000 steps | none, none |
| C | isotropic η=0.60 ×2 | disordered phase, S_loc 0.19 (isotropic baseline 0.15) | none, none |
| C | isotropic η=0.80 ×2 | S_loc 0.16 | none, none |

The standalone `neg_isotropic` (i.i.d. random orientations) is the same object as
Class A random_uniform / spatial_white_noise, so it is covered there rather than
duplicated.

**False positives (2/25), both Class B catalog mates, both screening only**
(secondary metrics re-derived in `p33_p37_panel_diagnostics.json`):
- *P17_collective_sensing*: D* = 0.0066 (screening gate 0.004). 50 agents binned
  onto the detector's 48×48 director grid leave most cells empty; the detector's
  own neighbour-fill smooths the few headings into large uniform patches (S_loc
  0.90), polar order is low (φ 0.09), |L| 0.06, and the patch seams register as
  ±1/2 windings. Not emerged (S_loc 0.83 from the first frame), not sustained (min
  half-defect density 0.0035), so it stops at screening.
- *P7_lane_formation*: D* = 0.0115. Counterflow lanes are a genuinely *apolar*
  ordered velocity field (two opposed streams, φ = 0.005), which is exactly what
  P33's polar gate cannot exclude. The 20×4 corridor also occupies only a strip of
  the detector's square grid (it reads a single box_size). Sustained, but not
  emerged (S_loc 0.52 → 0.48), so screening.

Both fire through the detector's particle-binning input path. The standalone
validation exercised that path only on dense Vicsek runs (400 particles in a
10-box, `--harden`), never on sparse or non-square particle systems. Lattice
director fields, the positive's own representation, were all rejected.
Continuous d = 7.86.

### P34 — adaptive-network fission
**Panel: 5 positives; A 10; B 5; C 8; 23 negatives. No false positives.**

| class | negative | what it is / why it belongs | result |
|---|---|---|---|
| A | 10 surrogates above | structureless / static graphs | all none |
| B mate | P18_voter | voter on its fixed periodic Moore lattice graph (1024 nodes), cell states as opinions; 40 frames | none |
| B mate | P21_hegselmann_krause | HK bounded-confidence *interaction* graph (|x_i − x_j| < ε), re-derived each frame, so it itself evolves and fragments as opinions cluster (9 frames: HK stops at convergence) | none |
| B' | p34_static_modular | standalone: fixed 2-block SBM (high Q, no emergence) | none |
| B' | p34_pre_fragmented | standalone: two disconnected halves from t=0 | none |
| B' | p34_ba_growth | Barabási–Albert growth (m=2): the canonical evolving network, no community mechanism | none |
| C | fixed_network_voter ×2 | standalone: no rewiring (coevolution off) | none, none |
| C | opinion_blind_rewiring ×2 | standalone null_random_rewire: rewiring ignores opinions | none, none |
| C | subcritical φ=0.50 ×2 | same model; ground truth: consensus on one connected component | none, none |
| C | subcritical φ=0.80 ×2 | same; ground truth: no fission for φ ≤ 0.80, fission at φ ≥ 0.90 | none, none |

Continuous d = 26.99. Every negative has modularity_gain 0, the gated metric's
value unless Q_early ≤ 0.22 and Q_late ≥ 0.30 both hold.

### P35 — entropy-driven crystallization
**Panel: 5 positives; A 9 (+1 skipped); B 9; C 8; 26 negatives.**

| class | negative | what it is / why it belongs | result |
|---|---|---|---|
| A | 9 surrogates above (permutation_shuffled skipped) | structureless point sets | all none |
| B mates | P2_abp, P5_vicsek, P6_dorsogna, P7_lane_formation, P17_collective_sensing, P19_informed_minority | continuous_2d particle systems, each in its own box | all none |
| B' | p35_static_crystal | standalone: final crystal held fixed (high ψ6, zero gain) | none |
| B' | p35_keller_segel | standalone (T2c system): chemotactic clumping, L=50 | none |
| B' | p35_agent_active_nematic | standalone (T2c system): self-propelled rods, nematic alignment, L=18 | none |
| C | dilute_gas η_end=0.20 ×2 | standalone neg_dilute_gas (same call) | none, none |
| C | fluid η_end=0.40 ×2, 0.55 ×2 | compressed only into the liquid range | none ×4 |
| C | fluid η_end=0.65 ×2 | below the hard-disk freezing range (liquid to η ≈ 0.700; Bernard & Krauth 2011) | none (seed 206), **confirmation (seed 207, FP)** |

**False positive (1/26):** fluid η_end=0.65, seed 207: ψ6 0.363 → 0.436, clearing
the hexatic floor (0.42) by 0.016 with gain 0.073 (confirmation gate 0.07). For
comparison the positives' ψ6 gains are 0.077–0.126. This regime is a judgment call
(decision 6). By the literature it is below freezing. In this minimal soft-disk
quench, though, the structural gap between η=0.65 and the positive's η=0.74 is
small. The model's own crystalline-particle fraction (local |ψ6| > 0.7) in the
final frame is 0.02–0.05 at η_end 0.20–0.55, 0.13 at 0.65, 0.15–0.16 at 0.70 and
0.22–0.23 at 0.74. Without the η=0.65 regime, P35 would be 24/24.
Continuous d = 6.66.

### P36 — Pareto wealth tail from heterogeneous exchange
**Panel: 5 positives; A 8 (+2 skipped); B 5; C 10; 23 negatives. No false positives.**

| class | negative | what it is / why it belongs | result |
|---|---|---|---|
| A | 8 scalar_wealth nulls (P28's: uniform, gaussian, binary, exponential ×2, constant, static gradient, binary) | static / resampled distributions | all none |
| B mate | P28_yard_sale | the catalog's P28 canonical positive (condensation) | none |
| B' | p36_yard_sale_condensation | standalone neg_condensation | none |
| B' | p36_lognormal_static | standalone neg_lognormal | none |
| B' | p36_bouchaud_mezard | multiplicative growth + weak redistribution (non-conserved route to power-law tails); parameters copied unchanged from the P28 Class C module | none |
| B' | p36_yard_sale_redistribution | yard-sale with wealth tax χ=0.05; P28 Class C parameters unchanged | none |
| C | homogeneous λ = 0.00 ×2 | standalone neg_boltzmann (exponential) | none, none |
| C | homogeneous λ = 0.50 ×2 | standalone neg_uniform_saving (Gamma) | none, none |
| C | homogeneous λ = 0.20 ×2, 0.80 ×2 | same model, other single λ (Gamma; Patriarca et al. 2004) | none ×4 |
| C | bounded heterogeneity λ_i ~ U[0,0.5) ×2 | heterogeneous, but the power-law tail needs λ density near 1 (Chatterjee et al. 2004; Mohanty 2006) | none, none |

Continuous d = 22.37. One caveat comes from the skipped-substrate diagnostic.
P36 has no emergence or dynamics check, so the frame-shuffled positive is
detected 5/5, and nothing in the detector asks how a tail arose. A static
heavy-tailed draw is not in this panel, because P28's scalar_wealth Class A has
none, so that case is untested. P36's clean result means it rejects every wealth
mechanism tested here. It does not show that the detector checks the
mechanism.

### P37 — field-mediated resource competition (Huisman–Weissing)
**Panel: 5 positives; A 8 (+2 skipped); B 4 (advisory); C 10; 22 negatives.**

| class | negative | what it is / why it belongs | result |
|---|---|---|---|
| A | **random_uniform** | i.i.d. U(0,1) abundances | **definitive (FP)**, osc_cv 0.60 |
| A | **random_gaussian** | i.i.d. N(1,0.3) | **definitive (FP)**, 0.30 |
| A | **random_binary** | i.i.d. {0.5,1.5} | **definitive (FP)**, 0.50 |
| A | **spatial_white_noise** | i.i.d. Exp(1) | **definitive (FP)**, 0.92 |
| A | **temporal_white_noise** | independent log-random walks | **definitive (FP)**, 0.15 |
| A | constant | all abundances 1 | none |
| A | **linear_gradient** | monotone linear drift, no oscillation | **definitive (FP)**, 0.12 |
| A | **checkerboard** | period-2 alternation | **definitive (FP)**, 0.50 |
| B mate | P11_lotka_volterra | lattice predator–prey; abundances = prey, predator densities; resource = empty-site fraction (space, which prey reproduce into) | none (2 species < 3) |
| B mate | **P12_rps** | lattice cyclic dominance (3 species); resource = empty-site fraction | **screening (FP)**, 0.25 |
| B' | **p37_may_leonard** | well-mixed 3-species cyclic competition with direct interaction (May & Leonard 1975, α=0.8, β=1.2, α+β=2 → sustained LV-type cycles); resource = free capacity 1−Σx | **screening (FP)**, 0.43 |
| B' | **p37_forced_stable_coexistence** | the standalone stable-coexistence chemostat (a fixed point unforced) with seasonal supply S(t)=S(1+0.3 sin 2πt/100) | **definitive (FP)**, 0.21 |
| C | competitive_exclusion ×2 | standalone (hierarchy → 1 survivor) | none, none |
| C | stable_coexistence ×2 | standalone (5 species at a fixed point) | none, none |
| C | single_resource ×2 | standalone (R* rule → 1 survivor) | none, none |
| C | hw_uniform_consumption ×2 | HW K and supply, intransitive consumption removed → exclusion to 1 | none, none |
| C | hw_low_supply ×2 | HW K and C, all S_j = 1 → exclusion to 1 | none, none |

**False positives (10/22).** All ten pass the same two conditions: ≥ 3 species
persisting and a median late-window coefficient of variation above 0.05. The
detector's discriminating metric is "several species persist and fluctuate". It
does not test whether the fluctuation is deterministic, endogenous, oscillatory
rather than a trend, or mediated by the resource field (the `resources` key is
only checked for presence and counted). Class C stays clean because every
non-forming regime of the HW chemostat ends at a fixed point or in exclusion,
which is what the detector's two conditions reject. The negatives' mean osc_cv
(0.18) is above the positives' (0.14–0.19), so the continuous d is −0.09: on this
panel the continuous metric does not separate the pattern at all. The standalone
d = 15.09 came from a panel whose negatives were all fixed points or exclusion.

## Decisions and judgment calls

1. **New formats for P33, P34, P35, P37; P36 reuses `scalar_wealth`.** P34 and P37
   read keys no existing format provides. P33's input contract also accepts the
   existing `particles` format, but its canonical positive is a lattice director
   field (`theta_field`). Particle nulls would test only the detector's
   particle-binning path, never the representation its positives use. So Class A
   is built in the positive's schema, and the particle path is tested by the
   continuous_2d catalog mates instead (which is where it failed). P35 reads its
   periodic box only from metadata. The `particles` nulls live in a 16-box
   without a box key, so a separate `positions` format carries each frame's box.
2. **Surrogate parameters were fixed by analogy to existing formats, not tuned.**
   Examples: the Gaussian director with std = half the range, as in the phases
   format; the timeseries convention of a temporal ramp for linear_gradient and
   an alternation for checkerboard, as in choice_timeseries; mean degree and
   sizes taken from the positives. Several of these fire on P37. That is the
   finding, and no generator was changed after seeing results.
3. **P35 is scored in each substrate's own box.** `run_p35` wraps the frozen
   detector function so that a substrate whose frames carry `box_size` is
   scored with that box. Everything else, the positives included, gets the
   canonical metadata (L=14), exactly as the standalone validation did with its
   `_with_box`. Without the wrapper the harness would wrap a 16- or 20-box
   system at L=14 and create spurious neighbours. Checked out-of-TNR: scoring
   all six catalog mates at L=14 instead changes no verdict.
4. **Class B typing and cross-type mates.** P33 → continuous_2d (6 mates). P35 →
   continuous_2d (6 mates). P36 → scalar_wealth (mate P28, which needed a catalog
   generator). P34 → `evolving_network`, the ontology value P34 was built to
   fill, so no catalog pattern has it; the closest catalog phenomena, opinion
   dynamics on interaction graphs, enter as explicit mates P18 voter and P21 HK.
   P22 SIR was not added: its 64×64 lattice is a 4096-node graph, and the P34
   detector keeps a dense float copy of every frame (≈134 MB per frame). P37 →
   `population_resource`, no catalog mate of that type; P11 Lotka–Volterra and
   P12 RPS enter as explicit mates, with empty-site fraction as the shared
   resource. LV yields only 2 species, so its rejection comes from the ≥3-species
   gate.
5. **Look-alikes.** The standalone look-alike negatives were reused unchanged. Added:
   Barabási–Albert growth (P34); Bouchaud–Mézard and yard-sale-with-tax (P36),
   copying the P28 Class C module's parameters unchanged; May–Leonard (P37, the
   "LV-type oscillation" look-alike) and a seasonally forced version of the
   standalone stable-coexistence control (P37). Forcing amplitude 0.3 and period
   100 are round a-priori values. Each look-alike runs once with a harness-drawn
   seed, as every other panel's supplements do.
6. **Class C regimes chosen and excluded on ground truth.**
   - *P33, excluded:* low activity from a random start (η 0–0.10). In this model
     the ±1/2 defects do not anneal: half-defect density is 0.012–0.145 even after
     4000 steps, so the observation contains the pattern's signature. Also
     excluded: the η 0.30–0.45 crossover (S_loc 0.27–0.35).
   - *P34, used:* φ = 0.50 and 0.80. Every run with φ ≤ 0.80 ends in consensus on
     one component; fission starts at φ = 0.90. That puts this implementation's
     transition between 0.80 and 0.90. The φ = 0.99 in the model docstring is the
     positive's operating point, above it.
   - *P37, excluded:* HW with high mortality m = 0.6, because all five species
     keep oscillating (late CV ≈ 1), so the pattern still forms.
   - *P35, kept:* η_end = 0.65. It is below freezing by the hard-disk literature,
     but in this model it sits between the clearly fluid runs and the positive on
     crystalline fraction (0.13 vs 0.22). It produced the one P35 false positive.
     Dropping it would give 24/24. It is kept because it was chosen a priori from
     the literature and the judgment that it is non-forming did not depend on the
     detector.
7. **Invariance flags a priori; skipped substrates reported** (table above).
8. **Factory defaults as frozen.** `make_p33_detector_fn` defaults to 49
   permutations; the standalone used 199. The permutation null only adds +0.1 to
   the definitive-tier confidence when p < 0.01, which is unreachable with 49.
   It changes no detection or tier, only the discrete score: positives score
   0.85 instead of up to 0.95.
9. **`run_phase2a_panel.py` `__main__` guard moved to the end.** It sat above the
   P33–P37 factories, so the command-line entry point could not have reached them.
10. **Provenance.** The code was committed first (8802759), then all five panels
    were run again from that commit. The second run is byte-identical to the
    first apart from `head_commit` and `elapsed_seconds`, so the panels are
    deterministic. `class_b_for_pattern` output for all pre-existing patterns is
    unchanged, and `tests/test_phase2a_panel.py` gives the same 7 failures before
    and after (87 pass). All 7 are pre-existing and unrelated: the P18
    network-override tests, two flag-table tests, and the v1.1-schema verdict
    test.

## Runtime
VPS (8 cores, shared load ~8–12 from other jobs), one process per panel,
BLAS single-threaded. Wall clock including positive generation, with the harness's
own `elapsed_seconds` (negatives + scoring) in parentheses: P33 19 s (13 s); P34
106 s (74 s); P35 336 s (208 s); P36 21 s (15 s); P37 132 s (91 s). Catalog
natives are cached under `analysis/outputs/phase2a_catalog_cache/` (gitignored).

## Reproduce
```
cd <worktree>; export PYTHONPATH=. OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
python analysis/validation_rebuild/p33_p37_class_c_groundtruth.py all   # Class C ground truth
python analysis/run_phase2a_panel.py p33-p37                            # or p33 ... p37
python analysis/validation_rebuild/p33_p37_panel_diagnostics.py         # out-of-TNR diagnostics
```

## Files
- Harness wiring: `epc/phase2a/synthetic_discovery.py` (new), `synthetic.py`, `panel.py`,
  `catalog.py` (`DISCOVERY_CLASS_B`, `_gen_p28_yard_sale`, five adapters),
  `discovery_lookalikes.py` (new) + registration in `structured.py`,
  `detector_invariance.py`, `failed_regimes/p33..p37_*.py` (new),
  `analysis/run_phase2a_panel.py` (`build_p3X_positives`, `run_p3X`, main).
- Outputs: `analysis/outputs/p33_active_nematic_phase2a_panel_v12.json`,
  `p34_adaptive_network_…_v12.json`, `p35_entropy_crystal_…_v12.json`,
  `p36_pareto_exchange_…_v12.json`, `p37_resource_competition_…_v12.json`,
  `p33_p37_class_c_groundtruth.json`, `p33_p37_panel_diagnostics.json`.

---

## Post-panel fix: P37 detector v2 (IN-SAMPLE)

**Everything above this line is the frozen record and is unchanged.** That
includes the P37 result: 12/22, definitive on structureless noise, continuous d
−0.09, in `analysis/outputs/p37_resource_competition_phase2a_panel_v12.json`,
produced by the v1 detector at 8802759. This section reports a fix to the P37
detector made **after the frozen panel had been seen**. The fixed panel result is
therefore **in-sample**: it shows that the v2 detector rejects these particular
negatives, not that it generalises. Only P37 was changed; the P33 and P35
detectors and their false positives stay as found.

- **Code:** 0a0b1c7. Files: `epc/detectors/p37_resource_competition.py`,
  `epc/metrics/competition_dynamics.py` (new `sustained_oscillation_test`; the v1
  `coexistence_oscillation` is kept unchanged), `epc/phase2a/continuous_metrics.py`,
  `epc/phase2a/detector_invariance.py`, and `analysis/run_phase2a_panel.py`.
- **Output:** `analysis/outputs/p37_resource_competition_phase2a_panel_v12_fixed.json`
  (`head_commit` 0a0b1c7).
- **Per-negative detail:** `analysis/validation_rebuild/p37_fix_diagnostics.py` →
  `analysis/outputs/p37_fix_diagnostics.json`.

### What changed and why
The detector's own docstring and its source (Huisman & Weissing 1999, *Nature*
402:407; 2001, *Ecology* 82:2682) describe the pattern as species competing for
shared depletable resources that **coexist through sustained oscillation**. The v1
detector measured something weaker: at least 3 species persist and their
late-window coefficient of variation exceeds 0.05. In other words, persistence
plus fluctuation, which white noise, drift and alternation all share. v2 tests
the documented signature directly, over the last half of the frames:

| test | v2 definition | what it is meant to exclude |
|---|---|---|
| coexistence | ≥ 3 species persist (late mean > 1% of global peak); ≥ 3 of them with detrended amplitude std/mean > 0.05 | exclusion; fixed points; a pure linear drift (zero detrended amplitude) |
| sustained oscillation | pooled over persisting species, the share of each linearly detrended series' variance in its dominant 3-bin periodogram band. The band needs ≥ 3 cycles per window and must sit below the Nyquist band (frequency < 0.45 per frame, so period-2 alternation is excluded). Tested against 199 surrogates that permute each series independently (amplitude distribution kept, temporal order destroyed): p_osc, z_osc | white noise (no dominant band); period-2 alternation (out of band); time-shuffled oscillation |
| coupling to the resource field | max \|corr\| between detrended abundances and detrended resources at lags ≤ 2 frames, against 199 multivariate phase-randomised resource surrogates (each resource keeps its spectrum and the resources keep their mutual relations; alignment with the abundances is destroyed): p_coupling | a "resources" series that is unrelated to the abundances, or constant |

Tiers:
- **screening:** coexistence, p_osc < 0.05 and p_coupling < 0.05.
- **confirmation:** additionally p_osc < 0.01, p_coupling < 0.01 and ≥ 4 species persisting.
- **definitive:** additionally median amplitude > 0.10, and ≥ 3 species whose
  second-half amplitude is at least half their first-half amplitude ("sustained").

Null draws are 199 for both tests (`make_p37_detector_fn(n_surrogates=199)`);
the minimum attainable p is 0.005.

**The canonical continuous metric had to change.** The v1 metric
(`coexistence_oscillation`, the late-window CV) is the quantity that does not
separate the pattern (d −0.09 on the frozen panel). `CONTINUOUS_METRIC["P37"]` is
now `sustained_oscillation`: z_osc when coexistence and coupling (p < 0.05) hold,
and 0 otherwise. The v1 score is still reported in
`primary_metric["coexistence_oscillation"]` for comparison.

**The invariance flag changed with the metric.** P37's `time_shuffle_invariant` is
now False, because v2 reads temporal order. So the time-shuffled positive is now
scored as a Class A negative (Class A evaluated: 9, versus 8 in the frozen
panel). `permutation_invariant` stays True: species relabelling leaves every v2
statistic unchanged.

The positive recordings are coarsely sampled: 37.5 time units per frame, with a
dominant period of about 3–4 frames. That is why the band reaches up to a
frequency of 0.45 per frame rather than stopping at 3 frames per cycle.

### What was seen before the v2 design was fixed (so the reader can judge tuning)
Scored with the **v1** detector before or while v2 was designed:
- **Smoke test** (v1): random_uniform, permutation_shuffled, time_shuffled,
  checkerboard, P11, P12, May–Leonard, forced stable coexistence, and two Class C
  regimes.
- **The full frozen panel**, run twice (v1).
- **`p33_p37_panel_diagnostics.py`** (v1): for every negative that fired, its
  osc_cv, persistence and tier (the frozen table above), plus the skipped
  permutation and time shuffles.
- **The Class C ground-truth scan** (`p33_p37_class_c_groundtruth.py`), which reads
  the models' own late-window CV and persistence on candidate non-forming regimes.
- **I wrote every Class A generator and look-alike myself**, so the generic null
  classes the design names (white noise, red noise, drift, alternation) were
  known from that.

Scored with **v2** before commit 0a0b1c7 fixed the design: **the five canonical
positives only**. No panel negative was scored with v2 until the design was
committed. The checks on positives are kept in
`docs/validation_rebuild/p37_fix_design_checks/`:
- `inspect_pos.py`: autocorrelations and periodogram peaks of the positives.
- `pos_test2.py`: the pooled oscillation and coupling statistics on the positives.

One design change was driven by the positives. The first version required ≥ 3
species *individually* significant (each p < 0.05). That detected only 3 of the 5
positives: seeds 0 and 4 had just 2 such species each, because a 40-frame window
of a near-aliased chaotic oscillation gives a single-species test little power.
It was replaced by the pooled community test above, which tests the same
hypothesis at the level of the community. With it, all positives reach
p_osc = 0.005 (z 3.8–10.3), and all have p_coupling = 0.005: every species is
strongly anticorrelated, −0.96 to −1.0 at lag 0, with one of the resources. The thresholds (0.05 / 0.01, amplitude floors 0.05 / 0.10, ≥ 3 cycles,
Nyquist band 0.45, lags ≤ 2) were set by convention or carried over from v1. None
was adjusted after any negative was scored, and nothing changed after the fixed
panel ran.

### Fixed-panel results (v2 detector, in-sample)
**Positives:** 5/5 **definitive** (sustained_oscillation z = 3.82, 6.12, 8.67,
10.27, 3.78; p_osc and p_coupling 0.005 for every seed).

| | n | TNR | false positives |
|---|---|---|---|
| Class A | 9 (permutation_shuffled skipped) | 9/9 = 1.000 | none |
| Class B | 4 (advisory, n < 5) | 3/4 = 0.750 | p37_forced_stable_coexistence (screening) |
| Class C | 10 | 10/10 = 1.000 | none |
| overall | 23 | **22/23 = 0.957** | 1 |

- Continuous d on the **new** metric (`sustained_oscillation`): **1.06**.
- Discrete (tier-confidence) d: 11.50.
- Harness verdict: **PASS**. Class B is advisory at n = 4, so it does not gate.
- Run time: harness 106 s, about 140 s wall clock including the positives.

Which test rejected each negative (from `p37_fix_diagnostics.json`):

| negative | v1 (frozen) | v2 | rejected by |
|---|---|---|---|
| random_uniform, random_gaussian, random_binary, spatial_white_noise | definitive ×4 | none ×4 | oscillation (p_osc 0.27–0.68) |
| temporal_white_noise (log-random walk) | definitive | none | **coupling only** (p_coupling 0.29); oscillation test *passed* (p_osc 0.005, z 6.1) |
| time_shuffled (now evaluated) | skipped | none | oscillation (p_osc 0.61) |
| constant, linear_gradient | none / definitive | none | amplitude (no detrended fluctuation) |
| checkerboard (period-2) | definitive | none | oscillation (out of band; p_osc 1.0) |
| P11_lotka_volterra | none | none | coexistence (1 species persists in the window) |
| P12_rps | screening | none | **coupling only** (p_coupling 0.22); oscillation passed (p_osc 0.005, z 6.5) |
| p37_may_leonard | screening | none | **coupling only** (p_coupling 0.50); oscillation passed (p_osc 0.005, z 13.9) |
| p37_forced_stable_coexistence | definitive | **screening (FP)** | passes oscillation (z 25.5) and coupling (p 0.03); stops at screening because p_coupling ≥ 0.01 |
| all 10 Class C | none | none | amplitude or coexistence (fixed points, exclusion) |
| permutation_shuffled (skipped, out of TNR, 5 seeds) | skipped | 5/5 detected | degenerate by construction: it is the positive with species relabelled |

### Caveats a reader should weigh
1. **The oscillation test does not reject red noise.** The log-random walk passes
   it, because its low-frequency power falls inside the ≥ 3-cycle band of a
   100-frame window, and the permutation null is white. It is rejected only by
   the coupling test. The design reasoning had expected the 3-cycle floor to be
   enough; it was not.
2. **Three Class B rejections rest on the coupling test, and that test has little
   power for narrowband signals.** The look-alikes are P12 RPS, May–Leonard and
   the random walk above. Phase-randomising a near-periodic resource series
   produces another near-periodic series that aligns with the abundances by
   chance about as often as the real one does. May–Leonard's "resource" is
   literally 1 − Σx, perfectly coupled by construction, and it is still not
   significant. So those rejections should not be read as the detector
   recognising field-mediated versus direct competition; they reflect the test's
   weak power on smooth periodic data. For the broadband chaotic positives the
   same test is decisive (p 0.005).
3. **The forced look-alike remains a false positive at screening.** Externally
   forced oscillation of a stable coexistence chemostat is oscillatory, resource
   coupled, and multi-species. Nothing in the v2 signature separates endogenous
   from exogenous oscillation.
4. **The continuous d of 1.06 is modest.** The forced look-alike scores
   z = 25.5, above every positive. Checkerboard scores z = −12.3: its alternation
   is perfectly "coupled", and its band statistic sits far below the null.
5. **In-sample.** Every negative here was seen under v1, and some generic null
   classes were reasoned about explicitly during the design. A fresh negative
   panel, especially red noise without a resource field, narrowband direct-
   interaction oscillators, and forced systems, would be needed before claiming
   the fix generalises.
6. **The battery's generic emergence channel was not changed.**
   `epc/phase2a/emergence.py` still uses the v1 `coexistence_oscillation` (the
   coexistence-oscillation channel), so that channel keeps v1's non-specificity.
   It is out of scope here and flagged for follow-up.
7. **Re-running `analysis/run_phase2a_panel.py p37` now uses v2 and writes the
   `_fixed` file.** To reproduce the frozen P37 record, check out 3163bcd. Note
   that `p33_p37_panel_diagnostics.py` reads the P37 output path, so a re-run
   would now describe the fixed panel.
