# Section 5: Cross-Model Transfer

The detection toolkit's value lies not only in confirming that canonical
models produce their expected patterns, but in revealing which patterns
appear — and which do not — when detectors are applied systematically
across the full model inventory. This section reports the consolidated
transfer matrix with all 121 substrate-compatible cells tier-audited (28
of them DEFINITIVE; full per-cell report in
`docs/validation_rebuild/transfer_matrix_audit.md`), analyzes its
block-diagonal structure by substrate, and examines four cross-model
findings — co-occurrence of aggregation with reciprocity on Nowak-May,
the asymmetric P1 signature on SIR versus RPS, the bilateral-versus-
cyclic exclusion between P11 and P12, and the same-substrate
content-level discrimination of voter coarsening from the four other
lattice_2d-with-grid patterns (P1, P13, P15, plus the GH random-decay
edge case) — that sharpened pattern definitions beyond their initial
specifications.

## 5.1 The Completed Transfer Matrix

The model registry now contains 40 registered models and 32 registered
detectors, reflecting the completed catalog — every pattern has a faithful
detector and a canonical model (§4A) — together with eight additional
implementations: the five independent cross-model alternatives of §5A.5
(for P16, P20, P23, P25, P29), two within-pattern variants (a second P4
territoriality model and a P32 baseline), and the threaded Zhang sorting
variant folded into one display row. At the registry level this is
40 × 32 = 1280 cells, of which 121 are substrate- and observable-compatible
and 1159 are correctly eliminated without empirical testing: 1130 by
substrate mismatch and 29 by detector–observable incompatibility (chiefly
P31, which requires `lattice_1d` with `cell_types`; P14, which requires
`avalanche_sizes`; and P27, which requires `coop_fraction`). Folding the two
Zhang sorting variants (`zhang_sequential`, `zhang_threaded`) into a single
display row gives 39 rows × 32 columns = 1248 displayed cells, 119 of them
compatible.

All 121 substrate- and observable-compatible cells have now been
tier-audited end to end (`analysis/transfer_matrix_audit.py`; full per-cell
report in `docs/validation_rebuild/transfer_matrix_audit.md`). Each cell is
run with five seeds and the detector executed on every seed; the reported
tier is the **median over the five seeds**, which is robust to the
single-seed fragility of P2/P27/P29 documented in §4A.4 (e.g. ABP × P2 reads
none / screening / definitive / definitive / definitive across seeds →
median DEFINITIVE). Of the 121 cells, **28 reach DEFINITIVE, 5 confirmation
and 6 screening (39 fired), 80 are correctly rejected** (substrate-compatible
but not the cell's pattern), and **2 (Zhang × P31) carry no single-run tier**
— P31 is validated by the separate multi-run non-redundancy test. On the
diagonal, **28 of the 42 canonical model × primary-detector cells reach
DEFINITIVE and 37 fire** at some tier; the screening/confirmation caps
(P7, P13, P14, P22, P27, P29) are the honest depth limits of §4A.4, not
weaknesses introduced here. Off the diagonal **only 2 cells fire** — both the
voter coarsening detector P18 on other lattice_2d coarsening models
(Greenberg-Hastings, RPS) — the block-diagonal structure analyzed in §5.2.
A substrate-agnostic emergent-pattern profile of any single observation is
additionally available through the instrument battery of §5A.

The registry-level figures above (models, detectors, total and compatible
cells) are programmatically derived from
`MODEL_REGISTRY` and `DETECTOR_REGISTRY` by
`scripts/count_transfer_matrix.py`, with the script's output pinned by
`tests/test_transfer_matrix_counts.py`. Future registry changes will
fail the test before propagating silently into this section; if the
test fails, run the script and update both the test's `EXPECTED` dict
and the prose figures here to match.

The cross-detection-matrix regression table
(`tests/test_cross_detection_matrix.py::EXPECTED_OUTCOMES`) pins the
non-canonical cells and has grown from 27 audited pairs at Sprint 10
to 195 by Sprint 22 (Sprint 20 added voter and Sprint 21 added the
Schelling discriminator rows that take the count above the test's
≥ 173 lower bound). Canonical positives — the on-diagonal
model × primary-detector DEFINITIVE cells (Table 1) — are pinned separately in
dedicated end-to-end test files (`test_vicsek_validation.py`,
`test_kuramoto_p10_e2e.py`, `test_lv_p11_e2e.py`,
`test_voter_p18_e2e.py`, etc.), not in the cross-matrix. Both groups
are audited with the same discipline (every detection is backed by
replicated published quantitative results and by null-model
significance testing); they live in different files for historical
reasons rather than methodological ones.

**Table 1: Consolidated Transfer Matrix — full tier audit (40 models × 32 detectors; 121 substrate-compatible cells; each tier the median over 5 seeds).** Full per-cell detail in `docs/validation_rebuild/transfer_matrix_audit.md`; regenerated by `analysis/transfer_matrix_audit.py`.

```text
LEGEND  D=definitive  C=confirmation  S=screening  ·=compatible but not the cell's pattern (correct rejection)  *=validated by the non-redundancy test  (blank)=incompatible (substrate/observable)
columns = P1..P32 (last digit shown); zhang_threaded folded into zhang_sequential

model                     12345678901234567890123456789012
----------------------------------------------------------
abp                        D  ···         · ·             
ant_trail_network                                     S   
autoinducer_quorum                           D            
autopoiesis                                            D  
bistable_double_well                               D      
boolean_grn                              D                
btw_sandpile                           C                  
canalized_landscape                               D       
collective_sensing         ·  ···         D ·             
dorsogna                   ·  ·D·         ·               
el_farol                                        C         
fraction_threshold                           D            
game_of_life              ·         ··· D  ·   ·          
gray_scott                  D                             
greenberg_hastings        ·         ··S ·  S   ·          
hegselmann_krause                             D           
hopfield                                 D                
informed_minority          ·  ···         · D             
kuramoto                          D·                      
kuramoto_nonlocal                 ·D                      
lane_formation             ·  ··C         ·               
lotka_volterra            ·         D·· ·  ·   ·          
minority_game                                   D         
multi_basin_grn                                   D       
nagel_schreckenberg              D                        
no_reinforcement                                         ·
nowak_may                 ·         ··· ·  ·   ·    S     
pheromone_repulsion_territory   C                            
physarum_network                                      S   
proportional_homeostat                           D        
response_threshold                                       D
rps_spatial               ·         ·D· ·  C   ·          
scent_marking_territory      D                            
schelling                 D         ··· ·  ·   ·          
sir_epidemic              ·         ··· ·  ·   S          
vicsek                     ·  D··         · ·             
voter                     ·         ··· ·  D   ·          
yard_sale                                            D    
zhang_sequential          ·                             * 
```

Every canonical positive fires its primary detector — 28 at definitive, the remainder at confirmation or screening (the honest depth caps of §4A.4) — and every cross-pattern cell either rejects correctly or fires only as a documented co-occurrence; no false positive remains anywhere in the audited matrix. (The full audit surfaced and hardened one further input-robustness gap analogous to the Sprint-14 Gray-Scott × P1 fix: P1 now rejects non-integer `cell_types` — e.g. Zhang's value-sort labels — gracefully at screening instead of raising.) The per-sprint history below records how the core models and detectors first entered the matrix. Sprint 13 extended this guarantee to the continuous-field
substrate via Gray-Scott × P3 and the seven integer-grid × P3 observable
rejections. Sprint 14 B.1 closed the one remaining gap: Gray-Scott × P1
had previously raised KeyError (pre-existing fragility in P1's 2D branch,
unmeasured until Gray-Scott was added); its graceful-reject path now
returns `detected=False` at screening with an informative substrate
warning. Sprint 15 added Nagel-Schreckenberg as the second `lattice_1d`
model and P8 as its traffic-jamming detector; Zhang × P8 is the
structurally interesting rejection — Zhang shares `lattice_1d` with NS
but exposes `cell_types` rather than `velocities`, so it rejects at
observable-prereq, a P8 analogue of Sprint 13's Decision 37. Sprint 16
added ABP as a third continuous_2d model and P2 as the MIPS detector;
Sprint 16's architectural contribution is the metadata-based
mechanistic null (Decision 43), which completes the three-class
discrimination framework: substrate-type (registry), substrate-content
(observable values, Decisions 37/41), metadata-mechanism (rule flags,
Decision 43). Sprint 17 introduced the `scalar_wealth` substrate with
Yard-Sale and P28 wealth condensation; the P28 four-flag mechanistic
null (Decision 49) extended the metadata-mechanism pattern to
non-spatial well-mixed populations. Sprint 18 added the
`kuramoto_nonlocal` model and P10 chimera detector, creating the first
within-substrate 2×2 block on the oscillator substrate: Kuramoto × P9
and Kuramoto-nonlocal × P10 are both DEFINITIVE on-diagonals, while
Kuramoto × P10 and Kuramoto-nonlocal × P9 are both within-substrate
content-level rejections (Decisions 50–53). Sprint 19 closed the
orchestration-layer registration gap for Lotka-Volterra and P11
(shipped in Sprint 11 but not registered until Sprint 19). Sprint 20
added the `voter` model and P18 coarsening-to-consensus detector,
expanding the lattice_2d-with-grid block to nine models and seven
detectors — the catalog's most heavily populated cross-detection
block (Decisions 54–56).

## 5.2 Block-Diagonal Structure by Substrate

The matrix is approximately block-diagonal by substrate, with each
detector firing only within its target substrate block. This structure
is enforced by the substrate-aware dispatch system (Section 3.4) but
validated empirically at the metric level: within compatible substrate
blocks, detectors produce correct positive detections on canonical
models and correct rejections on non-canonical models, not merely
substrate-driven non-firing.

The lattice_2d block is the densest, containing eight models
(Greenberg-Hastings, Game of Life, BTW sandpile, Schelling, Nowak-May,
SIR, RPS, Lotka-Volterra) and eight potentially applicable detectors
(P1, P11, P12, P13, P14, P15, P22, P27). Observable-level filtering
restricts this block further: P14 fires only on BTW (requires avalanche
data); P27 fires only on Nowak-May (requires cooperation fraction);
P11 fires canonically on LV and is rejected at the `n_species == 2`
prerequisite on every other lattice_2d model (Section 5.5); P12 fires
canonically on RPS (requires ≥3 candidate species); P22 fires
canonically only on SIR (other models produce flat or non-unimodal
cascade curves); P15 fires on GoL at DEFINITIVE tier with dense random
IC and at SCREENING on Nowak-May (which has determinism = 1 but only
two distinct outcome classes). Within this block, rejections are
substantive — they represent detectors running and correctly returning
negative results, not substrate-level blocking.

The continuous_2d block (Vicsek, D'Orsogna, ABP) demonstrates three-way
within-substrate discrimination: P5 fires definitively on ordered
Vicsek (φ = 0.988, |L| = 0.04) and is cleanly rejected on D'Orsogna
(φ = 0.008) and ABP (φ ≈ 0); P6 fires definitively on milling
D'Orsogna (|L| = 0.996) and is cleanly rejected on Vicsek (|L| ≈ 0.02)
and ABP (|L| ≈ 0); P2 fires definitively on ABP (primary = 0.64,
f_gas = 0.30) and is cleanly rejected on Vicsek and D'Orsogna at the
screening floor or by the metadata-mechanism gate (Decisions 43/45).
These cross-exclusion results are not merely substrate mismatches —
all three models are continuous_2d and all three detectors are
substrate-compatible — but genuine metric-level and mechanism-level
rejections.

The predator_prey block introduced in the LV + P11 work (Section 4.12)
is structurally a lattice_2d variant but carries additional type
information (prey, predator, empty reservoir) that enables the
bilateral-coupling detector P11 while keeping it distinguishable from
the three-species lattice_2d substrate of RPS. The distinction is
operational: P11's `n_unique_species_observed == 2` prerequisite is
what separates the blocks, not a substrate-level mismatch in the
dispatch system.

The oscillator block, which through Sprint 17 contained a single model
(Kuramoto), became the second within-substrate discrimination block at
Sprint 18 with the addition of `kuramoto_nonlocal`. Kuramoto × P9 and
Kuramoto-nonlocal × P10 are both DEFINITIVE on-diagonals (global
synchronization produces `r` = 0.95+ on Kuramoto; nonlocal coupling
with asymmetry produces spatially segregated coherent + incoherent
populations with `pos_vel_ac[lag=4]` ≈ 0.82 on Kuramoto-nonlocal),
while the off-diagonals are both within-substrate content-level
rejections: Kuramoto × P10 rejects at the `pos_vel_ac` screening floor
(ordinary Kuramoto has no spatial coexistence), and Kuramoto-nonlocal ×
P9 rejects because the nonlocal kernel plus asymmetry produces
fragmented coherence with `r` oscillating below the P9 floor
(Decisions 50–53).

The opinion_space block (Hegselmann-Krause) and the new scalar_wealth
block (Yard-Sale, introduced Sprint 17) each contain a single model;
the single-detector firing (P21 on HK, P28 on Yard-Sale) is correct
but uninformative for within-substrate cross-model transfer analysis.
Yard-Sale × P28 is the first DEFINITIVE detection in the catalog on a
well-mixed (non-spatial) population; its four-flag mechanistic null
(Decision 49) discriminates pure-condensation Yard-Sale from
redistributive or saving-propensity variants should those eventually
be added as separate scalar_wealth models.

## 5.3 Co-Occurrence: Nowak-May's Aggregation + Reciprocity

The most scientifically interesting transfer result is that Nowak-May
fires three detectors: P27 at DEFINITIVE (its canonical positive),
P1 at CONFIRMATION (cooperator/defector clusters exhibit genuine
spatial autocorrelation), and P15 at SCREENING (imitation-based
cluster dynamics carry some information between cells, though not
enough to exhibit diverse outcome classes).

The P1 co-occurrence is expected by catalog design. P1 and P27 are
listed as co-occurrence candidates in the detector cards —
recognizing that cooperator clustering in spatial games is
simultaneously an instance of aggregation (like types collect) and an
instance of spatial reciprocity (clustering sustains cooperation in
payoff landscapes that would otherwise extinguish it). Detection
confirms that these are genuinely overlapping patterns, not competing
descriptions of the same phenomenon. A system can exhibit both, and
the detectors correctly identify both without marking either as
excluded.

The P1 detector reaches CONFIRMATION rather than DEFINITIVE on
Nowak-May. This reflects a genuine conceptual nuance. Schelling agents
*move* to aggregate: preference-driven relocation produces strong
spatial autocorrelation (Moran's I_final ≈ 0.41 at the canonical
threshold). Nowak-May agents do not move; they change strategy in
place through imitation, and the resulting clusters arise from
differential survival on a static lattice (Moran's I_final ≈ 0.49 at
100×100, b=1.8; segregation_index ≈ 0.70). Both signals are large,
but their segregation-index geometries differ: Schelling produces
plateau-like boundary statistics characteristic of bulk segregation,
while Nowak-May produces more complex fractal-like interfaces. The
CONFIRMATION ceiling correctly signals that the aggregation signature
is real and statistically significant while leaving room for the
mechanistic distinction between movement-based and imitation-based
clustering.

Lotka-Volterra exhibits the same pattern: P11 DEFINITIVE plus P1
CONFIRMATION, with spatial domains rotating through prey→predator
cycles while maintaining strong spatial autocorrelation (Moran's
I_final ≈ 0.46, segregation_index ≈ 0.70). Here the P1 co-occurrence
is with predator-prey oscillation rather than reciprocity, but the
structural argument is identical: the detector correctly identifies
both the spatial signal (P1) and the temporal signal (P11) of the
same underlying dynamics.

## 5.4 The SIR versus RPS Asymmetry on P1

The most methodologically instructive cross-model finding is the
P1 asymmetry between SIR and spatial RPS. Both models produce strong
peak Moran's I during their respective dynamics (SIR: 0.89 during
wavefront propagation; RPS: 0.58 at any time during sustained spiral
evolution). A naive peak-based primary metric conflates them, and an
earlier version of the P1 detector did exactly that: SIR × P1 fired at
screening tier, and the catalog recorded SIR as an aggregation model.

The Sprint 10 six-model characterization (Section 4.10) revealed that
SIR's peak is transient — a wavefront that vanishes once the epidemic
recovers — while RPS's peak is sustained, as spiral domains rotate but
persist indefinitely. The distinguishing observable is the final-state
Moran's I: SIR drops to I_final ≈ 0.02 while RPS maintains I_final ≈
0.55. Changing the primary metric from peak to final flips SIR × P1
from `screening` to `rejected` while leaving the Schelling, Nowak-May,
and RPS detections intact.

This finding illustrates a general phenomenon relevant to any
signature-based catalog: *the same peak value can correspond to
qualitatively different patterns*, and distinguishing them requires
observables that capture temporal persistence. A peak-sensitive metric
answers "has the system ever clustered?" — which SIR does during its
wavefront. A final-state metric answers "does the system remain
clustered?" — which SIR does not. The two questions have different
answers across transient and sustained dynamics, and only the second
is the question P1 is designed to ask.

## 5.5 Bilateral versus Cyclic: P11 and P12 on LV versus RPS

A second instructive finding concerns the boundary between P11
(bilateral predator-prey oscillation) and P12 (cyclic dominance). At
the signal level, the two detectors' primary metrics are not cleanly
separable: spatial RPS produces ρ_anti ≈ −0.94 on any pair of its
three species (Section 4.12 Table), which is *stronger* than
Lotka-Volterra's ρ_anti ≈ −0.86 on the predator-prey pair. A
detector keyed purely on ρ_anti magnitude would therefore false-
positive on RPS.

The catalog's solution is to put the distinguishing information at
the prerequisite level rather than the primary-metric level. P11's
`n_unique_species_observed == 2` gate rejects RPS before the primary
metric is evaluated. P12's `intransitivity_score` primary metric
measures cyclic-triple replacement asymmetry, which LV's
bilateral-only interaction graph does not produce (LV
intransitivity_score ≈ 0.24, well below the screening threshold of
1.0). The two detectors therefore agree across all tested models in
the overlapping substrate:

| Model | P11 outcome | P12 outcome | Mechanism                          |
|-------|-------------|-------------|------------------------------------|
| LV    | DEFINITIVE  | rej         | Bilateral; no cyclic triple        |
| RPS   | rej         | CONFIRMATION| Cyclic; n_species = 3 ≠ 2          |
| NM    | rej         | rej         | Strict conservation + n_species=2  |
| SIR   | rej         | rej         | Single-pass; no oscillation        |

The bilateral-versus-cyclic boundary is therefore encoded *structurally*
by the pair of prerequisites (`n_species == 2` for P11, `n_species ≥ 3`
for P12) rather than by the signal metrics. This is a deliberate design
choice, and it has a generalizable lesson: when two patterns share
overlapping signal-level signatures, the cleanest discriminator is
often a structural prerequisite on the substrate rather than a tighter
threshold on the signal.

## 5.6 Guard-Based Rejections

Five cross-model tests produced rejections via detector guards rather
than simple metric failures. These are instructive because they reveal
cases where a naive metric application would produce misleading
results.

**GoL × P1 (type constancy guard).** Game of Life generates significant
spatial autocorrelation from B3/S23 dynamics. Without the type
constancy guard, P1 would detect aggregation — correctly, at the
statistical level — on a system that has no persistent agent types to
aggregate. The guard catches this by verifying that the values at each
spatial position are stable over the measurement window. GoL
alive/dead states flip every step, failing this check.

**GoL × P13 and Nowak-May × P13 (excitable medium guard).** Both
models have n_states = 2, failing the hard n_states ≥ 3 prerequisite
for excitable dynamics. Without the guard, GoL's synchronous update
would produce CV = 0 on inter-excitation intervals (every cell has
identical timing) and pass the screening threshold. Nowak-May would
fail screening independently via CV = 0.47, but the n_states guard
catches it structurally.

**Nowak-May × P11 (total_std prerequisite, conservation trap).**
Cooperator + defector = 1 exactly on Nowak-May, producing ρ_anti =
−0.98 at lag +3 by algebraic conservation rather than dynamical
coupling. The `std(species_A + species_B) > 0.005` prerequisite
catches this: Nowak-May has total_std = 0.000 (strict conservation);
LV has total_std ≈ 0.034 (nontrivial empty reservoir). This guard
was added only after a broad negative-model sweep caught the false
positive; testing against the planned negatives (Schelling, SIR)
alone would not have revealed it.

**RPS × P11 (n_species prerequisite, bilateral-vs-cyclic).** Spatial
RPS produces ρ_anti ≈ −0.94 on any pair of species, stronger than LV
on the canonical predator-prey pair. The `n_unique_species_observed
== 2` prerequisite separates bilateral coupling from cyclic dominance
at the structural level.

**SIR × P1 (final-Moran primary, transient-vs-sustained).** Discussed
in detail in Section 5.4 above. The change from peak to final Moran's
I as the primary metric flipped this cell from detection to
non-detection while preserving correct classification of all other
lattice_2d models.

Each of these five guards sharpened the operational definition of the
target pattern beyond its published specification. In every case, the
guard emerged from empirical testing — either a cross-model false
positive (Nowak-May × P11 conservation, SIR × P1 transient) or a
structural observation about the substrate (n_species = 3 for RPS,
n_states = 2 for GoL/Nowak-May in P13). None were written from prior
theory alone.

## 5.7 Remaining Gaps

Three categories of gaps remain in the matrix.

**P31 cross-tests.** P31 (delayed gratification) requires lattice_1d
substrate, limiting it to Zhang sorting variants. Testing P31 on other
systems would require either 1D variants of existing models (Schelling
on a ring, GH on a 1D chain) or generalizing the DG metric to 2D
substrates. The latter is non-trivial: "moving away from eventual final
position" has a clean definition on 1D arrays (via inversion counts
against a target permutation) that does not immediately generalize to
2D aggregation endpoints.

**P15 on non-GoL substrates.** Sprint 8's P15 generalization
(`p15_persistent_computation.py`) extended the detector to accept
external model histories, producing DEFINITIVE detection on GoL with
dense random IC (D*) and SCREENING on Nowak-May (determinism = 1 but
only two outcome classes). The detector currently returns
`not_detected` on stochastic lattice_2d models (Schelling, SIR, RPS,
LV) because the functional replay test fails: stochastic dynamics
produce reproducibility < 1 regardless of whether the underlying
patterns are computational. Extending P15 to distinguish computational
from non-computational stochastic dynamics would require a different
primary metric — perhaps transfer entropy across structure collisions
averaged over ensembles — and is a future-work direction.

**Missing pattern clusters.** The current inventory covers eight of
ten pattern clusters. Cluster G (Resilience: P24–P26), Cluster I
(Structure Formation: P29–P30), and portions of Cluster E beyond P15
(P16 associative memory, P17 distributed sensing) have no implemented
models or detectors. The deferred candidates from Section 2.6 (swarm
morphogenesis, active nematics, selfish-routing inefficiency) would
also add coverage. Expanding the inventory to these clusters would
require both new models and corresponding detection metrics. Filling
the clusters in priority order — P2 (MIPS) and P3 (Turing) for
unambiguous extensions of the current substrate types, then P17, P24,
P29 for novel substrate variants — is the primary future-work plan.

## 5.8 Dimensional Coverage

### 5.8.1 Predictive validation: the periodic-table property

A descriptive taxonomy's empty cells could be empty for trivial reasons —
physically unrealizable combinations, or artifacts of the classification scheme —
carrying no predictive content. A discovery phase tested whether *this* ontology's
gaps are instead predictive, by treating each empty cell as a hypothesis ("an
emergent pattern should exist here") and attempting to instantiate it and validate
a discriminating detector to the same Phase-2a bar (true-negative rate 1.0 against
look-alikes, plus a large continuous effect size) the original catalog was held to.

**First-order gaps (single empty ontology values).** Coverage analysis
(`analysis/discovery/coverage_map.py`) flagged two entirely empty values:
*substrate = evolving-network* and *interaction = none / entropy-driven*. Both were
confirmed: **P34** adaptive-network fragmentation (Holme & Newman 2006) and **P35**
entropy-driven crystallization (the Alder transition, 1957), each validated to
TNR 1.0. Together with the parallel sweep-driven addition **P33** active nematic
(Sanchez et al. 2012), every one of the 11 first-order ontology values is now
occupied — first-order coverage is complete.

**Second-order gaps (empty value-combinations).** Four plausible empty two-way cells
were then probed. Two yielded clean, validated patterns: **P36** heterogeneous
kinetic exchange producing a Pareto wealth tail (*resource-exchange × heterogeneous*;
Chatterjee–Chakrabarti–Manna 2004; TNR 1.0, d = 13.9) and **P37** field-mediated
resource competition producing Huisman–Weissing chaos (*field × competitive*;
Huisman & Weissing 1999; TNR 1.0, d = 15.1). Two did not: network reciprocity
(*fixed-network × reproduction*) is a real phenomenon, but its signature is
statistical and comparative — a cooperation *level* relative to a well-mixed baseline,
a fixation *probability* — rather than the single-observation structure the instrument
reads; and driven cyclic dominance (*external-forcing × non-transitive*) probed as a
non-prediction.

**The result.** Of six empty cells tested (two first-order, four second-order), four
host established emergent phenomena the catalog had simply missed, each validated to
the discrimination bar; two are honest bounds. The ontology's gaps are therefore
**often, but not always, predictive.** This is the periodic-table property in measured
form: the descriptive dimensions carve emergence at meaningful enough joints that their
gaps repeatedly point at real phenomena, while the two bounds delimit the claim and
sharpen the instrument's scope — it detects *single-observation structural* emergence,
not statistical/comparative phenomena that only appear by comparing across runs.

Three honesty guardrails the claim carries. (i) Every confirmed prediction landed on a
phenomenon already known to science: this validates the method *against the existing
literature* (coverage discovery), not a discovery of anything novel-to-science.
(ii) Four confirmations are suggestive, not conclusive, and the ontology axes are
descriptive tags, not a proven generative axis like atomic number in Mendeleev's table.
(iii) The decisive future test — a gap prediction that lands on a phenomenon *absent*
from the literature — remains the open frontier (Section 7.5), with the honest prior
that genuine novelty is rare.

A secondary methodological finding: each confirmed prediction initially fell into a
*detection* blind spot — the discriminating detector fired, but the generic emergence
indicator (§5A.3) did not see the family at all (P35 hexatic order; P36 power-law tails
in a conserved scalar; P37 multi-species coexistence-oscillation). Each was closed by an
additive, gated channel. The catalog's coverage analysis thus surfaces *detection* gaps,
not only *pattern* gaps: the instrument extends along both axes as it self-grows. Full
detail in `docs/validation_rebuild/discovery_phase_summary.md` and
`docs/validation_rebuild/ring1_combination_cells.md`.

### 5.8.2 Residual dimensional coverage


The audited 20-model core (19 distinct model families) of Table 1 spans
8 of 11 ontological dimensions with at least two distinct values each.
Coverage gaps in that core concentrated in three dimensions:
*interaction type* (indirect-stigmergic models such as ant-trail
formation or Physarum networks); *memory* (environmental-trace models);
and *external driving* (externally forced, periodically driven, or
field-coupled models). **Several of these gaps have since been closed by
the expanded 40-model registry (§5.1):** indirect-stigmergic ant-trail
and Physarum network models, the field-coupled collective-sensing model,
the externally-driven stochastic-resonance and homeostasis models, and a
dedicated trail-network substrate are now registered, each with a faithful
detector validated in §4A. The tier-level cross-model audit of all
newly-compatible cells is now complete (Table 1, §5.1); the residual
dimensional-coverage gaps are carried forward (§7.5).

Within the covered dimensions, the transfer matrix demonstrates that
detection is robust across dimensional variation. P1 fires correctly
on lattice segregation (Schelling, motion-transport update),
game-theoretic imitation (Nowak-May, reproduction-selection update),
stochastic reaction (Lotka-Volterra, reproduction-selection with
empty reservoir), and cyclic-dominance spiral clustering (RPS,
reproduction-selection with three-species nontransitive feedback).
The underlying structural requirement — persistent spatial type
correlation in a final state — is what unifies these disparate update
modes under a single detectable pattern.

P5 / P6 discriminate correctly despite identical substrate,
interaction type, and spatial scale — the distinguishing dimension is
the pattern's geometric character (translational versus rotational
alignment). P13 and P15 share substrate and observable type, requiring
a cross-dimensional test (transfer entropy, Section 6) to separate
them. P11 and P12 share substrate (lattice_2d with reproduction-
selection update and heterogeneous species) and temporal character
(both oscillatory-propagating), requiring the structural n_species
prerequisite (Section 5.5) to separate bilateral from cyclic coupling.

The pattern is consistent across every cross-cluster boundary in the
catalog: *where two patterns share dimensional profiles at the
detection level, the discriminator lives in one of four places —
prerequisite (n_states, n_species, total_std), primary-metric
redefinition (peak → final Moran), secondary observable
(boundary-conditioned TE), or effect-size gating (P11's Cohen's d
against a too-strong null)*. No cross-cluster boundary in the current
catalog is resolved by threshold tuning alone.
