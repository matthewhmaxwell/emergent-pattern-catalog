# Detector Specification Cards v0.6.0

Bridge document from pattern taxonomy (v0.4) to detection toolkit.

v0.6.0 (Sprint 16): P2 card completely rewritten to match implemented
detector (ABP + P2MIPSDetector). The pre-existing recipe based on
Hartigan dip on density histograms was shown empirically (Phase 1c/d)
to false-positive every discrete-count distribution — dip is a bad
primary for this substrate. Replaced with two_phase_coexistence_score
= min(f_gas, f_liquid) plus a three-part confirmation gate (r band,
CV_v, frac_stalled) and a metadata-driven mechanistic null for
DEFINITIVE. See ADR 43-46 and REPLICATION_NOTES.md Sprint 16.

v0.5.5 (Sprint 15): P8 card rewritten to match implemented detector
(Nagel-Schreckenberg + P8TrafficJammingDetector). Calibration table at
L=1000, v_max=5, p_slow=0.3 regime. Density-saturation false-positive
trap discriminator documented.

v0.5.4 (Sprint 13): P3 card rewritten to match implemented detector
(Gray-Scott + P3TuringWavelengthDetector). New lattice_2d_continuous
substrate documented.

v0.5.3 (Sprint 11): P11 card rewritten to match implemented detector
(primary metric is `rho_anti` at nonzero lag, not PSD Q-factor).
Empirical measurements from LV lattice canonical positive added.

---

## Detector Output Schema

Every detector returns a structured result, not just pass/fail:

```python
@dataclass
class DetectorResult:
    pattern_id: str              # e.g., "P5"
    detected: bool               # screening-tier pass
    tier: str                    # "screening" | "confirmation" | "definitive"
    confidence: float            # 0.0–1.0, capped by tier
    primary_metric: dict         # e.g., {"polarization": 0.84}
    secondary_metrics: dict      # supporting measurements
    effect_size: dict            # e.g., {"cohens_d": 1.2, "I_vs_null": 0.65}
    null_p_value: float          # p-value from primary null model
    null_type: str               # "shuffle" | "surrogate" | "mechanistic_intervention"
    exclusions_checked: list     # neighbor patterns tested
    exclusion_results: dict      # e.g., {"P6": "excluded", "P7": "excluded"}
    co_occurrence_candidates: list  # compatible patterns also detected
    metadata_available: bool     # whether model metadata was used
    warnings: list[str]          # e.g., ["N < 50", "run length < 10τ", "metadata missing for structural check"]
    notes: str                   # free text for edge cases
```

**`effect_size`** enables cross-model comparison. Recommended entries:
primary metric value minus null-model mean, divided by null-model SD (Cohen's
d equivalent). Include raw metric values for transparency.

**`null_type`** distinguishes null-model evidence strength:
- `"shuffle"`: label/timing permutation of observed data. Available from any state history.
- `"surrogate"`: synthetic data matching marginal statistics (e.g., phase-randomized surrogates). Available from state history.
- `"mechanistic_intervention"`: requires modifying model rules (e.g., zero-coupling, constant-speed, open-loop). Only available in simulation contexts.

**`warnings`** flags: low N, short run length relative to τ, boundary artifacts,
metadata missing for structural checks, approaching screening threshold, etc.

**Confidence score** is tier-capped:

| Tier | Max confidence | Base | Bonuses |
|---|---|---|---|
| Screening | 0.60 | 0.35 if primary passes | +0.15 secondaries, +0.10 shuffle null p < 0.01 |
| Confirmation | 0.85 | 0.55 if confirmation criteria met | +0.15 null p < 0.001, +0.10 effect size > 1.0, +0.05 all secondaries |
| Definitive | 1.00 | 0.75 if definitive criteria met | +0.10 all exclusions cleared, +0.10 both null types rejected, +0.05 finite-size/robustness checks |

Confidence never exceeds the cap for its tier. A screening result cannot
reach 0.85 regardless of how strong the primary signal is — it must pass
confirmation criteria to unlock the higher range.

---

## Detection Tier Definitions

| Tier | Purpose | Requirements |
|---|---|---|
| **Screening** | Exploratory flagging. Catches candidates at cost of some false positives. | Primary metric above relaxed threshold. |
| **Confirmation** | Standard detection. Suitable for cross-model comparison studies. | Primary + ≥1 secondary metric. Null-model p < 0.01. Scale-normalized persistence. |
| **Definitive** | Publication-grade. Suitable for claiming pattern presence in a paper. | All confirmation requirements + all nearest-neighbor exclusions passed + mechanistic and shuffle nulls both rejected. |

---

## Scale-Normalization Convention

Fixed timestep thresholds are fragile across implementations. All persistence
and duration requirements use system-intrinsic timescales:

| System type | Timescale unit | Symbol |
|---|---|---|
| Oscillators | Complete oscillation period | T_osc |
| Cellular automata | Grid width / wavefront speed | T_prop = L / v |
| Agent motion (continuous) | Domain crossing time = L / v_mean | T_cross |
| Agent motion (lattice) | Mean steps to traverse grid = L² / swap_rate | T_traverse |
| Sorting algorithms | Total convergence time (steps to sorted state) | T_sort |
| Selection/evolution | Generation time | T_gen |
| Perturbation-recovery | Relaxation time of unperturbed system | T_relax |

When a card says "≥ 10 τ", τ is the appropriate timescale from this table.

---

## Observable Scope Categories

Each card is labeled with one of:

- **State-history only:** All metrics computable from observed state sequences.
  Pattern can be detected in real-world or unknown-mechanism data.
- **Model-metadata assisted:** Primary detection works from state history;
  some secondary metrics or structural checks benefit from model parameters.
- **Model-metadata required:** Primary detection requires knowledge of model
  rules, parameters, or interaction structure that cannot be inferred from
  state history alone.

---

## Conventions

- Significance α = 0.01 unless noted.
- Permutation tests: ≥ 999 permutations.
- Power-law fitting: MLE + KS (Clauset et al. 2009), never log-log regression.
- Bimodality: Hartigan's dip test unless noted.
- "τ" in thresholds = system-appropriate timescale (see table above).

---

## Cluster A: Spatial Organization

### P1 — Similarity-driven aggregation

**Observable scope:** State-history only

**Required raw observables:**
- `positions[t]`: (N,) grid indices or (N, 2) (row, col)
- `types`: (N,) integer labels (constant, ≥ 2 distinct values required)
- `grid_dims`: (rows, cols) for 2D models; `n` for 1D chimeric arrays

**Preprocessing:**
1. Construct adjacency (1D: left/right; 2D: 4- or 8-connected neighbors).
2. Type-occupancy grid G(t).
3. Same-type neighbor fraction `s_i(t) = (# same-type neighbors) / (# total neighbors)`.

**Primary metric: Final-state Moran's I**

`I_final = (N/W) · Σ w_{ij}(x_i−x̄)(x_j−x̄) / Σ(x_i−x̄)²`

computed on the **terminal state** `state_history[-1]`.

Sprint 10 change: the primary was previously `max(peak, final)` over the
trajectory. That metric falsely classified models with transient
wavefronts (SIR epidemic, traveling excitation waves) as P1 co-
occurrences because peak Moran's I during the wavefront can be very
high (≈ 0.9) even when the final state is uniform (≈ 0.02 after the
epidemic dies out and every cell is recovered). Final-state Moran's I
correctly distinguishes:

- **Genuine sustained aggregation** (Schelling, Nowak-May b=1.8): final ≈ peak.
- **Rotating-but-persistent clustering** (RPS spatial spirals): final ≈ peak
  ≈ 0.55 (spirals rotate but the spatial structure persists at every step).
- **Transient wavefronts** (SIR epidemic): peak ≈ 0.89 but final ≈ 0.02 —
  correctly rejected.

See REPLICATION_NOTES.md (Sprint 10 section) for the 6-model empirical
characterization that motivated this change, and PROJECT_STATUS.md
Decision 32 for the architecture rationale.

**Retained diagnostic fields in `primary_metric`:**
- `morans_i_final` — primary, equal to `morans_i`
- `morans_i_peak` — max over sampled trajectory (diagnostic; shows transient gap)
- `morans_i_sustained` — mean over last 20% of trajectory (diagnostic; stability check)
- `expected_i` — E[I] under random spatial null = −1/(N−1)
- `n_unique_types` — hard prerequisite: must be ≥ 2 (Moran's I is degenerate with 1 type)

**Secondary metrics:**
- Segregation index S(t) = ⟨s_i(t)⟩ at final step
- Cluster count and mean cluster size (connected components)
- `sustained_i_mean` and `sustained_i_cv` over last 20% of trajectory
  (stability check — used at confirmation tier)
- Mean path length (sorting models only)

**Null model:**
- *Label-shuffle on final state*: permute type labels across occupied
  positions on the terminal-state grid, recompute Moran's I, build null
  distribution (default 999 permutations). Matches the primary-metric
  substrate (final state) so observed-vs-null comparison is coherent.

**Detection tiers:**
- **Screening:** `n_unique_types ≥ 2` AND `I_final > expected_i` AND
  `I_final ≥ 0.05` (magnitude floor). The 0.05 floor was added in Sprint 10
  to reject pure-noise false positives; without it a random-label grid
  passes screening at p ≈ 0.03 through sampling variance alone (observed
  I ≈ 0.015, expected I ≈ −0.0002).
- **Confirmation:** Screening + null p < 0.01 + segregation index > 0.4 +
  `sustained_i_cv` < 0.3 (stability check — transient wavefronts with
  collapsing final states fail this gate even if final I happened to be
  above the floor).
- **Definitive:** Confirmation + p < 0.001 + P2/P3/P30 exclusions all cleared.

**Common false positives:**
- P2 (MIPS): type-blind clustering. Excluded via type-label preservation check.
- P3 (Turing): periodic spacing. Excluded via 2D FFT dominant-peak check.
- Transient wavefronts (SIR-style): **no longer a false positive** after
  Sprint 10; rejected at screening by the final-state primary.
- Initial-condition artifacts; finite-size trapping (require N ≥ 50).

**Nearest-neighbor exclusions:**
- P2: Type-blind density analysis — clustering persists without type distinction → P2
- P3: 2D FFT dominant peak → P3
- P30: Topological closure → P30

**Co-occurrence:**
- *Allowed:* P31 (if non-redundancy passes), P27 (spatial reciprocity uses
  P1-like clustering of cooperators — co-detection expected).
- *Screening-only co-occurrence:* RPS spatial (P12 definitive + P1 screening)
  — spiral domains produce sustained non-block clustering that screens but
  not confirms, because cluster geometry is cyclically-evolving rather than
  Schelling-style static segregation.
- *Excluded:* P2 (same spatial signal, different mechanism), P3 (periodic
  vs irregular), P30 (closure vs open clusters).

**Implementation note:** The legacy `morans_i_peak` field is preserved
in the primary metric dict for backward compatibility with tests that
historically inspected it directly (e.g. `test_sir_p22_e2e.py` asserts
peak > 0.5 and final < 0.1 on SIR as a physics diagnostic).

---

### P2 — Activity-induced phase separation (MIPS)

**Status:** Implemented Sprint 16. Canonical positive: Fily-Marchetti
Active Brownian Particles at phi=0.5, Pe=100, N>=800.

**Observable scope:** Model-metadata-assisted. Empirical detection runs
entirely on state_history (positions + velocities). Metadata flags
promote CONFIRMATION → DEFINITIVE via the mechanistic-null gate.

**Required raw observables:**
- `positions[t]`: (N, 2) float continuous in [0, L)² with periodic BC
- `velocities[t]`: (N, 2) float — OR `speeds[t]`: (N,) magnitudes
- (optional) `local_density[t]`: (N,) per-particle counts / area; if
  absent the detector computes it via cKDTree on positions.

**Required metadata for DEFINITIVE:**
- `box_size`: float — periodic box side length L
- `rho_star`: float — density threshold at which v(rho) → 0
  (fallback default 4.0 if absent, matches ABP with r_cg=1.0=sigma)
- `has_density_dependent_speed`: bool — must be True for DEFINITIVE
- `has_attraction_rule`: bool — must be False for DEFINITIVE
- `has_alignment_rule`: bool — must be False for DEFINITIVE

**Preprocessing:**
1. Compute local density rho_i(t) at each particle via cKDTree disk
   query with coarse-graining radius r_cg (default 1.0 = particle
   diameter). Self-inclusion in count; divide by pi*r_cg^2.
2. Compute per-particle speed magnitudes from velocities.
3. Collect rho and speed across burn-in-skipped measurement frames.

**Primary metric: two_phase_coexistence_score = min(f_gas, f_liquid)**

  f_gas = fraction of particles with rho_local < rho_star/2
  f_liquid = fraction of particles with rho_local > rho_star

True MIPS shows BOTH phases present in steady state. Flocking shows
f_liquid ≈ 1 with f_gas ≈ 0; uniform gas shows f_gas ≈ 1 with
f_liquid ≈ 0; density-saturated (stuck) shows f_liquid ≈ 1 with
f_gas ≈ 0. Only genuine coexistence produces a nontrivial minimum.

Score range [0, 0.5]. Empirical calibration at Fily-Marchetti canonical
regime (phi=0.5, Pe=100, N=1000, 3000 post-burn steps):
  score = 0.34 (f_gas=0.63, f_liquid=0.34).

**Secondary metrics:**
- `density_speed_r` (Pearson r of rho and |v|): genuine MIPS shows
  -0.95 < r < -0.30. r <= -0.99 is the dilute-Poisson artifact
  (few discrete density values give spurious perfect correlation).
- `cv_v` (coefficient of variation of |v|): > 0.30 discriminates
  genuine density-dependent slowdown from thermal/flocking regimes
  where speed is uniform.
- `p90_p10_ratio` (density dynamic range): > 5 indicates clustering;
  not itself a tier gate but a diagnostic.
- `frac_stalled` (fraction with |v| < 5% of mean): in (0.0, 0.98)
  for MIPS; = 1.0 for fully stuck.
- `mean_v`, `mean_rho`: sanity diagnostics.

**Null models:**

1. *Shuffle null (primary confirmation gate):* Independently permute
   rho_local[i] ↔ speed[i] pairings across particles, keeping each
   marginal distribution fixed. Under H0 (no density-velocity
   coupling), r ≈ 0. Observed r much more negative rejects H0.
   At canonical regime: null mean r ≈ 0.00, null std ≈ 0.013,
   observed r ≈ −0.94, Cohen's d ≈ −71 (effectively infinite).
   Requires n_permutations >= 199 for p < 0.01 floor.

2. *Mechanistic null (DEFINITIVE gate, metadata-based):* Model must
   affirm via metadata flags: has_density_dependent_speed=True,
   has_attraction_rule=False, has_alignment_rule=False. This is
   the P2 analogue of Sprint 13's continuous-field substrate gate
   for P3 (Decision 37) and Sprint 15's integer-velocity substrate
   gate for P8 (Decision 41). See ADR 43.

3. *Mechanistic intervention (orthogonal validation, NOT gated by
   detector):* Running the same model with rho_star → ∞ (constant
   speed, no density feedback) produces two_phase_score ≈ 0.0,
   f_liquid = 0.0. Phase 1f verified this on ABP canonical. The
   detector does not re-run the model; this validation is left to
   follow-on analysis.

**Detection tiers:**

- *Screening:* prereqs + two_phase_score > 0.03
- *Confirmation:* screening + two_phase_score > 0.08 +
  (-0.99 < density_speed_r < -0.30) + cv_v > 0.30 +
  frac_stalled < 0.98 + null_p < 0.01
- *Definitive:* confirmation + two_phase_score > 0.15 +
  mechanistic null metadata affirms (has_density_dependent_speed=True,
  has_attraction_rule=False, has_alignment_rule=False)

**Prerequisite gates (cause screening_rejection_reason):**
- `empty_state_history`: no snapshots
- `substrate_mismatch`: no positions, wrong shape, no velocities/speeds
- `too_few_particles`: N < 50
- `run_too_short`: post-burn < 300 snapshots
- `below_two_phase_floor`: primary <= 0.03

**Three false-positive traps identified in Phase 1 characterization:**

1. *Thermal regime* (low Pe, e.g. phi=0.5 Pe=5): v(rho) is active
   but the overall speed distribution is narrow because few
   particles reach rho_star. two_phase_score ~ 0.05, r ~ -0.95,
   but cv_v ~ 0.3 — fails confirmation gate.

2. *Dilute regime* (low phi, e.g. phi=0.1 Pe=100): local density
   takes few discrete values, f_liquid ≈ 0. r ≈ -1.0 (Poisson
   artifact from discrete counts), cv_v ≈ 0. Primary below floor
   rejects at screening.

3. *Over-saturated (stuck) regime* (high phi long runtime, e.g.
   phi=0.85 Pe=100 >= 3000 steps): system coarsens to one dense
   cluster, f_gas << 0.03. Short-runtime traces of this regime
   show transient two-phase coarsening that CAN false-positive to
   DEFINITIVE — detector requires >= 3·T_rot post-burn to avoid
   this. See ADR 46.

**Common false positives (across-model):**
- Flocking (Vicsek ordered): all particles in a moving band; f_liquid
  dominant but f_gas near zero → below screening. Additionally
  cv_v = 0 exactly (constant speed) blocks confirmation anyway.
- Attractive clustering (D'Orsogna milling): tight flock + empty
  surround; f_liquid ~ 0.7 but f_gas ~ 0.05 → SCREENING only.
  Additionally has_attraction_rule=True blocks DEFINITIVE via the
  mechanistic-null metadata gate.

**Nearest-neighbor exclusions:**
- P1 (similarity aggregation): excluded_by_substrate (lattice only).
- P6 (milling): excluded when has_attraction_rule=False;
  not_excluded when has_attraction_rule=True; inconclusive without
  metadata.

**Co-occurrence:**
- *Allowed:* P5 (flocking + v(rho) would co-exist in a Vicsek+v(rho)
  hybrid model — not currently in catalog).
- *Excluded:* P1 (different substrate), P6 (different mechanism).

**Empirical calibration (Sprint 16 Phase 1):**

  Canonical MIPS positive (N=1000, phi=0.5, Pe=100, rho_star=4.0,
  r_cg=1.0, dt=0.05, 2000 burn + 3000 measurement):

    two_phase_score: 0.25 – 0.40   (seed-dependent finite-size)
    f_gas: 0.35 – 0.65
    f_liquid: 0.30 – 0.65
    r: −0.87 to −0.95
    cv_v: 0.8 – 2.5 (highly regime-dependent)
    frac_stalled: 0.3 – 0.7
    null_p: < 0.005 (= 1/201 floor at n_permutations=199)
    Cohen's d: −50 to −100 (null std ~ 0.01, observed r ~ −0.9)
    → DEFINITIVE at all tested seeds.

  Negative sweep (continuous_2d neighbours):
    Vicsek ordered (eta=0.1):     score=0.017, f_gas=0.02
                                  → rejected (below_two_phase_floor)
    Vicsek disordered (eta=2.5):  score=0.002, f_gas=0.003
                                  → rejected
    D'Orsogna milling:            score=0.056, f_gas=0.056
                                  → SCREENING only
                                  (primary above floor but below
                                  confirmation; metadata flag would
                                  block DEFINITIVE)

  Within-model false-positive traps (ABP, varying regime):
    thermal (phi=0.5 Pe=5):       score=0.05, cv_v=0.3 → SCREENING
    dilute (phi=0.1 Pe=100):      score=0.00, f_liq=0 → rej_screening
    stuck long (phi=0.85 3500):   score=0.07, f_gas=0.07 → SCREENING
    stuck short (phi=0.85 2000):  score=0.24 → transient CONFIRMATION
                                  (documented caveat, ADR 46)

**References:**
- Fily, Y. & Marchetti, M. C. (2012). *Athermal Phase Separation of
  Self-Propelled Particles with No Alignment.* Phys. Rev. Lett.
  108, 235702.
- Redner, G. S., Hagan, M. F. & Baskaran, A. (2013). *Structure and
  Dynamics of a Phase-Separating Active Colloidal Fluid.* Phys.
  Rev. Lett. 110, 055701.
- Cates, M. E. & Tailleur, J. (2015). *Motility-Induced Phase
  Separation.* Annu. Rev. Condens. Matter Phys. 6, 219.

---

### P3 — Turing pattern formation

**Observable scope:** State-history only (continuous scalar field on a 2D periodic lattice).

**Canonical positive:** Gray-Scott reaction-diffusion (`GrayScott`; Pearson 1993). N=128, feed_rate=0.037, kill_rate=0.060 (labyrinth regime), Du=0.16, Dv=0.08, seed=42, ≥4000 timesteps.

**Substrate:** `lattice_2d_continuous` — new in Sprint 13. Distinguished from the integer-labelled `lattice_2d` substrate by the `n_unique_values` content prerequisite (see Prerequisites below). Gray-Scott is currently the only model on this substrate.

**Required raw observables:**
- `state_history`: list of snapshots each carrying a `field` observable (2D float array). The detector intentionally does NOT fall back to `grid` — consuming an integer grid would be a silent substrate violation.
- Grid dimensions (periodic boundaries).

**Preprocessing:**
1. Extract the final `field` snapshot (v-field for Gray-Scott).
2. Subtract spatial mean (kills DC bin).
3. 2D FFT, azimuthal average → radial power spectrum `radial[k]` for k = 0, 1, ..., N/2 − 1.
4. Search for the peak in radial[k_min : k_max] with k_min = 2 (skip DC and first bin) and k_max = N/2 (full radial range, Sprint 13 characterization showed this gives the cleanest signal-to-noise).

**Prerequisites (all three must pass):**
- **n_unique_values ≥ 50** in the final field. The load-bearing Sprint 13 discriminator. RPS at mobility = 1e-4 produces a raw integer-grid radial-FFT peak-to-mean ≈ 23 — *numerically matching* Gray-Scott labyrinth. Discrimination is substrate-level (GS has ~16k unique float values; every integer-grid model has ≤ 10), not empirical threshold tuning. Rejects RPS, GH, Schelling, Nowak-May, GoL, SIR, LV cleanly without any parameter fitting.
- **field_std ≥ 0.01**. Rejects non-Turing Gray-Scott parameter regimes (uniform decay at F=0.10, k=0.10 gives v → 0; uniform high at F=0.01, k=0.02 gives constant v ≈ 0.30).
- **2 ≤ peak_k ≤ N/4**. Rejects DC/first-bin drift and grid-discretization artifacts near Nyquist.

**Primary metric: radial-FFT peak-to-mean**

$$\text{peak\_to\_mean} = \max_{k_\text{min} \le k < k_\text{max}} \text{radial}[k] \;\big/\; \text{mean}(\text{radial}[k_\text{min} : k_\text{max}])$$

**Secondary metrics:**
- `peak_k`, `wavelength_pixels = N / peak_k`
- `peak_k_cv`: std/mean of peak_k across the last N snapshots (default last 5 at stride 50) — stationarity check for the selected wavelength
- `peak_to_mean_min`: minimum peak-to-mean across late snapshots (weakness check)

**Null model: spatial shuffle (SHUFFLE)**

Shuffle the field values while preserving the marginal intensity distribution. Destroys all spatial structure including the Turing wavelength. This is a clean null for Turing-pattern detection — for Gray-Scott labyrinth the null peak-to-mean distribution is centered at ≈ 1.42 ± 0.16, giving Cohen's d ≈ 107 against the observed ≈ 18.75.

**Detection tiers:**
- *Screening:* prerequisites pass AND peak_to_mean > 5.0 AND peak_k ∈ [2, N/4].
- *Confirmation:* screening AND peak_to_mean > 10.0 AND Cohen's d > 10.0 AND peak_k_cv < 0.15 AND null p < 0.01.
- *Definitive:* confirmation AND peak_to_mean > 15.0 AND Cohen's d > 30.0 AND peak_k_cv < 0.05.

Requires `n_permutations ≥ 199` for confirmation (achieves p ≤ 0.005 at the floor). Default is 199.

**Canonical-positive measurements (Sprint 13):**

| Seed | N | T | peak_k | λ (px) | peak/mean | Cohen's d | null p | peak_k_cv | tier |
|------|---|---|--------|--------|-----------|-----------|--------|-----------|------|
| 42   | 128 | 4000 | 10 | 12.8 | 18.75 | 102.8 | 0.005 | 0.000 | DEFINITIVE |
| 7    | 128 | 4000 | 10 | 12.8 |  ≈18  |  ≈100 | 0.005 | 0.000 | DEFINITIVE |
| 123  | 128 | 4000 | 10 | 12.8 |  ≈18  |  ≈100 | 0.005 | 0.000 | DEFINITIVE |
| 42   |  64 | 3000 |  5 | 12.8 | 16.98 |  >70  | 0.005 | 0.000 | DEFINITIVE |
| 42   |  96 | 3000 |  7 | 13.7 | 22.84 |  >70  | 0.005 | 0.000 | DEFINITIVE |

Spots regime (F=0.030, k=0.062) at N=128, T=4000, seed=42: peak_k=11, λ=11.6 px, peak/mean=13.35, d=58.8 → **CONFIRMATION**.

**Wavelength invariance (Sprint 13 scaling test):**
peak_k scales linearly with N; wavelength in pixels is ≈ 12 and invariant. This is a physical property of the reaction-diffusion system (Du, Dv, F, k), not a grid artifact.

**Common false positives and how P3 rejects them:**
- **RPS (low mobility, raw grid)**: peak/mean ≈ 23 — higher than Gray-Scott labyrinth. Rejected at substrate prerequisite (no `field` observable on the state snapshots); also caught at the n_unique_values prerequisite if someone manually injects the grid as a `field`.
- **GH excitable spirals**: peak/mean ≈ 6.6 on raw grid. Rejected at substrate prerequisite.
- **Schelling segregation**: peak/mean ≈ 3.3 on raw grid. Rejected at substrate prerequisite (also fails peak/mean screening).
- **Non-Turing Gray-Scott regimes**: v-field decays to 0 or a spatially uniform value. Rejected at field_std prerequisite.
- **"Pearson spots" at small N**: the literature-standard F=0.062, k=0.0609 gives wavelength ≈ 64 pixels at N=128, which is a domain-size artifact (the pattern needs N ≥ 256 to resolve short-wavelength structure). Documented in `GrayScott` docstring; use F=0.030, k=0.062 as the short-wavelength "spots" alternative for N=128.

**Nearest-neighbor exclusions (Sprint 13 implementation):**
- **P1 (aggregation)**: canonical P1 positives (Schelling, Zhang sorting) live on integer-grid substrates; the `n_unique_values ≥ 50` prerequisite auto-clears P1. Marked `excluded` when P3 fires.
- **P13 (excitable waves)**: canonical P13 positive (GH) is integer-grid. Same auto-clear; marked `excluded`.

**Co-occurrence:**
- *Allowed:* none registered at Sprint 13 (Gray-Scott is the only model on `lattice_2d_continuous`).
- *Excluded:* P1, P13 (substrate-level).

**Implementation:** `epc.detectors.p3_turing_wavelength.P3TuringWavelengthDetector`. Metric module: `epc.metrics.turing_wavelength`.

**References:**
- Turing, A. M. (1952). "The Chemical Basis of Morphogenesis." Phil. Trans. R. Soc. Lond. B 237, 37-72.
- Gray, P. & Scott, S. K. (1983). "Autocatalytic reactions in the isothermal, continuous stirred tank reactor." Chem. Eng. Sci. 38, 29-43.
- Pearson, J. E. (1993). "Complex Patterns in a Simple System." Science 261, 189-192.

---

### P4 — Territoriality / exclusion boundaries

**Observable scope:** State-history only (scent field is part of state)

**Required raw observables:**
- `positions[t]`: (N, 2) per timestep
- Agent identity labels: (N,)
- `scent[t]`: (rows, cols, N) intensity per agent per cell

**Preprocessing:**
1. Home range per agent: KDE at 95% isopleth.
2. Pairwise overlap: O_{ij} = area(HR_i ∩ HR_j) / area(HR_i ∪ HR_j).
3. Boundary cells: where ≥ 2 agents' scent within 20% of each other.

**Primary metric: Home range overlap**
Mean ⟨O_{ij}⟩ trending to near-zero.

**Secondary metrics:**
- Boundary persistence over sliding window of 0.1 × T_cross
- Spatial cross-correlation (occupancy vs foreign scent)
- Territory size CV

**Null models:**
- *No-scent:* random walkers, no marking. High overlap.
- *Scent-no-avoidance:* mark but don't avoid. Tests avoidance necessity.

**Detection tiers:**
- *Screening:* ⟨O⟩ < 0.2
- *Confirmation:* ⟨O⟩ < 0.1 AND boundary persistence > 80% over last 50% AND below no-scent null (p < 0.01)
- *Definitive:* Confirmation + scent-no-avoidance null shows high overlap + P1/P29 exclusions

**Common false positives:**
- P1 (clustering, not exclusion). Hard-core crowding. Transient separation.

**Nearest-neighbor exclusions:**
- P1: Attraction → P1. Exclusion with boundary → P4
- P29: Shared trails → P29. Exclusive territories → P4

**Co-occurrence:**
- *Allowed:* P5 (territorial flockers), P32 (territorial specialists)
- *Excluded:* P1 (attraction vs exclusion), P29 (shared vs exclusive)

---

## Cluster B: Collective Motion

### P5 — Translational alignment (flocking)

**Observable scope:** State-history only

**Required raw observables:**
- `positions[t]`: (N, 2)
- `velocities[t]`: (N, 2)

**Preprocessing:**
1. Unit headings: v̂_i = v_i/|v_i|. Exclude zero-velocity agents.
2. COM velocity: V_cm = (1/N) Σ v_i.
3. Unwrap for periodic boundaries.

**Primary metric: Polarization order parameter**
φ(t) = |(1/N) Σ v̂_i(t)|

**Secondary metrics:**
- Group speed ratio R = |V_cm| / ⟨|v_i|⟩
- Angular momentum L = (1/N) Σ (r̂_i × v̂_i)
- Directional persistence (heading autocorrelation)

**Null models:**
- *Heading-shuffle:* permute directions, keep speeds. φ ≈ 1/√N.
- *Zero-coupling:* alignment strength = 0. Independent walks.

**Detection tiers:**
- *Screening:* φ_mean > 0.5 over ≥ 5 × T_cross
- *Confirmation:* φ_mean > 0.7 over ≥ 10 × T_cross AND R > 0.5 AND above shuffle null (p < 0.01)
- *Definitive:* Confirmation + P6/P7/P8 exclusions + zero-coupling null shows φ ≈ 1/√N

**Common false positives:**
- P6 (milling): high |L|, low R. P7 (lanes): bimodal antiparallel headings. External field. Transient alignment.

**Nearest-neighbor exclusions:**
- P6: |L| > 0.5 → milling
- P7: Bimodal antiparallel heading distribution → lanes
- P8: Backward-propagating density waves → jamming

**Co-occurrence:**
- *Allowed:* P17 (flocking + gradient sensing), P19 (flocking + leadership), P9 (flocking + synchronization)
- *Excluded:* P6 (translational vs rotational), P7 (uniform vs counterflow)

---

### P6 — Milling / vortex formation

**Observable scope:** State-history only

**Required raw observables:**
- `positions[t]`: (N, 2)
- `velocities[t]`: (N, 2)

**Preprocessing:**
1. COM: X_cm = (1/N) Σ r_i. Radial vectors: r̂_i = (r_i−X_cm)/|r_i−X_cm|.
2. Radial density profile ρ(r).

**Primary metric: Angular momentum**
L(t) = (1/N) Σ (r̂_i × v̂_i). |L| near 1 = milling.

**Secondary metrics:**
- Group speed ratio R < 0.3
- Ring density profile: hollowness ρ(0)/ρ(r*) < 0.5
- Rotation coherence (angular velocity variance)

**Null models:**
- *Heading-shuffle:* |L| ≈ 0.
- *Zero-coupling:* independent walkers.

**Detection tiers:**
- *Screening:* |L| > 0.3 over ≥ 5 × T_cross
- *Confirmation:* |L| > 0.5 AND R < 0.3 over ≥ 10 × T_cross
- *Definitive:* Confirmation + ring density confirmed + P5 exclusion

**Common false positives:**
- P5 (high R). Transient vortices. Boundary-confined circulation.

**Nearest-neighbor exclusions:**
- P5: R > 0.5 → flocking
- P12: Species spirals (nontransitive competition, not physical rotation)

**Co-occurrence:**
- *Allowed:* P9 (synchronized milling)
- *Excluded:* P5 (translational vs rotational)

---

### P7 — Lane formation in counterflow

**Observable scope:** State-history only (goal directions inferable from persistent heading)

**Required raw observables:**
- `positions[t]`: (N, 2)
- `velocities[t]`: (N, 2)
- `goal_dir`: (N,) target directions

**Preprocessing:**
1. Classify by goal direction. 2. Local directional segregation per agent. 3. Head-on collision count.

**Primary metric: Lane order parameter**
Ψ = ⟨same-direction neighbor fraction⟩ − f_random

**Secondary metrics:**
- Collision rate reduction vs mixed baseline
- Throughput gain
- Lane count

**Null models:**
- *Random-direction:* no counterflow structure.
- *No-avoidance:* counterflow but no steering.

**Detection tiers:**
- *Screening:* Ψ significantly above random (p < 0.05)
- *Confirmation:* Ψ significant (p < 0.01) AND collision rate < 50% of mixed baseline
- *Definitive:* Confirmation + domain width > 5× lane width + P5/P29 exclusions

**Common false positives:**
- P5 (unimodal heading). P29 (persistent infrastructure). Narrow corridors.

**Nearest-neighbor exclusions:**
- P5: Unimodal heading → P5. Bimodal antiparallel → P7.
- P29: Vanish when agents stop → P7. Persist → P29.

**Co-occurrence:**
- *Allowed:* P8 (lanes can jam)
- *Excluded:* P5 (uniform vs counterflow)

---

### P8 — Self-organized jamming / stop-go waves

**Observable scope:** State-history only

**Canonical positive model:** Nagel-Schreckenberg traffic CA (1D ring,
`lattice_1d` substrate), `L=1000`, `v_max=5`, `p_slow=0.3`, `ρ=0.15`.

**Required raw observables:**
- `velocities[t]`: (n_cars,) integer array with values in [0, v_max].
  Provided by Nagel-Schreckenberg; not emitted by any other catalog
  model (substrate-level discrimination).
- `gaps[t]`: (n_cars,) integer array. Used for secondary metrics.

**Preprocessing:**
1. Collect per-step velocity history across post-burn-in snapshots into a
   `(T, n_cars)` integer matrix.
2. Compute marginal stopped-fraction (primary).
3. For each car, extract consecutive-`v=0` run lengths → pooled lifetime
   distribution (secondary).

**Primary metric: Stopped fraction `P(v=0)`**
`stopped_fraction = <1[v_i(t) = 0]>_{t > burn_in, i}` — time-and-car-
averaged fraction of v=0 states. This is the Bette-Habel-Emig-
Schreckenberg (2017) order parameter for the NS jamming transition.

Choice rationale (Sprint 15 empirical finding): spacetime-FFT primaries
(the Sprint 6-era design) work but require a 2D spacetime diagram and
are more expensive; `stopped_fraction` has tiny seed-to-seed variance
(σ/mean < 2% at L=1000) and matches the published order parameter.

**Secondary metrics:**
- `jam_lifetime_p95`, `jam_lifetime_max`: 95th-percentile and maximum of
  per-car consecutive-v=0 run lengths. Heavy-tailed in true jamming
  (p95 ≈ 13, max > 50 at canonical positive); short and bounded in
  pigeonhole density saturation (p95 ≤ 4 even at ρ=0.80).
- `gap_cv`, `zero_gap_fraction`: coefficient of variation and zero-gap
  fraction of the gap distribution. Bimodal gap distribution in jammed
  regimes (bumper-to-bumper inside jams, large gaps between).
- `n_jam_events`: total number of stopped-runs across all cars.

**Null model:** Per-car independent temporal shuffle of `v(t)`, applied
to the jam_lifetime_p95 statistic. Preserves marginal v-distribution
per car (so stopped_fraction is unchanged under the null), destroys the
persistence structure that makes real NS lifetimes heavy-tailed. Under
the null, `p95` collapses to the geometric-distribution p95 at the
observed stopped marginal (≈ 2 steps), giving effectively infinite
Cohen's d at canonical parameters (null std ≈ 0 across 199 permutations).

n_permutations = 199 required for confirmation tier (p < 0.01 floor
= 1/200 = 0.005).

**Detection tiers (Sprint 15 empirical calibration):**
- *Screening:* substrate prereqs PASSED AND `stopped_fraction > 0.05`.
- *Confirmation:* screening + `jam_lifetime_p95 > 5` + null `p < 0.01`.
  The `p95` gate is the designed discriminator between emergent NS
  jamming (p95 > 10) and pigeonhole density saturation (p95 ≤ 4 at
  ρ=0.80, p=0).
- *Definitive:* confirmation + `stopped_fraction > 0.15` +
  `jam_lifetime_max > 20`.

**Substrate-level prerequisites (hard rejections):**
- `velocities` observable must be present as a 1D integer array
  (rejects all lattice_2d, continuous_2d, oscillator, opinion_space
  models; rejects Zhang sorting which is lattice_1d but has no
  velocities observable).
- Velocity values must lie in [0, 64] (rejects Vicsek/D'Orsogna's
  continuous 2D velocity vectors should substrate gating ever be
  bypassed).
- `n_cars >= 20` and post-burn-in run length `>= 100`.

**Common false positives (and how they are rejected):**
- **Pigeonhole density saturation** (NS at ρ > 1/(v_max+1) with p=0):
  stopped_fraction high, but `jam_lifetime_p95 ≤ 4`. Rejected at
  confirmation.
- **P14 (SOC scale-free avalanches):** 2D substrate, incompatible with
  `velocities` observable. Rejected at substrate registration.
- **Periodic forcing / external bottleneck:** not present in the NS
  canonical model family; would need a separate model variant to test.

**Nearest-neighbor exclusions:**
- None currently active. P8's substrate prerequisite (1D integer
  velocities) cleanly excludes all other EPC patterns by construction.
  Future pattern additions on lattice_1d would trigger an exclusion
  check update (e.g., a hypothetical P4 on 1D occupancy, or a P13
  variant on 1D excitable media).

**Co-occurrence:**
- None currently. NS is the only lattice_1d traffic model in the
  catalog. Future multi-lane or counterflow variants could allow P7
  co-occurrence.

**Sprint 15 calibration anchors (L=1000, v_max=5, p_slow=0.3):**

| Regime | ρ | stopped | lt_p95 | lt_max | Tier |
|---|---|---|---|---|---|
| Free flow | 0.05 | 0.000 | 0 | 0 | Screening (rej) |
| Free flow | 0.08 | 0.000 | 1 | ~1 | Screening (rej) |
| Onset | 0.10 | 0.003 | 0 | 0 | Screening (rej) |
| **Near-transition** | **0.12** | **0.082** | **12** | ~50 | **CONFIRMATION** |
| **Canonical positive** | **0.15** | **0.181** | **13** | ~60 | **DEFINITIVE** |
| Deep jam | 0.30 | 0.431 | 13 | ~85 | DEFINITIVE |
| Deterministic (p=0) | 0.15 | 0.000 | 0 | 0 | Screening (rej) |
| Density saturation (p=0) | 0.80 | 0.750 | 4 | 6 | Screening (confirmation rej) |

---

## Cluster C: Temporal Dynamics and Synchronization

### P9 — Temporal synchronization

**Observable scope:** State-history only (K_c sweep requires model access; primary detection does not)

**Required raw observables:**
- `θ[t]`: (N,) phases in [0, 2π)
- Coupling topology, coupling strength K

**Preprocessing:**
1. Instantaneous phase (Hilbert transform if needed). 2. Local r for subregions.

**Primary metric: Kuramoto order parameter**
r(t) = |(1/N) Σ exp(iθ_i)|

**Secondary metrics:**
- Frequency entrainment (σ of instantaneous frequencies → 0)
- Order parameter fluctuation variance

**Null models:**
- *Uncoupled (K=0):* r ≈ 1/√N.
- *Random-phase:* uniform [0, 2π).

**Detection tiers:**
- *Screening:* r > 0.7 over ≥ 10 T_osc
- *Confirmation:* r > 0.9 over ≥ 50 T_osc AND above uncoupled null (p < 0.01)
- *Definitive:* Confirmation + ≥ 100 T_osc sustained + P10 exclusion (uniform local r)

**Common false positives:**
- P10 (partial sync). External pacemaker. Transient coherence.

**Nearest-neighbor exclusions:**
- P10: Spatially heterogeneous local r → P10. Uniform → P9.

**Co-occurrence:**
- *Allowed:* P5 (synchronized flockers), P20 (quorum triggers sync)
- *Excluded:* P10 (full vs partial synchronization)

---

### P10 — Chimera states

**Observable scope:** Model-metadata assisted (coupling kernel shape useful but not required)

**Required raw observables:**
- `θ[t]`: (N,) with spatial ordering
- Coupling kernel specification

**Preprocessing:**
1. Partition into W spatial windows. 2. Local r_w per window. 3. Spatial r_w profile.

**Primary metric: Coherent/incoherent coexistence**
Bimodal distribution of {r_w}: some > 0.9, others < 0.5.

**Secondary metrics:**
- Spatial persistence of partition (temporal autocorrelation of r_w profile)
- Global r in (0.3, 0.8)

**Null models:**
- *Global-coupling:* all-to-all → full sync, no chimera.
- *Short-range:* nearest-neighbor only → no chimera.

**Detection tiers:**
- *Screening:* Dip test on {r_w} p < 0.05 AND global r < 0.9
- *Confirmation:* Dip p < 0.01 AND persists ≥ 50 T_osc AND global r in (0.3, 0.8)
- *Definitive:* Confirmation + ≥ 100 T_osc + global-coupling null eliminates chimera + P9 exclusion

**Common false positives:**
- P9 (uniform high r). Transient disorder. P21 (opinion fragmentation — different substrate).

**Nearest-neighbor exclusions:**
- P9: Uniform r_w → P9. Heterogeneous → P10.

**Co-occurrence:**
- *Allowed:* (chimera is typically exclusive within oscillator systems)
- *Excluded:* P9 (full vs partial sync)

---

### P11 — Predator-prey oscillations

**Observable scope:** State-history only (species labels on a grid, plus an empty-reservoir state)

**Canonical positive:** Lotka-Volterra lattice (`LotkaVolterraLattice`; Mobilia-Georgiev-Täuber 2007). L=100, predation_rate=4.0, prey_reproduction_rate=1.0, predator_death_rate=1.0, seed=42, ≥1200 generations.

**Required raw observables:**
- `state_history`: list of snapshots each with `grid` (2D int array; state 0 = empty, states 1+ = species)
- Exactly two non-empty species observed

**Preprocessing:**
1. Auto-detect the two most-populous non-zero grid states (or use explicit `species_state_hint` from metadata).
2. Extract fraction-vs-time for each species.
3. Drop first `burn_in` generations (default 100).

**Prerequisites (all three must pass — see Decision 35):**
- **n_unique_species_observed == 2**. Rejects RPS (3 species) even though RPS species pairs give stronger anti-correlation than LV on the primary metric.
- **std(species_A) > 0.005 AND std(species_B) > 0.005**. Rejects Schelling (per-agent identity is conserved, so species fractions never change).
- **std(species_A + species_B) > 0.005**. Rejects strictly-conserved 2-species systems like Nowak-May (coop + defect = 1 exactly → rho_anti = -1.0 by algebra, a conservation artifact). LV has a nontrivial empty reservoir: std(prey + pred) ≈ 0.03.

**Primary metric: rho_anti at nonzero lag**

$$\rho_{\text{anti}} = \min_{|\tau| \ge 5} \text{Pearson}\big(\text{prey}(t), \text{pred}(t + \tau)\big)$$

Rationale (Decision 34): the LV quarter-period phase shift plus strong anti-phase coupling means the minimum of the cross-correlogram is at |τ| ~ 10-20 rather than at the prettier τ = T/4 ~ 40. Empirically measured |tau_anti| ∈ [11, 18] for canonical LV. The |τ| ≥ 5 exclusion band rejects near-instantaneous anti-correlation (which would be a conservation artifact in systems like Nowak-May).

**Secondary metrics:**
- `fft_peak_to_mean`: predator density FFT peak-to-mean ratio (skip DC; period ≥ 3). LV empirical range: 18–32 across seeds.
- `|tau_anti|`: lag at which rho_anti is achieved.
- `rho_at_zero_lag`: contrast with lagged rho_anti.
- `rho_anti` on first half / second half with `halves_agree` flag (robustness).

**Null model: circular shift (SHUFFLE)**

At each permutation, circularly shift the predator time series by a uniformly random offset in [1, n-1]. Preserves each series' marginal distribution, autocorrelation, and FFT magnitude spectrum; destroys the cross-series phase relationship.

**Decision 36**: p-value is NOT a tier gate. The circular-shift null preserves the autocorrelation that makes LV's extreme rho_anti possible, so on LV the null frequently produces its own extreme anti-correlations (p ≈ 0.05-0.15 even when Cohen's d ≈ -2). Signal-vs-noise separation is achieved by rho_anti magnitude (|LV| ~ 0.8, |noise| ~ 0.1) and cohens_d. The null is retained for reporting the effect-size diagnostic.

**Detection tiers:**
- *Screening:* prerequisites pass AND rho_anti < −0.3 AND |tau_anti| ≥ 5 AND fft_peak_to_mean > 6.
- *Confirmation:* screening AND rho_anti ≤ −0.5 AND cohens_d ≤ −1.5 AND halves_agree.
- *Definitive:* confirmation AND rho_anti ≤ −0.7 AND fft_peak_to_mean > 12 AND cohens_d ≤ −1.5.

**Canonical-positive measurements (Sprint 11):**

| Seed | n_steps | rho_anti | tau_anti | fft_p2m | cohens_d | tier |
|------|---------|----------|----------|---------|----------|------|
| 42   | 1500    | −0.863   | −16      | 25.1    | −2.21    | DEFINITIVE |
| 42   | 1200    | −0.863   | −15      | 21.6    | −1.75    | DEFINITIVE |
| 7    | 1200    | −0.716   | −16      | 13.6    | −1.77    | DEFINITIVE |
| 123  | 1200    | −0.780   | −14      | 17.6    | −1.95    | DEFINITIVE |

**Common false positives and how P11 rejects them:**
- **RPS (3-species cyclic)**: rho_anti measured at −0.93 to −0.97 on species pairs, *stronger* than LV. Rejected via n_species prerequisite.
- **Nowak-May (strict 2-species conservation)**: rho_anti = −0.979 at lag +3 by algebraic conservation. Rejected via total_std prerequisite (= 0.000).
- **Schelling (static identity)**: species fractions never change. Rejected via species_std prerequisite.
- **SIR (transient wavefront)**: active phase is O(L) generations only; post-burn-in variance ~ 0. Rejected via species_std prerequisite.
- **White noise**: rho_anti ≈ −0.08 << −0.3. Rejected at screening.

**Nearest-neighbor exclusions (Sprint 11 implementation):**
- **P12**: requires ≥ 3 cyclically-competing species. P11's n_species == 2 prerequisite auto-clears P12; when fired, P11 marks P12 as `excluded`.
- **P9**: intra-species synchronization (phase coherence among identical oscillators). P11 is inter-species anti-phase; marked `excluded` when P11 fires with n_species == 2.

**Co-occurrence:**
- *Allowed:* P1 (spatial clustering of predator/prey populations). Measured: LV × P1 = CONFIRMATION with I_final ≈ 0.46, seg ≈ 0.70 at standard n_perm = 999.
- *Excluded:* P12 (bilateral vs nontransitive).

**Implementation:** `epc.detectors.p11_predator_prey_oscillation.P11PredatorPreyDetector`. Metric module: `epc.metrics.predator_prey_crosscorr`.

**References:**
- Mobilia, M., Georgiev, I.T. & Täuber, U.C. (2007). "Phase Transitions and Spatio-Temporal Fluctuations in Stochastic Lattice Lotka-Volterra Models." J. Stat. Phys. 128, 447-483. [arXiv: q-bio/0512039]
- Heiba, B., Chen, S. & Täuber, U.C. (2018). "Boundary effects on population dynamics in stochastic lattice Lotka-Volterra models." Physica A 491, 582-590. [arXiv: 1706.02567]
- Täuber, U.C. (2024). "Stochastic spatial Lotka-Volterra predator-prey models." arXiv: 2405.05006.

---

### P12 — Cyclic dominance / coexistence waves

**Observable scope:** State-history only

**Substrate:** lattice_2d (implemented). Generalizes in principle to any
discrete-state substrate with local neighborhoods (e.g., graph-based
RPS), but the current implementation assumes a 2D grid.

**Required raw observables:**
- `grid[t]`: (rows, cols) int array with discrete states. At least three
  distinct non-empty species must be observed across the trajectory. A
  state labelled `0` is treated as empty if the trajectory contains ≥ 4
  distinct states; otherwise all observed states are candidate species.
- `grid_dims`: (rows, cols)

**Preprocessing:**
1. Drop transient: analyze only the second half of the trajectory.
2. Identify candidate species set (observed non-empty states).
3. For each ordered pair `(X, Y)` with `X ≠ Y`, accumulate four counts
   across `(t, t+1)` pairs:
   - `n_X_with_Y_nbr` — cells in state X with ≥ 1 Y-neighbor at t
   - `n_X_without_Y_nbr` — cells in state X with no Y-neighbor at t
   - `n_became_Y_with` — subset of `n_X_with_Y_nbr` that transitioned to Y
   - `n_became_Y_without` — subset of `n_X_without_Y_nbr` that transitioned to Y
4. Compute the **neighbor-conditional replacement ratio**:
   ```
   ρ(X, Y) = (p_with + ε) / (p_without + ε)
   ```
   where `p_with = n_became_Y_with / n_X_with_Y_nbr`,
   `p_without = n_became_Y_without / n_X_without_Y_nbr`, and `ε = 1e-6`
   is a Laplace-smoothing constant.
5. For every 3-subset of candidate species and every cyclic ordering
   `(a, b, c)` (meaning "b replaces a, c replaces b, a replaces c"),
   compute `min(ρ(a,b), ρ(b,c), ρ(c,a))`. Pick the ordering that
   maximizes this minimum.

**Primary metric: log₁₀ of min forward-cycle ρ over best triple**

```
intransitivity_score = log10(max over triples (min forward ρ))
```

**Why this separates RPS from GH.** In RPS, transitions are entirely
neighbor-dependent: a cell changes state only when a specific-species
neighbor is selected by the Gillespie scheduler. For dominance edges,
`p_with ≫ p_without`, giving ρ ≈ 70–200 on canonical coexistence
parameters. In Greenberg-Hastings, transitions `k → (k+1) mod n` for
`k > 0` are purely clock-driven and do NOT depend on neighbors; the
only neighbor-dependent transition is `0 → 1`. Because the clock forces
state-k cells to transition whether or not they have a (k+1)-neighbor,
`p_with = p_without = 1` exactly, so `ρ = 1.0` for GH. This gives a
two-orders-of-magnitude separation between the two mechanisms.

**Secondary metrics:**
- `coexistence_fraction` — fraction of second-half snapshots where all
  three triple-members are present at ≥ 5% density
- `min_species_density` — minimum over time of min species fraction
- `direction_stable` — whether the identified cyclic triple from the
  third and fourth quarters of the trajectory agree up to rotation
- `identified_triple` — the cyclic direction, e.g. `[1, 3, 2]` meaning
  species 3 replaces 1, species 2 replaces 3, species 1 replaces 2

**Null model:**

*Spatial shuffle:* at each timestep in the analyzed window, globally
permute the cells of the grid. This preserves species marginals
(the proportion of each state) but destroys the neighbor-transition
correlation. Under this null, `p_with ≈ p_without` for every edge,
so `ρ → 1` and `intransitivity_score → 0`. `n_permutations = 199`
(floor p = 0.005) is the default; use 499+ to reach `p < 0.005`
needed for definitive tier.

For the canonical positive (RPS at M = 10⁻⁴, L = 30), null scores
cluster at `0.01 ± 0.01` against an observed score of 1.8, giving
Cohen's d > 200.

**Detection tiers** (implemented in `P12CyclicDominanceDetector`):

- *Screening:* `n_candidate_species ≥ 3` AND
  `intransitivity_score > 1.0` (i.e., min ρ > 10)
- *Confirmation:* Screening + `intransitivity_score > 1.3` (min ρ > 20)
  AND `coexistence_fraction ≥ 0.8` AND null `p < 0.01`
- *Definitive:* Confirmation + `intransitivity_score > 2.0`
  (min ρ > 100) AND `direction_stable = True` AND null `p < 0.005`

**Common false positives:**

- P13 (excitable waves): LOOKS similar — persistent spiral-like
  wavefronts on a lattice with ≥ 3 states — but clock-driven; produces
  ρ ≈ 1.0 which fails P12 screening.
- P11 (predator-prey oscillations): bilateral (2 species); fails the
  3-candidate-species prerequisite.
- P3 (Turing patterns): no neighbor-dependent species replacement;
  produces low ρ.

**Nearest-neighbor exclusions:**

- P13: RPS dominance edges yield ρ ≈ 70+. If `best_rho_min > 10`, mark
  P13 as `excluded` (rules out clock-driven mechanism). If
  `best_rho_min < 2`, mark P13 as `not_excluded`. Otherwise
  `inconclusive`.
- P22: RPS is re-entrant (species cycle indefinitely); P22 cascades
  are single-pass. If model metadata reports `model_class =
  "cyclic_competition"` or `model_name = "rps_spatial"`, mark P22
  as `excluded`. If metadata reports an epidemic-class model, mark
  P22 as `not_excluded`.

**Co-occurrence:**

- *Allowed:* P1 (cyclic spirals produce transient spatial clustering;
  RPS × P1 fires at screening with Moran's I ≈ 0.57 during wavefront
  passage, analogous to SIR × P1).
- *Excluded:* P13 (mechanistic alternative for the same visual
  signature), P22 (re-entrant vs. single-pass is definitional).

**Verified on:** RPS (Reichenbach 2007, M = 10⁻⁴, L = 30, 80 generations,
n_permutations = 199) × P12 → CONFIRMATION (score = 1.83, ρ_min = 67.6,
p = 0.005, Cohen's d > 200, coexistence = 1.0, direction stable, P13
and P22 both excluded). Reaches DEFINITIVE at M = 10⁻⁵, 150 generations,
n_permutations = 499 (score = 2.06, p = 0.002, confidence = 0.85). See
`tests/test_rps_p12_e2e.py::TestRPSDetectedByP12` and
`tests/test_rps_replication.py` for the underlying model validation.

**Verified rejections:** Greenberg-Hastings (n=3 and n=5, `ρ = 1.0`,
score = 0.0), SIR (no cyclic dominance, score = 0.0), Game of Life
(n_candidate_species = 2, fails prerequisite), Nowak-May (same
prerequisite failure). See `tests/test_rps_p12_e2e.py::TestP12CrossDetectionRejections`.

**P13 boundary note:** The Sprint 9 investigation documented that P13
already rejects RPS cleanly at screening via wavefront-speed CV
(CV = 0.59–0.68 on RPS vs. CV ≈ 0.05–0.15 on GH). This is an INDEPENDENT
discrimination mechanism from P12's neighbor-conditional ratio — P13
distinguishes clock-driven from stochastic-neighbor-driven wavefront
propagation via speed uniformity, while P12 distinguishes them via
per-cell transition logic. Both mechanisms agree on every tested case.
See `tests/test_rps_p13_boundary.py` for the pinned behavior.

---

## Cluster D: Wave Propagation and Excitable Dynamics

### P13 — Excitable spiral and target waves

**Observable scope:** State-history only

**Required raw observables:**
- `grid[t]`: (rows, cols) with discrete states (rest/excited/refractory)
- State transition rules

**Preprocessing:**
1. Wavefronts: rest → excited transitions. 2. Spiral tips via topological charge. 3. Propagation speed.

**Primary metric: Wavefront speed and spiral persistence**
Constant-speed propagation (CV < 0.15) with persistent wave sources.

**Secondary metrics:**
- Topological charge (net 0 in periodic BC)
- Refractory tail length consistency
- Wave source count

**Null models:**
- *No-refractory:* immediate re-excitation → uniform flash.
- *High-threshold:* no propagation → quiescent.

**Detection tiers:**
- *Screening:* Persistent wavefront for ≥ 5 × T_prop with speed CV < 0.2
- *Confirmation:* Speed CV < 0.15 AND (spiral tip ≥ 50 rotations OR target source ≥ 50 cycles)
- *Definitive:* Confirmation + ≥ 100 rotations/cycles + P15 discrimination via TE AND functional test

**Common false positives:**
- P15 (computation): propagating structures that compute. P12 (species spirals). Periodic forcing.

**Nearest-neighbor exclusions (revised P13/P15 boundary):**
- P15: Two-stage discriminator required:
  1. *TE test:* Compute Transfer Entropy across wave/structure collisions. TE ≈ 0 supports P13 (annihilative collisions). TE > 0 is necessary but not sufficient for P15.
  2. *Functional test:* Do collision outcomes vary systematically with input configuration? If collisions produce input-dependent, reproducible output variation (routing, gating, logic) → P15. If waves simply annihilate or pass through regardless of input configuration → P13.
  - TE ≈ 0 → P13 (definitive).
  - TE > 0 but no demonstrated functional transformation → classify P13, flag P15 candidate.
  - TE > 0 AND demonstrated functional transformation → P15.
- P12: Excitable states → P13. Nontransitive species → P12.

**Co-occurrence:**
- *Allowed:* P9 (synchronized excitable media)
- *Excluded:* P15 (propagation geometry vs computation — mutually exclusive for the same collision set, but same system could have both annihilative and computational collisions in different regions)

---

### P14 — Self-organized criticality

**Observable scope:** Model-metadata required (must verify self-tuning, not external parameter tuning)

**Required raw observables:**
- `avalanche_sizes`: list of ints
- `avalanche_durations`: list of ints
- `activity[t]`: topplings per driving step

**Preprocessing:**
1. Discard T_burn = max(1000, 10% of total). 2. Empirical survival function. 3. PSD via Welch.

**Primary metric: Power-law avalanche distribution**
MLE fit P(s) ∼ s^(−τ). Clauset et al. bootstrap p-value.

**Secondary metrics:**
- 1/f noise: β ∈ (0.5, 1.5)
- Finite-size scaling collapse
- Duration scaling T ∼ s^γ
- Alternative-distribution LR tests: Clauset et al. likelihood-ratio comparing power-law against log-normal, exponential, and stretched exponential. Power-law must not be rejected in favor of any alternative.

**Null models:**
- *Subcritical dissipative:* p_diss = 0.2 → exponential distribution.
- *Shuffled-cascade:* permute sizes in time → destroys temporal correlations.

**Detection tiers:**
- *Screening:* Power-law fit p > 0.05 AND τ in plausible range
- *Confirmation:* p > 0.1 AND τ in model-class range AND LR tests do not favor log-normal, exponential, or stretched exponential over power-law (all three LR p > 0.1). AND (β ∈ (0.5, 1.5) OR duration scaling consistent). Verified self-tuning (model metadata).
- *Definitive:* Confirmation + dissipative null exponential + finite-size scaling collapse + P8 exclusion

**Common false positives:**
- Log-normal masquerading. Truncated exponentials. P8 (characteristic scale). External tuning (not SOC).

**Nearest-neighbor exclusions:**
- P8: Characteristic wave speed → P8. Scale-free → P14.
- P15: If TE > 0 across structure collisions → flag P15 co-occurrence.

**Co-occurrence:**
- *Allowed:* P15 (SOC can enable computation — co-detection expected)
- *Excluded:* P8 (characteristic vs scale-free)

---

## Cluster E: Information Processing and Collective Cognition

### P15 — Persistent propagating computation

**Observable scope:** State-history only (logic gate demo may use specific ICs but detection is from observed history)

**Required raw observables:**
- `grid[t]`: (rows, cols)
- Cell state alphabet size

**Preprocessing:**
1. Structure detection: connected components, track across frames. Persist ≥ T_persist = max(10, 0.1 × T_prop).
2. Collision detection: bounding box overlap/merge.
3. ROI extraction around collision sites.

**Primary metric: Transfer Entropy across collisions (necessary condition)**
TE(input → output) per collision. KSG or plug-in estimator.

**Functional discriminator (required for P15, not just TE):**
Collision outcomes must vary systematically with input configuration — demonstrating routing, gating, or logic. Tested by: holding one input constant and varying the other; if output changes predictably → functional transformation confirmed.

**Secondary metrics:**
- Transient length: T_transient / N
- Structure census (distinct types)
- Logic gate truth table fidelity

**Null models:**
- *Class I/II CA:* no structures, TE ≈ 0.
- *Subcritical-density:* structures die.
- *Excitable-wave control (P13):* Greenberg-Hastings. Waves propagate, annihilate. TE ≈ 0.

**Detection tiers:**
- *Screening:* TE > 0 for > 10% of collisions (permutation p < 0.05) AND T_transient > N
- *Confirmation:* TE > 0 for > 25% (p < 0.01) AND T_transient > 10·N AND functional discriminator shows input-dependent output
- *Definitive:* Confirmation + logic gate demo with truth table fidelity > 0.9 + Greenberg-Hastings control shows TE ≈ 0

**Common false positives:**
- P13 (annihilative collisions, TE ≈ 0). Chaotic rules without structure. Transient artifacts.

**Nearest-neighbor exclusions (revised P13/P15 boundary):**
- P13: Two-stage test (see P13 card):
  - TE ≈ 0 → P13
  - TE > 0 without functional transformation → P13, P15 candidate flag
  - TE > 0 with functional transformation → P15
- P14: SOC is Layer 2A descriptor. Can co-occur.

**Co-occurrence:**
- *Allowed:* P14 (SOC enabling computation)
- *Excluded:* P13 (for the same collision set — see note at P13)

---

### P16 — Associative memory / pattern completion

**Observable scope:** Model-metadata required (stored patterns and weights needed for overlap computation)

**Required raw observables:**
- `s[t]`: (N,) node states
- `P`: (K, N) stored patterns
- `s[0]`: corrupted cue
- `W`: (N, N) weights

**Preprocessing:**
1. Overlap m_μ(t) = (1/N) Σ s_i · P_{μ,i}. 2. Identify target (max initial overlap). 3. Track all overlaps.

**Primary metric: Pattern completion accuracy**
Final overlap m_target(T_final).

**Secondary metrics:**
- Basin size (fraction of corruptions converging correctly)
- Storage capacity K/N
- Spurious attractor rate

**Null models:**
- *Random-weight:* spurious attractors.
- *Zero-weight:* no dynamics, accuracy = initial overlap.

**Detection tiers:**
- *Screening:* Accuracy > 0.8 from ≤ 50% corruption
- *Confirmation:* Accuracy > 0.9 AND above random-weight null (p < 0.01) AND K > 1 retrievable patterns
- *Definitive:* Confirmation + basin sizes measured for ≥ 3 patterns + P25 exclusion

**Common false positives:**
- Trivial attractors. Single dominant attractor. P25 (one target, many ICs).

**Nearest-neighbor exclusions:**
- P25: Single target → P25. Multiple selectively retrievable → P16.

**Co-occurrence:**
- *Allowed:* P25 (same system may show both equifinality and memory)
- *Excluded:* (no hard exclusions beyond P25 edge case)

---

### P17 — Distributed sensing / collective gradient detection

**Observable scope:** Model-metadata required (gradient field and noise level needed for CI computation)

**Required raw observables:**
- `positions[t]`: (N, 2)
- `velocities[t]`: (N, 2)
- Gradient field and direction ĝ
- σ_sense

**Preprocessing:**
1. CI_i = v̂_i · ĝ. 2. CI_group = V̂_cm · ĝ. 3. CI_group vs N across trials.

**Primary metric: Superlinear sensing**
CI_group scales with N above individual baseline.

**Secondary metrics:**
- Navigation accuracy (within 30° of ĝ) vs N
- Time-to-detection vs N
- Collective gain CI_group / CI_individual

**Null models:**
- *Zero-gradient:* CI ≈ 0.
- *No-interaction:* CI_group = CI_individual for all N.

**Detection tiers:**
- *Screening:* CI_group > CI_individual at N ≥ 10 (p < 0.05)
- *Confirmation:* Paired test p < 0.01 across group sizes AND monotonic increase with N
- *Definitive:* Confirmation + collective gain exceeds √N (genuine sensing, not noise averaging) + P5 exclusion

**Common false positives:**
- P5 (alignment without gradient response). Simple noise averaging (gain = √N).

**Nearest-neighbor exclusions:**
- P5: No gradient correlation → P5. Gradient-aligned → P17. Can co-occur.

**Co-occurrence:**
- *Allowed:* P5 (flocking + sensing is the canonical case), P19 (leadership + sensing)
- *Excluded:* (none — P17 is always layered on top of a motion pattern)

---

## Cluster F: Decision-Making and Social Dynamics

### P18 — Collective consensus / decision-making

**Observable scope:** Model-metadata assisted (option quality needed for accuracy measurement)

**Required raw observables:**
- `choice[t]`: (N,) in {1..K}
- `quality`: (K,) (if available)
- `f_k(t)`: commitment fractions

**Preprocessing:**
1. f_k(t) per option. 2. Consensus: max f_k > 0.9. 3. Accuracy if quality known.

**Primary metric: Convergence to consensus**
t_consensus, accuracy, final committed fraction.

**Secondary metrics:**
- Sensitivity (accuracy vs Δq)
- Wisdom-of-crowds gain

**Null models:**
- *Independent-choice:* random, accuracy = 1/K.
- *Majority-rule:* Condorcet baseline.

**Detection tiers:**
- *Screening:* max f_k > 0.8 in ≥ 70% of trials
- *Confirmation:* max f_k > 0.9 in ≥ 90% AND accuracy > independent null (p < 0.01)
- *Definitive:* Confirmation + P20/P21/P22 exclusions all cleared

**Common false positives:**
- P21 (multimodal). P22 (one-way spread). P20 (density trigger).

**Nearest-neighbor exclusions:**
- P21: Multimodal final → P21. Unimodal → P18.
- P22: Bidirectional influence → P18. One-way → P22.
- P20: Quality-triggered → P18. Density-triggered → P20.

**Co-occurrence:**
- *Allowed:* P19 (leadership-guided consensus), P17 (sensing-informed consensus)
- *Excluded:* P21 (consensus vs fragmentation, same decision variable), P22 (deliberation vs cascade)

---

### P19 — Emergent leadership / minority guidance

**Observable scope:** Model-metadata required (informed labels and preferred direction needed)

**Required raw observables:**
- `positions[t]`, `velocities[t]`
- `informed`: (N,) boolean
- `preferred_dir`: (2,) unit vector

**Preprocessing:**
1. Separate subpopulations. 2. Heading per subpopulation. 3. Cross-correlation and TE.

**Primary metric: Influence asymmetry**
TE(minority→majority) / TE(majority→minority) > 2.0. Guidance efficacy > 0.7 with < 20% informed.

**Secondary metrics:**
- Critical informed fraction
- Mechanism verification (no special signal)
- Majority response lag

**Null models:**
- *No-informed:* all naive (pure P5).
- *Labeled-but-uninformed:* no preferred direction.

**Detection tiers:**
- *Screening:* Guidance efficacy > 0.5 with informed < 30%
- *Confirmation:* TE ratio > 2.0 AND efficacy > 0.7 with informed < 20%
- *Definitive:* Confirmation + labeled-but-uninformed null shows no guidance + P5/P32 exclusions

**Common false positives:**
- P5 (no informed subset). P18 (symmetric pooling). P32 (persistent roles).

**Nearest-neighbor exclusions:**
- P5: No informed subset → P5. Disproportionate minority → P19.
- P32: Persistent role → P32. Context-dependent → P19.

**Co-occurrence:**
- *Allowed:* P5 (leadership implies flocking), P18 (guided consensus)
- *Excluded:* P32 (transient vs persistent roles for same agents)

---

### P20 — Quorum sensing / threshold-activated collective response

**Observable scope:** Model-metadata required (density/signal mechanism must be verified)

**Required raw observables:**
- `active[t]`: (N,) boolean
- `density[t]`
- `f_active(t)`

**Preprocessing:**
1. Sweep density, record f_active. 2. Upward + downward sweep. 3. Hill function fit.

**Primary metric: Sharp transition**
Hill coefficient n > 2.

**Secondary metrics:**
- Hysteresis width
- Response onset time
- Bistability range

**Null models:**
- *Independent-activation:* n = 1.
- *Linear-response:* f ∝ ρ.

**Detection tiers:**
- *Screening:* n > 1.5
- *Confirmation:* n > 2 AND (hysteresis detected OR bistability range > 0) AND above independent null (p < 0.01)
- *Definitive:* Confirmation + density-mechanism verified from model metadata + P18/P22 exclusions

**Common false positives:**
- P18 (quality trigger). P22 (irreversible spread). Continuous transition (n ≈ 1).

**Nearest-neighbor exclusions:**
- P18: Density trigger → P20. Quality trigger → P18.
- P22: Reversible bistable → P20. Irreversible → P22.

**Co-occurrence:**
- *Allowed:* P9 (quorum triggers synchronization)
- *Excluded:* P18 (different trigger mechanism for same collective switch), P22 (reversible vs irreversible)

---

### P21 — Polarization / fragmentation

**Observable scope:** Model-metadata assisted (interaction threshold ε useful but distributional analysis works without it)

**Required raw observables:**
- `opinion[t]`: (N,) continuous or discrete
- Interaction rule / ε (if bounded confidence)

**Preprocessing:**
1. Distribution histogram per timestep. 2. GMM fit (BIC). 3. Dip test. 4. Between-cluster distance.

**Primary metric: Persistent multimodality**
Dip test p < 0.01 at steady state.

**Secondary metrics:**
- Between-cluster distance > ε
- Time to fragmentation
- Final opinion variance

**Null models:**
- *Full-range (ε→∞):* consensus.
- *Random-opinion:* uniform, not clustered.

**Detection tiers:**
- *Screening:* Dip p < 0.05
- *Confirmation:* Dip p < 0.01 AND ≥ 2 clusters AND persistence ≥ 50 × mean interaction time
- *Definitive:* Confirmation + ≥ 100 × interaction time + full-range null converges + emerged from initially unimodal IC + P18 exclusion

**Common false positives:**
- P18 (unimodal). P1 (spatial, not opinion). Initial bimodality artifact.

**Nearest-neighbor exclusions:**
- P18: Unimodal → P18. Multimodal → P21.
- P22: Stable partition → P21. One-way spread → P22.

**Co-occurrence:**
- *Allowed:* P1 (opinion + spatial fragmentation can co-occur)
- *Excluded:* P18 (consensus vs fragmentation, same decision variable)

---

### P22 — Information cascade / social contagion

**Observable scope:** State-history only

**Substrate:** lattice_2d (implemented); generalizes to any substrate with
discrete states and local transition rules (e.g., network with adjacency).

**Required raw observables:**
- `grid[t]`: (rows, cols) int array with discrete states, including at
  least "susceptible" (state 0) and "infected/active" (state 1)
- `grid_dims`: (rows, cols)
- `newly_infected` per step (optional; derivable from grid deltas)

> **Note:** An earlier version of this card specified a network-centric
> API (`adopted[t]`, `t_adopt`, adjacency matrix, cascade-size distribution
> across trials). The implementation in `epc/detectors/p22_information_cascade.py`
> uses the lattice-based approach below. The network-based formulation is
> conceptually equivalent but operates on a different substrate; this card
> documents what the code actually measures.

**Preprocessing:**
1. Extract infection time map `T_inf[i, j] = t*` where `t*` is the first
   step at which cell `(i, j)` entered the infected state.
2. Time series: `i_count[t]`, `r_count[t]` (recovered count where
   applicable), `newly_infected[t]`.
3. Compute cascade reach at end of trajectory: fraction of cells that
   ever left the susceptible state.

**Primary metric: Moran's I on the infection TIME map**

Spatial autocorrelation of `T_inf`. Cells near the seed have small `t*`;
distant cells have large `t*`. A propagating wavefront produces strong
spatial clustering (I ≈ 0.98 for a circular wavefront). Random assignment
of infections to cells produces I ≈ 0.

This metric is preferred over Moran's I on the final binary state because
when the epidemic infects most cells the final state is nearly uniform
(I → 0 by construction) while the infection time map retains the full
wavefront structure regardless of final reach.

**Secondary metrics:**
- `cascade_reach_total` = (R∞ + I∞) / N — fraction ever infected
- `r0_estimate` — early-growth exponential fit to `newly_infected[t]`
- `peak_infected_fraction` — max(I[t]) / N
- `is_unimodal` — binary flag: single peak in I(t) curve
- `epidemic_died_out` — binary flag: I[-1] = 0
- `wavefront_velocity` — secondary estimate from radial expansion

**Null model:**

*Random-timing null:* preserve the NUMBER of cells that get newly infected
at each timestep, but randomly assign which cells get infected. This
destroys spatial correlation while preserving the overall epidemic curve
shape (same `i_count[t]` trajectory expected value). Re-compute Moran's I
on the null infection time map. n_permutations = 199 (floor p = 0.005).

For the canonical positive (SIR epidemic), null Moran's I ≈ 0 ± 0.01,
giving Cohen's d > 100 against observed.

**Detection tiers** (implemented in `P22CascadeDetector`):

- *Screening:* `cascade_reach_total ≥ 0.05` AND `moran_i_time ≥ 0.10`
- *Confirmation:* Screening + `p < 0.01` AND `is_unimodal` AND
  `r0_estimate > 1.0`
- *Definitive:* Confirmation + `cascade_reach_total ≥ 0.30` AND
  `epidemic_died_out = True` (single-pass wave, not persistent activity)

**Common false positives:**

- P13 (excitable waves): produces wavefronts but waves are persistent and
  re-entrant, not single-pass. Exclusion check: if activity does NOT die
  out in the last quarter of the trajectory, mark P13 inconclusive.
- P1 (aggregation): can co-occur with P22 (cascade creates transient
  clustering), but P1 alone does not imply cascade dynamics.

**Nearest-neighbor exclusions:**

- P13: persistent re-entrant waves (cells re-excited many times) vs.
  single-pass cascade (each cell infected at most once)
- P1: static spatial clustering (convergence to an aggregated state) vs.
  dynamic cascade (the spreading process itself)

**Co-occurrence:**

- *Allowed:* P1 (spatially-structured cascade produces transient
  aggregation; co-occurrence expected on SIR during the epidemic peak)
- *Excluded:* P13 (transient vs. persistent is definitional)

**Verified on:** SIR epidemic CA × P22 → DEFINITIVE (Moran I_time = 0.987,
Cohen's d = 109.5, p = 0.005, cascade_reach = 1.0, R0_est > 1, unimodal,
died out). See `tests/test_sir_p22_e2e.py::TestSIRP22Canonical`.

**Verified rejections:** GoL, Greenberg-Hastings, Nowak-May, Schelling all
correctly rejected by P22 screening (either low cascade reach or low
Moran's I on infection time map).

---

### P23 — Anti-coordination / emergent load balancing

**Observable scope:** State-history only

**Required raw observables:**
- `choice[t]`: (N,) in {1..K}
- `attendance_k(t)` per option
- `capacity_k`: (K,) option capacities (if applicable; default = N/K)

**Preprocessing:**
1. Attendance time series. 2. Variance. 3. Autocorrelation. 4. Mean absolute deviation from capacity: MAD_cap = ⟨|attendance_k − capacity_k|⟩.

**Primary metric: Reduced attendance variance**
Observed / random-choice variance ratio.

**Secondary metrics:**
- Negative lag-1 autocorrelation
- Nash equilibrium proximity
- **Capacity deviation:** MAD_cap / capacity_mean. Low ratio = balancing near capacity, not just low-variance around the wrong level. Must be < 0.2 for meaningful load balancing.

**Null models:**
- *Random-choice:* binomial variance.
- *Fixed-strategy:* maximum variance.

**Detection tiers:**
- *Screening:* Variance ratio < 0.7
- *Confirmation:* Ratio < 0.5 AND lag-1 autocorrelation significantly negative (p < 0.01) AND capacity deviation MAD_cap/capacity < 0.2
- *Definitive:* Confirmation + N ≥ 50 + P18/P32 exclusions

**Common false positives:**
- P18 (convergence to one option). P32 (fixed roles). Small-N artifacts.

**Nearest-neighbor exclusions:**
- P18: One option → P18. Distributed → P23.
- P32: Stable assignment → P32. Dynamic switching → P23.

**Co-occurrence:**
- *Allowed:* (typically standalone within decision-making cluster)
- *Excluded:* P18 (distributed vs concentrated), P32 (dynamic vs fixed)

---

## Cluster G: Resilience and Regulation

### P24 — Homeostatic regulation

**Observable scope:** Model-metadata required (perturbation protocol and set-point must be specified)

**Required raw observables:**
- `x[t]`: controlled variable
- `x_target`: set-point
- `p[t]`: perturbation signal

**Preprocessing:**
1. Deviation e(t) = x − x_target. 2. Perturbation onset/offset. 3. Recovery trajectory.

**Primary metric: Recovery time and accuracy**
t_rec < T_relax. Steady-state deviation < 10% of |p|.

**Secondary metrics:**
- Regulatory gain |p|/|e_ss|
- Perturbation robustness (t_rec vs |p|)
- Recovery completeness

**Null models:**
- *Open-loop:* feedback disabled → permanent displacement.
- *Random-response:* random corrections → deviation random-walks.

**Detection tiers:**
- *Screening:* System returns to within 20% of |p| within 2 × T_relax
- *Confirmation:* t_rec < T_relax AND deviation < 10% of |p| AND better than open-loop (p < 0.01)
- *Definitive:* Confirmation + active corrective mechanism verified + multiple perturbation magnitudes + P25 exclusion

**Common false positives:**
- P25 (IC convergence, not perturbation recovery). Passive dissipation. Transient recovery.

**Nearest-neighbor exclusions:**
- P25: Perturbation-recovery → P24. IC-diversity → P25.

**Co-occurrence:**
- *Allowed:* P25 (same system, different experimental designs), P30 (homeostatic autopoietic system)
- *Excluded:* (none hard — P24 and P25 test different things)

---

### P25 — Canalized restoration (equifinality)

**Observable scope:** Model-metadata required (target macrostate must be defined)

**Required raw observables:**
- `S[t]`: macrostate descriptor
- `S_target`: target
- K ≥ 20 diverse initial conditions
- Disruption protocol

**Preprocessing:**
1. Run from each IC. 2. d_k = dist(final, target). 3. Var({d_k}). 4. Disruption recovery.

**Primary metric: Low convergence variance**
Var({d_k}) small. > 90% of ICs converge within tolerance.

**Secondary metrics:**
- Basin volume vs tolerance
- Disruption recovery rate
- Path diversity (high trajectory variance, low endpoint variance)

**Null models:**
- *Random-dynamics:* uniform final states → high variance.
- *Small-basin:* most ICs don't converge.

**Detection tiers:**
- *Screening:* > 80% of ICs converge within tolerance
- *Confirmation:* > 90% AND variance < 10% of state-space AND better than random (p < 0.01)
- *Definitive:* Confirmation + disruption recovery demonstrated + target has nontrivial structure + P16/P24 exclusions

**Common false positives:**
- P24 (active regulation). Trivial attractor. P16 (multiple patterns, not one target).

**Nearest-neighbor exclusions:**
- P24: IC-diversity → P25. Perturbation-recovery → P24.
- P16: Edge case — if system has only one target macrostate from diverse ICs → P25 only. If it has multiple selectively retrievable patterns → P16 only. If both properties demonstrated in different experimental protocols (IC-diversity test shows equifinality; cued-recall test shows pattern selection) → both allowed.

**Co-occurrence:**
- *Allowed:* P24 (complementary tests), P16 (system with equifinal basin + memory)
- *Excluded:* (none hard)

---

### P26 — Stochastic resonance

**Observable scope:** Model-metadata required (subthreshold signal parameters A, f₀ needed)

**Required raw observables:**
- `y[t]`: output signal
- Signal: A, f₀
- Noise amplitude σ (swept)
- Performance measure

**Preprocessing:**
1. Run at ≥ 20 noise levels with ≥ 50 trials each. 2. Performance per σ.

**Primary metric: Inverted-U performance curve**
Peak at σ* > 0.

**Secondary metrics:**
- Peak SNR
- Performance gain at σ*
- Peak width

**Null models:**
- *No-signal (A=0):* monotonic decrease or flat.
- *Suprathreshold:* A above threshold → monotonic decrease.

**Detection tiers:**
- *Screening:* Performance at σ* > performance at σ=0 (p < 0.05)
- *Confirmation:* p < 0.01 AND inverted-U shape confirmed (negative 2nd derivative) AND σ* not at sweep boundary
- *Definitive:* Confirmation + no-signal null flat + suprathreshold null monotonic + ≥ 50 trials per level

**Common false positives:**
- Frequency resonance (sweep frequency, not noise). Coincidental peak (small sample). Non-monotonic from other mechanisms.

**Nearest-neighbor exclusions:**
- Unique niche. No direct neighbor confusion.

**Co-occurrence:**
- *Allowed:* P11 (noise-enhanced predator-prey oscillations), P18 (noise-enhanced consensus)
- *Excluded:* (none)

---

## Cluster H: Competition and Cooperation

### P27 — Spatial reciprocity / emergent cooperation

**Observable scope:** Model-metadata required (payoff matrix needed to verify PD structure)

**Required raw observables:**
- `strategy[t]`: (N,) C/D
- Lattice positions (fixed)
- Payoff matrix (T > R > P > S)

**Preprocessing:**
1. f_C(t). 2. Moran's I on cooperator indicator. 3. Well-mixed equilibrium.

**Primary metric: Cooperator survival above well-mixed baseline**
f_C > 0 long-term when well-mixed predicts f_C → 0.

**Secondary metrics:**
- Spatial autocorrelation of cooperation
- Cluster survival time
- Interface dynamics

**Null models:**
- *Well-mixed:* cooperators extinct.
- *Random-strategy:* f_C = 0.5, Moran's I ≈ 0.

**Detection tiers:**
- *Screening:* f_C > 0.02 at t > 100 × T_gen
- *Confirmation:* f_C > 0.05 at t > 1000 × T_gen AND above well-mixed (p < 0.01) AND Moran's I > 0 (p < 0.01)
- *Definitive:* Confirmation + PD structure verified (T>R>P>S) + N ≥ 100 + P1 exclusion

**Common false positives:**
- P1 (movement clustering, not selection). Non-PD payoffs. Finite-size drift.

**Nearest-neighbor exclusions:**
- P1: Movement → P1. Static lattice + selection → P27.
- P18: No payoff conflict → P18. Mixed-motive → P27.

**Co-occurrence:**
- *Signature overlap:* P1 — cooperator clusters produce positive spatial autocorrelation (P1-like Moran's I signal), but this is a byproduct of selection dynamics on a static lattice, not preference-driven movement. The P1 detector may fire on P27 output; this is a known false positive for P1, not genuine co-occurrence. Use the P1→P27 exclusion (movement vs selection) to resolve.
- *Allowed:* P12 (cyclic games + spatial reciprocity)
- *Excluded:* (none hard beyond P1 signature overlap)

---

### P28 — Wealth condensation / spontaneous inequality

**Observable scope:** State-history only

**Required raw observables:**
- `wealth[t]`: (N,) non-negative
- Exchange events
- W_total (conserved)

**Preprocessing:**
1. Gini G(t). 2. Wealth distribution; Pareto tail. 3. Oligarchic fraction.

**Primary metric: Gini toward 1.0**
Monotonically increasing. Spearman ρ(G, t) > 0.9.

**Secondary metrics:**
- Pareto exponent α
- Oligarchic fraction
- Wealth mobility (rank correlation)

**Null models:**
- *Random redistribution:* Gini fluctuates, no trend.
- *Non-interacting:* Gini constant.

**Detection tiers:**
- *Screening:* G_final > 0.6 AND Spearman > 0.7
- *Confirmation:* G > 0.8 AND Spearman > 0.9 AND Pareto tail (p > 0.1)
- *Definitive:* Confirmation + redistribution null flat + identical initial agents verified + P14/P32 exclusions

**Common false positives:**
- P14 (stationary power law). P32 (functional roles). Pre-existing heterogeneity.

**Nearest-neighbor exclusions:**
- P14: Stationary → P14. Trending to absorption → P28.
- P32: Functional roles → P32. Resource quantity → P28.

**Co-occurrence:**
- *Allowed:* P22 (cascade dynamics driving condensation)
- *Excluded:* P14 (stationary vs non-stationary), P32 (resource vs role)

---

## Cluster I: Emergent Structure Formation

### P29 — Trail / network formation

**Observable scope:** State-history only

**Required raw observables:**
- `positions[t]`: (N, 2)
- `trail[t]`: (rows, cols)
- `targets`: (M, 2)

**Preprocessing:**
1. Threshold trail → binary network. 2. Skeletonize → graph. 3. Efficiency, Steiner ratio, fault tolerance.

**Primary metric: Network efficiency**
Relative to optimal. Steiner ratio.

**Secondary metrics:**
- Fault tolerance (connectivity after 20% edge removal)
- Trail reinforcement dynamics

**Null models:**
- *Random-network (ER):* same nodes/edges.
- *Direct-path:* straight lines, high total length.

**Detection tiers:**
- *Screening:* Efficiency > 0.6 × optimal
- *Confirmation:* > 0.8 × optimal AND Steiner ratio < 1.5 AND fault tolerance > ER null (p < 0.01)
- *Definitive:* Confirmation + emerged on uniform substrate + P4/P7 exclusions

**Common false positives:**
- P4 (exclusive territories). Pre-existing channels. P7 (dynamic flow lanes).

**Nearest-neighbor exclusions:**
- P4: Shared → P29. Exclusive → P4.
- P7: Persistent after agents stop → P29. Vanishes → P7.

**Co-occurrence:**
- *Allowed:* P32 (trail network + specialized foragers)
- *Excluded:* P4 (shared vs exclusive), P7 (persistent vs dynamic)

---

### P30 — Spontaneous boundary formation (autopoiesis)

**Observable scope:** State-history only

**Required raw observables:**
- `positions[t]`, `types[t]`
- `conc[t]`: (rows, cols)

**Preprocessing:**
1. Identify closed boundary (Euler characteristic). 2. Interior/exterior. 3. Concentration gradient.

**Primary metric: Topological closure + gradient**
χ = 1. Interior/exterior ratio > 2.0.

**Secondary metrics:**
- Self-repair (after 10% removal)
- Semi-permeability
- Boundary stability duration

**Null models:**
- *No-boundary:* uniform concentration.
- *Random-placement:* no persistent closure.

**Detection tiers:**
- *Screening:* Topological closure detected for ≥ 5 × T_cross
- *Confirmation:* Closure ≥ 50 × T_cross AND concentration ratio > 2.0 AND above random (p < 0.01)
- *Definitive:* Confirmation + self-repair demonstrated + emerged from dispersed agents + P1 exclusion

**Common false positives:**
- P1 (clusters without closure). Static boundary. Transient enclosure.

**Nearest-neighbor exclusions:**
- P1: Closure → P30. Open clusters → P1.

**Co-occurrence:**
- *Allowed:* P24 (homeostatic autopoietic system), P32 (membrane specialists)
- *Excluded:* P1 (closure vs open)

---

## Cluster J: Emergent Agent-Level Competencies

### P31 — Delayed gratification (provisional)

**Observable scope:** State-history only

**Required raw observables:**
- `positions[t]`: (N,) grid indices
- `types`: (N,)
- `positions[T]`: final positions
- `grid_dims`

**Preprocessing:**
1. d_i(t) = Manhattan dist to final. 2. Δd_i(t) = d_i(t+1)−d_i(t). 3. DG event: Δd > 0. 4. DG_i = (#DG events) / (#moves). 5. Flag frozen agents.

**Primary metric: Population DG distribution**
⟨DG⟩, σ_DG, quartiles. Correlation with sorting efficiency. Per-algorithm profiles.

**Secondary metrics:**
- Temporal concentration (CV of inter-DG intervals)
- DG spatial clustering (Moran's I on DG indicator)
- Frozen-cell fraction

**Null models:**
- *Random-walk:* random swaps. DG ≈ 0.5.
- *Marginal-shuffle:* permute each agent's move sequence. Preserves marginals, destroys timing.

**Detection tiers:**
- *Screening:* ⟨DG⟩ significantly < 0.5 AND > 0.05 (z-test vs random walk, p < 0.05)
- *Confirmation:* p < 0.01 AND σ_DG > 0.05 AND efficiency correlates with DG (p < 0.01) AND sorting completion ≥ 90%
- *Definitive:* Confirmation + non-redundancy test passed (see below)

**Non-redundancy protocol (P31 survival criterion):**
1. *Baseline:* predict outcome from aggregation features only (Moran's I trajectory, cluster count/size, path length, segregation index, endpoint quality)
2. *Extended:* add DG features (mean, variance, quantiles, temporal concentration, spatial clustering)
3. *Ablation:* DG from shuffled timing (marginals preserved)
- 5-fold CV, paired t-test p < 0.05
- P31 survives iff Extended > Baseline AND Extended > Ablation

**Common false positives:**
- Random walk artifact. Aggregation redundancy. Endpoint bias (require completion ≥ 90%).

**Nearest-neighbor exclusions:**
- P1: Non-redundancy test directly evaluates. Extended ≤ Baseline → P31 is P1 redescription.

**Co-occurrence:**
- *Allowed:* P1 (DG occurs during aggregation — co-detection expected), P32 (specialized roles + DG)
- *Excluded:* (P31 is provisional; if non-redundancy fails it collapses to P1)

---

### P32 — Emergent specialization (division of labor)

**Observable scope:** State-history only (identical initial agents can be verified from state history if initial state recorded)

**Required raw observables:**
- `task_history[t]`: (N, K) cumulative counts
- `efficiency[t]`: global metric
- Agent identity

**Preprocessing:**
1. H_i(t) = −Σ p_{ik} log p_{ik}. 2. ⟨H⟩(t). 3. Switching frequency.

**Primary metric: Behavioral entropy reduction**
⟨H⟩ decreasing. ⟨H⟩_final < 50% of initial.

**Secondary metrics:**
- Efficiency gain > 10% over non-specialized
- Switching frequency decreasing
- Role persistence (autocorrelation)
- Role count (k-means silhouette)

**Null models:**
- *Random-task:* H = H_max.
- *Identical-response:* no differentiation.

**Detection tiers:**
- *Screening:* ⟨H⟩_final < 70% of initial
- *Confirmation:* < 50% AND efficiency gain > 10% AND above random null (p < 0.01)
- *Definitive:* Confirmation + identical initial agents verified + P19/P23/P28 exclusions

**Common false positives:**
- P19 (transient role). P28 (resource, not function). Pre-existing heterogeneity. P23 (dynamic switching).

**Nearest-neighbor exclusions:**
- P19: Transient → P19. Persistent → P32.
- P23: High switching → P23. Decreasing switching → P32.
- P28: Resource quantity → P28. Functional role → P32.

**Co-occurrence:**
- *Allowed:* P31 (specialized + DG), P29 (trail builders + specialists), P1 (specialists clustering)
- *Excluded:* P19 (transient vs persistent for same agents), P23 (dynamic vs stable)

---

## Cross-Reference: Exclusion Graph

```
Cluster A
  P1  ⊗ P2, P3, P30
  P2  ⊗ P1
  P3  ⊗ P1, P12
  P4  ⊗ P1, P29

Cluster B
  P5  ⊗ P6, P7
  P6  ⊗ P5
  P7  ⊗ P5, P29
  P8  ⊗ P14

Cluster C
  P9  ⊗ P10
  P10 ⊗ P9
  P11 ⊗ P12
  P12 ⊗ P11, P3

Cluster D
  P13 ⊗ P15 (two-stage: TE + functional test)
  P14 ⊕ P15 (co-occurrence allowed)

Cluster E
  P15 ⊗ P13 (two-stage)
  P16 ⊗/⊕ P25 (edge case: if only one target → P25 only; if multiple selective → P16 only; if both properties in different protocols → both allowed)
  P17 ⊕ P5 (co-occurrence expected)

Cluster F
  P18 ⊗ P21 (same variable), P22
  P19 ⊗ P32 (same agents)
  P20 ⊗ P18, P22
  P21 ⊗ P18 (same variable)
  P22 ⊗ P20, P21
  P23 ⊗ P18, P32

Cluster G
  P24 ⊕ P25 (complementary)
  P25 ⊗/⊕ P16 (see P16 edge case above)
  P26 — (unique niche)

Cluster H
  P27 ⊗ P1 (signature overlap — P1 fires as false positive on P27 cooperator clusters), P18
  P28 ⊗ P14, P32

Cluster I
  P29 ⊗ P4, P7
  P30 ⊗ P1

Cluster J
  P31 ⊗ P1 (pending non-redundancy)
  P32 ⊗ P19, P23, P28
```

Legend: ⊗ = must exclude before declaring. ⊕ = co-occurrence expected/allowed.

---

## Observable Scope Summary

| Scope | Patterns |
|---|---|
| **State-history only** | P1, P4, P5, P6, P7, P8, P13, P15, P22, P23, P28, P29, P30, P31, P32 |
| **Model-metadata assisted** | P3, P9, P10, P11, P18, P21 |
| **Model-metadata required** | P2, P12, P14, P16, P17, P19, P20, P24, P25, P26, P27 |

This matters for real-world applicability: state-history-only detectors can
be applied to empirical data from unknown systems. Model-metadata-required
detectors are limited to simulation contexts where the rule set is known.

---

## Methodological Notes

### On threshold calibration
All thresholds are initial targets. After implementation, validate against
canonical positive (expected model) and negative (null model) cases. Target:
>95% sensitivity on canonical positives, >95% specificity on nulls.

### On multi-pattern detection
The exclusion graph identifies patterns to rule out, not patterns that
cannot co-occur. A system can legitimately trigger multiple detectors.
The co-occurrence fields make allowed combinations explicit.

### On null model philosophy
Each card specifies ≥ 2 nulls: mechanistic (removes mechanism, tests
necessity) and shuffle (preserves marginals, destroys structure). Both
required for definitive-tier detection.

### On observable naming
Standard state-history keys for cross-model metric application:
`positions`, `velocities`, `types`, `grid`, `wealth`, `opinion`,
`strategy`, `trail`, `scent`, `active`, `θ` (phases).
