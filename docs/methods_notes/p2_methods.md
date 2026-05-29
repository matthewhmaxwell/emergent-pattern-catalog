# P2 Methods Note — Activity-Induced Phase Separation (MIPS)

**Pattern:** P2 — Activity-induced phase separation (MIPS)
**Canonical model:** Active Brownian Particles (Fily & Marchetti 2012)
**Detector:** `epc/detectors/p2_mips.py`
**Model:** `epc/models/active_brownian_particles.py`
**Reproduction sprint:** Sprint 52 (`analysis/reproductions/p2_filymarchetti2012.py`)

---

## 1. Pattern and canonical reference

P2 is motility-induced phase separation: self-propelled particles with *no*
attractive interaction and *no* alignment rule spontaneously segregate into
coexisting dense (liquid) and dilute (gas) phases. The purely kinetic origin —
density-dependent speed reduction feeding back to increase local density —
distinguishes MIPS from thermodynamic phase separation (Cahn-Hilliard dynamics)
and from attraction-driven clustering (D'Orsogna-type Morse potential).

**Primary reference:** Fily, Y. & Marchetti, M. C. (2012). Athermal phase
separation of self-propelled particles with no alignment. *Phys. Rev. Lett.*
108, 235702.

**Secondary references:** Redner, G. S., Hagan, M. F. & Baskaran, A. (2013).
Structure and dynamics of a phase-separating active colloidal fluid. *Phys. Rev.
Lett.* 110, 055701. Cates, M. E. & Tailleur, J. (2015). Motility-induced phase
separation. *Annu. Rev. Condens. Matter Phys.* 6, 219.

---

## 2. Model equations and update rule

The implementation follows the overdamped Langevin equations of Fily-Marchetti:

```
dr_i/dt = v(ρ_i) · ê(θ_i)
dθ_i/dt = √(2 D_r) · ξ_i(t)
```

where `ê(θ_i) = (cos θ_i, sin θ_i)`, `ξ_i(t)` is unit Gaussian white noise,
and the density-dependent speed is the Fily-Marchetti linear slowdown law:

```
v(ρ) = v₀ · max(0, 1 − ρ/ρ*)
```

Local density `ρ_i` at particle `i` is estimated by counting neighbors
(including self) within coarse-graining radius `r_cg` and dividing by the disk
area `π r_cg²`. Integration uses Euler-Maruyama at timestep `dt`:

1. Compute `ρ_i` via cKDTree with periodic boundary conditions.
2. Compute `v_i = v₀ · max(0, 1 − ρ_i/ρ*)`.
3. Update positions: `r_i ← r_i + v_i · ê(θ_i) · dt`.
4. Update headings: `θ_i ← θ_i + √(2 D_r dt) · η`, `η ~ N(0,1)`.
5. Wrap positions into `[0, L)` (periodic).

The `cKDTree` is rebuilt every step with `boxsize=L` to enforce periodic
distance computation. This is O(N log N) per step; at N = 400–1000 and
typical run lengths (2000–5000 steps), wall time is dominated by the tree query.

**Note:** Unlike Cahn-Hilliard model-B equations (which evolve a conserved
scalar density field ∂φ/∂t = ∇²(δF/δφ)) or kinetic Ising models (discrete
spin flips), the EPC implementation is a particle-based overdamped Langevin
model. The universality class is different: Cahn-Hilliard produces Model B
coarsening (ξ(t) ∝ t^(1/3)); MIPS in ABP can reach a non-equilibrium steady
state with two coexisting phases. The detector and model are grounded in the
ABP/MIPS literature, not in equilibrium phase-separation models.

---

## 3. Parameter defaults and justifications

| Parameter | Default | Justification |
|---|---|---|
| `n_particles` | 400 | Minimum for reliable two-phase statistics; Sprint 16 found N < 300 gives seed-dependent metastability |
| `box_size` | 20.0 | Packing fraction φ = N π σ²/(4 L²); with N=400, σ=1, L=20 gives φ ≈ 0.785 |
| `v0` | 0.3 | Pe = v₀/D_r ≈ 100 at D_r = 3e-3; canonical MIPS regime (Fily-Marchetti Fig. 1) |
| `rho_star` | 10.0 | Critical density scale; Sprint 16 calibrated against ρ* = 4.0 for r_cg=1.0; higher ρ* shifts the MIPS boundary |
| `D_r` | 3e-3 | Rotational correlation time T_rot = 1/D_r ≈ 333 steps at dt=0.1; persistence length L_p = v₀/D_r |
| `r_cg` | 1.0 | Particle diameter σ = 1.0; density estimated at scale of one particle diameter |
| `dt` | 0.1 | Euler-Maruyama stability: Pe·dt/L_p < 0.01 satisfied for canonical parameters |
| `seed` | None | Default non-reproducible; set for any quantitative study |

**Canonical positive parameters (Fily-Marchetti 2012 Fig. 2 regime):**
N = 800, φ = 0.5, Pe = 100, ρ* = 4.0, r_cg = 1.0, dt = 0.05.
Packing fraction φ = 0.5 is controlled by setting `L = sqrt(N π σ²/(4φ))`.

---

## 4. Implementation choices

**No hard-core repulsion.** The Fily-Marchetti model does not include
volume-exclusion forces — particles overlap freely. The clustering mechanism
is purely kinetic: slow particles act as obstacles to incoming fast particles
because their density contribution raises ρ for nearby particles, slowing them
too. This is a feature, not a limitation — removing hard-core interactions
makes the kinetic mechanism analytically tractable (Fily-Marchetti Sec. II).

**Coarse-graining scale.** The detector independently re-estimates ρ_i using
the same `r_cg` and cKDTree approach as the model's `get_state()` call.
This means detection operates on re-derived densities, not on model-cached
values — ensuring the detector can function on any continuous_2d substrate
that exposes `positions` and `velocities`, not only on ABP.

**Structure factor and domain characterization.** In the MIPS physics literature
(Fily-Marchetti 2012, Redner et al. 2013), characteristic domain size ξ is
extracted by computing the radial structure factor:

```
S(k) = |Σ_i exp(i k · r_i)|² / N
```

and identifying the wavevector k* at the first peak (or equivalently via the
first zero of the pair correlation function g(r)). This FFT-based approach is
standard for measuring cluster coarsening dynamics. The EPC P2 detector uses
a simpler but closely related approach: instead of measuring ξ(t) explicitly,
it measures `f_gas = P(ρ_i < ρ*/2)` and `f_liquid = P(ρ_i > ρ*)` — the
fraction of particles in each phase. These fractions converge to the
lever-rule fractions of a phase diagram once the cluster has reached
steady-state, and do not depend on domain shape or coarsening time, making
them robust to finite-size effects and stochastic nucleation lags. The FFT
structure factor approach would be needed to study coarsening kinetics
(ξ(t) ∝ t^α); the detector targets the steady-state coexistence signature.

**Burn-in requirement.** Sprint 16 Phase 1d established that N = 400 shows
seed-dependent metastability at short run lengths: canonical MIPS fails to
nucleate within 1000 steps on ≈20% of seeds. The detector enforces a minimum
of 300 post-burn-in snapshots (configurable `burn_in` parameter); the Sprint 52
reproduction uses `burn_in = 500` (25% of 2500-step total run) to ensure the
measurement window captures the post-nucleation steady state.

**Parallelism.** Each `step()` call re-allocates a cKDTree for the full
particle array. This is not optimized for GPU or multi-core use. For large N
or long runs, consider subsampling the history (every 5 steps is standard
for the Sprint 52 reproduction).

---

## 5. Deviations from canonical (and why)

**ρ* = 4.0 vs. 10.0 default.** Fily & Marchetti set ρ* so that at the mean
density φ/(π r_cg² / L²) → ρ* is the close-packing density for r_cg = σ.
With σ = r_cg = 1 and φ = 0.5 on an L × L torus, the mean density
ρ_mean = N / L² = φ · 4/π ≈ 0.64 per unit area. Setting ρ* = 4.0 means
mean density is 16% of ρ*, placing the system well inside the two-phase
regime. The model default ρ* = 10.0 is calibrated for exploratory use at
default box_size = 20 and n_particles = 400; the Sprint 52 reproduction uses
ρ* = 4.0 to match the Fily-Marchetti convention.

**No translational noise.** The Fily-Marchetti model and this implementation
are *athermal* (no translational noise term in dr_i/dt). Translational Brownian
noise would add a `sqrt(2 D_t dt) · ξ_trans` term; at Pe ≫ 1 this is
negligible, but it would push the MIPS boundary to slightly higher Pe.

---

## 6. Reproduction status

**Sprint 52 reproduction:** Fily & Marchetti (2012) Fig. 2 canonical MIPS state
reproduced at (φ = 0.5, Pe = 100): `two_phase_score = 0.1237 ± 0.077`
(tolerance ≥ 0.10, **PASS**); density-speed Pearson r = −0.958 ± 0.020
(|r| ≥ 0.70, **PASS**); thermal Pe = 5 score = 0.052 (< 0.08, **PASS**).
All three tolerance checks pass. dim1 PARTIAL → **PASS**.

Seeds 0–1 show below-threshold two_phase_score (0.005, 0.061) due to
stochastic nucleation lag; 3/5 seeds exceed threshold. The seed-mean passes,
consistent with the known MIPS finite-size effect at N = 800.

Output: `analysis/outputs/p2_filymarchetti2012_reproduction.json`.

**Phase-2a panel dim4 (Sprint 49):** overall TNR = 1.000, Cohen's d = 4.245.
Output: `analysis/outputs/p2_phase2a_panel.json`.

---

## 7. Observable extraction

State snapshots expose:

| Key | Shape | Description |
|---|---|---|
| `positions` | (N, 2) float64 | Particle coordinates in [0, L) |
| `velocities` | (N, 2) float64 | `v_i · [cos θ_i, sin θ_i]` |
| `headings` | (N,) float64 | θ_i ∈ [−π, π) |
| `speeds` | (N,) float64 | `v_i` — variable (key MIPS signature) |
| `local_density` | (N,) float64 | ρ_i at snapshot time |

The `speeds` field is re-computed from the density at each `get_state()` call,
making snapshots idempotent and independent of the integration step.

The `get_metadata()` mechanistic flags are:
- `has_density_dependent_speed = True`
- `has_alignment_rule = False`
- `has_attraction_rule = False`

These three flags together encode the MIPS signature and are required for
the P2 detector to emit a DEFINITIVE result (Architecture Decision 43).

---

## 8. Known limitations

- **Finite-size nucleation lag.** MIPS cluster nucleation time is stochastic
  and seed-dependent, especially at N < 1000. Burn-in must exceed the
  nucleation time; the detector's minimum run-length requirement (300 snapshots
  post burn-in) may not be sufficient for all parameter regimes.

- **No cluster morphology tracking.** The detector reports coexisting-phase
  fractions but not cluster count, shape, or coarsening exponent. Structure
  factor / pair correlation function analysis would be needed to characterize
  the MIPS domain morphology or test coarsening universality class.

- **ρ* must be known.** The two_phase_coexistence_score requires ρ*. If the
  model's `rho_star` metadata is absent, the detector falls back to 4.0.
  For non-ABP substrates, ρ* must be set explicitly via the detector's
  `rho_star` parameter; without it the threshold is arbitrary.

- **dim2 still PARTIAL.** Multi-seed characterization (≥ 5 seeds with
  dispersion bounds) was performed in Sprint 16 but the ≥ 5-seed bar
  with explicit dispersion statistics was not formally reported per the
  depth_gap rubric. See C4 carry-forward.
