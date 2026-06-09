# P20 Methods Note — Quorum Sensing / Threshold-Activated Collective Response

## 1. Pattern definition

**P20 — Quorum sensing** describes a binary collective behavioral switch
triggered when agent density crosses a critical threshold. The system is OFF
below the threshold and ON above it, with a sharp (step-like) transition.
Positive feedback (e.g., autoinducer-enhanced production) creates hysteresis:
the activation density is higher than the deactivation density.

**Canonical reference:** Waters & Bassler (2005) — bacterial quorum sensing;
Nealson, Platt & Hastings (1970) — *Vibrio fischeri* luminescence.

**Distinctness from catalog neighbours:**
- P18 (consensus): P18 is graded convergence among competing alternatives;
  P20 is a binary on/off toggle with no competing alternatives.
- P14 (SOC): P14 is avalanche criticality with power-law statistics; P20 is a
  single population-level switch at a fixed density threshold.
- P1 (aggregation): P1 measures spatial clustering; P20 measures a density-
  dependent collective state transition with hysteresis.

## 2. Model implementations

### 2.1 AutoinducerQuorum (canonical)

Mean-field bacterial quorum-sensing model. Population-level ODE for
autoinducer concentration C:

  dC/dt = density × [f_off × r_base + f_on × (r_base + r_enh)] − γ × C

where f_on is the fraction of ON agents, r_base = 1.0 and r_enh = 4.0
are production rates, and γ = 1.0 is the degradation rate.

Agent switching:
- OFF → ON when C > activation_threshold (1.5)
- ON → OFF when C < deactivation_threshold (1.0)
- Per-agent Gaussian noise (σ = 0.02) on effective concentration

At equilibrium with all OFF: C_eq = density × r_base / γ = density.
Activation occurs at density_ON = activation_threshold = 1.5.

At equilibrium with all ON: C_eq = density × (r_base + r_enh) / γ = 5 × density.
Deactivation occurs at density_OFF = deactivation_threshold / 5 = 0.2.

Hysteresis width ≈ 1.3 (theoretical), measured ≈ 1.19.

The model sweeps density from 0.1 to 3.0 (40 levels up, 40 down), with 300
equilibration steps per density level. Only the last 50% of steps are used
for the equilibrium curve.

### 2.2 FractionThresholdModel (T1b cross-model)

Simple fraction-threshold model without biochemical dynamics. Each agent
switches ON when the density-amplified effective ON fraction exceeds
activation_fraction (0.5), and OFF when it drops below deactivation_fraction
(0.2). A density-dependent seed term (0.35 × density) enables bootstrapping
the collective switch at high density without requiring prior ON agents.

This provides an independent implementation that exhibits the same P20
signature (sharp threshold + hysteresis) via a different mechanism.

### 2.3 GradedResponseModel (negative control)

Linear response: fraction_on = slope × density + noise. No sharp threshold,
no hysteresis. Used as a negative control — the P20 detector correctly
rejects this model at screening (no OFF→ON transition).

## 3. Detector design

### 3.1 T1a observation contract

The detector reads a density-sweep observation bundle via
`extract_observation_bundle(history)`. Required keys: `density`,
`concentration`, `collective_state`, `fraction_on`, `density_idx`,
`sweep_direction`. Any model or external system producing these keys
can be evaluated.

### 3.2 Primary metric: step-function R²

For each sweep direction, compute the equilibrium fraction_on at each
density level (last 50% of steps). Fit the best step function
f(d) = f_low if d < d_c, f_high if d ≥ d_c by exhaustive search over
all possible split points. The R² of this fit measures transition
sharpness: R² close to 1.0 for a perfect step, lower for gradual
responses.

### 3.3 Null model: density-shuffle

Permute the equilibrium fraction_on values across density levels and
re-fit the best step function. Under H0 (no density-dependent switch),
the shuffled R² should be modest (random step functions don't fit well).
The observed R² should far exceed shuffled R² for genuine quorum sensing.

p-value: fraction of null R² values ≥ observed R², floored at
1/(n_permutations + 1).

### 3.4 Hysteresis analysis

Compute the density at which fraction_on crosses 0.5 in the up-sweep
(activation density) and down-sweep (deactivation density). Hysteresis
is present when activation density > deactivation density by more than
5% of the total density range.

### 3.5 Tier criteria

| Tier | Criteria |
|------|----------|
| Screening | OFF→ON transition observed (f_high − f_low > 0.5) |
| Confirmation | Step R² > 0.7, relative transition width < 0.4, null p < 0.05 |
| Definitive | Hysteresis present, Cohen's d > 1.0, null p ≤ 0.005, metadata has_threshold_activation |

### 3.6 Design decisions

- **Step-function R² over max_slope:** Max slope of shuffled (randomly
  permuted) fraction curves is artificially high due to discontinuous
  jumps between non-adjacent values. Step-function R² correctly captures
  the global structure of a sharp transition and is lower for shuffled
  data. (ADR: step-function R² > raw max-slope for density-shuffle null.)

- **Equilibrium last-50% averaging:** The model equilibrates over
  n_steps_per_density timesteps. Using only the last 50% avoids transient
  dynamics contaminating the equilibrium curve. This is analogous to
  burn-in removal in MCMC.

- **Hysteresis via crossing detection:** Linear interpolation between
  consecutive density levels to find the 0.5-crossing point. For the
  down-sweep, the curve is already sorted ascending by density; the
  crossing point gives the deactivation density directly.

## 4. Reproduction results

### dim1 (Sprint 83)
- Step-function R² = 1.000 (tolerance > 0.9: PASS)
- Critical density = 1.401
- Hysteresis width = 1.190 (tolerance > 0.1: PASS)
- Activation density = 1.401, deactivation density = 0.212
- Detector: DEFINITIVE (confidence 0.90)
- All 3 tolerance checks PASS.

### dim2 (Sprint 83)
- 20-seed campaign: R² = 1.000 ± 0.000 (CV ≈ 0%)
- Critical density = 1.401 ± 0.000
- Hysteresis width = 1.190 ± 0.000
- All 20 seeds DEFINITIVE
- CV ≈ 0%: the model operates in a high-SNR regime where stochastic
  noise (σ = 0.02) is negligible relative to threshold gaps (0.5).
  This is physically correct — bacterial QS is robust to noise.

### T1b (Sprint 83)
- FractionThresholdModel: DEFINITIVE (confidence 0.90)
- Step R² = 1.000, hysteresis width = 0.892
- Independent mechanism (fraction-threshold vs autoinducer ODE)
  produces the same P20 signature.
