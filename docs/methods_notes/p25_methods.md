# P25 — Canalized Restoration / Equifinality: Methods Note

## Pattern definition

**P25** detects equifinality: a system that converges to the same target
macrostate from a wide range of initial conditions, beyond what a simple
dominant attractor predicts. Canalization implies that diverse starting
points are funneled into the same developmental outcome via genuine
dynamical attractors, not trivial collapse.

**Distinctness from P24 (homeostatic regulation):** P24 is cybernetic —
continuous active regulation vs. sustained perturbation (deviation integral,
growth ratio). P25 is geometric — convergence from arbitrary starting
points (IC variance → final-state variance ratio).

**Distinctness from P16 (associative memory):** P16 has MULTIPLE
retrievable targets selectively recalled by different cues. P25 funnels
diverse ICs to ONE canalized target.

## Model implementations

### CanalizedLandscape (`epc/models/canalization.py`)

Multi-dimensional gradient-flow on a combined linear + quartic potential:

```
dx/dt = -α·(x - x*) - 0.4·|x - x*|²/d²·(x - x*) + noise
```

where α = `basin_strength`, x* = target attractor, d = `n_dims`.

The quartic term creates a genuine basin structure: far-away ICs
experience stronger restoring forces. Convergence is NOT just exponential
decay — trajectories exhibit nonlinear relaxation with distance-dependent
timescales.

Parameters: `n_dims=10`, `basin_strength=2.0`, `ic_spread=5.0`,
`n_ics=20`, `n_steps=200`, `dt=0.05`.

ICs are sampled uniformly in a ball of radius `ic_spread` centered on
the target, ensuring diverse initial conditions spanning the basin.

### MultiBasinGRN (`epc/models/canalization.py`)

Continuous-valued gene regulatory network with sigmoidal activation:

```
dx_i/dt = -decay·x_i + σ(W·x + bias)
```

where W is a Hebbian weight matrix constructed from a target pattern,
and bias = 3.0 × signed_pattern ensures a deep dominant basin.
State variables x_i ∈ [0, 1] represent gene expression levels.

This provides a biologically-motivated T1b cross-model for P25,
distinct from the gradient-flow landscape.

### Negative controls

- **DiffusiveDynamics**: pure random walk (dx = noise·dW). ICs diverge;
  final-state variance ≥ IC variance. Detector correctly rejects.
- **TrivialCollapse**: instant map (x → constant at step 0). All "finals"
  identical but relaxation time = 0. Detector rejects via the
  trivial-collapse gate (requires relaxation_time > 1).

## Detector design

### T1a observation contract

Input: canalization observation bundle — trajectories from multiple ICs
with keys: `state`, `step`, `trial`, `ic`, `target`,
`distance_to_target`, `converged`.

The adapter `extract_observation_bundle()` groups records by trial and
extracts per-trial IC, final state, convergence status, and
steps-to-convergence.

### Primary metric

**Convergence variance ratio** = Var(finals) / Var(ICs), computed as
the ratio of total variances (sum of per-dimension variances) across
all trials. Canalized systems have ratio ≪ 1.

### Secondary metrics

- **Basin volume**: fraction of ICs that converge (distance < 0.1).
- **Relaxation time**: median steps-to-convergence across converged trials.
- **Restoration accuracy** (if perturbation applied): fraction of
  perturbed trajectories that return to target.

### Null model

IC-distribution surrogate: generate "final states" by sampling from
the IC distribution (multivariate Gaussian with IC mean and covariance).
If the system is NOT canalized, finals should be as spread as ICs,
giving ratio ≈ 1.0. The p-value is the fraction of surrogate ratios
≤ observed ratio.

### Tier criteria

| Tier | Criteria |
|------|----------|
| Screening | convergence_variance_ratio < 0.1 |
| Confirmation | ratio < 0.1 AND basin_volume ≥ 0.8 AND null p < 0.01 |
| Definitive | confirmation AND null p ≤ 0.005 AND Cohen's d > 1.0 AND relaxation_time > 1 AND (if perturbation: restoration ≥ 0.8) |

### Trivial-collapse gate

The detector explicitly rejects systems with relaxation_time ≤ 1 step.
This prevents false positives from trivial constant maps (x → constant)
where convergence ratio is 0 but no dynamics-driven convergence occurs.

### Exclusion checks

- **P16**: not applicable (P16 = multiple retrievable targets, P25 = single target)
- **P24**: not applicable (P24 = cybernetic regulation under perturbation, P25 = geometric IC convergence)

## Reproduction results

### dim1: Convergence variance ratio (Sprint 81)

CanalizedLandscape at canonical parameters (n_dims=10, basin_strength=2.0,
ic_spread=5.0, 20 ICs, 200 steps):
- IC variance: 20.25
- Final variance: ≈ 0 (machine precision)
- Convergence ratio: ≈ 0 (< 0.10 tolerance PASS)
- Basin volume: 1.0 (≥ 0.80 tolerance PASS)
- Detector tier: DEFINITIVE (confidence 0.90)

### dim2: 20-seed campaign (Sprint 81)

20 seeds (0–19), same parameters:
- Convergence ratio: 0.000 ± 0.000 (CV ≈ 0%)
- Basin volume: 1.000 ± 0.000
- All 20 seeds DEFINITIVE (confidence 0.90)
- Detection fraction: 100%

### T1b: MultiBasinGRN cross-model (Sprint 81)

MultiBasinGRN (n_genes=10, n_ics=20, n_steps=400, seed=42):
- Convergence ratio: ≈ 0 (< 0.10 PASS)
- Basin volume: 1.0
- Detector tier: DEFINITIVE (confidence 0.90)

Both models exhibit equifinality and are correctly detected. Both
negative controls (diffusive, trivial collapse) are correctly rejected.
