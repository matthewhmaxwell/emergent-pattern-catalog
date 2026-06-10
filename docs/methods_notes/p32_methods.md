# P32 Emergent Specialization / Division of Labor — Methods Note

## Pattern and canonical reference

**P32 — Emergent specialization (division of labor).** Initially IDENTICAL
agents spontaneously differentiate into distinct functional roles that
increase collective efficiency, without pre-assigned task allocation.
Canonical minimal model: the **response-threshold model** (Bonabeau,
Theraulaz & Deneubourg, 1996).

**Distinctness from P23:** P32 is role-FIXED specialization — agents settle
into stable task identities with low switching. P23 (anti-coordination) is
continuous re-balancing where agents keep changing choices. The P32 detector
excludes P23 via low late-window task-switching frequency (< 0.3).

**Distinctness from P19:** P32 yields stable task roles across the
population; P19 leadership is a transient, individual-level phenomenon.

## Model update rules and state encoding

### ResponseThresholdModel (canonical)

**State:** N agents, each with a vector of per-task response thresholds
θ_{i,j} ∈ [0.01, 1.0]. M task types, each with stimulus s_j ∈ [0, 1].

**Response probability (Bonabeau et al. 1996, Eq. 1):**

    P(agent i responds to task j) = s_j² / (s_j² + θ_{i,j}²)

Each agent evaluates all tasks independently. If it responds to multiple
tasks, it performs the one with highest response probability.

**Threshold reinforcement:**
- Performing task j: θ_{i,j} ← max(0.01, θ_{i,j} - ξ)
- Not performing task j: θ_{i,j} ← min(1.0, θ_{i,j} + φ)

where ξ is the reinforcement rate and φ is the forgetting rate.

**Stimulus dynamics:** Stimulus accumulates at a base rate and decreases
proportional to the number of workers:

    s_j ← clip(s_j + r - 0.1 · workers_j, 0, 1)

### NoReinforcementModel (negative control)

Identical dynamics but with ξ = 0, φ = 0. Thresholds never change,
so agents remain generalists.

### Multiplicative-threshold variant (T1b cross-model)

Same architecture but with multiplicative threshold update:
- Performing: θ_{i,j} ← max(0.01, θ_{i,j} × (1 - rate))
- Not performing: θ_{i,j} ← min(1.0, θ_{i,j} × (1 + rate))

Structurally different reinforcement scheme verifying OOD-readiness.

## Parameter defaults and justifications

| Parameter | Default | Justification |
|-----------|---------|---------------|
| n_agents | 20 | Sufficient for role differentiation across 3 tasks |
| n_tasks | 3 | Minimum for non-trivial division (Bonabeau 1996 uses 2–5) |
| reinforcement_rate (ξ) | 0.05 | Moderate learning rate; 0.15 for faster convergence in reproduction |
| forgetting_rate (φ) | 0.01 | Slow threshold drift for non-performed tasks |
| stimulus_rate | 0.1 | Base demand accumulation |
| initial_threshold | 0.5 | All agents start identical (symmetric initial conditions) |
| n_steps | 500–1000 | Sufficient for threshold convergence at ξ=0.05 |

## Detection metrics

### Primary: entropy decline

Per-agent task entropy computed in early and late windows. Agents
specializing will show declining entropy (from near-maximal toward 0).

    H_i = -Σ_j p_{i,j} log₂(p_{i,j})

where p_{i,j} is the fraction of time agent i spends on task j within
the window. Screening threshold: decline > 0.1 bits.

### Secondary: role diversity

Number of distinct dominant tasks across agents, normalized by n_tasks.
Must be ≥ 0.5 (at least half the tasks have dedicated specialists) and
NOT single-task collapse (all agents on the same task).

### Secondary: efficiency

Fraction of timesteps where all tasks have at least one worker assigned.
Late-window efficiency should exceed early-window efficiency (specialization
improves collective performance).

### Secondary: switching frequency

Mean per-agent task-switching frequency: fraction of consecutive timestep
pairs where an active agent changes tasks. Low late switching (< 0.3)
indicates stable specialization, excluding P23 anti-coordination.

## Null model

**Time-shuffle null (NullType.SHUFFLE):** For each agent independently,
randomly permute its assignment time series. This preserves the marginal
task-frequency distribution but destroys any temporal trend in
specialization. Under the null, early and late entropy should be
statistically identical.

p-value: fraction of null entropy-decline values ≥ observed decline.
Floor: 1/(n_permutations + 1).

## Three-tier gates

| Tier | Gates |
|------|-------|
| Screening | entropy_decline > 0.1 bits |
| Confirmation | null p < 0.01 AND role_diversity ≥ 0.5 AND NOT single_task_collapse AND late_coverage ≥ 0.8 AND late_entropy < 1.0 |
| Definitive | null p < 0.005 AND entropy_decline > 0.3 AND (late_efficiency > 0.3 OR efficiency_gain > 0) AND role_diversity ≥ 0.67 |

## Observation bundle (T1a)

| Key | Type | Description |
|-----|------|-------------|
| `task_assignments` | `ndarray(T, N)` int | Task index per agent per step (-1 = idle) |
| `n_agents` | `int` | Number of agents |
| `n_tasks` | `int` | Number of task types |
| `steps` | `ndarray(T,)` int | Step numbers |

## Multi-seed characterization (dim2)

20 seeds at ξ=0.05, φ=0.01, n_steps=1000. Final per-agent entropy
mean ≈ 1.21 ± 0.21 (CV = 0.17). Efficiency gain mean ≈ 0.04 ± 0.04.
4/20 seeds reach confirmation; specialization strength varies by seed
due to stochastic task encounters.
