# P24 Methods Note — Homeostatic Regulation

**Pattern:** P24 — Homeostatic regulation
**Canonical model:** Proportional negative-feedback homeostat (Ashby 1956)
**Detector:** `epc/detectors/p24_homeostasis.py`
**Model:** `epc/models/homeostasis.py`
**Reproduction sprint:** Sprint 72 (`analysis/reproductions/p24_homeostasis.py`)

---

## 1. Pattern and canonical reference

P24 is the maintenance of an internal variable near a set-point despite
sustained external perturbation, through active corrective (negative)
feedback. The defining signatures are:

1. **Bounded deviation**: the regulated variable does not drift under
   sustained perturbation (deviation growth ratio ≈ 1.0).
2. **Small deviation integral**: ∫|x − setpoint| dt is orders of magnitude
   smaller than what an uncontrolled system would produce.
3. **Active feedback**: the bounded deviation is caused by a corrective
   mechanism, not passive damping or equilibrium convergence.

**Distinctness from P25 (equifinality):** P24 is cybernetic — continuous
active regulation maintaining a set-point under perturbation. P25 is
geometric — convergence from arbitrary initial conditions to an attractor.
P24 requires ongoing perturbation and ongoing correction; P25 requires
diverse ICs and convergence to the same macrostate.

**Primary reference:** Ashby, W. R. (1956). *An Introduction to Cybernetics*.
Chapman & Hall.

**Secondary references:**
- Wiener, N. (1948). *Cybernetics*. MIT Press.
- Cannon, W. B. (1932). *The Wisdom of the Body*. W. W. Norton.

---

## 2. Model implementation

### 2.1 Dynamics

Two controller variants are provided:

**Proportional controller (canonical):**
```
dx/dt = -gain × (x - setpoint) + perturbation(t) + noise
```
Steady-state deviation under sustained perturbation: x_ss = setpoint + pert/gain.

**Integral controller (T1b cross-model):**
```
dx/dt = perturbation(t) + control(t) + noise
control(t) = -gain × ∫(x - setpoint) dt
```
The integral controller accumulates error and achieves zero steady-state error
for step perturbations (exact regulation), providing an independent implementation
for cross-model generalization testing.

### 2.2 Perturbation schedule

Perturbations are parameterized by onset time, optional offset time, and amplitude.
Sustained perturbation (no offset) is the primary test case; pulse perturbation
(onset + offset) tests recovery.

### 2.3 Integration

Euler method with timestep dt. Process noise is Gaussian with scale noise_std × √dt.
Deterministic (noise_std=0) for dim1 reproduction; stochastic for dim2 multi-seed.

### 2.4 Parameters

| Parameter | Canonical value | Role |
|-----------|----------------|------|
| setpoint  | 10.0           | Target regulated value |
| gain      | 5.0            | Feedback strength (0 = no regulation) |
| dt        | 0.1            | Integration timestep |
| amplitude | 5.0            | Perturbation magnitude |
| onset     | 50.0           | Perturbation start time |
| noise_std | 0.0 (dim1) / 0.5 (dim2) | Process noise |
| n_steps   | 2000           | Simulation length |

---

## 3. Detector design

### 3.1 Observation contract (T1a)

The detector reads a **scalar-regulated-variable observation bundle** via
`extract_observation_bundle()`. The bundle contains four aligned time series:
`time`, `x` (regulated variable), `setpoint`, and `perturbation`. Any system
producing these four fields can be tested by the P24 detector without
modification.

### 3.2 Null model

The null model is a **surrogate uncontrolled trajectory**: integrate
`dx = perturbation(t) × dt + noise` starting from the observed x at
perturbation onset, with no corrective feedback. Noise scale is estimated
from pre-onset trajectory residuals. This tests whether bounded deviation
could arise without active regulation.

For a deterministic system with sustained perturbation, all surrogates
produce linearly growing deviation, giving deviation integrals 2-3 orders
of magnitude larger than the regulated system.

### 3.3 Tier criteria

- **Screening:** Deviation growth ratio ≤ 2.0 (deviation is bounded,
  not growing without bound).
- **Confirmation:** Deviation integral significantly less than surrogate
  null (p < 0.01) AND deviation ratio < 0.5 vs null mean.
- **Definitive:** Deviation ratio < 0.3 AND Cohen's d > 1.0 AND p ≤ 0.005
  AND metadata (if available) confirms active feedback.

### 3.4 Key design decisions

1. **Surrogate null over shuffle null:** A time-shuffle null is ineffective
   for regulated systems because the marginal distribution of x values is
   narrow (all near setpoint + pert/gain), so shuffled trajectories have
   similar deviation integrals. The surrogate uncontrolled null explicitly
   tests "what if there were no feedback?" — the core homeostasis claim.

2. **Growth ratio screening:** Rather than checking recovery time (which
   requires pulse perturbation), the growth ratio test works for both
   sustained and pulse perturbation: regulated systems show ratio ≈ 1.0,
   uncontrolled show ratio >> 1.0.

3. **Metadata gate at definitive:** If model metadata explicitly indicates
   `has_active_feedback=False`, definitive tier is denied even if the
   statistical criteria pass. This prevents false positives from passive
   systems that happen to have bounded deviation (e.g., friction-damped).

---

## 4. Reproduction results

### 4.1 dim1 (Sprint 72)

Proportional controller (gain=5.0) under sustained perturbation (amplitude=5.0):
- Controlled deviation integral: 149.75
- Uncontrolled deviation integral: 56,175.03
- **Deviation ratio: 0.0027** (tolerance: < 0.30) — **PASS**
- Detector tier: **DEFINITIVE** (confidence=0.90) — **PASS**

### 4.2 dim2 (Sprint 72)

20-seed campaign with process noise (noise_std=0.5):
- Deviation integral: 149.59 ± 1.20 (CV = 0.8%)
- Deviation ratio: 0.0027 ± 0.00003 (CV = 1.2%)
- Steady-state deviation: 0.999 ± 0.015
- All 20 seeds: DEFINITIVE tier

### 4.3 T1b cross-model (Sprint 72)

Integral controller (gain=2.0) tested with same detector:
- DEFINITIVE detection (p ≤ 0.005)
- Confirms detector recognizes the *phenomenon* (homeostatic regulation),
  not the specific implementation (proportional vs integral control).
