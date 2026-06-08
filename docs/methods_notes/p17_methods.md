# P17 Methods Note — Distributed Sensing / Collective Gradient Detection

**Pattern:** P17 — Distributed sensing / collective gradient detection
**Canonical model:** Collective sensing (Berdahl et al. 2013)
**Detector:** `epc/detectors/p17_collective_sensing.py`
**Model:** `epc/models/collective_sensing.py`
**Reproduction sprint:** Sprint 67 (`analysis/reproductions/p17_berdahl2013.py`)

---

## 1. Pattern and canonical reference

P17 is the emergence of group-level environmental sensing that no individual
can achieve alone. A group of mobile agents in a noisy scalar field collectively
navigates toward the field optimum (gradient climbing) even though each individual's
sensing is too noisy to detect the gradient reliably. The defining signature is
that navigation accuracy (chemotactic index, CI) **rises with group size N**
while isolated agents remain at chance.

**Primary reference:** Berdahl, A., Torney, C. J., Ioannou, C. C., Faria, J. J.,
& Couzin, I. D. (2013). Emergent sensing of complex environments by mobile
animal groups. *Science* 339(6119), 574–576.

**Secondary references:**
- Camley, B. A., Zimmermann, J., Levine, H., & Rappel, W.-J. (2016). Collective
  gradient sensing and chemotaxis: modeling and recent developments. *J. Phys.:
  Condens. Matter* 28, 253001.
- Torney, C. J., Neufeld, Z., & Couzin, I. D. (2009). Context-dependent
  interaction leads to emergent search behaviour in social aggregates. *PNAS*
  106(52), 22055–22060.

---

## 2. Model equations and update rule

The implementation uses a speed-modulation mechanism in a periodic 2D domain with
a Gaussian scalar field:

```
Field:       F(x) = A * exp(-|x - x_peak|²_periodic / (2σ²))
Sensing:     q_i = F(x_i) + N(0, σ_sense)
Speed:       v_i = v_max * (1 - α * clamp(q_i / A, 0, 1))
Social:      torque_i = s * angle_to_centroid_of_neighbors
Heading:     θ_i(t+1) = θ_i(t) + N(0, σ_turn) + torque_i
Position:    x_i(t+1) = x_i(t) + v_i * [cos(θ_i), sin(θ_i)] * dt   (mod L)
```

**Mechanism:** Agents that happen to be near the field peak sense higher values
(on average, despite noise) and slow down. This causes the group's center of
mass to drift toward the peak. Social attraction toward the centroid ensures
the group stays cohesive. The key insight (Berdahl 2013): because each agent's
sensing noise is independent, the centroid-attraction mechanism effectively
averages N independent samples, reducing the effective noise by √N. Hence
CI scales with group size.

---

## 3. Parameters and canonical regime

| Parameter | Canonical value | Role |
|-----------|----------------|------|
| N | 50 | Group size (sweep target) |
| box_size | 20.0 | Periodic domain side length |
| v_max | 0.4 | Maximum speed (at zero field) |
| turn_noise | 0.3 rad | Heading noise per step |
| sensing_noise | 0.8 | Additive Gaussian noise on field sensing |
| α | 0.95 | Speed-modulation depth |
| social_strength | 0.2 | Attraction toward group centroid |
| field_sigma | 5.0 | Gaussian field width |
| field_amplitude | 1.0 | Peak field value |
| dt | 1.0 | Integration timestep |
| n_steps | 1000 | Simulation duration |
| init_mode | "offset" | Start clustered ~L/3 from peak |

**SNR analysis:** At the offset start position (~6.6 units from peak), the true
field is ~exp(-(6.6)²/(2×25)) ≈ 0.42. The gradient signal (difference from peak
to start) is ~0.58. With sensing_noise=0.8, individual SNR ≈ 0.58/0.8 = 0.73
(unreliable). For N=50, effective SNR ≈ 0.73×√50 ≈ 5.2 (strong). This is
exactly the Berdahl scaling mechanism.

---

## 4. Detection methodology

**Primary metric:** Chemotactic Index (CI) = (d_initial - d_final) / d_initial,
where d is the periodic distance from group center of mass to the field peak.
CI = 1 means perfect navigation; CI = 0 means random walk.

**Group-size sweep:** The detector runs the model at N ∈ {1, 5, 10, 25, 50}
(5 seeds each) and computes the slope of CI vs log(N).

**Null model:** α = 0 (no speed modulation). Preserves social structure and
random walks but removes the mechanism linking field quality to movement speed.
Without speed modulation, no group size produces positive CI.

**Tier criteria:**
- Screening (conf ≤ 0.60): CI(N_max) exceeds CI(N=1) + 0.05 AND CI(N_max) > 0.1
- Confirmation (conf ≤ 0.85): slope > 0.02 AND null p < 0.05
- Definitive (conf ≤ 1.00): monotonic CI-vs-N AND Cohen's d > 2.0 AND p < 0.01
  AND CI(N_max) > 0.2

---

## 5. Distinctness from P5 (flocking)

P5 detects collective motion (alignment, group polarization) without requiring
any environmental field or navigation task. P17 requires:
1. An external scalar field that agents sense
2. Speed (or persistence) modulation correlated with field value
3. Group-level navigation accuracy that scales with N

A Vicsek flock has high P5 but zero P17 (no field, no gradient climbing). A
collective-sensing group has moderate P5 (cohesive motion) AND high P17 (N-scaling
CI). The two patterns can co-occur but are independently detectable.

---

## 6. Architecture decisions

- **ADR 57:** Speed-modulation mechanism chosen over turning-rate modulation.
  Both produce the Berdahl effect, but speed modulation gives a stronger and
  more robust N-scaling signal in this parameter regime.
- **ADR 58:** α=0 null chosen over spatial-shuffle null. Breaking speed-field
  coupling directly tests the mechanism; shuffling positions would also destroy
  the social structure, making the null less informative.
- **ADR 59:** CI computed from periodic center-of-mass (circular mean) to handle
  domain wrapping correctly. Early-window (first 10%) vs late-window (last 25%)
  for robust estimation.
