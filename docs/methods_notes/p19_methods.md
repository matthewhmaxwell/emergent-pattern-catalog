# P19 Methods Note — Emergent Leadership / Minority Guidance

**Pattern:** P19 — Emergent leadership / minority guidance
**Canonical model:** Informed minority flocking (Couzin et al. 2005)
**Detector:** `epc/detectors/p19_emergent_leadership.py`
**Model:** `epc/models/informed_minority.py`
**Reproduction sprint:** Sprint 69 (`analysis/reproductions/p19_couzin2005.py`)

---

## 1. Pattern and canonical reference

P19 is the emergence of group-level directional decision-making driven by a
small informed minority. A subset of agents (fraction ρ) has a weak bias
toward a preferred direction, while the majority aligns via local social
interactions (Vicsek-style alignment). The whole group converges to move
in the preferred direction, even when the informed fraction is very small
(ρ ~ 0.05). The defining signatures are:

1. **Accuracy increases with ρ** (diminishing returns curve)
2. **Disproportionate influence**: small ρ achieves high group accuracy
3. **Asymmetric influence**: informed agents' headings are closer to the
   preferred direction than naive agents', detectable via label-shuffle null

**Primary reference:** Couzin, I. D., Krause, J., Franks, N. R., & Levin, S. A.
(2005). Effective leadership and decision-making in animal groups on the move.
*Nature* 433(7025), 513–516.

**Secondary references:**
- Couzin, I. D. (2009). Collective cognition in animal groups. *Trends Cogn.
  Sci.* 13(1), 36–43.
- Berdahl, A. et al. (2013). Emergent sensing of complex environments by
  mobile animal groups. *Science* 339(6119), 574–576. (Related: P17.)

---

## 2. Model implementation

### 2.1 Dynamics

The model is a Vicsek-style SPP system with an informed minority.
Synchronous update on N agents in a periodic [0, L)² box:

**Naive agents (fraction 1-ρ):**
$$θ_i(t+1) = \text{arg}[\langle e^{iθ_j} \rangle_{j \in N_i}] + η ξ_i$$

**Informed agents (fraction ρ):**
$$θ_i(t+1) = \text{arg}[(1-ω)\langle e^{iθ_j} \rangle_{j \in N_i} + ω e^{iθ_{\text{pref}}}] + η ξ_i$$

where $ω$ is the bias weight (strength of preference for $θ_{\text{pref}}$),
$N_i$ is the set of neighbors within metric radius $r$ (including self),
$η$ is the angular noise strength, and $ξ_i \sim U[-1/2, 1/2]$.

Positions update as:
$$x_i(t+1) = x_i(t) + v_0 \cos(θ_i(t+1)) Δt$$

### 2.2 Canonical parameters

| Parameter | Symbol | Canonical value | Notes |
|-----------|--------|-----------------|-------|
| N agents | N | 200 | Couzin 2005 |
| Box size | L | 10.0 | ρ_area = 2.0 |
| Speed | v₀ | 0.03 | |
| Noise | η | 0.1 | Low noise → strong alignment |
| Interaction radius | r | 1.0 | Metric, including self |
| Informed fraction | ρ | 0.1 | dim1 sweep variable |
| Bias weight | ω | 0.3 | Moderate bias |
| Preferred direction | θ_pref | 0.0 | East |

### 2.3 Informed agent assignment

The first `n_informed = round(N × ρ)` agents are marked as informed.
Assignment is deterministic from the initial setup (not shuffled each step).
This matches Couzin (2005) where informed individuals maintain their
knowledge throughout the trial.

---

## 3. Detector design

### 3.1 Three-tier detection

**Screening (confidence ≤ 0.60):**
Group mean heading aligns with the preferred direction with accuracy > 0.3.
Accuracy = cos(θ_group - θ_pref), where θ_group is the circular mean of
all headings. Measured over the last 50% of history (steady state).

**Confirmation (confidence ≤ 0.85):**
Influence asymmetry: informed agents' mean heading is closer to the preferred
direction than naive agents' mean heading, tested via label-shuffle permutation
null. The "pull" metric is:

$$\text{pull}(t) = \cos(\bar{θ}_{\text{inf}}(t) - θ_{\text{pref}}) - \cos(\bar{θ}_{\text{naive}}(t) - θ_{\text{pref}})$$

where $\bar{θ}$ denotes the circular mean. Mean pull > 0 with p < 0.05
(label-shuffle null, ≥ 99 permutations) confirms that the informed agents
are specifically leading toward the preferred direction.

**Definitive (confidence ≤ 1.00):**
Guidance efficacy = accuracy / ρ > 2.0 (disproportionate influence) AND
ρ ≤ 0.25 (genuine minority) AND pull p < 0.01 (≥ 199 permutations).

### 3.2 Null model

Label-shuffle permutation: randomly reassign which agents are "informed"
(preserving the count) and recompute the pull. Under the null hypothesis
(no informed-naive distinction matters), shuffled labels produce the same
expected pull. This is more direct than a time-shuffle null because it
tests whether the *specific agents labeled as informed* are systematically
biased toward the preferred direction.

### 3.3 Architecture decision: directional cross-correlation over TE

The brief specifies "transfer entropy (or directional cross-correlation)"
for the influence-asymmetry metric. In practice, TE on mean-heading time
series produces zero signal in the converged Vicsek regime because both
informed and naive group means are constant (the KSG estimator requires
temporal variability). The label-shuffle directional pull metric captures
the same information asymmetry more robustly: it detects that specific
agents (the informed ones) are consistently biased toward the preferred
direction, which is the mechanistic signature of minority guidance.

---

## 4. Distinctness from related patterns

| Pattern | Mechanism | How P19 differs |
|---------|-----------|-----------------|
| P5 (flocking) | Symmetric alignment → polarization | P5 has no preferred direction; P19 requires asymmetric influence |
| P17 (collective sensing) | Group averages noisy field samples → gradient climbing | P17 uses environmental sensing; P19 uses social information only |
| P18 (consensus) | Symmetric pooling → agreement | P18 is symmetric; P19 is disproportionately driven by a minority |
| P32 (division of labor) | Stable specialized roles | P32 roles are persistent; P19 "leadership" is assigned externally |

P19 is distinguished from P5 by the presence of `informed_mask` and the
directional bias of informed agents. A pure Vicsek flock (ρ=0) will
polarize but not toward any particular direction — P19 requires
directional accuracy toward a specific preferred heading.

---

## 5. Reproduction notes

### 5.1 dim1: accuracy vs ρ (Couzin 2005 Fig. 2a)

Sweep over ρ ∈ {0, 0.025, 0.05, 0.1, 0.15, 0.2, 0.3, 0.5}, 10 seeds each.
Accuracy rises monotonically with ρ (Spearman ρ = 1.0, p < 0.001).
ρ=0 is chance-level (accuracy ≈ 0.13 ± 0.76). Even ρ=0.025 achieves
near-perfect accuracy (0.977 ± 0.040). This is more abrupt than Couzin
(2005) Fig. 2a because our noise (η=0.1) is low relative to the alignment
strength; the paper uses higher effective noise.

### 5.2 dim2: multi-seed (20 seeds at ρ=0.1)

All 20 seeds achieve accuracy > 0.999, confirming robust guidance at the
canonical informed fraction. CV = ~0% for accuracy. Influence pull positive
in all seeds.

### 5.3 Convergence speed

With canonical parameters (η=0.1, ρ=0.1, N=200, L=10), the group converges
to the preferred direction within ~100 steps. This fast convergence means
the steady-state accuracy is nearly 1.0, which matches the theoretical
prediction: at low noise, even a weak bias ω=0.3 in 10% of agents is
sufficient for perfect group alignment because the Vicsek alignment
mechanism amplifies any directional bias.
