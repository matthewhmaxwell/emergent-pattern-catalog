# P16 Associative Memory / Pattern Completion — Methods Note

## Pattern and canonical reference

**P16 — Associative memory / pattern completion.** A system stores and
retrieves distributed patterns via attractor dynamics, completing partial
inputs to stored templates. The canonical minimal model is the **Hopfield
network** (Hopfield, 1982); the detector also generalizes to Boolean
gene-regulatory networks with multiple fixed-point attractors (T1b).

The storage capacity transition at α_c ≈ 0.138 was established by
Amit, Gutfreund, & Sompolinsky (1985).

## Model update rules and state encoding

### Hopfield network

**State:** N binary neurons s_i ∈ {−1, +1}.

**Pattern storage (Hebbian learning):**
Given P random patterns {ξ^μ}, μ = 1..P, each ξ^μ_i ∈ {−1, +1}:

    w_ij = (1/N) Σ_μ ξ^μ_i ξ^μ_j,    w_ii = 0

**Dynamics:** Asynchronous sign-update. At each sweep, neurons are updated
in random order:

    h_i = Σ_j w_ij s_j
    s_i → sign(h_i)    (unchanged if h_i = 0)

Convergence: the network reaches a fixed point (no further state changes)
within O(N) sweeps for typical initial conditions.

**Retrieval:** Start from a corrupted cue (stored pattern with fraction
`corruption` of bits flipped). If the network converges to a state with
high overlap m = (1/N) Σ_i ξ^μ_i s_i with the target pattern, retrieval
is successful.

### Boolean GRN (T1b cross-model)

Same binary-state, Hebbian-weight, asynchronous-update structure but
interpreted as a gene-regulatory network where gene expression states
(on/off) settle into designed attractor configurations. Structurally
identical dynamics ensure the T1a adapter reads both models identically.

## Parameter defaults and justifications

| Parameter | Default | Justification |
|-----------|---------|---------------|
| N | 100 | Sufficient for clean retrieval at low load; N=500 for capacity reproduction |
| P | 5 | α = 0.05, well below α_c ≈ 0.138 |
| corruption | 0.2 | 20% bit-flip: standard retrieval challenge |
| max_steps | 100 | Convergence typically in 1–3 sweeps at low load |
| seed | 42 | Reproducibility |

## Canonical positive parameters

Hopfield N=100, P=5, corruption=0.2, seed=42: all 5 cue patterns retrieved
with overlap > 0.95 and selective recall accuracy = 1.0. Detector reaches
DEFINITIVE (confidence=0.90).

## Detection pipeline

### T1a observation contract

The detector reads an **attractor-network state-vector trajectory** via
`extract_observation_bundle()`:

- `trials`: list of per-trial dicts, each with `states` (list of state vectors),
  `overlaps`, `cue_pattern_idx`, `converged`
- `stored_patterns`: (P, N) array of the stored templates
- `N`, `P`: network size and pattern count

Any system producing history dicts with keys `state`, `step`, `trial`,
`cue_pattern_idx`, `overlap`, `stored_patterns`, `converged` is compatible.

### Three-tier detection

**SCREENING** (confidence cap 0.60):
- Mean completion accuracy (best overlap with any stored pattern, averaged
  across trials) exceeds max(0.5, chance_overlap + 0.2) where
  chance_overlap ≈ 2/√N by CLT.

**CONFIRMATION** (confidence cap 0.85):
- Completion accuracy > 0.7
- Selective recall accuracy > 0.5 (more than half the cues retrieve the
  correct target pattern)
- Random-weights null p < 0.01

**DEFINITIVE** (confidence cap 1.00):
- All confirmation criteria
- n_distinct patterns selectively recalled ≥ 2 (genuine content-addressable
  memory, not convergence to a single attractor)
- Recall accuracy ≥ 0.8
- Completion accuracy > 0.8
- Null p < 0.01
- Metadata check: `has_multiple_attractors` must not be False

## Null model specification

**Random-weights surrogate.** For each null run:
1. Generate P random ±1 patterns (independent of the stored patterns)
2. Compute Hebbian weights from these random patterns
3. Attempt to retrieve the ORIGINAL stored patterns from corrupted cues
   using the random weights
4. Record mean best-overlap (completion accuracy)

The null tests whether the observed completion accuracy could arise from
a network with no relationship between its weights and the target patterns.
At N=100, the null distribution has mean ≈ 0.17, std ≈ 0.065 — cleanly
separated from the observed ≈ 0.96 (Cohen's d ≈ 12).

**Rare null outlier note:** With probability ≈ 1/200, a random-weights
network produces a spurious fixed point that happens to align well with
a stored pattern. This is a genuine rare event (not a null model failure);
the definitive tier uses p < 0.01 (not ≤ 0.005) to accommodate this.

## Cross-detection matrix outcomes

| Model × P16 | Expected | Rationale |
|---|---|---|
| hopfield (canonical) | DEFINITIVE | Low-load regime, perfect selective recall |
| boolean_grn (T1b) | DEFINITIVE | Same attractor structure, independent implementation |

## dim1 reproduction: AGS1985 storage capacity

Amit, Gutfreund & Sompolinsky (1985): at N=500, sweep α ∈ [0.02, 0.25]:
- α < 0.10: retrieval overlap ≈ 1.0 (PASS)
- Transition midpoint: α ≈ 0.173 (finite-size shifted from N→∞ value 0.138;
  tolerance [0.10, 0.20] PASS)
- α = 0.20: overlap ≈ 0.37 (< 0.50 threshold, PASS)

All three tolerance checks pass. The finite-size shift is a documented
property of the Hopfield model (the AGS result is exact only at N→∞).

## dim2 multiseed

20 seeds at α = 0.05 (N=100, P=5): completion accuracy = 1.000 ± 0.000
(CV = 0.0%). All 20 seeds reach DEFINITIVE. The zero variance reflects
the deep sub-capacity regime where retrieval is essentially deterministic.
