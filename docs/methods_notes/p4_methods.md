# P4 Methods Note — Territoriality / Exclusion Boundaries

## Pattern definition

Mobile agents using transient environmental marks create stable, mutually
exclusive spatial domains with persistent boundaries. Outcome is spatial
EXCLUSION + boundary maintenance, not clustering (distinct from P1
aggregation) and not connective transport networks (distinct from P29
stigmergic transport).

## Canonical reference

Giuggioli, L., Potts, J. R., & Harris, S. (2011). Animal interactions
and the emergence of territoriality. *PLoS Computational Biology*, 7(3),
e1002008.

## Model implementations

### ScentMarkingModel (canonical)

N agents on an L×L torus. Each agent deposits own-ID scent (a decaying
scalar field); scent decays exponentially at rate γ per step. Movement
uses a two-stage rule:

1. **Avoidance (Stage 1):** Among von Neumann neighbors + stay, exclude
   cells where total foreign scent exceeds `repulsion_strength × deposition_rate`.
   This creates hard territorial boundaries where active foreign marking
   is present, while stale marks (decayed below threshold) become available
   for territory expansion.

2. **Home attraction (Stage 2):** Among allowed cells, select with
   probability ∝ (1 + home_attraction × own_scent)^(1/temperature),
   biasing movement toward familiar territory.

If all neighbors are blocked (boundary deadlock), agent stays in place.

Default parameters: N=4, L=48, deposition=0.1, decay=0.03,
repulsion_strength=2.0, home_attraction=2.0, temperature=0.5, steps=20000.

### PheromoneRepulsionModel (T1b cross-model)

Independent implementation with simplified mechanics: agents deposit
pheromone and avoid cells where foreign pheromone exceeds a hard
threshold. Among allowed cells, selection is uniform random (no
home-attraction bias). Provides an independent implementation mechanism
(hard threshold vs softmax) for T1b cross-model generalization.

Default parameters: N=4, L=48, deposition=5.0, decay=0.005,
avoidance_threshold=0.3, steps=5000.

### PlainRandomWalkModel (negative control)

Agents deposit scent identically but do not respond to it. Movement is
uniform random among von Neumann neighbors + stay. Used as the
non-territorial baseline: home ranges overlap freely.

## Detector design

### T1a observation bundle

The detector reads inputs via `extract_observation_bundle(history)` which
extracts:
- `positions`: list of ndarray(N, 2)
- `scent_fields`: list of ndarray(N, L, L)
- `occupancy`: list of ndarray(N, L, L) — cumulative visit counts
- `steps`, `n_agents`, `grid_size`

Any system emitting this schema can be consumed by P4.

### Primary metrics

1. **Exclusivity index:** For each visited cell, fraction of total visits
   from the dominant agent. Averaged over all visited cells. Perfect
   territories → 1.0; N agents sharing uniformly → 1/N.

2. **Pairwise home-range overlap:** Sum of min(normed_occ_i, normed_occ_j)
   for each agent pair. Perfect exclusion → ~0; complete overlap → ~1.

3. **Boundary persistence:** Assign each cell to the dominant agent at
   two time windows (early=quarter, late=final). Fraction of cells with
   same owner = persistence. Stable territories → ~1.0.

4. **Occupancy-scent correlation:** Pearson r between agent's own
   occupancy and total foreign scent. Negative = scent-mediated exclusion.

### Null model

Cell-level multinomial shuffle. For each visited cell with T total visits,
draw agent assignments from Multinomial(T, [1/N, ..., 1/N]) and recompute
exclusivity. This tests whether per-cell agent dominance exceeds what
random (non-territorial) assignment would produce. 199 permutations,
one-sided p-value.

### Tier criteria

| Tier | Requirements |
|------|-------------|
| Screening | Exclusivity > 0.80 AND p < 0.10 |
| Confirmation | + boundary persistence > 0.6, p < 0.05, Cohen's d > 1.0 |
| Definitive | + occ-scent correlation < −0.005, p ≤ 0.005, metadata confirms foreign avoidance |

### Confidence computation

Uses `compute_confidence(tier, bonuses)` from `epc/detector_result.py`:
- Screening: base 0.35, bonuses for secondaries_pass, shuffle_null_p_001
- Confirmation: base 0.55, bonuses for null_p_0001, effect_size_gt_1, all_secondaries
- Definitive: base 0.75, bonuses for all_exclusions_cleared, both_null_types_rejected, finite_size_robustness

## Reproduction results (Sprint 86)

Single-seed (seed=42) canonical ScentMarkingModel:
- Exclusivity = 0.902 (tolerance: > 0.85)
- Pairwise overlap = 0.034 (tolerance: < 0.10)
- Boundary persistence = 0.865 (tolerance: > 0.70)
- Detector tier = definitive (confidence = 0.90)
- Cohen's d = 157.5
- All tolerances PASS
