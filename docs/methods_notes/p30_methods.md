# P30 Methods Note — Spontaneous Boundary Formation (Autopoiesis)

## Pattern definition

P30 detects spontaneous boundary formation: agents organize into a closed,
semi-permeable topology that maintains an internal micro-environment distinct
from the exterior. Distinct from P1 (similarity-driven aggregation) in that P30
requires topological closure with a functional inside/outside distinction, not
merely clustering. Distinct from P4 (territorial exclusion) in that P4 creates
exclusive domains between agents; P30 creates one self-maintaining boundary.

## Model: SCL-style autopoiesis (`epc/models/autopoiesis.py`)

Three particle types in a 2D periodic domain:

| Type | Label | Role |
|------|-------|------|
| 0 | Substrate (S) | Free-diffusing resource |
| 1 | Catalyst (C) | Converts nearby S → L; near-stationary |
| 2 | Link (L) | Membrane particle; attracted to equilibrium radius around C |

### Update rules (per timestep)

1. **Production:** For each substrate particle within `production_radius` of any
   catalyst, convert to link with probability `production_rate`.

2. **Decay:** Each link particle reverts to substrate with probability
   `decay_rate`.

3. **Forces:**
   - Catalyst-link radial spring: links attracted toward
     `membrane_equilibrium_radius` from nearest catalyst
     (spring constant = `catalyst_link_attraction`).
   - Link-link tangential attraction: nearby links attract
     (strength = `link_attraction`; cohesion for chain formation).
   - Short-range steric repulsion: all pairs within `repulsion_radius`.

4. **Diffusion:** Each particle type has its own diffusion coefficient:
   substrate (high) > link (low) > catalyst (very low).

5. **Periodic boundary conditions.**

### Parameter defaults and justifications

| Parameter | Default | Justification |
|-----------|---------|---------------|
| n_substrate | 100 | Sufficient substrate pool for membrane maintenance |
| n_catalyst | 3 | Small catalyst cluster; fewer = cleaner signal |
| box_size | 20.0 | Large enough for clear interior/exterior separation |
| production_rate | 0.15 | Moderate conversion; avoids overwhelming link count |
| decay_rate | 0.01 | Low decay → stable membrane; production/decay ratio ≈ 15 |
| production_radius | 3.0 | Zone around catalyst where membrane forms |
| membrane_equilibrium_radius | 3.0 | Preferred membrane distance from catalyst |
| catalyst_link_attraction | 0.5 | Keeps membrane at preferred radius |
| link_attraction | 0.3 | Chain cohesion between neighboring links |
| substrate_diffusion | 0.3 | High: substrate flows freely |
| link_diffusion | 0.01 | Low: membrane particles stay in place |
| catalyst_diffusion | 0.002 | Near-stationary production center |

### Negative controls

1. **NonBondingParticleModel:** All particles are type 0 (substrate), Brownian
   diffusion only. No production, no bonds → closure_fraction = 0.

2. **DenseClusterModel:** All particles attract each other (P1-like aggregation).
   No type differentiation, no production zone → association_score ≈ 0.

## Detection metrics

### Primary: association_score

**Definition:** Fraction of link particles within `association_radius` of any
catalyst, divided by the expected fraction under complete spatial randomness
(CSR).

```
frac_observed = n_links_near_catalyst / n_links_total
frac_CSR = min(n_cat × π × r² / L²)
association_score = frac_observed / frac_CSR
```

For autopoiesis: score >> 1 (links concentrate near catalysts).
For random: score ≈ 1.

### Secondary: closure_fraction

Angular coverage test per catalyst: link particles around catalyst have maximum
angular gap < π → catalyst is enclosed. Fraction of enclosed catalysts.

### Secondary: enrichment_ratio

Catalyst density inside the 90th-percentile link radius vs uniform expectation.
Measures catalyst confinement by the membrane.

## Null model

**Type-shuffle permutation:** Keep all positions fixed, randomly permute type
labels (S/C/L) among particles. Test statistic: association_score.

Under null: "catalysts" are at random positions among the particle cloud →
association_score ≈ 1.0.
Under autopoiesis: links concentrate around real catalyst positions →
association_score >> 1.

p-value: fraction of null permutations with association_score ≥ observed.

## Three-tier gates

| Tier | Gates |
|------|-------|
| Screening | association_score > 1.5 AND closure_fraction > 0.5 AND n_links ≥ 3 |
| Confirmation | + enrichment_ratio > 1.2 + null p < 0.01 + persistence > 0.5 |
| Definitive | + closure > 0.7 + enrichment > 2.0 + persistence > 0.8 + link_cv < 0.3 |

## Observation bundle (T1a)

| Key | Type | Description |
|-----|------|-------------|
| `positions` | `list[ndarray(N, 2)]` | Particle positions at each snapshot |
| `types` | `list[ndarray(N,)]` int | Particle type (0=S, 1=C, 2=L) |
| `bonds` | `list[list[tuple]]` | Bonded link pairs (optional) |
| `steps` | `ndarray(T,)` int | Step numbers |
| `box_size` | `float` | Domain side length |
| `n_particles` | `int` | Total particle count |

## Limitations

1. **Topological closure metric:** The angular-coverage test (max gap < π) is a
   necessary but not sufficient condition for topological closure. With many link
   particles, angular coverage can be trivially satisfied. The association_score
   metric (spatial co-location) provides the primary discrimination.

2. **Self-repair:** Self-repair after breach is mechanistically present (production
   regenerates lost links) but the recovery rate depends on decay/production
   balance. Full recovery may require extended simulation time beyond the standard
   detection window.

3. **Escape clause (Sprint 92 brief):** P30 is acknowledged as the hardest pattern.
   Closure + gradient reproduce reliably at CONFIRMATION tier. Self-repair is
   partial and is not required for the CONFIRMATION gate.

## Multi-seed characterization (dim2)

See `analysis/outputs/p30_multiseed.json` for 20-seed campaign results.
