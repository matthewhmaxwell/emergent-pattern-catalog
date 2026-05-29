# P22 Methods Note — Information Cascade / Social Contagion (SIR)

**Pattern:** P22 — Information cascade / social contagion
**Canonical model:** SIR epidemic cellular automaton on a 2D lattice
**Detector:** `epc/detectors/p22_information_cascade.py`
**Model:** `epc/models/sir_epidemic.py`
**Reproduction sprint:** Sprint 51 (`analysis/reproductions/p22_dattaacharyya2005.py`)

---

## 1. Pattern and canonical reference

P22 is an information cascade: a local perturbation propagates outward through
a network or medium, converting susceptible nodes into a new state. The spatial
SIR epidemic CA is the canonical implementation: a single infected cell at the
lattice center spreads infection outward in an expanding circular wavefront
that eventually burns through all susceptible cells and extinguishes.

**Primary reference:** Datta, A. & Acharyya, M. (2021/2022). Modelling the
spread of an epidemic in presence of vaccination using cellular automata.
*International Journal of Modern Physics C* 33, 2250094 (arXiv:2104.10456).

**Related references:**
- Rousseau, G., Giorgini, G., Livi, R. & Chaté, H. (1997). Dynamical phases
  in a cellular automaton model for epidemic propagation. *Physica D* 103,
  554–563.
- Grassberger, P. (1983). On the critical behavior of the general epidemic
  process and dynamical percolation. *Math. Biosci.* 63, 157–172.
- Hoya White, S., Martín del Rey, A. & Rodríguez Sánchez, G. (2007).
  Modeling epidemics using cellular automata. *Appl. Math. Comput.* 186,
  193–202.

**Note on Fuks & Lawniczak (2002):** An earlier version of this catalog cited
Fuks & Lawniczak as the primary SIR reference. Their paper describes a
Lattice Gas CA where individuals *move* between sites — a fundamentally
different model. Datta & Acharyya (2021) uses the same fixed-site SIR CA
as this implementation. The citation was corrected in Sprint 8.

---

## 2. Model update rules and state encoding

The SIR CA uses three cell states with integer encoding:

| State | Integer | Meaning |
|---|---|---|
| Susceptible | 0 | Can become infected |
| Infected | 1 | Active infection, can spread |
| Recovered | 2 | Permanently immune, no transitions |

Update is **synchronous**: all cells transition simultaneously at each step.
Two update rules apply:

**Infection (S → I):** A Susceptible cell with ≥ 1 Infected neighbour is
infected with probability:

```
P(S → I) = 1 − (1 − p)^{n_infected_neighbours}   [independent-neighbours model]
```

This is the standard independent-neighbours infection model (Datta-Acharyya
2021, Hoya White et al. 2007). Each infected neighbour independently attempts
infection; `n_infected_neighbours` infected neighbours each fail to infect
with probability `(1−p)`, so the probability of remaining susceptible is
`(1−p)^{n_infected_neighbours}`. An alternative threshold model
(`independent_neighbors=False`) uses `P(S → I) = p` if any neighbour is
infected; this is a less standard variant not used in the canonical tests.

**Recovery (I → R):** An Infected cell recovers to Recovered with probability
`q` per timestep. Recovery is independent of neighbourhood. No re-infection
is modelled in the default configuration.

**Boundary conditions:** Periodic (default) or fixed-boundary (dead border).
The canonical single-seed experiment uses periodic boundaries, so the
wavefront wraps around for very large epidemics.

---

## 3. Irreversibility enforcement (Sprint 41 prerequisite)

The defining structural feature of SIR dynamics is **irreversibility**:
transitions are strictly S(0) → I(1) → R(2). Once a cell enters the
Recovered state, it never returns to Susceptible or Infected. This is the
spatial analog of herd immunity accumulation: the recovered region behind
the wavefront acts as a firewall that prevents re-infection.

The `P22CascadeDetector.detect()` enforces irreversibility as a hard
prerequisite gate (`_check_irreversibility_prereq()`): it scans all
consecutive frame pairs and rejects if any cell shows `curr_state < prev_state`
(a state decrease, indicating backward transition). This gate was added in
Sprint 41, grounded in three published references:

1. **Datta & Acharyya (2021):** SIR CA explicitly encodes permanent immunity.
2. **Mobilia, Georgiev & Täuber (2007):** Lotka-Volterra lattice allows
   backward transitions (predator death: state 2 → 0); correctly rejected.
3. **Reichenbach, Mobilia & Frey (2007):** Spatial RPS allows species
   elimination through state 0 (empty) cycles; correctly rejected.

Before Sprint 41, LV and RPS were both false positives on P22 at SCREENING
tier (their persistent spatial activity satisfies cascade-size prerequisites).
The irreversibility gate eliminates both: LV's `PREDATOR(2) → EMPTY(0)` and
RPS's cyclic `→ EMPTY(0)` each trigger the guard immediately. After Sprint 41,
Phase-2a panel overall TNR = 1.000 (Cohen's d = +∞).

---

## 4. Parameter defaults and justifications

| Parameter | Default | Justification |
|---|---|---|
| `rows`, `cols` | 100 × 100 | Sufficient for single-seed wavefront to propagate before boundary effects; Sprint 51 uses 120×120 |
| `infection_prob` | 0.3 | Well above Moore-neighbourhood p_c ≈ 0.038 (at q=0.1); produces a spreading epidemic |
| `recovery_prob` | 0.1 | Standard canonical value across SIR CA literature |
| `neighborhood` | `"moore"` | 8-connected; lower p_c than 4-connected VN, giving broader parameter regime for visible cascades |
| `boundary` | `"periodic"` | Torus avoids fixed-boundary effects on wavefront shape |
| `init_mode` | `"single_seed"` | One infected cell at center; produces clean circular wavefront for P22 cascade detection |
| `independent_neighbors` | `True` | Standard independent-infection model (Datta-Acharyya 2021) |
| `seed` | 42 | Fixed for reproducibility |

**Canonical positive parameters (DEFINITIVE P22 detection):**
80 × 80, `infection_prob=0.30`, `recovery_prob=0.20`, `single_seed`.
At 64 × 64, the same parameters reach only SCREENING tier due to
reduced cascade size.

---

## 5. Percolation threshold context

The SIR CA on a 2D lattice exhibits a sharp **percolation-type phase transition**
at critical infection probability p_c (Grassberger 1983):

| Neighbourhood | p_c (q = 0.1) |
|---|---|
| Moore (8-connected) | ≈ 0.038 |
| Von Neumann (4-connected) | ≈ 0.10 |

Below p_c, the epidemic dies locally (only a finite number of cells ever
become infected). Above p_c, it percolates across the entire grid. The EPC
implementation uses `infection_prob` values well above the Moore-neighbourhood
p_c (default 0.3 >> 0.038), ensuring reliable spreading cascades.

**Important calibration note:** Phase-2a panel Class C (Sprint 39) used
`infection_prob ∈ linspace(0.05, 0.18, 10)` as "failed-regime" negatives,
expecting these to be below p_c. In fact, all 10 values are above the Moore-
neighbourhood p_c ≈ 0.038. The epidemic genuinely spreads from the single
seed at these parameters; P22 correctly identifies the cascade at SCREENING
tier. The label "below percolation threshold ~0.2" referred to an effective
CONFIRMATION threshold, not the physical percolation point. This
misidentification was documented in Sprint 39 and is a known calibration
limitation (Class C ground-truth annotation issue).

The mean-field R₀ approximation `R₀ = (1 − (1−p)^n_neighbours) / q` substantially
overestimates the true lattice reproductive number because it ignores spatial
correlations (depletion of susceptibles around the wavefront). At the measured
p_c values, R₀_approx ≈ 2.7–3.4 rather than the mean-field prediction R₀_c = 1.
This discrepancy is a well-known feature of spatial epidemic models
(Grassberger 1983, Datta-Acharyya 2021 §3.1.1).

---

## 6. Reproduction status

**Sprint 51 reproduction:** Datta & Acharyya (2021) §3.1.1/Fig. 11 wavefront
speed reproduced. The Sprint 51 reproduction script
(`analysis/reproductions/p22_dattaacharyya2005.py`) implements a modified
SIR CA (fixed-duration infection t_τ = 4 timesteps + re-infection probability
p2 = 0.10) inline rather than via `epc.models.sir_epidemic` — the paper's
variant differs from the default model in using a deterministic recovery time
instead of stochastic recovery probability per step.

| Observable | Published (Datta-Acharyya 2021 Fig. 11) | Measured | Tolerance | Status |
|---|---|---|---|---|
| Wavefront speed v | 0.4405 ± 0.0008 | 0.4612 ± 0.0164 | rel. error < 15% | **PASS** |

Relative error = (0.4612 − 0.4405) / 0.4405 = 4.7%, within the ≤ 15% tolerance.

dim1 PARTIAL → **PASS**.
Output: `analysis/outputs/p22_phase2a_panel.json` (irreversibility prereq).

**Phase-2a panel dim4 (Sprint 41):** overall TNR = 1.000, Cohen's d = +∞.
Output: `analysis/outputs/p22_phase2a_panel.json`.

---

## 7. Observable extraction

State snapshots expose:

| Key | Shape/Type | Description |
|---|---|---|
| `grid` | (rows, cols) int | Cell states {0=S, 1=I, 2=R} |
| `grid_dims` | tuple (rows, cols) | Grid dimensions |
| `n_states` | int = 3 | Always 3 (S, I, R) |
| `step` | int | Timestep |
| `s_count`, `i_count`, `r_count` | int | Cell counts by state |
| `newly_infected`, `newly_recovered` | int | Transitions this step |
| `activity_density` | float | `i_count / (rows × cols)` |
| `s_fraction`, `i_fraction`, `r_fraction` | float | Normalized fractions |

The primary P22 observable is the **infection time-map** T(x,y) = the
timestep at which cell (x,y) first becomes infected. For a cascade with
a single seed, this map records the wavefront's spatial propagation pattern.
The P22 detector's primary metric is `Moran's I` computed on the time-map
across the infected region: cells infected simultaneously (same wavefront
passage) are spatially clustered, producing high spatial autocorrelation.
Spatial-shuffle null permutes the time values among infected cells
(n_permutations = 199).

Early termination: `run()` stops when `i_count == 0` (epidemic over), which
is the expected outcome for all supercritical SIR runs on finite grids.

---

## 8. Known limitations

- **Single-pass dynamics.** SIR's wavefront cannot re-enter the recovered
  region. There is no spiral formation, no sustained activity, and no
  oscillation. This is by design (SIR = epidemiological "burn-through") but
  means P22 only fires during the active epidemic. Once the cascade
  completes, the SIR state history is a dead record.

- **P13 false-positive boundary.** SIR has n_states = 3 and passes P13's
  hard prerequisite guard. However, P13 correctly rejects SIR because
  the excitable-wave screening metric (wavefront speed stability CV < 0.15)
  fails: SIR's single-pass dynamics cause activity to die out before the
  measurement window accumulates sufficient speed samples. P13/SIR
  discrimination is confirmed in `tests/test_sir_p13_boundary.py`.

- **Wavefront speed model difference.** The Sprint 51 reproduction uses a
  deterministic fixed-duration infection (t_τ = 4) with re-infection (p2 = 0.10),
  which differs from `epc.models.sir_epidemic`'s stochastic recovery (probability
  q per step). The two formulations produce quantitatively similar wavefront
  speeds but are not identical. Datta-Acharyya (2021)'s specific variant was
  implemented inline in the reproduction script rather than modifying the
  canonical EPC model.

- **dim2 still PARTIAL.** Per-N basin fraction (fraction of parameter space
  giving spreading cascades vs. dying epidemics across multiple system
  sizes N = L²) has not been systematically reported with ≥ 5-seed dispersion
  statistics per the depth_gap rubric.
