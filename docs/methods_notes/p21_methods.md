# P21 Methods Note — Polarization / Fragmentation (Hegselmann-Krause)

**Pattern:** P21 — Polarization / fragmentation
**Canonical model:** Hegselmann-Krause bounded-confidence opinion dynamics (2002)
**Detector:** `epc/detectors/p21_polarization.py`
**Model:** `epc/models/hegselmann_krause.py`
**Reproduction sprint:** Sprint 53 (`analysis/reproductions/p21_hegselmann2002.py`)

---

## 1. Pattern and canonical reference

P21 is the emergence of opinion polarization or fragmentation in a population
of agents who average only with neighbours who are sufficiently similar to
themselves. Below a critical confidence bound ε_c, the initial uniform
distribution self-organizes into two or more stable opinion clusters. At
ε ≈ 0.1 the result is fragmentation (many small clusters); at ε ≈ 0.2 it
is polarization (two clusters); above ε_c ≈ 0.25–0.27 it is consensus
(one cluster).

**Primary reference:** Hegselmann, R. & Krause, U. (2002). Opinion dynamics
and bounded confidence: models, analysis, and simulation. *Journal of
Artificial Societies and Social Simulation*, 5(3), 2.

---

## 2. Model equations and update rule

The HK model places N agents with continuous opinions x_i(t) ∈ [0,1]. At
each round, every agent simultaneously averages over all agents within
confidence bound ε:

```
x_i(t+1) = (1 / |N_i(t)|) Σ_{j ∈ N_i(t)} x_j(t)
```

where `N_i(t) = {j : |x_i(t) − x_j(t)| ≤ ε}` (always includes i itself).

**Update mode: synchronous.** All agents read the opinions at time t and
write new opinions simultaneously. This is the canonical HK model as defined
in Hegselmann-Krause (2002). An alternative asynchronous variant exists in
the literature (random sequential update) and produces slightly different
quantitative ε_c values; the EPC implementation uses exclusively synchronous
updates to match the paper's replication target.

**Convergence detection.** After each synchronous round, the maximum opinion
change is computed:

```
max_change = max_i |x_i(t+1) − x_i(t)|   (L∞ norm)
```

If `max_change < convergence_tol`, the system is flagged as converged and
subsequent `step()` calls are no-ops (returning the frozen state). The model
default is `convergence_tol = 1e-8`. The Sprint 53 reproduction uses
`convergence_tol = 1e-6` (matching the paper's description of numerical
convergence). For practical purposes both values yield the same cluster
structure — the L∞ norm at actual convergence is numerically zero once
clusters are separated by gaps > ε, since isolated clusters have no
cross-boundary neighbours.

**Cluster counting.** The model computes `n_clusters` at each step via
gap detection: sort opinions, compute consecutive differences, count gaps
`> gap_threshold` (default `gap = ε/2`). This is the standard method
in the HK literature and matches the algorithm in Hegselmann-Krause (2002).
The gap threshold `ε/2` is proportional to the confidence bound: a gap larger
than half the interaction range implies two opinions would not average with
each other even if initialized at the gap boundaries.

**Vectorized inner loop.** For N ≤ 5000, the pairwise difference matrix
`|x_i − x_j|` is computed as a dense O(N²) operation using NumPy broadcasting.
This is fast enough for the catalog's typical N = 100–500 usage. For N > 5000,
a sorted-array sweep with binary search would reduce cost to O(N log N) per
step; this optimization is not implemented.

---

## 3. Parameter defaults and justifications

| Parameter | Default | Justification |
|---|---|---|
| `n_agents` | 500 | Sprint 5 canonical; Sprint 53 reproduction uses N=100 to match HK (2002) Fig. 2 |
| `epsilon` | 0.2 | Canonical polarization regime: produces 2 stable clusters from uniform IC |
| `init_mode` | `"uniform"` | U[0,1] initial opinions, matches HK (2002) §3 |
| `convergence_tol` | 1e-8 | Stricter than paper's nominal 1e-6; results identical in practice |
| `seed` | 42 | Fixed for reproducibility |

**Init modes available:** `"uniform"` (U[0,1], canonical), `"gaussian"`
(N(0.5, 0.2) clipped), `"bimodal"` (mixture of two Gaussians at 0.3 and 0.7).
The bimodal mode is used for Phase-2a panel negative-control tests to verify
that the detector does not spuriously promote pre-existing structure.

---

## 4. Implementation choices

**Synchronous vs. asynchronous.** The synchronous rule (all agents update
simultaneously) is the defining rule of the HK model and is used in all
Hegselmann-Krause (2002) simulations. Asynchronous variants (random activation)
break the symmetry between agents at different opinion positions and tend to
produce smoother cluster merging dynamics but do not change the qualitative
steady-state cluster count. The EPC implementation is exclusively synchronous
for replication fidelity.

**Convergence and early stopping.** The `run()` method stops early on
convergence to avoid running unnecessary no-op steps. State history length
therefore varies by ε: low ε (many clusters, fast convergence) produces
shorter histories than high ε (consensus after slow drift). Detectors must
handle variable-length histories; `P21PolarizationDetector` uses the full
history up to convergence as the measurement window.

**No noise term.** The model is deterministic conditional on the initial
condition. The only stochasticity is in the random initialization (seeded
by `seed`). Adding noise (e.g., re-sampling `ε` each step) would prevent
convergence; this extension is not in the catalog. The dip-test bootstrap
in the detector introduces Monte Carlo variability, but the model itself
is deterministic once initialized.

**Opinion domain [0, 1].** Opinions are bounded to the unit interval. Boundary
effects are mild because initial conditions are U[0,1] and the averaging rule
never creates opinions outside the range of existing opinions. The steady-state
clusters near 0 and 1 at low ε are genuine boundary effects: agents at the
extreme ends have fewer neighbours and may form singleton clusters.

---

## 5. Deviations from canonical (and why)

**N = 100 vs. N = 500 default.** The Sprint 53 reproduction uses N = 100
to match Hegselmann-Krause (2002) Fig. 2, which uses N = 100. The model's
default N = 500 was set in Sprint 5 for detection robustness (more agents
→ sharper dip statistic). The methods note distinguishes the two uses:
N = 100 for reproduction fidelity, N = 500 for detection power.

**ε_c = 0.25–0.27 boundary zone.** At ε = 0.25, finite-size effects cause
stochastic outcomes: with N = 100, 14/20 seeds reach consensus (1 cluster)
and 6/20 reach polarization (2 clusters). The Sprint 53 reproduction widens
the published tolerance to [1, 3] at ε = 0.25, consistent with the paper's
transition-zone discussion. This boundary behaviour is a fundamental property
of the HK model, not an implementation artifact.

**gap_threshold = ε/2 vs. 0.05.** The `_count_clusters()` function in the
model uses `gap = ε/2` (proportional to confidence bound). The Sprint 53
reproduction script and REPLICATION_NOTES use a fixed threshold of 0.05
(matching the paper's description of "gaps between the groups"). Both produce
identical cluster counts at convergence because HK steady-state clusters are
separated by gaps >> ε; the distinction only matters during transient merging
steps.

---

## 6. Reproduction status

**Sprint 53 reproduction:** Hegselmann & Krause (2002) Fig. 2 cluster-count
vs. ε curve reproduced across 8 ε points (0.10 to 0.50), 20 seeds each,
N = 100, uniform initial conditions, synchronous update, convergence
tolerance 1e-6 or T = 10000 steps.

| ε | Published range | Measured median | Status |
|------|----------------|-----------------|--------|
| 0.10 | [4, 7] | 4 | PASS |
| 0.15 | [3, 5] | 3 | PASS |
| 0.20 | [2, 4] | 2 | PASS |
| 0.25 | [1, 3]† | 1 | PASS |
| 0.27 | [1, 2] | 1 | PASS |
| 0.30 | [1, 1] | 1 | PASS |
| 0.50 | [1, 1] | 1 | PASS |

†ε = 0.25 is in the 2→1 transition zone (ε_c ≈ 0.24–0.27).

dim1 PARTIAL → **PASS**.
Output: `analysis/outputs/p21_hegselmann2002_reproduction.json`.

**Phase-2a panel dim4 (Sprint 49):** overall TNR = 0.955, Cohen's d = 5.487.
`time_shuffled` FP at CONFIRMATION (0.850) — carry-forward C-p21-time-shuffled-fp.
Output: `analysis/outputs/p21_phase2a_panel.json`.

---

## 7. Observable extraction

State snapshots expose:

| Key | Shape/Type | Description |
|---|---|---|
| `opinions` | (N,) float64 | Agent opinions in [0, 1] |
| `step` | int | Interaction round |
| `n_clusters` | int | Cluster count via gap detection (gap = ε/2) |
| `variance` | float | Opinion variance `var(x)` |
| `converged` | bool | True if L∞ norm < convergence_tol |

The `opinions` array is the primary observable for the P21 detector: the
Hartigan dip test on the final-state opinion distribution tests for
multimodality (bimodal → polarization, multimodal → fragmentation).
The `n_clusters` field provides a fast diagnostic but is not used directly
by the detector's primary test.

Metadata exposes:

| Key | Value |
|---|---|
| `model` | `"hegselmann_krause"` |
| `epsilon` | confidence bound |
| `n_agents` | N |
| `init_mode` | initialization mode |
| `substrate` | `"opinion_space"` |

---

## 8. Known limitations

- **P21 detector: dip test requires ≥ 2 clusters to be present.** Hartigan's
  dip test on a unimodal final-state distribution (ε ≥ ε_c) correctly
  returns `detected=False`. But at ε slightly below ε_c (stochastic boundary),
  the initial steps are unimodal and the final state is bimodal; the time-
  shuffled Class A substrate (which mixes steps from all timepoints) can still
  show a weak bimodal signal from early-convergence snapshots, causing the
  known C-p21-time-shuffled-fp carry-forward.

- **dim2 still PARTIAL.** The ε-sweep study was mostly single-seed; multi-seed
  characterization with ≥ 5 seeds and explicit dispersion statistics on the
  cluster-count distribution across seeds has not been formally reported
  per the depth_gap rubric.

- **N-dependence of ε_c.** The critical confidence bound ε_c depends on N.
  At N = 100 it is approximately 0.24–0.27; at N → ∞ the theoretical ε_c
  for uniform initial conditions converges to ≈ 0.265 (Hegselmann-Krause
  2002 §5). The detection parameters are calibrated for N = 100–500.

- **No spatial embedding.** The HK model has no spatial structure — agents
  interact purely in opinion space. The `opinion_space` substrate type in
  the orchestration registry reflects this; spatial detectors (P1, P12,
  P13, etc.) are correctly rejected at the substrate-type gate.
