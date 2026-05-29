# Sprint 57 Return — P2 + P21 + P22 dim3 closure: methods-note authoring

**Date:** 2026-05-29
**Base HEAD (sprint start):** `7cd824c`
**Sprint goal:** Author dim3 methods notes for P2 (MIPS/ABP), P21 (Hegselmann-Krause opinion dynamics), P22 (SIR epidemic CA). Close dim3 depth gap on all three.
**Tag:** `v0.57.0`
**Sprint type:** Pure documentation — no code, model, or detector changes.

---

## Part A — Existing methods-note template

No prior `docs/methods_notes/` directory existed. The sprint created it and established the standard methods-note structure from the brief:

- Pattern + canonical reference
- Model equations / update rule
- Parameter defaults + their justifications
- Implementation choices (synchronous vs async, boundary conditions, etc.)
- Deviations from canonical (if any), and why
- Reproduction status (cross-link to Sprint X reproduction)
- Observable extraction
- Known limitations

Referenced prior examples in depth_gap.md and REPLICATION_NOTES.md to calibrate scope and depth for "substantive" dim3 PASS.

---

## Part B — Methods notes authored

### P2 — `docs/methods_notes/p2_methods.md`

**Implementation:** Fily-Marchetti (2012) ABP model (not Cahn-Hilliard or Ising). Overdamped Langevin equations: `dr_i/dt = v(ρ_i) ê(θ_i)`, `dθ_i/dt = √(2 D_r) ξ_i`, with density-dependent speed `v(ρ) = v₀ max(0, 1 − ρ/ρ*)`. Euler-Maruyama integration; periodic cKDTree for local density.

**Key content:**
- ABP vs. Cahn-Hilliard/Ising distinction: ABP is a particle-based non-equilibrium model; MIPS universality class differs from Model-B coarsening
- Primary metric: `two_phase_coexistence_score = min(f_gas, f_liquid)` where `f_gas = P(ρ_i < ρ*/2)`, `f_liquid = P(ρ_i > ρ*)` — documents why this is chosen over FFT/structure-factor approaches
- FFT structure-factor context: radial S(k) is the standard MIPS domain-size characterization in Fily-Marchetti / Redner et al.; EPC uses phase fractions as a robust steady-state proxy that doesn't require coarsening-time tracking
- Hartigan dip failure (ADR 44): integer-count discrete distributions universally reject at floor p=0.005 regardless of regime
- Mechanistic-null metadata flags (ADR 43): `has_density_dependent_speed`, `has_alignment_rule=False`, `has_attraction_rule=False`
- Burn-in / nucleation-lag requirements; three false-positive traps (thermal, dilute, over-saturated)
- Sprint 52 reproduction: PASS at (φ=0.5, Pe=100): two_phase_score=0.1237±0.077, r=−0.958±0.020, thermal score=0.052
- Known limitations: finite-size nucleation lag, no cluster morphology tracking, ρ* must be known, dim2 still PARTIAL

**~8 sections, ~650 words.**

### P21 — `docs/methods_notes/p21_methods.md`

**Implementation:** Hegselmann-Krause (2002) synchronous bounded-confidence averaging. Update rule: `x_i(t+1) = (1/|N_i(t)|) Σ_{j: |x_i - x_j| ≤ ε} x_j(t)`. All agents update simultaneously.

**Key content:**
- Synchronous vs. asynchronous distinction: HK (2002) uses synchronous; asynchronous produces slightly different ε_c but same qualitative behaviour
- Convergence detection: `||x(t+1) − x(t)||_∞ < tol`; model default tol=1e-8; Sprint 53 reproduction uses 1e-6 (paper convention); in practice identical at convergence since isolated clusters have no cross-boundary neighbours
- Cluster counting: sorted-gap detection with gap = ε/2 (model internal) vs. 0.05 (reproduction script); produce identical results at convergence
- N=100 (HK 2002 Fig. 2 replication) vs. N=500 (detection power) distinction
- ε_c ≈ 0.24–0.27 finite-size boundary zone; ε=0.25 stochastic (14/20 consensus, 6/20 two-cluster at N=100)
- No noise: model deterministic conditional on IC; bootstrap variability comes from dip-test
- C-p21-time-shuffled-fp carry-forward context: pre-convergence unimodal steps in shuffled history
- Sprint 53 reproduction: PASS, all 8 ε points within HK (2002) Fig. 2 tolerance
- Known limitations: dim2 PARTIAL, N-dependence of ε_c, O(N²) dense inner loop

**~8 sections, ~600 words.**

### P22 — `docs/methods_notes/p22_methods.md`

**Implementation:** Probabilistic SIR CA on 2D lattice. States: S=0, I=1, R=2. Synchronous update. Infection: `P(S→I) = 1 − (1−p)^{n_infected_nbrs}` (independent-neighbours). Recovery: `P(I→R) = q` per step.

**Key content:**
- S/I/R encoding (0/1/2); integer state constants SUSCEPTIBLE=0, INFECTED=1, RECOVERED=2
- Synchronous update; independent-neighbours vs. threshold-model alternative
- Irreversibility prerequisite (Sprint 41): `_check_irreversibility_prereq()` rejects substrates with `curr_state < prev_state`; eliminates LV (PREDATOR→EMPTY=2→0) and RPS (cyclic→EMPTY=0); grounds in Datta-Acharyya (2021) + Mobilia et al. (2007) + Reichenbach (2007)
- Percolation threshold: Moore p_c ≈ 0.038, VN p_c ≈ 0.10 at q=0.1; default p=0.3 well above threshold
- Class C calibration issue: infection_prob 0.05–0.18 is above Moore p_c; "below percolation threshold ~0.2" referred to CONFIRMATION threshold, not physical p_c
- Mean-field R₀ overestimation and spatial-correlation explanation (Grassberger 1983)
- Model difference between Sprint 51 reproduction (fixed t_τ=4 + re-infection p2=0.10, implemented inline) vs. `epc.models.sir_epidemic` (stochastic geometric recovery)
- Sprint 51 reproduction: PASS, wavefront speed 0.4612 ± 0.0164 vs. published 0.4405 ± 0.0008 (4.7% rel. error)
- Known limitations: single-pass dynamics, P13 false-positive boundary (correctly handled), dim2 PARTIAL

**~8 sections, ~620 words.**

---

## Part C — REPLICATION_NOTES + depth_gap updates

**REPLICATION_NOTES.md changes:**

| Section | Added |
|---|---|
| Sprint 16 / Sprint 52 (P2) | `## Dim3 Methods Note — Sprint 57` with coverage summary + `docs/methods_notes/p2_methods.md` reference. dim3 PARTIAL→PASS. |
| Sprint 5 / Sprint 53 (P21) | `## Dim3 Methods Note — Sprint 57` with coverage summary + `docs/methods_notes/p21_methods.md` reference. dim3 PARTIAL→PASS. |
| Sprint 7–8 / Sprint 51 (P22) | `## Dim3 Methods Note — Sprint 57` with coverage summary + `docs/methods_notes/p22_methods.md` reference. dim3 PARTIAL→PASS. |

**depth_gap.md changes:**

| Field | Before | After |
|---|---|---|
| P2 dim3 | PARTIAL | **PASS** |
| P2 grade | GAP | GAP (dim2 still PARTIAL) |
| P2 effort | M | M |
| P21 dim3 | PARTIAL | **PASS** |
| P21 grade | GAP | GAP (dim2 still PARTIAL) |
| P21 effort | L | **M** (dim3 PASS reduces remaining work to dim2 only) |
| P22 dim3 | PARTIAL | **PASS** |
| P22 grade | GAP | GAP (dim2 still PARTIAL) |
| P22 effort | L | **M** (dim3 PASS reduces remaining work to dim2 only) |
| C5 carry-forward | Open | **CLOSED** |
| Sprint 57 finding | — | Added |
| AT-DEPTH count | 13/19 | 13/19 (unchanged — all three still have dim2 PARTIAL) |

---

## Part D — Paper sync

**§3.4 (Substrate-Aware Dispatch):** Added 3-sentence paragraph after the block-diagonal structure sentence. Documents the three new methods notes and their core content, with file links.

**§4.2 (Greenberg-Hastings / P13):** Added cross-reference to `p22_methods.md §5` for the P13/SIR boundary test (both n_states=3, but P13 requires refractory→resting cycle absent in SIR).

**§4.9 (HK / P21):** Added paragraph after Sprint 53 reproduction table cross-referencing `p21_methods.md`.

**§4.10 (SIR / P22):** Added paragraph after Sprint 51 numerical reproduction result cross-referencing `p22_methods.md`.

**§4.15 (ABP / P2):** Added paragraph after Sprint 52 numerical reproduction result cross-referencing `p2_methods.md`.

**§6.11 aggregate:** Sprint 57 paragraph added after Sprint 53 paragraph. Documents P2+P21+P22 dim3 closures, C5 CLOSED, AT-DEPTH count unchanged at 13/19.

**paper_CHANGELOG.md:** Sprint 57 entry added at top. Lists all 3 methods-note files, 3 REPLICATION_NOTES updates, depth_gap changes, and paper section cross-references.

---

## Part E — Post-flight

```
python3 -m pytest tests/test_orchestration.py tests/test_cross_detection_matrix.py -m "not slow" -q --tb=no

87 passed in 92.15s
```

No code modified — pure documentation sprint. No regressions expected or observed on the key registry tests.

---

## Acceptance criteria

- [x] 3 new methods-note files exist with ~6 sections each (`docs/methods_notes/p2_methods.md`, `p21_methods.md`, `p22_methods.md`)
- [x] REPLICATION_NOTES P2/P21/P22 reference the new files (Dim3 Methods Note sections added)
- [x] depth_gap P2 + P21 + P22 dim3 → PASS
- [x] Paper §3.4 + §4.2/§4.9/§4.10/§4.15 + §6 + CHANGELOG synced
- [x] `pytest tests/test_orchestration.py tests/test_cross_detection_matrix.py -m "not slow"` — 87 passed, 0 failed
- [x] Commit + tag `v0.57.0`

---

## Final commit hash and tag

**Commit:** (post-sprint commit)
**Tag:** `v0.57.0`

---

**Decision: GO**
