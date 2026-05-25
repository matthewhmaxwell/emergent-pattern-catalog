# Sprint 46 Return — continuous_2d dim4 batch: P2 ABP + P5 Vicsek + P6 D'Orsogna

**Date:** 2026-05-25
**Base HEAD (sprint start):** `b25c18e`
**Sprint goal:** Run v1.2 panel against all 3 continuous_2d patterns (P5, P2, P6). No detector/model/spec changes.

---

## Part A — P5 panel (Vicsek / flocking)

**Canonical positive:** VicsekModel(N=300, box_size=7.0, speed=0.03, noise=0.1, n_steps=5000), 5 seeds (0–4). T_cross = 233 steps; n_steps=5000 → measurement window ≈ 10.7 T_cross (satisfies ≥10 T_cross DEFINITIVE prerequisite).

**Class C failed regimes:** `epc/phase2a/failed_regimes/p5_vicsek.py` — 10 high-noise regimes, noise ∈ linspace(0.70, 1.50, 10), N=300, box_size=7.0, n_steps=200. All regimes well above order-disorder transition (~noise=0.4 for these parameters); φ ≈ 0 in all cases → P5 rejects at screening.

```
pattern  overall    syn    cat    fai      d verdict
P5         0.957  0.889  1.000  1.000  4.987 PASS-with-weakness
```

- All 5 positives: DEFINITIVE, confidence=0.850
- Class A synthetic TNR = 0.889 (8/9 evaluated; `permutation_shuffled` SKIPPED — `permutation_invariant=True`; `time_shuffled` FP at DEFINITIVE — carry-forward C-p5-time-shuffle-fp)
- Class B catalog TNR = 1.000 (advisory: 2 mates + 2 supplements; P2_abp and P6_dorsogna rejected via unaligned headings; random-walk supplements rejected)
- Class C failed regimes TNR = 1.000 (10/10 correctly rejected)
- Cohen's d = 4.987

**Invariance flag check (P5):** `epc/phase2a/detector_invariance.py` confirmed `P5: (True, False)` — `permutation_invariant=True, time_shuffle_invariant=False`. Flag was already present; no auto-flip performed.

**Files added:**
- `epc/phase2a/failed_regimes/p5_vicsek.py`
- `analysis/outputs/p5_phase2a_panel.json`

**Files modified:**
- `analysis/run_phase2a_panel.py`: `build_p5_positives()`, `make_p5_detector_fn()`, `run_p5()` added; `p5` case in `main()`
- `epc/phase2a/synthetic.py`: `"particles"` format branches added to all 10 Class A generators + `_particles_history_from_random()` helper
- `epc/phase2a/catalog.py`: `_adapt_to_particles()` adapter + `"particles"` case in `load_catalog_substrate_for_format()`
- `epc/phase2a/failed_regimes/__init__.py`: p5_vicsek registered

---

## Part B — P2 panel (Active Brownian Particles / MIPS)

**Canonical positive:** ActiveBrownianParticles(N=800, phi_packing=0.5, Pe=100, v0=1.0, D_r=0.01, box≈35.4, rho_star=4.0, r_cg=1.0, dt=0.05, n_steps=2500), 5 seeds (0–4).

**Class C failed regimes:** `epc/phase2a/failed_regimes/p2_active_brownian.py` — 10 low-Pe regimes, Pe ∈ linspace(0.50, 10.00, 10), N=400, phi=0.5, n_steps=600. All well below MIPS onset (Pe_c ≈ 40–60 for these parameters); uniform density → two_phase_score ≈ 0 → P2 rejects at screening.

```
pattern  overall    syn    cat    fai      d verdict
P2         0.958  0.900  1.000  1.000  3.401 PASS
```

- 3/5 positives: DEFINITIVE (confidence=0.950), 1/5 screening (0.600), 1/5 screening (0.000) — MIPS nucleation is seed-dependent; burn-in at seeds 0–1 insufficient within 2500 steps
- Class A synthetic TNR = 0.900 (9/10 evaluated; `permutation_shuffled` FP at screening — two_phase_score is spatial-distribution invariant; carry-forward C-p2-perm-shuffled-fp; flag NOT auto-set per brief)
- Class B catalog TNR = 1.000 (advisory: P5_vicsek and P6_dorsogna + 2 supplements)
- Class C failed regimes TNR = 1.000 (10/10 correctly rejected)
- Cohen's d = 3.401

**Invariance flag check (P2):** `epc/phase2a/detector_invariance.py` — P2 NOT in dict (defaults to `(False, False)` via `_DEFAULT`). Brief instructs: DO NOT auto-flip; flag in carry-forward. Carry-forward C-p2-perm-shuffled-fp opened.

**Files added:**
- `epc/phase2a/failed_regimes/p2_active_brownian.py`
- `analysis/outputs/p2_phase2a_panel.json`

**Files modified:**
- `analysis/run_phase2a_panel.py`: `build_p2_positives()`, `make_p2_detector_fn()`, `run_p2()` added; `p2` case in `main()`
- `epc/phase2a/failed_regimes/__init__.py`: p2_active_brownian registered

*(synthetic.py and catalog.py changes shared with Part A above)*

---

## Part C — P6 panel (D'Orsogna SPP / milling)

**Canonical positive:** DOrsognaSPPModel(N=100, C_a=0.5, C_r=1.0, l_a=3.0, l_r=0.5, alpha=1.0, beta=0.5, dt=0.05, init_mode="ring", init_radius=5.0, n_steps=3000), 5 seeds (0–4). Ring init ensures immediate milling onset.

**Class C failed regimes:** `epc/phase2a/failed_regimes/p6_dorsogna.py` — 10 mismatched-radii regimes, l_a ∈ linspace(0.10, 0.49, 10) with l_r=0.5 fixed. Throughout, l_a ≤ l_r → attraction range ≤ repulsion range → no stable orbit forms → |L| ≈ 0 → P6 rejects at screening.

```
pattern  overall    syn    cat    fai      d verdict
P6         0.958  0.900  1.000  1.000  5.087 PASS
```

- All 5 positives: DEFINITIVE, confidence=0.850
- Class A synthetic TNR = 0.900 (9/10; `time_shuffled` FP at DEFINITIVE — milled frames retain |L|>0 independent of temporal order; carry-forward C-p6-time-shuffle-fp)
- Class B catalog TNR = 1.000 (advisory: P5_vicsek and P2_abp + 2 supplements)
- Class C failed regimes TNR = 1.000 (10/10 correctly rejected)
- Cohen's d = 5.087

**Invariance flag check (P6):** Confirmed `P6: (False, False)` in `detector_invariance.py`. Not permutation-invariant (angular momentum requires positions+velocities to be coherent). No flag change needed.

**Files added:**
- `epc/phase2a/failed_regimes/p6_dorsogna.py`
- `analysis/outputs/p6_phase2a_panel.json`

**Files modified:**
- `analysis/run_phase2a_panel.py`: `build_p6_positives()`, `make_p6_detector_fn()`, `run_p6()` added; `p6` case in `main()`
- `epc/phase2a/failed_regimes/__init__.py`: p6_dorsogna registered

---

## Part D — depth_gap.md

- **P5 row:** dim4 PARTIAL → **PASS**; all 4 dims now PASS → grade **GAP → AT-DEPTH**.
- **P2 row:** dim4 PARTIAL → **PASS**; dims 1–3 still PARTIAL; grade remains GAP.
- **P6 row:** dim4 PARTIAL → **PASS**; dim2 still PARTIAL; grade remains GAP.
- **AT-DEPTH count: 8 → 9 / 19** (P3, P5, P9, P13, P15, P18, P27, P28, P31).
- Aggregate findings block updated: Sprint 46 finding added; dim4 resolved patterns list updated; AT-DEPTH/Gap counts updated.

---

## Part E — REPLICATION_NOTES

New sections appended:
- `## Phase-2a Panel Result (v1.2) — Sprint 46 (P5 Vicsek / flocking, PASS-with-weakness)`
- `## Phase-2a Panel Result (v1.2) — Sprint 46 (P2 ABP / MIPS, PASS)`
- `## Phase-2a Panel Result (v1.2) — Sprint 46 (P6 D'Orsogna / milling, PASS)`

Sprint 30 rule notation (PASS → depth_gap update; no detector/model changes).

---

## Part F — Paper sync

- `docs/paper_section4_draft.md` §4.24 added: Sprint 46 P5/P2/P6 panel results.
- `docs/paper_section6_draft.md`: Sprint 46 AT-DEPTH narrative added (count 9/19; P5 AT-DEPTH, P2/P6 dim4 PASS).
- `docs/paper_CHANGELOG.md`: Sprint 46 entry added.

---

## Part G — Carry-forward review

**Carry-forwards closed this sprint:** None.

**Carry-forwards opened this sprint:**
- **C-p5-time-shuffle-fp (Sprint 46):** P5 `time_shuffled` Class A FP. φ = |⟨e^iθ⟩| is a per-frame alignment metric; every frame in a flocked trajectory has high φ regardless of temporal order. Resolution requires `time_shuffle_invariant=True` flag in `detector_invariance.py` (spec decision out of scope Sprint 30 rule) or Class A substrate redesign. Informational only; panel PASS-with-weakness.
- **C-p2-perm-shuffled-fp (Sprint 46):** P2 `permutation_shuffled` Class A FP. two_phase_score is a spatial density statistic computed on coarse grid; invariant to particle-index permutation. Resolution requires `permutation_invariant=True` flag in `detector_invariance.py`. Brief explicitly instructs DO NOT auto-flip; flag for spec review. When flag is set, permutation_shuffled will be SKIPPED and syn TNR will rise to 1.000.
- **C-p6-time-shuffle-fp (Sprint 46):** P6 `time_shuffled` Class A FP. |L| = |Σ r_i × v_i| / N is a per-frame angular momentum metric; every frame in a milled trajectory has high |L|. Analogous to C-p5-time-shuffle-fp above. Resolution: same spec decision path.

**Active carry-forwards from prior sprints (unchanged):**
- C-p1-time-shuffle-fp (Sprint 43): P1 `time_shuffled` Class A FP.
- C-p1-linear-gradient-fp (Sprint 43): P1 `linear_gradient` Class A FP.
- C-p1-class-c-subthreshold-fp (Sprint 43): P1 Class C calibration issue at density=0.9.
- C3 / C-p12-wavelength-scaling (Sprint 9 / Sprint 26): P12 λ ∝ √M wavelength scaling not replicated. Sprint 48 target per brief.
- C-class-a-constant-field-trivial-sync (Sprint 35): P9 constant_field degenerate substrate.
- C-p13-positive-tier-screening (Sprint 44): P13 canonicals reach SCREENING (0.500) not CONFIRMATION; informational only.

---

## Pre-flight / post-flight

- Pre-flight: HEAD verified at `b25c18e`; `pytest tests/test_orchestration.py tests/test_phase2a_panel.py` confirmed 154 tests passing (11.6s).
- Infrastructure added: `synthetic.py` particles format, `catalog.py` `_adapt_to_particles()` adapter, 3 failed-regime modules, 3 panel runner functions.
- Panel runs executed; JSONs written.
- Post-flight: `test_orchestration.py` + `test_phase2a_panel.py` + `test_cross_detection_matrix.py` passing (178 tests in 102s); transfer matrix unchanged (no model/detector/orchestration changes in this sprint).

---

**Decision: GO**

P5: PASS-with-weakness (TNR=0.957, d=4.987); all 5 positives DEFINITIVE. dim4 PARTIAL→PASS; P5 advances to AT-DEPTH (9/19).
P2: PASS (TNR=0.958, d=3.401); 3/5 positives DEFINITIVE. dim4 PARTIAL→PASS; grade remains GAP.
P6: PASS (TNR=0.958, d=5.087); all 5 positives DEFINITIVE. dim4 PARTIAL→PASS; grade remains GAP.
No new carry-forwards that block chain. Three time_shuffle/permutation FPs opened as informational carry-forwards. AT-DEPTH count: **9 / 19**. Chain may proceed to Sprint 47 (mixed dim4 — P8 lattice_1d Nagel-Schreckenberg + P10 chimera oscillator).
