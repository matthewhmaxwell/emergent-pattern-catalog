# Sprint 44 Return — lattice_2d dim4 batch 3: P12 RPS + P13 GH

**Date:** 2026-05-25
**Base HEAD (sprint start):** `1f12266`
**Sprint goal:** Run v1.2 Phase-2a panel against P12 RPS (cyclic dominance) and P13 GH (excitable waves). No detector/model/spec changes.

---

## Part A — P12 panel (spatial RPS, cyclic dominance)

**Canonical positive:** RPSSpatialModel (rows=50, cols=50, mobility=1e-4, n_steps=200), 5 seeds.

Mobility=1e-4 is below M_c ≈ 4.5×10⁻⁴ (Reichenbach et al. 2007): three-species coexistence maintained, cyclic-dominance spirals form.

**Class C failed regimes:** `epc/phase2a/failed_regimes/p12_rps.py` — 10 high-mobility extinction regimes, mobility ∈ linspace(5e-3, 5e-2, 10). All values ≥ 11× M_c; cyclic coexistence collapses to monoculture; P12's three-species persistence gate fails.

```
pattern  overall    syn    cat    fai      d verdict
P12        1.000  1.000  1.000  1.000    inf PASS
```

- All 5 positives: CONFIRMATION, confidence=0.700 (log₁₀(min ρ) > 2.0; null p < 0.005)
- Class A synthetic TNR = 1.000 (10/10 — both perm_shuffled + time_shuffled tested; flags=(False, False))
- Class B catalog TNR = 1.000 (7/7 lattice_2d mates: P11_LV, P13_GH, P14_BTW, P15_GoL, P1_Schelling, P22_SIR, P27_NowakMay — all correctly rejected at screening)
- Class C failed regimes TNR = 1.000 (10/10 extinction regimes correctly rejected)
- Cohen's d = +∞

**Files added:**
- `epc/phase2a/failed_regimes/p12_rps.py`
- `analysis/outputs/p12_phase2a_panel.json`

**Files modified:**
- `analysis/run_phase2a_panel.py`: `build_p12_positives()`, `make_p12_detector_fn()`, `run_p12()` added
- `epc/phase2a/failed_regimes/__init__.py`: p12_rps registered

---

## Part B — P13 panel (Greenberg-Hastings, excitable waves)

**Canonical positive:** GreenbergHastings (rows=50, cols=50, n_states=8, threshold=1, moore, random, density=0.3, n_steps=300), 5 seeds.

n_states=8, threshold=1, Moore neighbourhood is the canonical spiral-forming regime (Greenberg & Hastings 1978; Fisch, Gravner & Griffeath 1991). n_steps=300 ensures ≥5 T_prop timescales for wavefront-speed CV metric to stabilise.

**Class C failed regimes:** `epc/phase2a/failed_regimes/p13_gh.py` — 10 low-density init regimes, density ∈ linspace(0.01, 0.10, 10). Below excitation threshold for persistent wavefronts: sparse isolated excited clusters burn out before sustaining spiral nucleation; P13 screening prerequisite (persistent wavefront ≥ 5×T_prop) not satisfied → rejected.

```
pattern  overall    syn    cat    fai      d verdict
P13        1.000  1.000  1.000  1.000    inf PASS
```

- All 5 positives: SCREENING, confidence=0.500 (persistent wavefront CV < 0.20 satisfied; CONFIRMATION 50-rotation gate not reached within 300 steps at 50×50)
- Class A synthetic TNR = 1.000 (10/10)
- Class B catalog TNR = 1.000 (7/7: P11_LV, P12_RPS, P14_BTW, P15_GoL, P1_Schelling, P22_SIR, P27_NowakMay — all correctly rejected)
- Class C failed regimes TNR = 1.000 (10/10 low-density regimes correctly rejected)
- Cohen's d = +∞ (positives 0.500 vs all negatives 0.000 — sharp class separation)

**Note on positive tier:** P13 positives reach SCREENING (not CONFIRMATION) with these parameters. The panel PASS criterion is satisfied via TNR + Cohen's d; the positive tier is informational. CONFIRMATION requires ≥50 spiral-tip rotations; 300-step trajectories on 50×50 accumulate sufficient wavefront statistics for the CV gate but insufficient rotation count. This is a documentation carry-forward, not a panel failure.

**Files added:**
- `epc/phase2a/failed_regimes/p13_gh.py`
- `analysis/outputs/p13_phase2a_panel.json`

**Files modified:**
- `analysis/run_phase2a_panel.py`: `build_p13_positives()`, `make_p13_detector_fn()`, `run_p13()` added
- `epc/phase2a/failed_regimes/__init__.py`: p13_gh registered

---

## Part C — REPLICATION_NOTES

Two new sections appended:
- `## Phase-2a Panel Result (v1.2) — Sprint 44 (P12 spatial RPS, PASS)`
- `## Phase-2a Panel Result (v1.2) — Sprint 44 (P13 Greenberg-Hastings, PASS)`

Sprint 30 rule notation (PASS → depth_gap update; no detector/model changes).

---

## Part D — depth_gap.md

- **P12 row:** dim4 PARTIAL → **PASS**; grade remains **GAP** (dim1: λ ∝ √M wavelength scaling not replicated, C3 carry-forward, Sprint 48 target; dim2: mostly single-seed characterisation).
- **P13 row:** dim4 PARTIAL → **PASS**; all 4 dims now PASS → grade **GAP → AT-DEPTH**.
- **AT-DEPTH count: 7 → 8 / 19** (P13 added; P3, P9, P15, P18, P27, P28, P31 unchanged).
- Aggregate findings block updated.

---

## Part E — Paper sync

- `docs/paper_section4_draft.md` §4.22 added: Sprint 44 P12 + P13 panel results.
- `docs/paper_section6_draft.md`: AT-DEPTH count 7→8; Sprint 44 narrative added.
- `docs/paper_CHANGELOG.md`: Sprint 44 entry added.

---

## Part F — Carry-forward review

**Carry-forwards closed this sprint:** None.

**Active carry-forwards from prior sprints (unchanged):**
- C-p1-time-shuffle-fp (Sprint 43): P1 `time_shuffled` Class A FP.
- C-p1-linear-gradient-fp (Sprint 43): P1 `linear_gradient` Class A FP.
- C-p1-class-c-subthreshold-fp (Sprint 43): P1 Class C calibration issue at density=0.9.
- C3 / C-p12-wavelength-scaling (Sprint 9 / Sprint 26): P12 λ ∝ √M wavelength scaling not replicated. Sprint 48 target per brief.
- C-class-a-constant-field-trivial-sync (Sprint 35): P9 constant_field degenerate substrate.

**Carry-forwards opened this sprint:**
- **C-p13-positive-tier-screening**: P13 canonical positives (GH 50×50, n_states=8, n_steps=300) reach SCREENING (confidence=0.500) but not CONFIRMATION. CONFIRMATION requires ≥50 spiral-tip rotations; longer runs or larger grids would likely close this gap. Informational only — panel PASS not affected. Deferred to a future cleanup sprint.

---

## Pre-flight / post-flight

- Pre-flight: HEAD verified at `1f12266`; `pytest tests/ -m "not slow"` confirmed passing (orchestration 87, phase2a panel 91, P12/P13 specific tests 50 = all pass).
- Panel runs: P12 and P13 panels executed and JSONs written.
- Post-flight: `test_orchestration.py` + `test_cross_detection_matrix.py` passing (87 tests); `test_phase2a_panel.py` passing (91 tests); transfer matrix unchanged (no model/detector/orchestration changes).

---

**Decision: GO**

P12: PASS (TNR=1.000, d=+∞). dim4 PARTIAL→PASS. Grade remains GAP (dim1+dim2 PARTIAL).
P13: PASS (TNR=1.000, d=+∞). dim4 PARTIAL→PASS; all dims PASS → AT-DEPTH. AT-DEPTH count: 8/19.
New carry-forward C-p13-positive-tier-screening (informational; panel not affected). Chain may proceed to Sprint 45 (final lattice_2d dim4 — P11 LV singleton).
