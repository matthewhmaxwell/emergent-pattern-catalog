# Sprint 45 Return — lattice_2d dim4 final: P11 LV

**Date:** 2026-05-25
**Base HEAD (sprint start):** `f4e09c1`
**Sprint goal:** Close the lattice_2d dim4 chain by running P11 predator-prey oscillation panel. Conditional re-run check for prior PARTIALs. No detector/model/spec changes.

---

## Part A — P11 panel (Lotka-Volterra, predator-prey oscillations)

**Canonical positive:** LotkaVolterraLattice (rows=100, cols=100, predation_rate=4.0, prey_reproduction_rate=1.0, predator_death_rate=1.0, n_steps=1200), 5 seeds (0–4).

Predation_rate=4.0 is the Sprint 11 canonical positive (not the Mobilia 2007 default of 2.0); established in `tests/test_lv_p11_e2e.py`. n_steps=1200 satisfies the ≥1200-generation DEFINITIVE prerequisite documented in carry-forward #8 (Sprint 11).

**Class C failed regimes:** `epc/phase2a/failed_regimes/p11_lotka_volterra.py` — 10 predator-extinction regimes at predator_death_rate ∈ linspace(2.0, 5.0, 10), predation_rate=2.0, n_steps=400 on 50×50. All regimes at μ ≥ λ = 2.0; predators die at or faster than they reproduce. Prey-only absorbing state reached before step 400; P11's n_species prerequisite (exactly 2 non-zero species) fails → P11 rejects at screening.

```
pattern  overall    syn    cat    fai      d verdict
P11        1.000  1.000  1.000  1.000    inf PASS
```

- All 5 positives: DEFINITIVE, confidence=0.900 (rho_anti < −0.7, fft_peak_to_mean > 12, cohens_d < −1.5)
- Class A synthetic TNR = 1.000 (9/9 evaluated; `time_shuffled` SKIPPED — P11 primary metric rho_anti is time_shuffle_invariant: measures inter-species lag structure, not temporal order)
- Class B catalog TNR = 1.000 (7/7 lattice_2d mates: P1_Schelling, P12_RPS, P13_GH, P14_BTW, P15_GoL, P22_SIR, P27_NowakMay — all correctly rejected at screening)
- Class C failed regimes TNR = 1.000 (10/10 predator-extinction regimes correctly rejected)
- Cohen's d = +∞ (positives score=0.900, all negatives score=0.000)

**Files added:**
- `epc/phase2a/failed_regimes/p11_lotka_volterra.py`
- `analysis/outputs/p11_phase2a_panel.json`

**Files modified:**
- `analysis/run_phase2a_panel.py`: `build_p11_positives()`, `make_p11_detector_fn()`, `run_p11()` added; `p11` case in `main()`; docstring updated
- `epc/phase2a/failed_regimes/__init__.py`: p11_lotka_volterra registered

---

## Part B — Conditional re-runs

**P1 PARTIAL (Sprint 43) examined:** The three open P1 carry-forwards are:
- C-p1-time-shuffle-fp: `time_shuffled` substrate fires because Schelling segregation is temporal-order-independent (segregated structure present in each frame). Resolution requires detector spec revision (invariance flag) or Class A spec revision — not chain-resolved.
- C-p1-linear-gradient-fp: `linear_gradient` fires because Moran's I responds to gradient structure. Resolution requires Class A substrate redesign or P1 screening gate adjustment — not chain-resolved.
- C-p1-class-c-subthreshold-fp: Class C calibration (accidental clustering at density=0.9, small grids). Resolution requires separate calibration work — not chain-resolved.

Per brief §B, no autonomous re-run performed. Carry-forwards remain open.

No other PARTIALs in the specified set (p3, p12, p13, p22, p27 all showed PASS).

---

## Part C — REPLICATION_NOTES

New section appended:
- `## Phase-2a Panel Result (v1.2) — Sprint 45 (P11 Lotka-Volterra, PASS)`

Sprint 30 rule notation (PASS → depth_gap update; no detector/model changes).

---

## Part D — depth_gap.md

- **P11 row:** dim4 PARTIAL → **PASS**; grade remains **GAP** (dim1 PARTIAL: Mobilia-Georgiev-Täuber 2007 cited but no specific Fig/table reproduced with stated tolerance).
- **AT-DEPTH count unchanged: 8 / 19** (P11 dim4 PASS does not advance grade because dim1 is still PARTIAL).
- Aggregate findings block updated: Sprint 45 finding added; dim4 resolved count updated (P11 added to list of patterns with Phase-2a panel PASS).

---

## Part E — Paper sync

- `docs/paper_section4_draft.md` §4.23 added: Sprint 45 P11 panel result.
- `docs/paper_section6_draft.md`: Sprint 45 AT-DEPTH narrative added (count unchanged 8/19; P11 dim4 PASS noted).
- `docs/paper_CHANGELOG.md`: Sprint 45 entry added.

---

## Part F — Carry-forward review

**Carry-forwards closed this sprint:** None.

**Active carry-forwards from prior sprints (unchanged):**
- C-p1-time-shuffle-fp (Sprint 43): P1 `time_shuffled` Class A FP.
- C-p1-linear-gradient-fp (Sprint 43): P1 `linear_gradient` Class A FP.
- C-p1-class-c-subthreshold-fp (Sprint 43): P1 Class C calibration issue at density=0.9.
- C3 / C-p12-wavelength-scaling (Sprint 9 / Sprint 26): P12 λ ∝ √M wavelength scaling not replicated. Sprint 48 target per brief.
- C-class-a-constant-field-trivial-sync (Sprint 35): P9 constant_field degenerate substrate.
- C-p13-positive-tier-screening (Sprint 44): P13 canonicals reach SCREENING (0.500) not CONFIRMATION; informational only.

**Carry-forwards opened this sprint:** None.

---

## Pre-flight / post-flight

- Pre-flight: HEAD verified at `f4e09c1`; `pytest tests/test_orchestration.py tests/test_phase2a_panel.py` confirmed 154 tests passing (11.6s).
- Files created/modified; panel run executed; JSONs written.
- Post-flight: `test_orchestration.py` + `test_phase2a_panel.py` + `test_cross_detection_matrix.py` passing (178 tests in 104s); transfer matrix unchanged (no model/detector/orchestration changes in this sprint).

---

**Decision: GO**

P11: PASS (TNR=1.000, d=+∞). All 5 positives DEFINITIVE. dim4 PARTIAL→PASS. Grade remains GAP (dim1 unresolved). AT-DEPTH count unchanged 8/19. No new carry-forwards. Chain may proceed to Sprint 46 (continuous_2d dim4 batch — P2 + P5 + P6, now possible with Sprint 37 generators).
