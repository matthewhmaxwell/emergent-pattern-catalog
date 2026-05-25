# Sprint 47 Return — lattice_1d + oscillator dim4 batch: P8 NS + P10 chimera

**Date:** 2026-05-25
**Base HEAD (sprint start):** `44313f9`
**Sprint goal:** Run v1.2 panel against P8 (Nagel-Schreckenberg, lattice_1d) and P10 (chimera, oscillator). No detector/model/spec changes (Sprint 30 rule).
**Tag:** `v0.47.0`
**Commit:** `1f06afc`

---

## Part A — P8 panel (Nagel-Schreckenberg / traffic jamming)

**Canonical positive:** NagelSchreckenberg(L=1000, density=0.30, v_max=5, p_slow=0.3, n_steps=2500, seeds 0–4). Detector uses burn_in=1000 → 1500 measurement steps. stopped_fraction ≈ 0.18–0.22 at rho=0.30.

**Class C failed regimes:** `epc/phase2a/failed_regimes/p8_nagel_schreckenberg.py` — 10 low-density regimes, rho ∈ linspace(0.05, 0.20, 10), L=1000, p_slow=0.3, n_steps=2500, seeds 200–209.

```
pattern  overall    syn    cat    fai      d verdict
P8         0.652  0.800  1.000  0.400  1.751 PARTIAL
```

- All 5 positives: DEFINITIVE, confidence=0.900
- Class A synthetic TNR = 0.800 (8/10; 2 FPs)
  - `permutation_shuffled` FP at SCREENING (score=0.500) — stopped_fraction is spatial-order-invariant
  - `time_shuffled` FP at SCREENING (score=0.500) — stopped_fraction is time-average-invariant
  - P8 absent from `epc/phase2a/detector_invariance.py`; flags not auto-set per brief
- Class B catalog + supplements TNR = 1.000 (P31_zhang_sorting + 2 supps)
- Class C failed regimes TNR = 0.400 (4/10 correctly rejected; 6 FPs at rho ≥ 0.1167)
  - rho=0.1167: CONFIRMATION (0.700); rho ∈ {0.1333–0.2000}: DEFINITIVE (0.900)
  - NS jamming onset at p=0.3 is rho ≈ 0.12 (Nagel & Schreckenberg 1992)
- Cohen's d = 1.751

**Infrastructure added (Sprint 30 scope: panel plumbing only):**
- `epc/phase2a/synthetic.py`: "sequence" format branches added to all 8 Class A generators + `permutation_shuffled_positive`. Added `SEQUENCE_DEFAULT_N = 200` and `_sequence_history_from_array()` helper.
- `epc/phase2a/catalog.py`: `_adapt_to_sequence()` function + "sequence" case in `load_catalog_substrate_for_format`.
- `epc/phase2a/panel.py`: "sequence" case in `class_a_kwargs` (sets `n = target_n`).

**Decision (PARTIAL → escalate):** Three carry-forwards opened:
- **C-p8-perm-shuffled-fp**: `permutation_shuffled` Class A fires at SCREENING. To close: set `permutation_invariant=True` in `detector_invariance.py` (requires invariance analysis sprint, not in scope here).
- **C-p8-time-shuffle-fp**: `time_shuffled` Class A fires at SCREENING. Same resolution path.
- **C-p8-class-c-near-onset**: Restrict next Class C sweep to rho ∈ linspace(0.01, 0.09, 10) to stay cleanly below the NS jamming transition at rho ≈ 0.12 (p=0.3). Recommended effort: S (parameter change + one-line re-run).

**Files added/modified:**
- `epc/phase2a/failed_regimes/p8_nagel_schreckenberg.py` (new)
- `epc/phase2a/synthetic.py` (modified: "sequence" format)
- `epc/phase2a/catalog.py` (modified: _adapt_to_sequence)
- `epc/phase2a/panel.py` (modified: sequence class_a_kwargs)
- `analysis/run_phase2a_panel.py` (modified: build_p8_positives, make_p8_detector_fn, run_p8, main case)
- `analysis/outputs/p8_phase2a_panel.json` (new, gitignored)

---

## Part B — P10 panel (chimera states / non-local Kuramoto ring)

**Canonical positive:** KuramotoNonlocal(N=128, A=0.995, beta=0.05, init_mode="asymmetric_gaussian", seeds 0–4), n_frames=50. Canonical Abrams-Strogatz (2004) parameters. metadata["has_nonlocal_coupling"]=True set by get_metadata().

**Class C failed regimes:** `epc/phase2a/failed_regimes/p10_chimera.py` — 10 ordinary all-to-all Kuramoto regimes, K ∈ linspace(1.5·K_c, 4.0·K_c) = linspace(1.5, 4.0, 10), gamma=0.5, N=128, dt=0.05, n_steps=6000, record_every=10, seeds 300–309.

```
pattern  overall    syn    cat    fai      d verdict
P10        0.957  0.900  1.000  1.000  9.679 PASS
```

- All 5 positives: DEFINITIVE, confidence=0.950
- Class A synthetic TNR = 0.900 (9/10; 1 FP)
  - `permutation_shuffled` FP at SCREENING (score=0.500) — permuting oscillator indices preserves global velocity distribution → pos_vel_ac reads similarly to chimera positive
  - P10 present in invariance dict as `(False, False)`; carry-forward C-p10-perm-shuffled-fp per brief
- Class B catalog + supplements TNR = 1.000 (P9_kuramoto + 2 supps: incoherent_phases, subcritical_kuramoto)
- Class C failed regimes TNR = 1.000 (10/10 correctly rejected)
  - All-to-all Kuramoto above K_c produces full synchronisation: no coexistence, pos_vel_ac[lag=4] << 0.55 → rejected at screening
  - Additionally has_nonlocal_coupling absent from Class C metadata → DEFINITIVE gate locked out
- Cohen's d = 9.679

**Decision (PASS):** P10 dim4 advances PARTIAL→PASS. All four dimensions PASS → **P10 advances to AT-DEPTH**. AT-DEPTH count: **10 / 19**.

Carry-forward opened:
- **C-p10-perm-shuffled-fp**: `permutation_shuffled` Class A fires at SCREENING. P10 is in invariance dict as (False, False) but permutation of oscillator spatial labels does preserve some chimera statistics. To close properly requires either setting `permutation_invariant=True` (would SKIP the substrate, not penalize) or hardening the coexistence gate to be permutation-sensitive.

**Files added/modified:**
- `epc/phase2a/failed_regimes/p10_chimera.py` (new)
- `analysis/run_phase2a_panel.py` (modified: build_p10_positives, make_p10_detector_fn, run_p10, main case)
- `analysis/outputs/p10_phase2a_panel.json` (new, gitignored)

---

## Part C — depth_gap.md update

| Pattern | dim4 before | dim4 after | grade before | grade after | AT-DEPTH Δ |
|---------|-------------|------------|--------------|-------------|------------|
| P8      | PARTIAL     | PARTIAL    | GAP          | GAP         | 0          |
| P10     | PARTIAL     | PASS       | GAP          | AT-DEPTH    | +1         |

AT-DEPTH count: 9 → **10 / 19**. Gap count: 10 → **9 / 19**.

---

## Part D — documentation sync

- `REPLICATION_NOTES.md`: Two new Phase-2a panel result sections (P8 PARTIAL, P10 PASS) appended.
- `docs/paper_section4_draft.md`: §4.25 added (Sprint 47 P8/P10 panel results).
- `docs/paper_section6_draft.md`: Sprint 47 AT-DEPTH +1 paragraph appended.
- `docs/paper_CHANGELOG.md`: Sprint 47 entry appended.

---

## Part E — post-flight

```
PYTHONPATH=. python3.12 -m pytest tests/test_orchestration.py tests/test_phase2a_panel.py -x -q
154 passed in 9.41s
```

No regressions. Test count unchanged (no new tests added; sprint 30 rule in force).

---

## Part F — open carry-forwards (sprint 47)

| ID | Pattern | Description | Priority | Recommended effort |
|----|---------|-------------|----------|--------------------|
| C-p8-perm-shuffled-fp | P8 | `permutation_shuffled` Class A FP at SCREENING; stopped_fraction spatial-order-invariant | Low | M (invariance analysis sprint) |
| C-p8-time-shuffle-fp | P8 | `time_shuffled` Class A FP at SCREENING; stopped_fraction time-average-invariant | Low | M (same as above) |
| C-p8-class-c-near-onset | P8 | 6/10 Class C regimes at rho ≥ 0.1167 fire — overlaps NS jamming onset; restrict next sweep to linspace(0.01, 0.09, 10) | Medium | S (one-sprint parameter change + re-run) |
| C-p10-perm-shuffled-fp | P10 | `permutation_shuffled` Class A FP at SCREENING; oscillator spatial relabelling preserves velocity distribution | Low | M (hardening coexistence gate or invariance flag) |

---

## Decision line

- **P8:** PARTIAL → escalate (carry-forwards C-p8-class-c-near-onset, C-p8-perm-shuffled-fp, C-p8-time-shuffle-fp). Dim4 remains PARTIAL. Recommended next action: restrict Class C density range and re-run (S effort).
- **P10:** PASS → AT-DEPTH. Sprint goal met. Carry-forward C-p10-perm-shuffled-fp is low priority.
- **Overall:** AT-DEPTH count advances 9→10. Sprint 47 is GO with C-p8-class-c-near-onset flagged for escalation.
