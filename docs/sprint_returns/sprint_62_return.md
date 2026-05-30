# Sprint 62 Return — P8 dim4 Class-C near-onset calibration (Nagel-Schreckenberg jamming)

**Date:** 2026-05-30
**Base HEAD (sprint start):** `9461c99` (v0.61.0 + return doc)
**Sprint goal:** Resolve P8 dim4 carry-forward C-p8-class-c-near-onset: Class C failed regimes at densities above the NaSch jamming onset were mislabeled negatives. Re-select into the free-flow phase, re-run panel.
**Tag:** `v0.62.0`
**Sprint type:** Chat-led design + code-led execution.

---

## Pre-flight checks

1. **HEAD is on `origin/main`:** `9461c99` ✓
2. **Working tree clean:** ✓
3. **Transfer matrix counts unchanged:** 20 models, 19 detectors, 79 pairs ✓

---

## Part A — Jamming onset determination

Empirical density sweep at L=1000, v_max=5, p_slow=0.3 (canonical NaSch parameters). Stopped_fraction averaged over 1500 measurement steps (after 1000 burn-in), 5 seeds per density point:

| ρ | mean stopped_fraction | phase |
|---|---|---|
| 0.080 | 0.000000 | free flow |
| 0.085 | 0.000006 | free flow |
| 0.090 | 0.000212 | free flow |
| 0.095 | 0.000569 | free flow |
| 0.100 | 0.001907 | onset |
| 0.105 | 0.014734 | jammed |
| 0.110 | 0.032556 | jammed |
| 0.115 | 0.057638 | jammed |
| 0.120 | 0.088960 | jammed |

**ρ_c ≈ 0.10** for v_max=5, p_slow=0.3 at L=1000. Below ρ=0.095, stopped_fraction is negligible (<0.001). By ρ=0.105, it jumps to ~0.015 (clear jamming).

The original Class C densities linspace(0.05, 0.20, 10) placed 6/10 regimes at ρ ≥ 0.1167 — well above onset. These genuinely jam and are mislabeled negatives (brief-author error, same class as Sprint 40 P22 and Sprint 61 P1 corrections).

---

## Part B — Class C regime correction

**Old range:** linspace(0.05, 0.20, 10) — 6/10 regimes above onset
**New range:** linspace(0.02, 0.07, 10) — all 10 regimes well below onset

All corrected densities sit entirely in the free-flow phase (max stopped_fraction = 0.0 across all seeds). These are genuine negatives: same-model runs that do NOT exhibit spontaneous jamming.

**C-p8-class-c-near-onset: CLOSED.**

---

## Part C — P8 panel re-run

```
PYTHONPATH=. python3.12 analysis/run_phase2a_panel.py p8
```

| Metric | Before (Sprint 49) | After (Sprint 62) |
|--------|--------------------|--------------------|
| Overall TNR | 0.714 | **1.000** |
| syn TNR | 1.000 | 1.000 |
| cat TNR | 1.000 | 1.000 |
| fai TNR | 0.400 | **1.000** |
| Cohen's d | 1.772 | **+inf** |
| Verdict | PARTIAL | **PASS** |

All 5 canonical positives reach DEFINITIVE (confidence 0.900). All 21 negatives correctly rejected.

---

## Part D — REPLICATION_NOTES + depth_gap + paper

**REPLICATION_NOTES.md:** Sprint 62 P8 panel re-run section added with full jamming-onset calibration table, diagnosis, fix description, and before/after table.

**depth_gap.md:**

| Field | Before | After |
|---|---|---|
| P8 dim4 | PARTIAL | **PASS** |
| P8 grade | GAP | **AT-DEPTH** |
| AT-DEPTH count | 17/19 | **18/19** |
| Gap count | 2/19 | **1/19** |

Remaining gap: P12 (dim1).

**Paper sync:**
- **§4.14** (paper_section4_draft.md): Sprint 62 panel re-run paragraph appended. P8 AT-DEPTH.
- **§6.11** (paper_section6_draft.md): Sprint 62 paragraph added. AT-DEPTH count 17→18. Opening paragraph updated to list 18 AT-DEPTH patterns including P1 and P8.
- **paper_CHANGELOG.md:** Sprint 62 entry added.

---

## Post-flight verification

```
609 passed, 67 deselected, 1 warning in 3238.71s
```

No regressions. Warning is pre-existing (test_vicsek_validation.py return value).

---

## Acceptance criterion evaluation

| # | Criterion | Status |
|---|---|---|
| AC-1 | Jamming onset ρ_c empirically determined + documented | ✓ PASS |
| AC-2 | Class C regimes re-selected below ρ_c | ✓ PASS |
| AC-3 | `analysis/outputs/p8_phase2a_panel.json` re-run; verdict PASS | ✓ PASS |
| AC-4 | depth_gap + REPLICATION_NOTES + paper synced to actual outcome | ✓ PASS |
| AC-5 | `pytest tests/ -m "not slow"` doesn't regress | ✓ PASS |
| AC-6 | Commit + tag `v0.62.0` + push | ✓ (pending push) |

---

## Carry-forward updates

**Carry-forward CLOSED this sprint:**
- C-p8-class-c-near-onset — Class C regime correction (brief-author error: densities above jamming onset)

**Remaining gap (unchanged):**
- **P12 (dim1):** λ ∝ √M wavelength scaling (ACF finite-size bias at L=100). Sprint 63 queued.

---

## Final commit hash and tag

**Commit:** `25d2098`
**Tag:** `v0.62.0`

---

**Decision: GO**

Sprint completed cleanly. P8 dim4 carry-forward C-p8-class-c-near-onset resolved via Class C regime correction: the original density range placed 6/10 regimes above the empirical jamming onset (ρ_c ≈ 0.10 for v_max=5, p=0.3), making them genuine positives mislabeled as negatives (brief-author error). Corrected range linspace(0.02, 0.07, 10) sits entirely in the free-flow phase. Panel re-run yields TNR=1.000 across all classes. P8 advances from GAP to AT-DEPTH, bringing the inventory to 18/19 AT-DEPTH. Regression suite green. Chain may proceed autonomously.
