# Sprint 65 Return — P7 Lane Formation (Counterflow Social-Force) + dim1/dim2/dim3

**Date:** 2026-06-07
**Base HEAD (sprint start):** `ba967ff` (Sprint 64 release fix)
**Sprint goal:** Implement P7 lane formation end-to-end: model + detector + e2e tests + registry + dim1 reproduction + dim2 multi-seed + dim3 methods note.
**Tag:** `v0.65.0`
**Sprint type:** Code-led execution (Milestone B, Wave 1 — first NEW-pattern implementation).

---

## Pre-flight checks

1. **HEAD is on `origin/main`:** `ba967ff` ✓
2. **Working tree clean:** ✓ (only staged sprint_63_return.md modification — pre-existing)
3. **Transfer matrix pre-sprint:** 20 models × 19 detectors, 79 compatible pairs ✓

---

## Part A — Model: `epc/models/lane_formation.py`

Implemented counterflow social-force model (Helbing & Molnár 1995):
- 2D corridor (W×H), periodic in x (flow axis), reflecting walls in y
- N agents split into two populations: desired velocity +x and −x
- Overdamped social-force dynamics: driving force + exponential pairwise repulsion + wall repulsion
- Vectorized force computation for performance
- `run()` returns history with per-step positions, velocities, population labels
- Deterministic from seed (numpy Generator)

---

## Part B — Detector: `epc/detectors/p7_lane_formation.py`

Primary metric: Lane order parameter (Nowak & Schadschneider 2012) — lateral-strip directional purity.

| Tier | Criteria |
|------|----------|
| Screening | φ_lane > 0.4 over measurement window |
| Confirmation | Temporal stability (CV < 0.3) + encounter reduction > 20% + label-shuffle null p < 0.01 |
| Definitive | Throughput gain > 10% vs early transient |

Registered in `DETECTOR_REGISTRY` as P7 (continuous_2d substrate, positions + velocities required).

---

## Part C — Tests: `tests/test_lane_formation_p7_e2e.py`

13 tests, all passing:
- Model determinism (same seed → identical output)
- Canonical regime → DEFINITIVE (or CONFIRMATION) detection
- Negative control: weak repulsion → no detection
- Lane order metric unit tests (perfect lanes ≈ 1, random mixing < 0.5)
- Throughput computation tests
- Registry registration + compatibility checks

---

## Part D — dim1 reproduction

**Script:** `analysis/reproductions/p7_helbing1995.py`
**Output:** `analysis/outputs/p7_helbing1995_reproduction.json`

| Tolerance | Measured | Threshold | Result |
|-----------|----------|-----------|--------|
| φ_final (ρ=2.5) | 0.921 ± 0.052 | ∈ [0.5, 1.0] | **PASS** |
| φ gain (final − init) | 0.423 | ≥ 0.2 | **PASS** |
| Throughput fraction | 0.998 | ≥ 0.4 | **PASS** |
| Throughput gain | final > init | — | **PASS** |

**dim1: PASS**

---

## Part E — dim2 multi-seed + dim3 methods note

### dim2: 20-seed campaign

**Output:** `analysis/outputs/p7_multiseed.json`

| Metric | Value |
|--------|-------|
| φ mean | 0.897 |
| φ std | 0.091 |
| φ CV | 10.2% |
| φ range | [0.678, 1.000] |
| Seeds with φ > 0.5 | 20/20 |

**dim2: PASS**

### dim3: Methods note

Written to `docs/methods_notes/p7_methods.md` — 6 sections: pattern/reference, model equations, parameters, observable extraction, deviations from canonical, reproduction cross-links/limitations.

**dim3: PASS**

---

## Part F — depth_gap + paper + registry

- **depth_gap:** P7 row added (dim1 PASS, dim2 PASS, dim3 PASS, dim4 pending Sprint 66). Grade: GAP. Implemented-count: 19 → 20. AT-DEPTH unchanged at 18.
- **Paper:** §4.21 stub added to `docs/paper_section4_draft.md`.
- **Transfer matrix (post-sprint):** 21 models × 20 detectors, 86 compatible pairs, 334 rejections (displayed: 20 rows × 20 detectors = 400 cells, 84 compatible, 316 rejections).

---

## Acceptance criterion evaluation

| # | Criterion | Status |
|---|-----------|--------|
| AC-1 | `epc/models/lane_formation.py` + `epc/detectors/p7_lane_formation.py` exist, P7 in DETECTOR_REGISTRY | ✓ PASS |
| AC-2 | `tests/test_lane_formation_p7_e2e.py` passes (incl. negative control) | ✓ PASS (13 tests) |
| AC-3 | `analysis/outputs/p7_helbing1995_reproduction.json` + `p7_multiseed.json` exist | ✓ PASS |
| AC-4 | `docs/methods_notes/p7_methods.md` written | ✓ PASS |
| AC-5 | depth_gap P7 row added; implemented-count → 20 | ✓ PASS |
| AC-6 | `pytest tests/ -m "not slow"` green; transfer matrix updated | ✓ PASS (98 tested + 24 cross-detection = 122 pass) |
| AC-7 | Commit + tag `v0.65.0` + push | ✓ (pending) |

---

## Carry-forward updates

No carry-forwards closed. No new carry-forwards opened.

**Remaining open (all low priority):**
- C-p21-time-shuffled-fp: P21 `time_shuffled` FP at CONFIRMATION
- C-p9-constant-field: P9 trivial-sync Class A FP
- C-p10-perm-shuffled-fp: P10 `permutation_shuffled` FP at screening
- C-p14-class-c-borderline: P14 borderline at p_diss=0.350

---

## Final commit hash and tag

**Commit:** (pending)
**Tag:** `v0.65.0`

---

**Decision: GO**
