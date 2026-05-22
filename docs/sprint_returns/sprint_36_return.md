# Sprint 36 Return Summary

**Sprint goal:** Clear 9 sprints' accumulated paper-debt by adding §3.7 (Phase-2a
standard negative panel methodology) and syncing §4 per-pattern subsections + §6
aggregate counts to reflect work landed in Sprints 27–35. **Status: complete.**

**Sprint 37 cue:** Chain begins. Sprint 37 brief (already queued) starts the
lattice_2d dim4-closure batch: P22 SIR + P27 Nowak-May. Orchestrator self-dispatches
as soon as this return doc lands.

---

## Pre-flight verification

- Base HEAD: `442b968` (Sprint 36 follow-up — orchestrator/ source commit, one ahead
  of the brief's stated 741a2d4; the follow-up added only `orchestrator/` files,
  non-conflicting with this sprint's paper-only scope). ✓
- `python3 scripts/count_transfer_matrix.py` (with PYTHONPATH=.): **20 / 19 / 79 /
  274 / 27 / 19 / 361 / 77 / 284** ✓ (unchanged).
- pytest pre-flight: skipped per chain-aware skip rule (Sprint 35 post-flight at
  HEAD 741a2d4 passed 585 / 0). ✓

---

## Part A — §3.7 added to `docs/paper_section3_draft.md`

New section "Phase-2a Standard Negative Panel" appended after §3.6 Statistical
Power Requirements. Five subsections:

1. **Motivation** (~150 words): cites Sprint 28 audit finding that dim4 was the
   dominant gap across 15/19 patterns; distinguishes cross-detection matrix
   (tests detectors against catalog neighbors only) from the Phase-2a panel
   (tests against substrate-diverse non-positives).

2. **Panel composition** (~200 words): 30 substrates across Class A (synthetic
   null, fixed across all patterns), Class B (catalog-derived non-positives,
   substrate-typed per v1.1), Class C (failed-regime, with N/A escape hatch).

3. **PASS criterion + invariance flags** (~200 words): overall TNR ≥ 0.95,
   per-class TNR ≥ 0.90 (advisory weak-class), Cohen's d ≥ 1.0; v1.2
   `primary_metric_permutation_invariant` and
   `primary_metric_time_shuffle_invariant` flags for degenerate Class A skips.

4. **Spec evolution as methodology contribution** (~200 words): honest narrative
   of v1.0 (Sprint 28, cross-format adapter contamination + substrate-type
   conflation) → v1.1 (Sprint 30, both fixed) → v1.2 (Sprint 34, invariance
   flags for degenerate Class A substrates); framed as emergent methodology, not
   designed up-front.

5. **Reproducibility note** (~100 words): spec at `docs/phase2a_panel_spec.md`,
   harness at `epc/phase2a/`, archived results, per-pattern JSONs.

Acceptance check: `grep -c "^## 3.7"` → 1 ✓; `grep -c "Phase-2a standard
negative panel"` → 2 ✓.

---

## Part B — §4 per-pattern Phase-2a paragraphs + Sprint 27 deferred §4.19 update

### §4.3 (P15 Game of Life)

**Phase-2a panel (v1.1, Sprint 33): PASS.** Overall TNR = 1.000, Cohen's d =
8.282. All three classes clean. Panel result provides independent confirmation
of AT-DEPTH grade. Paragraph appended before `## 4.4`.

### §4.5 (P9 Kuramoto / P9 synchronization)

**Phase-2a panel (v1.0 PARTIAL → v1.1 PARTIAL → v1.2 PASS-with-weakness).**
v1.2 result: overall TNR = 0.952, Cohen's d = 4.781. Class A TNR = 0.875
(weak — `constant_field` carry-forward). Class B = 1.000. Class C = 1.000.
Multi-paragraph history of spec evolution appended before `## 4.6`, serving
as the §3.7 worked example.

### §4.7 (P14 BTW Sandpile)

**Phase-2a panel (v1.2, Sprint 35): PASS.** Overall TNR = 0.960, Cohen's d =
10.585. Class A clean, Class B (lattice_2d) 7/7 clean, Class C (dissipative)
0.900 with one borderline at p_diss=0.350 (carry-forward C-p14-class-c-
borderline). Paragraph appended before `## 4.8`.

### §4.19 (P10 chimera / Sprint 27 deferred Phase 1n update)

Sprint 27's multi-seed (A, β) phase-boundary scan (Phase 1n) was deferred from
paper text in Sprint 27 itself. Added as "Sprint 27 update" block appended to
§4.19 before `## 4.20`:

- Phase 1n parameters: N=128, 5 seeds × 6 A values × 5 β values = 150 runs.
- Key table of basin fractions (β=0.05: 30/30; β=0.10: 29/30; β=0.18: 18/30;
  β=0.20: 12/30; β=0.22: 1/30).
- Key finding: Sprint 26's apparent sharp transition was a single-seed artifact;
  the true boundary is a smooth basin-volume gradient.
- Reconciliation with Sprint 26 Phase 1k at β=0.18.
- References `analysis/outputs/p10_phase_boundary_multiseed.json` and
  `analysis/outputs/p10_basin_volume_multiseed.png`.

Note: the time-units methods note (Sprint 27 carry-forward #34) was already
present in §4.19 as "Honest caveat on time-unit convention" from a prior draft
pass; not duplicated.

### §4.20 (P18 voter model)

**Phase-2a panel (v1.1, Sprint 33): PASS.** Overall TNR = 1.000, Cohen's d =
+∞ (no overlap between canonical-positive and pooled-negative score
distributions). All three classes clean. Panel result provides second
independent confirmation of AT-DEPTH grade (alongside §6.10 metric-level
discriminators). Paragraph appended at end of file.

Acceptance check: `grep -c "Phase-2a panel" docs/paper_section4_draft.md` → 7
(≥5 required) ✓.

---

## Part C — §6 aggregate updates

New subsection `## 6.11 Aggregate Grading Status` appended at end of
`docs/paper_section6_draft.md`:

- **AT-DEPTH: 5 / 19** (P9, P15, P18, P28, P31).
- **GAP: 14 / 19**.
- One-line note: 4 of 5 AT-DEPTH patterns (P9, P15, P18) reached that grade via
  Phase-2a panel; P28 and P31 were AT-DEPTH pre-panel via mechanistic-null gates.
- Transfer-matrix figures confirmed locked: **20 / 19 / 79 / 274 / 27 / 19 / 361
  / 77 / 284** (unchanged from Sprint 26).

Acceptance checks: `grep -c "AT-DEPTH"` → 4 (≥1) ✓; `grep -c "5 / 19"` → 1 ✓.

---

## Part D — `docs/paper_CHANGELOG.md` created

New file at `docs/paper_CHANGELOG.md` with:
- Header comment explaining purpose and per-sprint format.
- Sprint 36 entry listing all 7 changes (§3.7, §4.3, §4.5, §4.7, §4.19, §4.20,
  §6.11 additions).
- Sprint 37 placeholder entry.

---

## Post-flight verification

All acceptance criteria met:

| Check | Result |
|---|---|
| `grep -c "^## 3.7" docs/paper_section3_draft.md` | **1** ✓ |
| `grep -c "Phase-2a standard negative panel" docs/paper_section3_draft.md` | **2** ✓ |
| `grep -c "Phase-2a panel" docs/paper_section4_draft.md` | **7** (≥5) ✓ |
| `grep -c "AT-DEPTH" docs/paper_section6_draft.md` | **4** (≥1) ✓ |
| `"5 / 19"` present in `docs/paper_section6_draft.md` | **1** ✓ |
| `test -f docs/paper_CHANGELOG.md` | **EXISTS** ✓ |
| `python3 scripts/count_transfer_matrix.py` | **20/19/79/274/27/19/361/77/284** ✓ |
| `pytest tests/ -m "not slow"` | **585 passed, 0 failed, 65 deselected** ✓ |

`git diff --stat` shows 4 files: `docs/paper_section3_draft.md` (+116),
`docs/paper_section4_draft.md` (+115), `docs/paper_section6_draft.md` (+24),
`docs/paper_CHANGELOG.md` (new, +20). No other files modified. ✓

---

## Deviations and judgment calls

### Deviation 1 — Base HEAD is 442b968, not 741a2d4

The brief's stated base HEAD was `741a2d4` (Sprint 35, v0.35.0). Actual HEAD was
`442b968` (Sprint 36 follow-up: commit orchestrator/ as version-tracked source of
truth). The follow-up commit added only `orchestrator/` files — no paper, code, or
test changes — so it is fully non-conflicting with this sprint's paper-only scope.
No action required.

### Deviation 2 — §4.19 time-units methods note already present

Sprint 27's carry-forward #34 ("add a one-paragraph methods note to §4.19 capturing
the time-units reconciliation") was found to already be present in §4.19 as the
"Honest caveat on time-unit convention" paragraph. Not duplicated. Only the Phase 1n
multi-seed phase boundary content (carry-forward #31) was missing and was added.

### Deviation 3 — Test environment requires explicit PYTHONPATH and local pip install

System Python3 on this machine lacks numpy/scipy/pytest by default. Installed pytest,
numpy, scipy via `pip3 install --break-system-packages`. Test run required
`PYTHONPATH=/home/matthewhmaxwell/emergent-pattern-catalog`. Test suite completed
in ~60 minutes wall time (585 tests, no slow tests, but many tests are
compute-intensive even without the `slow` mark). All 585 passed.

---

## Carry-forward summary (unchanged from Sprint 35)

- **C-class-a-constant-field-trivial-sync (Sprint 35):** P9's residual false positive
  from `constant_field` substrate. Uniformity-sensitive detectors only. v1.3
  candidate for chat-led review; not in scope for Sprint 37 lattice_2d batch
  (grid-based detectors look for spatial structure, not uniformity).
- **C-p27-time-shuffle-invariance (Sprint 34):** Still provisional. Sprint 37's P27
  panel run will validate or invalidate.
- **C-p14-class-c-borderline (Sprint 33):** Persists at p_diss=0.350. P14 still PASS.
- **C-supplements (Sprint 31):** Still OPEN for lattice_1d / lattice_2d_continuous /
  scalar_wealth / opinion_space. Not needed for Sprint 37's P22+P27 batch.
- **C-pyproject-pin (Sprint 29):** Still OPEN. 1-line pyproject.toml change deferred.

---

## HEAD commit hash and tag at end of sprint

- **Commit:** `e005295`
- **Tag:** `v0.36.0` (pushed to origin)
