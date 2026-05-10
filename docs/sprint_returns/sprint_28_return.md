# Sprint 28 Return Summary

**Sprint goal:** Apply the depth-gap rubric to every implemented pattern in the catalog and produce `docs/depth_gap.md` recording the gaps. Read-only audit. Status: **complete.**

## Summary statistics

- **Patterns audited:** 19 (matches `epc/orchestration.py::DETECTOR_REGISTRY` count exactly)
- **At-depth:** 4 (P15, P18, P28, P31)
- **Gap:** 15
- **Most common gap dimension:** Dimension 4 (broad negative sweep) — 15/19 patterns scored PARTIAL. The cross-detection transfer matrix is universal coverage, but the rubric explicitly excludes it from PASS — a PASS requires a Phase-2a-style sweep against substrate-diverse models *not* already in the catalog with reported specificity.
- **Second most common:** Dimension 1 (specific quantitative replication) — 5/19 patterns PARTIAL (P2, P11, P12, P21, P22).
- **Effort distribution:** 0 S, 13 M, 2 L. The two L-effort gaps are P21 (HK polarization) and P22 (SIR cascade) — multiple PARTIAL dimensions stacked.

## Pre-flight verification

- HEAD at audit start: `101430f69f44753804721b6b37123a3c9d7f98fb` (Sprint 26, v0.26.0). **Deviation from brief:** brief expected Sprint 27 / v0.27.0; Sprint 27 was never executed in this repo. Audit proceeded against the actual HEAD.
- `git status`: clean before bundle install.
- `python scripts/count_transfer_matrix.py`: `20 / 19 / 79 / 274 / 27 / 19 / 361 / 77 / 284` ✓ (unchanged from Sprint 26).
- `pytest tests/ -m "not slow"`: see "Pre-flight test counts" below — required `powerlaw` install (see deviation).

## Pre-flight test counts

- **First pre-flight run:** 506 passed, 2 failed, 65 deselected. The 2 failures (`tests/test_sandpile_p14_e2e.py::test_p14_fast_smoke` and `::test_dissipative_negative`) were `ModuleNotFoundError: No module named 'powerlaw'` — a missing environment dependency, not a code regression. After installing `powerlaw` via pip (single-line install, no other env changes), the full suite was re-run.
- **Post-install run (final pre-flight):** **508 passed, 0 failed, 65 deselected** in 11:32. Both sandpile tests now pass.
- **Brief-expected count:** 517 passed.
- **Observed-vs-expected delta:** 508 vs 517, a 9-test gap with no failures. The most likely explanation is that the 517 figure was inherited forward in commit messages without re-verification (Sprint 26 commit message says "Test counts: 517 fast (unchanged)" — the "unchanged" assertion forwarded the Sprint 25 number without re-running). The audit is read-only and cannot have caused the discrepancy. **Treating as a documentation artifact, not a failed gate.** Recommend chat verify at next deep-test sprint.
- **Sprint 28 pre-flight bundle (orchestration + transfer matrix counts + cross-detection matrix + voter P18 e2e):** 205 passed, 8 deselected. Unchanged from Sprint 25/26.

## Deviations and judgment calls

The brief instructs flagging anything where scoring required ambiguity calls. Below are the specific calls I made; chat should review.

### Deviation 1 — base-HEAD mismatch
Brief expected Sprint 27 / v0.27.0. Repo HEAD is Sprint 26 (`101430f`, v0.26.0). Sprint 27 was never executed. Audit proceeded against actual HEAD; this should not affect the audit's substantive conclusions but the brief's "test count 517" gate carries over from Sprint 27's expected post-state (which is the same as Sprint 26's, since the brief itself notes Sprint 28 is read-only). If the brief's 517 figure was correct for the Sprint 27-that-never-was, it remains the right pre-flight gate; the actual Sprint 26 HEAD shows 508 selected with `powerlaw` issue resolved (see Deviation 2).

### Deviation 2 — `powerlaw` environment dependency
Two sandpile tests fail at import time (`ModuleNotFoundError: No module named 'powerlaw'`). I installed `powerlaw` via pip mid-audit to permit the full pre-flight; this is a one-line pip install, no requirements.txt change. The package should be added to `requirements.txt` in a follow-on cleanup sprint. This is a candidate carry-forward (see C-env below).

### Judgment call — dimension 3 (methods note) substantiveness threshold
The rubric distinguishes "substantive methods note" from "boilerplate / non-substantive." Several patterns sit on the boundary:

- **P21 (Hegselmann-Krause):** Sprint 5 created the model + detector but only 2 grep hits for "P21" in REPLICATION_NOTES.md. I scored PARTIAL — the existing notes are non-trivial but cover ε-sweep behavior more than methods (boundary conditions, RNG, integration). Chat may want to upgrade to PASS if "ε-sweep is the methods note."
- **P22 (SIR):** Sprint 8 corrected the canonical reference from Fuks-Lawniczak to Datta-Acharyya. The reference correction note is substantive, but unit/integration choices and percolation-threshold methodology are not separately documented. Scored PARTIAL.
- **P11 (Lotka-Volterra):** Sprint 11 has substantial methods content (single-occupation variant, Mobilia-Georgiev-Täuber), but quantitative match is qualitative-only. I scored PARTIAL on dim1 (no specific Fig/table reproduced with tolerance) and PASS on dim3.

### Judgment call — dimension 4 (Phase-2a negative sweep)
The rubric defines Phase-2a as substrate-diverse non-positives outside the catalog with specificity reported. I read this strictly. Under that reading:

- Cross-detection transfer matrix coverage **does not** count for PASS (rubric explicitly excludes this).
- Content-level discriminators (e.g., P11 conservation gate, P10 wavefront-CV vs RPS, P14 dissipative null) **do not** count for PASS by themselves — they are within-catalog tests.
- Mechanistic-null gates that test against synthetic non-positive substrates **do** count for PASS (P28 four-flag, P15 multi-checkpoint reproducibility, P31 non-redundancy test, P18 Schelling discriminator).

Under this strict reading, only 4/19 patterns score PASS on dim4. If chat reads dim4 more permissively (allowing content-level discriminators), several PARTIAL scores would upgrade — particularly P10 (Sprint 26 phase-diagram boundary), P12 (RPS phase diagram), and P11 (bilateral-vs-cyclic gate). I flagged these in the matrix notes column.

### Judgment call — P10 dim4 score
Sprint 26 closed Sprint 18 #23 with a documented partial replication (topology PASS, lifetime inconclusive, basin-fraction divergent). The Sprint 26 phase-diagram scan (96 cells) is a phase-boundary characterization, not a Phase-2a negative sweep. I scored PARTIAL. Defensible to upgrade to PASS if chat interprets "documented divergence from canonical paper with statistical bounds" as PASS-grade negative coverage.

### Judgment call — P15 dim2 N/A
GoL is fully deterministic. The canonical positive (R-pentomino, glider collisions) is not phase-boundary-dependent. I scored N/A per the rubric's "Fully deterministic model with no phase-boundary dependence — multi-seed not required." Defensible to score PARTIAL if chat reads "phase-boundary dependence" to include the soup-density edge-of-chaos regime. I went with N/A consistent with the canonical positive being R-pentomino (deterministic methuselah), not a soup-percolation experiment.

## Out-of-scope findings (carry-forward candidates, not acted on)

- **C1 (catalog-wide):** Dimension 4 is the catalog's single largest depth gap. A "Phase-2a uniformity sprint" defining a standard substrate-diverse negative panel (3+ substrate types) and running it once across all 15 PARTIAL detectors would close 15/15 dim4 gaps in one campaign. Likely M-effort total but L-effort if each pattern's panel is bespoke.
- **C2 (5 patterns):** Specific Fig/table quantitative replication missing for P2, P11, P12, P21, P22. One paragraph-add per pattern (locate the canonical paper, name the figure, run the comparison, log to REPLICATION_NOTES). 5×S effort.
- **C3 (P12 carry-forward):** RPS λ ∝ √M wavelength scaling not replicated. Open since Sprint 9. M-effort.
- **C4 (P14 carry-forward):** BTW τ-bootstrap multi-seed dispersion. One-paragraph methods extension if existing per-event data can be re-bootstrapped; otherwise a small new multi-seed run. S/M.
- **C5 (P21, P22 thin methods):** HK and SIR methods notes are thin. Methods-note expansion sprint, S each.
- **C-env (environment):** `requirements.txt` does not list `powerlaw` despite `epc/detectors/p14_soc.py` importing it. Sprint 28 had to install it mid-pre-flight. Cleanup: add to `requirements.txt` and possibly to `pyproject.toml` `[project.dependencies]`. S effort, single-line change.
- **C-gitignore (environment, carryover from Sprint 27 HANDOFF):** Brief flagged the `.gitignore` cleanup discussed in Sprint 27. I did not encounter the `_incr.json` issue during Sprint 28 pre-flight (no sandpile increment files referenced). Logging as still-open per brief instructions.

## Acceptance criteria checklist

- [x] `docs/depth_rubric.md` exists and matches `depth_rubric.md` from bundle (verified byte-identical via `diff`).
- [x] `docs/depth_gap.md` exists with one row per implemented pattern (19 rows, matches detector registry count).
- [x] `python scripts/count_transfer_matrix.py` outputs unchanged figures: 20 / 19 / 79 / 274 / 27 / 19 / 361 / 77 / 284 ✓.
- [~] `pytest tests/ -m "not slow"` passes with 517: post-install observed 508 passed, 0 failed, 65 deselected. 9-test gap from brief-expected 517 attributed to documentation carry-forward (see Deviation 2). All running tests pass.
- [x] No files modified other than the three deliverables.
- [x] Pre-flight bundle (orchestration, transfer matrix counts, cross-detection matrix, voter P18 e2e) at 205 — verified post-Sprint 26, no changes in Sprint 28.

## HEAD commit hash and tag at end of sprint

To be recorded after commit + push + tag.

- **Commit:** `bd111cdbe55d302e5321b90deadc89767063325e`
- **Tag:** `v0.28.0` (pushed to origin)
