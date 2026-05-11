# Sprint 34 Return Summary

**Sprint goal:** Land Phase-2a panel spec v1.2 (primary-metric invariance flags + skip-degenerate-by-construction substrates), archive v1.1 spec and the four v1.1 panel JSONs, add invariance flags to the panel-harness configuration, update the harness to read flags + skip + extend the JSON schema, and ratify P15's canonical positive change. **No new panel runs; no detector logic changes.** Status: **complete.** Closes Sprint 32's C-class-a-oscillator-degenerate / Sprint 33's C-class-a-permutation-degenerate (same finding, cross-format confirmed).

## Pre-flight verification

- Base HEAD: `c2c5774` (Sprint 33, v0.33.0). Working tree clean. ✓
- `python3.12 scripts/count_transfer_matrix.py`: **20 / 19 / 79 / 274 / 27 / 19 / 361 / 77 / 284** ✓ (unchanged).
- Pre-flight fast tests: 566 passed (Sprint 33 ground truth).

## Part A — v1.1 archived (5 `git mv` operations)

- `docs/phase2a_panel_spec.md` → `docs/archive/phase2a_panel_spec_v1_1.md`
- `analysis/outputs/p9_phase2a_panel.json` → `analysis/outputs/archive/v1_1/p9_phase2a_panel.json`
- `analysis/outputs/p14_phase2a_panel.json` → `analysis/outputs/archive/v1_1/p14_phase2a_panel.json`
- `analysis/outputs/p15_phase2a_panel.json` → `analysis/outputs/archive/v1_1/p15_phase2a_panel.json`
- `analysis/outputs/p18_phase2a_panel.json` → `analysis/outputs/archive/v1_1/p18_phase2a_panel.json`

## Part B — v1.2 spec landed

- `docs/phase2a_panel_spec.md` (byte-identical to bundle's `phase2a_panel_spec_v1_2.md`).

## Part C — Detector invariance flags: **Approach 2 (centralized file)**

**Approach chosen:** Approach 2 — single centralized file `epc/phase2a/detector_invariance.py`.

**Rationale.**
- The detector modules in `epc/detectors/` have heterogeneous structures: some are `BaseDetector` subclasses (P1, P18, P28), some are bare functions (P9 `detect_p9`, P14 `detect_p14`), some are classes that include their own `step_fn` machinery (P15). Adding module-level constants of the form `PRIMARY_METRIC_PERMUTATION_INVARIANT = …` would be a new convention not present elsewhere in those modules and would require 19 independent edits.
- The invariance flags are **panel-harness metadata only** — `epc/phase2a/panel.py` is the only consumer. They are not part of the detector's own contract or its return type, and callers of the detector outside the panel (e.g., the e2e tests in `tests/test_*_e2e.py`) never need them.
- A single centralized file is **greppable from one place** when audits or future v1.x revisions retune the table. It's also the natural home for the `InvarianceFlags` dataclass, helper lookups, and the "P27 provisional" caveat carried as a comment.

**Module:** `epc/phase2a/detector_invariance.py`. Exports `InvarianceFlags` dataclass, `DETECTOR_INVARIANCE_FLAGS: Dict[str, InvarianceFlags]`, and three lookup helpers: `get_flags(pattern_id)`, `is_permutation_invariant(pattern_id)`, `is_time_shuffle_invariant(pattern_id)`. Default for unknown patterns is `(False, False)` — the conservative choice that runs every substrate.

**Flag values applied** (per v1.2 spec §"Change 2" table):

| Pattern | perm_inv | time_inv | Notes |
|---|---|---|---|
| P1  | F | F | |
| P3  | F | F | |
| P5  | T | F | Heading order parameter |⟨e^iθ⟩| is aggregate over headings. |
| P6  | F | F | |
| P9  | T | T | Kuramoto r is aggregate over phases, final-state. |
| P10 | F | F | |
| P11 | F | T | Trajectory order matters; spatial doesn't (well-mixed). |
| P12 | F | F | |
| P13 | F | F | |
| P14 | T | T | Avalanche-size power-law exponent is order-free. |
| P15 | F | F | TE depends on adjacency + ordering. |
| P17 | F | F | (not implemented; declarative-only) |
| P18 | T | T | Consensus fraction max f_k is aggregate, final-state. |
| P19 | F | F | (not implemented; declarative-only) |
| P21 | T | T | Dip test on opinion distribution; distributional + final-state. |
| P22 | F | F | |
| P27 | F | T | **Provisional** per spec Change 2 — see carry-forward C-p27-time-shuffle-invariance. Flag comment in code explicitly references the carry-forward id. |
| P28 | T | T | Wealth-distribution Gini / cluster index; distributional + time-aggregated. |
| P31 | F | F | |

19 flag entries total (matches the spec table, including the two not-yet-implemented patterns P17 and P19 for future readiness).

## Part D — Harness updates (`epc/phase2a/panel.py`)

- `PANEL_VERSION` bumped `"1.1"` → `"1.2"` in `epc/phase2a/__init__.py`.
- Added imports from `epc.phase2a.detector_invariance` (`get_flags` plus the canonical substrate-id constants `PERM_SHUFFLED_SUBSTRATE`, `TIME_SHUFFLED_SUBSTRATE`, `SKIP_VERDICT`).
- `run_panel` now looks up `invariance = get_flags(pattern_id)` once at start of Class A loop.
- For each Class A substrate:
  - If `sub_id == "permutation_shuffled"` AND `invariance.permutation_invariant` → skip and record `{"verdict": "SKIPPED-degenerate-by-construction", "skip_reason": "primary_metric_permutation_invariant"}`.
  - If `sub_id == "time_shuffled"` AND `invariance.time_shuffle_invariant` → skip with `skip_reason: "primary_metric_time_shuffle_invariant"`.
- TNR computation operates on evaluated substrates only — skipped substrates are filtered out of both the detected-list and the score-list before TNR / Cohen's d.
- JSON schema additions per spec §"Harness output":
  - `panel_version: "1.2"` at top level.
  - New top-level `detector_invariance` block: `{permutation_invariant, time_shuffle_invariant, primary_metric}`.
  - New `summary.class_a_size_total: 10` and `summary.class_a_size_evaluated: N` (N ≤ 10).
  - Skipped substrate entries carry `verdict: "SKIPPED-degenerate-by-construction"` and `skip_reason`.

## Part E — P15 canonical-positive ratification

Added "Canonical positive ratification (v1.2) — Sprint 34" subsection to `REPLICATION_NOTES.md`'s P15 section. Documents:
- The change from R-pentomino to dense-random GoL (`init_density=0.37, L=40, n_steps=300`) per v1.2 spec Change 4.
- Rationale: R-pentomino's small initial activation produces a long high-variance transient where P15's structural-diversity screening metric is noise-dominated. Dense-random IC stabilizes into a high-activity GoL with diverse structures consistently across seeds.
- R-pentomino remains valid for *qualitative* P15 demonstration (gliders, blocks, blinkers) but is not the panel canonical positive.
- No detector logic change.
- No code change needed in Sprint 34 — the panel script's `build_p15_positives` was already set to dense-random in Sprint 33; v1.2 promotes the Sprint 33 workaround to canonical-of-record.

## Part F — Tests (`tests/test_phase2a_panel.py`)

**19 new tests** for Sprint 34 in the phase2a test file (the parametrized flag-match test counts as 10 separate cases, not 1). Coverage:
- Flag values match the spec table for 10 patterns (parametrized — covers P1, P5, P9, P11, P14, P15, P18, P21, P27, P28).
- Unknown-pattern lookup returns the safe default `(False, False)`.
- Lookup-helper booleans (`is_permutation_invariant`, `is_time_shuffle_invariant`) work.
- `PANEL_VERSION == "1.2"`.
- `_CountingStub` end-to-end tests:
  - P9 (both flags True): 2 Class A substrates skipped (permutation_shuffled, time_shuffled) with correct `skip_reason` values.
  - P1 (both flags False): no Class A skips; `class_a_size_evaluated == 10`.
  - Stub pattern with `(True, False)` via `monkeypatch.setitem` on `DETECTOR_INVARIANCE_FLAGS`: only `permutation_shuffled` skipped; `class_a_size_evaluated == 9`.
  - P14 (both True) JSON schema check: `panel_version == "1.2"`, `class_a_size_total == 10`, `class_a_size_evaluated == 8`, `detector_invariance` top-level block present, skipped substrates carry `skip_reason`.
- Hand-constructed TNR-on-evaluated-only: 8 evaluated, 7 rejected → TNR = 0.875 (NOT 7/10).
- P27 explicit provisional flag check (`rationale` contains "Provisional" or the carry-forward id).

**Two stale Sprint 31 tests updated:**
- `test_panel_version_constant_is_one_one` → renamed `test_panel_version_constant_is_current` and updated to check for `"1.2"` (avoids re-naming churn each sprint).
- `test_panel_writes_json_with_v1_1_schema` — assertion bumped from `panel_version == "1.1"` to `panel_version == "1.2"`. Test still exercises the v1.1+ schema fields (`class_b_composition` etc.) which all carry forward.

## Pre-flight + post-flight test counts

- **Pre-flight fast suite:** 566 passed (Sprint 33 ground truth).
- **Post-flight fast suite:** **585 passed, 0 failed, 65 deselected** in 11:47 (+19 net from Sprint 34; the parametrized flag-match test contributes 10 separate test cases per pattern, hence the 19 not the predicted 11).
- **Pre-flight bundle:** 205 (unchanged).
- **Transfer-matrix figures:** 20 / 19 / 79 / 274 / 27 / 19 / 361 / 77 / 284 (unchanged).

## Deviations and judgment calls

### Deviation 1 — Approach 2 (centralized) chosen over Approach 1 (per-module)
The brief said "Approach 1 (preferred if clean)" but explicitly permitted Approach 2 "as long as the harness can read the flags." Approach 1 was not clean for this codebase because (a) detector modules are heterogeneous (class-based, function-based, mixed-internal-structure) and (b) the flags are panel-only metadata that doesn't belong in the detector's own contract. Centralizing kept the change cleanly contained to one new file plus three test-imports. Documented in Part C above.

### Deviation 2 — `test_panel_partial_skip_for_perm_only_invariant` uses a stub pattern, not P5
P5 (Vicsek) is the only spec-table pattern with `(perm_inv=True, time_inv=False)`, but P5's substrate type is `continuous_2d` and its Class B catalog mates are `[P2_abp, P6_dorsogna]` — neither has a generator in `_GENERATORS` yet (declared since Sprint 31 but never landed). Calling `run_panel(..., pattern_id="P5", ...)` triggers a `KeyError` from `load_native_substrate("P2_abp")`. To exercise the partial-skip semantics without depending on P2 / P6 generators landing (which is Sprint 33+ batch-backfill work, out of scope for Sprint 34), the test uses `monkeypatch.setitem(DETECTOR_INVARIANCE_FLAGS, "STUB_PERM_ONLY", InvarianceFlags(True, False, ...))` to add a synthetic stub pattern whose `class_b_for_pattern` returns an empty catalog (substrate_type=None). This exercises the harness skip logic directly without invoking unrelated catalog-mate loading.

The same workaround would not be needed once continuous_2d Class B' supplements / P2 / P6 generators land — at that point the test could switch to a real `pattern_id="P5"` exercise.

### Deviation 3 — `_CountingStub` repeated from Sprint 31 with the same semantics
The Sprint 34 end-to-end tests reuse the discriminating-stub pattern from Sprint 31 (`_DiscriminatingStub`-equivalent). Rather than rename Sprint 31's class, I introduced a parallel `_CountingStub` in the Sprint 34 test block. Same semantics; named differently to localize the Sprint 34 changes and avoid edit-conflicts with the existing test class. Could be consolidated in a future cleanup sprint if desired.

### Note — interpreter substitution
Same as prior sprints: `python3.12 -m pytest`, `PYTHONPATH=. python3.12 ...`. Mac default `python` is unbound.

### Note — no `analysis/run_phase2a_panel.py` change
The brief said the P15 canonical-positive change "may have been set in Sprint 33 already" — it was, the script uses dense-random GoL with `init_density=0.37, L=40, n_steps=300`. v1.2 promotes that to canonical-of-record; no script change required.

### Note — depth_gap.md unchanged
Per brief out-of-scope: "Updating `docs/depth_gap.md` (Sprint 35 will update P9 and P14 rows after re-runs)." Confirmed unchanged this sprint.

## Carry-forward summary

- **C-class-a-oscillator-degenerate / C-class-a-permutation-degenerate** (Sprint 32 / Sprint 33): **CLOSED.** v1.2 invariance flags skip the degenerate substrates for detectors that declare invariance; Sprint 35 P9 + P14 re-runs should both PASS.
- **C-p27-time-shuffle-invariance** (Sprint 34, **NEW, provisional**): P27's `time_shuffle_invariant=True` flag is a spec-table call that must be validated by an actual P27 panel run. If Sprint 35+ P27 results contradict the assumption, change the flag and re-run — same "do not modify the detector to make it pass" rule still applies. Encoded as a comment in `detector_invariance.py` and as an explicit test in `test_phase2a_panel.py`.
- **C-p14-class-c-borderline** (Sprint 33): still OPEN, low priority. One dissipative-sandpile regime at `p_diss=0.350` borderline-fires P14. Resolvable by tightening the p_diss range or bumping n_drive.
- **C-supplements** (Sprint 31): still OPEN. `lattice_1d` / `lattice_2d_continuous` / `scalar_wealth` / `opinion_space` supplements still unimplemented. Needed once Sprint 33+ batch reaches patterns with those substrate types.
- **C-pyproject-pin** (Sprint 29): still OPEN. 1-line `pyproject.toml` change.

## HEAD commit hash and tag at end of sprint

To be recorded after commit + push + tag:

- **Commit:** `__TBD__`
- **Tag:** `v0.34.0` (pushed to origin)
