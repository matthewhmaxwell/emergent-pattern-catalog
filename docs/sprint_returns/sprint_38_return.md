# Sprint 38 Return Summary

**Sprint goal:** Add `_gen_p8_nagel_schreckenberg` native generator + minimum-viable
B' supplement builders for `lattice_1d` substrate type. Final infrastructure piece
before Sprint 44's mixed dim4 batch (P8 + P10). **Status: code complete; commit/push
blocked by transient OS process exhaustion (see Environment Issue section below).**

**Sprint 39 cue:** Dim4 closure chain begins. P22 SIR + P27 Nowak-May (lattice_2d
batch 1).

---

## Pre-flight verification

- Base HEAD: `4e36f37` (Sprint 37 follow-up — sprint_37_return.md added; matches
  expected Sprint 37 post-commit). ✓
- `PYTHONPATH=. python3.12 scripts/count_transfer_matrix.py`: **20 / 19 / 79 / 274 /
  27 / 19 / 361 / 77 / 284** ✓ (unchanged, as expected for infrastructure-only sprint).
- pytest pre-flight: chain-aware skip per brief.

---

## Part A — `_gen_p8_nagel_schreckenberg` added to `epc/phase2a/catalog.py`

**`SUBSTRATE_PARAMS["P8_nagel_schreckenberg"]`** added with canonical jamming-regime
parameters per Nagel-Schreckenberg 1992 and Bette et al. 2017:
- `road_length=100` (panel-scale; L=1000 canonical but qualitative jam pattern present)
- `density=0.3` (deep jam: stopped ≈ 0.43 at L=1000, p_slow=0.3 per model docstring)
- `v_max=5, p_slow=0.3` (NS 1992 illustrative parameter choice)
- `n_steps=200, seed=0`

**`_gen_p8_nagel_schreckenberg(p)`** instantiates `NagelSchreckenberg(L=road_length,
density=..., v_max=..., p_slow=..., seed=...)`, calls `model.run(n_steps=n_steps)`,
and converts each snapshot's `positions` array to a binary cell-occupancy array of
length `road_length`. Returns `{"kind": "sequence", "arrays": (T, road_length)}` in
int8 dtype, matching `_gen_p31_zhang_sorting`'s shape contract.

**Registered** in `_GENERATORS["P8_nagel_schreckenberg"]`.

**`PATTERN_TO_SUBSTRATE_ID["P8"]`** comment updated from
`"declarative; generator NOT yet implemented"` to `"generator added Sprint 38"`.

---

## Part B — `class_b_for_pattern("P31")` now includes `P8_nagel_schreckenberg`

`P8` was already in `PATTERN_TO_SUBSTRATE_ID` and `SUBSTRATE_TYPE_BY_PATTERN["P8"]`
already resolves to `"lattice_1d"` (NS model is registered in MODEL_REGISTRY with
`substrate_type="lattice_1d"`, `primary_patterns=["P8"]`). Adding the generator in
Part A completes the functional path. `class_b_for_pattern("P31")["catalog_mates"]`
now includes `"P8_nagel_schreckenberg"` (and vice versa for P8). ✓

---

## Part C — 2 B' supplement builders for `lattice_1d` in `epc/phase2a/structured.py`

### `independent_lane_traffic`

N cars on a ring of `road_length` cells, each advancing 1 cell per step with no
gap lookahead and no randomization. No jam formation: `stopped_fraction = 0` by
construction. Returns `list[dict]` with `{"array": binary occupancy (int8, length
road_length), "step": t}` per step. Default: `road_length=100`, `n_cars=30`,
`n_steps=200`.

### `reverse_sorted_sequence`

Array of `n` integers initialized to `[n, n-1, ..., 1]` (fully reverse-sorted). Each
step adds i.i.d. uniform noise in `[-noise_scale, +noise_scale]` (rounded) — no swap
or sorting dynamics. Tests that emergent-ordering detectors (P31 delayed gratification)
do not fire on pre-imposed order. Returns `list[dict]` with `{"array": int32 array of
length n, "step": t}` per step. Default: `n=64`, `n_steps=200`, `noise_scale=0.5`.

Both builders:
- Added to `SUPPLEMENTS_BY_SUBSTRATE_TYPE["lattice_1d"]`
- Added to `SUPPLEMENT_BUILDERS` dict

---

## Part D — 4 new tests in `tests/test_phase2a_panel.py`

| Test | What it verifies |
|---|---|
| `test_gen_p8_nagel_schreckenberg_deterministic` | Same seed → byte-identical NS occupancy trajectory |
| `test_gen_p8_nagel_schreckenberg_output_shape` | kind="sequence", ndim=2, shape[1]==road_length |
| `test_class_b_p31_includes_p8_nagel_schreckenberg` | substrate_type="lattice_1d", P8_nagel_schreckenberg in catalog_mates |
| `test_lattice_1d_supplements_registered` | both builders in SUPPLEMENTS_BY_SUBSTRATE_TYPE and SUPPLEMENT_BUILDERS |

All 4 tests verified passing. 150 panel + orchestration tests verified passing (87
test_phase2a_panel.py + 63 test_orchestration.py). Full suite count not confirmed due
to environment issue (see below).

---

## Part E — `docs/paper_CHANGELOG.md` updated

Sprint 38 entry prepended with the 3 required bullet points. Mechanical-only sprint;
no §3–§8 prose changes.

---

## Part F — Carry-forward review

Per brief: C-supplements now closed for both continuous_2d (Sprint 37) and
lattice_1d (Sprint 38).

| Substrate type | C-supplements status |
|---|---|
| `continuous_2d` | **CLOSED** (Sprint 37) |
| `lattice_1d` | **CLOSED** (this sprint) |
| `lattice_2d_continuous` | OPEN |
| `scalar_wealth` | OPEN |
| `opinion_space` | OPEN |

All other open carry-forwards unchanged from Sprint 37:
- **C-class-a-constant-field-trivial-sync (Sprint 35):** Open. P9 residual false
  positive from constant_field substrate. v1.3 candidate; not in scope.
- **C-p27-time-shuffle-invariance (Sprint 34):** Open/provisional. Flag=True per spec.
- **C-p14-class-c-borderline (Sprint 33):** Persists at p_diss=0.350; P14 still PASS.
- **C-pyproject-pin (Sprint 29):** Open. 1-line pyproject.toml change deferred.

---

## Environment issue — OS process exhaustion blocking commit/push

During post-flight verification, three background pytest tasks were launched and left
running simultaneously. This exhausted the OS fork limit on the host machine. All
subsequent `fork()` calls from the shell return `EAGAIN` ("Resource temporarily
unavailable"), making it impossible to run any shell command including `git commit`,
`git tag`, or `git push`.

**Impact:** The 4 modified files are staged (uncommitted) in the working tree.
The commit, tag (`v0.38.0`), and push have NOT been executed.

**Manual recovery steps (operator required):**

```bash
# 1. Verify tests pass (should be 595 fast tests)
cd /home/matthewhmaxwell/emergent-pattern-catalog
python3.12 -m pytest tests/ -m "not slow" -q --no-header

# 2. Commit
git add epc/phase2a/catalog.py epc/phase2a/structured.py \
        tests/test_phase2a_panel.py docs/paper_CHANGELOG.md \
        docs/sprint_returns/sprint_38_return.md
git commit -m "Sprint 38: lattice_1d catalog generator (P8 NS) + B' supplements

- epc/phase2a/catalog.py: add _gen_p8_nagel_schreckenberg; NS jamming-regime
  params (density=0.3, v_max=5, p_slow=0.3, road_length=100, n_steps=200, seed=0);
  register in _GENERATORS; update PATTERN_TO_SUBSTRATE_ID P8 comment.
- epc/phase2a/structured.py: add independent_lane_traffic + reverse_sorted_sequence
  lattice_1d B' supplements; register in SUPPLEMENTS_BY_SUBSTRATE_TYPE and
  SUPPLEMENT_BUILDERS.
- tests/test_phase2a_panel.py: 4 new Sprint 38 tests (deterministic, shape,
  class_b P31 includes P8_NS, lattice_1d supplements registered).
- docs/paper_CHANGELOG.md: Sprint 38 entry appended.

Co-Authored-By: Claude Sonnet 4.6 <noreply@anthropic.com>"

# 3. Tag and push
git tag v0.38.0
git push origin main
git push origin v0.38.0
```

---

## Post-flight verification (partial)

Acceptance criteria verified before fork exhaustion:

| Check | Result |
|---|---|
| `'P8_nagel_schreckenberg' in _GENERATORS` | **PASS** ✓ |
| `class_b_for_pattern('P31')["catalog_mates"]` includes `P8_nagel_schreckenberg` | **PASS** ✓ |
| `SUPPLEMENTS_BY_SUBSTRATE_TYPE["lattice_1d"]` has 2 entries | **PASS** ✓ |
| 4 new Sprint 38 tests pass | **PASS** ✓ (verified via pytest -k "sprint_38 or p8_nagel or lattice_1d") |
| 150 panel + orchestration tests pass | **PASS** ✓ |
| Transfer matrix: **20/19/79/274/27/19/361/77/284** | **PASS** ✓ |
| `pytest tests/ -m "not slow"` → 595 passed | **BLOCKED** (OS fork exhaustion) |
| Commit + tag `v0.38.0` + push | **BLOCKED** (OS fork exhaustion) |

`git diff --stat` (vs Sprint 37 HEAD): 4 files changed (before return doc):
`epc/phase2a/catalog.py`, `epc/phase2a/structured.py`,
`tests/test_phase2a_panel.py`, `docs/paper_CHANGELOG.md`. No other files modified. ✓

---

## Deviations and judgment calls

None. All deliverables implemented exactly as specified. Generator shape contract
matches `_gen_p31_zhang_sorting` (kind="sequence", int8 arrays, (T, road_length)).
Supplement format follows structured.py convention of `list[dict]` with array/step
keys, matching the lattice_1d sequence idiom.

---

## HEAD commit hash and tag at end of sprint

- **Commit:** NOT YET COMMITTED (see Environment Issue above)
- **Tag:** NOT YET TAGGED

**Decision: GO-LIMITED** — All Sprint 38 code is complete and correct; all verifiable
acceptance criteria pass; transfer matrix unchanged. Chain is blocked ONLY by a
transient OS process-exhaustion incident (too many background pytest processes) that
prevents `git commit` / `git tag` / `git push`. Operator must run the manual recovery
steps above, verify 595 fast tests pass, then commit + tag + push. Once committed,
Sprint 39 may proceed autonomously.
