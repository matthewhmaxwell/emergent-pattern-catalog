# Sprint 29 Return Summary

**Sprint goal:** Apply the prepared-but-unpushed Sprint 27 bundle, close Sprint 28 carry-forward C-env (`powerlaw` in `requirements.txt`), close Sprint 27 carry-forward C-gitignore (`*_incr.json` scratchfile rule), and move Sprint 28's root-level `SPRINT_RETURN.md` under `docs/sprint_returns/` to establish the convention. Mechanical cleanup only — no detector / model / test changes. Status: **complete.**

## Pre-flight verification

- Base HEAD: `4250e2c` (Sprint 28, v0.28.0). Working tree clean. ✓
- `python3.12 scripts/count_transfer_matrix.py` (with `PYTHONPATH=.`): **20 / 19 / 79 / 274 / 27 / 19 / 361 / 77 / 284** ✓ (unchanged).
- `python3.12 -m pytest tests/ -m "not slow"`: **508 passed, 0 failed, 65 deselected** in 11:14. ✓ (matches Sprint 28's corrected ground truth.)

## Part A — Sprint 27 bundle landed

Five files copied from the prepared-but-unpushed Sprint 27 bundle:

| Bundle file | Repo path | Op |
|---|---|---|
| `REPLICATION_NOTES.md` | `REPLICATION_NOTES.md` | modified (+145 lines, "## Sprint 27 — close #31" section) |
| `p10_phase_boundary_multiseed.py` | `analysis/p10_phase_boundary_multiseed.py` | new (197 lines) |
| `p10_make_basin_figure.py` | `analysis/p10_make_basin_figure.py` | new (72 lines) |
| `p10_phase_boundary_multiseed.json` | `analysis/outputs/p10_phase_boundary_multiseed.json` | new (1442 lines) |
| `p10_basin_volume_multiseed.png` | `analysis/outputs/p10_basin_volume_multiseed.png` | new (99,663 bytes) |

Verifications:
- `grep -c "## Sprint 27 — close #31" REPLICATION_NOTES.md` → 1 ✓
- All five paths present ✓

## Part B — `powerlaw` in `requirements.txt`

`powerlaw` was already present in `requirements.txt` from the initial commit (`a856c35`) at `>=1.5.0`. The Sprint 28 C-env carry-forward predates verification of the requirements file: the package is listed there but was missing from the Sprint 28 thread's local environment, hence the mid-pre-flight pip install. The dependency line itself was never absent from the file.

Per the brief's pinning instruction (pin to `>=` of installed version), bumped from `powerlaw>=1.5.0` to `powerlaw>=2.0.0` to match the installed version (`pip show powerlaw` → 2.0.0).

`pyproject.toml` still carries the original `powerlaw>=1.5.0` line under `[project] dependencies`. Tightening it was out of scope for this sprint (brief mentioned `requirements.txt` only) — flagging as a 1-line follow-up if version-floor consistency between the two files matters.

Verification: `grep -c "^powerlaw" requirements.txt` → 1 ✓

## Part C — `docs/sprint_returns/` convention established

- `mkdir -p docs/sprint_returns/` ✓
- `git mv SPRINT_RETURN.md docs/sprint_returns/sprint_28_return.md` ✓ (rename preserved, no content change)
- This sprint's return summary lives at `docs/sprint_returns/sprint_29_return.md` (NOT at the repo root).

Verifications:
- `test -f docs/sprint_returns/sprint_28_return.md` ✓
- `! test -f SPRINT_RETURN.md` ✓ (root cleared)
- `test -f docs/sprint_returns/sprint_29_return.md` ✓

## Part D — `.gitignore` cleanup

The existing `.gitignore` (line 33: `outputs/`) already matches `analysis/outputs/*_incr.json` because the bare `outputs/` rule applies to any directory named `outputs/` at any depth. No `.gitignore` modification was required to close C-gitignore — the scratchfile cannot be committed accidentally under the current rules.

A side effect of this rule: the new `p10_phase_boundary_multiseed.json` and `p10_basin_volume_multiseed.png` bundle outputs are also caught by `outputs/` and required `git add -f` to stage. (The pre-existing `analysis/outputs/p10_*.json|.npz|.png` files in this sprint's HEAD were similarly force-added in Sprint 26.) This is consistent with the existing repo convention.

Verification: `git check-ignore -v analysis/outputs/p10_phase_boundary_multiseed_incr.json` → `.gitignore:33:outputs/` ✓

## Post-flight verification

- `python3.12 scripts/count_transfer_matrix.py`: **20 / 19 / 79 / 274 / 27 / 19 / 361 / 77 / 284** ✓ (unchanged).
- `python3.12 -m pytest tests/ -m "not slow"`: **508 passed, 0 failed, 65 deselected** ✓ (matches expected ground truth).

## Acceptance criteria checklist

- [x] `grep -c "## Sprint 27 — close #31" REPLICATION_NOTES.md` → 1
- [x] `analysis/p10_phase_boundary_multiseed.py` exists
- [x] `analysis/p10_make_basin_figure.py` exists
- [x] `analysis/outputs/p10_phase_boundary_multiseed.json` exists
- [x] `analysis/outputs/p10_basin_volume_multiseed.png` exists
- [x] `grep -c "^powerlaw" requirements.txt` → 1
- [x] `docs/sprint_returns/sprint_28_return.md` exists
- [x] `SPRINT_RETURN.md` no longer at repo root
- [x] `docs/sprint_returns/sprint_29_return.md` exists (this file)
- [x] `git check-ignore -v analysis/outputs/p10_phase_boundary_multiseed_incr.json` returns a match (line 33: `outputs/`)
- [x] `python3.12 scripts/count_transfer_matrix.py` → 20 / 19 / 79 / 274 / 27 / 19 / 361 / 77 / 284
- [x] `python3.12 -m pytest tests/ -m "not slow"` → 508 passed
- [x] `git diff --stat`: 1 mod (REPLICATION_NOTES.md), 1 mod (requirements.txt), 1 rename (SPRINT_RETURN.md → docs/sprint_returns/sprint_28_return.md), 5 new files (4 Sprint 27 deliverables + sprint_29_return.md). `.gitignore` not modified — see Part D.

## Deviations and judgment calls

### Deviation 1 — `powerlaw` already in `requirements.txt`
The Sprint 28 C-env carry-forward was logged based on a missing local install, not a missing line in `requirements.txt`. The line `powerlaw>=1.5.0` has been present since initial commit. Acted on the brief's pinning instruction by bumping the floor to `>=2.0.0` (matches the version installed during Sprint 28's mid-pre-flight `pip install`). Pyproject.toml's parallel `powerlaw>=1.5.0` line was left alone (out of scope per brief).

### Deviation 2 — `.gitignore` not modified
Acceptance criterion item is satisfied by the pre-existing `outputs/` rule (line 33), which transitively covers `analysis/outputs/*_incr.json`. Adding a more-specific `analysis/outputs/*_incr.json` line would be redundant. The brief's "if not covered, add" clause anticipated this case. `git diff --stat` therefore shows no `.gitignore` modification.

### Deviation 3 — Force-add required for output files
The existing `outputs/` ignore rule meant `git add` would not pick up the two new files in `analysis/outputs/`. Used `git add -f` (consistent with how Sprint 26's `analysis/outputs/p10_*` files must have been added). Not flagged in the brief; logged here for transparency.

### Note — interpreter
The Mac default `python` is unbound; project tests run under `python3.12` (`/usr/local/bin/python3`, pytest 9.0.2, powerlaw 2.0.0). Where the brief reads `python` / `pytest`, I substituted `python3.12 -m pytest` and `PYTHONPATH=. python3.12 scripts/count_transfer_matrix.py`. No effect on test counts or matrix figures.

## Carry-forward status

- **C-env:** CLOSED. `powerlaw` floor pinned in `requirements.txt`.
- **C-gitignore:** CLOSED. Scratchfile already covered by existing `outputs/` rule; verified with `git check-ignore`.
- **Optional follow-up (1 line):** Bump `powerlaw` floor in `pyproject.toml` from `>=1.5.0` to `>=2.0.0` for consistency with `requirements.txt`. Not load-bearing.

## HEAD commit hash and tag at end of sprint

To be recorded after commit + push + tag:

- **Commit:** `__TBD__`
- **Tag:** `v0.29.0` (pushed to origin)
