# Sprint 31 Return Summary

**Sprint goal:** Land Phase-2a panel spec v1.1 (substrate-typed Class B + Class C N/A escape hatch), archive v1.0 spec and prototype JSONs, and update the harness Class B logic. **No new panel runs; no detector changes.** Status: **complete.** Closes Sprint 30 carry-forward C-panel-spec.

## Pre-flight verification

- Base HEAD: `43ef2b1` (Sprint 30, v0.30.0). Working tree clean. ✓
- `python3.12 scripts/count_transfer_matrix.py`: **20 / 19 / 79 / 274 / 27 / 19 / 361 / 77 / 284** ✓ (unchanged).
- Pre-flight fast tests: 543 passed (Sprint 30 ground truth, not re-run pre-flight to save time — re-verified post-flight).

## File inventory

**Part A — v1.0 archived (3 git mv operations):**
- `docs/phase2a_panel_spec.md` → `docs/archive/phase2a_panel_spec_v1_0.md`
- `analysis/outputs/p18_phase2a_panel.json` → `analysis/outputs/archive/v1_0/p18_phase2a_panel.json`
- `analysis/outputs/p9_phase2a_panel.json` → `analysis/outputs/archive/v1_0/p9_phase2a_panel.json`

**Part B — v1.1 spec landed:**
- `docs/phase2a_panel_spec.md` (byte-identical to bundle's `phase2a_panel_spec_v1_1.md`).

**Part C — Harness updates:**
- `epc/phase2a/__init__.py` — `PANEL_VERSION` bumped `"1.0"` → `"1.1"`.
- `epc/phase2a/structured.py` — **NEW** module with 4 Class B' supplement builders:
  - `incoherent_phases`, `subcritical_kuramoto` (oscillator)
  - `random_graph_evolution`, `network_random_walks` (network)
  - `SUPPLEMENTS_BY_SUBSTRATE_TYPE` registry maps `oscillator` → 2 builders, `network` → 2 builders.
- `epc/phase2a/catalog.py`:
  - `PATTERN_TO_SUBSTRATE_ID` — declares all 19 implemented patterns. Six (P1, P3, P5, P9, P10, P12, P14, P15, P18, P27, P31) have native generators landed in Sprint 30; the rest (P2, P6, P8, P11, P13, P21, P22, P28) are declarative-only — Class-B-membership computation works for them but loading the substrate would raise (Sprint 32+ task).
  - `_build_substrate_type_by_pattern()` — derives the pattern → substrate_type map from `epc.orchestration.MODEL_REGISTRY`, with explicit v1.1 spec overrides for `P18` (`lattice_2d` → `network`) and `P21` (`opinion_space` → `network`).
  - `SUBSTRATE_TYPE_BY_PATTERN` — module-level constant.
  - `class_b_for_pattern(pattern_id)` — returns `{"catalog_mates": [...], "synthetic_supplements": [...], "substrate_type": ...}`. Substrate-typed selection from PATTERN_TO_SUBSTRATE_ID, capped at 10 mates. Adds Class B' supplements when mates < 3.
  - The legacy v1.0 `catalog_ids_for_pattern` is retained but documented as superseded.
- `epc/phase2a/panel.py`:
  - `run_panel` rewritten to use `class_b_for_pattern` for Class B + B'.
  - Class C N/A handling: reads `failed_regime_module.CONFIG.get("status")`; when `"N/A"`, Class C is skipped, the pooled-negative pool excludes it, and `class_c_status` + `class_c_n_a_reason` are written to the JSON.
  - JSON schema additions per spec §"Harness output (v1.1)": `class_b_composition` (`substrate_type`, `catalog_mates`, `synthetic_supplements`), `class_c_status`, `class_c_n_a_reason`. Per-class TNRs report a sub-object `{n, tnr, advisory}` (with the flat `..._tnr` fields preserved for backward compatibility).
  - `overall_verdict` rewritten for v1.1: PASS / PASS-with-weakness / PARTIAL / FAIL with explicit Cohen's d gating (≥1.0 for PASS, ≥0.5 for PARTIAL, <0.5 → FAIL).
  - Per-class TNR with ≥5-substrate gating threshold (`ADVISORY_CLASS_SIZE = 5`); classes below that report TNR but do not gate the verdict.
- `epc/phase2a/failed_regimes/p18_voter.py` — replaced with a Class C N/A declaration. `CONFIG = {"status": "N/A", "n_a_reason": ..., "regimes": []}`. `build_substrate` now raises `RuntimeError` defensively (it must not be reached when status is N/A).
- `epc/phase2a/failed_regimes/p9_kuramoto.py` — unchanged (Class C populated for P9; sub-K_c regimes were already clean in Sprint 30).

**Part D — Tests (`tests/test_phase2a_panel.py`):**
- Net **46 tests, all passing in 0.17 s** (vs 35 in Sprint 30). Test count delta: **+11 net**.
- New tests: `class_b_for_pattern` for P9 / P18 / P21 / P1 / unknown; substrate-type override sanity; supplement determinism; supplement format; v1.1 verdict labels (PASS / PASS-with-weakness / PARTIAL / FAIL plus the TNR-passes-but-d-fails edge case); Class C N/A skipping; JSON schema additions; v1.1 always-fire-stub FAIL semantics.
- Updated tests: `test_panel_version_constant_is_one_one` (was `_one_zero`); P18 failed-regime check now asserts `status == "N/A"` rather than counting 10 regimes; v1.0 verdict tests rewritten as v1.1 verdict tests with explicit `cohens_d_value` argument; stub end-to-end tests use a discriminating-stub class so the PASS path can be exercised under Cohen's-d-gated verdict logic.

**Part E — This file:**
- `docs/sprint_returns/sprint_31_return.md`

## Substrate-type taxonomy used

| pattern_id | substrate_type (registry) | substrate_type (effective in v1.1) | Class B mates | Class B' supplements |
|---|---|---|---|---|
| P1  | lattice_2d | lattice_2d | 7 (P11, P12, P13, P14, P15, P22, P27) | 0 |
| P2  | continuous_2d | continuous_2d | 2 (P5, P6) | 0 (no continuous_2d supplements defined yet) |
| P3  | lattice_2d_continuous | lattice_2d_continuous | 0 | 0 (no lattice_2d_continuous supplements defined yet) |
| P5  | continuous_2d | continuous_2d | 2 (P2, P6) | 0 |
| P6  | continuous_2d | continuous_2d | 2 (P2, P5) | 0 |
| P8  | lattice_1d | lattice_1d | 1 (P31) | 0 (no lattice_1d supplements defined yet) |
| P9  | oscillator | oscillator | 1 (P10) | **2** (incoherent_phases, subcritical_kuramoto) |
| P10 | oscillator | oscillator | 1 (P9) | **2** |
| P11 | lattice_2d | lattice_2d | 7 | 0 |
| P12 | lattice_2d | lattice_2d | 7 | 0 |
| P13 | lattice_2d | lattice_2d | 7 | 0 |
| P14 | lattice_2d | lattice_2d | 7 | 0 |
| P15 | lattice_2d | lattice_2d | 7 | 0 |
| P18 | lattice_2d | **network** ← v1.1 override | 1 (P21) | **2** (random_graph_evolution, network_random_walks) |
| P21 | opinion_space | **network** ← v1.1 override | 1 (P18) | **2** |
| P22 | lattice_2d | lattice_2d | 7 | 0 |
| P27 | lattice_2d | lattice_2d | 7 | 0 |
| P28 | scalar_wealth | scalar_wealth | 0 | 0 (no scalar_wealth supplements defined yet) |
| P31 | lattice_1d | lattice_1d | 1 (P8) | 0 |

## Deviations and judgment calls

### Deviation 1 — substrate-type overrides only for P18 and P21
The v1.1 spec table introduces a **`network`** substrate type that does not exist in the orchestrator registry (`epc.orchestration.MODEL_REGISTRY` uses `lattice_2d`, `lattice_2d_continuous`, `continuous_2d`, `lattice_1d`, `oscillator`, `opinion_space`, `scalar_wealth`). The brief tells us to "build the constant by inspecting the registry. ... If the substrate-type assignment for any pattern differs from v1.1 spec's illustrative table, document the chosen assignment and rationale." Per the brief test prescriptions (P18 → network with [P21] mate + 2 supps; P21 → network with [P18] mate + 2 supps), I applied the `network` override for **only P18 and P21**, leaving every other pattern's substrate type at its registry value. The override is implemented as an explicit dict mutation in `_build_substrate_type_by_pattern()` so it is greppable and easy to revisit in Sprint 32+.

The spec table also lists P3 as `lattice_2d`, P11 as `continuous_2d`, P28 as `lattice_2d` and adds non-registered patterns (P17 Berdahl, P19 Couzin) — those are illustrative, not authoritative. Where the registry disagrees with the spec table for a pattern other than P18/P21, the registry wins.

### Deviation 2 — only oscillator and network B' supplements implemented this sprint
The spec table flags supplements needed for: oscillator (P9, P10), network (P18, P21), lattice_1d (P31, plus P8), lattice_2d_continuous (P3), scalar_wealth (P28), opinion_space (only P21 if overridden away from network). The brief explicitly named only oscillator and network supplements. I implemented those four supplement builders. The other substrate types currently return `synthetic_supplements: []` from `class_b_for_pattern`. Adding lattice_1d / lattice_2d_continuous / scalar_wealth / opinion_space supplements is **deferred** — none of the three Sprint 32 prototype re-run targets (P9, P18) needs them, and Sprint 33+ batch runs against P31, P3, P28 will need them defined before those panels can run cleanly. **Logged as carry-forward C-supplements.**

### Deviation 3 — declarative-only catalog entries (no native generators yet)
`PATTERN_TO_SUBSTRATE_ID` declares all 19 registry patterns so that `class_b_for_pattern` can compute substrate-type mates correctly. Eight of those patterns (P2, P6, P8, P11, P13, P21, P22, P28) do **not** have native catalog generators in `_GENERATORS` / `SUBSTRATE_PARAMS`. They appear correctly as catalog-mates for other patterns, but if a v1.1 panel run actually tries to *load* one of them (e.g., P1's panel attempting to load P22_sir_epidemic), the call raises `KeyError`. Sprint 32 will need to add at least the generators required by P9 and P18 panel re-runs (the brief notes for Sprint 32 only mention P9 and P18, both of whose Class B mates are already implemented: P10_chimera for P9, and P21_hegselmann_krause for P18 — the latter still needs a generator). **Concretely, Sprint 32 must implement `_gen_p21_hegselmann_krause` before re-running the P18 panel.** Logged as carry-forward C-p21-generator.

### Deviation 4 — `np.var` floating-point noise required a `cohens_d` clamp
While running the rewritten end-to-end stub tests I discovered that `np.var([0.9]*20, ddof=1)` returns ~5e-32 instead of exactly 0 due to ULP-level accumulation in the variance algorithm. Without a degeneracy clamp, the always-fire stub (positives and negatives all scoring 0.9) computed `cohens_d` ≈ 1.0 from `(0.9 - 0.9000000000000002) / sqrt(5e-32)` — a meaningless artifact that masked the true "no signal" case. Fixed by treating `pooled_var ≤ 1e-12` as degenerate in `cohens_d`, returning ±inf or 0 from the mean difference rather than dividing by tiny noise. Also tightened the inf handling in `overall_verdict` so that legitimately-degenerate-perfect detectors (positives all fire, negatives all reject → ±inf d) still gate correctly.

### Note — interpreter substitution
Same as Sprint 29/30: `python3.12 -m pytest`, `PYTHONPATH=. python3.12 ...`. Mac default `python` is unbound.

### Note — analysis/run_phase2a_panel.py left at v1.0 wiring
The one-shot run script was not updated this sprint (the brief says no new panel runs). It still uses the v1.0 `n_steps`/`record_every` defaults appropriate for Kuramoto positives. Sprint 32 will revise the script as part of the v1.1 re-run effort.

## Pre-flight + post-flight test counts

- **Pre-flight fast suite:** 543 passed (Sprint 30 ground truth).
- **Post-flight fast suite:** **554 passed, 0 failed, 65 deselected** in 11:06 (matches predicted 543 + 11).
- **Pre-flight bundle (transfer matrix counts + cross-detection matrix + voter P18 e2e + orchestration):** 205 (unchanged).
- **Transfer-matrix figures:** 20 / 19 / 79 / 274 / 27 / 19 / 361 / 77 / 284 (unchanged).

## Carry-forward summary

- **C-panel-spec (Sprint 30):** **CLOSED.** v1.1 spec landed; harness updated; v1.0 spec + prototype JSONs archived under `docs/archive/` and `analysis/outputs/archive/v1_0/`.
- **C-supplements (Sprint 31, NEW):** Implement Class B' supplement builders for `lattice_1d`, `lattice_2d_continuous`, `scalar_wealth`, and (if not overriding) `opinion_space` substrate types. Required before Sprint 33+ batch runs that include P31 / P3 / P28. Brief constrained Sprint 31 to oscillator + network supplements.
- **C-p21-generator (Sprint 31, NEW):** Implement `_gen_p21_hegselmann_krause` in `epc/phase2a/catalog.py`. Required before Sprint 32 P18 panel re-run, since P18's Class B catalog-mate under the network override is `P21_hegselmann_krause`. The HK model itself exists at `epc/models/hegselmann_krause.py`; the generator is a small new function that calls it.
- **C-pyproject-pin (Sprint 29):** Bump `powerlaw>=1.5.0` to `>=2.0.0` in `pyproject.toml`. Still open. 1-line change.

## HEAD commit hash and tag at end of sprint

To be recorded after commit + push + tag:

- **Commit:** `__TBD__`
- **Tag:** `v0.31.0` (pushed to origin)
