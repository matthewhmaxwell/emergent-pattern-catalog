# Phase-2a Standard Negative Panel — Specification v1.1

This is a revision of v1.0 (`docs/phase2a_panel_spec.md`, dated Sprint 30). v1.1 supersedes v1.0; the v1.0 file should be archived to `docs/archive/phase2a_panel_spec_v1_0.md` rather than deleted.

The Sprint 30 prototype runs (P18 voter, P9 Kuramoto) both returned PARTIAL under v1.0. Investigation revealed two structural issues with the v1.0 spec rather than detector quality issues. v1.1 fixes both.

---

## What changed and why

### Change 1 — Class B is now substrate-typed

**v1.0 problem:** Class B specified "10 catalog-derived non-positives... run for the same duration and grid size as the detector under test (where compatible)." The "where compatible" was doing too much work. Vicsek headings adapted to Kuramoto phases retain alignment; a binary GoL grid mapped to phases via `value × π` gives near-uniform phases. In each case the adapted substrate genuinely *is* the pattern, so the detector firing isn't a false positive — it's correctly detecting structure that the cross-format adaptation introduced.

**v1.1 fix:** Each pattern's Class B contains only catalog-derived substrates that share its native substrate type. The catalog already maintains four substrate types in the orchestrator's substrate-aware dispatch system (`lattice_1d`, `lattice_2d`, `continuous_2d`, `oscillator`); v1.1 inherits that taxonomy. A pattern's Class B substrates are drawn from the canonical positives of *other* patterns whose models target the same substrate type.

This shrinks Class B per pattern (from a fixed 10 to a variable 0–N depending on substrate type membership). For patterns whose substrate type has fewer than 3 catalog-mates, Class B is supplemented with synthetic structured substrates (see Class B' below). For patterns whose substrate type has no catalog-mates at all, Class B is reported as N/A and the panel runs on Classes A and C only.

### Change 2 — Class C has an N/A escape hatch

**v1.0 problem:** Class C ("the detector's own model run in 10 parameter regimes that do NOT produce the target pattern") assumes every pattern has a parameter regime that suppresses it. Some patterns don't. The voter model has no parameter that suppresses consensus; running it long enough always reaches consensus from any initial condition that isn't trivially balanced. The Sprint 30 P18 prototype used extreme initial conditions (0.93–0.999) as a proxy, which trips the detector ~60% of the time — but those substrates start in or near consensus, so the detector firing is a true positive in disguise, not a failure of specificity.

**v1.1 fix:** Class C is N/A for patterns that lack a parameter regime suppressing the pattern. The catalog enumerates these patterns up front (see "Class C N/A list" below). For these patterns the panel runs on Classes A and B only. The PASS criterion adjusts proportionally — see "PASS criterion" below.

### Change 3 — PASS criterion adjusts for class N/A and per-class minimum size

**v1.0 problem:** The 30-substrate panel was a fixed assumption. With Class B varying by substrate type and Class C sometimes N/A, panel sizes range from ~15 (one class N/A) to 30+ (full). The PASS criterion in v1.0 was "≥95% TNR overall." This stays. But per-class size now matters for whether per-class TNR is meaningful.

**v1.1 fix:** PASS criterion stays at ≥95% TNR overall and per-class TNR is reported only when that class has ≥5 substrates. Classes with <5 substrates are reported as "N≤4 — TNR: x/N (advisory only)" and do not gate PASS.

---

## Composition (v1.1)

The panel for pattern P_i contains substrates from up to three classes. The composition for P_i is determined by P_i's substrate type and Class C eligibility.

### Class A — Synthetic null substrates (10, fixed across all patterns)

Unchanged from v1.0. The 10 synthetic substrates listed in v1.0 §Class A apply identically.

### Class B — Substrate-typed catalog-derived non-positives (variable, 0–N)

**Substrate-type taxonomy (inherited from ADR #25):**

| Substrate type | Patterns whose canonical positive uses this substrate type |
|---|---|
| `lattice_1d` | P31 (Zhang sorting) |
| `lattice_2d` | P1 (Schelling), P3 (Gray-Scott), P12 (RPS), P13 (GH excitable), P14 (BTW sandpile), P15 (GoL), P22 (SIR-on-grid), P27 (Nowak-May), P28 (Yard-Sale-on-grid), and any other lattice-based positives |
| `continuous_2d` | P5 (Vicsek), P6 (D'Orsogna mill), P11 (LV continuous), P17 (Berdahl), P19 (Couzin), and any other particle/continuous positives |
| `oscillator` | P9 (Kuramoto), P10 (chimera) |
| `network` | P18 (voter on graph), P21 (Hegselmann-Krause) — opinion dynamics on graphs |

If the audit reveals additional substrate types in implemented patterns, this table extends. Default rule: if a pattern's substrate type is ambiguous or hybrid, assign it to its primary substrate type and note the hybrid in the panel run output.

**Class B membership rule:** for pattern P_i with substrate type T, Class B contains the canonical positives of all *other* implemented patterns with substrate type T, up to a maximum of 10. If fewer than 3 such patterns exist, Class B is supplemented from Class B' below.

**Class B' — Synthetic structured non-pattern substrates (substrate-typed):** when Class B has fewer than 3 catalog-mates, supplement with substrate-typed structured-but-non-pattern substrates. For example: a `continuous_2d` pattern with no catalog-mates gets a Class B' of randomly-walking particles (no alignment), uniform random positions with random walks, etc. — substrates that are "live" in the right substrate but lack the target pattern. Document the chosen Class B' substrates in the panel run JSON.

**Class B membership table (computed from current catalog):**

| Pattern | Substrate type | Class B catalog-mates (count) | Class B' supplement needed? |
|---|---|---|---|
| P1 (Schelling) | lattice_2d | 8 | no |
| P3 (Gray-Scott) | lattice_2d | 8 | no |
| P5 (Vicsek) | continuous_2d | 4 | no |
| P6 (D'Orsogna mill) | continuous_2d | 4 | no |
| P9 (Kuramoto) | oscillator | 1 (P10) | yes (need ≥2 supplements) |
| P10 (chimera) | oscillator | 1 (P9) | yes (need ≥2 supplements) |
| P11 (LV) | continuous_2d | 4 | no |
| P12 (RPS) | lattice_2d | 8 | no |
| P13 (GH) | lattice_2d | 8 | no |
| P14 (BTW) | lattice_2d | 8 | no |
| P15 (GoL) | lattice_2d | 8 | no |
| P17 (Berdahl) | continuous_2d | 4 | no |
| P18 (voter) | network | 1 (P21) | yes (need ≥2 supplements) |
| P19 (Couzin) | continuous_2d | 4 | no |
| P21 (HK) | network | 1 (P18) | yes (need ≥2 supplements) |
| P22 (SIR) | lattice_2d | 8 | no |
| P27 (Nowak-May) | lattice_2d | 8 | no |
| P28 (Yard-Sale) | lattice_2d | 8 | no |
| P31 (Zhang) | lattice_1d | 0 | yes (need 3 — full Class B') |

The membership counts assume all 19 implemented patterns are in their listed substrate types. Code should compute this table programmatically from the registry rather than copying from this spec — the spec table is illustrative.

### Class C — Failed-regime substrates (10, pattern-specific) OR N/A

Unchanged from v1.0 in form: 10 failed-regime substrates per pattern, specified by the pattern-specific config in `epc/phase2a/failed_regimes/<pattern_id>.py`.

**Class C N/A list:** patterns for which Class C is N/A because no parameter regime suppresses the pattern:

| Pattern | Reason for N/A |
|---|---|
| P15 (GoL) | Deterministic; canonical positive (R-pentomino) is a fixed initial condition. No parameter to vary. |
| P18 (voter) | No parameter regime suppresses consensus given non-trivial initial conditions. |
| P31 (Zhang sorting) | Algorithm always sorts; no parameter regime that prevents convergence to sorted state. |

This list may extend as further patterns are characterized. A pattern's Class C is N/A when both of these hold: (a) the model has no parameter whose value gates whether the pattern emerges, and (b) all reasonable initial conditions reach the pattern. Default to required if uncertain — the burden is on showing that no parameter regime exists.

For patterns where Class C is N/A, the panel runs on Classes A and B (or A only, if both B and C are N/A — which would itself be a finding worth flagging).

---

## PASS criterion (v1.1)

A detector passes Dim 4 against the v1.1 panel when:

1. **Overall TNR ≥ 95%** across all substrates in the panel for that pattern. (At most 1 false positive per 20 substrates.)
2. **Per-class TNR is reported** for every class with ≥5 substrates. Classes with <5 are reported as "TNR: x/N (advisory)".
3. **Per-class TNR ≥ 90% is the soft expectation** for classes with ≥5 substrates. PASS-with-weakness if any class drops below 90%.
4. **Effect size on the primary metric: Cohen's d ≥ 1.0** between the canonical positive and the pooled non-N/A panel.

Verdicts:
- **PASS** — meets (1)–(4) cleanly.
- **PASS-with-weakness** — meets (1) and (4) but a class with ≥5 substrates is below 90%.
- **PARTIAL** — fails (1) but Cohen's d ≥ 0.5 (detector has signal but specificity is below the bar).
- **FAIL** — fails (1) and Cohen's d < 0.5 (detector lacks signal against this panel).

---

## Harness output (v1.1)

The JSON schema gains two fields:

- `panel_version`: `"1.1"` (was `"1.0"`).
- `class_b_composition`: `{"catalog_mates": [...], "synthetic_supplements": [...]}` — explicit record of which substrates were drawn from the catalog vs. synthesized as Class B' supplements.
- `class_c_status`: `"populated"` or `"N/A"` with `n_a_reason` if N/A.

Per-class TNR fields handle the variable-size case: `synthetic_tnr`, `catalog_tnr` (or `catalog_tnr_advisory` if N<5), `failed_regime_tnr` (or `null` if N/A).

The `verdict` field uses the v1.1 verdict labels above.

---

## Migration path

1. Sprint 32 (code-led): Apply v1.1 spec. Re-run the panel against P18 and P9 under v1.1. Both should PASS or PASS-with-weakness; if either does not, the spec is wrong, not the detectors. Same "do not modify the detector to make it pass" rule from Sprint 30.
2. After P18 and P9 both pass v1.1 cleanly, Sprint 33+ runs v1.1 against the remaining 13 PARTIAL detectors in batches.
3. v1.0 results in `analysis/outputs/p<i>_phase2a_panel.json` from Sprint 30 are archived to `analysis/outputs/archive/v1_0/` rather than overwritten. The v1.1 results overwrite the active panel files.
