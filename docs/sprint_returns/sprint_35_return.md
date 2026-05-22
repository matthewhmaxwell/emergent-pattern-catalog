# Sprint 35 Return Summary

**Sprint goal:** Re-run the Phase-2a panel against P9 (Kuramoto) and P14 (BTW sandpile) under v1.2 spec, validate that both PASS now that degenerate-by-construction substrates are skipped, and update `docs/depth_gap.md` rows for those two patterns only. **Status: complete.**

**Sprint 36 decision: GO.** Both patterns PASS (P9 PASS-with-weakness, P14 PASS clean). The v1.2 invariance-flag fix lands cleanly. Sprint 36 can start the lattice_2d batch.

## Pre-flight verification

- Base HEAD: `b0affd0` (Sprint 34, v0.34.0). Working tree clean. ✓
- `python3.12 scripts/count_transfer_matrix.py`: **20 / 19 / 79 / 274 / 27 / 19 / 361 / 77 / 284** ✓ (unchanged).
- Pre-flight fast tests: 585 passed (Sprint 34 ground truth).

## Part A — P9 v1.2 panel result

**Verdict: PASS-with-weakness.** Overall TNR = **0.952**, Cohen's d = **4.781**.

| Class | TNR | n_eval / n_total | Notes |
|---|---|---|---|
| Synthetic (Class A) | 0.875 | 7 / 8 | weak (< 0.90); only false positive: `constant_field` |
| Catalog (substrate-typed: oscillator) | 1.000 | 3 / 3 | advisory (n<5) |
| Failed-regime sub-K_c | **1.000** | 10 / 10 | gating |
| **Overall** | **0.952** | 20 / 21 | crosses 0.95 PASS gate |

- **Class A skip evidence:** `permutation_shuffled` and `time_shuffled` both appear in the JSON with `verdict="SKIPPED-degenerate-by-construction"`; `class_a_size_total=10`, `class_a_size_evaluated=8`. `detector_invariance` block: `{permutation_invariant: True, time_shuffle_invariant: True, primary_metric: "Kuramoto order parameter r"}`.
- **`constant_field` still trips** (as the brief predicted). All phases ≡ 0 → r = 1 trivially → P9 (correctly given its primary metric) fires. This is the same family of issue as the permutation degenerates but is **not** a within-substrate transformation of the positive, so the v1.2 invariance flags don't catch it. **New carry-forward C-class-a-constant-field-trivial-sync** logged; do not modify spec or detector this sprint per the brief.
- v1.1 → v1.2 delta: 0.913 → 0.952 (crosses PASS gate). Output: `analysis/outputs/p9_phase2a_panel.json`.

## Part B — P14 v1.2 panel result

**Verdict: PASS.** Overall TNR = **0.960**, Cohen's d = **10.585**.

| Class | TNR | n_eval / n_total | Notes |
|---|---|---|---|
| Synthetic (Class A) | **1.000** | 8 / 8 | gating; clean |
| Catalog (substrate-typed: lattice_2d) | 1.000 | 7 / 7 | gating |
| Failed-regime (dissipative sandpile) | 0.900 | 9 / 10 | gating; 1 borderline at `p_diss=0.350` (C-p14-class-c-borderline persists, low priority) |
| **Overall** | **0.960** | 24 / 25 | crosses 0.95 PASS gate |

- Class A skip evidence: same two substrates SKIPPED, evaluated size = 8, **all 8 reject cleanly**.
- Class C borderline at `p_diss=0.350` persists from Sprint 33 (mid-range dissipation retains heavy-tailed structure). Failed-regime TNR = 0.900 is **exactly at** the weak-class threshold (≥0.90 = not weak) → not flagged as weakness.
- Cohen's d = 10.585 (positives at CONFIRMATION/0.700, pooled negatives mostly 0.0 with one borderline at 0.350).
- v1.1 → v1.2 delta: 0.889 → 0.960 (crosses PASS gate). Output: `analysis/outputs/p14_phase2a_panel.json`.

## Part C — `docs/depth_gap.md` row updates

| Pattern | dim4 | Grade change | Notes |
|---|---|---|---|
| P9  | PARTIAL → **PASS** | GAP → **AT-DEPTH** | dim1/dim2/dim3 already PASS; v1.2 panel PASS-with-weakness closes dim4 |
| P14 | PARTIAL → **PASS** | GAP (narrowed, still GAP) | dim2 still PARTIAL (single-run τ, no bootstrap); v1.2 panel closes dim4 only |

**Header counts updated:** AT-DEPTH **4 → 5** (P9 joins the AT-DEPTH set with P15, P18, P28, P31). GAP **15 → 14**. The aggregate-findings summary in the lower section of `depth_gap.md` was also updated for consistency.

No other rows modified.

## Part D — Sprint 36 GO/GO-LIMITED/NO-GO decision

**Decision: GO.**

Both target patterns PASS under v1.2 — P9 PASS-with-weakness, P14 PASS clean. The v1.2 invariance-flag mechanism lands as designed and closes C-class-a-permutation-degenerate cleanly (no regressions in other classes; Class B and C results are positively clean across both patterns).

Per the brief's GO taxonomy:

> "Both PASS / PASS-with-weakness: **GO**. Sprint 36 starts batch."

**Suggested initial batch for Sprint 36 (per the brief's note):** P22 (SIR), P27 (Nowak-May), P28 (Yard-Sale). Considerations:

- **P22 SIR** — generator landed in Sprint 33 (lattice_2d, `_gen_p22_sir_epidemic`). Brief mentioned this. Class B = 7 lattice_2d mates (all generators landed). Expected to follow P14's pattern (clean cat + Class A + Class C; likely PASS).
- **P27 Nowak-May** — generator already in catalog (Sprint 30). `time_shuffle_invariant=True` is **provisional** per C-p27-time-shuffle-invariance (Sprint 34 carry-forward); Sprint 36's P27 panel run is the natural validation. Class B = 7 lattice_2d mates. Expected PASS if the provisional flag holds; if not, Class A 7/8 evaluated with `time_shuffled` (incorrectly) reinstated could PARTIAL on a real FP — that finding would close the C-p27-time-shuffle-invariance carry-forward by flipping the flag to False.
- **P28 Yard-Sale** — already AT-DEPTH per audit; would serve as a sanity check (like P15 in Sprint 33). substrate_type is `scalar_wealth` per the registry, which doesn't have B' supplements yet (C-supplements) and lacks lattice_2d catalog mates, so this would be the first scalar_wealth panel run and may surface new infrastructure gaps. **Recommend deferring P28 to a later sprint after Sprint 36 surfaces what scalar_wealth needs.**

**Recommended Sprint 36 batch (revised):** **P22 SIR + P27 Nowak-May** (both lattice_2d, both PARTIAL on dim4, full Class B coverage already landed). P28 deferred. If both PASS, Sprint 37 takes another lattice_2d pair (P3 GS, P12 RPS, P13 GH — choose two).

**Carry-forwards to watch in Sprint 36:**
- C-p27-time-shuffle-invariance — Sprint 36 P27 panel will validate or invalidate.
- C-class-a-constant-field-trivial-sync — will affect any future detector whose primary metric is uniformity-sensitive (Kuramoto-like). Not in scope for Sprint 36; flagged for the lattice_2d patterns to confirm `constant_field` rejects there (it should, since grid detectors look for spatial structure or dynamics).

## Pre-flight + post-flight test counts

- **Pre-flight fast suite:** 585 passed (Sprint 34 ground truth).
- **Post-flight fast suite:** **585 passed, 0 failed, 65 deselected** in 11:13 (unchanged from Sprint 34, exactly as expected; no new tests this sprint).
- **Pre-flight bundle:** 205 (unchanged).
- **Transfer-matrix figures:** 20 / 19 / 79 / 274 / 27 / 19 / 361 / 77 / 284 (unchanged).

## Deviations and judgment calls

### Deviation 1 — P9 came in slightly better than predicted (PASS-with-weakness, not PARTIAL)
The brief predicted P9 was "likely to PARTIAL again under v1.2" because `constant_field` was expected to trip and give 8/9 syn = 0.889 overall. Actual: 7/8 evaluated syn = 0.875 (still weak), but the smaller denominator + perfect cat + perfect fai produced overall 20/21 = **0.952**, just above the PASS gate. The `constant_field` failure mode is the same one the brief anticipated and is logged as a new carry-forward; the verdict label is just slightly more favourable than predicted because of the math on the smaller evaluated pool.

### Deviation 2 — P14 fai = 0.900 is **at** the threshold, not below
The v1.2 weak-class threshold is `< 0.90`. P14's fai = 0.900 sits exactly at the boundary → not flagged as weak → verdict is plain PASS, not PASS-with-weakness. The C-p14-class-c-borderline carry-forward from Sprint 33 (`p_diss=0.350` cell) is unchanged.

### Deviation 3 — Sprint 36 batch composition (P28 deferred)
The brief's suggested batch was P22 + P27 + P28. After re-reading the depth-gap context (P28 substrate is `scalar_wealth`, with no Class B' supplements yet per C-supplements and no scalar_wealth catalog mates), I recommend P28 be deferred until at least one scalar_wealth supplement builder lands. The brief explicitly tagged P28 as already AT-DEPTH, so this is purely about avoiding infrastructure thrash. **Revised batch: P22 + P27** (two PARTIAL lattice_2d patterns, full Class B coverage).

### Deviation 4 — depth_gap.md P14 row marked **GAP-narrowed** rather than AT-DEPTH
The brief said "if PASS/PASS-with-weakness change dim4 to PASS … recompute grade (AT-DEPTH if all dimensions PASS)". P14's dim2 was already PARTIAL per the Sprint 28 audit (`τ=1.247 from a single 100k-event run; ≥5-seed bootstrap dispersion not reported`). Closing dim4 alone is not sufficient to flip P14 to AT-DEPTH — that needs C4-style multi-seed bootstrap work (still in open carry-forwards from the Sprint 28 audit). P14's row notes are updated to reflect that dim4 closes but dim2 remains the open gap. Header AT-DEPTH count therefore goes 4 → 5 from P9 only, not 4 → 6.

### Note — interpreter substitution
Same as prior sprints: `python3.12 -m pytest`, `PYTHONPATH=. python3.12 ...`. Mac default `python` is unbound.

## Carry-forward summary

- **C-class-a-permutation-degenerate (Sprint 33):** Validated CLOSED in this sprint. P9 and P14 both crossed the PASS gate under v1.2; the SKIPPED-degenerate-by-construction mechanism works.
- **C-class-a-constant-field-trivial-sync (Sprint 35, NEW):** P9's residual false positive is `constant_field`. Same family as permutation degenerates (mathematically the substrate produces the canonical positive's metric value trivially) but not a within-substrate transform of the positive, so the v1.2 invariance flags don't catch it. Affects uniformity-sensitive detectors (Kuramoto-like). **v1.3 candidate** — chat-led review. Do not modify spec or detector this sprint.
- **C-p14-class-c-borderline (Sprint 33):** Persists at `p_diss=0.350`. Class C TNR 0.900 is at the threshold but not below, so doesn't gate P14's verdict. Low priority.
- **C-p27-time-shuffle-invariance (Sprint 34):** Still provisional. **Sprint 36 P27 panel run is the natural validation.**
- **C-supplements (Sprint 31):** Still OPEN for lattice_1d / lattice_2d_continuous / scalar_wealth / opinion_space. Not needed for Sprint 36's recommended P22+P27 batch (both lattice_2d).
- **C-pyproject-pin (Sprint 29):** Still OPEN. 1-line `pyproject.toml` change.

## HEAD commit hash and tag at end of sprint

To be recorded after commit + push + tag:

- **Commit:** `__TBD__`
- **Tag:** `v0.35.0` (pushed to origin)
