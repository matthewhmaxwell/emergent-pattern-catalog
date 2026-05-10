# Sprint 32 Return Summary

**Sprint goal:** Add the P21 (Hegselmann-Krause) native generator that closes carry-forward C-p21-generator, then re-run the Phase-2a panel against P18 and P9 under v1.1, updating `docs/depth_gap.md` only for those two patterns based on the v1.1 outcomes. **Status: complete, with a deviation — P9 returned PARTIAL under v1.1, not PASS.**

## Pre-flight verification

- Base HEAD: `de79d6b` (Sprint 31, v0.31.0). Working tree clean. ✓
- `python3.12 scripts/count_transfer_matrix.py`: **20 / 19 / 79 / 274 / 27 / 19 / 361 / 77 / 284** ✓ (unchanged).
- Pre-flight fast tests: 554 passed (Sprint 31 ground truth, re-verified post-flight).

## Part A — P21 (Hegselmann-Krause) native generator

- **File:** `epc/phase2a/catalog.py`. New `_gen_p21_hegselmann_krause(p)` thin wrapper around `epc.models.hegselmann_krause.HegselmannKrauseModel` (the existing canonical-positive model from Sprint 5; **wrapped, not inlined**).
- **Parameters:** `n_agents=400, epsilon=0.2, init_mode="uniform", n_steps=100, seed=0`.
  - ε = 0.2 produces the canonical fragmented outcome (~2 stable opinion clusters from uniform IC) per Sprint 5 replication notes; ε ≥ 0.5 → consensus, ε ≤ 0.1 → many small clusters.
  - N = 400 chosen so the grid adapter can reshape cleanly to a 20×20 occupancy field for P18's grid-format detector under the v1.1 network override.
- **Native "kind" added:** `"opinions"`. Both `_adapt_to_grid` (binarize at 0.5, reshape to 20×20) and `_adapt_to_phases` (multiply by 2π) gain a branch for it.
- **Tests added:** 4 in `tests/test_phase2a_panel.py`:
  - `test_p21_generator_deterministic` — same seed → byte-identical opinion trajectory.
  - `test_p21_canonical_positive_is_fragmented` — confirms ≥ 2 clusters and a non-trivial opinion spread at ε=0.2.
  - `test_p21_native_kind_is_opinions` — output `kind` and metadata fields.
  - `test_p18_class_b_now_includes_p21_generator` — `class_b_for_pattern("P18")["catalog_mates"]` now contains `P21_hegselmann_krause` AND that id is backed by a real generator in `_GENERATORS`.

## Part B — Re-run P18 and P9 under v1.1

### P18 voter — verdict **PASS** (overall TNR = 1.000, Cohen's d = +inf)

| Class | TNR | n | Notes |
|---|---|---|---|
| Synthetic (Class A) | 1.000 | 10 | gating |
| Catalog (substrate-typed: network) | 1.000 | 3 (P21 + 2 supps) | advisory only (n<5) |
| Failed-regime (Class C) | **N/A** | 0 | Sprint 31 declaration; no parameter regime suppresses voter consensus |
| **Overall** | **1.000** | 13 | |

- Class B composition: catalog mate `P21_hegselmann_krause`; B' supplements `random_graph_evolution`, `network_random_walks`. All 3 correctly rejected.
- Cohen's d = +inf (positives uniformly score 0.5, negatives uniformly score 0.0; the degenerate-perfect path in `cohens_d` returns +inf rather than dividing by zero).
- v1.0 → v1.1 delta: catalog 0.900 → 1.000 (the v1.0 P3_gray_scott firing is no longer in P18's substrate-typed Class B); Class C correctly recognised as N/A (eliminating v1.0's largest source of false positives — the init-fraction proxy regimes, which were true positives in disguise).
- Output: `analysis/outputs/p18_phase2a_panel.json`.

### P9 Kuramoto — verdict **PARTIAL** (overall TNR = 0.913, Cohen's d = 3.445)

| Class | TNR | n | Notes |
|---|---|---|---|
| Synthetic (Class A) | 0.800 | 10 | gating; 2 false positives (constant phases, permutation-shuffled positive) |
| Catalog (substrate-typed: oscillator) | 1.000 | 3 (P10 + 2 supps) | advisory only (n<5) |
| Failed-regime (sub-K_c Class C) | **1.000** | 10 | gating |
| **Overall** | **0.913** | 23 | |

- Class B composition: catalog mate `P10_chimera`; B' supplements `incoherent_phases`, `subcritical_kuramoto`. All 3 correctly rejected.
- Cohen's d = 3.445 (positives at 0.350, negatives largely at 0.0 with 2 hits at 0.350).
- v1.0 → v1.1 delta: catalog 0.400 → **1.000** (substrate-typed fix worked exactly as predicted); failed_regime stayed at 1.000; overall 0.733 → 0.913. **Did not reach the 0.95 PASS gate.**
- The 2 residual false positives are degenerate-sync substrates: `constant_field` (all phases at the same value → r=1 trivially; this *is* mathematically synchronization) and `permutation_shuffled_positive` (cell-permuted positive trajectory → r is permutation-invariant so unchanged from canonical). The detector's behavior on both is arguably correct given how Kuramoto's order parameter is defined; the panel's Class A semantics — "no emergent pattern" — is what fails to apply cleanly to oscillator detectors.
- Output: `analysis/outputs/p9_phase2a_panel.json`.
- Per Sprint 30 rule, the detector was **not** modified to engineer a PASS. Logged as **carry-forward C-class-a-oscillator-degenerate** for chat-led panel-spec revision.

## Part C — `docs/depth_gap.md` row updates

- **P18 row:** dim4 stays PASS (was already PASS via content-level discriminators); grade stays AT-DEPTH; notes updated to record the v1.1 panel result as positive confirmation.
- **P9 row:** dim4 stays PARTIAL (panel did not move it to PASS); grade stays GAP; notes updated to record the v1.1 panel still-PARTIAL outcome and the C-class-a-oscillator-degenerate carry-forward.
- **Header counts unchanged:** AT-DEPTH 4, Gap 15. P18 was already AT-DEPTH; P9 did not move out of GAP.

No other rows touched.

## Part D — Substrate-type ground-truth note appended to spec

Appended a new "Note on substrate-type ground truth" subsection at the end of `docs/phase2a_panel_spec.md`. Clarifies that:
1. The substrate-type taxonomy is computed from `epc.orchestration.MODEL_REGISTRY`, not from the spec table.
2. The spec table's illustrative entries (including non-registered P17/P19 and the P3/P11/P28 disagreements) are forward-looking guidance, not overrides.
3. Any divergence between the spec table and the registry is resolved in favor of the registry, with two deliberate exceptions (P18/P21 → `network`) that are encoded as explicit dict assignments in `_build_substrate_type_by_pattern`.

This is an **addendum**, not a version bump (still v1.1).

## Pre-flight + post-flight test counts

- **Pre-flight fast suite:** 554 passed (Sprint 31 ground truth).
- **Post-flight fast suite:** **558 passed, 0 failed, 65 deselected** in 10:43 (matches predicted 554 + 4 P21 tests).
- **Pre-flight bundle:** 205 (unchanged).
- **Transfer-matrix figures:** 20 / 19 / 79 / 274 / 27 / 19 / 361 / 77 / 284 (unchanged).

## Deviations and judgment calls

### Deviation 1 — P9 returned PARTIAL under v1.1, not PASS
The brief and the v1.1 spec migration path expected both P18 and P9 to PASS or PASS-with-weakness under v1.1. P18 PASSed cleanly (TNR=1.000). P9 came in at TNR=0.913 — within striking distance of the 0.95 gate but not over it. Per Sprint 30 strict rule (re-cited in this brief), the detector was **not** modified, and the panel composition was **not** altered beyond what v1.1 spec dictates. The 2 residual false positives are both Class A synthetic substrates that produce mathematically-synchronized phase distributions:
- **constant_field** in phases format: all `theta = 0` → `r = |⟨e^{iθ}⟩| = 1`. This is trivially synchronized.
- **permutation_shuffled_positive** in phases format: shuffles cells of the positive trajectory's final phase array. Since `r` is invariant under permutation of phases, the shuffled state still has `r ≈ r_positive`.
The detector firing on these is arguably correct: they are synchronization states, just degenerate ones. The panel's Class A is supposed to test against substrates that do not exhibit the pattern, but for oscillator detectors the boundary between "no pattern" and "trivially-perfect pattern" is where the synthetic generators bite.

This is a **panel-design question**, not a detector quality question. **Logged as carry-forward C-class-a-oscillator-degenerate** for chat-led review. Suggested resolution paths for v1.2: (a) skip `constant` and `permutation_shuffled` in oscillator-format Class A; (b) replace them with permutations that destroy phase coherence (e.g., add per-phase noise after permutation); (c) accept that oscillator detectors have a documented blind spot at degenerate sync and require a non-degeneracy precondition before screening fires.

### Judgment call — P21 generator wraps the existing model
The brief preferred wrapping over inlining if a model exists. `epc/models/hegselmann_krause.py` was already present, so I wrapped it. The generator is a 12-line thin adapter; no model changes.

### Judgment call — n_agents=400 for HK canonical positive
Chosen to give a clean 20×20 reshape under the grid adapter (`int(sqrt(400)) == 20` exactly, no truncation needed). N=500 (HK module default) would reshape to 22×22 with 16 truncated agents — workable but messier. The canonical fragmented outcome at ε=0.2 is N-insensitive in the relevant range (Sprint 5 used N=500 and reported 2 clusters; my N=400 also produces 2 clusters per the new generator test).

### Judgment call — `opinions` kind grid adapter binarizes at 0.5
Opinions ∈ [0, 1] with HK's bimodal cluster structure typically straddling 0.5 for fragmented IC. Binarizing at 0.5 captures the cluster polarization on a binary grid — exactly what a grid-reading detector like P18 sees natively. Other reasonable choices (binarize at the per-step median, ternary discretization at 0.33/0.67) were considered; 0.5 chosen for simplicity and because HK fragmentation is empirically symmetric around 0.5 with uniform IC.

### Note — interpreter substitution
Same as Sprint 29/30/31: `python3.12 -m pytest`, `PYTHONPATH=. python3.12 ...`. Mac default `python` is unbound.

### Note — `analysis/outputs/p{9,18}_phase2a_panel.json` force-added
The same `outputs/` gitignore line that needed force-add in Sprints 30 and 31 still applies. Used `git add -f` consistently.

## Carry-forward summary

- **C-p21-generator (Sprint 31):** **CLOSED.** `_gen_p21_hegselmann_krause` landed in `epc/phase2a/catalog.py`; P18's Class B catalog-mate is now loadable; P18 v1.1 panel ran cleanly to PASS.
- **C-class-a-oscillator-degenerate (Sprint 32, NEW):** Class A's `constant_field` and `permutation_shuffled_positive` produce mathematically-synchronized phase substrates that any specific Kuramoto-style detector will (correctly) fire on. The panel's Class A semantic ("no emergent pattern") doesn't cleanly apply to oscillator detectors at these substrates. Resolution requires panel-spec revision (likely v1.2). Affects: P9 (this sprint), and will affect P10 (chimera) when its panel is run in Sprint 33+. **Until resolved, oscillator-detector v1.1 panels will plateau at ~0.91 overall TNR and remain PARTIAL.**
- **C-supplements (Sprint 31):** Still OPEN. lattice_1d / lattice_2d_continuous / scalar_wealth supplements still unimplemented. Not needed for any pattern run yet.
- **C-pyproject-pin (Sprint 29):** Still OPEN. 1-line `pyproject.toml` change.

## HEAD commit hash and tag at end of sprint

To be recorded after commit + push + tag:

- **Commit:** `__TBD__`
- **Tag:** `v0.32.0` (pushed to origin)
