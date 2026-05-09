# Sprint 24 Phase 1 — Baseline characterization for #20b recovery

**Sprint:** 24 (science sprint, focused on carry-forward #20b)
**Phase:** 1 (baseline + grading) and 2 (implementation + verification)
**HEAD at characterization:** `3285f3cf918a2f5174fe9b36ac0fb79c44bda2eb` (Sprint 23 D6)
**Status:** Both phases complete. C6 implemented; #20b closed.

## Purpose

Sprint 21 carry-forward #20b documented that Schelling at threshold = 0.5
reaches P18 DEFINITIVE on all 5 characterized seeds (a metric-level false
positive). The transfer prompt's Option A workflow specifies a four-step
characterization-first process: lock in a baseline, evaluate the four
candidate fixes against that baseline, then implement the chosen
candidate. This document records the baseline and the candidate-fix
grading.

**Discipline:** No detector code, no orchestration code, no registry
was modified during Phase 1. Pure read-only characterization. All numbers
derive from the saved baseline JSON; the candidate-fix grading is a
synthetic dry-run (apply the candidate's threshold rule to already-
collected metrics, predict the new tier).

## What was characterized

**Voter (canonical positive):** L ∈ {64, 128, 256} × seeds {0, 1, 42,
200, 500} = 15 runs. n_steps follows the existing
`TestSprint20SlowReplication` convention: {64: 400, 128: 400, 256: 300}.

**Schelling (negative across parameter regimes):** thresholds {0.30,
0.375, 0.43, 0.5} × seeds {0, 1, 2, 3, 4} = 20 runs at L = 64, density =
0.9, n_steps = 300. Each Schelling run is graded twice (once with
realistic Schelling metadata `{'threshold': τ, 'density': 0.9}`, once
with `model_metadata=None`); the two paths produce identical tier
outcomes (confirming the metadata path doesn't alter the false
positive). Recommendation grading uses the with_metadata path
exclusively.

Detector: `P18ConsensusDetector(n_permutations=199, seed=0)` at HEAD.

Total runs: 35 model runs, 50 detector evaluations. Wall-clock cost ~5
minutes.

## Headline findings

### Finding A — voter holds DEFINITIVE at every L and seed

All 15/15 voter runs reach DEFINITIVE with `null_p ≤ 0.005` and
`P1=excluded`. The canonical positive baseline is unchanged from
Sprint 20 / Sprint 23 expectations. No regression.

### Finding B — Schelling thr = 0.43 is *also* a false positive (NEW)

Sprint 21's enumeration tested {0.30, 0.375, 0.5} only. Phase 1 added
threshold = 0.43 and discovered: **threshold = 0.43 produces tier
outcomes bit-for-bit identical to threshold = 0.5** across all 5 seeds.
The reason is geometric: at full neighborhood (8 neighbors), the
possible same-fraction values are {0/8, 1/8, ..., 8/8}. None falls in
the half-open interval [0.43, 0.5). So thresholds 0.43 and 0.5 are
dynamically equivalent at full-neighborhood positions. This generalizes
the false positive: it is not isolated to threshold = 0.5 but extends
across the full half-open interval **τ ∈ (0.375, 0.5]** because of
Schelling's discrete neighbor counting.

This generalization strengthens the case for a real fix; the false
positive is broader than the Sprint 21 carry-forward suggested.

### Finding C — confirmed canonical-threshold separation at 0.375

At threshold = 0.375 (Schelling's original), 4 of 5 seeds fail
screening (`moran_final_qtr` ∈ [0.239, 0.285], below the 0.30
screening floor); 1 seed (seed=2) passes screening with
`moran_final_qtr` = 0.301 but stops at SCREENING because
`wall_final_qtr` = 0.356 fails the 0.30 confirmation ceiling. Exactly
as Sprint 21 documented. No regression.

### Finding D — threshold = 0.30 fails screening cleanly

All 5 seeds at threshold = 0.30 fail screening with
`moran_final_qtr` ∈ [0.203, 0.259]. Below the screening floor. No
discrimination concern in this regime.

### Finding E — voter wall_final creep with L

Voter `wall_final_qtr_mean` increases with L:
- L = 64:  range [0.153, 0.201], mean 0.175
- L = 128: range [0.187, 0.238], mean 0.214
- L = 256: range [0.213, 0.234], mean 0.222

This is consistent with longer time-to-plateau at larger L (the run
length scales modestly with L; the Sprint 20 finite-size convention
uses L = 256 with only 300 sweeps). Schelling thr ∈ {0.43, 0.5} has
`wall_final_qtr` ≈ 0.275. The gap between voter L = 128 max (0.238)
and Schelling (0.275) is only 0.037 — meaning **a wall-ceiling-based
discriminator (candidate C1) has thin headroom at L ≥ 128.**

By contrast, the gap on `moran_final_qtr_mean` is much cleaner:
- Voter min across all 15 runs: 0.499 (at L=64 seed=0)
- Schelling thr ∈ {0.43, 0.5} max across all 10 runs: 0.410 (at thr=0.43 seed=3)
- Separation: **0.089**, with clean midpoint at 0.45.

This suggests `moran_final` discrimination is far more robust than
`wall_final` discrimination for this pair.

## Candidate-fix evaluation

Six candidates evaluated (the original Sprint 21 four plus two
combinations surfaced by the empirical baseline):

| Cand | Description | Voter | Sch thr=0.43 | Sch thr=0.5 | Notes |
|------|-------------|-------|--------------|-------------|-------|
| C1 | wall ceiling 0.30 → 0.25 | 15/15 DEF ✓ | 5/5 SCR | 5/5 SCR | thin margin (+0.012 at L=128 max) |
| C2 | P1-aware DEF downgrade | 15/15 DEF ✓ | 5/5 CONF | 5/5 CONF | architectural; doesn't change confidence-tier outcome below DEF |
| C3 | Schelling 'update' token alone | 15/15 DEF ✓ | 5/5 DEF ✗ | 5/5 DEF ✗ | NO EFFECT (definitive doesn't consult exclusions) |
| C2+C3 | C3 plumbing + C2 gate | 15/15 DEF ✓ | 5/5 CONF | 5/5 CONF | same outcome as C2 alone |
| C4 = C1+C2 | wall + P1-aware | 15/15 DEF ✓ | 5/5 SCR | 5/5 SCR | thin margin (C1 component) |
| **C5** | moran floor 0.30 → 0.45 | **15/15 DEF ✓** | 5/5 CONF | 5/5 CONF | **clean +0.049 margin both sides** |
| **C6 = C5+C2** | moran floor + P1-aware | **15/15 DEF ✓** | 5/5 CONF | 5/5 CONF | **clean margin + architectural defense in depth** |

Key observations:

1. **C3 standalone is a no-op.** Adding a registry 'update' token to
   Schelling without changing the tier-determination logic does
   nothing because `_check_definitive` does not consult exclusion
   results. C3 only matters as plumbing for C2.

2. **C1 (and C4) tighten a metric ceiling with thin headroom.**
   Voter L = 128 seed = 42 has `wall_final_qtr_mean` = 0.238, just
   0.012 below a 0.25 ceiling. Across more seeds at higher L, voter
   could plausibly breach 0.25; the candidate trades a known false
   positive for a probable future false negative.

3. **C5 is the cleanest single-metric fix.** A 0.089 separation
   between voter and Schelling on `moran_final_qtr_mean`, with ~0.05
   margin on both sides of a 0.45 floor. Empirically much safer than
   C1.

4. **C2 (and C6) honestly demote to CONFIRMATION rather than to
   SCREENING.** C1/C4 demote Schelling to SCREENING, claiming the
   confirmation gates fail. But Schelling's wall-decay and Moran
   growth ARE genuinely there at thr ∈ {0.43, 0.5} — calling the
   confirmation gates failed is, strictly, wrong. The honest demotion
   is to CONFIRMATION (metrics consistent, but neighbors not
   excluded). The three-tier system was designed for exactly this
   gradation.

5. **C2 has a hidden architectural payoff.** The current code at
   `epc/base_detector.py:333` writes
   `bonuses["all_exclusions_cleared"] = True` *before* exclusions are
   actually checked, with a comment "Updated after exclusion check"
   that is currently a lie (no later code overwrites it). C2 forces
   `_check_definitive` to consult the exclusion result, which makes
   that bonus assertion truthful and tightens the contract. This is a
   latent bug fix beyond just the false positive.

## Recommendation

**Implement C6 (= C5 + C2):** raise `DEFINITIVE_MORAN_FINAL_MIN` from
0.30 to 0.45, AND modify `_check_definitive` to require
`exclusion_results['P1'] == 'excluded'` for the DEFINITIVE tier.

Rationale:
- Empirically clean margins on voter (≥ +0.049 on `moran_final`).
- Honest demotion of Schelling thr ∈ {0.43, 0.5} to CONFIRMATION
  (metrics-consistent, neighbors-not-excluded), which is the tier
  designed for exactly this case.
- Architecturally sound: closes the latent
  `all_exclusions_cleared` bug at the same time.
- Avoids the C1 thin-margin risk at L ≥ 128.
- Robust against future Schelling parameters (any τ that produces
  `moran_final < 0.45` would be metric-rejected by C5; any τ that
  somehow produced `moran_final ≥ 0.45` without P1 exclusion would
  be architecturally rejected by C2).

Rejected alternatives:
- **C3 alone**: no-op against the false positive.
- **C1 alone**: thin headroom at L ≥ 128; trades known FP for likely
  future FN.
- **C2 alone**: only demotes to CONFIRMATION, but does not close the
  metric gap; if a future change touches the exclusion logic, the
  fix could regress.
- **C4 (C1+C2)**: same C1 thin-margin risk + dishonest "fails
  confirmation" outcome (Schelling's confirmation gates DO pass at
  baseline thresholds; demoting to SCREENING claims they don't).
- **C5 alone**: empirically equivalent to C6 on this characterization
  ensemble, but lacks the architectural defense in depth. C6's
  marginal cost over C5 is ~5 lines of code and is worth it.

The transfer prompt's option 4 ("some combination of the above")
explicitly anticipated this kind of recommendation; C6 is a principled
combination that the original three-candidate framing did not name but
which the empirical baseline supports.

## Note on `DEFINITIVE_MORAN_FINAL_MAX = 0.75`

C5 raises the lower bound but does NOT change the upper bound. The
current 0.75 ceiling rejects GH broken_wave (`moran_final` ≈ 0.87,
Sprint 20 §4.20) and is unchanged. The new window [0.45, 0.75] is
narrower than the old [0.30, 0.75] but still covers voter's empirical
range [0.499, 0.663] with margin on both sides.

## Note on the `_all_secondaries_pass` confirmation bonus

`P18ConsensusDetector._all_secondaries_pass` at line 564 does NOT
consume `DEFINITIVE_MORAN_FINAL_MIN`; it only checks wall_spearman,
wall_final, and minority. C5 affects only the DEFINITIVE gate, so the
confirmation tier is unchanged. This is the correct compositional
behavior.

## Phase 2 work — COMPLETED

Phase 2 implemented C6 in code and verified the predictions. All
seven items from the original Phase 2 plan completed:

1. ✓ `epc/detectors/p18_consensus.py`:
   - `DEFINITIVE_MORAN_FINAL_MIN = 0.30` → `0.45` (Decision 57)
   - `_check_definitive` now calls `_check_exclusions` and requires
     all three exclusions to return `"excluded"` (Decision 58)
   - Top-of-file detection-tiers docstring + `_check_definitive`
     docstring updated
2. ✓ Re-ran Phase 1 baseline against modified detector. Outcomes
   match dry-run predictions exactly: voter 15/15 DEFINITIVE,
   Schelling thr ∈ {0.43, 0.5} 5/5 CONFIRMATION (was DEFINITIVE).
3. ✓ Added `TestSprint24Schelling0p5Regression` to
   `tests/test_voter_p18_e2e.py` — 30 parametrized tests across
   three assertions per (threshold × seed) pair.
4. ✓ `docs/detector_cards.md` v0.6.3: tier-gate spec, discriminator
   table (added thr=0.43 row, updated thr=0.5 row), narrative
   rewrite of #20b paragraph, Decisions 57 + 58, P1 exclusion
   bullet.
5. ✓ `docs/paper_section4_draft.md` §4.20: "Sprint 24 update: #20b
   closed via combined gate fix (C6)" paragraph after the Sprint
   21 caveat.
6. ✓ `REPLICATION_NOTES.md`: full Sprint 24 section appended
   (~250 lines) documenting Phase 1 baseline, Findings A–E,
   candidate-fix grading, Phase 2 implementation, verification,
   and Sprint 24 newly surfaced findings (#27 latent contract
   bug, #28 metric-choice principle, #29 dry-run grading workflow).
7. ✓ Pre-flight bundle: 175 → 205 fast tests pass under modified
   detector. No regressions.

## Files generated by Phase 1 + Phase 2

Repo deliverables (land in the named paths via Claude Code):

| Path | Purpose |
|---|---|
| `epc/detectors/p18_consensus.py` | Modified — Decisions 57, 58, docstrings |
| `tests/test_voter_p18_e2e.py` | Modified — +`TestSprint24Schelling0p5Regression` |
| `docs/detector_cards.md` | Modified — v0.6.3, P18 card updates |
| `docs/paper_section4_draft.md` | Modified — §4.20 closure paragraph |
| `REPLICATION_NOTES.md` | Modified — Sprint 24 section appended |
| `docs/sprint24/phase1_baseline.md` | New — this document |
| `docs/sprint24/baseline_voter_schelling.json` | New — 35-run baseline data archive |
| `docs/sprint24/candidate_grades.json` | New — dry-run grading detail |
| `scripts/sprint24_baseline.py` | New — reproducibility: characterization driver |
| `scripts/sprint24_grade_candidates.py` | New — reproducibility: candidate dry-run grader |
