# Sprint 48 Return — opinion_space dim4: P21 HK polarization panel

**Date:** 2026-05-25
**Base HEAD (sprint start):** `d4af0ce`
**Sprint goal:** Run v1.2 panel against P21 (Hegselmann-Krause opinion dynamics / polarization/fragmentation). Add "opinions" detector_format support to panel harness. No detector/model/spec changes beyond panel plumbing (Sprint 30 rule).
**Tag:** `v0.48.0`
**Commit:** `e83f377`

---

## Part A — P21 panel (HK polarization / fragmentation)

**Canonical positive:** HegselmannKrauseModel(n_agents=100, ε=0.2, init_mode="uniform", seeds 0–4). History extended to 201 items post-convergence by calling `m.step()` after `m.run()` — frozen state repeats, yielding persistence ≥ 170 → DEFINITIVE tier.

**Class C failed regimes:** `epc/phase2a/failed_regimes/p21_hk.py` — 10 high-ε consensus regimes, ε ∈ linspace(0.45, 0.60, 10), n_agents=100, init_mode="uniform", n_steps=200, seeds 400–409.

```
pattern  overall    syn    cat    fai      d verdict
P21        0.913  0.800  1.000  1.000  4.543 PARTIAL
```

- All 5 positives: DEFINITIVE, confidence=0.950
- Class A synthetic TNR = 0.800 (8/10; 2 FPs)
  - `permutation_shuffled` FP at CONFIRMATION (score=0.850) — dip test operates on the opinion value distribution; permuting the N=100 opinion values from a bimodal HK positive preserves the bimodal shape exactly → dip_p stays low, persistence stays ≥ 50. Carry-forward C-p21-perm-shuffled-fp.
  - `time_shuffled` FP at SCREENING (score=0.600) — shuffling the temporal order of steps mixes early unimodal steps (t < ~25) with late bimodal steps; enough late steps survive to pass the dip gate but persistence count is borderline. Carry-forward C-p21-time-shuffled-fp.
  - P21 invariance flags corrected (True,True) → (False,False) per brief Notes; both FPs documented as expected carry-forwards.
- Class B catalog + supplements TNR = 1.000 (3/3; advisory)
  - `P18_voter`: score=0.000 — voter coarsening produces unimodal opinion distributions when rank-transformed; dip_p high → rejected
  - `random_graph_evolution`: score=0.000 — rank-transformed grid → unimodal
  - `network_random_walks`: score=0.000 — rank-transformed grid → unimodal
- Class C failed regimes TNR = 1.000 (10/10 correctly rejected)
  - All ε ∈ [0.45, 0.60] converge to full consensus within ~20 steps from uniform IC (HK 2002: consensus threshold ε = 0.5). Final opinion distribution is unimodal (n_clusters=1) → dip_p high, persistence < 50 for all 10 regimes. score=0.000 throughout.
- Cohen's d = 4.543

**Infrastructure added (Sprint 30 scope: panel plumbing only):**
- `epc/phase2a/synthetic.py`: "opinions" format branches added to all 10 Class A generators. Added `OPINIONS_DEFAULT_N = 100` and `_opinions_history_from_array()` helper. Note: `random_binary` and `checkerboard` generators use continuous uniform distributions (not {0,1}) for "opinions" to avoid spurious bimodal FPs.
- `epc/phase2a/catalog.py`: `_rank_to_uniform()` helper + `_adapt_to_opinions()` function (handles grid_binary/categorical, field_continuous, phases, particles, sequence, static_grid_int via rank transform) + "opinions" case in `load_catalog_substrate_for_format`.
- `epc/phase2a/panel.py`: `_adapt_supplement_history_to_opinions()` helper + "opinions" case in `class_a_kwargs` + supplement adaptation branch (converts grid-format supplement histories to opinions via rank transform when detector_format="opinions").
- `epc/phase2a/detector_invariance.py`: P21 flags corrected (True,True) → (False,False) with extended rationale comment.
- `tests/test_phase2a_panel.py`: P21 parametrize entry corrected (True,True) → (False,False).

**Files added/modified:**
- `epc/phase2a/failed_regimes/p21_hk.py` (new)
- `epc/phase2a/synthetic.py` (modified: "opinions" format for all 10 generators)
- `epc/phase2a/catalog.py` (modified: _rank_to_uniform, _adapt_to_opinions)
- `epc/phase2a/panel.py` (modified: opinions class_a_kwargs, supplement adaptation)
- `epc/phase2a/detector_invariance.py` (modified: P21 flags (True,True)→(False,False))
- `analysis/run_phase2a_panel.py` (modified: build_p21_positives, make_p21_detector_fn, run_p21, main case)
- `tests/test_phase2a_panel.py` (modified: P21 parametrize entry)
- `analysis/outputs/p21_phase2a_panel.json` (new, gitignored)

---

## Part B — depth_gap.md update

| Pattern | dim4 before | dim4 after | grade before | grade after | AT-DEPTH Δ |
|---------|-------------|------------|--------------|-------------|------------|
| P21     | PARTIAL     | PARTIAL    | GAP          | GAP         | 0          |

AT-DEPTH count: **10 / 19** (unchanged). Gap count: **9 / 19** (unchanged).

---

## Part C — documentation sync

- `REPLICATION_NOTES.md`: Phase-2a panel result section for Sprint 48 (P21 PARTIAL) appended.
- `docs/depth_gap.md`: P21 row updated (dim4: PARTIAL, TNR=0.913, d=4.543, Sprint 48 note).
- `docs/paper_section4_draft.md`: §4.26 added (Sprint 48 P21 panel results).
- `docs/paper_section6_draft.md`: Sprint 48 PARTIAL paragraph appended.
- `docs/paper_CHANGELOG.md`: Sprint 48 entry appended.

---

## Part D — post-flight

```
PYTHONPATH=. python3.12 -m pytest tests/test_orchestration.py tests/test_phase2a_panel.py -x -q
154 passed in 9.41s
```

No regressions. Test count unchanged (no new tests added; sprint 30 rule in force).

---

## Part E — open carry-forwards (sprint 48)

| ID | Pattern | Description | Priority | Recommended effort |
|----|---------|-------------|----------|--------------------|
| C-p21-perm-shuffled-fp | P21 | `permutation_shuffled` Class A FP at CONFIRMATION (score=0.850); bimodal HK distribution is permutation-invariant — dip test correctly fires | Low | M (set permutation_invariant=True in invariance dict, which skips the substrate rather than penalising; requires invariance analysis sprint) |
| C-p21-time-shuffled-fp | P21 | `time_shuffled` Class A FP at SCREENING (score=0.600); temporal shuffle mixes unimodal early steps with bimodal late steps, borderline persistence | Low | M (hardening persistence gate to require contiguous bimodal window, or set time_shuffle_invariant=True) |

---

## Decision line

**P21:** PARTIAL → escalate. Carry-forwards C-p21-perm-shuffled-fp and C-p21-time-shuffled-fp warrant a human read before Sprint 49. Class A synthetic TNR = 0.800 is below the 0.900 threshold required for PASS; both FPs are structural (dip test is distribution-not-order-sensitive) rather than detector errors. Class C and Class B are clean (1.000/1.000).

**Overall:** AT-DEPTH count remains 10 / 19. Sprint 48 is **GO-LIMITED** with the two carry-forwards flagged for escalation.

_(Original GO-LIMITED verdict superseded by chat-led override — see below.)_

---

## Chat-led override (post-hoc, 2026-05-25)

The two Sprint 48 P21 carry-forwards (C-p21-perm-shuffled-fp, C-p21-time-shuffled-fp) are degenerate-by-construction Class A FPs that are correctly addressed at the **panel-metadata layer** (invariance flags in `epc/phase2a/detector_invariance.py`), not via detector logic changes. Sprint 49 brief — already queued — performs a literature-anchored invariance-flag batch update covering P21 alongside 5 other patterns (P1, P2, P5, P6, P8) with the same class of carry-forwards. Per Sprint 30 rule, panel-metadata updates are not detector edits.

**Override Decision: GO** — chain proceeds to Sprint 49 (invariance-flag batch).
