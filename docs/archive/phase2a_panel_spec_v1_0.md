# Phase-2a Standard Negative Panel — Specification

This document specifies the standard substrate-diverse negative panel used to score Dimension 4 (broad negative sweep / specificity) of the depth-gap rubric. Every detector in the catalog is tested against this panel as part of its Phase-2a validation.

The panel is the methodological backbone that closes the catalog's largest depth gap (audit C1: 15/19 patterns PARTIAL on Dim 4).

---

## Composition

The panel contains **30 negative substrates** across three classes of 10 each. The synthetic and catalog classes are fixed across all patterns; the failed-regime class is pattern-specific (the substrates in it depend on which detector is being tested).

### Class A — Synthetic null substrates (10, fixed)

Substrates with no emergent pattern structure, generated from null processes. These test the most basic specificity claim: the detector should not fire on noise.

The 10 synthetic substrates are:

1. **Random uniform field** — values drawn i.i.d. uniform on the detector's natural range, no spatial or temporal structure
2. **Random Gaussian field** — values i.i.d. Gaussian, mean and variance matched to the canonical positive's metric range
3. **Random binary field** — i.i.d. Bernoulli(0.5) on a grid, the simplest binary substrate
4. **Spatial white noise time series** — independent random fields at each time step (no temporal autocorrelation)
5. **Temporal white noise per cell** — each cell evolves as an independent random walk (no spatial structure)
6. **Permutation-shuffled positive** — take the canonical positive's final state and randomly permute cells, destroying spatial structure but preserving the marginal value distribution
7. **Time-shuffled positive** — take the canonical positive's full trajectory and shuffle time steps, destroying temporal structure
8. **Constant field** — all cells fixed to a single value (tests against detectors that fire on triviality)
9. **Linear gradient field** — smooth monotonic spatial gradient, no emergent structure but non-trivial spatial correlation
10. **Periodic boundary checkerboard** — alternating values, deterministic structure but not the target pattern

These 10 are identical for every detector tested against the panel.

### Class B — Catalog-derived non-positives (10, fixed)

The canonical positive states of 10 *other* patterns from the catalog, run for the same duration and grid size as the detector under test (where compatible). For a detector D_i testing pattern P_i, Class B contains the canonical positives of 10 other patterns — these are *not* P_i but they are emergent patterns, so the test is whether D_i discriminates between its own pattern and other emergent phenomena.

The 10 catalog-derived substrates are the canonical positives of: **P1 (Schelling), P3 (Gray-Scott), P5 (Vicsek), P9 (Kuramoto), P10 (chimera), P14 (BTW sandpile), P15 (GoL), P18 (voter), P27 (Nowak-May), P31 (Zhang sorting)**. These were chosen to span all 11 ontological dimensions of the catalog with maximal coverage at minimum count.

When a detector D_i is testing pattern P_i and P_i is in this list, that entry is replaced with the canonical positive of **P12 (RPS)** as a fallback (so the test always uses 10 non-self positives).

### Class C — Failed-regime substrates (10, pattern-specific)

The detector's *own* model run in 10 parameter regimes that do NOT produce the target pattern. This is the strongest test: same model, same substrate type, but a regime where the pattern fails to emerge.

The 10 failed regimes are pattern-specific and must be specified per detector when the panel is applied. Each pattern's panel-application sprint specifies its 10 failed regimes (e.g., for Vicsek: 10 high-noise η regimes above the order-disorder transition; for Schelling: 10 low-tolerance regimes below the segregation threshold; for Kuramoto: 10 low-coupling K regimes below K_c).

**Default if pattern-specific regimes are not yet known:** use 10 evenly-spaced parameter values within the documented "no pattern" region of the model's phase diagram, drawn from a single dimension if the phase diagram is 1D, or from a Latin Hypercube sample if the phase diagram is multi-dimensional.

---

## PASS criterion

A detector passes Dim 4 against the panel when:

1. **Overall TNR ≥ 95%** across all 30 substrates. (At most 1 false positive out of 30.)
2. **Per-class TNR is reported** explicitly: synthetic TNR, catalog TNR, failed-regime TNR.
3. **Per-class TNR ≥ 90% is the soft expectation** (no class shows ≤80% in particular). If any class drops below 90%, the panel result is reported as PASS-with-weakness with the weak class flagged in the detector card. Soft because some patterns may have a class where 100% rejection is easy and another where 90% is a real achievement; we want the data, not a binary gate.
4. **Effect size on the primary metric**: Cohen's d ≥ 1.0 between the canonical positive and the pooled panel, reported alongside TNR.

A detector that meets (1) and (2) and reports (3) and (4) is scored PASS on Dim 4. A detector that meets (1) but has a per-class weakness is PASS-with-weakness. A detector that does not meet (1) is PARTIAL on Dim 4 (or FAIL if no panel run exists at all).

---

## Harness output

When the panel is run for a detector D_i, the harness produces:

1. A JSON results file at `analysis/outputs/p<i>_phase2a_panel.json` with structure:
   ```
   {
     "pattern_id": "P9",
     "head_commit": "...",
     "panel_version": "1.0",
     "canonical_positive": {"score": ..., "verdict": "DEFINITIVE"},
     "synthetic": [{"substrate": "random_uniform", "score": ..., "verdict": "..."}, ...],
     "catalog": [{"substrate": "P1_schelling", "score": ..., "verdict": "..."}, ...],
     "failed_regime": [{"substrate": "K=0.5_below_Kc", "params": {...}, "score": ..., "verdict": "..."}, ...],
     "summary": {
       "overall_tnr": 0.967,
       "synthetic_tnr": 1.0,
       "catalog_tnr": 1.0,
       "failed_regime_tnr": 0.9,
       "cohens_d_positive_vs_panel": 2.34,
       "verdict": "PASS"
     }
   }
   ```

2. A summary section appended to the pattern's entry in `REPLICATION_NOTES.md` with the verdict, TNRs, and Cohen's d.

3. An update to the pattern's row in `docs/depth_gap.md` flipping Dim 4 from PARTIAL to PASS (or PASS-with-weakness), with a one-line note pointing to the panel JSON.

---

## Panel versioning

This panel is **v1.0**. Subsequent revisions (adding substrates, removing substrates, changing PASS criteria) get a version bump and a note in the panel spec changelog. Detector results record the panel version they were tested against, so re-running becomes a tractable depth-gap question rather than a silent invalidation.

---

## Implementation location

- Harness module: `epc/phase2a/panel.py`
- Synthetic substrate generators: `epc/phase2a/synthetic.py`
- Catalog-derived substrate runner: `epc/phase2a/catalog.py`
- Pattern-specific failed-regime configs: `epc/phase2a/failed_regimes/<pattern_id>.py`
- Test harness: `tests/test_phase2a_panel.py`

The harness is reusable across patterns; only the failed-regime config and the detector being tested change per pattern.
