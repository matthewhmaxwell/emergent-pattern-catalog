# Depth-Gap Audit Rubric

This rubric defines the depth standard every implemented pattern (P1–P32) must meet for the Emergent Pattern Catalog to be publication-ready. Each pattern is scored on four dimensions, each scored as **PASS / PARTIAL / FAIL / N/A**.

The rubric is the standard. The Sprint 29 audit applies it to the currently-implemented patterns and produces a depth-gap matrix. Subsequent sprints close the gaps it surfaces.

---

## Dimension 1 — Quantitative replication of a published canonical paper

**Required for:** every pattern.

A pattern passes this dimension when the repository contains a canonical positive model that quantitatively reproduces a specific result from a peer-reviewed publication, with the comparison documented in `REPLICATION_NOTES.md`.

| Score | Criterion |
|---|---|
| PASS | A specific published paper is named. A specific quantitative result (figure, table, or named statistic) from that paper is reproduced within stated tolerance. The reproduction is documented in `REPLICATION_NOTES.md` with the canonical citation, the figure/table reference, the tolerance, and the comparison values. |
| PARTIAL | A canonical paper is named and the model qualitatively matches its behavior, but no specific quantitative result is reproduced, OR the comparison is not documented in `REPLICATION_NOTES.md`. |
| FAIL | No canonical paper is named, OR the model does not match the paper's behavior. |
| N/A | Not applicable. (No exceptions: every pattern needs a canonical paper.) |

---

## Dimension 2 — Multi-seed validation

**Required for:** stochastic models and any pattern whose canonical positive depends on a phase-boundary determination. Optional for fully deterministic models with no phase-boundary dependence.

A pattern passes this dimension when, where required, validation runs use ≥5 seeds and the result is reported with appropriate dispersion (variance, confidence intervals, or basin-fraction).

| Score | Criterion |
|---|---|
| PASS | The canonical positive (and any phase-boundary determination) was validated across ≥5 seeds, and dispersion is reported in `REPLICATION_NOTES.md` (CI, variance, or basin fraction). |
| PARTIAL | Multi-seed validation exists (≥2 seeds) but does not meet the ≥5-seed bar, OR dispersion is not reported. |
| FAIL | Single-seed only, where multi-seed is required by the model's stochasticity or phase-boundary dependence. |
| N/A | Fully deterministic model with no phase-boundary dependence — multi-seed not required. |

**How to decide N/A:** the auditor checks the model source. If the model uses any randomness (initial conditions, transition rules, parameter draws) OR if the canonical positive's classification depends on parameter values near a phase boundary, multi-seed is required. Otherwise N/A is appropriate. Default to required if uncertain.

---

## Dimension 3 — Methods note

**Required for:** every pattern.

A pattern passes this dimension when `REPLICATION_NOTES.md` contains a written methods note explaining the specific implementation choices made relative to the canonical paper. Examples: unit conventions, discretization scheme, boundary conditions, time-step, integration method, RNG choice. The methods note must be substantive enough that a reader could re-derive the implementation decisions.

| Score | Criterion |
|---|---|
| PASS | A methods note exists in `REPLICATION_NOTES.md` covering: units/conventions, discretization or rule choice, any deviation from the canonical paper, and the rationale. The note is specific and substantive (not boilerplate). |
| PARTIAL | A methods note exists but is incomplete (e.g., covers some choices but not all material ones), OR is boilerplate / non-substantive. |
| FAIL | No methods note, OR only inline code comments without a written note in `REPLICATION_NOTES.md`. |
| N/A | Not applicable. (No exceptions: every pattern needs a methods note.) |

---

## Dimension 4 — Broad negative sweep

**Required for:** every pattern.

A pattern passes this dimension when the detector has been validated against a broad negative set (Phase 2a-style sweep) of substrate-diverse non-positives, with results documented in `REPLICATION_NOTES.md` or `docs/detector_cards.md`. The cross-detection transfer matrix is necessary but not sufficient — a Phase 2a-style sweep specifically tests detector specificity against models *not* in the catalog.

| Score | Criterion |
|---|---|
| PASS | A documented Phase 2a-style negative sweep exists. The sweep covers substrate-diverse non-positives (≥3 substrate types where applicable). Specificity (≥95% true negative rate, or stated rationale for lower) is reported. |
| PARTIAL | Negative testing exists but is limited to the cross-detection transfer matrix, OR the sweep covers fewer than 3 substrate types, OR specificity is not reported. |
| FAIL | No negative testing beyond the canonical positive's complement. |
| N/A | Not applicable. (No exceptions: every pattern needs broad negative testing.) |

---

## Overall pattern grade

A pattern is **at depth** when all four dimensions are PASS (or N/A where allowed). Anything else is a depth gap.

The audit produces a matrix with rows = patterns, columns = dimensions, cells = scores, plus a column for the overall grade and a notes column for any auditor commentary.

---

## Audit procedure

For each implemented pattern P_i:

1. Locate the detector card (`docs/detector_cards.md`).
2. Locate the canonical positive model in the codebase.
3. Locate the pattern's section in `REPLICATION_NOTES.md`.
4. Locate the pattern's tests in `tests/`.
5. Score each of the four dimensions per the rubric above.
6. Record findings in `docs/depth_gap.md` as a row of the audit matrix, with:
   - pattern ID and name
   - score for each dimension (PASS / PARTIAL / FAIL / N/A)
   - overall grade (AT-DEPTH / GAP)
   - one-sentence note per non-PASS dimension explaining what's missing
   - estimated effort to close the gap (S / M / L)

The audit does not modify any model, detector, or test. It only produces the matrix.
