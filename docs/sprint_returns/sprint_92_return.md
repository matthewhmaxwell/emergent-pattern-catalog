# Sprint 92 Return — P30 Autopoiesis (ALL 32 PATTERNS IMPLEMENTED)

**Date:** 2026-06-10
**Base HEAD (sprint start):** `f2cdb7b` (Sprint 91 follow-up)
**Tag:** `v0.92.0`
**Sprint type:** Chat-led design + code-led execution.

## Part A — AutopoiesisModel (`epc/models/autopoiesis.py`)

SCL-style autopoiesis with three particle types in 2D periodic domain:

- **Substrate (S, type=0):** Free-diffusing resource (diffusion=0.3).
- **Catalyst (C, type=1):** Near-stationary production center (diffusion=0.002).
- **Link (L, type=2):** Membrane particle (diffusion=0.01).

Update rules: production (S → L near C), decay (L → S), radial spring
(L attracted to membrane_equilibrium_radius=3.0 from nearest C), tangential
L-L attraction (chain cohesion), steric repulsion. Total particle count
conserved. Catalyst count fixed.

Key design decision: radial spring mechanism (membrane_equilibrium_radius +
catalyst_link_attraction) produces a thin-shell membrane that concentrates
links at a preferred distance from catalysts, giving strong discrimination
under type-shuffle null.

Two negative controls:
- **NonBondingParticleModel:** All substrate, Brownian diffusion only.
- **DenseClusterModel:** P1-like universal attraction, no type differentiation.

## Part B — P30AutopoiesisDetector (`epc/detectors/p30_autopoiesis.py`)

Three-tier detection:
- **Screening:** association_score > 1.5 AND closure_fraction > 0.5 AND
  mean_n_links ≥ 3.
- **Confirmation:** + enrichment_ratio > 1.2 + null_p < 0.01 + persistence > 0.5.
- **Definitive:** + closure > 0.7 + enrichment > 2.0 + persistence > 0.8 +
  link_cv < 0.3.

Primary metric: **association_score** — fraction of link particles within
association_radius of any catalyst, divided by CSR expectation.

Null model: type-shuffle permutation (keep all positions, randomize type
labels S/C/L). Under null, association_score ≈ 1.0; under autopoiesis,
association_score >> 1.

P1 exclusion: association_score measures spatial co-location of *specific
types* (links near catalysts), not generic clustering.

T1a adapter: `extract_observation_bundle()` produces positions, types, bonds,
steps, box_size, n_particles.

## Part C — E2E tests (`tests/test_autopoiesis_p30_e2e.py`)

25 tests, all passing (31.4s):
- TestP30CanonicalDetection: 8 tests (detected, CONFIRMATION tier, closure,
  association, enrichment, null p, confidence, P1 excluded)
- TestP30NegativeControls: 2 tests (non-bonding, dense cluster rejected)
- TestP30Registration: 5 tests (orchestration registry)
- TestP30ObservationBundle: 2 tests (T1a adapter keys + shapes)
- TestP30Determinism: 2 tests (seeded reproducibility)
- TestP30ModelState: 6 tests (type validity, catalyst conservation, particle
  count, bounds, link production, metadata)

## Part D — dim1 reproduction

`analysis/outputs/p30_scl_reproduction.json`. Varela-Maturana-Uribe 1974.

| Sub-signature | Metric | Value | Threshold | Result |
|---------------|--------|-------|-----------|--------|
| Closure | closure_fraction | 1.000 | > 0.5 | **PASS** |
| Gradient | enrichment_ratio | 1.997 | > 1.2 | **PASS** |
| Self-repair | recovery_fraction | 1.102 | > 0.7 | **PASS** |

Detection: CONFIRMATION (association_score=2.301, p=0.005, confidence=0.70).

## Part E — dim2 multiseed + dim3 methods note

**dim2:** 20-seed campaign at canonical regime:

| Metric | Mean | Std | CV |
|--------|------|-----|-----|
| association_score | 2.211 | 0.087 | 3.9% |
| closure_fraction | 0.996 | 0.013 | 1.3% |
| enrichment_ratio | 1.687 | 0.586 | 34.7% |

Detection: 20/20 detected, 14/20 confirmation+, 5/20 definitive.

**dim3:** `docs/methods_notes/p30_methods.md` — SCL particle dynamics,
three detection metrics, type-shuffle null, three-tier gates, T1a bundle,
limitations (angular-coverage necessary-not-sufficient, partial self-repair
per escape clause).

## Part F — Documentation updates

- `epc/orchestration.py`: autopoiesis model + P30 detector registered.
  Substrate: particle_membrane (17th). Registry: 42 models × 32 detectors.
- `docs/observation_schema.md`: Particle-membrane bundle (P30) added.
- `docs/depth_gap.md`: P30 row added (dim1–dim3 PASS, dim4 pending → GAP).
  Patterns audited: 32. AT-DEPTH: 30/32. Gap: 2 (P12, P30).
- `docs/paper_section4_draft.md`: §4.30 added.
- `docs/paper_section6_draft.md`: Sprint 92 entry added.
- `docs/paper_CHANGELOG.md`: Sprint 92 entry added.
- `REPLICATION_NOTES.md`: Sprint 92 section added.

## Files changed

**New (6):**
- `epc/models/autopoiesis.py` — AutopoiesisModel + 2 negative controls
- `epc/detectors/p30_autopoiesis.py` — P30 detector + T1a adapter
- `tests/test_autopoiesis_p30_e2e.py` — 25 e2e tests
- `analysis/reproductions/p30_scl.py` — dim1 reproduction script
- `analysis/reproductions/p30_multiseed.py` — dim2 campaign script
- `docs/methods_notes/p30_methods.md` — dim3 methods note

**Modified (7):**
- `epc/orchestration.py` — autopoiesis + P30 registration
- `docs/observation_schema.md` — P30 bundle
- `docs/depth_gap.md` — P30 row + counts
- `docs/paper_section4_draft.md` — §4.30
- `docs/paper_section6_draft.md` — Sprint 92 entry
- `docs/paper_CHANGELOG.md` — Sprint 92 entry
- `REPLICATION_NOTES.md` — Sprint 92 section

## Carry-forwards

- **C-p30-dim4**: P30 dim4 (Phase-2a panel) pending. Effort: S.
- **C-p12-dim1**: P12 (RPS λ ∝ √M scaling) remains dim1 GAP.
  Documented finite-size measurement limitation (Sprints 54/58/59/63).
- **C-p30-enrichment-cv**: Enrichment ratio CV=34.7% across seeds — highest
  variance metric. Driven by stochastic catalyst positioning relative to
  membrane. Association_score CV=3.9% is robust.

## Summary

| Metric | Value |
|--------|-------|
| New files | 6 |
| Modified files | 7 |
| Sprint-specific tests | 25/25 PASS |
| dim1 reproduction | PASS (3/3 sub-signatures) |
| dim2 detection rate | 20/20 (100%) |
| dim2 confirmation+ rate | 14/20 (70%) |
| AT-DEPTH count | 30 / 32 |
| Gap count | 2 (P12, P30) |
| **Patterns implemented** | **32 / 32 (100%)** |

**Decision: GO** — P30 implemented at CONFIRMATION tier with robust 20/20
detection across seeds (association_score CV=3.9%). All three dim1 sub-signatures
PASS (closure=1.0, enrichment=2.0, self-repair=1.10). 25/25 e2e tests PASS.
Both negative controls correctly rejected. dim4 (Phase-2a panel) deferred to
next sprint. **ALL 32 PATTERNS NOW IMPLEMENTED.** Chain may proceed to Sprint 93
(P30 dim4 panel).
