# Sprint 3 GitHub Update — Collective Motion (P5/P6)

## What's new

### Models (2 new)
- `epc/models/vicsek.py` — Standard Vicsek Model (1995)
  - 2D periodic box, constant speed, angular noise
  - arctan2 circular mean + sparse-matrix vectorized neighbor averaging
  - 1.2 ms/step at N=300; 4 init modes
  - Phase transition replicated: 6/6 published claims match

- `epc/models/dorsogna_spp.py` — D'Orsogna SPP with Morse potential (2006)
  - Second-order Newtonian, Rayleigh self-propulsion, RK4 integration
  - Pairwise Morse attraction-repulsion, open boundary
  - Published milling parameters: |L|=0.996, ring density confirmed

### Metrics (1 new module, 5 metrics)
- `epc/metrics/collective_motion.py`
  - PolarizationMetric (φ) — P5 primary
  - GroupSpeedRatioMetric (R) — P5 secondary
  - AngularMomentumMetric (L) — P6 primary, P5 exclusion
  - HeadingAutocorrelationMetric — directional persistence
  - HeadingDistributionMetric — nematic order for P7 exclusion

### Detectors (2 new)
- `epc/detectors/p5_flocking.py` — Full 3-tier P5 detector
  - Heading-shuffle null (199 perms), P6/P7 exclusions
  - Vicsek η=0.5 → DEFINITIVE; η=5.0 → not detected

- `epc/detectors/p6_milling.py` — Full 3-tier P6 detector
  - Ring density profile, P5 exclusion
  - D'Orsogna → DEFINITIVE; Vicsek → not detected

### Tests (2 new)
- `tests/test_vicsek_validation.py` — 8 unit + 6 physics tests
- `tests/test_p5_detector_validation.py` — 8 metric + 4 detector tests

### Docs (3 updated)
- `PROJECT_STATUS.md` — Sprint 2 complete, Sprint 3 progress, 23 arch decisions
- `REPLICATION_NOTES.md` — Added Vicsek + D'Orsogna + P5/P6 detector notes
- `CLAUDE.md` — Updated code layout and status

## Cross-detection transfer matrix

|             | Vicsek (η=0.5)      | D'Orsogna (milling) |
|-------------|----------------------|---------------------|
| P5 flocking | ✓ DEFINITIVE (φ=.99)| ✗ none (φ=.008)     |
| P6 milling  | ✗ none (|L|=.017)   | ✓ DEFINITIVE (|L|=1)|

## Files to copy into repo

```
# New files (create these)
epc/models/vicsek.py
epc/models/dorsogna_spp.py
epc/metrics/collective_motion.py
epc/detectors/p5_flocking.py
epc/detectors/p6_milling.py
tests/test_vicsek_validation.py
tests/test_p5_detector_validation.py

# Updated files (replace these)
CLAUDE.md
PROJECT_STATUS.md
REPLICATION_NOTES.md

# Ensure these exist (may already)
epc/__init__.py
epc/models/__init__.py
epc/metrics/__init__.py
epc/detectors/__init__.py
```

## Suggested commit message

```
Sprint 3: Vicsek + D'Orsogna models, P5/P6 detectors, collective motion metrics

New models:
- Standard Vicsek Model (1995): phase transition replicated, 6/6 claims match
- D'Orsogna SPP with Morse potential (2006): milling at published params

New detectors:
- P5 flocking: 3-tier, shuffle null, P6/P7 exclusions → definitive on Vicsek
- P6 milling: 3-tier, ring density, P5 exclusion → definitive on D'Orsogna
- Perfect cross-discrimination: no false positives in 2×2 transfer matrix

New metrics: polarization, group speed ratio, angular momentum,
heading autocorrelation, nematic order (5 metrics in collective_motion.py)

Design decisions: nematic order for antiparallel detection (#21),
random-uniform null instead of permutation (#20), measured T_cross
for open-space models (#22)
```
