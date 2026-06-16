# Battery Confusion Matrix — Instrument Self-Recognition (T2c-flavored)

**Date:** 2026-06-16. **HEAD:** 592ad2a. Each of the 31 canonical positives
(P1–P30, P32; single seed) was run through the full battery via
`analysis/battery_profile.profile_observation`. Top-1 = the pattern whose
detector fired at the highest tier (det=True first), calibrated confidence breaking ties.

**Diagonal self-recognition: 30/31.**

| truth | top-1 | tier | verdict | emergence |
|---|---|---|---|---|
| P1 | P1 | definitive | MATCH | 0.7769 |
| P2 | P2 | none | EMERGENT-UNCLASSIFIED | 0.5368 |
| P3 | P3 | definitive | MATCH | 0.7551 |
| P4 | P4 | definitive | MATCH | 0.0105 |
| P5 | P5 | definitive | MATCH | 0.9811 |
| P6 | P6 | definitive | MATCH | 0.0003 |
| P7 | P7 | confirmation | MATCH | 0.0 |
| P8 | P8 | definitive | MATCH | 0.0 |
| P9 | P9 | definitive | MATCH | 0.7513 |
| P10 | P10 | definitive | MATCH | 0.728 |
| P11 | P11 | definitive | MATCH | 0.7602 |
| P12 | P12 | confirmation | MATCH | 0.963 |
| P13 ⚠ | P18 | screening | MATCH | 0.7624 |
| P14 | P14 | confirmation | MATCH | 0.0 |
| P15 | P15 | definitive | MATCH | 0.7837 |
| P16 | P16 | definitive | MATCH | 0.0892 |
| P17 | P17 | definitive | MATCH | 0.1034 |
| P18 | P18 | definitive | MATCH | 0.7777 |
| P19 | P19 | definitive | MATCH | 0.8502 |
| P20 | P20 | definitive | MATCH | 0.0894 |
| P21 | P21 | definitive | MATCH | 0.0985 |
| P22 | P22 | screening | MATCH | 0.0598 |
| P23 | P23 | definitive | MATCH | 0.0894 |
| P24 | P24 | definitive | MATCH | 0.0894 |
| P25 | P25 | definitive | MATCH | 0.7504 |
| P26 | P26 | definitive | MATCH | 0.0894 |
| P27 | P27 | screening | MATCH | 0.0604 |
| P28 | P28 | definitive | MATCH | 0.1142 |
| P29 | P29 | screening | MATCH | 0.0 |
| P30 | P30 | definitive | MATCH | 0.1892 |
| P32 | P32 | definitive | MATCH | 0.0903 |

## Residuals (honest)
- **P13 → P18 (miss):** on the excitable-waves grid the consensus detector (P18)
  fires at screening and edges P13, which reached only screening on this single
  seed. Both are grid-coarsening patterns. Candidate fix: P18-excludes-P13, or
  strengthen the P13 single-seed positive. One cell.
- **P12 → P18 (off-diagonal confirmation+):** RPS and consensus share grid
  coarsening; P12 still ranks first (its own confirmation), P18 is the only
  other pattern reaching confirmation+ on any positive. One cell.
- Single-seed fragility (P22 transient, P27/P29) yields screening-tier or
  tier-none self-recognition on seed 0; multi-seed would lift these.