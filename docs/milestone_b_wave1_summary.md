# Milestone B — Wave 1 Summary

**Date:** 2026-06-08
**Sprints:** 65–70
**Scope:** First three new-pattern implementations of Milestone B (self-propelled-particle family): P7 lane formation, P17 collective sensing, P19 emergent leadership.

---

## Outcome Table

| | P7 Lane Formation | P17 Collective Sensing | P19 Emergent Leadership |
|---|---|---|---|
| **Model file** | `epc/models/lane_formation.py` | `epc/models/collective_sensing.py` | `epc/models/informed_minority.py` |
| **Detector file** | `epc/detectors/p7_lane_formation.py` | `epc/detectors/p17_collective_sensing.py` | `epc/detectors/p19_emergent_leadership.py` |
| **Reference** | Helbing & Molnár (1995) | Berdahl et al. (2013) | Couzin et al. (2005) |
| **dim1 result** | PASS — φ_final=0.921±0.052 ∈ [0.5,1.0]; φ gain=0.423 ≥ 0.2; throughput fraction=0.998 ≥ 0.4 | PASS — CI slope vs log(N)=0.133 (>0.02); Spearman ρ=0.90 (p=0.037); CI(N=1)=−0.167 (chance); CI(N=50)=0.396 (>0.15) | PASS — accuracy-vs-ρ Spearman ρ=1.0; acc(ρ=0)=0.125 (chance); acc(ρ=0.05)=1.0 (>0.2); acc(ρ=0.5)=1.0 (>0.8) |
| **dim1 tolerance** | All 4 checks PASS | All 4 checks PASS | All 5 checks PASS |
| **dim1 published anchor** | Helbing & Molnár (1995) — lane order parameter in bidirectional pedestrian flow | Berdahl et al. (2013) Fig. 1 — CI rises with group size N (√N scaling) | Couzin et al. (2005) Fig. 2a — accuracy vs informed fraction ρ |
| **dim4 panel TNR** | 0.955 (PASS-with-weakness) | 1.000 (PASS) | 0.960 (PASS) |
| **dim4 Cohen's d** | 6.932 | 11.117 | 5.418 |
| **AT-DEPTH?** | Yes (Sprint 66) | Yes (Sprint 68) | Yes (Sprint 70) |

---

## Content Prerequisites Added During Panels

| Pattern | Prerequisite | Literature Grounding |
|---------|-------------|---------------------|
| P7 | Counterflow requires ≥2 populations with minority ≥10% | Helbing & Molnár (1995): lane formation is defined as spontaneous segregation of **bidirectional** pedestrian streams |
| P17 | (1) `field_samples` present in history; (2) individual SNR ≤ 3.0; (3) social cohesion ratio ≤ 0.20 | Berdahl et al. (2013): collective sensing requires noisy individual sensing (SNR < 3) + social coupling; group advantage is the mechanism |
| P19 | Early-window informed→naive leadership gap (convergence phase 10–40% of trajectory) | Couzin et al. (2005): informed minority must **lead** the alignment process — persistent directional preference steers the group |

---

## Open Carry-Forwards from Wave 1

| ID | Pattern | Description | Priority |
|----|---------|-------------|----------|
| C-p7-time-shuffled-fp | P7 | `time_shuffled` FP at screening — each frame independently preserves lane structure (high φ_lane) | Low |
| C-p19-bias-zero-chance-alignment | P19 | 1/5 bias_zero Class C regimes (seed=410) reaches confirmation by chance alignment toward θ_pref=0.0 | Low (4% FP rate, within TNR ≥ 0.95 tolerance) |
| C-p19-abrupt-saturation | P19 | accuracy saturates to ~1.0 at ρ=0.025, more abruptly than Couzin (2005) Fig. 2a; root cause: η=0.1 is low | Not a validation failure |
| C-p19-te-vs-pull | P19 | TE (KSG) on mean heading collapses in converged regime; label-shuffle pull used instead | Documented architectural decision |

---

## Confirmation: Each dim1 Reproduction Hits a Real Published Anchor

- **P7:** Helbing & Molnár (1995) — the canonical social-force pedestrian dynamics paper. Lane order parameter φ and throughput gain reproduced against the published bidirectional-flow lane formation phenomenology. Explicit tolerance thresholds defined and met.
- **P17:** Berdahl et al. (2013) Fig. 1 — CI (chemotactic index) scaling with group size N. The √N improvement in collective accuracy is the paper's central result. Four quantitative tolerance checks all pass.
- **P19:** Couzin et al. (2005) Fig. 2a — accuracy vs informed fraction ρ. The defining result of the paper (small informed minority guides the whole group). Five quantitative tolerance checks all pass. Note: saturation is more abrupt than published due to low noise parameter; documented, not a validation gap.

---

## Aggregate Counts (Post-Wave 1)

- **Implemented patterns:** 22
- **AT-DEPTH:** 21 / 22
- **Remaining gap:** P12 (dim1 — documented finite-size measurement limitation)
- **Wave 1 result:** 3/3 new patterns reached AT-DEPTH
