# Ring 2 — Comprehensive hunt (frozen 21-lens battery)

The definitive novelty hunt after the comprehensive lens build-out ("breadth + depth first, hunt
once"). Battery FROZEN at 21 admitted coordinate lenses (the original 6 + the 15 build-out lenses),
integrity-confirmed (particle_life seed 3000 N=120 → 0 leads with the full battery, identical to the
6-lens run: same configs → same leads; coordinate lenses cannot move the tripwire).

## Sweep
440 configs across 4 generative families with the full descriptor (full 21-lens fingerprints saved):
particle_life 120 (seed 3000), reaction_diffusion 120 (4000), kuramoto 100 (6000), lenia 100 (7000).

## Result 1 — tripwire frontier: 0 genuine-novelty leads / 440
| family | leads |
|---|---|
| particle_life | 0 / 120 |
| reaction_diffusion | 0 / 120 |
| kuramoto | 0 / 100 |
| lenia | 0 / 100 |
The em/MPR-C complex+unclassified frontier is empty. The prior broadened-pass conclusion holds
under the comprehensive battery (the tripwire is unchanged; the richer lenses do not move leads).

## Result 2 — lens-space outlier scan (the discovery mechanism the battery enables)
`lens_space_outliers.py`: per-family RANK-transform of every numeric lens coordinate (robust to
heavy tails / scale — no single coord can dominate) → kNN outlier score. Finds configs structurally
unusual across ANY of the 21 axes, which the em/C tripwire alone would miss.

Findings: the top outliers are MILD (z ≈ 2.0–3.5; expected from a homogeneous population of ~100 —
no config stands dramatically apart) AND every one is `tripped=False` with a RECOGNISED `em_kind`,
i.e. CLASSIFIED by the named channels. They are known emergent states at the edge of their family's
parameter range, not unclassified novelty:
- kuramoto #15 — em_kind local-phase-order, circulation rank 0.97 → rotating / traveling phase wave
  (the known TWISTED/WAVE family, already closed via round-6 net-widening).
- particle_life #95 — em_kind orientation, polarization rank 1.0 → a polar flock.
- reaction_diffusion #2 — em_kind field, LOW on coarsen/loops/symmetry → near-homogeneous/transient field.
- lenia #41 — em_kind field → a structured creature field.

(A first MAD-z pass was discarded: it was dominated by scale artifacts in a few heavy-tailed coords
— sk_peak z≈2376, em_score z≈−641 — which is a units issue, not novelty. The rank transform fixes it.)

## Conclusion
**0 genuine novelty-to-science across both the tripwire frontier and the lens-space outlier frontier.**
The prior result held under the raised ceiling. The comprehensive 21-lens battery did not surface any
unclassified structural novelty; its outliers are all recognised states. The honest deliverable is
the AUDITED FRONTIER + the comprehensive interpretable instrument itself (21 validated lenses, the
differentiation from ASAL's single opaque descriptor). No catalog pattern is added (the catalog
remains 37). Next, if pursued: the documented second-tier missing axes (traveling-wave dispersion,
glassy χ4, self-replication, …) and/or new generative families beyond the four swept.
