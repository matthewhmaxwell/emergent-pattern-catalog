# EPC Emergence Census — post-run registered analyses

Generated 2026-10-25 05:56 on the frozen census code. These are the analyses promised in PREREG_census_v1.md (S3, S4.5, S5, S10) that were not run before the first results report.

## S5 False-flag audit (placebo on fresh census programs)

- 3000 programs (0 errors), 5972 interaction-free views compared with an independent interaction-free run: **0 false flags** -> rate <= 0.050% (95%, rule of three)

## S10 Minimal programs re-checked on 10 seeds; S5 rule knock-outs

Criterion: PASS = flagged under the census rule with 10 seeds (>=2 seeds emergent AND interaction-driven by the paired knock-out) AND the behaviour named/shown by the namer, AND emergent in >= half of the 10 seeds

| behaviour | grammar | census minimum (bits) | confirmed minimum (bits) | confirmed program | knock-outs |
|---|---|---|---|---|---|
| cyclic waves and spirals | g1 | 13 | 13 | `C kC=3 | C.CYCLE(th=3,p=1.0)` | sole rule -> interaction-free -> behaviour absent by construction (necessary) |
| cyclic waves and spirals | g2 | 18 | 18 | `C kC=3 | C[any] NEXT_GE(th=1) -> COPY_RAND() p=1.0` |  |
| domain coarsening | g1 | 11 | 11 | `C kC=2 | C.COPY(p=0.1)` | sole rule -> interaction-free -> behaviour absent by construction (necessary) |
| domain coarsening | g2 | 14 | 14 | `C kC=2 | C[any] ALWAYS() -> COPY_RAND() p=0.1` |  |
| flocking | g1 | 14 | 14 | `A kA=1 v=1.0 eta=0.3 | A.ALIGN(r=1.0)` | sole rule -> interaction-free -> behaviour absent by construction (necessary) |
| flocking | g2 | 17 | 17 | `A kA=1 v=1.0 eta=0.3 | A[any] ALWAYS() -> HEAD_MEAN(b=any,r=1.0) p=1.0` |  |
| lattice phase locking | g1 | 9 | 9 | `C kC=1 | C.KURA(K=1.0)` | sole rule -> interaction-free -> behaviour absent by construction (necessary) |
| lattice phase locking | g2 | 11 | 11 | `C kC=1 | C[any] ALWAYS() -> PHASE_COUPLE(K=1.0) p=1.0` |  |
| network consensus | g1 | 10 | 10 | `N kN=3 | N.COPY(p=1.0)` | sole rule -> interaction-free -> behaviour absent by construction (necessary) |
| network consensus | g2 | 14 | 14 | `N kN=2 | N[any] ALWAYS() -> COPY_RAND() p=1.0` |  |
| network fragmentation | g1 | 11 | 11 | `N kN=2 | N.REWIRE(p=1.0)` | sole rule -> interaction-free -> behaviour absent by construction (necessary) |
| network fragmentation | g2 | 14 | 14 | `N kN=2 | N[any] ALWAYS() -> REWIRE_SAME() p=1.0` |  |
| network synchronization | g1 | 9 | 9 | `N kN=1 | N.KURA(K=1.0)` | sole rule -> interaction-free -> behaviour absent by construction (necessary) |
| network synchronization | g2 | 11 | 11 | `N kN=1 | N[any] ALWAYS() -> PHASE_COUPLE(K=1.0) p=1.0` |  |
| power-law avalanches | g1 | 15 | 15 | `C kC=4 | C.SAND(p=0.001)` | sole rule -> interaction-free -> behaviour absent by construction (necessary) |
| power-law avalanches | g2 | 18 | 18 | `C kC=4 | C[any] ALWAYS() -> ADD_GRAIN() p=0.001` |  |
| segregation with vacancies | g1 | 11 | 11 | `C kC=3 | C.SCHELL(th=0.5)` | sole rule -> interaction-free -> behaviour absent by construction (necessary) |
| segregation with vacancies | g2 | 18 | 18 | `C kC=3 | C[any] NEXT_GE(th=3) -> MOVE_EMPTY() p=0.1` |  |
| spatial game chaos | g1 | 15 | 15 | `CF kC=2 fields=1 D=[0.2] | C.GAME(g=(1.9, 0.0))` | sole rule -> interaction-free -> behaviour absent by construction (necessary) |
| spatial game chaos | g2 | 16 | 16 | `C kC=2 | C[any] ALWAYS() -> IMITATE_BEST(g=(1.9, 0.0)) p=1.0` |  |

## S3 Avalanche confirmation (5,000 steps, confirmation tier)

- [g1] 15b `C kC=4 | C.SAND(p=0.001)` -> tiers ['confirmation', 'confirmation', 'confirmation'] confirmed=True
- [g1] 18b `CF kC=4 fields=1 D=[0.2] | C.SAND(p=0.001)` -> tiers ['confirmation', 'confirmation', 'none'] confirmed=True
- [g1] 18b `CF kC=4 fields=1 D=[0.05] | C.SAND(p=0.001)` -> tiers ['confirmation', 'none', 'confirmation'] confirmed=True
- [g2] 18b `C kC=4 | C[any] ALWAYS() -> ADD_GRAIN() p=0.001` -> tiers ['confirmation', 'confirmation', 'confirmation'] confirmed=True
- [g2] 18b `C kC=4 | C[3] ALWAYS() -> ADD_GRAIN() p=0.003` -> tiers ['none', 'none', 'none'] confirmed=False
- [g2] 19b `C kC=4 | C[any] SAME_LT(th=0.2) -> ADD_GRAIN() p=0.1` -> tiers ['confirmation', 'confirmation', 'confirmation'] confirmed=True

## S5 Robustness of each minimal program

fraction of all census programs with the same rule structure as the confirmed minimal program (substrates, type counts, field count, rule templates; all rule/agent/diffusion parameters free) whose run is flagged and shows the behaviour

| behaviour | grammar | minimal program | programs with that structure | showing the behaviour | fraction |
|---|---|---|---|---|---|
| cyclic waves and spirals | g1 | `C kC=3 | C.CYCLE(th=3,p=1.0)` | 112 | 7 | 0.062 |
| cyclic waves and spirals | g2 | `C kC=3 | C[any] NEXT_GE(th=1) -> COPY_RAND() p=1.0` | 48 | 6 | 0.125 |
| domain coarsening | g1 | `C kC=2 | C.COPY(p=0.1)` | 14 | 8 | 0.571 |
| domain coarsening | g2 | `C kC=2 | C[any] ALWAYS() -> COPY_RAND() p=0.1` | 14 | 8 | 0.571 |
| flocking | g1 | `A kA=1 v=1.0 eta=0.3 | A.ALIGN(r=1.0)` | 216 | 186 | 0.861 |
| flocking | g2 | `A kA=1 v=1.0 eta=0.3 | A[any] ALWAYS() -> HEAD_MEAN(b=any,r=1.0) p=1.0` | 1488 | 756 | 0.508 |
| lattice phase locking | g1 | `C kC=1 | C.KURA(K=1.0)` | 14 | 4 | 0.286 |
| lattice phase locking | g2 | `C kC=1 | C[any] ALWAYS() -> PHASE_COUPLE(K=1.0) p=1.0` | 196 | 30 | 0.153 |
| network consensus | g1 | `N kN=3 | N.COPY(p=1.0)` | 14 | 5 | 0.357 |
| network consensus | g2 | `N kN=2 | N[any] ALWAYS() -> COPY_RAND() p=1.0` | 14 | 4 | 0.286 |
| network fragmentation | g1 | `N kN=2 | N.REWIRE(p=1.0)` | 14 | 8 | 0.571 |
| network fragmentation | g2 | `N kN=2 | N[any] ALWAYS() -> REWIRE_SAME() p=1.0` | 14 | 8 | 0.571 |
| network synchronization | g1 | `N kN=1 | N.KURA(K=1.0)` | 14 | 4 | 0.286 |
| network synchronization | g2 | `N kN=1 | N[any] ALWAYS() -> PHASE_COUPLE(K=1.0) p=1.0` | 196 | 24 | 0.122 |
| power-law avalanches | g1 | `C kC=4 | C.SAND(p=0.001)` | 14 | 1 | 0.071 |
| power-law avalanches | g2 | `C kC=4 | C[any] ALWAYS() -> ADD_GRAIN() p=0.001` | 14 | 1 | 0.071 |
| segregation with vacancies | g1 | `C kC=3 | C.SCHELL(th=0.5)` | 6 | 6 | 1.0 |
| segregation with vacancies | g2 | `C kC=3 | C[any] NEXT_GE(th=3) -> MOVE_EMPTY() p=0.1` | 48 | 12 | 0.25 |
| spatial game chaos | g1 | `CF kC=2 fields=1 D=[0.2] | C.GAME(g=(1.9, 0.0))` | 36 | 4 | 0.111 |
| spatial game chaos | g2 | `C kC=2 | C[any] ALWAYS() -> IMITATE_BEST(g=(1.9, 0.0)) p=1.0` | 84 | 19 | 0.226 |

## S4.5 Known-behaviour-atypical-look views, reviewed like unnamed groups

- 37860 atypical-look views -> 155 fingerprint groups -> 73 distinct phenomena; 73 literature-checked
- verdicts: {'KNOWN': 69, 'KNOWN-VARIANT': 4}
- **literature-silent: 0**

