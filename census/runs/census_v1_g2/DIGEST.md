# Emergence Census digest — 2026-10-24 05:55

**98,082 / 98,082** programs done (100.0%), up to 21 bits. Mean 5.6 s per program. Flagged: 21,172. Errors: 0. Fingerprint errors: 0 (must be 0).

## Tier 4 — vetted findings
_none (vetting has not run)_

## Tier 3 — novelty candidates
_none (no class has passed the literature check as literature-silent)_

## Tier 2 — catalog expansions
_none yet_

## Tier 1 — unnamed behaviour classes awaiting a literature check (61)
- [11 bits, 2725 programs] `C kC=1 | C[any] ALWAYS() -> PHASE_COUPLE(K=0.1) p=0.1`
- [11 bits, 647 programs] `N kN=1 | N[any] ALWAYS() -> PHASE_COUPLE(K=0.1) p=1.0`
- [11 bits, 501 programs] `C kC=1 | C[any] ALWAYS() -> PHASE_COUPLE(K=0.1) p=1.0`
- [13 bits, 285 programs] `C kC=1 | C[any] ALWAYS() -> PHASE_COUPLE(K=0.03) p=1.0`
- [13 bits, 268 programs] `N kN=1 | N[any] ALWAYS() -> PHASE_COUPLE(K=1.0) p=0.01`
- [14 bits, 79 programs] `N kN=2 | N[any] ALWAYS() -> COPY_RAND() p=0.1`
- [14 bits, 46 programs] `C kC=2 | C[any] ALWAYS() -> COPY_RAND() p=1.0`
- [15 bits, 2000 programs] `N kN=3 | N[0] ALWAYS() -> COPY_RAND() p=0.1`
- [15 bits, 1195 programs] `C kC=3 | C[0] ALWAYS() -> COPY_RAND() p=1.0`
- [15 bits, 628 programs] `C kC=3 | C[1] ALWAYS() -> MOVE_EMPTY() p=1.0`
- [15 bits, 431 programs] `N kN=4 | N[0] ALWAYS() -> COPY_RAND() p=0.1`
- [15 bits, 353 programs] `C kC=2 | C[0] ALWAYS() -> COPY_RAND() p=0.1`
- [15 bits, 18 programs] `C kC=3 | C[any] ALWAYS() -> COPY_RAND() p=1.0`
- [15 bits, 12 programs] `N kN=4 | N[any] DISCORD() -> SET_NEXT() p=1.0`
- [15 bits, 4 programs] `C kC=3 | C[any] ALWAYS() -> COPY_MAJ() p=1.0`

## Review — known behaviour, unusual look (4024)

## Behaviour map so far

| behaviour | catalog | bits | simplest program | programs | routes |
|---|---|---|---|---|---|
| network synchronization | P9 | 11 | `N kN=1 | N[any] ALWAYS() -> PHASE_COUPLE(K=1.0) p=1.0` | 1205 | 6 |
| lattice phase locking | P9 | 11 | `C kC=1 | C[any] ALWAYS() -> PHASE_COUPLE(K=1.0) p=1.0` | 2236 | 10 |
| domain coarsening | P18 | 14 | `C kC=2 | C[any] ALWAYS() -> COPY_RAND() p=0.1` | 1885 | 46 |
| network fragmentation | P34 | 14 | `N kN=2 | N[any] ALWAYS() -> REWIRE_SAME() p=1.0` | 237 | 5 |
| network consensus | P18 | 14 | `N kN=2 | N[any] ALWAYS() -> COPY_RAND() p=1.0` | 1552 | 23 |
| spatial game chaos | P27 | 16 | `C kC=2 | C[any] ALWAYS() -> IMITATE_BEST(g=(1.9, 0.0)) p=1.0` | 103 | 6 |
| flocking | P5 | 17 | `A kA=1 v=1.0 eta=0.3 | A[any] ALWAYS() -> HEAD_MEAN(b=any,r=1.0) p=1.0` | 843 | 4 |
| cyclic waves and spirals | — | 18 | `C kC=3 | C[any] NEXT_GE(th=1) -> COPY_RAND() p=1.0` | 19 | 3 |
| power-law avalanches | P14 | 18 | `C kC=4 | C[any] ALWAYS() -> ADD_GRAIN() p=0.001` | 24 | 5 |
| segregation with vacancies | P1 | 18 | `C kC=3 | C[any] NEXT_GE(th=3) -> MOVE_EMPTY() p=0.1` | 108 | 3 |
