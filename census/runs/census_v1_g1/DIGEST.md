# Emergence Census digest — 2026-10-24 05:53

**289,798 / 289,798** programs done (100.0%), up to 22 bits. Mean 11.6 s per program. Flagged: 141,585. Errors: 0. Fingerprint errors: 0 (must be 0).

## Tier 4 — vetted findings
_none (vetting has not run)_

## Tier 3 — novelty candidates
_none (no class has passed the literature check as literature-silent)_

## Tier 2 — catalog expansions
_none yet_

## Tier 1 — unnamed behaviour classes awaiting a literature check (70)
- [9 bits, 9051 programs] `C kC=1 | C.KURA(K=0.1)`
- [9 bits, 1219 programs] `N kN=1 | N.KURA(K=0.1)`
- [11 bits, 2643 programs] `N kN=2 | N.GAME(g=(1.9, 0.0))`
- [11 bits, 2160 programs] `C kC=2 | C.COPY(p=1.0)`
- [11 bits, 494 programs] `C kC=4 | C.COPY(p=1.0)`
- [12 bits, 429 programs] `C kC=2 | C.GAME(g=(1.5, -0.5))`
- [12 bits, 64 programs] `C kC=2 | C.GAME(g=(1.3, 0.0))`
- [13 bits, 6631 programs] `C kC=2 | C.COPY(p=0.003)`
- [13 bits, 493 programs] `N kN=3 | N.CYCLE(th=3,p=0.1)`
- [13 bits, 486 programs] `N kN=3 | N.COPY(p=0.05)`
- [13 bits, 383 programs] `C kC=4 | C.CYCLE(th=1,p=1.0)`
- [13 bits, 344 programs] `C kC=3 | C.CYCLE(th=1,p=1.0)`
- [13 bits, 198 programs] `N kN=2 | N.REWIRE(p=0.01)`
- [13 bits, 195 programs] `N kN=3 | N.REWIRE(p=0.003)`
- [13 bits, 139 programs] `N kN=4 | N.CYCLE(th=1,p=1.0)`

## Review — known behaviour, unusual look (31126)

## Behaviour map so far

| behaviour | catalog | bits | simplest program | programs | routes |
|---|---|---|---|---|---|
| network synchronization | P9 | 9 | `N kN=1 | N.KURA(K=1.0)` | 4002 | 10 |
| lattice phase locking | P9 | 9 | `C kC=1 | C.KURA(K=1.0)` | 9356 | 24 |
| network consensus | P18 | 10 | `N kN=3 | N.COPY(p=1.0)` | 6198 | 32 |
| segregation with vacancies | P1 | 11 | `C kC=3 | C.SCHELL(th=0.5)` | 5608 | 18 |
| domain coarsening | P18 | 11 | `C kC=2 | C.COPY(p=0.1)` | 23542 | 68 |
| network fragmentation | P34 | 11 | `N kN=2 | N.REWIRE(p=1.0)` | 2751 | 9 |
| cyclic waves and spirals | — | 13 | `C kC=3 | C.CYCLE(th=3,p=1.0)` | 1309 | 12 |
| flocking | P5 | 14 | `A kA=1 v=1.0 eta=0.3 | A.ALIGN(r=1.0)` | 19080 | 21 |
| power-law avalanches | P14 | 15 | `C kC=4 | C.SAND(p=0.001)` | 151 | 8 |
| spatial game chaos | P27 | 15 | `CF kC=2 fields=1 D=[0.2] | C.GAME(g=(1.9, 0.0))` | 688 | 17 |
