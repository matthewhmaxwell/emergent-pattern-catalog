# Census triage — 2026-10-24 05:53

289,798 programs; 141,585 with a flagged view; errors 0; fingerprint errors 0 (must be 0).
Flagged views: 145,471 = 39,114 named + 31,126 known-behaviour-atypical-look (review) + 75,231 unnamed in 70 classes (literature check).

## Behaviour map (shortest program per behaviour = the emergence-complexity table)

| behaviour | catalog | bits | shortest program | programs | routes |
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

## Routes (every distinct mechanism that produces each behaviour)

### network synchronization — 10 routes; ingredient common to all: N.KURA

| layers | rule templates | programs | bits | shortest example |
|---|---|---|---|---|
| N | N.KURA | 681 | 9 | `N kN=1 | N.KURA(K=1.0)` |
| N | N.COPY, N.KURA | 2358 | 13 | `N kN=1 | N.COPY(p=1.0) ; N.KURA(K=1.0)` |
| N | N.KURA, N.MAJ | 240 | 17 | `N kN=2 | N.MAJ(p=1.0) ; N.KURA(K=1.0)` |
| N | N.KURA, N.REWIRE | 137 | 17 | `N kN=2 | N.KURA(K=1.0) ; N.REWIRE(p=1.0)` |
| N | N.GAME, N.KURA | 48 | 17 | `N kN=2 | N.GAME(g=(1.9, 0.0)) ; N.KURA(K=1.0)` |
| N | N.CYCLE, N.KURA | 384 | 19 | `N kN=2 | N.CYCLE(th=1,p=1.0) ; N.KURA(K=1.0)` |
| N | N.FLIP, N.KURA | 60 | 19 | `N kN=2 | N.FLIP(a=0,b=1,p=1.0) ; N.KURA(K=1.0)` |
| N | N.CONTAG, N.KURA | 24 | 21 | `N kN=2 | N.CONTAG(a=0,b=1,th=1,p=1.0) ; N.KURA(K=1.0)` |
| N | N.COPY, N.KURA, N.MAJ | 48 | 22 | `N kN=3 | N.COPY(p=1.0) ; N.MAJ(p=1.0) ; N.KURA(K=1.0)` |
| N | N.COPY, N.KURA, N.REWIRE | 22 | 22 | `N kN=3 | N.COPY(p=1.0) ; N.KURA(K=1.0) ; N.REWIRE(p=0.1)` |

### lattice phase locking — 24 routes; ingredient common to all: C.KURA

| layers | rule templates | programs | bits | shortest example |
|---|---|---|---|---|
| C | C.KURA | 678 | 9 | `C kC=1 | C.KURA(K=1.0)` |
| C | C.COPY, C.KURA | 2232 | 13 | `C kC=1 | C.COPY(p=1.0) ; C.KURA(K=1.0)` |
| CF | C.KURA | 1050 | 13 | `CF kC=1 fields=1 D=[0.2] | C.KURA(K=1.0)` |
| CA | C.KURA | 1920 | 17 | `CA kC=1 kA=1 v=1.0 eta=0.3 | C.KURA(K=1.0)` |
| C | C.KURA, C.MAJ | 208 | 17 | `C kC=2 | C.MAJ(p=1.0) ; C.KURA(K=1.0)` |
| C | C.KURA, C.SWAP | 208 | 17 | `C kC=2 | C.SWAP(p=1.0) ; C.KURA(K=1.0)` |
| C | C.KURA, C.SCHELL | 96 | 17 | `C kC=3 | C.SCHELL(th=0.5) ; C.KURA(K=1.0)` |
| C | C.GAME, C.KURA | 48 | 17 | `C kC=2 | C.GAME(g=(1.9, 0.0)) ; C.KURA(K=1.0)` |
| CF | C.COPY, C.KURA | 520 | 18 | `CF kC=1 fields=1 D=[0.2] | C.COPY(p=1.0) ; C.KURA(K=1.0)` |
| CF | C.EMIT, C.KURA | 288 | 18 | `CF kC=1 fields=1 D=[0.2] | C.KURA(K=1.0) ; C.EMIT(b=0,f=0,p=1.0)` |
| C | C.CYCLE, C.KURA | 312 | 19 | `C kC=2 | C.CYCLE(th=1,p=1.0) ; C.KURA(K=1.0)` |
| CF | C.KURA, F.DECAY | 200 | 19 | `CF kC=1 fields=1 D=[0.2] | C.KURA(K=1.0) ; F.DECAY(f=0,d=1.0)` |
| CF | C.KURA, F.FEED | 200 | 19 | `CF kC=1 fields=1 D=[0.2] | C.KURA(K=1.0) ; F.FEED(f=0,F=1.0)` |
| C | C.FLIP, C.KURA | 60 | 19 | `C kC=2 | C.FLIP(a=0,b=1,p=1.0) ; C.KURA(K=1.0)` |
| C | C.KURA, C.SAND | 40 | 19 | `C kC=4 | C.SAND(p=1.0) ; C.KURA(K=1.0)` |
| … | 9 more routes in triage.json | | | |

### network consensus — 32 routes; ingredient common to all: none

| layers | rule templates | programs | bits | shortest example |
|---|---|---|---|---|
| N | N.COPY | 424 | 10 | `N kN=3 | N.COPY(p=1.0)` |
| N | N.MAJ | 369 | 11 | `N kN=2 | N.MAJ(p=1.0)` |
| N | N.GAME | 32 | 11 | `N kN=2 | N.GAME(g=(1.5, 0.5))` |
| N | N.CYCLE | 20 | 13 | `N kN=2 | N.CYCLE(th=3,p=0.1)` |
| N | N.CONTAG | 32 | 15 | `N kN=2 | N.CONTAG(a=0,b=1,th=1,p=1.0)` |
| N | N.COPY, N.MAJ | 929 | 16 | `N kN=3 | N.COPY(p=1.0) ; N.MAJ(p=1.0)` |
| N | N.COPY, N.KURA | 430 | 16 | `N kN=3 | N.COPY(p=1.0) ; N.KURA(K=1.0)` |
| N | N.COPY, N.REWIRE | 275 | 16 | `N kN=3 | N.COPY(p=1.0) ; N.REWIRE(p=0.1)` |
| N | N.KURA, N.MAJ | 483 | 17 | `N kN=2 | N.MAJ(p=1.0) ; N.KURA(K=1.0)` |
| N | N.GAME, N.REWIRE | 154 | 17 | `N kN=2 | N.GAME(g=(1.9, 0.0)) ; N.REWIRE(p=1.0)` |
| N | N.MAJ, N.REWIRE | 149 | 17 | `N kN=2 | N.MAJ(p=1.0) ; N.REWIRE(p=0.1)` |
| N | N.GAME, N.MAJ | 144 | 17 | `N kN=2 | N.MAJ(p=1.0) ; N.GAME(g=(1.9, 0.0))` |
| N | N.COPY, N.GAME | 135 | 17 | `N kN=2 | N.COPY(p=1.0) ; N.GAME(g=(1.9, 0.0))` |
| N | N.GAME, N.KURA | 132 | 17 | `N kN=2 | N.GAME(g=(1.5, 0.5)) ; N.KURA(K=1.0)` |
| N | N.COPY, N.CYCLE | 906 | 18 | `N kN=3 | N.COPY(p=1.0) ; N.CYCLE(th=3,p=1.0)` |
| … | 17 more routes in triage.json | | | |

### segregation with vacancies — 18 routes; ingredient common to all: C.SCHELL

| layers | rule templates | programs | bits | shortest example |
|---|---|---|---|---|
| C | C.SCHELL | 81 | 11 | `C kC=3 | C.SCHELL(th=0.5)` |
| CF | C.SCHELL | 440 | 15 | `CF kC=3 fields=1 D=[0.2] | C.SCHELL(th=0.5)` |
| C | C.MAJ, C.SCHELL | 326 | 17 | `C kC=3 | C.MAJ(p=1.0) ; C.SCHELL(th=0.3)` |
| C | C.COPY, C.SCHELL | 320 | 17 | `C kC=3 | C.COPY(p=1.0) ; C.SCHELL(th=0.5)` |
| C | C.SCHELL, C.SWAP | 304 | 17 | `C kC=3 | C.SWAP(p=1.0) ; C.SCHELL(th=0.5)` |
| C | C.KURA, C.SCHELL | 281 | 17 | `C kC=3 | C.SCHELL(th=0.5) ; C.KURA(K=1.0)` |
| CA | C.SCHELL | 1458 | 18 | `CA kC=3 kA=1 v=1.0 eta=0.3 | C.SCHELL(th=0.5)` |
| C | C.CYCLE, C.SCHELL | 690 | 19 | `C kC=3 | C.CYCLE(th=1,p=1.0) ; C.SCHELL(th=0.5)` |
| C | C.FLIP, C.SCHELL | 136 | 20 | `C kC=3 | C.FLIP(a=0,b=1,p=0.1) ; C.SCHELL(th=0.5)` |
| CAF | C.SCHELL | 1096 | 21 | `CAF kC=3 kA=1 v=1.0 eta=0.3 fields=1 D=[0.2] | C.SCHELL(th=0.5)` |
| CF | C.COPY, C.SCHELL | 160 | 21 | `CF kC=3 fields=1 D=[0.2] | C.COPY(p=1.0) ; C.SCHELL(th=0.5)` |
| CF | C.MAJ, C.SCHELL | 137 | 21 | `CF kC=3 fields=1 D=[0.2] | C.MAJ(p=1.0) ; C.SCHELL(th=0.3)` |
| C | C.SAND, C.SCHELL | 23 | 21 | `C kC=4 | C.SCHELL(th=0.5) ; C.SAND(p=0.001)` |
| CF | C.SCHELL, C.SWAP | 32 | 22 | `CF kC=3 fields=1 D=[0.2] | C.SWAP(p=1.0) ; C.SCHELL(th=0.5)` |
| CF | C.KURA, C.SCHELL | 32 | 22 | `CF kC=3 fields=1 D=[0.2] | C.SCHELL(th=0.5) ; C.KURA(K=1.0)` |
| … | 3 more routes in triage.json | | | |

### domain coarsening — 68 routes; ingredient common to all: none

| layers | rule templates | programs | bits | shortest example |
|---|---|---|---|---|
| C | C.MAJ | 422 | 11 | `C kC=2 | C.MAJ(p=0.1)` |
| C | C.COPY | 340 | 11 | `C kC=2 | C.COPY(p=0.1)` |
| C | C.SCHELL | 13 | 12 | `C kC=3 | C.SCHELL(th=0.7)` |
| C | C.CYCLE | 137 | 13 | `C kC=2 | C.CYCLE(th=3,p=0.1)` |
| CF | C.MAJ | 1819 | 14 | `CF kC=2 fields=1 D=[0.2] | C.MAJ(p=1.0)` |
| CF | C.COPY | 1234 | 14 | `CF kC=2 fields=1 D=[0.2] | C.COPY(p=0.1)` |
| CF | C.CYCLE | 1101 | 16 | `CF kC=2 fields=1 D=[0.2] | C.CYCLE(th=3,p=0.1)` |
| CF | C.SCHELL | 84 | 16 | `CF kC=3 fields=1 D=[0.2] | C.SCHELL(th=0.7)` |
| CA | C.MAJ | 4428 | 17 | `CA kC=2 kA=1 v=1.0 eta=0.3 | C.MAJ(p=1.0)` |
| CA | C.COPY | 3168 | 17 | `CA kC=2 kA=1 v=1.0 eta=0.3 | C.COPY(p=0.1)` |
| C | C.COPY, C.MAJ | 752 | 17 | `C kC=2 | C.COPY(p=1.0) ; C.MAJ(p=0.1)` |
| C | C.MAJ, C.SWAP | 628 | 17 | `C kC=2 | C.MAJ(p=1.0) ; C.SWAP(p=1.0)` |
| C | C.KURA, C.MAJ | 624 | 17 | `C kC=2 | C.MAJ(p=1.0) ; C.KURA(K=0.1)` |
| C | C.COPY, C.KURA | 437 | 17 | `C kC=2 | C.COPY(p=0.1) ; C.KURA(K=1.0)` |
| C | C.MAJ, C.SCHELL | 304 | 17 | `C kC=3 | C.MAJ(p=1.0) ; C.SCHELL(th=0.3)` |
| … | 53 more routes in triage.json | | | |

### network fragmentation — 9 routes; ingredient common to all: N.REWIRE

| layers | rule templates | programs | bits | shortest example |
|---|---|---|---|---|
| N | N.REWIRE | 376 | 11 | `N kN=2 | N.REWIRE(p=1.0)` |
| N | N.COPY, N.REWIRE | 521 | 16 | `N kN=3 | N.COPY(p=1.0) ; N.REWIRE(p=1.0)` |
| N | N.KURA, N.REWIRE | 488 | 17 | `N kN=2 | N.KURA(K=1.0) ; N.REWIRE(p=1.0)` |
| N | N.MAJ, N.REWIRE | 475 | 17 | `N kN=2 | N.MAJ(p=1.0) ; N.REWIRE(p=1.0)` |
| N | N.CYCLE, N.REWIRE | 662 | 19 | `N kN=2 | N.CYCLE(th=1,p=0.1) ; N.REWIRE(p=1.0)` |
| N | N.FLIP, N.REWIRE | 47 | 20 | `N kN=3 | N.FLIP(a=0,b=1,p=1.0) ; N.REWIRE(p=1.0)` |
| N | N.CONTAG, N.REWIRE | 42 | 21 | `N kN=2 | N.CONTAG(a=0,b=1,th=1,p=0.1) ; N.REWIRE(p=1.0)` |
| N | N.COPY, N.KURA, N.REWIRE | 72 | 22 | `N kN=3 | N.COPY(p=1.0) ; N.KURA(K=1.0) ; N.REWIRE(p=1.0)` |
| N | N.COPY, N.MAJ, N.REWIRE | 68 | 22 | `N kN=3 | N.COPY(p=1.0) ; N.MAJ(p=1.0) ; N.REWIRE(p=1.0)` |

### cyclic waves and spirals — 12 routes; ingredient common to all: none

| layers | rule templates | programs | bits | shortest example |
|---|---|---|---|---|
| C | C.CYCLE | 56 | 13 | `C kC=3 | C.CYCLE(th=3,p=1.0)` |
| CF | C.CYCLE | 301 | 16 | `CF kC=3 fields=1 D=[0.2] | C.CYCLE(th=3,p=1.0)` |
| C | C.CYCLE, C.MAJ | 125 | 19 | `C kC=3 | C.MAJ(p=1.0) ; C.CYCLE(th=1,p=1.0)` |
| C | C.COPY, C.CYCLE | 115 | 19 | `C kC=3 | C.COPY(p=1.0) ; C.CYCLE(th=3,p=1.0)` |
| C | C.CYCLE, C.SWAP | 92 | 19 | `C kC=3 | C.CYCLE(th=3,p=1.0) ; C.SWAP(p=1.0)` |
| C | C.CYCLE, C.KURA | 71 | 19 | `C kC=3 | C.CYCLE(th=3,p=1.0) ; C.KURA(K=1.0)` |
| C | C.CYCLE, C.SCHELL | 25 | 19 | `C kC=3 | C.CYCLE(th=3,p=1.0) ; C.SCHELL(th=0.3)` |
| CA | C.CYCLE | 485 | 20 | `CA kC=3 kA=1 v=1.0 eta=0.3 | C.CYCLE(th=3,p=1.0)` |
| CF | C.CYCLE, C.MAJ | 16 | 22 | `CF kC=3 fields=1 D=[0.2] | C.MAJ(p=1.0) ; C.CYCLE(th=1,p=1.0)` |
| CF | C.COPY, C.CYCLE | 15 | 22 | `CF kC=3 fields=1 D=[0.2] | C.COPY(p=1.0) ; C.CYCLE(th=3,p=1.0)` |
| C | C.CYCLE, C.FLIP | 4 | 22 | `C kC=3 | C.CYCLE(th=1,p=0.1) ; C.FLIP(a=0,b=1,p=0.1)` |
| CA | A.WRITECELL, C.MAJ | 4 | 22 | `CA kC=3 kA=1 v=1.0 eta=0.3 | C.MAJ(p=0.1) ; A.WRITECELL()` |

### flocking — 21 routes; ingredient common to all: none

| layers | rule templates | programs | bits | shortest example |
|---|---|---|---|---|
| A | A.ALIGN | 4774 | 14 | `A kA=1 v=1.0 eta=0.3 | A.ALIGN(r=1.0)` |
| A | A.REPEL | 83 | 16 | `A kA=1 v=1.0 eta=0.0 | A.REPEL(b=any,r=0.5)` |
| A | A.ATTRACT | 61 | 16 | `A kA=1 v=1.0 eta=0.0 | A.ATTRACT(b=any,r=0.5)` |
| CA | A.ALIGN | 2153 | 17 | `CA kC=1 kA=1 v=1.0 eta=0.3 | A.ALIGN(r=1.0)` |
| AF | A.ALIGN | 7698 | 18 | `AF kA=1 v=1.0 eta=0.3 fields=1 D=[0.2] | A.ALIGN(r=1.0)` |
| CA | A.REPEL | 73 | 18 | `CA kC=1 kA=1 v=1.0 eta=0.0 | A.REPEL(b=any,r=1.0)` |
| A | A.ALIGN, A.ATTRACT | 1778 | 19 | `A kA=1 v=1.0 eta=0.3 | A.ATTRACT(b=any,r=1.0) ; A.ALIGN(r=1.0)` |
| AF | A.REPEL | 95 | 19 | `AF kA=1 v=1.0 eta=0.0 fields=1 D=[0.2] | A.REPEL(b=any,r=1.0)` |
| CA | A.ATTRACT | 54 | 19 | `CA kC=1 kA=1 v=1.0 eta=0.0 | A.ATTRACT(b=any,r=0.5)` |
| A | A.ALIGN, A.REPEL | 288 | 20 | `A kA=1 v=1.0 eta=0.3 | A.ALIGN(r=3.0) ; A.REPEL(b=any,r=0.5)` |
| AF | A.ATTRACT | 49 | 20 | `AF kA=1 v=1.0 eta=0.0 fields=1 D=[0.2] | A.ATTRACT(b=any,r=0.5)` |
| A | A.ALIGN, A.QUORUM | 784 | 21 | `A kA=1 v=1.0 eta=0.3 | A.ALIGN(r=1.0) ; A.QUORUM(th=1,r=1.0)` |
| CAF | A.ALIGN | 775 | 21 | `CAF kC=1 kA=1 v=1.0 eta=0.3 fields=1 D=[0.2] | A.ALIGN(r=1.0)` |
| CA | A.ALIGN, A.WRITECELL | 121 | 21 | `CA kC=1 kA=1 v=1.0 eta=0.3 | A.ALIGN(r=1.0) ; A.WRITECELL()` |
| A | A.ATTRACT, A.REPEL | 48 | 21 | `A kA=1 v=1.0 eta=0.0 | A.ATTRACT(b=any,r=1.0) ; A.REPEL(b=any,r=0.5)` |
| … | 6 more routes in triage.json | | | |

### power-law avalanches — 8 routes; ingredient common to all: C.SAND

| layers | rule templates | programs | bits | shortest example |
|---|---|---|---|---|
| C | C.SAND | 8 | 15 | `C kC=4 | C.SAND(p=0.001)` |
| CF | C.SAND | 44 | 18 | `CF kC=4 fields=1 D=[0.2] | C.SAND(p=0.001)` |
| CA | C.SAND | 81 | 21 | `CA kC=4 kA=1 v=1.0 eta=0.3 | C.SAND(p=0.001)` |
| C | C.SAND, C.SWAP | 6 | 21 | `C kC=4 | C.SWAP(p=1.0) ; C.SAND(p=0.001)` |
| C | C.COPY, C.SAND | 4 | 21 | `C kC=4 | C.COPY(p=1.0) ; C.SAND(p=0.001)` |
| C | C.KURA, C.SAND | 4 | 22 | `C kC=4 | C.SAND(p=0.001) ; C.KURA(K=1.0)` |
| C | C.FLIP, C.SAND | 2 | 22 | `C kC=4 | C.FLIP(a=1,b=0,p=1.0) ; C.SAND(p=0.1)` |
| C | C.SAND, C.SCHELL | 2 | 22 | `C kC=4 | C.SCHELL(th=0.7) ; C.SAND(p=0.001)` |

### spatial game chaos — 17 routes; ingredient common to all: C.GAME

| layers | rule templates | programs | bits | shortest example |
|---|---|---|---|---|
| CF | C.GAME | 24 | 15 | `CF kC=2 fields=1 D=[0.2] | C.GAME(g=(1.9, 0.0))` |
| C | C.GAME, C.SWAP | 46 | 17 | `C kC=2 | C.SWAP(p=1.0) ; C.GAME(g=(1.5, 0.5))` |
| C | C.COPY, C.GAME | 43 | 17 | `C kC=2 | C.COPY(p=1.0) ; C.GAME(g=(1.5, 0.5))` |
| C | C.GAME, C.MAJ | 35 | 17 | `C kC=2 | C.MAJ(p=0.1) ; C.GAME(g=(1.9, 0.0))` |
| C | C.GAME, C.KURA | 19 | 17 | `C kC=2 | C.GAME(g=(1.9, 0.0)) ; C.KURA(K=1.0)` |
| C | C.GAME | 2 | 17 | `C kC=2 | C.GAME(g=(1.9, 0.0)) ; C.GAME(g=(1.9, 0.0))` |
| CA | C.GAME | 92 | 18 | `CA kC=2 kA=1 v=1.0 eta=0.3 | C.GAME(g=(1.9, 0.0))` |
| C | C.CYCLE, C.GAME | 128 | 19 | `C kC=2 | C.CYCLE(th=1,p=1.0) ; C.GAME(g=(1.5, 0.5))` |
| C | C.FLIP, C.GAME | 85 | 19 | `C kC=2 | C.FLIP(a=0,b=1,p=0.1) ; C.GAME(g=(1.5, 0.5))` |
| CAF | C.GAME | 123 | 21 | `CAF kC=2 kA=1 v=1.0 eta=0.3 fields=1 D=[0.2] | C.GAME(g=(1.9, 0.0))` |
| CF | C.COPY, C.GAME | 30 | 21 | `CF kC=2 fields=1 D=[0.2] | C.COPY(p=1.0) ; C.GAME(g=(1.5, 0.5))` |
| C | C.CONTAG, C.GAME | 24 | 21 | `C kC=2 | C.CONTAG(a=0,b=1,th=1,p=0.1) ; C.GAME(g=(1.5, 0.5))` |
| CF | C.GAME, C.MAJ | 15 | 21 | `CF kC=2 fields=1 D=[0.2] | C.MAJ(p=0.1) ; C.GAME(g=(1.9, 0.0))` |
| CF | C.GAME, C.SWAP | 12 | 22 | `CF kC=2 fields=1 D=[0.2] | C.SWAP(p=1.0) ; C.GAME(g=(1.5, 0.5))` |
| CF | C.GAME, C.KURA | 4 | 22 | `CF kC=2 fields=1 D=[0.2] | C.GAME(g=(1.9, 0.0)) ; C.KURA(K=1.0)` |
| … | 2 more routes in triage.json | | | |

## Review — known behaviour, atypical look (31126)

- [10 bits] `N kN=4 | N.COPY(p=1.0)` view N — shows: network consensus
- [10 bits] `N kN=4 | N.COPY(p=0.1)` view N — shows: network consensus
- [11 bits] `C kC=3 | C.MAJ(p=1.0)` view C — shows: domain coarsening
- [11 bits] `N kN=3 | N.MAJ(p=0.1)` view N — shows: network consensus
- [11 bits] `N kN=4 | N.MAJ(p=1.0)` view N — shows: network consensus
- [11 bits] `N kN=4 | N.REWIRE(p=0.1)` view N — shows: network fragmentation
- [11 bits] `N kN=3 | N.REWIRE(p=1.0)` view N — shows: network fragmentation
- [11 bits] `N kN=3 | N.MAJ(p=1.0)` view N — shows: network consensus
- [11 bits] `N kN=3 | N.REWIRE(p=0.1)` view N — shows: network fragmentation
- [11 bits] `N kN=4 | N.REWIRE(p=1.0)` view N — shows: network fragmentation
- [12 bits] `C kC=1 | C.KURA(K=0.2)` view C.phase — shows: lattice phase locking
- [12 bits] `C kC=3 | C.SCHELL(th=0.2)` view C — shows: segregation with vacancies
- [12 bits] `N kN=1 | N.KURA(K=0.2)` view N.phase — shows: network synchronization
- [12 bits] `N kN=2 | N.GAME(g=(1.3, 0.0))` view N — shows: network consensus
- [12 bits] `N kN=4 | N.COPY(p=0.3)` view N — shows: network consensus
- [13 bits] `C kC=2 | C.MAJ(p=0.01)` view C — shows: domain coarsening
- [13 bits] `N kN=3 | N.REWIRE(p=0.3)` view N — shows: network fragmentation
- [13 bits] `N kN=4 | N.COPY(p=0.2)` view N — shows: network consensus
- [13 bits] `N kN=4 | N.REWIRE(p=0.03)` view N — shows: network fragmentation
- [13 bits] `C kC=2 | C.MAJ(p=0.003)` view C — shows: domain coarsening
- [13 bits] `N kN=3 | N.COPY(p=0.2)` view N — shows: network consensus
- [13 bits] `N kN=3 | N.REWIRE(p=0.03)` view N — shows: network fragmentation
- [13 bits] `C kC=3 | C.CYCLE(th=3,p=0.1)` view C — shows: domain coarsening
- [13 bits] `C kC=4 | C.CYCLE(th=3,p=0.1)` view C — shows: domain coarsening
- [13 bits] `N kN=1 | N.KURA(K=1.0) ; N.KURA(K=1.0)` view N.phase — shows: network synchronization

## Unnamed behaviour classes -> literature check (70)

| class | bits | shortest program | size |
|---|---|---|---|
| lattice-phase-1 | 9 | `C kC=1 | C.KURA(K=0.1)` | 9051 |
| network-phase-1 | 9 | `N kN=1 | N.KURA(K=0.1)` | 1219 |
| network-1 | 11 | `N kN=2 | N.GAME(g=(1.9, 0.0))` | 2643 |
| grid-1 | 11 | `C kC=2 | C.COPY(p=1.0)` | 2160 |
| grid-22 | 11 | `C kC=4 | C.COPY(p=1.0)` | 494 |
| grid-2 | 12 | `C kC=2 | C.GAME(g=(1.5, -0.5))` | 429 |
| grid-20 | 12 | `C kC=2 | C.GAME(g=(1.3, 0.0))` | 64 |
| grid-5 | 13 | `C kC=2 | C.COPY(p=0.003)` | 6631 |
| network-6 | 13 | `N kN=3 | N.CYCLE(th=3,p=0.1)` | 493 |
| network-15 | 13 | `N kN=3 | N.COPY(p=0.05)` | 486 |
| grid-4 | 13 | `C kC=4 | C.CYCLE(th=1,p=1.0)` | 383 |
| grid-3 | 13 | `C kC=3 | C.CYCLE(th=1,p=1.0)` | 344 |
| network-12 | 13 | `N kN=2 | N.REWIRE(p=0.01)` | 198 |
| network-4 | 13 | `N kN=3 | N.REWIRE(p=0.003)` | 195 |
| network-5 | 13 | `N kN=4 | N.CYCLE(th=1,p=1.0)` | 139 |
| network-2 | 13 | `N kN=3 | N.CYCLE(th=1,p=1.0)` | 68 |
| network-3 | 13 | `N kN=3 | N.REWIRE(p=0.01)` | 44 |
| agents-1 | 14 | `A kA=1 v=1.0 eta=0.3 | A.ATTRACT(b=any,r=1.0)` | 25261 |
| agents-2 | 14 | `A kA=1 v=1.0 eta=0.3 | A.REPEL(b=any,r=3.0)` | 4493 |
| agents-12 | 15 | `A kA=1 v=2.0 eta=0.3 | A.REPEL(b=any,r=3.0)` | 3962 |
| agents-6 | 15 | `A kA=1 v=1.0 eta=0.3 | A.ATTRACT(b=any,r=8.0)` | 2882 |
| agents-4 | 15 | `A kA=1 v=0.1 eta=0.3 | A.REPEL(b=any,r=1.0)` | 2287 |
| agents-3 | 15 | `A kA=1 v=1.0 eta=3.14159 | A.ALIGN(r=1.0)` | 2123 |
| grid-10 | 15 | `C kC=3 | C.CYCLE(th=1,p=0.01)` | 942 |
| agents-5 | 15 | `A kA=1 v=0.0 eta=1.0 | A.ALIGN(r=1.0)` | 528 |
| grid-8 | 15 | `C kC=3 | C.CYCLE(th=5,p=0.1)` | 369 |
| grid-7 | 15 | `C kC=2 | C.CYCLE(th=4,p=1.0)` | 186 |
| network-phase-4 | 15 | `N kN=1 | N.KURA(K=0.1) ; N.KURA(K=0.03)` | 103 |
| grid-6 | 15 | `C kC=2 | C.CONTAG(a=0,b=1,th=1,p=0.1)` | 97 |
| grid-28 | 15 | `C kC=2 | C.CONTAG(a=0,b=1,th=3,p=0.1)` | 40 |
| agents-8 | 16 | `AF kA=1 v=1.0 eta=0.3 fields=1 D=[0.05] | A.CLIMB(f=0)` | 1556 |
| agents-9 | 16 | `A kA=1 v=0.3 eta=3.14159 | A.ALIGN(r=2.0)` | 661 |
| agents-7 | 16 | `A kA=1 v=0.1 eta=1.0 | A.ALIGN(r=0.5)` | 479 |
| grid-9 | 16 | `C kC=2 | C.CYCLE(th=7,p=0.01)` | 116 |
| agents-16 | 17 | `AF kA=1 v=0.3 eta=0.3 fields=1 D=[0.15] | A.CLIMB(f=0)` | 1225 |
| agents-15 | 17 | `AF kA=1 v=0.3 eta=0.0 fields=1 D=[0.2] | A.CLIMB(f=0)` | 521 |
| agents-11 | 17 | `AF kA=1 v=0.1 eta=0.3 fields=1 D=[0.2] | A.CLIMB(f=0)` | 489 |
| agents-10 | 17 | `CA kC=1 kA=2 v=1.0 eta=1.0 | A.ALIGN(r=1.0)` | 396 |
| grid-11 | 17 | `C kC=2 | C.COPY(p=1.0) ; C.MAJ(p=1.0)` | 182 |
| network-phase-3 | 17 | `N kN=1 | N.KURA(K=0.3) ; N.KURA(K=0.03)` | 117 |
