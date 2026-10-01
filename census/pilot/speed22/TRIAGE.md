# Census triage — 2026-10-01 21:04

600 programs; 297 with a flagged view; errors 0; fingerprint errors 0 (must be 0).
Flagged views: 303 = 86 named + 64 known-behaviour-atypical-look (review) + 153 unnamed in 15 classes (literature check).

## Behaviour map (shortest program per behaviour = the emergence-complexity table)

| behaviour | catalog | bits | shortest program | programs | routes |
|---|---|---|---|---|---|
| domain coarsening | P18 | 17 | `C kC=4 | C.MAJ(p=0.1) ; C.SCHELL(th=0.3)` | 53 | 19 |
| segregation with vacancies | P1 | 17 | `C kC=4 | C.MAJ(p=0.1) ; C.SCHELL(th=0.3)` | 12 | 5 |
| flocking | P5 | 18 | `CA kC=2 kA=3 v=1.0 eta=0.3 | A.ALIGN(r=3.0)` | 41 | 9 |
| network consensus | P18 | 18 | `N kN=4 | N.COPY(p=1.0) ; N.COPY(p=0.002)` | 13 | 8 |
| lattice phase locking | P9 | 19 | `C kC=1 | C.KURA(K=1.0) ; C.KURA(K=0.03) ; C.COPY(p=1.0)` | 21 | 10 |
| network synchronization | P9 | 19 | `N kN=1 | N.KURA(K=0.1) ; N.COPY(p=0.1) ; N.KURA(K=0.3)` | 6 | 2 |
| network fragmentation | P34 | 19 | `N kN=4 | N.COPY(p=0.5) ; N.REWIRE(p=1.0)` | 8 | 5 |
| cyclic waves and spirals | — | 21 | `C kC=3 | C.CYCLE(th=2,p=0.1) ; C.COPY(p=0.1)` | 2 | 2 |
| spatial game chaos | P27 | 22 | `C kC=2 | C.GAME(g=(1.9, 0.0)) ; C.CYCLE(th=3,p=0.05)` | 1 | 1 |

## Routes (every distinct mechanism that produces each behaviour)

### domain coarsening — 19 routes; ingredient common to all: none

| layers | rule templates | programs | bits | shortest example |
|---|---|---|---|---|
| CF | C.MAJ | 8 | 17 | `CF kC=2 fields=1 D=[0.15] | C.MAJ(p=0.01)` |
| C | C.MAJ, C.SCHELL | 3 | 17 | `C kC=4 | C.MAJ(p=0.1) ; C.SCHELL(th=0.3)` |
| CA | C.COPY | 8 | 18 | `CA kC=4 kA=4 v=1.0 eta=0.3 | C.COPY(p=1.0)` |
| C | C.KURA, C.MAJ | 3 | 18 | `C kC=4 | C.KURA(K=0.1) ; C.MAJ(p=0.1)` |
| CF | C.COPY | 2 | 19 | `CF kC=4 fields=2 D=[0.2, 0.2] | C.COPY(p=0.5)` |
| C | C.MAJ, C.SWAP | 1 | 19 | `C kC=3 | C.MAJ(p=0.1) ; C.SWAP(p=0.03)` |
| CF | C.CYCLE | 5 | 20 | `CF kC=2 fields=1 D=[0.15] | C.CYCLE(th=3,p=0.002)` |
| C | C.CYCLE, C.MAJ | 2 | 20 | `C kC=4 | C.CYCLE(th=8,p=0.1) ; C.MAJ(p=0.1)` |
| C | C.COPY | 1 | 20 | `C kC=4 | C.COPY(p=0.0003) ; C.COPY(p=0.1)` |
| CA | C.MAJ | 8 | 21 | `CA kC=2 kA=4 v=0.1 eta=1.0 | C.MAJ(p=0.01)` |
| CAF | C.COPY | 2 | 21 | `CAF kC=2 kA=3 v=1.0 eta=0.3 fields=1 D=[0.2] | C.COPY(p=0.1)` |
| C | C.COPY, C.CYCLE | 1 | 21 | `C kC=3 | C.CYCLE(th=2,p=0.1) ; C.COPY(p=0.1)` |
| C | C.CYCLE, C.SCHELL | 1 | 21 | `C kC=4 | C.SCHELL(th=0.3) ; C.CYCLE(th=4,p=0.1)` |
| CF | C.MAJ, C.SCHELL | 2 | 22 | `CF kC=4 fields=1 D=[0.02] | C.MAJ(p=0.1) ; C.SCHELL(th=0.3)` |
| CAF | C.MAJ | 2 | 22 | `CAF kC=2 kA=3 v=0.0 eta=0.3 fields=1 D=[0.2] | C.MAJ(p=1.0)` |
| … | 4 more routes in triage.json | | | |

### segregation with vacancies — 5 routes; ingredient common to all: C.SCHELL

| layers | rule templates | programs | bits | shortest example |
|---|---|---|---|---|
| C | C.MAJ, C.SCHELL | 3 | 17 | `C kC=4 | C.MAJ(p=0.1) ; C.SCHELL(th=0.3)` |
| CA | C.SCHELL | 3 | 20 | `CA kC=3 kA=3 v=0.0 eta=0.1 | C.SCHELL(th=0.3)` |
| CAF | C.SCHELL | 3 | 21 | `CAF kC=3 kA=3 v=1.0 eta=0.3 fields=1 D=[0.2] | C.SCHELL(th=0.5)` |
| C | C.CYCLE, C.SCHELL | 1 | 21 | `C kC=4 | C.SCHELL(th=0.3) ; C.CYCLE(th=4,p=0.1)` |
| CF | C.MAJ, C.SCHELL | 2 | 22 | `CF kC=4 fields=1 D=[0.02] | C.MAJ(p=0.1) ; C.SCHELL(th=0.3)` |

### flocking — 9 routes; ingredient common to all: none

| layers | rule templates | programs | bits | shortest example |
|---|---|---|---|---|
| CA | A.ALIGN | 6 | 18 | `CA kC=2 kA=3 v=1.0 eta=0.3 | A.ALIGN(r=3.0)` |
| AF | A.ALIGN | 15 | 19 | `AF kA=3 v=2.0 eta=0.3 fields=1 D=[0.2] | A.ALIGN(r=1.0)` |
| A | A.ALIGN, A.REPEL | 2 | 20 | `A kA=1 v=1.0 eta=1.0 | A.ALIGN(r=3.0) ; A.REPEL(b=any,r=0.5)` |
| CA | A.REPEL | 1 | 20 | `CA kC=4 kA=1 v=1.0 eta=0.0 | A.REPEL(b=any,r=0.5)` |
| A | A.ALIGN | 8 | 21 | `A kA=2 v=1.0 eta=0.3 | A.ALIGN(r=2.0) ; A.ALIGN(r=0.5)` |
| A | A.ALIGN, A.QUORUM | 3 | 21 | `A kA=1 v=0.3 eta=0.3 | A.ALIGN(r=3.0) ; A.QUORUM(th=3,r=3.0)` |
| A | A.ALIGN, A.ATTRACT | 3 | 22 | `A kA=1 v=1.0 eta=0.0 | A.ALIGN(r=0.5) ; A.ATTRACT(b=any,r=0.5)` |
| A | A.ALIGN, A.CYCLE | 2 | 22 | `A kA=2 v=0.3 eta=0.3 | A.CYCLE(th=3,r=3.0) ; A.ALIGN(r=3.0)` |
| CAF | A.ALIGN | 1 | 22 | `CAF kC=3 kA=1 v=0.3 eta=0.3 fields=1 D=[0.05] | A.ALIGN(r=2.0)` |

### network consensus — 8 routes; ingredient common to all: none

| layers | rule templates | programs | bits | shortest example |
|---|---|---|---|---|
| N | N.COPY | 2 | 18 | `N kN=4 | N.COPY(p=1.0) ; N.COPY(p=0.002)` |
| N | N.COPY, N.REWIRE | 1 | 18 | `N kN=4 | N.COPY(p=0.3) ; N.REWIRE(p=0.1)` |
| N | N.KURA, N.MAJ | 1 | 19 | `N kN=4 | N.MAJ(p=1.0) ; N.KURA(K=0.01)` |
| N | N.COPY, N.KURA | 1 | 20 | `N kN=2 | N.COPY(p=1.0) ; N.KURA(K=0.05)` |
| N | N.COPY, N.MAJ | 1 | 21 | `N kN=4 | N.COPY(p=0.0003) ; N.MAJ(p=0.3)` |
| N | N.COPY, N.CYCLE | 3 | 22 | `N kN=2 | N.COPY(p=1.0) ; N.CYCLE(th=3,p=0.02)` |
| N | N.CYCLE, N.GAME | 3 | 22 | `N kN=2 | N.CYCLE(th=1,p=0.05) ; N.GAME(g=(1.5, 0.5))` |
| N | N.CYCLE, N.MAJ | 1 | 22 | `N kN=4 | N.MAJ(p=0.5) ; N.CYCLE(th=3,p=0.1)` |

### lattice phase locking — 10 routes; ingredient common to all: C.KURA

| layers | rule templates | programs | bits | shortest example |
|---|---|---|---|---|
| C | C.COPY, C.KURA | 5 | 19 | `C kC=1 | C.KURA(K=1.0) ; C.KURA(K=0.03) ; C.COPY(p=1.0)` |
| CF | C.EMIT, C.KURA | 2 | 19 | `CF kC=1 fields=1 D=[0.1] | C.KURA(K=1.0) ; C.EMIT(b=0,f=0,p=1.0)` |
| C | C.KURA, C.MAJ | 3 | 20 | `C kC=4 | C.KURA(K=1.0) ; C.MAJ(p=0.3)` |
| CA | C.KURA | 3 | 20 | `CA kC=3 kA=1 v=1.0 eta=0.3 | C.KURA(K=0.3)` |
| CF | C.KURA | 2 | 21 | `CF kC=2 fields=2 D=[0.05, 0.02] | C.KURA(K=0.2)` |
| CAF | C.KURA | 2 | 21 | `CAF kC=1 kA=4 v=1.0 eta=1.0 fields=1 D=[0.05] | C.KURA(K=1.0)` |
| C | C.KURA | 1 | 22 | `C kC=1 | C.KURA(K=1.0) ; C.KURA(K=0.01) ; C.KURA(K=0.005)` |
| C | C.KURA, C.SWAP | 1 | 22 | `C kC=3 | C.SWAP(p=0.01) ; C.KURA(K=0.5)` |
| CF | C.KURA, F.FEED | 1 | 22 | `CF kC=1 fields=1 D=[0.05] | F.FEED(f=0,F=1.0) ; C.KURA(K=0.5)` |
| CF | C.COPY, C.KURA | 1 | 22 | `CF kC=3 fields=1 D=[0.02] | C.COPY(p=1.0) ; C.KURA(K=1.0)` |

### network synchronization — 2 routes; ingredient common to all: N.KURA

| layers | rule templates | programs | bits | shortest example |
|---|---|---|---|---|
| N | N.COPY, N.KURA | 4 | 19 | `N kN=1 | N.KURA(K=0.1) ; N.COPY(p=0.1) ; N.KURA(K=0.3)` |
| N | N.CYCLE, N.KURA | 2 | 22 | `N kN=3 | N.CYCLE(th=3,p=0.1) ; N.KURA(K=0.2)` |

### network fragmentation — 5 routes; ingredient common to all: N.REWIRE

| layers | rule templates | programs | bits | shortest example |
|---|---|---|---|---|
| N | N.COPY, N.REWIRE | 2 | 19 | `N kN=4 | N.COPY(p=0.5) ; N.REWIRE(p=1.0)` |
| N | N.MAJ, N.REWIRE | 3 | 20 | `N kN=4 | N.MAJ(p=0.02) ; N.REWIRE(p=0.1)` |
| N | N.CYCLE, N.REWIRE | 1 | 21 | `N kN=2 | N.REWIRE(p=0.3) ; N.CYCLE(th=3,p=1.0)` |
| N | N.CONTAG, N.REWIRE | 1 | 22 | `N kN=3 | N.CONTAG(a=0,b=1,th=1,p=1.0) ; N.REWIRE(p=0.1)` |
| N | N.FLIP, N.REWIRE | 1 | 22 | `N kN=3 | N.REWIRE(p=1.0) ; N.FLIP(a=0,b=1,p=0.01)` |

### cyclic waves and spirals — 2 routes; ingredient common to all: C.CYCLE

| layers | rule templates | programs | bits | shortest example |
|---|---|---|---|---|
| C | C.COPY, C.CYCLE | 1 | 21 | `C kC=3 | C.CYCLE(th=2,p=0.1) ; C.COPY(p=0.1)` |
| C | C.CYCLE, C.SWAP | 1 | 21 | `C kC=4 | C.SWAP(p=0.1) ; C.CYCLE(th=1,p=0.3)` |

### spatial game chaos — 1 routes; ingredient common to all: C.CYCLE, C.GAME

| layers | rule templates | programs | bits | shortest example |
|---|---|---|---|---|
| C | C.CYCLE, C.GAME | 1 | 22 | `C kC=2 | C.GAME(g=(1.9, 0.0)) ; C.CYCLE(th=3,p=0.05)` |

## Review — known behaviour, atypical look (64)

- [17 bits] `CF kC=2 fields=1 D=[0.15] | C.MAJ(p=0.01)` view C — shows: domain coarsening
- [18 bits] `N kN=4 | N.COPY(p=0.3) ; N.REWIRE(p=0.1)` view N — shows: network consensus
- [18 bits] `CA kC=2 kA=3 v=1.0 eta=0.3 | A.ALIGN(r=3.0)` view A — shows: flocking
- [18 bits] `N kN=4 | N.COPY(p=1.0) ; N.COPY(p=0.002)` view N — shows: network consensus
- [18 bits] `CA kC=4 kA=4 v=1.0 eta=0.3 | C.COPY(p=1.0)` view C — shows: domain coarsening
- [18 bits] `CF kC=3 fields=1 D=[0.02] | C.MAJ(p=0.002)` view C — shows: domain coarsening
- [19 bits] `CA kC=1 kA=1 v=0.3 eta=0.0 | A.ALIGN(r=8.0)` view A — shows: flocking
- [19 bits] `N kN=4 | N.COPY(p=0.5) ; N.REWIRE(p=1.0)` view N — shows: network fragmentation
- [19 bits] `CA kC=1 kA=4 v=1.0 eta=0.0 | A.ALIGN(r=0.5)` view A — shows: flocking
- [19 bits] `N kN=1 | N.KURA(K=0.1) ; N.COPY(p=0.1) ; N.KURA(K=0.3)` view N.phase — shows: network synchronization
- [19 bits] `CA kC=3 kA=3 v=1.0 eta=1.0 | A.ALIGN(r=5.0)` view A — shows: flocking
- [19 bits] `CF kC=4 fields=2 D=[0.1, 0.05] | C.MAJ(p=0.003)` view C — shows: domain coarsening
- [19 bits] `N kN=4 | N.MAJ(p=1.0) ; N.KURA(K=0.01)` view N — shows: network consensus
- [19 bits] `AF kA=3 v=2.0 eta=0.3 fields=1 D=[0.2] | A.ALIGN(r=1.0)` view A — shows: flocking
- [20 bits] `AF kA=1 v=0.3 eta=1.0 fields=1 D=[0.1] | A.ALIGN(r=8.0)` view A — shows: flocking
- [20 bits] `A kA=1 v=1.0 eta=1.0 | A.ALIGN(r=3.0) ; A.REPEL(b=any,r=0.5)` view A — shows: flocking
- [20 bits] `CA kC=2 kA=1 v=2.0 eta=0.3 | A.ALIGN(r=2.0)` view A — shows: flocking
- [20 bits] `CA kC=4 kA=1 v=1.0 eta=0.0 | A.REPEL(b=any,r=0.5)` view A — shows: flocking
- [20 bits] `N kN=4 | N.MAJ(p=0.02) ; N.REWIRE(p=0.1)` view N — shows: network fragmentation
- [20 bits] `CA kC=4 kA=3 v=2.0 eta=1.0 | A.ALIGN(r=8.0)` view A — shows: flocking
- [20 bits] `AF kA=1 v=1.0 eta=0.3 fields=2 D=[0.2, 0.05] | A.ALIGN(r=3.0)` view A — shows: flocking
- [20 bits] `AF kA=1 v=0.3 eta=0.1 fields=1 D=[0.2] | A.ALIGN(r=5.0)` view A — shows: flocking
- [21 bits] `N kN=4 | N.COPY(p=0.0003) ; N.MAJ(p=0.3)` view N — shows: network consensus
- [21 bits] `AF kA=4 v=2.0 eta=1.0 fields=1 D=[0.15] | A.ALIGN(r=5.0)` view A — shows: flocking
- [21 bits] `C kC=4 | C.SCHELL(th=0.3) ; C.CYCLE(th=4,p=0.1)` view C — shows: domain coarsening, segregation with vacancies

## Unnamed behaviour classes -> literature check (15)

| class | bits | shortest program | size |
|---|---|---|---|
| network-1 | 14 | `N kN=3 | N.MAJ(p=0.0003)` | 7 |
| agents-3 | 17 | `A kA=1 v=1.0 eta=0.3 | A.QUORUM(th=1,r=0.5)` | 26 |
| grid-1 | 17 | `C kC=2 | C.SWAP(p=1.0) ; C.COPY(p=0.1)` | 21 |
| agents-4 | 18 | `A kA=2 v=0.3 eta=0.1 | A.ATTRACT(b=0,r=5.0)` | 23 |
| agents-7 | 18 | `A kA=2 v=2.0 eta=0.1 | A.ATTRACT(b=any,r=5.0)` | 20 |
| grid-4 | 18 | `CA kC=2 kA=4 v=0.3 eta=1.0 | C.COPY(p=1.0)` | 8 |
| agents-1 | 19 | `AF kA=3 v=0.1 eta=0.3 fields=1 D=[0.1] | A.CLIMB(f=0)` | 7 |
| grid-3 | 19 | `CA kC=2 kA=2 v=0.0 eta=1.0 | C.GAME(g=(1.5, 0.5))` | 3 |
| grid-2 | 19 | `CF kC=4 fields=1 D=[0.15] | C.SAND(p=0.0003)` | 2 |
| lattice-phase-1 | 20 | `C kC=1 | C.KURA(K=0.1) ; C.KURA(K=0.1) ; C.KURA(K=0.005)` | 8 |
| agents-5 | 21 | `A kA=2 v=0.1 eta=0.3 | A.ALIGN(r=1.0) ; A.ATTRACT(b=any,r=1.0)` | 12 |
| agents-2 | 21 | `AF kA=2 v=1.0 eta=3.14159 fields=1 D=[0.02] | A.ALIGN(r=2.0)` | 9 |
| grid-5 | 21 | `C kC=3 | C.CONTAG(a=0,b=1,th=2,p=0.005)` | 1 |
| agents-6 | 22 | `A kA=1 v=1.0 eta=0.0 | A.REPEL(b=any,r=0.5) ; A.REPEL(b=any,r=8.0)` | 4 |
| network-phase-1 | 22 | `N kN=1 | N.KURA(K=0.1) ; N.COPY(p=0.001) ; N.COPY(p=0.3)` | 2 |
