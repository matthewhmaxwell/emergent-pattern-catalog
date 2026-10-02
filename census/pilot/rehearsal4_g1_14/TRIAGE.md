# Census triage — 2026-10-02 14:23

914 programs; 407 with a flagged view; errors 0; fingerprint errors 0 (must be 0).
Flagged views: 407 = 191 named + 57 known-behaviour-atypical-look (review) + 159 unnamed in 18 classes (literature check).

## Behaviour map (shortest program per behaviour = the emergence-complexity table)

| behaviour | catalog | bits | shortest program | programs | routes |
|---|---|---|---|---|---|
| network synchronization | P9 | 9 | `N kN=1 | N.KURA(K=1.0)` | 23 | 2 |
| lattice phase locking | P9 | 9 | `C kC=1 | C.KURA(K=1.0)` | 27 | 3 |
| network consensus | P18 | 10 | `N kN=3 | N.COPY(p=1.0)` | 43 | 4 |
| segregation with vacancies | P1 | 11 | `C kC=3 | C.SCHELL(th=0.5)` | 11 | 1 |
| network fragmentation | P34 | 11 | `N kN=2 | N.REWIRE(p=1.0)` | 24 | 1 |
| domain coarsening | P18 | 11 | `C kC=2 | C.COPY(p=0.1)` | 88 | 6 |
| cyclic waves and spirals | — | 13 | `C kC=3 | C.CYCLE(th=3,p=1.0)` | 3 | 1 |
| flocking | P5 | 14 | `A kA=1 v=1.0 eta=0.3 | A.ALIGN(r=1.0)` | 32 | 1 |

## Routes (every distinct mechanism that produces each behaviour)

### network synchronization — 2 routes; ingredient common to all: N.KURA

| layers | rule templates | programs | bits | shortest example |
|---|---|---|---|---|
| N | N.KURA | 19 | 9 | `N kN=1 | N.KURA(K=1.0)` |
| N | N.COPY, N.KURA | 4 | 13 | `N kN=1 | N.COPY(p=1.0) ; N.KURA(K=1.0)` |

### lattice phase locking — 3 routes; ingredient common to all: C.KURA

| layers | rule templates | programs | bits | shortest example |
|---|---|---|---|---|
| C | C.KURA | 17 | 9 | `C kC=1 | C.KURA(K=1.0)` |
| CF | C.KURA | 6 | 13 | `CF kC=1 fields=1 D=[0.2] | C.KURA(K=1.0)` |
| C | C.COPY, C.KURA | 4 | 13 | `C kC=1 | C.COPY(p=1.0) ; C.KURA(K=1.0)` |

### network consensus — 4 routes; ingredient common to all: none

| layers | rule templates | programs | bits | shortest example |
|---|---|---|---|---|
| N | N.COPY | 15 | 10 | `N kN=3 | N.COPY(p=1.0)` |
| N | N.MAJ | 23 | 11 | `N kN=2 | N.MAJ(p=1.0)` |
| N | N.GAME | 4 | 11 | `N kN=2 | N.GAME(g=(1.5, 0.5))` |
| N | N.CYCLE | 1 | 13 | `N kN=2 | N.CYCLE(th=3,p=0.1)` |

### segregation with vacancies — 1 routes; ingredient common to all: C.SCHELL

| layers | rule templates | programs | bits | shortest example |
|---|---|---|---|---|
| C | C.SCHELL | 11 | 11 | `C kC=3 | C.SCHELL(th=0.5)` |

### network fragmentation — 1 routes; ingredient common to all: N.REWIRE

| layers | rule templates | programs | bits | shortest example |
|---|---|---|---|---|
| N | N.REWIRE | 24 | 11 | `N kN=2 | N.REWIRE(p=1.0)` |

### domain coarsening — 6 routes; ingredient common to all: none

| layers | rule templates | programs | bits | shortest example |
|---|---|---|---|---|
| C | C.MAJ | 38 | 11 | `C kC=2 | C.MAJ(p=0.1)` |
| C | C.COPY | 26 | 11 | `C kC=2 | C.COPY(p=0.1)` |
| C | C.SCHELL | 2 | 12 | `C kC=3 | C.SCHELL(th=0.7)` |
| C | C.CYCLE | 4 | 13 | `C kC=2 | C.CYCLE(th=3,p=0.1)` |
| CF | C.MAJ | 11 | 14 | `CF kC=2 fields=1 D=[0.2] | C.MAJ(p=1.0)` |
| CF | C.COPY | 7 | 14 | `CF kC=2 fields=1 D=[0.2] | C.COPY(p=0.1)` |

### cyclic waves and spirals — 1 routes; ingredient common to all: C.CYCLE

| layers | rule templates | programs | bits | shortest example |
|---|---|---|---|---|
| C | C.CYCLE | 3 | 13 | `C kC=3 | C.CYCLE(th=3,p=1.0)` |

### flocking — 1 routes; ingredient common to all: A.ALIGN

| layers | rule templates | programs | bits | shortest example |
|---|---|---|---|---|
| A | A.ALIGN | 32 | 14 | `A kA=1 v=1.0 eta=0.3 | A.ALIGN(r=1.0)` |

## Review — known behaviour, atypical look (57)

- [10 bits] `N kN=4 | N.COPY(p=1.0)` view N — shows: network consensus
- [10 bits] `N kN=4 | N.COPY(p=0.1)` view N — shows: network consensus
- [11 bits] `N kN=3 | N.REWIRE(p=1.0)` view N — shows: network fragmentation
- [11 bits] `N kN=4 | N.REWIRE(p=1.0)` view N — shows: network fragmentation
- [11 bits] `N kN=3 | N.REWIRE(p=0.1)` view N — shows: network fragmentation
- [11 bits] `N kN=4 | N.REWIRE(p=0.1)` view N — shows: network fragmentation
- [11 bits] `N kN=3 | N.MAJ(p=1.0)` view N — shows: network consensus
- [11 bits] `N kN=4 | N.MAJ(p=1.0)` view N — shows: network consensus
- [11 bits] `N kN=3 | N.MAJ(p=0.1)` view N — shows: network consensus
- [11 bits] `C kC=3 | C.MAJ(p=1.0)` view C — shows: domain coarsening
- [12 bits] `C kC=3 | C.SCHELL(th=0.2)` view C — shows: segregation with vacancies
- [12 bits] `N kN=4 | N.COPY(p=0.3)` view N — shows: network consensus
- [12 bits] `C kC=1 | C.KURA(K=0.2)` view C.phase — shows: lattice phase locking
- [12 bits] `N kN=2 | N.GAME(g=(1.3, 0.0))` view N — shows: network consensus
- [12 bits] `N kN=1 | N.KURA(K=0.2)` view N.phase — shows: network synchronization
- [13 bits] `C kC=3 | C.MAJ(p=0.003)` view C — shows: domain coarsening
- [13 bits] `C kC=4 | C.CYCLE(th=3,p=0.1)` view C — shows: domain coarsening
- [13 bits] `C kC=1 | C.KURA(K=1.0) ; C.KURA(K=1.0)` view C.phase — shows: lattice phase locking
- [13 bits] `C kC=2 | C.MAJ(p=0.01)` view C — shows: domain coarsening
- [13 bits] `N kN=4 | N.COPY(p=0.5)` view N — shows: network consensus
- [13 bits] `C kC=2 | C.MAJ(p=0.003)` view C — shows: domain coarsening
- [13 bits] `C kC=4 | C.MAJ(p=0.003)` view C — shows: domain coarsening
- [13 bits] `N kN=3 | N.COPY(p=0.2)` view N — shows: network consensus
- [13 bits] `N kN=4 | N.COPY(p=0.2)` view N — shows: network consensus
- [13 bits] `N kN=1 | N.KURA(K=1.0) ; N.KURA(K=1.0)` view N.phase — shows: network synchronization

## Unnamed behaviour classes -> literature check (18)

| class | bits | shortest program | size |
|---|---|---|---|
| lattice-phase-1 | 9 | `C kC=1 | C.KURA(K=0.1)` | 15 |
| network-phase-1 | 9 | `N kN=1 | N.KURA(K=0.1)` | 8 |
| network-1 | 11 | `N kN=2 | N.GAME(g=(1.9, 0.0))` | 47 |
| lattice-phase-2 | 11 | `C kC=1 | C.KURA(K=0.03)` | 11 |
| grid-1 | 11 | `C kC=2 | C.COPY(p=1.0)` | 9 |
| grid-3 | 11 | `C kC=2 | C.GAME(g=(1.5, 0.5))` | 2 |
| network-5 | 12 | `N kN=2 | N.GAME(g=(1.5, -0.5))` | 11 |
| grid-2 | 12 | `C kC=2 | C.GAME(g=(1.5, -0.5))` | 1 |
| grid-4 | 13 | `C kC=2 | C.COPY(p=0.003)` | 26 |
| network-4 | 13 | `N kN=2 | N.REWIRE(p=0.01)` | 10 |
| network-2 | 13 | `N kN=3 | N.CYCLE(th=1,p=1.0)` | 2 |
| grid-5 | 13 | `C kC=3 | C.CYCLE(th=1,p=1.0)` | 1 |
| grid-6 | 13 | `C kC=4 | C.CYCLE(th=1,p=1.0)` | 1 |
| network-3 | 13 | `N kN=3 | N.REWIRE(p=0.01)` | 1 |
| network-6 | 13 | `N kN=4 | N.CYCLE(th=1,p=1.0)` | 1 |
| agents-1 | 14 | `A kA=1 v=1.0 eta=0.3 | A.ATTRACT(b=any,r=1.0)` | 8 |
| agents-2 | 14 | `A kA=1 v=1.0 eta=0.3 | A.REPEL(b=any,r=3.0)` | 4 |
| grid-7 | 14 | `C kC=2 | C.CYCLE(th=7,p=0.1)` | 1 |
