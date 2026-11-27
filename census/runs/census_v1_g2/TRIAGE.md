# Census triage — 2026-10-24 05:55

98,082 programs; 21,172 with a flagged view; errors 0; fingerprint errors 0 (must be 0).
Flagged views: 21,277 = 4,122 named + 4,024 known-behaviour-atypical-look (review) + 13,131 unnamed in 61 classes (literature check).

## Behaviour map (shortest program per behaviour = the emergence-complexity table)

| behaviour | catalog | bits | shortest program | programs | routes |
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

## Routes (every distinct mechanism that produces each behaviour)

### network synchronization — 6 routes; ingredient common to all: none

| layers | rule templates | programs | bits | shortest example |
|---|---|---|---|---|
| N | N:ALWAYS>PHASE_COUPLE | 458 | 11 | `N kN=1 | N[any] ALWAYS() -> PHASE_COUPLE(K=1.0) p=1.0` |
| N | N:ALWAYS>COPY_RAND, N:ALWAYS>PHASE_COUPLE | 675 | 15 | `N kN=1 | N[any] ALWAYS() -> COPY_RAND() p=1.0 ; N[any] ALWAYS() -> PHASE_COUPLE(K=1.0) p=1.0` |
| N | N:DISCORD>PHASE_COUPLE | 33 | 17 | `N kN=2 | N[any] DISCORD() -> PHASE_COUPLE(K=1.0) p=1.0` |
| N | N:CNT_LE>PHASE_COUPLE | 21 | 19 | `N kN=2 | N[any] CNT_LE(b=0,th=3) -> PHASE_COUPLE(K=1.0) p=0.1` |
| N | N:CNT_GE>PHASE_COUPLE | 10 | 19 | `N kN=2 | N[any] CNT_GE(b=0,th=1) -> PHASE_COUPLE(K=1.0) p=1.0` |
| N | N:NEXT_GE>PHASE_COUPLE | 8 | 19 | `N kN=2 | N[any] NEXT_GE(th=1) -> PHASE_COUPLE(K=1.0) p=1.0` |

### lattice phase locking — 10 routes; ingredient common to all: none

| layers | rule templates | programs | bits | shortest example |
|---|---|---|---|---|
| C | C:ALWAYS>PHASE_COUPLE | 490 | 11 | `C kC=1 | C[any] ALWAYS() -> PHASE_COUPLE(K=1.0) p=1.0` |
| C | C:ALWAYS>COPY_RAND, C:ALWAYS>PHASE_COUPLE | 748 | 15 | `C kC=1 | C[any] ALWAYS() -> COPY_RAND() p=1.0 ; C[any] ALWAYS() -> PHASE_COUPLE(K=1.0) p=1.0` |
| CF | C:ALWAYS>PHASE_COUPLE | 206 | 17 | `CF kC=1 fields=1 D=[0.2] | C[any] ALWAYS() -> PHASE_COUPLE(K=1.0) p=1.0` |
| CA | C:ALWAYS>PHASE_COUPLE | 368 | 19 | `CA kC=1 kA=1 v=1.0 eta=0.3 | C[any] ALWAYS() -> PHASE_COUPLE(K=1.0) p=1.0` |
| CA | C:AGENT_ON>PHASE_COUPLE | 263 | 19 | `CA kC=1 kA=1 v=1.0 eta=0.3 | C[any] AGENT_ON() -> PHASE_COUPLE(K=1.0) p=1.0` |
| CF | C:FIELD_GT>PHASE_COUPLE | 92 | 19 | `CF kC=1 fields=1 D=[0.2] | C[any] FIELD_GT(f=0,th=0.5) -> PHASE_COUPLE(K=1.0) p=1.0` |
| C | C:CNT_GE>PHASE_COUPLE | 19 | 19 | `C kC=2 | C[any] CNT_GE(b=0,th=1) -> PHASE_COUPLE(K=1.0) p=1.0` |
| C | C:NEXT_GE>PHASE_COUPLE | 19 | 19 | `C kC=2 | C[any] NEXT_GE(th=1) -> PHASE_COUPLE(K=1.0) p=1.0` |
| C | C:SAME_LT>PHASE_COUPLE | 16 | 20 | `C kC=2 | C[any] SAME_LT(th=0.7) -> PHASE_COUPLE(K=1.0) p=1.0` |
| C | C:CNT_LE>PHASE_COUPLE | 15 | 20 | `C kC=2 | C[any] CNT_LE(b=0,th=7) -> PHASE_COUPLE(K=1.0) p=1.0` |

### domain coarsening — 46 routes; ingredient common to all: none

| layers | rule templates | programs | bits | shortest example |
|---|---|---|---|---|
| C | C:ALWAYS>COPY_MAJ | 40 | 14 | `C kC=2 | C[any] ALWAYS() -> COPY_MAJ() p=1.0` |
| C | C:ALWAYS>COPY_RAND | 28 | 14 | `C kC=2 | C[any] ALWAYS() -> COPY_RAND() p=0.1` |
| C | C:SAME_LT>COPY_MAJ | 160 | 17 | `C kC=2 | C[any] SAME_LT(th=0.5) -> COPY_MAJ() p=0.1` |
| C | C:SAME_LT>COPY_RAND | 157 | 17 | `C kC=2 | C[any] SAME_LT(th=0.5) -> COPY_RAND() p=0.1` |
| C | C:SAME_LT>SET_NEXT | 136 | 17 | `C kC=2 | C[any] SAME_LT(th=0.5) -> SET_NEXT() p=0.1` |
| C | C:CNT_LE>COPY_MAJ | 119 | 17 | `C kC=2 | C[any] CNT_LE(b=0,th=3) -> COPY_MAJ() p=1.0` |
| C | C:NEXT_GE>COPY_MAJ | 103 | 17 | `C kC=2 | C[any] NEXT_GE(th=1) -> COPY_MAJ() p=1.0` |
| C | C:NEXT_GE>COPY_RAND | 86 | 17 | `C kC=2 | C[any] NEXT_GE(th=1) -> COPY_RAND() p=1.0` |
| C | C:CNT_GE>COPY_MAJ | 73 | 17 | `C kC=2 | C[any] CNT_GE(b=0,th=1) -> COPY_MAJ() p=1.0` |
| C | C:SAME_LT>SET | 70 | 17 | `C kC=2 | C[any] SAME_LT(th=0.5) -> SET(b=0) p=0.1` |
| C | C:CNT_GE>COPY_RAND | 57 | 17 | `C kC=2 | C[any] CNT_GE(b=0,th=1) -> COPY_RAND() p=1.0` |
| C | C:NEXT_GE>SET_NEXT | 56 | 17 | `C kC=2 | C[any] NEXT_GE(th=3) -> SET_NEXT() p=0.1` |
| C | C:ALWAYS>IMITATE_BEST | 48 | 17 | `C kC=2 | C[any] ALWAYS() -> IMITATE_BEST(g=(1.5, -0.5)) p=0.1` |
| C | C:CNT_LE>SET | 40 | 17 | `C kC=2 | C[any] CNT_LE(b=0,th=3) -> SET(b=1) p=1.0` |
| C | C:SAME_LT>SWAP_RAND | 40 | 17 | `C kC=2 | C[any] SAME_LT(th=0.5) -> SWAP_RAND() p=1.0` |
| … | 31 more routes in triage.json | | | |

### network fragmentation — 5 routes; ingredient common to all: none

| layers | rule templates | programs | bits | shortest example |
|---|---|---|---|---|
| N | N:ALWAYS>REWIRE_SAME | 42 | 14 | `N kN=2 | N[any] ALWAYS() -> REWIRE_SAME() p=1.0` |
| N | N:DISCORD>REWIRE_SAME | 23 | 15 | `N kN=2 | N[any] DISCORD() -> REWIRE_SAME() p=1.0` |
| N | N:CNT_GE>REWIRE_SAME | 43 | 17 | `N kN=2 | N[any] CNT_GE(b=0,th=1) -> REWIRE_SAME() p=1.0` |
| N | N:NEXT_GE>REWIRE_SAME | 28 | 17 | `N kN=2 | N[any] NEXT_GE(th=1) -> REWIRE_SAME() p=1.0` |
| N | N:CNT_LE>REWIRE_SAME | 101 | 18 | `N kN=2 | N[any] CNT_LE(b=0,th=7) -> REWIRE_SAME() p=1.0` |

### network consensus — 23 routes; ingredient common to all: none

| layers | rule templates | programs | bits | shortest example |
|---|---|---|---|---|
| N | N:ALWAYS>COPY_MAJ | 36 | 14 | `N kN=2 | N[any] ALWAYS() -> COPY_MAJ() p=1.0` |
| N | N:ALWAYS>COPY_RAND | 25 | 14 | `N kN=2 | N[any] ALWAYS() -> COPY_RAND() p=1.0` |
| N | N:DISCORD>SET | 48 | 15 | `N kN=2 | N[any] DISCORD() -> SET(b=0) p=1.0` |
| N | N:DISCORD>COPY_RAND | 34 | 15 | `N kN=2 | N[any] DISCORD() -> COPY_RAND() p=1.0` |
| N | N:DISCORD>COPY_MAJ | 30 | 15 | `N kN=2 | N[any] DISCORD() -> COPY_MAJ() p=1.0` |
| N | N:DISCORD>SET_NEXT | 18 | 15 | `N kN=2 | N[any] DISCORD() -> SET_NEXT() p=1.0` |
| N | N:ALWAYS>IMITATE_BEST | 115 | 16 | `N kN=2 | N[any] ALWAYS() -> IMITATE_BEST(g=(1.5, 0.5)) p=1.0` |
| N | N:CNT_LE>SET | 257 | 17 | `N kN=2 | N[any] CNT_LE(b=0,th=1) -> SET(b=1) p=0.1` |
| N | N:CNT_LE>COPY_MAJ | 154 | 17 | `N kN=2 | N[any] CNT_LE(b=0,th=3) -> COPY_MAJ() p=1.0` |
| N | N:CNT_LE>COPY_RAND | 120 | 17 | `N kN=2 | N[any] CNT_LE(b=0,th=3) -> COPY_RAND() p=1.0` |
| N | N:DISCORD>IMITATE_BEST | 84 | 17 | `N kN=2 | N[any] DISCORD() -> IMITATE_BEST(g=(1.5, 0.5)) p=1.0` |
| N | N:CNT_GE>SET | 74 | 17 | `N kN=2 | N[any] CNT_GE(b=0,th=1) -> SET(b=0) p=1.0` |
| N | N:CNT_GE>COPY_MAJ | 67 | 17 | `N kN=2 | N[any] CNT_GE(b=0,th=1) -> COPY_MAJ() p=1.0` |
| N | N:NEXT_GE>SET | 60 | 17 | `N kN=2 | N[any] NEXT_GE(th=1) -> SET(b=0) p=1.0` |
| N | N:CNT_GE>COPY_RAND | 51 | 17 | `N kN=2 | N[any] CNT_GE(b=0,th=1) -> COPY_RAND() p=1.0` |
| … | 8 more routes in triage.json | | | |

### spatial game chaos — 6 routes; ingredient common to all: none

| layers | rule templates | programs | bits | shortest example |
|---|---|---|---|---|
| C | C:ALWAYS>IMITATE_BEST | 19 | 16 | `C kC=2 | C[any] ALWAYS() -> IMITATE_BEST(g=(1.9, 0.0)) p=1.0` |
| C | C:CNT_GE>IMITATE_BEST | 29 | 19 | `C kC=2 | C[any] CNT_GE(b=0,th=1) -> IMITATE_BEST(g=(1.5, 0.5)) p=0.1` |
| C | C:CNT_LE>IMITATE_BEST | 25 | 19 | `C kC=2 | C[any] CNT_LE(b=0,th=3) -> IMITATE_BEST(g=(1.9, 0.0)) p=0.1` |
| C | C:NEXT_GE>IMITATE_BEST | 14 | 19 | `C kC=2 | C[any] NEXT_GE(th=1) -> IMITATE_BEST(g=(1.9, 0.0)) p=1.0` |
| CF | C:ALWAYS>IMITATE_BEST | 13 | 20 | `CF kC=2 fields=1 D=[0.2] | C[any] ALWAYS() -> IMITATE_BEST(g=(1.9, 0.0)) p=1.0` |
| C | C:SAME_LT>IMITATE_BEST | 3 | 20 | `C kC=2 | C[any] SAME_LT(th=0.7) -> IMITATE_BEST(g=(1.9, 0.0)) p=1.0` |

### flocking — 4 routes; ingredient common to all: none

| layers | rule templates | programs | bits | shortest example |
|---|---|---|---|---|
| A | A:ALWAYS>HEAD_MEAN | 771 | 17 | `A kA=1 v=1.0 eta=0.3 | A[any] ALWAYS() -> HEAD_MEAN(b=any,r=1.0) p=1.0` |
| A | A:ALWAYS>HEAD_TO | 16 | 19 | `A kA=1 v=1.0 eta=0.0 | A[any] ALWAYS() -> HEAD_TO(b=any,r=0.5) p=1.0` |
| A | A:ALWAYS>HEAD_AWAY | 13 | 20 | `A kA=1 v=1.0 eta=0.0 | A[any] ALWAYS() -> HEAD_AWAY(b=any,r=0.5) p=1.0` |
| CA | A:ALWAYS>HEAD_MEAN | 43 | 21 | `CA kC=1 kA=1 v=1.0 eta=0.3 | A[any] ALWAYS() -> HEAD_MEAN(b=any,r=1.0) p=1.0` |

### cyclic waves and spirals — 3 routes; ingredient common to all: none

| layers | rule templates | programs | bits | shortest example |
|---|---|---|---|---|
| C | C:NEXT_GE>SET_NEXT | 11 | 18 | `C kC=3 | C[any] NEXT_GE(th=3) -> SET_NEXT() p=1.0` |
| C | C:NEXT_GE>COPY_RAND | 7 | 18 | `C kC=3 | C[any] NEXT_GE(th=1) -> COPY_RAND() p=1.0` |
| C | C:SAME_LT>SET_NEXT | 1 | 19 | `C kC=4 | C[any] SAME_LT(th=0.6) -> SET_NEXT() p=0.1` |

### power-law avalanches — 5 routes; ingredient common to all: none

| layers | rule templates | programs | bits | shortest example |
|---|---|---|---|---|
| C | C:ALWAYS>ADD_GRAIN | 2 | 18 | `C kC=4 | C[any] ALWAYS() -> ADD_GRAIN() p=0.001` |
| C | C:SAME_LT>ADD_GRAIN | 11 | 19 | `C kC=4 | C[any] SAME_LT(th=0.2) -> ADD_GRAIN() p=0.1` |
| C | C:NEXT_GE>ADD_GRAIN | 5 | 20 | `C kC=4 | C[any] NEXT_GE(th=1) -> ADD_GRAIN() p=0.003` |
| C | C:CNT_GE>ADD_GRAIN | 3 | 21 | `C kC=4 | C[any] CNT_GE(b=0,th=1) -> ADD_GRAIN() p=0.003` |
| C | C:CNT_LE>ADD_GRAIN | 3 | 21 | `C kC=4 | C[any] CNT_LE(b=1,th=1) -> ADD_GRAIN() p=0.003` |

### segregation with vacancies — 3 routes; ingredient common to all: none

| layers | rule templates | programs | bits | shortest example |
|---|---|---|---|---|
| C | C:SAME_LT>MOVE_EMPTY | 83 | 18 | `C kC=3 | C[any] SAME_LT(th=0.5) -> MOVE_EMPTY() p=1.0` |
| C | C:NEXT_GE>MOVE_EMPTY | 24 | 18 | `C kC=3 | C[any] NEXT_GE(th=3) -> MOVE_EMPTY() p=0.1` |
| C | C:CNT_LE>MOVE_EMPTY | 1 | 21 | `C kC=3 | C[1] CNT_LE(b=1,th=2) -> MOVE_EMPTY() p=0.1` |

## Review — known behaviour, atypical look (4024)

- [11 bits] `C kC=1 | C[any] ALWAYS() -> PHASE_COUPLE(K=1.0) p=0.1` view C.phase — shows: lattice phase locking
- [13 bits] `C kC=1 | C[any] ALWAYS() -> PHASE_COUPLE(K=1.0) p=0.03` view C.phase — shows: lattice phase locking
- [13 bits] `C kC=1 | C[any] ALWAYS() -> PHASE_COUPLE(K=1.0) p=0.01` view C.phase — shows: lattice phase locking
- [13 bits] `N kN=1 | N[any] ALWAYS() -> PHASE_COUPLE(K=0.3) p=1.0` view N.phase — shows: network synchronization
- [13 bits] `C kC=1 | C[any] ALWAYS() -> PHASE_COUPLE(K=0.3) p=0.1` view C.phase — shows: lattice phase locking
- [13 bits] `N kN=1 | N[any] ALWAYS() -> PHASE_COUPLE(K=0.3) p=0.1` view N.phase — shows: network synchronization
- [14 bits] `C kC=1 | C[any] ALWAYS() -> PHASE_COUPLE(K=1.0) p=0.05` view C.phase — shows: lattice phase locking
- [14 bits] `C kC=1 | C[any] ALWAYS() -> PHASE_COUPLE(K=0.5) p=0.1` view C.phase — shows: lattice phase locking
- [14 bits] `C kC=1 | C[any] ALWAYS() -> PHASE_COUPLE(K=1.0) p=0.02` view C.phase — shows: lattice phase locking
- [14 bits] `C kC=1 | C[any] ALWAYS() -> PHASE_COUPLE(K=0.2) p=1.0` view C.phase — shows: lattice phase locking
- [14 bits] `C kC=1 | C[any] ALWAYS() -> PHASE_COUPLE(K=0.2) p=0.1` view C.phase — shows: lattice phase locking
- [14 bits] `N kN=1 | N[any] ALWAYS() -> PHASE_COUPLE(K=1.0) p=0.02` view N.phase — shows: network synchronization
- [14 bits] `N kN=1 | N[any] ALWAYS() -> PHASE_COUPLE(K=0.2) p=1.0` view N.phase — shows: network synchronization
- [14 bits] `N kN=1 | N[any] ALWAYS() -> PHASE_COUPLE(K=0.2) p=0.1` view N.phase — shows: network synchronization
- [15 bits] `C kC=1 | C[any] ALWAYS() -> PHASE_COUPLE(K=1.0) p=0.1 ; C[any] ALWAYS() -> COPY_RAND() p=1.0` view C.phase — shows: lattice phase locking
- [15 bits] `N kN=3 | N[any] ALWAYS() -> REWIRE_SAME() p=1.0` view N — shows: network fragmentation
- [15 bits] `C kC=1 | C[any] ALWAYS() -> PHASE_COUPLE(K=1.0) p=0.1 ; C[any] ALWAYS() -> COPY_RAND() p=0.1` view C.phase — shows: lattice phase locking
- [15 bits] `C kC=1 | C[any] ALWAYS() -> PHASE_COUPLE(K=0.3) p=0.3` view C.phase — shows: lattice phase locking
- [15 bits] `N kN=3 | N[any] ALWAYS() -> REWIRE_SAME() p=0.1` view N — shows: network fragmentation
- [15 bits] `N kN=4 | N[any] ALWAYS() -> COPY_MAJ() p=1.0` view N — shows: network consensus
- [15 bits] `N kN=4 | N[0] ALWAYS() -> REWIRE_SAME() p=1.0` view N — shows: network fragmentation
- [15 bits] `C kC=1 | C[any] ALWAYS() -> COPY_RAND() p=0.1 ; C[any] ALWAYS() -> PHASE_COUPLE(K=1.0) p=0.1` view C.phase — shows: lattice phase locking
- [15 bits] `C kC=1 | C[any] ALWAYS() -> PHASE_COUPLE(K=0.3) p=0.03` view C.phase — shows: lattice phase locking
- [15 bits] `N kN=2 | N[0] ALWAYS() -> REWIRE_SAME() p=1.0` view N — shows: network fragmentation
- [15 bits] `N kN=4 | N[0] ALWAYS() -> REWIRE_SAME() p=0.1` view N — shows: network fragmentation

## Unnamed behaviour classes -> literature check (61)

| class | bits | shortest program | size |
|---|---|---|---|
| lattice-phase-2 | 11 | `C kC=1 | C[any] ALWAYS() -> PHASE_COUPLE(K=0.1) p=0.1` | 2725 |
| network-phase-1 | 11 | `N kN=1 | N[any] ALWAYS() -> PHASE_COUPLE(K=0.1) p=1.0` | 647 |
| lattice-phase-1 | 11 | `C kC=1 | C[any] ALWAYS() -> PHASE_COUPLE(K=0.1) p=1.0` | 501 |
| lattice-phase-3 | 13 | `C kC=1 | C[any] ALWAYS() -> PHASE_COUPLE(K=0.03) p=1.0` | 285 |
| network-phase-2 | 13 | `N kN=1 | N[any] ALWAYS() -> PHASE_COUPLE(K=1.0) p=0.01` | 268 |
| network-1 | 14 | `N kN=2 | N[any] ALWAYS() -> COPY_RAND() p=0.1` | 79 |
| grid-1 | 14 | `C kC=2 | C[any] ALWAYS() -> COPY_RAND() p=1.0` | 46 |
| network-2 | 15 | `N kN=3 | N[0] ALWAYS() -> COPY_RAND() p=0.1` | 2000 |
| grid-4 | 15 | `C kC=3 | C[0] ALWAYS() -> COPY_RAND() p=1.0` | 1195 |
| grid-5 | 15 | `C kC=3 | C[1] ALWAYS() -> MOVE_EMPTY() p=1.0` | 628 |
| network-7 | 15 | `N kN=4 | N[0] ALWAYS() -> COPY_RAND() p=0.1` | 431 |
| grid-2 | 15 | `C kC=2 | C[0] ALWAYS() -> COPY_RAND() p=0.1` | 353 |
| grid-25 | 15 | `C kC=3 | C[any] ALWAYS() -> COPY_RAND() p=1.0` | 18 |
| network-3 | 15 | `N kN=4 | N[any] DISCORD() -> SET_NEXT() p=1.0` | 12 |
| grid-3 | 15 | `C kC=3 | C[any] ALWAYS() -> COPY_MAJ() p=1.0` | 4 |
| grid-6 | 16 | `C kC=2 | C[any] ALWAYS() -> COPY_RAND() p=0.003` | 203 |
| network-8 | 16 | `N kN=3 | N[1] DISCORD() -> SET(b=0) p=1.0` | 183 |
| network-4 | 16 | `N kN=2 | N[any] ALWAYS() -> REWIRE_SAME() p=0.01` | 130 |
| network-5 | 16 | `N kN=2 | N[0] DISCORD() -> REWIRE_SAME() p=0.1` | 103 |
| grid-7 | 16 | `C kC=2 | C[any] ALWAYS() -> IMITATE_BEST(g=(1.5, 0.5)) p=1.0` | 72 |
| grid-15 | 16 | `C kC=2 | C[any] ALWAYS() -> IMITATE_BEST(g=(1.9, 0.0)) p=0.1` | 27 |
| grid-14 | 17 | `C kC=4 | C[1] ALWAYS() -> ADD_GRAIN() p=0.01` | 496 |
| agents-1 | 17 | `A kA=1 v=1.0 eta=0.3 | A[any] ALWAYS() -> HEAD_TO(b=any,r=1.0) p=1.0` | 438 |
| grid-26 | 17 | `C kC=4 | C[any] ALWAYS() -> SWAP_RAND() p=0.01` | 291 |
| grid-21 | 17 | `C kC=2 | C[any] ALWAYS() -> COPY_RAND() p=0.0003` | 233 |
| grid-13 | 17 | `C kC=3 | C[0] ALWAYS() -> COPY_RAND() p=0.01` | 219 |
| grid-20 | 17 | `C kC=2 | C[any] CNT_LE(b=0,th=3) -> SET_NEXT() p=1.0` | 124 |
| grid-32 | 17 | `C kC=2 | C[any] ALWAYS() -> COPY_MAJ() p=0.0003` | 111 |
| grid-12 | 17 | `C kC=2 | C[0] ALWAYS() -> COPY_RAND() p=0.003` | 85 |
| agents-2 | 17 | `A kA=1 v=0.3 eta=1.0 | A[any] ALWAYS() -> HEAD_MEAN(b=any,r=1.0) p=0.1` | 83 |
| grid-11 | 17 | `C kC=2 | C[any] SAME_LT(th=0.3) -> SET_NEXT() p=1.0` | 83 |
| grid-9 | 17 | `C kC=2 | C[any] CNT_GE(b=0,th=3) -> SET_NEXT() p=1.0` | 72 |
| grid-29 | 17 | `C kC=3 | C[0] ALWAYS() -> COPY_MAJ() p=0.003` | 55 |
| grid-8 | 17 | `C kC=2 | C[any] ALWAYS() -> IMITATE_BEST(g=(1.5, -0.5)) p=1.0` | 39 |
| grid-17 | 17 | `C kC=2 | C[any] NEXT_GE(th=3) -> SWAP_RAND() p=1.0` | 33 |
| grid-22 | 17 | `C kC=2 | C[any] SAME_LT(th=0.3) -> SWAP_RAND() p=1.0` | 24 |
| grid-24 | 17 | `C kC=2 | C[any] CNT_GE(b=0,th=3) -> COPY_MAJ() p=0.1` | 16 |
| grid-10 | 17 | `C kC=2 | C[any] CNT_LE(b=0,th=3) -> COPY_RAND() p=1.0` | 14 |
| network-6 | 17 | `N kN=2 | N[any] CNT_LE(b=0,th=3) -> REWIRE_SAME() p=1.0` | 12 |
| agents-5 | 18 | `A kA=1 v=1.0 eta=0.3 | A[any] ALWAYS() -> HEAD_AWAY(b=any,r=3.0) p=1.0` | 151 |
