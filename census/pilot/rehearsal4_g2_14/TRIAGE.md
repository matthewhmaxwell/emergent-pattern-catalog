# Census triage — 2026-10-02 14:28

564 programs; 73 with a flagged view; errors 0; fingerprint errors 0 (must be 0).
Flagged views: 73 = 23 named + 14 known-behaviour-atypical-look (review) + 36 unnamed in 6 classes (literature check).

## Behaviour map (shortest program per behaviour = the emergence-complexity table)

| behaviour | catalog | bits | shortest program | programs | routes |
|---|---|---|---|---|---|
| network synchronization | P9 | 11 | `N kN=1 | N[any] ALWAYS() -> PHASE_COUPLE(K=1.0) p=1.0` | 14 | 1 |
| lattice phase locking | P9 | 11 | `C kC=1 | C[any] ALWAYS() -> PHASE_COUPLE(K=1.0) p=1.0` | 15 | 1 |
| network consensus | P18 | 14 | `N kN=2 | N[any] ALWAYS() -> COPY_RAND() p=1.0` | 3 | 2 |
| domain coarsening | P18 | 14 | `C kC=2 | C[any] ALWAYS() -> COPY_RAND() p=0.1` | 3 | 2 |
| network fragmentation | P34 | 14 | `N kN=2 | N[any] ALWAYS() -> REWIRE_SAME() p=1.0` | 2 | 1 |

## Routes (every distinct mechanism that produces each behaviour)

### network synchronization — 1 routes; ingredient common to all: N:ALWAYS>PHASE_COUPLE

| layers | rule templates | programs | bits | shortest example |
|---|---|---|---|---|
| N | N:ALWAYS>PHASE_COUPLE | 14 | 11 | `N kN=1 | N[any] ALWAYS() -> PHASE_COUPLE(K=1.0) p=1.0` |

### lattice phase locking — 1 routes; ingredient common to all: C:ALWAYS>PHASE_COUPLE

| layers | rule templates | programs | bits | shortest example |
|---|---|---|---|---|
| C | C:ALWAYS>PHASE_COUPLE | 15 | 11 | `C kC=1 | C[any] ALWAYS() -> PHASE_COUPLE(K=1.0) p=1.0` |

### network consensus — 2 routes; ingredient common to all: none

| layers | rule templates | programs | bits | shortest example |
|---|---|---|---|---|
| N | N:ALWAYS>COPY_MAJ | 2 | 14 | `N kN=2 | N[any] ALWAYS() -> COPY_MAJ() p=1.0` |
| N | N:ALWAYS>COPY_RAND | 1 | 14 | `N kN=2 | N[any] ALWAYS() -> COPY_RAND() p=1.0` |

### domain coarsening — 2 routes; ingredient common to all: none

| layers | rule templates | programs | bits | shortest example |
|---|---|---|---|---|
| C | C:ALWAYS>COPY_MAJ | 2 | 14 | `C kC=2 | C[any] ALWAYS() -> COPY_MAJ() p=1.0` |
| C | C:ALWAYS>COPY_RAND | 1 | 14 | `C kC=2 | C[any] ALWAYS() -> COPY_RAND() p=0.1` |

### network fragmentation — 1 routes; ingredient common to all: N:ALWAYS>REWIRE_SAME

| layers | rule templates | programs | bits | shortest example |
|---|---|---|---|---|
| N | N:ALWAYS>REWIRE_SAME | 2 | 14 | `N kN=2 | N[any] ALWAYS() -> REWIRE_SAME() p=1.0` |

## Review — known behaviour, atypical look (14)

- [11 bits] `C kC=1 | C[any] ALWAYS() -> PHASE_COUPLE(K=1.0) p=0.1` view C.phase — shows: lattice phase locking
- [13 bits] `N kN=1 | N[any] ALWAYS() -> PHASE_COUPLE(K=0.3) p=1.0` view N.phase — shows: network synchronization
- [13 bits] `N kN=1 | N[any] ALWAYS() -> PHASE_COUPLE(K=0.3) p=0.1` view N.phase — shows: network synchronization
- [13 bits] `C kC=1 | C[any] ALWAYS() -> PHASE_COUPLE(K=1.0) p=0.03` view C.phase — shows: lattice phase locking
- [13 bits] `C kC=1 | C[any] ALWAYS() -> PHASE_COUPLE(K=1.0) p=0.01` view C.phase — shows: lattice phase locking
- [13 bits] `C kC=1 | C[any] ALWAYS() -> PHASE_COUPLE(K=0.3) p=0.1` view C.phase — shows: lattice phase locking
- [14 bits] `C kC=1 | C[any] ALWAYS() -> PHASE_COUPLE(K=0.2) p=1.0` view C.phase — shows: lattice phase locking
- [14 bits] `N kN=1 | N[any] ALWAYS() -> PHASE_COUPLE(K=0.2) p=1.0` view N.phase — shows: network synchronization
- [14 bits] `C kC=1 | C[any] ALWAYS() -> PHASE_COUPLE(K=0.2) p=0.1` view C.phase — shows: lattice phase locking
- [14 bits] `N kN=1 | N[any] ALWAYS() -> PHASE_COUPLE(K=0.2) p=0.1` view N.phase — shows: network synchronization
- [14 bits] `C kC=1 | C[any] ALWAYS() -> PHASE_COUPLE(K=1.0) p=0.05` view C.phase — shows: lattice phase locking
- [14 bits] `C kC=1 | C[any] ALWAYS() -> PHASE_COUPLE(K=1.0) p=0.02` view C.phase — shows: lattice phase locking
- [14 bits] `N kN=1 | N[any] ALWAYS() -> PHASE_COUPLE(K=1.0) p=0.02` view N.phase — shows: network synchronization
- [14 bits] `C kC=1 | C[any] ALWAYS() -> PHASE_COUPLE(K=0.5) p=0.1` view C.phase — shows: lattice phase locking

## Unnamed behaviour classes -> literature check (6)

| class | bits | shortest program | size |
|---|---|---|---|
| lattice-phase-2 | 11 | `C kC=1 | C[any] ALWAYS() -> PHASE_COUPLE(K=0.1) p=0.1` | 23 |
| network-phase-2 | 11 | `N kN=1 | N[any] ALWAYS() -> PHASE_COUPLE(K=0.1) p=0.1` | 8 |
| lattice-phase-1 | 11 | `C kC=1 | C[any] ALWAYS() -> PHASE_COUPLE(K=0.1) p=1.0` | 2 |
| network-phase-1 | 11 | `N kN=1 | N[any] ALWAYS() -> PHASE_COUPLE(K=0.1) p=1.0` | 1 |
| grid-1 | 14 | `C kC=2 | C[any] ALWAYS() -> COPY_RAND() p=1.0` | 1 |
| network-1 | 14 | `N kN=2 | N[any] ALWAYS() -> COPY_RAND() p=0.1` | 1 |
