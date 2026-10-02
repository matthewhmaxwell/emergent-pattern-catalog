# Placebo test — 2026-10-02 10:51

1510 programs (0 errors), 2444 interaction-free views compared with an independent interaction-free run of the same program. Any flag is false by construction.

**False flags: 23 / 2444**

| view family | views | FALSE FLAGS | screen emergent in >= 2 of 3 seeds | per seed-run: generic >= 0.5 | model-free complexity | consensus gain | avalanche |
|---|---|---|---|---|---|---|---|
| agents | 399 | 0 | 18 | 119 (5.0%) | 16 (0.7%) | 0 (0.0%) | 0 (0.0%) |
| field | 912 | 0 | 598 | 3592 (65.6%) | 680 (12.4%) | 0 (0.0%) | 0 (0.0%) |
| grid | 641 | 0 | 3 | 32 (0.8%) | 0 (0.0%) | 0 (0.0%) | 0 (0.0%) |
| lattice-phase | 118 | 12 | 48 | 0 (0.0%) | 280 (39.5%) | 0 (0.0%) | 0 (0.0%) |
| network | 291 | 0 | 9 | 52 (3.0%) | 3 (0.2%) | 52 (3.0%) | 0 (0.0%) |
| network-phase | 83 | 11 | 36 | 0 (0.0%) | 185 (37.1%) | 0 (0.0%) | 0 (0.0%) |

## False flags

- `CF kC=1 fields=2 D=[0.05, 0.05] | C.KURA(K=0.002)` view C.phase — null run: `CF kC=1 fields=2 D=[0.05, 0.05] | C.KURA(K=0.0)`
- `C kC=4 | C.MAJ(p=0.5) ; C.KURA(K=1.0)` view C.phase — null run: `C kC=4 | C.KURA(K=0.0)`
- `CA kC=1 kA=2 v=0.1 eta=0.3 | C.KURA(K=0.03)` view C.phase — null run: `CA kC=1 kA=2 v=0.1 eta=0.3 | C.KURA(K=0.0)`
- `C kC=3 | C.COPY(p=0.001) ; C.KURA(K=1.0)` view C.phase — null run: `C kC=3 | C.KURA(K=0.0)`
- `C kC=1 | C.COPY(p=0.01) ; C.KURA(K=0.01) ; C.COPY(p=1.0)` view C.phase — null run: `C kC=1 | C.KURA(K=0.0)`
- `C kC=2 | C.KURA(K=0.02) ; C.CYCLE(th=3,p=1.0)` view C.phase — null run: `C kC=2 | C.KURA(K=0.0)`
- `C kC=1 | C.KURA(K=0.002)` view C.phase — null run: `C kC=1 | C.KURA(K=0.0)`
- `N kN=1 | N.KURA(K=1.0) ; N.COPY(p=1.0)` view N.phase — null run: `N kN=1 | N.KURA(K=0.0)`
- `N kN=1 | N.KURA(K=0.1) ; N.COPY(p=0.1)` view N.phase — null run: `N kN=1 | N.KURA(K=0.0)`
- `N kN=4 | N.KURA(K=0.3)` view N.phase — null run: `N kN=4 | N.KURA(K=0.0)`
- `C kC=1 | C.KURA(K=0.1) ; C.COPY(p=1.0)` view C.phase — null run: `C kC=1 | C.KURA(K=0.0)`
- `N kN=3 | N.KURA(K=0.05)` view N.phase — null run: `N kN=3 | N.KURA(K=0.0)`
- `N kN=1 | N.KURA(K=1.0)` view N.phase — null run: `N kN=1 | N.KURA(K=0.0)`
- `N kN=2 | N.KURA(K=0.1)` view N.phase — null run: `N kN=2 | N.KURA(K=0.0)`
- `CF kC=1 fields=1 D=[0.02] | C.KURA(K=0.1)` view C.phase — null run: `CF kC=1 fields=1 D=[0.02] | C.KURA(K=0.0)`
- `C kC=3 | C.KURA(K=0.5)` view C.phase — null run: `C kC=3 | C.KURA(K=0.0)`
- `C kC=3 | C.KURA(K=0.002)` view C.phase — null run: `C kC=3 | C.KURA(K=0.0)`
- `N kN=4 | N.KURA(K=0.5)` view N.phase — null run: `N kN=4 | N.KURA(K=0.0)`
- `N kN=1 | N.KURA(K=0.1)` view N.phase — null run: `N kN=1 | N.KURA(K=0.0)`
- `C kC=2 | C.KURA(K=0.1)` view C.phase — null run: `C kC=2 | C.KURA(K=0.0)`
- `N kN=2 | N.KURA(K=0.3)` view N.phase — null run: `N kN=2 | N.KURA(K=0.0)`
- `N kN=3 | N.KURA(K=0.01)` view N.phase — null run: `N kN=3 | N.KURA(K=0.0)`
- `N kN=2 | N.KURA(K=0.01)` view N.phase — null run: `N kN=2 | N.KURA(K=0.0)`

(1279 s; runs: census/pilot/speed22, census/pilot/rehearsal_g1_14; salt '')
