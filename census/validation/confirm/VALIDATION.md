# Validation gate — round 9 (multi-label, verified names) — 2026-09-28 05:49

| behaviour | verified | flagged | primary = own behaviour | primary = another (co-present) behaviour | withheld -> primary name (LOCO) |
|---|---|---|---|---|---|
| Schelling segregation | 54 | 54 | 47 | 0 | 0 |
| co-evolving network fragmentation | 37 | 37 | 34 | 0 | 0 |
| cyclic-CA waves/spirals | 41 | 41 | 28 | 0 | 0 |
| flocking | 108 | 107 | 95 | 0 | 0 |
| lattice phase locking | 36 | 36 | 34 | 0 | 0 |
| majority-rule coarsening | 72 | 70 | 59 | 0 | 0 |
| network synchronization | 36 | 36 | 34 | 0 | 0 |
| network voter consensus | 33 | 33 | 26 | 0 | 0 |
| spatial PD chaos (Nowak-May) | 3 | 3 | 0 | 0 | 0 |
| voter coarsening | 60 | 60 | 52 | 0 | 0 |

Status of flagged verified examples: {'KNOWN-BEHAVIOUR-ATYPICAL-LOOK': 68, 'NAMED': 409}
Also-shows pairs (behaviour -> also shows): majority-rule coarsening -> voter coarsening: 70; voter coarsening -> majority-rule coarsening: 60; cyclic-CA waves/spirals -> majority-rule coarsening: 24; cyclic-CA waves/spirals -> voter coarsening: 24; cyclic-CA waves/spirals -> cyclic-CA waves/spirals: 13; Schelling segregation -> majority-rule coarsening: 12; Schelling segregation -> voter coarsening: 12; flocking -> flocking: 12; majority-rule coarsening -> majority-rule coarsening: 11; voter coarsening -> voter coarsening: 8; Schelling segregation -> Schelling segregation: 7; network voter consensus -> network voter consensus: 7

## Criteria

- PASS — F0 fingerprint errors = 0
- PASS — P1 verified examples flagged = 477 / 480 (need >= 90%)
- PASS — P3 names the run does not show = 0 (implementation check)
- PASS — C1 names whose behaviour survives the interaction knock-out = 0 / 477
- PASS — N1 verified trivial/disordered interacting negatives flagged = 0 / 24
- PASS — N2 interaction-free flagged = 0 / 400 (consistency check)
- (reported) primary name = own behaviour: 409 / 477 (86%; target >= 80%)
- (reported) P4 novelty-risk: withheld-class examples still given a primary name = 0 / 477 (named behaviour genuinely co-present in 0)
- negatives dropped as NOT verified trivial/disordered (not scored): ['align at max noise r=1.0 v=1.0', 'align at max noise r=1.0 v=0.3']

(info) textbook measures passing on UNFLAGGED negatives (never named): 8

**GATE: PASS**
