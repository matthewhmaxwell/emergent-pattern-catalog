# Validation gate — round 9 (multi-label, verified names) — 2026-09-28 05:49

| behaviour | verified | flagged | primary = own behaviour | primary = another (co-present) behaviour | withheld -> primary name (LOCO) |
|---|---|---|---|---|---|
| Schelling segregation | 54 | 54 | 47 | 0 | 0 |
| co-evolving network fragmentation | 38 | 38 | 33 | 1 | 1 |
| cyclic-CA waves/spirals | 40 | 40 | 27 | 0 | 0 |
| flocking | 108 | 107 | 97 | 0 | 0 |
| lattice phase locking | 36 | 36 | 34 | 0 | 0 |
| majority-rule coarsening | 72 | 67 | 59 | 0 | 0 |
| network synchronization | 36 | 36 | 33 | 0 | 0 |
| network voter consensus | 31 | 31 | 25 | 0 | 0 |
| spatial PD chaos (Nowak-May) | 2 | 2 | 0 | 0 | 0 |
| voter coarsening | 56 | 56 | 49 | 0 | 0 |

Status of flagged verified examples: {'NAMED': 405, 'KNOWN-BEHAVIOUR-ATYPICAL-LOOK': 62}
Also-shows pairs (behaviour -> also shows): majority-rule coarsening -> voter coarsening: 67; voter coarsening -> majority-rule coarsening: 56; cyclic-CA waves/spirals -> majority-rule coarsening: 22; cyclic-CA waves/spirals -> voter coarsening: 22; cyclic-CA waves/spirals -> cyclic-CA waves/spirals: 13; Schelling segregation -> majority-rule coarsening: 12; Schelling segregation -> voter coarsening: 12; flocking -> flocking: 10; majority-rule coarsening -> majority-rule coarsening: 8; voter coarsening -> voter coarsening: 7; Schelling segregation -> Schelling segregation: 7; network voter consensus -> network voter consensus: 6

## Criteria

- PASS — F0 fingerprint errors = 0
- PASS — P1 verified examples flagged = 467 / 473 (need >= 90%)
- PASS — P3 names the run does not show = 0 (implementation check)
- PASS — C1 names whose behaviour survives the interaction knock-out = 0 / 467
- PASS — N1 verified trivial/disordered interacting negatives flagged = 0 / 24
- PASS — N2 interaction-free flagged = 0 / 400 (consistency check)
- (reported) primary name = own behaviour: 404 / 467 (87%; target >= 80%)
- (reported) P4 novelty-risk: withheld-class examples still given a primary name = 1 / 467 (named behaviour genuinely co-present in 1)
- negatives dropped as NOT verified trivial/disordered (not scored): ['align at max noise r=1.0 v=1.0', 'align at max noise r=1.0 v=0.3']

### Primary name = another behaviour (verified co-present)
- co-evolving network fragmentation -> network voter consensus: 1

(info) textbook measures passing on UNFLAGGED negatives (never named): 5

### Withheld behaviour -> primary name (LOCO)
- co-evolving network fragmentation -> network voter consensus: 1

**GATE: PASS**
