# Validation gate — round 9 (multi-label, verified names) — 2026-09-28 06:27

| behaviour | verified | flagged | primary = own behaviour | primary = another (co-present) behaviour | withheld -> primary name (LOCO) |
|---|---|---|---|---|---|
| Schelling segregation | 54 | 54 | 50 | 0 | 0 |
| co-evolving network fragmentation | 37 | 37 | 34 | 0 | 1 |
| cyclic-CA waves/spirals | 42 | 42 | 29 | 0 | 0 |
| flocking | 108 | 108 | 98 | 0 | 0 |
| lattice phase locking | 36 | 36 | 34 | 0 | 0 |
| majority-rule coarsening | 71 | 67 | 60 | 0 | 0 |
| network synchronization | 36 | 36 | 34 | 0 | 0 |
| network voter consensus | 28 | 28 | 21 | 0 | 0 |
| spatial PD chaos (Nowak-May) | 3 | 3 | 0 | 0 | 0 |
| voter coarsening | 57 | 57 | 52 | 0 | 0 |

Status of flagged verified examples: {'NAMED': 412, 'KNOWN-BEHAVIOUR-ATYPICAL-LOOK': 56}
Also-shows pairs (behaviour -> also shows): majority-rule coarsening -> voter coarsening: 67; voter coarsening -> majority-rule coarsening: 57; cyclic-CA waves/spirals -> majority-rule coarsening: 24; cyclic-CA waves/spirals -> voter coarsening: 24; cyclic-CA waves/spirals -> cyclic-CA waves/spirals: 13; Schelling segregation -> majority-rule coarsening: 12; Schelling segregation -> voter coarsening: 12; flocking -> flocking: 10; majority-rule coarsening -> majority-rule coarsening: 7; network voter consensus -> network voter consensus: 7; voter coarsening -> voter coarsening: 5; Schelling segregation -> Schelling segregation: 4

## Criteria

- PASS — F0 fingerprint errors = 0
- PASS — P1 verified examples flagged = 468 / 472 (need >= 90%)
- PASS — P3 names the run does not show = 0 (implementation check)
- PASS — C1 names whose behaviour survives the interaction knock-out = 0 / 468
- PASS — N1 verified trivial/disordered interacting negatives flagged = 0 / 312
- PASS — N2 interaction-free flagged = 0 / 400 (consistency check)
- (reported) primary name = own behaviour: 412 / 468 (88%; target >= 80%)
- (reported) P4 novelty-risk: withheld-class examples still given a primary name = 1 / 468 (named behaviour genuinely co-present in 1)
- negatives dropped as NOT verified trivial/disordered (not scored): ['swamped C k=4 C.SAND #2 after', 'swamped C k=4 C.SAND #1 after', 'noise-drowned N copy p=0.01 flip q=0.5', 'noise-drowned N copy p=0.03 flip q=0.5', 'align at max noise r=1.0 v=1.0', 'align at max noise r=1.0 v=0.3']

(info) textbook measures passing on UNFLAGGED negatives (never named): 141

### Withheld behaviour -> primary name (LOCO)
- co-evolving network fragmentation -> network voter consensus: 1

**GATE: PASS**
