# Validation gate — round 9 (multi-label, verified names) — 2026-10-02 11:51

| behaviour | verified | flagged | primary = own behaviour | primary = another (co-present) behaviour | withheld -> primary name (LOCO) |
|---|---|---|---|---|---|
| Schelling segregation | 54 | 54 | 45 | 0 | 0 |
| co-evolving network fragmentation | 37 | 37 | 33 | 1 | 1 |
| cyclic-CA waves/spirals | 41 | 41 | 26 | 0 | 0 |
| flocking | 108 | 108 | 98 | 0 | 0 |
| lattice phase locking | 36 | 36 | 33 | 0 | 0 |
| majority-rule coarsening | 72 | 63 | 56 | 0 | 0 |
| network synchronization | 36 | 36 | 33 | 0 | 0 |
| network voter consensus | 30 | 30 | 20 | 0 | 0 |
| spatial PD chaos (Nowak-May) | 2 | 2 | 0 | 0 | 0 |
| voter coarsening | 58 | 58 | 52 | 0 | 0 |

BEHAVIOUR level (classes sharing a textbook measure merged; presentation only): primary names the example's own behaviour in 396 / 465 (85%).
Also-shows at behaviour level (own behaviour -> other behaviour): cyclic waves and spirals -> domain coarsening: 22; segregation with vacancies -> domain coarsening: 12

Status of flagged verified examples: {'KNOWN-BEHAVIOUR-ATYPICAL-LOOK': 68, 'NAMED': 397}
Also-shows pairs (behaviour -> also shows): majority-rule coarsening -> voter coarsening: 63; voter coarsening -> majority-rule coarsening: 58; cyclic-CA waves/spirals -> majority-rule coarsening: 22; cyclic-CA waves/spirals -> voter coarsening: 22; cyclic-CA waves/spirals -> cyclic-CA waves/spirals: 15; Schelling segregation -> majority-rule coarsening: 12; Schelling segregation -> voter coarsening: 12; network voter consensus -> network voter consensus: 10; flocking -> flocking: 10; Schelling segregation -> Schelling segregation: 9; majority-rule coarsening -> majority-rule coarsening: 7; voter coarsening -> voter coarsening: 6

## Criteria

- PASS — F0 fingerprint errors = 0
- PASS — P1 verified examples flagged = 465 / 474 (need >= 90%)
- PASS — P3 names the run does not show = 0 (implementation check)
- PASS — C1 names whose behaviour survives the interaction knock-out = 0 / 465
- PASS — N1 verified trivial/disordered interacting negatives flagged = 0 / 311
- PASS — N2 interaction-free flagged = 0 / 400 (consistency check)
- (reported) primary name = own behaviour: 396 / 465 (85%; target >= 80%)
- (reported) P4 novelty-risk: withheld-class examples still given a primary name = 1 / 465 (named behaviour genuinely co-present in 1)
- negatives dropped as NOT verified trivial/disordered (not scored): ['swamped C k=4 C.SAND #1 after', 'swamped C k=4 C.SAND #2 after', 'noise-drowned N copy p=0.003 flip q=0.5', 'noise-drowned N copy p=0.01 flip q=0.5', 'noise-drowned N copy p=0.03 flip q=0.5', 'align at max noise r=1.0 v=1.0', 'align at max noise r=1.0 v=0.3']

### Primary name = another behaviour (verified co-present)
- co-evolving network fragmentation -> network voter consensus: 1

(info) textbook measures passing on UNFLAGGED negatives (never named): 141

### Withheld behaviour -> primary name (LOCO)
- co-evolving network fragmentation -> network voter consensus: 1

**GATE: PASS**
