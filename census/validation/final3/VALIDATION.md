# Validation gate — round 9 (multi-label, verified names) — 2026-10-02 10:25

| behaviour | verified | flagged | primary = own behaviour | primary = another (co-present) behaviour | withheld -> primary name (LOCO) |
|---|---|---|---|---|---|
| Schelling segregation | 54 | 54 | 47 | 0 | 0 |
| co-evolving network fragmentation | 37 | 37 | 34 | 0 | 1 |
| cyclic-CA waves/spirals | 41 | 41 | 27 | 0 | 0 |
| flocking | 108 | 108 | 96 | 0 | 0 |
| lattice phase locking | 36 | 36 | 34 | 0 | 0 |
| majority-rule coarsening | 72 | 67 | 61 | 0 | 0 |
| network synchronization | 36 | 36 | 35 | 0 | 0 |
| network voter consensus | 31 | 31 | 27 | 0 | 0 |
| spatial PD chaos (Nowak-May) | 3 | 3 | 0 | 0 | 0 |
| voter coarsening | 58 | 58 | 52 | 0 | 0 |

BEHAVIOUR level (classes sharing a textbook measure merged; presentation only): primary names the example's own behaviour in 413 / 471 (88%).
Also-shows at behaviour level (own behaviour -> other behaviour): cyclic waves and spirals -> domain coarsening: 23; segregation with vacancies -> domain coarsening: 12; network fragmentation -> network consensus: 1

Status of flagged verified examples: {'KNOWN-BEHAVIOUR-ATYPICAL-LOOK': 58, 'NAMED': 413}
Also-shows pairs (behaviour -> also shows): majority-rule coarsening -> voter coarsening: 67; voter coarsening -> majority-rule coarsening: 58; cyclic-CA waves/spirals -> majority-rule coarsening: 23; cyclic-CA waves/spirals -> voter coarsening: 23; cyclic-CA waves/spirals -> cyclic-CA waves/spirals: 14; Schelling segregation -> majority-rule coarsening: 12; Schelling segregation -> voter coarsening: 12; flocking -> flocking: 12; Schelling segregation -> Schelling segregation: 7; voter coarsening -> voter coarsening: 6; majority-rule coarsening -> majority-rule coarsening: 6; network voter consensus -> network voter consensus: 4

## Criteria

- PASS — F0 fingerprint errors = 0
- PASS — P1 verified examples flagged = 471 / 476 (need >= 90%)
- PASS — P3 names the run does not show = 0 (implementation check)
- PASS — C1 names whose behaviour survives the interaction knock-out = 0 / 471
- FAIL — N1 verified trivial/disordered interacting negatives flagged = 1 / 311
- PASS — N2 interaction-free flagged = 0 / 400 (consistency check)
- (reported) primary name = own behaviour: 413 / 471 (88%; target >= 80%)
- (reported) P4 novelty-risk: withheld-class examples still given a primary name = 1 / 471 (named behaviour genuinely co-present in 1)
- negatives dropped as NOT verified trivial/disordered (not scored): ['swamped C k=4 C.SAND #1 after', 'swamped C k=4 C.SAND #2 after', 'noise-drowned N copy p=0.003 flip q=0.5', 'noise-drowned N copy p=0.01 flip q=0.5', 'noise-drowned N copy p=0.03 flip q=0.5', 'align at max noise r=1.0 v=1.0', 'align at max noise r=1.0 v=0.3']

(info) textbook measures passing on UNFLAGGED negatives (never named): 141

### Withheld behaviour -> primary name (LOCO)
- co-evolving network fragmentation -> network voter consensus: 1

### Negatives flagged
- lattice oscillators k=2 K=0.002 ['C.phase']

**GATE: FAIL**
