# Simulation verification: summary

214 phenomena. Stage 3 candidates: 84.

| class | census | review | qd | all |
|---|---|---|---|---|
| ARTIFACT | 14 | 10 | 13 | 37 |
| DECOMPOSES | 15 | 12 | 39 | 66 |
| COMBINATION EFFECT | 12 | 28 | 27 | 67 |
| SINGLE MECHANISM | 19 | 23 | 2 | 44 |
| NOT RUNNABLE | 0 | 0 | 0 | 0 |

## Flagged (verdict text notes a mismatch or artifact) against class

| class | flagged | not flagged |
|---|---|---|
| ARTIFACT | 31 | 6 |
| DECOMPOSES | 40 | 26 |
| COMBINATION EFFECT | 52 | 15 |
| SINGLE MECHANISM | 17 | 27 |
| NOT RUNNABLE | 0 | 0 |

## Literature verdict against class

| verdict | ARTIFACT | DECOMPOSES | COMBINATION EFFECT | SINGLE MECHANISM |
|---|---|---|---|---|
| COVERED-BY-TWIN | 1 | 0 | 0 | 5 |
| KNOWN | 33 | 65 | 65 | 36 |
| KNOWN-VARIANT | 3 | 1 | 2 | 3 |

Tested features: 1035; robust: 514 (of which 10 by deviation V2).
Card replication (reported only): 913 of 1008 tested features with a card value replicate it.

| class | replicate / with card value |
|---|---|
| ARTIFACT | 72 / 82 |
| DECOMPOSES | 299 / 335 |
| COMBINATION EFFECT | 344 / 363 |
| SINGLE MECHANISM | 198 / 228 |

## Stage 3 candidates by batch and class

| batch | COMBINATION EFFECT | SINGLE MECHANISM (flagged) |
|---|---|---|
| census | 12 | 5 |
| review | 28 | 12 |
| qd | 27 | 0 |
