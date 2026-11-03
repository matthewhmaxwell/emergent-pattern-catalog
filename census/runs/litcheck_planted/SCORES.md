# Planted positive controls: scores (PLANTED_CONTROLS_PLAN.md)

| arm | runs | valid | correct (read) | 95% CI | close calls | correct (field) | leaks | invalid | cites target | verdicts |
|---|---|---|---|---|---|---|---|---|---|---|
| N_blocked_r1 | 20 | 11 | 5 | 0.17 to 0.77 | 8 | 0 | 9 | 0 | 9 | {'KNOWN-VARIANT': 10, 'KNOWN': 10} |
| N_blocked_r2 | 20 | 14 | 5 | 0.13 to 0.65 | 5 | 0 | 6 | 0 | 6 | {'KNOWN-VARIANT': 9, 'KNOWN': 11} |
| N_blocked_r3 | 20 | 12 | 3 | 0.05 to 0.57 | 4 | 0 | 8 | 0 | 8 | {'KNOWN-VARIANT': 8, 'KNOWN': 12} |
| N_open | 20 | 19 | 0 | 0.00 to 0.18 | 0 | 0 | 0 | 1 | 19 | {'KNOWN': 19, 'PARSE_FAILED': 1} |
| N_neutral | 20 | 11 | 7 | 0.31 to 0.89 | 4 | 0 | 8 | 1 | 8 | {'KNOWN-VARIANT': 10, 'KNOWN': 9, 'PARSE_FAILED': 1} |
| N_strictall | 20 | 11 | 4 | 0.11 to 0.69 | 6 | 1 | 9 | 0 | 9 | {'KNOWN-VARIANT': 9, 'KNOWN': 11} |
| D_open | 10 | 10 | 10 | 0.69 to 1.00 | 0 | 0 | 0 | 0 | 9 | {'KNOWN': 10} |

Primary (N blocked): 13 of 37 valid runs correct; by card majority over the 3 runs, 4 of 15 cards correct (exact 95% CI 0.08 to 0.55); same score in every valid run: 10.
Superseded field rule: 0 of 15 cards. Agreement over the 11 cards with three valid runs: same score 6, same verdict 7.
Pending reading: 0 []
