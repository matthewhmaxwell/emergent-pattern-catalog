# Planted positive controls: scores (PLANTED_CONTROLS_PLAN.md)

| arm | runs | valid | correct | 95% CI | leaks | invalid | cites target | verdicts |
|---|---|---|---|---|---|---|---|---|
| N_blocked_r1 | 20 | 11 | 0 | 0.00 to 0.28 | 9 | 0 | 9 | {'KNOWN-VARIANT': 10, 'KNOWN': 10} |
| N_blocked_r2 | 20 | 14 | 0 | 0.00 to 0.23 | 6 | 0 | 6 | {'KNOWN-VARIANT': 9, 'KNOWN': 11} |
| N_blocked_r3 | 20 | 12 | 0 | 0.00 to 0.26 | 8 | 0 | 8 | {'KNOWN-VARIANT': 8, 'KNOWN': 12} |
| N_open | 20 | 19 | 0 | 0.00 to 0.18 | 0 | 1 | 19 | {'KNOWN': 19, 'PARSE_FAILED': 1} |
| N_neutral | 20 | 11 | 5 | 0.17 to 0.77 | 8 | 1 | 8 | {'KNOWN-VARIANT': 10, 'KNOWN': 9, 'PARSE_FAILED': 1} |
| N_strictall | 20 | 11 | 1 | 0.00 to 0.41 | 9 | 0 | 9 | {'KNOWN-VARIANT': 9, 'KNOWN': 11} |
| D_open | 10 | 10 | 10 | 0.69 to 1.00 | 0 | 0 | 9 | {'KNOWN': 10} |

Primary (N blocked): 0 of 37 valid runs correct; by card majority over the 3 runs, 0 of 15 cards correct (exact 95% CI 0.00 to 0.22); unanimous across runs: 15.
Pending reading: 0 []
