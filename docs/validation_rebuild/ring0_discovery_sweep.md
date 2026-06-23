# Discovery Sweep — Ring 0 (out-of-catalog emergent systems)

**Seeds/system:** 3. **OOD gate:** match_min_tier=`confirmation`. Candidates are canonical emergent phenomena absent from the 32-pattern catalog; the target read is EMERGENT-UNCLASSIFIED (emergence detected, no false MATCH) with a named firing channel. Null controls must read NO-EMERGENCE.

## Outcome summary

- Candidate verdicts: `{'EMERGENT-UNCLASSIFIED': 16, 'MATCH': 6, 'NO-EMERGENCE': 3}`
- MATCH to a GRADUATED catalog pattern (correct recognition): **3** [('active_nematic', 'P33', 0), ('active_nematic', 'P33', 1), ('active_nematic', 'P33', 2)]
- EMERGENT-UNCLASSIFIED (correct abstention): **16**
- false MATCH (over-claim — instrument finding): **3** [('cahn_hilliard_spinodal', 'P3', 0), ('cahn_hilliard_spinodal', 'P3', 1), ('gray_scott_selfrep', 'P3', 1)]
- NO-EMERGENCE (missed — coverage gap): **3** [('swarmalators', 1), ('swarmalators', 2), ('kauffman_rbn_critical', 1)]
- Null controls: `{'NO-EMERGENCE': 12}` — leaks: []

## Candidates

| system | family | verdict(s) | mean em | channel(s) | closest | ref |
|---|---|---|---|---|---|---|
| dla_fractal | fractal aggregation | EMERGENT-UNCLASSIFIED×3 | 0.822 | synergy(psi_ce:parity), field | P18/screening, P18/none | Witten & Sander 1981, PRL 47:1400 |
| keller_segel_chemotaxis | chemotactic collapse | EMERGENT-UNCLASSIFIED×3 | 1.0 | positions | P11/none | Keller & Segel 1970, J Theor Biol 26:399 |
| active_nematic | nematic / topological defects | MATCH×3 | 0.834 | orientation | P33/definitive | Doostmohammadi et al. 2018, Nat Commun 9:3246 |
| eden_kpz_interface | interface roughening (KPZ) | EMERGENT-UNCLASSIFIED×3 | 0.887 | field | P18/screening | Eden 1961; Kardar-Parisi-Zhang 1986, PRL 56:889 |
| cahn_hilliard_spinodal | conserved phase separation | MATCH×2, EMERGENT-UNCLASSIFIED | 0.861 | field | P3/confirmation, P3/definitive | Cahn & Hilliard 1958, J Chem Phys 28:258 |
| gray_scott_selfrep | self-replicating spots | EMERGENT-UNCLASSIFIED×2, MATCH | 0.714 | field | P3/none, P3/definitive | Pearson 1993, Science 261:189 |
| swarmalators | space-phase coupling (sync+swarm) | NO-EMERGENCE×2, EMERGENT-UNCLASSIFIED | 0.2 | positions, temporal(spectral-peak) | P11/none | O'Keeffe, Hong & Strogatz 2017, Nat Commun 8:1504 |
| langton_ant | emergent order from chaos (highway) | EMERGENT-UNCLASSIFIED | 0.794 | field | P18/none | Langton 1986, Physica D 22:120 |
| kauffman_rbn_critical | edge-of-chaos Boolean dynamics | EMERGENT-UNCLASSIFIED×2, NO-EMERGENCE | 0.565 | vector, temporal(spectral-peak) | P11/none | Kauffman 1969, J Theor Biol 22:437 |

## Null controls

| system | family | verdict(s) | mean em | channel(s) | closest | ref |
|---|---|---|---|---|---|---|
| null_spatial_noise | null | NO-EMERGENCE×3 | 0.065 | positions | P11/none | — |
| null_random_walk | null | NO-EMERGENCE×3 | 0.077 | positions | P11/none | — |
| null_uncoupled_phases | null | NO-EMERGENCE×3 | 0.131 | phases | P11/none | — |
| null_frozen_noise | null | NO-EMERGENCE×3 | 0.124 | positions | P11/none | — |
