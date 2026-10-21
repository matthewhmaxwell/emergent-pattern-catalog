# EPC Emergence Census — deviations from PREREG_census_v1.md

Every change to code, thresholds, library or procedure after registration is logged here with its date, the reason,
and whether it was made before or after seeing the affected data. Catalog expansions (§4.3) are planned procedure,
not deviations, but each one is logged with its gate and placebo results.

Registered: 2026-10-11 (OSF, embargoed). Frozen code: commit 52cee0f. Plan hash (SHA-256):
7a398029068d04e14d4be3253269799aee906017460d36b6d74900db72ab8427

| # | date | change | reason | before / after seeing affected data |
|---|---|---|---|---|
| 1 | 2026-10-21 | `census/sim.py`: the G2 `ADD_GRAIN` op now records avalanche **durations** (`dur += top.any(...)`, `self.aval_dur.append(dur)`), exactly as the G1 `_C_SAND` op already did; then re-ran only the affected programs on their original seeds. | During the run, all 1,206 G2 sandpile (`ADD_GRAIN`) programs crashed (`ValueError: need at least one array to stack` at sim.py:517): the op recorded avalanche sizes but never durations, so stacking the empty duration list threw. The ≤14-bit rehearsal contained no such program, so the gate never hit it. Without the fix these 1,206 programs are absent from the G2 map (its avalanche routes); G1's avalanche map was complete and unaffected (0 errors in all 289,798). | **After** seeing the affected data (the 1,206 errors). The change is confined to the `ADD_GRAIN` branch, which is G2-only; verified result-identical on 16 already-run programs (12 G2 non-ADD_GRAIN + 4 G1 → 0 changed) and that the 1,206 then run with 0 errors. No other program's result changes. The novelty arm ran entirely under the frozen code 52cee0f; it is not re-run. Fixed code committed after 52cee0f on branch `emergence-census`. |
