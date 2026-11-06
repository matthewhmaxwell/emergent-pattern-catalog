# Rule-by-rule review of plug-ins N01, N12, N13 (agent report, 2026-11-03)

Agent "Review plug-ins N01 N12 N13" finished

I edited nothing and ran no pipeline code. The three modules follow their papers rule by rule, and I found no equation, sign or update-order errors. None needs a MUST FIX. Each has a few should-fix items, and one N13 item needs your decision before the freeze.

**Limits on what I could check:** I read the papers through a web-fetch tool that returns summaries and short quotes, not full text. I could not view the figures themselves, so I couldn't directly confirm the N13 marker positions read off Fig. 6.

## N01 (2609.01668), pulse-coupled oscillators: OK with notes

What I confirmed:
- **Rules match Sec. 2–3.** Charging is dx/dt = 1 − γx at γ = 0, with threshold 1 and reset 0. The pulse is k/N, initial voltages are Uniform(0,1), and absorption is single-pass with no avalanche. That matches the paper's recurrence "b_j′ = b_j − b_f + N + s_f", where the pulse is the firer's own cluster size.
- **Exact counts reproduce eq. (4).** Enumeration gives 227/256, 3124/3125 and 42275/46656.
- **The simulate() path agrees with the fidelity path.** `_evolve` and `_batch_outcome` gave 0 mismatches on all words for N = 4 and 5, and on 4,000 sampled words for N = 6.
- **Other checks pass.** A reduced fidelity run passed (Monte Carlo z = −0.17, −1.5, −0.15). Knock-out and `keep` handling, sizes 200/800, step length and `check_out` (census and area4) are all fine.

Findings:
1. **Should fix: RULE_TEXT (line 67) leaves N undefined.** The edit "N identical" → "200 identical" removed the definition, but rule (3) still says "kicks every other voltage up by k/N". Fix: "World: N = 200 identical oscillators …".
2. **Note: RULE_TEXT (line 73) changes the measurement clause.** "record whether all N end as one group" became "record the voltages", which goes beyond size, length or setting. N12 and N13 keep their "Measure …" clauses. Dropping the sizes, many starts and enumeration is justified; for the rest, a smaller edit would be "…from one random start and record whether all 200 end as one group."
3. **Should fix: `fidelity()` (lines 189–216) never runs `_evolve`, the code simulate() uses.** Both the enumeration and the Monte Carlo go through the separate `_batch_outcome`. I showed they agree, but the recorded check should cover the census code path. Fix: also count the enumeration through `_evolve`; it takes under 1 s for N ≤ 6.
4. **Note: the standard error at line 206 uses the enumerated value instead of the paper's A/N^N.** This is harmless while the two are equal.
5. **Note, a consequence of the plan, not a module bug: both runs look frozen.** With one census step equal to the free period, the knock-out frames all equal frame 0 exactly. All three real runs were fully synchronised and static from step 5. So the phase view compares two static configurations; this is worth recording before the freeze.

## N12 (2609.39267, with its Ref. [42] = arXiv 2411.01811): OK with notes

What I confirmed:
- **Couplings and boundaries match Ref. [42].** Couplings are c/r^4 with "∑_j J_ij = 4", each pair counted once, using minimum image on the torus. That is the model the paper says it uses for Fig. 4 ("uniform model with J = b J_LR", β_c = 0.7315 from Ref. [42]).
- **The normalisation matches the paper's results.** My Wolff run at L = 128 gave M² = 0.5194 ± 0.0024, against the Table 3 fit value of 0.5187. At L = 64 it gave 0.540 against 0.539.
- **The dynamics are exact.**
  - The block-wise FFT heat-bath matches a naive exact sequential heat-bath on identical random draws to 1e-14.
  - Heat-bath and Wolff agree at L = 16 (z = −0.4).
  - The von Mises conditional, the Wolff bond probability and the reflection formula are correct.
- **Other checks pass.** Knock-out at J = 0 gives uniform angles each sweep. One step is one random-order sweep, sizes are 64/128, and the format is valid. A census run takes about 30 s and reaches M² ≈ 0.51–0.54 within about 200 sweeps.

Findings:
1. **Should fix: fidelity criterion (c) (lines 274–282) is dominated by noise.** It fits η_ℓ from ln ln L with L0 fixed at 0.21, and the lever arm is tiny: ln ln(L/0.21) moves only 0.19 between L = 128 and 512. In my reduced run (sizes 64 and 128), (b) passed with a deviation of only −0.005, yet (c) gave η = 0.44 and failed. The fallback to L ≥ 64 also uses sizes below the paper's L_min of 128. Fix: pass if |η_fit − 0.34| ≤ max(0.07, 3σ_fit), or drop (c) and keep (b), which already tests the paper's measure point by point.
2. **Note: docstring (lines 9–14) could spell out the mapping.** In the notation of eqs. 27–29, with the sum over r &gt; 1, the module corresponds to J_LR = 1/(1+b). This is justified because β_c is quoted from Ref. [42], and confirmed by the M² match. Also, the Table 2–3 captions don't name the model; it is inferred from Fig. 4 and Sec. IV.3.
3. **Note: the normalisation was cross-checked at L = 128 and 256 before the tolerance was written.** The model choice comes from Ref. [42]'s stated normalisation, and the 0.01 tolerance is 10× looser than the match found, so this is not tuning.
4. **Note: the RULE_TEXT edits are needed and accurate.** "cluster" became "single-spin", and "from random angles" and "1000 sweeps" were added. No leak terms.

## N13 (2609.39858), noisy voter model with shocks: OK with notes (one item needs your decision)

What I confirmed:
- **Rules match the paper.**
  - Shocks follow eq. (10): with probability r per step, a shock replaces that step's update.
  - Noise is a flip, as in eq. (27); copying is from another agent, as in eqs. (3)–(4); time steps are Δt = 1/N.
  - The δ-shock and the start are both at m_0 = m_r = 0.3 (Fig. 6, N = 500).
- **The agent code matches the exact chain.** It agrees with the exact stationary distribution with copying on and off (KS 0.002 and 0.001). The per-sweep variance confirms N updates per sweep (90.9 against 91.5 expected).
- **Other checks pass.** Knock-out and `keep` handling and the format (census and area4) are fine. A reduced fidelity run passed (KS 0.007).
- **The figure-read values are robust.** η = 10^-3.5 and r = 10^-4.5 come from marker d) in Fig. 6; the text gives no numbers. In the exact chain at N = 500 the point sits inside a trimodal region about two decades wide in both r and η, so a half-decade reading error would not move it out. The module's marker b) at (1e-2, 1e-4) comes out as a single peak at m_r, as "shock-dominated" should.

Findings:
1. **Should fix, needs your decision before the freeze: r is held fixed per step while η is scaled with N** (HEADLINE line 69, `_params` lines 86–87). The paper's own stationary solutions (appendix eqs. 34–35, e.g. "Nr/(2+Nr)", with r a rate per unit time) depend on the shock rate only through N·r per unit time, which is N² × the per-step r. Holding the per-step r fixed moves the setting away from the paper's (shares are of total probability; heights are relative to the tallest peak):

   | Setting | N | Mass at consensus | Middle-peak height |
   |---|---|---|---|
   | Paper | 500 | 0.036 | 0.175 |
   | Census, module's fixed r | 200 | 0.24 | 0.03 |
   | Census, N²·r held at the paper's value | 200 | 0.042 | 0.37 |
   | Area4, module's fixed r | 800 | 0.005 | 1.0 (the −1 consensus peak drops to 0.055) |
   | Area4, N²·r held at the paper's value | 800 | 0.033 | 0.12 |

   Under the module's choice, each census seed also sees only about 6 shocks, against about 40 with N²·r held fixed. The verification stage's 4×-node condition then runs a visibly different shape. Fix: r = 10^-4.5 × (500/N)², i.e. 1.98e-4 at the census size and 1.2e-5 at area4, with RULE_TEXT updated to match. The plan's wording ("parameters the paper defines per population size") arguably supports either reading, which is why it's your call.
2. **Should fix: fidelity part (a) (lines 228–230, `_state_peaks`) almost always passes.** For any η &lt; 1/N and r &gt; 0, shock injection creates a cusp at n_r, so three state-level maxima appear across almost the whole grid I scanned. Part (b), with binned prominence and KS, carries the check. Fix: apply (b)'s binned-prominence criterion to the exact distribution in (a).
3. **Note: the module's peak rule treats consensus spikes as always prominent.** As a result, the paper's asymmetric bimodal regimes e) and f) never appear in the exact chain under that rule. The headline is trimodal under any reasonable rule.
4. **Note: the RULE_TEXT edits are setting-only and accurate.** "with probability r = 0.0000316 per step" correctly replaces "at rate". No coined terms, surnames or arXiv id in any of the three texts; the only title words present ("oscillators", "opinion") were already in the frozen cards.

All modules are in /Users/matthewmaxwell/VPS-staging/projects/epc-emergence-census/census/planted_models/ (n01.py, n12.py, n13.py).
