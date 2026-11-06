# Rule-by-rule review of plug-ins N02, N07, N19 (agent report, 2026-11-03)

Agent "Review plug-ins N02 N07 N19" finished

I found no MUST FIX in any of the three modules. Each one implements its paper's stated equations and update rules correctly. N02's fidelity check will fail on the c = 0.49 leg, but that comes from an inconsistency inside the paper, not from the module. I edited no files. My test scripts are in `/private/tmp/claude-501/-Users-matthewmaxwell-Downloads/ed5d6e8e-67aa-43e5-b2b6-30638ce75e8d/scratchpad/rv/`, all run at small sizes with seeds other than the census seeds.

## N02 (`census/planted_models/n02.py`): OK with notes

**What the paper's model really is**, and why c = 0.49 cannot come out uniform:

- **The model is Eqs. 10-11 (repeated as 112-113), and the module implements it exactly.** I checked `rhs_hat` against an analytic evaluation of Eqs. 10-11 on a smooth field: error 2e-13. The chemotaxis signs and the divergence form are right. I also re-derived the paper's Eq. 14/15 from Eqs. 10-11, and they match, so the paper's own linear theory uses this sign convention.
- **The sign convention is as written.** The prose after Eq. 9 says ξ &gt; 0 makes the current run up the gradient, but Eq. 8 makes it run down. The paper's analysis and its single-coupling simulations (Figs. 9-16) both follow Eqs. 8-11. Flipping the signs would contradict Figs. 9-16, and by Eq. 15 it would also remove the pattern for every c above about 0.08, including the sharp spots at 0.17.
- **The paper's own theory says the uniform state is unstable at every c &gt; 0 on this branch.** With c1 = c, c2 = -c, Eq. 22 becomes Det(k) = 100 − (9+30c)k² + (0.1+c²)k⁴. Its instability condition reduces to 41 + 540c + 500c² &gt; 0, which always holds. At c = 0.49 the growth rate is 0.79 per time unit. The unstable band (k² from 4.5 to 65) exists under every reading of the domain (a 2π box with any grid, or L = 128 at spacing 1) and also for the D1 = 0.1 row. The paper's statement that k_c "ultimately vanishes" is only asymptotic (k_c² ≈ 15/c).
- **The paper's mean-fixing does not explain it.** Pinning only the k = 0 mode leaves every other growth rate unchanged. I tested it at 64×64: pinned c = 0.49 still gives a pattern (amplitude ≈ 25), and pinned c = 0.17 blew up at t ≈ 58, which confirms the author's report.
- **The authors' simulation probably did not integrate Eqs. 10-11.** Appendix A's Eqs. 115-117 are linearised yet called "nonlinear PDEs". They have ∇²v where ∇²u belongs in the ψ chemotaxis term, −2ub/(b²+a) instead of −2b²u/(b²+a), and implicit factors 1+Δt k² and 1+Δt D k² with no D1. The paper also says its simulated pattern spacing does not depend on ξ, against its own Eq. 14. In the module the spacing shifts as Eq. 14 predicts (peak k² moves from 36 to 20). The most direct candidate I tried, chemotaxis acting only on the fluctuations, blows up at 0.49 rather than going uniform. The authors' actual code cannot be reconstructed from the paper. Version 2 changes none of this.

**Findings**

1. **Should fix** (latent bug, `_Field.substep`, line 143). The preconditioner's cross terms use the full K², while `rhs_hat` (lines 113-114) zeroes the Nyquist first derivatives. On the Nyquist lines the two disagree, and that produces spurious checkerboard growth whenever c1 &lt; 0 with c2 ≠ 0: 27 per time unit at dt = 0.01 and 10 at dt = 0.005, where the true rate is −81. A flipped run at c = 0.49 reached amplitude 6.6 in one time unit, and `advance()` accepted it because the fields stayed finite and under 1e6. No planned run triggers it: headline, knock-out, sub-programs and fidelity all have c1 ≥ 0 and c2 ≤ 0. For those, the scheme matches the exact linear growth rates at 64 and 128, and a census-size headline run had zero power near the Nyquist frequency. Fix: compute `Kc = self.iKX.imag**2 + self.iKY.imag**2` once and use `n12 = dt*c1*ub*Kc`, `n21 = dt*c2*vb*Kc`.
2. **Note** (fidelity, third leg). It will fail for the reason above, not because of the module. Keep the check frozen, and log a note before the freeze, as was done for P1, so that "implementation not confirmed" is read correctly. Dropping the mean-fixing is the right call: it does not change the 0.49 result and it makes the headline blow up.
3. **Note** (plan level). At D2/D1 = 1000 the pattern is a Turing pattern that exists without chemotaxis. The knock-out removes the reactions too, so it kills that pattern whatever chemotaxis does. The paper's new behaviour (amplitude rising then falling with c, spots turning to stripes) only shows across a sweep, never in one run at c = 0.17.
4. **Note** (minor):
   - `_run` line 179: after a failed step it records the restored start-of-step state, whereas census `sim.run` records the blown-up state.
   - `fidelity` line 213: the S = 1 fallback does not re-estimate the cost, and there is no deadline guard.
5. **OK**:
   - Parameters match the plan row.
   - Grid spacing 2π/128 is well supported, because Appendix B's stripe count does not change with grid shape.
   - One time unit per census step; output format passes `check_out`; RULE_TEXT is clean.
   - Census-size headline run takes 63 s for 3 seeds, amplitude ≈ 13.1, steady, presence = True.

## N07 (`census/planted_models/n07.py`): OK with notes

1. **OK** (rules, Sec. II Eqs. 2-8). The sender keeps λ1 + ξ1(1−λ1). The receiver pools what it gets with its own mass and keeps λ2 + ξ2(1−λ2) of the total. All four left-turn targets match Eqs. 2 and 8.
2. **OK** (parameters). λ1 = 0.4 and λ2 = 0.6 are the fractions kept, as in Figs. 4-6. The Fig. 2 case with λ = 0 only makes sense if λ is the fraction kept.
3. **OK** (knock-out). Rule 2 off sets λ2 = 1, which makes the odd diffusion and the off-diagonal mobility zero. That gives the paper's delta-correlated achiral model (Sec. VI.3). Measured achiral variance is 0.669, against 2/3 in theory. The knock-out also lowers D from 0.14 to 0.075; that is inherent to the model.
4. **OK** (batching and time). The batched update is bit-identical to a plain sequential loop: 0.0 difference at W = 12 (three chunk sizes, both versions) and at W = 32. Mass is conserved exactly. One sweep is one unit of the paper's time. The steady variance is 0.416 against 0.415 from `theory_C`, which itself matches Eqs. 51-52, 57, 60 and 67-68. It is very fast (0.2 s for 200 sweeps at 64×64 with 3 seeds).
5. **Note** (`presence`). At census size the verdict is borderline. One run gave η = 3.48 ± 0.21, just under the 3.5 threshold, so presence = False, while the amplitude with η fixed at 4 was 0.117 ± 0.004 (the paper's value is 0.11). The achiral run correctly shows nothing. I suggest also reporting the fixed-η amplitude and a fit against `theory_C(64)`; the criterion itself stays as is.
6. **Note**:
   - The fit uses the cos 4θ ring coefficient rather than Fig. 6's axis correlations. That is documented and sound.
   - The uniform start m = 1 is a choice, since the paper does not give one.
   - RULE_TEXT contains the common title words "mass" and "correlations", but both come from the frozen card.

## N19 (`census/planted_models/n19.py`, 05:42 version): OK with notes

1. **OK** (Eq. 1a-c). The cell flux (Scharfetter-Gummel scheme) has the correct form and sign. Bump tests confirm b moves up both the attractant and the nutrient gradient, with sensitivity κ b/(1+b²). The discrete operator agrees with Eq. 1a to second order (relative error 5e-4). Cell number is conserved to 1e-16 and b stays ≥ 0. The attractant is produced at rate b and decays implicitly; the nutrient is consumed at rate β b n and resupplied at rate β; diffusion ratios are 1, 2.5 and 0.75. I re-derived Eq. 2 from the module's Jacobian.
2. **OK** (point K, read independently). I rendered page 4 at 300 dpi and measured the ticks and marker myself: κa = 11.10 ± 0.19, β = 0.0839 ± 0.003. That agrees with the module. The axis calibration holds: along β = 0.084, Eq. 2 puts the region boundaries at 7.5, 10.1 and 14.6, and the figure's curves cross at 7.4, 10.0-10.2 and 14.7. K lies in the mixed region, with max growth 0.21 at q = 0.59, as in Fig. 3K. The Fig. 2 caption values are stable (max −1.0e-3 for q ≥ 2π/50). Fig. 1's caption uses δa = 2; with that value the caption point is only marginally unstable (+4e-4), and of the travelling-wave type, not the mixed region. So P1 is justified.
3. **Should fix** (`presence`). At census size, point K gives two or three aggregates in moving bands. "At least one aggregate" holds in 100% of late frame pairs, but "exactly one" holds in only 0-18%, so presence reports absent. The census domain is 64 against the paper's 50. Add a `present_n_agg_ge_1` field so a loss caused by scale is not read as absence of the behaviour.
4. **Note**:
   - The P1 log says β ≈ 0.083, while the module uses 0.084 and docstring line 80 describes a run at 0.083; these should be made to agree.
   - `linear_growth` reports the slow q = 0.01 mode for stable cases; scan q ≥ 2π/L instead.
   - `advance` line 178 silently cuts a census step short if it exceeds 200,000 sub-steps; it should flag the run instead.
   - Grid spacing dx = 1 is a stated assumption.
5. **OK**:
   - Knock-out and sub-programs remove exactly rules 1 and 2.
   - One time unit per census step; output passes `check_out`.
   - RULE_TEXT changes only the size and the P1 numbers.
   - A census-size run takes 63 s for 3 seeds.
