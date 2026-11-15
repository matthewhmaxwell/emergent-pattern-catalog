# Literature-reproduction specs, group C (P26, P12, P32, P17, P28)

Status: FROZEN PLAN, written 2026-11-15 before any simulation. No model code was run to write this file. The only computations were closed-form formulas from the sources (Kramers rate, Eq. 2.8, Eq. 2.7a, Eq. 61) and pixel digitization of published figures.

Standard: a reproduction counts only when a specific published number, or an analytic result stated in the source, is compared with our model output under the tolerance fixed below. A diagnostic is reported but never turns a FAIL into a PASS. If a primary check fails, the result is recorded as FAIL. Parameters, estimators and tolerances are not changed after the run.

Quotes are kept under 15 words (shorter than the 40-word limit in the brief) because of a copyright constraint on quoting. Each target therefore also gives an exact locator (equation, figure, page) so a reviewer can check it. Equations and numbers are reproduced in full.

Figure reads come from pixel digitization of the source PDF: frame and tick positions were detected, then marker centroids. The reading uncertainty for each value is stated with it.

## Summary

| Pattern | Primary target | Verdict | Model change needed | Rough CPU |
|---|---|---|---|---|
| P26 stochastic resonance | RMP 1998 Fig. 7: Floquet spectral amplification η(D) at Ω=0.1, A=0.1. Peak η = 12.23 at D = 0.150. Two-state cross-check: Eq. 2.8 gives D_SR = 0.148 | **Target found, feasible** | Parameters only (a=b=1, A=0.1, f=Ω/2π), plus a vectorized driver with no rule change | < 1 CPU-h |
| P12 cyclic dominance | Nature 2007 Fig. 2b: P_ext(M) at t=N for L=40 and L=100 (figure reads). M_c = (4.5±0.5)×10⁻⁴ is a reading of those curves | **Target found, feasible with a Numba kernel; moderate risk (time-unit and directed-selection convention)** | Run length must be converted to rate time units. An accelerated inner loop must be verified bit-identical | ~8–10 CPU-h with Numba (infeasible in pure Python) |
| P32 specialization | None applies to our model. The 1996 paper is closed-access here, and our model is neither its fixed-threshold model nor the 1998 reinforcement model | **Drop** for the current model. Weak alternative: TBD 1998 Fig. 2a, but only with a new model implementation | A new model class is needed | n/a |
| P17 collective sensing | None. The authors state that the simulation curve's quantitative values depend on model parameters | **Drop** (keep as qualitative) | n/a | n/a |
| P28 wealth condensation | Boghosian et al. 2017 Eq. (61): steady-state Gini with redistribution, G∞ = 1/(1+2.3682 τ̄^0.5943). Secondary: Boghosian 2014 slope of 1/a(t) = 2.04046×10⁻⁵ | **Target found, feasible** | Parameters only (χ, redistribute_every, f), plus a calling convention | ~4 CPU-h |

Record corrections found while doing this work. These are for the catalog owners and do not gate anything.
1. P12: `rps_spatial.py` and REPLICATION_NOTES give the arXiv ID of Reichenbach, Mobilia & Frey 2007 as `q-bio/0702032`. That ID is a different paper (Danchin, on bacteria and Lamarckism). The correct ID is **arXiv:0709.0217**. The companion JTB 2008 paper is arXiv:0801.1798.
2. P32: the catalog cites Bonabeau, Theraulaz & Deneubourg 1996 (the **fixed**-threshold model) for a **threshold-reinforcement** model. Reinforcement was introduced in Theraulaz, Bonabeau & Deneubourg 1998, Proc. R. Soc. B 265:327–332, doi:10.1098/rspb.1998.0299. The P32 panel positive also uses hand-picked "specializing" seeds (43–46, 49).
3. P28: paper §4.18 claims the runs reproduce "published Gini trajectories" from Boghosian 2014. arXiv:1212.6300v2 never mentions the Gini coefficient (zero occurrences in the full text). Its published simulation quantity is the 1/a(t) fit (see P28-B).
4. P17: the docstring of `collective_sensing.py` says noise averages as σ/√N through centroid attraction. Berdahl et al. argue the opposite: the gain comes from emergent speed differentials across the group, not from pooling estimates, and they say the key determinant is the group's spatial extent.

---

## P26: stochastic resonance in a periodically forced bistable well

### 1. Published target

**Source.** Gammaitoni, L., Hänggi, P., Jung, P. & Marchesoni, F. (1998). Stochastic resonance. *Rev. Mod. Phys.* 70, 223–287. doi:10.1103/RevModPhys.70.223. The figure is credited "After Jung and Hänggi (1991a)", which is *Phys. Rev. A* 44, 8032, per the RMP reference list. Copy read: Augsburg OPUS repository PDF.

**Primary target: RMP Fig. 7, p. 236.** It plots the "spectral amplification η versus the noise intensity D at a fixed modulation frequency" Ω = 0.1, for A = 0.1, 0.3, 0.5, 0.8, plus the linear-response curve. The curves were computed with nonadiabatic Floquet theory for the Fokker–Planck equation (Fig. 6 caption, p. 235; method in Sec. IV.A). They are numerically exact for the continuous model, not a two-state approximation. The model (Eqs. 4.5–4.8) is ẋ = x − x³ + A cos Ωt + ξ with ⟨ξ(t)ξ(0)⟩ = 2Dδ(t) (Eqs. 2.1–2.3), so ΔV = 1/4 and x_m = 1. The amplification is η = (2|M₁|/A)² (Eq. 4.21), where M₁ is the Ω Fourier component of the asymptotic mean ⟨x(t)⟩_as.

Digitized values for the A = 0.1 curve:

| Quantity | Read value | Reading uncertainty |
|---|---|---|
| η_max | **12.23** | ±0.15 |
| D at maximum (flat top spans 0.145–0.155) | **0.150** | ±0.008 |
| η(D = 0.20) | **10.8** | ±0.3 (overlaps the dotted LRT line) |
| η(D = 0.30) | **6.5** | ±0.3 |
| For reference: A = 0.3 peak | 9.3 at D ≈ 0.125 | ±0.2 / ±0.01 |
| For reference: LRT (dotted) peak | ≈13.1 at D ≈ 0.145 | ±0.2 |

**Analytic cross-checks.** These come from formulas stated in the source and were evaluated here at Ω = 0.1 and ΔV = 1/4.
- Two-state theory, Eq. (2.8): 4r_K²(D_SR) = Ω²(ΔV/D_SR − 1), with r_K = (√2π)⁻¹ exp(−ΔV/D) (Eq. 2.4). This gives **D_SR = 0.148**, matching the Floquet peak.
- Eq. (2.7a) with the exact stationary ⟨x²⟩₀ (computed by quadrature) peaks at η = 13.2, D = 0.144. That matches the dotted LRT curve, which suggests the figure reading is sound.
- Eq. (2.15): the SNR maximum is at D̄_SR = ΔV/2 = **0.125**. The source states this "if the prefactor of the Kramers rate is independent of D". It also notes (p. 239, Fig. 8) that in the full model the SNR diverges as D → 0 because of intrawell motion. This value is therefore **reported only, not gated**.

### 2. Model match (`epc/models/stochastic_resonance.py`, `BistableDoubleWell`)

| Item | Source (Fig. 7) | Ours | Fix |
|---|---|---|---|
| Equation | ẋ = ax − bx³ + A cos Ωt + √(2D)ξ | Same SDE with A sin(2πft). Euler–Maruyama | Same model. The phase is irrelevant to η |
| a, b | 1, 1 (dimensionless) | canonical 4, 1 | **Parameter**: a = b = 1. The canonical frame maps exactly: x = 2x̃, t = t̃/4, D = 16D̃, A = 8Ã, Ω = 4Ω̃, η = η̃/16 |
| Amplitude | A = 0.1 | canonical 1.0 (equivalent to Ã = 0.125, which matches no published curve) | **Parameter**: A = 0.1 |
| Frequency | Ω = 0.1 | f = 0.005 | **Parameter**: f = 0.1/2π = 0.0159155 |
| Solver | Floquet (exact) | Euler–Maruyama, dt = 0.01 | dt check (see Protocol) |
| Output | asymptotic ⟨x(t)⟩ Fourier amplitude | `simulate()` writes one dict per step per trial | **Code addition, no rule change**: a vectorized driver that integrates all trials together and stores only Fourier projections. At 3.8×10⁵ steps × 400 trials × 15 D, `simulate()` would create about 2×10⁹ dicts |
| Initial condition | asymptotic state | alternating wells | Burn-in removes it |

Equivalence test for the new driver, to run before the production run: feed the driver and `simulate()` the same pre-drawn normal increments for 4 trials × 2000 steps at D = 0.15. Trajectories must agree to ≤1e-12.

### 3. Protocol

- Frame: a = b = 1, A = 0.1, Ω = 0.1, dt = 0.01.
- D grid: {0.06, 0.08, 0.10, 0.12, 0.13, 0.14, 0.15, 0.16, 0.17, 0.18, 0.20, 0.25, 0.30, 0.40, 0.50}.
- 10 seeds (1001–1010) per D, each with 400 independent trials.
- Burn-in: 10 drive periods (628.3 time units). That is at least 4 relaxation times 1/(2r_K) at D = 0.06 and more than 12 at D ≥ 0.08.
- Measurement: 50 whole periods (3141.6 time units, 314,159 steps).
- Estimator, per trial: c_k = (2/T)Σ x(t)e^{−iΩt}Δt over the measurement window. Per seed: M̂ = mean_k c_k and s² = sample variance of the c_k. Bias-corrected η̂ = (|M̂|² − s²/n)/A². Report the mean ± SE across the 10 seeds.
- Peak location: fit a quadratic to the seed-mean η̂(D) on D ∈ [0.12, 0.18] (7 points). D_max is the vertex and η_max the vertex value. Get the SE from 1000 bootstrap resamples of seeds.
- dt check, required: repeat D ∈ {0.10, 0.15, 0.30} at dt = 0.005 with the same seeds. If any |Δη|/η exceeds 2%, the gated values are taken from dt = 0.005, rerun over the full D window [0.12, 0.18]. This rule is fixed now.
- Frame check, required: repeat D ∈ {0.10, 0.15, 0.30} in the canonical frame (a=4, b=1, A=0.8, f=0.063662, D×16, dt=0.0025) and compare 16·η. It must agree with the a=b=1 result within 2 combined SE. If not, record a FAIL of the scaling and investigate before gating.
- Reported, not gated: the SNR local maximum (background estimated from the residual periodogram in a ±5-bin annulus around Ω) against ΔV/2 = 0.125.

### 4. Tolerance and pass rule (fixed)

PASS only if all four hold:
1. η_max ∈ 12.23 × (1 ± 0.05) = [11.62, 12.84].
2. |D_max − 0.150| ≤ 0.015.
3. η(0.20) ∈ 10.8 × (1 ± 0.08).
4. η(0.30) ∈ 6.5 × (1 ± 0.08).

Justification:
- Error sources for η_max: the figure read (1.2%), our statistical SE (about 1.5% per point, from 4000 trials with a per-trial Fourier SD of about 0.11 against x̄ ≈ 0.35), and residual Euler bias bounded by the dt check (≤ 2%). Their root sum of squares is about 3%.
- ±5% is wide enough for those errors and still excludes the A = 0.3 curve (9.3) and the bare two-state x_m formula (18.7).
- The D tolerance covers the reading error (±0.008) plus the bootstrap SE (~0.005). It still excludes the A = 0.3 peak (0.125).
- The two curve points have larger reading error because the A = 0.1 curve and the dotted LRT curve overlap there, hence ±8%.

### 5. Compute and risks

- Compute: 3.8×10⁵ steps per (D, seed), each step a numpy update of a 400-element array (about 15–20 µs), so 6–8 s per (D, seed). 150 main runs take about 20 min. The dt and frame checks add about 15 min. Under 1 CPU-h in total. Run it on the VPS.
- Risks not caused by the model:
  - Digitization of a 1998 journal figure (mitigated by two independent analytic cross-checks: Eq. 2.8 and Eq. 2.7a).
  - The Floquet curve for A = 0.1 includes nonlinear response, so a linear-response solver must not be substituted.
  - The canonical P26 panel parameters are not the parameters used here. A pass validates the code and the SR physics, not the panel's detector settings.

---

## P12: cyclic dominance, spatial rock-paper-scissors

### 1. Published target

**Source.** Reichenbach, T., Mobilia, M. & Frey, E. (2007). Mobility promotes and jeopardizes biodiversity in rock–paper–scissors games. *Nature* 448, 1046–1049. doi:10.1038/nature06095. arXiv:0709.0217. Companion: Reichenbach, Mobilia & Frey (2008) *J. Theor. Biol.* 254, 368–383, arXiv:0801.1798.

**What M_c is.** The Fig. 2 caption (arXiv p. 12) says the transition "sharpens at a critical mobility Mc ≈ (4.5 ± 0.5) × 10−4". M_c is a reading of where P_ext curves for growing N saturate at 1, at σ = µ = 1. The paper gives no finite-size-scaling estimator for it. In the digitized panel, the drawn M_c line sits at 5.0×10⁻⁴. The reproducible quantities are therefore the **P_ext(M) curves themselves**.
- Definition (Fig. 2b caption and SI "Extensivity"): P_ext is the probability that, starting from random initial conditions, the system has reached an absorbing state (one species left) after waiting time t = N, where N = L² sites.
- Methods: 500–2000 realizations per point, L from 20 to 200.
- The SI also gives λ_c = 0.8 ± 0.05 in system-length units. That value is not used here (see the four failed wavelength attempts).

**Fig. 2b digitized values.** M grid = 10⁻⁴ × 2^{k/2}. Marker spacing is 0.1506 decades, which is √2. Reading uncertainty ±0.02, covering marker centroid and overlap.

| M | 7.07e-5 | 1.00e-4 | 1.41e-4 | 2.00e-4 | 2.83e-4 | 4.00e-4 | 5.66e-4 |
|---|---|---|---|---|---|---|---|
| L = 40 (magenta) | 0.427 | 0.599 | 0.776 | 0.908 | 0.945 | 0.976 | 0.987 |
| L = 100 (blue) | 0.247 | 0.407 | 0.507 | 0.611 | 0.734 | 0.911 | 0.959 |
| L = 200 (black, optional) | 0.104 | 0.167 | 0.210 | 0.309 | 0.534 | 0.830 | 0.964 |

### 2. Model match (`epc/models/rps_spatial.py`)

| Item | Source | Ours | Fix |
|---|---|---|---|
| Lattice | square, periodic, von Neumann | same | none |
| States, reactions | A, B, C, ∅; selection σ, reproduction µ, exchange/hop ε | same | none |
| Mobility | M = 2ε/N | ε = M·N/2 | none |
| Initial condition | i.i.d., each of 4 states w.p. 1/4 (JTB p. 4, at σ = µ) | `init_mode='random'` (uniform over 0..3) | none |
| Absorbing state | one species left | `run()` stops when ≤1 species remains | none |
| Reaction choice per step | pair chosen; reaction w.p. ∝ rate; done "if the move is allowed" (JTB p. 4) | ordered pair (site, neighbor); selection only if site dominates neighbor; reproduction only into an empty neighbor | **Convention risk**, see below |
| **Time unit** | rates set the time unit (σ = µ = 1). JTB Eq. (49), p. 18: transition rate W = Nσac for CA→C. JTB p. 4: D = ε/(2N) | one "generation" = N elementary steps with probabilities normalized by σ+µ+ε | **Protocol conversion** (no code change). For our directed scheme, matching both W = Nσac and D = ε/(2N) requires N(σ+µ+ε) elementary steps per unit time, so **t = N corresponds to G = N(σ+µ+ε) = N(2 + MN/2) of our generations**. At L = 100 that is 23,536–48,285 generations, not 10⁴ |
| Speed | 500 realizations at L = 100 | pure-Python inner loop (~27 ms per generation at L = 100, per notes) | **Code addition**: a Numba kernel that applies the same pre-drawn numpy arrays (`site_rows`, `site_cols`, `nbr_idx`, `reaction_roll`) in the same order, so results are bit-identical. Equivalence test: 3 seeds × 50 generations at L = 30, M ∈ {1e-5, 1e-3}; grids must be identical at every generation |

**Convention risk.** It is stated now, not discovered later.
- The Nature Methods text defines the time unit as one generation in which every individual has reacted on average once. Read literally, that is G = N, which conflicts with fixing time by σ = µ = 1. The JTB rate formula and the diffusion constant fix the conversion above for a directed scheme.
- If the original code applied selection and reproduction regardless of pair order, but still defined time by W = Nσac, it is equivalent to our code at **M/2**.

These alternatives are run as **diagnostics only**:
- (D1) G = N.
- (D2) our code at M/2 with G = N(2 + (M/2)N/2).

They are reported next to the primary result and cannot turn a primary FAIL into a PASS.

### 3. Protocol

- σ = µ = 1, von Neumann neighborhood, random initial condition.
- Primary sizes L = 40 and L = 100. Optional L = 200 at M ∈ {2.83e-4, 4.00e-4} only, with 200 realizations each, run only if the budget allows.
- M grid: the 7 values in the table above.
- 500 independent realizations per (L, M). Seeds: base 120000 + 1000·index(L) + 10·index(M) + r for realization r, so streams never overlap.
- Run each realization until absorption or until G_max = N(2 + MN/2) generations. Record a binary outcome: absorbed by G_max or not. Keep the absorption generation for the extinction-time histogram (reported).
- Use `record_every` larger than G_max so only the initial state is kept. Memory stays O(L²).
- Estimator: P̂_ext = k/500 with a Wilson 95% CI.
- Finite size: comparisons are made at the source's own L values, so no extrapolation is needed. Report ordering checks between L = 40 and L = 100.

### 4. Tolerance and pass rule (fixed)

PASS only if all of the following hold:
1. At **L = 100**, every one of the 7 points satisfies |P̂ − P_fig| ≤ 0.10, and the mean signed deviation over the 7 points is within ±0.05.
2. The same two conditions hold at **L = 40**.
3. Sharpening: P̂(L=100) < P̂(L=40) at each M ≤ 2.83e-4, where the figure gap is at least 0.17. This is a sign test at each of the 5 points.

Justification:
- Binomial SE is ≤ 0.022 for our 500 realizations and ≤ 0.022 for the source's (at least 500). Figure reading adds ±0.02. The SD of the difference is about 0.035, so ±0.10 is about 3σ per point.
- The mean-deviation bound catches a systematic shift, such as a time-unit error. Diagnostics D1 and D2 are expected to fail this bound if the primary convention is right.
- The "M_c = 4.5e-4" statement is treated as **reproduced** only if the L = 100 and optional L = 200 points at M ≥ 4e-4 are within tolerance. It is not fitted separately.

### 5. Compute and risks

- With Numba (RNG pre-draw dominates, ~30 ns per elementary step):
  - L = 40: 5–6×10⁶ steps per run, about 0.2 s. 3500 runs take about 12 min.
  - L = 100: 2.4–4.8×10⁸ steps per run, 7–15 s at full length. Many high-M runs stop early. About 8 CPU-h for 3500 runs.
  - L = 200 (optional): 1.2–1.6×10¹⁰ steps per run, 6–8 min. About 30–45 CPU-h for 400 runs.
  - Diagnostics D1 and D2 at L = 100 add about 5 CPU-h.
- Pure Python is about 100× slower: around 800 CPU-h at L = 100, which is **not feasible**.
- Run it on the VPS with process-level parallelism.
- Risks not caused by the model:
  - The convention ambiguity above.
  - Figure digitization.
  - The source uses an E-exchanges-then-reaction Gillespie variant (JTB Eq. 3). It is statistically equivalent to per-step reaction selection but not identical in event ordering.
  - Numba becomes a new dependency.

---

## P32: emergent specialization, response-threshold division of labor

### 1. Published target

**No clean published number applies to our model.**
- **Bonabeau, Theraulaz & Deneubourg (1996)**, *Proc. R. Soc. B* 263, 1565–1569, doi:10.1098/rspb.1996.0229.
  - Access: the publisher returned 403, there is no Europe PMC or PubMed record, and Semantic Scholar marks it CLOSED. It could not be read here, so no quote or number from it is given.
  - Content: per the 1998 paper's summary, it is a **fixed**-threshold model. Castes with different, constant thresholds respond with probability s²/(s²+θ²), and the paper compares the model with Wilson's (1984) *Pheidole* caste-ratio data.
  - A full-text copy is needed (institutional access) before it can anchor anything. It still would not match our model.
- **Theraulaz, Bonabeau & Deneubourg (1998)**, *Proc. R. Soc. B* 265, 327–332, doi:10.1098/rspb.1998.0299 (PMC1688885; IAS archive copy read) is the actual threshold-reinforcement source.
  - It gives quantitative results for **its own** continuous-time model, Eqs. (1)–(7).
  - Fig. 2a (p. 330): numbers of specialists N₁ and N₂ against the forgetting rate φ, with N = 100, m = 2, φ + ξ = 11, α = 3, δ = 1, p = 0.2, σ = 0.1, θ(0) = 500, x(0) = 0.1, specialist meaning θ < 100.
  - Text, p. 330: all individuals are specialists below φ ≈ 0.4, and "For φ > 2, forgetting is quick and no specialization is observed". The comparison symbols are glyph-encoded in the PDF and were checked against the rendered page.
  - Fig. 2a digitized values (±3): N_j ≈ 100 for φ ≤ 0.4, 96 at 0.5, 85–87 at 0.6, 70–72 at 1.0, 52–56 at 1.4–1.6, 46–51 at 1.8–1.9, and 0 for φ ≥ 2.0.
  - Fig. 2b: all specialists for p < 0.04 and again for p > 0.42.

### 2. Model match (`epc/models/division_of_labor.py`, `ResponseThresholdModel`)

The response function s²/(s²+θ²) matches. Nothing else does:

| Item | TBD 1998 | BTD 1996 (as summarized in 1998) | Ours |
|---|---|---|---|
| Thresholds | reinforced, ξ (learn) and φ (forget) per unit time, bounded in [θ_min, θ_max] (Fig. 1c: [1, 1000]) | fixed, caste-specific | reinforced ±0.05/±0.01 per step, clipped to [0.01, 1] |
| Activity | continuous time-fraction x_ij, Eq. (5), quitting rate p, noise σ | engage w.p. T, quit w.p. p | one discrete task per step chosen by argmax; no persistence or quitting |
| Stimulus | ∂s/∂t = δ − (α/N)Σx, unbounded | similar | s += 0.1 − 0.1·(workers), clipped to [0, 1] |
| Castes | identical start | two castes | identical start |

No parameter setting maps our model onto either source. Changing the code means writing a new model.

### 3–5. Recommendation

**Drop P32 from the quantitative-reproduction set for the current model.** Label it qualitative and fix the citation (Record correction 2).

If a quantitative anchor is wanted, the weak alternative is:
- Implement TBD 1998 Eqs. (1)–(7) as a new class.
- Run 10 seeds at φ ∈ {0.2, 0.3, 0.6, 1.0, 1.5, 2.2, 2.5, 3.0} with the Fig. 2a parameters.
- Gate on: N₁ = N₂ = 100 at φ ≤ 0.3 in ≥ 9/10 seeds; N₁ = N₂ = 0 at φ ≥ 2.2 in all seeds; mean N_j at φ = 1.0 within 71 ± 12.

It is weak for three reasons:
1. The source does not state the integration step, the noise discretization, the θ bounds used in Fig. 2, or the convergence cap.
2. The figure comes from single runs with no uncertainty.
3. It would validate a new implementation, not the catalog's model.

Compute would be trivial (minutes).

---

## P17: collective sensing (Berdahl et al. 2013)

### 1. Published target

**Source.** Berdahl, A., Torney, C. J., Ioannou, C. C., Faria, J. J. & Couzin, I. D. (2013). Emergent sensing of complex environments by mobile animal groups. *Science* 339, 574–576. doi:10.1126/science.1225883. Konstanz KOPS copy read. The supplement could not be retrieved (Science returned 403).

**No clean published number for the simulation model.**
- The model output appears only as the red line in Fig. 1B. Its performance Ψ rises from about 1.17 at N = 1 to about 1.5 at N = 256 (coarse read, ±0.03).
- The caption states: "The quantitative values of this curve are dependent on model parameters" (Fig. 1 caption, p. 575). It adds that only the qualitative trend is robust.
- The other numbers in the paper are experimental: fish data, and Ψ_max = 3.2 ± 0.5 for the experimental arena.

### 2. Model match (`epc/models/collective_sensing.py`)

The models differ in every component that sets the numbers:
- Environment: the source uses projected moving noisy light fields; ours is a static Gaussian peak.
- Interactions: the source uses local zonal interactions (refs 17–18, Couzin-type); ours is global attraction to the centroid by heading torque.
- Metric: the source uses Ψ, the time-averaged darkness normalized by a null value; ours is a chemotactic index computed from centroid distance.
- Mechanism claim: the source says the key determinant is group spatial extent relative to the environment length scale; our docstring claims σ/√N averaging.

Matching would mean adopting their released code and light fields (icouzin.princeton.edu). That reproduces their code, not ours, and the authors themselves disclaim the curve's values.

### 3–5. Recommendation

**Drop P17 from the quantitative-reproduction set and keep it as qualitative.** The trend can be stated, but the anchor should not be labeled "Quantitative" as the catalog currently does.

There is a possible non-P17 anchor, not recommended. The single-agent kinesis baseline predicts a stationary density ρ(x) ∝ 1/v_eff(x) for speed-modulated active walkers in the slow-variation limit (Schnitzer 1993, *PRE* 48:2553; Tailleur & Cates 2008, *PRL* 100:218103; neither retrieved in this pass). Our persistence length (~v/D_r ≈ 6) is comparable to field_sigma (4–5), so the limit is not met. It also tests only N = 1, not the emergent group effect.

---

## P28: wealth condensation, yard-sale model

### 1. Published targets

**P28-A (primary): steady-state Gini with redistribution.**
- Source: Boghosian, B. M., Devitt-Lee, A., Johnson, M., Marcq, J. A. & Wang, H. (2017). Oligarchy as a phase transition: the effect of wealth-attained advantage in a Fokker–Planck description of asset exchange. *Physica A* 476, 15–37. arXiv:1511.00770v2 read.
- Sec. III.D, Eq. (61), p. 14: G∞ ≈ 1/(1 + 2.3682 τ̄^0.5943). The paper calls it "the solid curve is the remarkably accurate fit to these points". The points are numerical steady-state solutions of the redistributive Fokker–Planck equation, Eq. (33), at W/N = 1 (Fig. 4, τ̄ ∈ (0, 1.5]).
- Definitions:
  - Eq. (8): transaction Δw = √(γΔt)·min(w, x)·η with η = ±1, so the stake fraction is f = √(γΔt).
  - Sec. III.A: a wealth tax at rate τ per unit time, returned equally, so the per-transaction-time drift is τΔt(W/N − w).
  - τ̄ = τ/γ is the "redistribution per transaction" (Sec. III.B).
- No error bar is given. In Fig. 4 the points lie on the fit to within about 0.005.

| τ̄ | 0.2 | 0.5 | 1.0 |
|---|---|---|---|
| G∞ from Eq. (61) | **0.5236** | **0.3893** | **0.2969** |
| For reference, G∞ at 2τ̄ | 0.4213 | 0.2969 | 0.2186 |

**P28-B (secondary, pure YSM, closest to the canonical positive).**
- Source: Boghosian, B. M. (2014). Kinetics of wealth and the Pareto law. *Phys. Rev. E* 89, 042804 (arXiv:1212.6300v2). Simulation section, pp. 10–12, Fig. 5.
- Setup: N = 5×10⁴ agents with 100 units each, fixed stake fraction α = 0.25 with a fair coin (η(β) from Eq. 22), up to 10⁹ transactions.
- Result: 1/a(t) ≈ 3.93264 + 0.0000204046 t, by least squares over t ≤ 8×10⁷. Here a(t) comes from w_min(t) through Eq. (37), and w_min(t) is the point where the straight-line fit to the log-linear part of A(w, t) reaches A = 1.
- The paper notes the slope is close to 1/N.

### 2. Model match (`epc/models/yard_sale.py`)

| Item | Source | Ours | Fix |
|---|---|---|---|
| Transaction | random pair, stake f·min(w_i, w_j), fair coin | identical (λ_save = 0) | none |
| Time mapping | t̄ = γt, one transaction per agent per Δt | pairs; each agent is in 2 transactions per N transactions | **Protocol**: t̄ = 2f²·n_tx/N |
| Redistribution | continuous OU tax τ(W/N − w) | every R transactions: w ← (1−χ)w + χ·mean | **Parameter**: R = N/2, χ = τ̄f². Per agent-transaction this is a tax of τ̄f², which matches τΔt = τ̄γΔt = τ̄f². Discretization error is O(χ) ≤ 0.0025 |
| **Batching gotcha** | n/a | `step(n)` applies all of a batch's redistributions **at the end** of the batch | **Calling convention**: always call `step(R)` (or `run(..., record_interval=R)`) so each tax follows its own R transactions |
| Limits | Fokker–Planck: β → 0, continuum, random-agent approximation | finite f, finite N | Check f at two values; N = 5000 (O(1/N) Gini bias ≈ 2×10⁻⁴) |
| P28-B wealth scale | W/N = 100, infinite precision | float64 | w0 = 100. The poorest agents reach ln w ≈ −160 by 8×10⁷ transactions, which is within float64 range; the rounding loss on rich agents is ~10⁻¹⁶ relative |

### 3. Protocol

**P28-A**
- N = 5000, f = 0.05, λ_save = 0, R = 2500, χ = τ̄·f² (5×10⁻⁴, 1.25×10⁻³, 2.5×10⁻³), `init_mode='equal'`.
- 10 seeds (2801–2810) per τ̄.
- Burn-in: t̄_burn = max(20/τ̄, 20), which is 100, 40 and 20.
- Measurement: t̄ = 50, sampling the Gini after every tax application.
- Estimator: time-average of G per seed, then mean ± SE across seeds.
- Stationarity, fixed now: the two halves of the measurement window must agree within 0.005. If they do not, double the burn-in once and rerun.
- Finite-f check, required: τ̄ = 0.5 at f = 0.025 (χ = 3.125×10⁻⁴), 10 seeds. If |G(0.05) − G(0.025)| > 0.01, the gated value at every τ̄ becomes the f² → 0 linear extrapolation from f ∈ {0.05, 0.025}, which requires running f = 0.025 at all three τ̄.
- Reported, not gated: τ̄ = 0.5 started from `init_mode='exponential'` (G₀ ≈ 0.5, above G∞), to bracket the steady state from both sides.

**P28-B**
- N = 5×10⁴, f = 0.25, χ = 0, w0 = 100, 10 seeds (2821–2830), run to t = 8×10⁷ transactions.
- At 16 equally spaced t, compute A(w, t) = fraction with wealth > w. Fit A against ln w by least squares on A ∈ [0.2, 0.8] (the fit window is fixed now; the source does not specify one). Set w_min where the fit reaches A = 1, then solve Eq. (37) w_min = (W/2a)·csch(N/2a)·exp(−N/2a) for a.
- Least-squares fit of 1/a against t gives the slope and intercept per seed.

### 4. Tolerance and pass rule (fixed)

- **P28-A PASS:** |Ḡ_sim(τ̄) − G_fit(τ̄)| ≤ 0.02 at each of τ̄ = 0.2, 0.5 and 1.0.
  - Error sources: fit error (~0.005), seed SE (~0.002), and finite-f and finite-N effects (≤ 0.01; the Kramers–Moyal truncation is O(f²) because odd moments vanish).
  - A factor-2 error in the τ̄ mapping shifts G by 0.09–0.10, so the test can tell the mapping is right.
- **P28-B PASS:** the mean slope across seeds is within ±15% of 2.04046×10⁻⁵, i.e. [1.73, 2.35]×10⁻⁵. The intercept is reported only.
  - The wide band is because the fit-window choice is ours, not the source's.
  - A slope near 1/N is the paper's stated reading.
- The P28 reproduction is reported as PASS only if P28-A passes. P28-B is listed separately.

### 5. Compute and risks

- The pure-Python loop takes about 1.5–2 µs per transaction.
- P28-A at f = 0.05: n_tx = t̄·N/(2f²) gives 3.0, 1.8 and 1.4×10⁸ per seed, scaled for N = 5000. About 1.5 CPU-h for 30 runs.
- The f = 0.025 check adds about 1.7 CPU-h, or about 5 more if extrapolation is triggered.
- P28-B: 8×10⁷ transactions per seed, about 2–3 min each. About 0.5 CPU-h.
- Gini sampling cost is negligible (O(N log N) per sample).
- Risks not caused by the model:
  - Eq. (61) is an empirical fit to Fokker–Planck solutions, not a closed form, and has no stated error.
  - The τ̄ mapping depends on reading Δt as one transaction per agent. Eqs. (8) and (13) support this; the test above can detect an error.
  - The P28-B estimator is underspecified in the source.
  - The canonical P28 panel uses f = 0.5 with no tax, so P28-A validates the code path and physics, not the panel regime. P28-B, at f = 0.25 with no tax, is closer.

---

## Sources read for this spec

- Gammaitoni et al. 1998 RMP (Augsburg OPUS PDF): Eqs. 2.1–2.15, 4.4–4.8, 4.21, 4.48–4.54; Figs. 6–8.
- Reichenbach, Mobilia & Frey 2007, arXiv:0709.0217v2 (main text, Methods, SI, Fig. 2b); Reichenbach, Mobilia & Frey 2008 JTB, arXiv:0801.1798 (simulation algorithm and D, p. 4; Eq. 49, p. 18).
- Theraulaz, Bonabeau & Deneubourg 1998 Proc. R. Soc. B (IAS archive PDF of PMC1688885): Eqs. 1–7, Figs. 1–2.
- Berdahl et al. 2013 Science (KOPS PDF): full text and Fig. 1B.
- Boghosian 2014 PRE, arXiv:1212.6300v2 (Secs. II–III, Figs. 3–5, Eqs. 22, 32–39); Boghosian et al. 2017 Physica A, arXiv:1511.00770v2 (Eqs. 8, 13, 31–33, 51, 61; Fig. 4).
- Not accessible: Bonabeau, Theraulaz & Deneubourg 1996 (publisher returned 403; closed access); Berdahl 2013 supplement (403).
