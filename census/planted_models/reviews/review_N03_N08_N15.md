# Rule-by-rule review of plug-ins N03, N08, N15 (agent report, 2026-11-03)

Agent "Review plug-ins N03 N08 N15" finished

## Verdicts

| Module | Verdict | What decides it |
|---|---|---|
| N03 | **OK with notes** | The code matches the printed model term by term. The phase-diagram mismatch is inside the paper, not the code. |
| N08 | **MUST FIX** (one item) | The census step uses the Lennard-Jones time unit; the paper declares 1/Dr as its time unit. Everything else is correct. |
| N15 | **OK** (notes only) | Steering, cone, inertia and integrator all match the paper. |

I edited nothing and called no forbidden census modules. I read all three papers from arXiv, viewed N03 Figs. 2–3, N08 Fig. 2 and N15 Fig. 2, and ran only small tests.

## N03 (`census/planted_models/n03.py`)

**On the author's concern: there is no code error, and no stated parameter closes the gap.**
- Each term matches the paper:
  - J_i (Eq. 4) and the image term (separation √(r² + (2h)²), in-plane part only).
  - V_self = +v0·p and V_int = +2·v0·J.
  - Ω = −(9/4)·v0·μr·(p × J), with p × J = cosθ·J_y − sinθ·J_x.
  - The WCA force of Eq. (22), Pe → Dr = 0.02, and Λc = 1 for the no-flux wall.
- Tests:
  - The current from one neighbour points toward it, with the expected magnitude (0.0828).
  - The torque turns the heading away from J and matches the formula to four digits (7.910e-4 against 7.907e-4).
  - The output format check passes, and the knock-out keeps only rules 1 and 4.
- Prefactors: I re-derived Eqs. (2)–(3) from the appendix's dimensional forms (Eqs. 13, 14, 20, 21, with ⟨A⟩ = A/2 and v0 = |⟨μ⟩|A/(4D)). They are self-consistent and physically correct. There is no hidden 1/(4π) or 1/b.
- The printed model cannot give the paper's gas. With printed Eq. 4 and Λc = 1, an edge disk whose heading points outward escapes only if its inward current is below 0.5 (from v0 &gt; 2·v0·J). At κ = 0.28 the static edge current is:
  - 0.80 for a 7-disk cluster up to 1.30 for 1459 disks, with disks 2.0 apart;
  - 0.95 to 1.58 with disks 1.83 apart.

  So any cluster of 4–7 or more disks is bound for every κ up to 0.5. The same holds for the paper's Λc = 0 gas near κ ≈ 0.1 (Fig. 3f): the edge current there is 0.66–1.19 at κ = 0.16.
- The unstated choices (viscosity, contact range) change the disk spacing between about 1.83 and 2.0, which moves the edge current by about 20%. Closing the gap would need a factor of about 2.5–3.
- I also tried reading r_ij in Eq. 4 as a unit vector. That overshoots: the edge current is below 0.5 at every κ, so there is no boundary near 0.25 either. This rules out the mirror-image, sign, torque-coupling and units explanations.
- So "implementation not confirmed" is the correct plan outcome if the fidelity check fails.

**Findings**
1. **Should fix (needs the owner's ruling) – time unit.** `STEP` (line 92) and `steps_per_census_step = 1000` use the simulation unit. The paper never declares a time unit and reports all times as t/τ (τ = 1/Dr = 50 units). The census run is therefore 20 τ, which ends before the paper's boiling state forms (Fig. 3a: S reaches 1 at about 25 τ; Fig. 3b: Q levels off around 30 τ).
   - If the ruling is τ per step, the 60-minute compute cap applies.
   - Measured cost is 0.05 s per time unit for the real run with clusters and 0.009 s for the knock-out, so the cap gives about 20 time units per census step, a 400 τ run.
   - In this code, clusters form by about 4 τ (S = 0.94 at t = 200), so the choice mainly matters for consistency across modules.
2. **Note – cavity test (`cavities`, lines 236–248).** A region counts as enclosed when its mouth is narrower than 4b centre to centre, which is looser than the paper's bound distance of 2.1b. The minimum area of π b², counted in space more than 2b from every disk, already implies a hole of roughly 4 or more missing disks, so it is acceptable. A stricter version would seal at 1.05b.
3. **Note – citations in the docstring.** The WCA potential is the paper's Eq. (22), not "B1", and the table is "Table 1". The appendix forms are Eqs. 13, 14, 20 and 21.
4. **Note – rule text (lines 94–100).** Besides the allowed numbers, it adds "(square of side 158.5)" and "(a free disk moves 2 radii per time unit)". These leak nothing, but they go beyond the plan's numbers-only edits. There are no coined terms, title words beyond "screening" (already in the frozen card), surnames or arXiv id.
5. **Note.** `simulate` passes only κ and μr; every other value comes from `_Sys` defaults that duplicate `HEADLINE`. They agree today but could drift if one is edited.
6. **Confirmed correct:**
   - The headline values (κ = 0.20, μr = 0.30): in Fig. 2, read off the figure, boiling appears at κ = 0.20 for μr from 0.15 to 0.40.
   - The 25b cutoff and the per-disk displacement cap of 0.01 (it rarely binds at dt = 0.001).
   - Positions are recorded in [0, L), and the heading is the propulsion direction.
   - `presence()` uses only the output arrays.
   - Fidelity: N = 1000, tolerance fixed in advance, and S, σS, Q measured exactly as Eq. 7 defines them.

## N08 (`census/planted_models/n08.py`)

**Findings**
1. **Must fix (or log as a deviation with the owner's ruling before the freeze) – census step.** Lines 27–30 and line 85 (`STEP`, `steps_per_census_step = 500`) use the Lennard-Jones time unit. Appendix A states "(m, σ, D_r⁻¹) units of mass, lengths and time" and "Using D_r⁻¹ as timescale". By the plan, one census step is therefore 1/Dr = 0.5 Lennard-Jones units, which is 250 sub-steps at dt = 0.002 and a 500-unit run. If changed, "Run 1,000 time units" in the rule text needs to name the unit.
2. **Resolved – Prandtl number.** The "Pr = 200" in Fig. 2's caption is 1/Pr: 1/0.005 = 200, which is the coefficient of the second time-derivative of r in the dimensionless Eq. A2 ("large inertia"). The values γ = 0.01, m = 1, Dr = 2, Pe = 5 → Fa = 0.1 and Re = 1000 are all consistent.
3. **Note – FENE spring.** The paper prints K·R0²·ln[1 − (r/R0)²], which has the wrong sign. The module uses the standard −½·K·R0²·ln form, with WCA also acting inside the bond. That is right for K = 30, R0 = 1.5 in LAMMPS, which the paper uses; the bond rests at 0.961.
4. **Confirmed correct:**
   - The Langevin noise is √(2γT/dt), and the velocity-Verlet integrator applies friction at the half-step velocity.
   - The torque turns every dumbbell clockwise. In a test with near-zero noise, one dumbbell spun up as ω∞(1 − e^(−γt/m)): −1.01, −1.97 and −2.88 at 5, 10 and 15 units.
   - The bond stays at 0.961–0.97 and the run is stable.
   - The knock-out removes only the bead repulsion between dumbbells and keeps the spring's internal WCA.
   - The heading is the dumbbell axis, and the output format check passes.
   - Fig. 2 places φ = 0.65, Pe = 5 inside the spinodal (a "stripe" symbol). The binodal [0.55, 0.73] is stated in the text (Appendix A.2), not only read off the figure.
   - Fidelity uses the paper's local area-fraction measure from bead positions on 20σ cells.
5. **Note – fidelity budget (`fidelity`, lines 316–327).** The size chooser does not count the start-up time (packing and relaxation, about 20 s or more per run, growing with size). The Pe = 0 run's reserve can therefore come out a little short of 300 units. Measure the start-up time and add it to the reserve.
6. **Note.** Headings are sampled every 2 census steps; fast-spinning dumbbells in dilute regions will alias. This is inherent to the census sampling.

## N15 (`census/planted_models/n15.py`)

**Findings**
1. **Confirmed correct – the three things you asked about:**
   - Steering: Eq. 2 with the 1/n_i normalisation, sin(α_ij − φ_i), and no steering when no neighbour is visible.
   - Vision cone: θ is the half-angle with a strict inequality (Eqs. 3–4), seeing is one-sided, and minimum image works. Tests matched 3·sin(9.46°) = 0.4932 exactly, including a pair across the periodic boundary.
   - Inertia and friction: m·r̈ = −γ·ṙ + γ·v0·e(φ) with m = 1, γ = 50, v0 = 0.1, Dr = 0.08, Ω = 3, dt = 0.005.
2. **Confirmed correct – plan items:**
   - Time is 200 steps per census step, as the plan's exception requires.
   - Box 20 at density 1; headline θ = π/8, Rv = 5 from the Fig. 1(c) caption. In Fig. 2 that point is foam with D ≈ 1.4, read off the heatmap.
   - The knock-out is Ω = 0, and the output format check passes.
   - The rule text changes only numbers, with no coined terms.
3. **Note – Eq. 6.** It is printed as log C / log r, but the text defines D as the log–log slope. The module correctly uses the slope; the docstring should say so.
4. **Note.** The fit range [1, min(2Rv, L/4)] is the author's choice (the paper gives none) and was fixed in advance. The non-coarsening test (length ratio ≤ 1.25) is also the author's, since the paper has no quantitative measure. Both are fine.
5. **Note.** Cost is 0.066 s per census step from a uniform start.

My test scripts are in `/private/tmp/claude-501/-Users-matthewmaxwell-Downloads/ed5d6e8e-67aa-43e5-b2b6-30638ce75e8d/scratchpad/`. `jedge.py` and `jedge2.py` hold the edge-current estimate; `t03.py`, `t03b.py`, `t08.py`, `t08b.py` and `t15.py` hold the module tests.
