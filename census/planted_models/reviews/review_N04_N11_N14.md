# Rule-by-rule review of plug-ins N04, N11, N14 (agent report, 2026-11-03)

Agent "Review plug-ins N04 N11 N14" finished

# Review of N04, N11 and N14 against their papers

I edited no files. I didn't import or call any census pipeline code, and I computed no scores or fingerprints. Test scripts are in `/private/tmp/claude-501/-Users-matthewmaxwell-Downloads/ed5d6e8e-67aa-43e5-b2b6-30638ce75e8d/scratchpad/rev/`. The papers, the authors' N11 code and their IBM cache are in the sibling `papers/` and `hoi/` folders.

| Module | Verdict | Blocking issue |
|---|---|---|
| N04 | **MUST FIX** | Transverse-force sigma contradicts the paper's Table 1; passive melting point about 3x the paper's, cause not resolvable |
| N11 | **MUST FIX** (a plan deviation, not code) | The plan's headline value 3.14e-5 is the paper's gap pattern, not stripes |
| N14 | **OK with notes** | None |

All three pass `base.check_out` for the real program, the knock-out and the sub-programs. `keep` switches exactly the plan's interaction rules (N04: 1, 2; N11: 3, 4; N14: 1, 2). Every `presence()` reads only output arrays. Every `RULE_TEXT` is the frozen card text edited only for size, length and setting, with no coined terms, title words or author names.

## N04: what causes the reversed headline

**What is faithful:**
- The forces match Eqs. 16–18 exactly: a finite-difference check of force and torque against the HOOMD-doc energy gives relative error 3e-10.
- Things the caller listed that are not the cause:
  - **Gay-Berne form:** HOOMD's form has no mu/nu exponents and no epsilon anisotropy, and neither does the module.
  - **LJ prefactors:** both are standard, so there is no factor-4 or factor-24 slip.
  - **Half-axes:** 0.55 and 0.65 are semi-axes, giving contact distances 1.1 to 1.3.
  - **Density and friction:** number density 0.87, Gamma = 10 and Gamma_r = 0.1 all match.
- At T = 0 the module's elastic constants match Table 1 in e_p units: mu = 42 (paper 40), kappa = 18 (paper 16), B = 875 (paper 1400).
- The KTHNY estimate k_B T_c = Y a^2/(16 pi) is 4.3 from the module's moduli and 4.1 from Table 1. Both agree with Fig. 1(a): at e_a/e_p = 0.1 the paper shows solid at T = 4, hexatic at 4.25 and liquid from T = 5.

**Single changes tested, at 400 particles for 10 time units:**

| Change | Passive melting | Active result | Conflict |
|---|---|---|---|
| Isotropic in-plane core (sigma = 1.1 for every pair; HOOMD's Gay-Berne with a body-z director in 2D) | Melts between 4.25 and 6, matching the paper | No instability at (1, 1), (2, 0.05) or (5, 1) | Zero torque contradicts kappa = 16 |
| Frozen, aligned axes | Still solid at T = 6 | — | — |
| Integrator heating (Gaussian or HOOMD-style uniform noise) | Solid at a stable dt = 5e-4; the paper's dt = 1e-3 blows up at T &gt;= 3 | — | — |
| Uniform energy rescale, e_p = 1/3 | Solid at 4, liquid at 6 | Matches (0.2, 3), (2, 0.05), (9, 1), (33, 4), (27, 6): 8 of 9 points in all | Headline (33, 6) still disordered; no ordered state above passive melting at activity 33–45, T = 5–6; no basis in the paper |
| Transverse sigma = 1.16 | Unchanged, still solid at T = 6 (psi6 0.83) | Restored order at (33, 6): psi6 0.87, no defects | — |

The e_p rescale works because the module's whole phase diagram is the paper's with both the T and e_a axes multiplied by about 3. Table 1 cannot rule this out, since moduli quoted in units of e_p are scale-free. Neither can anything else in the paper.

Bottom line: no single change reconciles both the passive melting point and the active headline. The passive melting temperature depends only on the Gay-Berne form, the density and T, and all three are stated. So the factor of 3 cannot come from anything the paper leaves unstated (transverse sigma, cut-offs, dt, friction ratio, initial condition). One likely contributor is thermal randomisation of the axes at T of about 4 or more (nematic order about 0.01): it raises the effective contact distance at a = 1.15, and the r^-12 core makes the melting point scale as about (1.2/1.1)^12, roughly 3. This is consistent with the isotropic-core test but not proven. The paper's own diagram also puts the headline in its last column, with activity of about 27 already liquid at T = 6, so the effect sits in a thin strip. At the plan's headline the module shows activity-driven defect proliferation over a passive crystal. That is the paper's other new behaviour, not restored order.

**N04 findings:**
1. **Must fix:** `HEADLINE["sigma_a"]=1.0` (line 105; used in `forces`, lines 189–192). With sigma = 1 the transverse force crosses zero at r = 1.122, right at the neighbour spacing a = 1.15. The resulting odd moduli are A = 24–35 e_a and K^o = 13–19 e_a, against Table 1's A = 400 e_a and K^o = 140 e_a. Fix: sigma_a of about 1.16 matches both (A about 330, K^o about 140), and any factor-2 convention still places it between 1.1 and 1.2. That needs dt of about 1e-5 (stable for 4 time units at the headline), so STEP must drop to about 0.02 time units to stay under the compute cap. Record this in the docstring.
2. **Must record:** the passive melting point is about 14 against the paper's about 5, and this is unresolvable from the paper (details above). The plan's fidelity check will fail with or without fix 1, because the passive system at T = 6 stays a crystal. Report N04 as "implementation not confirmed" and log the mismatch in `PLANTED_PIPELINE_DEVIATIONS.md` before the freeze.
3. **Note:** STEP is 0.2 time units, but the cap allows about 0.3 at the measured 1.08 ms per step; the plan says "the largest that fits". Also, one random generator (line 245) is shared by all seeds, so a seed's trajectory depends on S. Finally, the 400-particle square box imposes a 6% lattice strain.
4. **Note:** dt = 1e-3 blows up in the module at T &gt;= 3, so the paper must also have used a smaller step in that range; its value is not stated.

## N11

1. **Must fix (plan deviation).** Fig. 1 shows 3.27e-6 as spots, 2.35e-5 as stripes, 3.14e-5 as gaps and 3.93e-5 as uniform. The authors' README table agrees (stripes 2.356e-5, gaps 3.1414e-5). I recomputed the authors' cache `all_coherence_IBM.npy`: anisotropy A = 0.025, 0.916, 0.169 and 0.100 respectively. The planted card's labels are shifted by one panel, and the plan inherited that error. The plan's check "stripes at the headline" therefore fails for a correct implementation (`fidelity`, line 338). Before the freeze, log a P1-style entry choosing one of:
   - move the headline to 2.356e-5, consistent with the plan's label and check, although stripes are still forming at t = 0.6 and the census run lasts 1 time unit; or
   - keep 3.14e-5 and restate the check as gaps against spots, using the authors' A values.
   Either way, update `HEADLINE["gamma"]` (line 92), `RULE_TEXT` (line 105) and `_AUTHORS_A` (line 108) to match.
2. **OK:** the rules match the authors' `01_run_ibm_sweep.py`:
   - birth at 1920 with Gaussian sd 0.05 per coordinate;
   - death at 192 + C_p/(pi 1.25^2) + gamma (C_ho/(pi 0.0625^2))^2, counting neighbours within the radius (&lt;=) and excluding self;
   - offspring inherit the parent's rate until the refresh every 1e-4.
   The module's exact branching between refreshes replaces the paper's tau-leap, which carries about a 1% bias, and is documented. The grid pair count is accurate to 0.09% rms against exact counts, and the three-body count is exact.
3. **Should fix:** the domain (side 10, about 160,000 individuals) departs from the plan's literal rule. It is the only reading that keeps the behaviour observable: about 400 individuals at the steady density gives a side of about 0.5, smaller than the competition radius 1.25, and "close packing" has no meaning for point individuals. Log it in the deviations file, not only the docstring.
4. **Should fix (owner decision):** the 50x blow-up rule. The real program settles at 40.6x the start, but the pair-only sub-program settles at 47–51x and crossed the threshold at step 13 in my test, so its verification outcome is a coin flip. Decide and log before any result is read.
5. **OK:** the compute-cap reduction to 1e-3 time units per step is justified. The plan's 1/192 would take about 3–4 CPU-hours; I measured 0.077 s per refresh, giving 39 minutes at 1e-3.
6. **Note:** the area4 verification (about 650,000 individuals, 10 seeds, 3x length) will likely exceed 24 CPU-hours and be recorded as not run.

## N14

1. **OK:**
   - The overlap test agrees with brute force: 0 disagreements in 20,000 random pairs.
   - Bonds follow the Kern–Frenkel rule on the vertex osculating circles: distance &lt;= 2 r_v + 0.17, 45 degree cones on both patches, -1 per bonded patch pair. The model paper's formula r_a = a^2/b is swapped; the module correctly uses b^2/a = 0.177.
   - No three-bond patch was found in a targeted search of 400,000 configurations.
   - A perfect kagome gives Gamma_6 = 0.974, n_L = 1, 56 six-loops and 112 triangle motifs.
   - Insertion/deletion acceptances are correct for mu* = beta mu with equal proposal probabilities. The ideal-gas knock-out gives &lt;N&gt; = 21.7 against the expected 20.6.
   - The static-linking cluster moves satisfy detailed balance for square-well bonds: the failed-link factors cancel the Boltzmann factor, and the rotation test uses distance from the centroid.
   - T = 0.145 and mu = -4 sit at Fig. 5(n)'s maximum (Gamma_6 about 0.79, n_L about 0.97, rho about 0.42).
2. **Note:** `_insert`/`_delete` (lines 280, 286) assume a thermal wavelength of 1 and no 2 pi orientation factor, which the paper doesn't state. A ln(2 pi) shift of about 1.84 could move mu = -4 to about -5.8, near the gas boundary.
3. **Should fix (log):** STEP uses the paper's cycle of 1296 moves, but the plan defines a step as "one sweep (one attempted update per entity)". These differ while N grows from 5 to 545. Separately, no 50x blow-up rule is applied, although the real program grows about 100x by design. Both are reasonable choices but should go in the deviations file.
4. **Note:** a 1,000-cycle census run is about 0.05% of the paper's 2e6-cycle start-up (the development run reached Gamma_6 = 0.03), so presence will likely be negative. The fidelity run is also unlikely to reach the paper's state within 2 CPU-hours.
