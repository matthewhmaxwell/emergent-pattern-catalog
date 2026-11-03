# Fidelity checks of the planted-model plug-ins (PLANTED_PIPELINE_PLAN.md; read before the code freeze)

Each check reproduces one observation from the model's paper with the paper's own measure (files in fidelity/).
Budgets were 2 hours of wall-clock time per check on a shared 8-core machine that was oversubscribed (load 11 to 16)
during the runs, so several checks ran smaller or shorter than the budget would otherwise allow.

| Card | Result | What was measured |
|---|---|---|
| N01 | confirmed | exact enumeration (two code paths) 227/256, 3124/3125, 42275/46656 = paper eq. 4; Monte Carlo z -0.86, 0.09, 0.30 |
| N02 | not confirmed | amplitude larger at c = 0.17 than at 0.01 (passes), but a strong pattern at c = 0.49 where the paper reports none; the paper's own Eq. 22 predicts instability there (paper inconsistency, logged before the run) |
| N03 | not confirmed | 1,000 disks, 2,460 time units: one cluster (S = 1, sigma_S = 0), contact number 5.14, no cavities; the printed equations bind small clusters at every screening value (reviewer), so the paper's states cannot be reached |
| N04 | not confirmed | headline ordered (psi6 0.87, no defects) but the passive system at T = 6 is also ordered (0.82); the implementation's passive crystal melts near T = 14 against the paper's ~5 |
| N07 | confirmed | 800 x 800 lattice: chiral correlations decay as r^-eta with eta in [3.5, 4.5]; the achiral version fails the same test |
| N08 | not confirmed | 2,048 dumbbells for 5,847 time units at Pe = 5: one density peak (phi 0.648, sd 0.008), no coexistence; the paper notes small systems do not segregate |
| N11 | confirmed (P2 rule) | anisotropy 0.107 at g = 3.14e-5 (authors 0.169) and 0.016 at 3.27e-6 (authors 0.025), spot psi6 0.89; runs reached t = 0.43 and 0.21 (paper 3); the plan's original "stripes" check fails, as expected for the authors' gap pattern |
| N12 | confirmed | M^2(L) at beta = 1 within 0.01 of the paper's fit for L >= 128; census heat-bath agrees with the cluster sampler |
| N13 | confirmed | exact and simulated stationary distributions: three peaks at -1, 0.3, +1; KS within 0.02 |
| N14 | not confirmed | 2,034 cycles (paper: >= 2 million start-up): Gamma6 0.10, fraction in loops 0.60, a branched network, not the paper's ring lattice |
| N15 | confirmed | 10,000 particles, 210 time units: correlation dimension 1.42, structure length not growing (ratio 0.98) |
| N17 | not testable | the paper's lanes need boxes of side ~1,200 (plan) |
| N19 | confirmed | point K (deviation P1), 50 x 50, 3 seeds: one dense aggregate in travelling bands in the majority of late frames |

Confirmed 7 of 12 plug-ins; not confirmed 5 (N02, N03, N04, N08, N14); N17 not testable. Under the plan, all 13
in-scope models run through the pipeline; results for the five unconfirmed ones carry "implementation not confirmed".
