# Simulation verification: every phenomenon and tested feature

Medians over 10 fresh seeds at census settings. e = median(real) - median(knock-out); s = robust spread; e/s in brackets. 4x and 3x are the effects at 4x area and 3x run length. 'ko absent' = the knock-out never produces the view (deviation V2). card = value on the literature card; rep = replicates the card (reported only). R = robust (interaction-driven and persistent).

## census:g1:lattice-phase-1  (SINGLE MECHANISM)

`C kC=1 | C.KURA(K=0.1)` · view C.phase · flagged False · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| p_local_r | 0.6704 | 0.4502 | 0.2202 [9.8] | 0.2176 | 0.2225 | 0.591 | yes | R |
| p_freq_spread | 0.3203 | 0.5009 | -0.1805 [-11.0] | -0.1784 | -0.7314 | 0.403 | yes | R |
| p_vortex_density | 0.1016 | 0.3359 | -0.2344 [-26.6] | -0.2264 | -0.229 | 0.178 | no | R |
| p_r | 0.0521 | 0.0143 | 0.0377 [8.6] | 0.018 | 0.0409 | 0.035 | yes |  |
| p_r_std | 0.0216 | 0.0068 | 0.0148 [7.4] | 0.0057 | 0.0182 | 0.014 | yes |  |

## census:g1:network-phase-1  (SINGLE MECHANISM)

`N kN=1 | N.KURA(K=0.1)` · view N.phase · flagged False · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| p_freq_spread | 0.3404 | 0.504 | -0.1635 [-6.3] | -0.1547 | -0.5985 | 0.351 | yes | R |
| p_r | 0.1742 | 0.0577 | 0.1165 [8.8] | 0.0725 | 0.1347 | 0.233 | yes | R |
| p_r_std | 0.0674 | 0.0304 | 0.037 [3.7] | 0.0235 | 0.0616 | 0.08 | yes | R |

## census:g1:network-1  (SINGLE MECHANISM)

`N kN=2 | N.GAME(g=(1.9, 0.0))` · view N · flagged False · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| n_components | 1.498 | 1.498 | 0 [0.0] | 0 | 0 | 1.68 | yes |  |
| n_giant | 0.9875 | 0.9875 | 0 [0.0] | 0 | 0 | 0.981 | yes |  |
| n_type_entropy | 0.6879 | 0.6913 | -0.0034 [-0.1] | -0.0063 | -0.0034 | 0.944 | no |  |
| n_oscillation | 0.134 | 0 | 0.134 [49.6] | 0.2474 | 0.0531 | 0.559 | no |  |
| n_degree_cv | 0.4765 | 0.4765 | 0 [0.0] | 0 | 0 | 0.489 | yes |  |
| n_type_modularity | 0.1412 | 0.0249 | 0.1164 [4.0] | 0.1626 | 0.1189 | 0.185 | yes | R |

## census:g1:grid-1  (SINGLE MECHANISM)

`C kC=2 | C.COPY(p=1.0)` · view C · flagged False · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_n_domains | 5.879 | 6.396 | -0.5171 [-2.1] | -0.5223 | -0.8792 | 5.514 | yes |  |
| g_n_types | 2 | 2 | 0 [0.0] | 0 | 0 | 2.32 | yes |  |
| g_corrlen | 1 | 1 | 0 [0.0] | 1 | 0 | 1.573 | no |  |
| g_reversals | 1 | 0 | 1 [50.0] | 1 | 1 | 0.863 | yes | R |
| g_type_entropy | 0.9493 | 1 | -0.0507 [-1.3] | -0.0067 | -0.3346 | 0.821 | yes |  |
| g_largest_domain | 0.5736 | 0.0441 | 0.5295 [12.7] | 0.2789 | 0.7131 | 0.547 | yes | R |

## census:g2:lattice-phase-2  (SINGLE MECHANISM)

`C kC=1 | C[any] ALWAYS() -> PHASE_COUPLE(K=0.1) p=0.1` · view C.phase · flagged False · verdict COVERED-BY-TWIN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| p_local_r | 0.6623 | 0.451 | 0.2113 [9.5] | 0.2117 | 0.221 | 0.578 | yes | R |
| p_vortex_density | 0.1099 | 0.3284 | -0.2185 [-24.8] | -0.2222 | -0.2334 | 0.194 | no | R |
| p_r | 0.0379 | 0.0165 | 0.0214 [2.5] | 0.0197 | 0.0322 | 0.029 | yes |  |
| p_freq_spread | 0.0373 | 0.4974 | -0.46 [-43.0] | -0.4625 | -1.312 | 0.027 | no | R |
| p_r_std | 0.0051 | 0.0063 | -0.0012 [-0.6] | -0.001 | 0.0063 | 0.003 | yes |  |

## census:g2:network-phase-1  (SINGLE MECHANISM)

`N kN=1 | N[any] ALWAYS() -> PHASE_COUPLE(K=0.1) p=1.0` · view N.phase · flagged False · verdict COVERED-BY-TWIN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| p_r | 0.1941 | 0.0657 | 0.1284 [6.6] | 0.0672 | 0.1679 | 0.207 | yes | R |
| p_freq_spread | 0.3352 | 0.4925 | -0.1574 [-9.5] | -0.1522 | -0.6072 | 0.142 | no | R |
| p_r_std | 0.0784 | 0.0332 | 0.0453 [4.2] | 0.0272 | 0.0533 | 0.04 | yes | R |

## census:g1:grid-2  (SINGLE MECHANISM)

`C kC=2 | C.GAME(g=(1.5, -0.5))` · view C · flagged False · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_corrlen | 3 | 1 | 2 [25.0] | 2 | 2 | 2.761 | yes | R |
| g_n_domains | 1.386 | 6.375 | -4.989 [-32.1] | -5.275 | -4.989 | 2.406 | yes | R |
| g_n_types | 1 | 2 | -1 [-16.7] | -1 | -1 | 1.567 | no | R |
| g_largest_domain | 0.9929 | 0.0459 | 0.947 [45.5] | 0.9754 | 0.947 | 0.958 | yes | R |
| g_moran | 0.6282 | 0.0002 | 0.628 [49.8] | 0.6412 | 0.628 | 0.549 | yes | R |
| g_reversals | 0 | 0 | 0 [0.0] | 0 | 0 | 0.254 | no |  |

## census:g1:network-6  (SINGLE MECHANISM)

`N kN=3 | N.CYCLE(th=3,p=0.1)` · view N · flagged False · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| n_components | 1.701 | 1.701 | 0 [0.0] | 0 | 0 | 1.682 | yes |  |
| n_type_entropy | 1.085 | 1.094 | -0.0098 [-0.2] | -0.0038 | -0.0098 | 1.168 | yes |  |
| n_giant | 0.9825 | 0.9825 | 0 [0.0] | 0 | 0 | 0.98 | yes |  |
| n_cyclic_changes | 0.448 | -0.5 | 0.948 [49.9] | 0.9454 | 0.8901 | 0.535 | yes | R |
| n_degree_cv | 0.4897 | 0.4897 | 0 [0.0] | 0 | 0 | 0.492 | yes |  |
| n_oscillation | 0.1083 | 0 | 0.1083 [14.2] | 0.1874 | 0.0509 | 0.402 | no |  |

## census:g1:network-15  (SINGLE MECHANISM)

`N kN=3 | N.COPY(p=0.05)` · view N · flagged False · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| n_components | 1.099 | 1.099 | 0 [0.0] | 0 | 0 | 1.905 | no |  |
| n_giant | 0.995 | 0.995 | 0 [0.0] | 0 | 0 | 0.974 | yes |  |
| n_type_entropy | 0.8043 | 1.095 | -0.2905 [-7.6] | -0.0883 | -0.7291 | 0.712 | yes |  |
| n_oscillation | 0.6884 | 0 | 0.6884 [49.9] | 0.6072 | 0.6983 | 0.64 | yes | R |
| n_degree_cv | 0.4989 | 0.4989 | 0 [0.0] | 0 | 0 | 0.507 | yes |  |
| n_consensus_gain | 0.275 | 0 | 0.275 [24.8] | 0.1731 | 0.54 | 0.337 | yes | R |

## census:g1:grid-4  (SINGLE MECHANISM)

`C kC=4 | C.CYCLE(th=1,p=1.0)` · view C · flagged False · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_n_domains | 7.55 | 7.666 | -0.1156 [-0.4] | -0.11 | -0.1156 | 7.297 | yes |  |
| g_n_types | 4 | 4 | 0 [0.0] | 0 | 0 | 3.694 | yes |  |
| g_corrlen | 1 | 1 | 0 [0.0] | 0 | 0 | 1.083 | yes |  |
| g_type_entropy | 1 | 0.9997 | 0.0003 [0.0] | 0 | 0.0003 | 0.99 | yes |  |
| g_activity | 1 | 0 | 1 [50.0] | 1 | 1 | 0.848 | yes | R |
| g_interface | 0.7134 | 0.7484 | -0.035 [-1.2] | -0.0327 | -0.035 | 0.643 | yes |  |

## census:g1:network-12  (SINGLE MECHANISM)

`N kN=2 | N.REWIRE(p=0.01)` · view N · flagged False · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| n_components | 1.386 | 1.946 | -0.5596 [-1.4] | -0.4527 | -0.3365 | 1.344 | yes |  |
| n_giant | 0.99 | 0.975 | 0.015 [0.4] | 0.0063 | -0.4525 | 0.989 | yes |  |
| n_type_entropy | 0.692 | 0.692 | 0 [0.0] | 0 | 0 | 0.799 | yes |  |
| n_type_modularity | 0.4904 | -0.0049 | 0.4953 [50.0] | 0.4934 | 0.5003 | 0.485 | yes | R |
| n_edge_turnover | 0.4941 | 0 | 0.4941 [49.9] | 0.4931 | 0.5006 | 0.471 | yes | R |
| n_degree_cv | 0.4502 | 0.494 | -0.0438 [-2.2] | -0.0374 | -0.0489 | 0.454 | yes |  |

## census:g1:agents-1  (SINGLE MECHANISM)

`A kA=1 v=1.0 eta=0.3 | A.ATTRACT(b=any,r=1.0)` · view A · flagged True · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| a_density_cv | 2.209 | 0.8089 | 1.401 [20.2] | 0.7113 | 1.449 | 3.762 | no | R |
| a_step | 0.2751 | 6.601 | -6.326 [-46.0] | -19.48 | -7.033 | 0.395 | no | R |
| a_nn_ratio | 0.1752 | 0.99 | -0.8148 [-35.0] | -0.8191 | -0.8327 | 0.142 | yes | R |
| a_local_align | -0.2056 | -0.0393 | -0.1663 [-23.8] | -0.1679 | -0.1586 | -0.102 | no | R |
| a_polar | 0.0055 | 0.0476 | -0.0421 [-38.3] | -0.018 | -0.03 | 0.013 | no |  |
| a_ang_mom | 0.0052 | 0.0277 | -0.0224 [-2.1] | -0.0122 | -0.0177 | 0.011 | yes |  |

## census:g1:agents-2  (SINGLE MECHANISM)

`A kA=1 v=1.0 eta=0.3 | A.REPEL(b=any,r=3.0)` · view A · flagged True · verdict KNOWN-VARIANT

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| a_step | 0.784 | 6.505 | -5.721 [-39.2] | -19.13 | -6.438 | 0.651 | yes | R |
| a_ang_mom | 0.0214 | 0.0272 | -0.0058 [-0.4] | 0.0023 | 0.0218 | 0.055 | yes |  |

## census:g2:network-1  (SINGLE MECHANISM)

`N kN=2 | N[any] ALWAYS() -> COPY_RAND() p=0.1` · view N · flagged False · verdict COVERED-BY-TWIN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| n_components | 1.609 | 1.609 | 0 [0.0] | 0 | 0 | 1.666 | yes |  |
| n_giant | 0.985 | 0.985 | 0 [0.0] | 0 | 0 | 0.981 | yes |  |
| n_type_entropy | 0.367 | 0.6919 | -0.3249 [-3.2] | -0.047 | -0.6604 | 0.831 | yes |  |
| n_oscillation | 0.5057 | 0 | 0.5057 [15.1] | 0.5885 | 0.5902 | 0.524 | yes | R |
| n_degree_cv | 0.4947 | 0.4947 | 0 [0.0] | 0 | 0 | 0.49 | yes |  |
| n_activity | 0.1075 | 0 | 0.1075 [7.6] | 0.1933 | 0 | 0.168 | yes |  |

## census:g2:grid-1  (SINGLE MECHANISM)

`C kC=2 | C[any] ALWAYS() -> COPY_RAND() p=1.0` · view C · flagged False · verdict COVERED-BY-TWIN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_n_domains | 5.756 | 6.381 | -0.6252 [-2.6] | -0.5309 | -0.9656 | 6.035 | yes |  |
| g_n_types | 2 | 2 | 0 [0.0] | 0 | 0 | 2.52 | yes |  |
| g_corrlen | 1 | 1 | 0 [0.0] | 1 | 0 | 1.893 | no |  |
| g_type_entropy | 0.9363 | 0.9999 | -0.0636 [-1.6] | -0.0057 | -0.1745 | 0.93 | yes |  |
| g_reversals | 1 | 0 | 1 [50.0] | 1 | 1 | 0.868 | yes | R |
| g_largest_domain | 0.6685 | 0.05 | 0.6184 [26.3] | 0.3197 | 0.7811 | 0.384 | no | R |

## census:g1:agents-3  (SINGLE MECHANISM)

`A kA=1 v=1.0 eta=3.14159 | A.ALIGN(r=1.0)` · view A · flagged True · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| a_step | 6.558 | 7.183 | -0.6244 [-2.3] | -1.002 | 0.1047 | 3.9 | no |  |
| a_density_cv | 1.054 | 0.7899 | 0.2646 [7.2] | 0.297 | 0.2398 | 1.067 | yes | R |
| a_nn_ratio | 0.9286 | 0.9869 | -0.0583 [-1.5] | -0.0934 | -0.0798 | 0.848 | yes |  |
| a_type_entropy | undefined | | | | | 0.54 |  | |
| a_type_segregation | undefined | | | | | 0.225 |  | |
| a_local_align | 0.1717 | -0.0462 | 0.218 [15.8] | 0.1949 | 0.2155 | 0.215 | yes | R |

## census:g1:grid-6  (SINGLE MECHANISM)

`C kC=2 | C.CONTAG(a=0,b=1,th=1,p=0.1)` · view C.transient · flagged True · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_largest_domain | 1 | ko absent | | ko absent | ko absent | 1 | yes | R |
| g_n_domains | 0.6931 | ko absent | | ko absent | ko absent | 0.8 | yes | R |

## census:g1:agents-8  (SINGLE MECHANISM)

`AF kA=1 v=1.0 eta=0.3 fields=1 D=[0.05] | A.CLIMB(f=0)` · view A · flagged True · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| a_density_cv | 13.77 | 0.7849 | 12.98 [7.7] | 6.269 | 13.99 | 13.02 | yes |  |
| a_step | 0.7902 | 19.75 | -18.96 [-46.2] | -19 | -22.25 | 0.576 | yes | R |
| a_type_entropy | undefined | | | | | 0.437 |  | |
| a_type_segregation | undefined | | | | | 0.146 |  | |
| a_polar | 0.0491 | 0.0378 | 0.0113 [1.1] | 0.0097 | -0.0087 | 0.084 | yes |  |
| a_ang_mom | 0.0059 | 0.0209 | -0.015 [-1.6] | -0.0069 | 0.0163 | 0.079 | no |  |

## census:g1:grid-11  (ARTIFACT, flagged view did not recur on fresh seeds)

`C kC=2 | C.COPY(p=1.0) ; C.MAJ(p=1.0)` · view C.transient · flagged True · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_activity_trend | view absent | | | | | -3.231 |  | |
| g_n_types | view absent | | | | | 1.606 |  | |
| g_largest_domain | view absent | | | | | 0.999 |  | |
| g_activity_cv | view absent | | | | | 0.949 |  | |
| g_reversals | view absent | | | | | 0.921 |  | |
| g_n_domains | view absent | | | | | 0.866 |  | |

## census:g1:network-phase-2  (COMBINATION EFFECT)

`N kN=3 | N.KURA(K=1.0) ; N.REWIRE(p=1.0)` · view N.phase · flagged True · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| p_r | 0.8471 | 0.0638 | 0.7833 [15.6] | 0.3826 | 0.4177 | 0.523 | yes |  |
| p_r_std | 0.1157 | 0.0289 | 0.0868 [6.1] | 0.0864 | 0.1611 | 0.176 | yes | R |
| p_freq_spread | 0.0552 | 0.4987 | -0.4435 [-19.2] | -0.4565 | -1.256 | 0.066 | yes | R |

## census:g1:grid-17  (ARTIFACT, flagged view did not recur on fresh seeds)

`C kC=4 | C.MAJ(p=1.0) ; C.MAJ(p=0.1)` · view C.transient · flagged True · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_activity_trend | view absent | | | | | -2.232 |  | |
| g_n_types | view absent | | | | | 2.089 |  | |
| g_largest_domain | view absent | | | | | 0.999 |  | |
| g_n_domains | view absent | | | | | 0.9 |  | |
| g_corrlen | view absent | | | | | 0.821 |  | |
| g_reversals | view absent | | | | | 0.731 |  | |

## census:g2:agents-1  (SINGLE MECHANISM)

`A kA=1 v=1.0 eta=0.3 | A[any] ALWAYS() -> HEAD_TO(b=any,r=1.0) p=1.0` · view A · flagged False · verdict COVERED-BY-TWIN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| a_density_cv | 2.18 | 0.8089 | 1.371 [22.9] | 0.708 | 1.502 | 3.388 | no | R |
| a_step | 0.2803 | 6.536 | -6.255 [-45.9] | -19.48 | -6.907 | 0.909 | no | R |
| a_nn_ratio | 0.1849 | 0.9949 | -0.81 [-34.3] | -0.812 | -0.8296 | 0.243 | yes | R |
| a_local_align | -0.2085 | -0.0461 | -0.1624 [-16.2] | -0.16 | -0.147 | -0.111 | no | R |
| a_polar | 0.0063 | 0.0493 | -0.0429 [-18.7] | -0.0161 | -0.0341 | 0.019 | no |  |
| a_ang_mom | 0.0066 | 0.0201 | -0.0135 [-1.8] | -0.0108 | -0.0076 | 0.015 | yes |  |

## census:novelty:network-2  (COMBINATION EFFECT)

`N kN=2 | N.GAME(g=(1.9, 0.0)) ; N.MAJ(p=0.1)` · view N · flagged False · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| n_components | 1.792 | 1.792 | 0 [0.0] | 0 | 0 | 1.86 | yes |  |
| n_giant | 0.98 | 0.98 | 0 [0.0] | 0 | 0 | 0.975 | yes |  |
| n_type_entropy | 0.657 | 0.6928 | -0.0358 [-1.3] | -0.0267 | -0.0272 | 0.678 | yes |  |
| n_degree_cv | 0.4609 | 0.4609 | 0 [0.0] | 0 | 0 | 0.5 | yes |  |
| n_activity | 0.2895 | 0 | 0.2895 [49.9] | 0.3294 | 0.3777 | 0.329 | yes | R |
| n_type_modularity | 0.167 | -0.0031 | 0.17 [12.4] | 0.1657 | 0.1415 | 0.157 | yes | R |

## census:novelty:grid-1  (DECOMPOSES)

`C kC=2 | C.GAME(g=(1.5, 0.5)) ; C.KURA(K=1.0)` · view C · flagged False · verdict KNOWN · reproduced by rules [1]

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_n_domains | 4.771 | 6.371 | -1.6 [-7.2] | -1.598 | -1.596 | 4.746 | yes | R |
| g_corrlen | 2 | 1 | 1 [16.7] | 1 | 1 | 2 | yes | R |
| g_n_types | 2 | 2 | 0 [0.0] | 0 | 0 | 2 | yes |  |
| g_activity_trend | -0.7781 | 0 | -0.7781 [-40.7] | -0.6662 | -0.0011 | -1.633 | yes |  |
| g_reversals | 1 | 0 | 1 [50.0] | 1 | 1 | 1 | yes | R |
| g_largest_domain | 0.8179 | 0.0481 | 0.7698 [44.5] | 0.7987 | 0.7703 | 0.811 | yes | R |

## census:novelty:network-phase-1  (DECOMPOSES)

`N kN=2 | N.GAME(g=(1.9, 0.0)) ; N.KURA(K=0.1)` · view N.phase · flagged False · verdict KNOWN · reproduced by rules [2]

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| p_freq_spread | 0.347 | 0.5007 | -0.1537 [-3.4] | -0.1448 | -0.5827 | 0.338 | yes | R |
| p_r | 0.235 | 0.0609 | 0.1742 [10.2] | 0.0673 | 0.1739 | 0.25 | yes |  |
| p_r_std | 0.082 | 0.0295 | 0.0525 [5.2] | 0.0187 | 0.0599 | 0.068 | yes |  |

## census:g1:network-7  (ARTIFACT)

`N kN=3 | N.COPY(p=0.1) ; N.REWIRE(p=0.03)` · view N · flagged True · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| n_components | 1.099 | 1.386 | -0.2877 [-0.8] | -1.393 | -0.4904 | 1.106 | yes |  |

## census:novelty:grid-3  (ARTIFACT)

`CA kC=3 kA=3 v=0.0 eta=3.14159 | A.WRITECELL()` · view C · flagged True · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_cyclic_changes | 0.5 | -0.5 | 1 [50.0] | 1 | 0 | 0.5 | yes |  |

## census:g1:field-1  (ARTIFACT)

`F fields=2 D=[0.01, 0.2] | F.AUTOCAT(f=0,g=1,p=0.05)` · view F0 · flagged True · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| f_std | 0.033 | 0.0171 | 0.0159 [4.3] | 0.0021 | -0.0089 | 6.682e+09 | no |  |
| f_mean | 0.9123 | 0.9836 | -0.0713 [-1.9] | -0.0059 | -0.9836 | 1.342e+09 | no |  |
| f_moran | 0.9941 | 0.9852 | 0.0089 [0.2] | 0.0025 | -0.9944 | 0.942 | yes |  |
| f_n_blobs | 1.609 | 1.946 | -0.3365 [-1.3] | 0 | -0.5108 | 0.894 | yes |  |

## census:g1:network-14  (ARTIFACT)

`N kN=2 | N.REWIRE(p=0.1) ; N.CYCLE(th=1,p=0.1)` · view N · flagged True · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| n_components | 1.099 | 1.609 | -0.5108 [-0.8] | -1.153 | -0.9163 | 1.428 | yes |  |
| n_giant | 0.995 | 0.985 | 0.01 [0.3] | 0.0156 | 0.015 | 0.986 | yes |  |

## census:g1:grid-25  (COMBINATION EFFECT)

`C kC=4 | C.SWAP(p=1.0) ; C.CYCLE(th=3,p=0.1)` · view C · flagged True · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_n_domains | 6.504 | 7.669 | -1.165 [-4.1] | -1.177 | -1.303 | 6.347 | yes | R |
| g_n_types | 4 | 4 | 0 [0.0] | 0 | -2 | 3.206 | yes |  |
| g_corrlen | 1 | 1 | 0 [0.0] | 0 | 0 | 1.118 | yes |  |
| g_reversals | 0.7457 | 0 | 0.7457 [50.0] | 0.7084 | 0.8615 | 0.704 | yes | R |
| g_type_entropy | 0.6692 | 0.9998 | -0.3305 [-9.9] | -0.038 | -0.5014 | 0.647 | yes |  |
| g_activity_cv | 0.4376 | 0 | 0.4376 [49.7] | 0.5215 | 0.4455 | 0.42 | yes | R |

## census:g1:network-10  (DECOMPOSES)

`N kN=4 | N.REWIRE(p=0.01) ; N.MAJ(p=1.0)` · view N · flagged False · verdict KNOWN · reproduced by rules [2]

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| n_components | 1.701 | 1.792 | -0.0912 [-0.4] | -0.1038 | -0.0912 | 1.766 | yes |  |
| n_giant | 0.98 | 0.98 | 0 [0.0] | 0 | -0.0025 | 0.892 | yes |  |
| n_type_entropy | 0.0962 | 1.379 | -1.282 [-43.5] | -1.3 | -1.282 | 0.559 | no | R |
| n_degree_cv | 0.4904 | 0.4917 | -0.0013 [-0.1] | -0.0084 | -0.0013 | 0.478 | yes |  |
| n_type_modularity | 0 | -0.0119 | 0.0119 [1.0] | 0.0052 | 0.0119 | 0.339 | no |  |
| n_consensus_gain | 0.7025 | 0 | 0.7025 [50.2] | 0.715 | 0.7025 | 0.261 | no | R |

## census:g2:field-1  (ARTIFACT)

`F fields=2 D=[0.01, 0.2] | F[any] ALWAYS() -> AUTOCAT(f=0,g=1,p=0.05) p=1.0` · view F0 · flagged False · verdict COVERED-BY-TWIN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| f_std | 0.033 | 0.0171 | 0.0159 [4.3] | 0.0021 | -0.0089 | 6.682e+09 | no |  |
| f_mean | 0.9123 | 0.9836 | -0.0713 [-1.9] | -0.0059 | -0.9836 | 1.342e+09 | no |  |
| f_change | 0.0016 | 0.0001 | 0.0014 [1.4] | 0.0001 | -0.0001 | 6.708e+08 | no |  |
| f_moran | 0.9941 | 0.9852 | 0.0089 [0.2] | 0.0025 | -0.9944 | 0.942 | yes |  |
| f_n_blobs | 1.609 | 1.946 | -0.3365 [-1.3] | 0 | -0.5108 | 0.894 | yes |  |
| f_peak_sharp | 0.82 | 0.4789 | 0.3411 [3.5] | -0.0095 | -0.7592 | 0.426 | yes |  |

## census:novelty:network-3  (DECOMPOSES)

`N kN=2 | N.REWIRE(p=0.003) ; N.GAME(g=(1.9, 0.0))` · view N · flagged True · verdict KNOWN · reproduced by rules [1]

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| n_components | 1.099 | 1.609 | -0.5108 [-3.8] | -0.4492 | -0.5108 | 1.38 | yes | R |
| n_giant | 0.995 | 0.9825 | 0.0125 [0.3] | 0.0075 | 0.0125 | 0.988 | yes |  |
| n_type_entropy | 0.5268 | 0.6913 | -0.1645 [-6.7] | -0.0027 | -0.4068 | 0.569 | yes |  |
| n_degree_cv | 0.5254 | 0.508 | 0.0174 [0.8] | 0.0254 | 0.0416 | 0.509 | yes |  |

## census:g1:grid-19  (COMBINATION EFFECT)

`C kC=2 | C.MAJ(p=0.003) ; C.GAME(g=(1.5, -0.5))` · view C · flagged True · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_corrlen | 5.5 | 1 | 4.5 [12.1] | 5 | 8 | 5.081 | yes | R |
| g_activity_cv | 2.285 | 0 | 2.285 [31.8] | 1.206 | 1.189 | 2.367 | yes | R |
| g_largest_domain | 0.9585 | 0.0419 | 0.9166 [45.8] | 0.9383 | 0.7355 | 0.968 | yes | R |

## census:novelty:grid-4  (COMBINATION EFFECT)

`C kC=2 | C.GAME(g=(1.3, 0.0)) ; C.COPY(p=0.03)` · view C.transient · flagged True · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_n_types | 1 | ko absent | | ko absent | ko absent | 1 | yes | R |
| g_largest_domain | 1 | ko absent | | ko absent | ko absent | 1 | yes | R |

## census:novelty:lattice-phase-1  (ARTIFACT)

`C kC=2 | C.KURA(K=0.01) ; C.GAME(g=(0.5, -0.5))` · view C.phase · flagged False · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| p_local_r | 0.4618 | 0.451 | 0.0108 [0.6] | 0.0118 | 0.0108 | 0.495 | yes |  |
| p_freq_spread | 0.497 | 0.4982 | -0.0013 [-0.1] | -0.0013 | -0.0301 | 0.482 | yes |  |
| p_vortex_density | 0.3169 | 0.3318 | -0.0149 [-1.1] | -0.0143 | -0.0242 | 0.267 | yes |  |
| p_r | 0.0169 | 0.0162 | 0.0007 [0.3] | 0.0003 | 0.0018 | 0.019 | yes |  |
| p_r_std | 0.008 | 0.0078 | 0.0002 [0.1] | 0.0003 | 0.0006 | 0.009 | yes |  |

## census:g1:field-2  (ARTIFACT)

`F fields=2 D=[0.2, 0.05] | F.FEED(f=0,F=1.0) ; F.AUTOCAT(f=0,g=1,p=1.0)` · view F0 · flagged True · verdict KNOWN-VARIANT

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| f_std | 2.631e+08 | 0 | 2.631e+08 [2.6] | 2.318e+07 | 2.631e+08 | 7.193e+06 | yes |  |
| f_mean | -1.697e+07 | 1 | -1.697e+07 [-3.3] | -1.722e+06 | -1.697e+07 | -5.672e+05 | yes |  |
| f_change | 1.697e+07 | 0 | 1.697e+07 [3.3] | 1.722e+06 | 1.697e+07 | 5.166e+05 | yes |  |

## census:g1:network-13  (DECOMPOSES)

`N kN=2 | N.KURA(K=0.03) ; N.REWIRE(p=0.01)` · view N · flagged False · verdict KNOWN · reproduced by rules [2]

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| n_components | 1.386 | 1.701 | -0.3143 [-1.3] | -0.463 | -0.0912 | 1.332 | yes |  |
| n_type_entropy | 0.6881 | 0.6881 | 0 [0.0] | 0 | 0 | 0.936 | no |  |
| n_giant | 0.99 | 0.9825 | 0.0075 [0.2] | 0.0075 | -0.435 | 0.86 | yes |  |
| n_type_modularity | 0.4884 | -0.0193 | 0.5077 [31.0] | 0.4852 | 0.5149 | 0.575 | yes | R |
| n_edge_turnover | 0.5007 | 0 | 0.5007 [50.1] | 0.4843 | 0.5148 | 0.5 | yes | R |
| n_degree_cv | 0.4341 | 0.479 | -0.045 [-1.3] | -0.0408 | -0.0476 | 0.447 | yes |  |

## census:g1:grid-29  (ARTIFACT, flagged view did not recur on fresh seeds)

`C kC=4 | C.CYCLE(th=3,p=0.3) ; C.MAJ(p=1.0)` · view C.transient · flagged True · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_activity_trend | view absent | | | | | -4.742 |  | |

## census:g1:grid-26  (ARTIFACT)

`C kC=4 | C.SAND(p=0.002) ; C.SCHELL(th=0.5)` · view C · flagged True · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_n_domains | 7.261 | 7.665 | -0.4043 [-1.4] | -0.332 | -0.4128 | 7.403 | yes |  |

## census:novelty:agents-6  (ARTIFACT)

`AF kA=1 v=1.0 eta=1.0 fields=1 D=[0.05] | F.FEED(f=0,F=0.1) ; A.CLIMB(f=0)` · view A · flagged False · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| a_step | 6.686 | 17.46 | -10.77 [-22.3] | -3.029 | -19.26 | 8.409 | yes |  |
| a_density_cv | 3.005 | 0.7899 | 2.215 [15.6] | 0.3799 | 2.116 | 2.304 | yes |  |
| a_nn_ratio | 0.3428 | 0.9953 | -0.6525 [-21.3] | -0.1885 | -0.6862 | 0.547 | no |  |
| a_type_entropy | undefined | | | | | 0.111 |  | |
| a_polar | 0.0423 | 0.0435 | -0.0012 [-0.3] | -0.0023 | -0.0007 | 0.045 | yes |  |
| a_local_align | -0.0199 | -0.0322 | 0.0123 [1.2] | -0.0043 | 0.0048 | -0.025 | yes |  |

## census:novelty:agents-17  (COMBINATION EFFECT)

`AF kA=1 v=1.0 eta=1.0 fields=1 D=[0.05] | A.CLIMB(f=0) ; A.REPEL(b=any,r=2.0)` · view A · flagged False · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| a_step | 2.342 | 17.44 | -15.1 [-38.2] | -15 | -23.26 | 3.09 | yes | R |
| a_density_cv | 1.179 | 0.792 | 0.3865 [9.8] | 0.6708 | 0.4286 | 1.637 | no | R |
| a_nn_ratio | 0.8208 | 0.9915 | -0.1707 [-4.7] | -0.203 | -0.2331 | 0.802 | yes | R |
| a_type_entropy | undefined | | | | | 0.711 |  | |
| a_type_segregation | undefined | | | | | 0.216 |  | |
| a_local_align | -0.1074 | -0.0265 | -0.0809 [-4.1] | -0.0703 | -0.0788 | -0.056 | yes | R |

## census:novelty:network-5  (DECOMPOSES)

`N kN=2 | N.KURA(K=0.1) ; N.COPY(p=1.0) ; N.GAME(g=(1.5, 0.5))` · view N · flagged False · verdict KNOWN · reproduced by rules [2, 3]

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| n_components | 1.386 | 1.386 | 0 [0.0] | 0 | 0 | 1.69 | yes |  |
| n_giant | 0.99 | 0.99 | 0 [0.0] | 0 | 0 | 0.98 | yes |  |
| n_degree_cv | 0.4711 | 0.4711 | 0 [0.0] | 0 | 0 | 0.483 | yes |  |
| n_type_entropy | 0.615 | 0.6929 | -0.0779 [-3.0] | -0.0673 | -0.0953 | 0.409 | yes |  |
| n_activity | 0.4017 | 0 | 0.4017 [36.5] | 0.4107 | 0.4016 | 0.272 | yes | R |
| n_consensus_gain | 0.1775 | 0 | 0.1775 [19.1] | 0.1612 | 0.2075 | 0.258 | yes | R |

## census:novelty:agents-8  (COMBINATION EFFECT)

`AF kA=1 v=1.0 eta=1.0 fields=1 D=[0.05] | A.CLIMB(f=0) ; A.ALIGN(r=2.0)` · view A · flagged True · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| a_density_cv | 10.97 | 0.8149 | 10.16 [38.8] | 6.338 | 11.66 | 8.021 | no | R |
| a_step | 0.7435 | 17.48 | -16.73 [-45.9] | -16.73 | -26.34 | 1.193 | no | R |
| a_local_align | 0.8728 | -0.0362 | 0.9091 [50.0] | 0.92 | 0.929 | 0.827 | yes | R |
| a_ang_mom | 0.1631 | 0.0302 | 0.133 [2.6] | 0.2469 | -0.0078 | 0.52 | yes |  |
| a_polar | 0.7327 | 0.0473 | 0.6854 [20.5] | 0.3246 | 0.9094 | 0.391 | yes |  |
| a_type_entropy | undefined | | | | | 0.187 |  | |

## census:novelty:network-4  (DECOMPOSES)

`N kN=2 | N.COPY(p=1.0) ; N.REWIRE(p=0.3) ; N.GAME(g=(1.9, 0.0))` · view N · flagged True · verdict KNOWN · reproduced by rules [1, 2]

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| n_components | 1.386 | 1.609 | -0.2231 [-3.7] | -0.3925 | -0.2231 | 1.456 | yes | R |
| n_giant | 0.99 | 0.985 | 0.005 [0.1] | 0.0069 | 0.0025 | 0.986 | yes |  |
| n_edge_turnover | 0.9812 | 0 | 0.9812 [50.1] | 0.9949 | 0.9861 | 0.541 | no | R |
| n_degree_cv | 0.4736 | 0.4858 | -0.0122 [-0.4] | -0.0077 | -0.0123 | 0.511 | yes |  |
| n_oscillation | 0.3107 | 0 | 0.3107 [28.5] | 0.2495 | 0.2978 | 0.387 | yes | R |
| n_type_entropy | 0.5805 | 0.6909 | -0.1104 [-2.3] | -0.0902 | -0.1678 | 0.36 | yes |  |

## census:novelty:field-1  (COMBINATION EFFECT)

`AF kA=1 v=1.0 eta=1.0 fields=1 D=[0.02] | A.ALIGN(r=8.0) ; A.DEPOSIT(f=0,p=0.1)` · view F0 · flagged False · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| f_mean | 10.75 | 10.75 | 0 [0.0] | 0 | -0 | 28.89 | no |  |
| f_std | 0.6339 | 0.261 | 0.3729 [11.9] | 1.049 | 0.3167 | 17.78 | no | R |
| f_n_blobs | 2.833 | 3.998 | -1.165 [-4.8] | -2.035 | -0.8343 | 1.49 | no | R |
| f_moran | 0.9802 | 0.9107 | 0.0696 [1.8] | 0.0683 | 0.0477 | 0.992 | yes |  |
| f_peak_sharp | 0.8488 | 0.395 | 0.4538 [4.1] | 0.4287 | 0.3433 | 0.772 | yes | R |
| f_change | 0.1042 | 0.1034 | 0.0007 [0.2] | 0.0011 | 0.0007 | 0.285 | no |  |

## census:novelty:agents-14  (DECOMPOSES)

`AF kA=2 v=0.5 eta=1.0 fields=1 D=[0.05] | A.ALIGN(r=1.0) ; A.DEPOSIT(f=0,p=1.0)` · view A · flagged True · verdict KNOWN · reproduced by rules [1]

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| a_step | 8.743 | 8.733 | 0.0104 [0.0] | -0.0751 | 0.2209 | 11.91 | no |  |
| a_density_cv | 1.688 | 0.7909 | 0.8967 [18.1] | 0.7489 | 0.8777 | 2.545 | no | R |
| a_type_entropy | 0.6925 | 0.6925 | 0 [0.0] | 0 | 0 | 1.095 | no |  |
| a_nn_ratio | 0.5712 | 1.008 | -0.4366 [-13.8] | -0.4188 | -0.4691 | 0.607 | yes | R |
| a_type_segregation | 0.4972 | 0.4969 | 0.0003 [0.0] | 0.0005 | -0.006 | 0.348 | no |  |
| a_local_align | 0.3438 | -0.0322 | 0.376 [22.2] | 0.3961 | 0.4032 | 0.204 | no | R |

## census:novelty:field-2  (ARTIFACT)

`AF kA=1 v=1.0 eta=1.0 fields=1 D=[0.05] | A.ATTRACT(b=any,r=0.5) ; A.DEPOSIT(f=0,p=0.2)` · view F0 · flagged False · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| f_mean | 20.51 | 20.51 | 0 [0.0] | -0 | -0 | 18.89 | yes |  |
| f_n_blobs | 3.689 | 3.511 | 0.1774 [0.7] | -0.2635 | -0.2782 | 3.037 | yes |  |
| f_moran | 0.9385 | 0.9504 | -0.0119 [-0.3] | 0.0108 | 0.0082 | 0.96 | yes |  |
| f_peak_sharp | 0.4747 | 0.5647 | -0.0901 [-0.5] | 0.0683 | 0.0974 | 0.743 | yes |  |
| f_std | 0.3828 | 0.4151 | -0.0323 [-0.6] | 0.1022 | 0.0562 | 0.714 | no |  |
| f_change | 0.2029 | 0.2026 | 0.0002 [0.0] | 0.0001 | 0 | 0.189 | yes |  |

## census:novelty:lattice-phase-2  (DECOMPOSES)

`C kC=2 | C.CYCLE(th=4,p=0.01) ; C.GAME(g=(1.5, 0.5)) ; C.KURA(K=0.1)` · view C.phase · flagged False · verdict KNOWN · reproduced by rules [3]

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| p_local_r | 0.6682 | 0.4489 | 0.2193 [9.8] | 0.2226 | 0.225 | 0.668 | yes | R |
| p_freq_spread | 0.3178 | 0.4995 | -0.1816 [-11.1] | -0.1782 | -0.7271 | 0.322 | yes | R |
| p_vortex_density | 0.1045 | 0.3354 | -0.231 [-26.2] | -0.2294 | -0.228 | 0.108 | yes | R |
| p_r | 0.0534 | 0.0141 | 0.0393 [5.6] | 0.0212 | 0.0346 | 0.055 | yes | R |
| p_r_std | 0.0152 | 0.0073 | 0.0079 [5.3] | 0.0034 | 0.0164 | 0.022 | no |  |

## census:novelty:agents-16  (COMBINATION EFFECT)

`CAF kC=1 kA=4 v=1.0 eta=0.3 fields=2 D=[0.02, 0.02] | A.CLIMB(f=1) ; A.TURNCELL()` · view A · flagged True · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| a_density_cv | 1.843 | 0.8 | 1.043 [19.7] | 1.065 | 1.077 | 1.598 | yes | R |

## census:novelty:grid-7  (COMBINATION EFFECT)

`CAF kC=2 kA=1 v=1.0 eta=1.0 fields=2 D=[0.02, 0.1] | C.MAJ(p=0.1) ; A.WRITECELL()` · view C · flagged True · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_n_domains | 6.099 | 6.35 | -0.2508 [-1.0] | -0.2661 | -0.2679 | 6.059 | yes |  |
| g_n_types | 2 | 2 | 0 [0.0] | 0 | 0 | 2 | yes |  |
| g_corrlen | 1 | 1 | 0 [0.0] | 0 | 0 | 1 | yes |  |
| g_reversals | 1 | 0 | 1 [50.0] | 1 | 1 | 1 | yes | R |
| g_type_entropy | 1 | 0.9999 | 0.0001 [0.0] | 0 | 0.0001 | 1 | yes |  |
| g_interface | 0.4315 | 0.498 | -0.0665 [-3.6] | -0.0682 | -0.0676 | 0.427 | yes | R |

## census:novelty:grid-8  (COMBINATION EFFECT)

`CAF kC=2 kA=1 v=1.0 eta=1.0 fields=2 D=[0.02, 0.1] | A.WRITECELL() ; A.CLIMB(f=0)` · view C · flagged True · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_n_domains | 6.377 | 6.394 | -0.0178 [-0.1] | -0.0173 | -0.0212 | 6.175 | yes |  |
| g_n_types | 2 | 2 | 0 [0.0] | 0 | 0 | 2 | yes |  |
| g_activity_trend | -3.288 | 0 | -3.288 [-50.0] | -3.544 | -3.434 | -1.206 | no | R |
| g_corrlen | 1 | 1 | 0 [0.0] | 0 | 0 | 1.111 | yes |  |
| g_reversals | 1 | 0 | 1 [50.0] | 1 | 1 | 1 | yes | R |
| g_type_entropy | 1 | 0.9999 | 0.0001 [0.0] | 0 | 0.0001 | 0.998 | yes |  |

## census:novelty:grid-9  (DECOMPOSES)

`C kC=2 | C.CYCLE(th=4,p=0.01) ; C.GAME(g=(1.5, -0.5)) ; C.KURA(K=0.0003)` · view C · flagged False · verdict KNOWN · reproduced by rules [1, 2]

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_corrlen | 8 | 1 | 7 [9.4] | 10.5 | 1 | 9.667 | yes |  |
| g_n_domains | 2.832 | 6.375 | -3.543 [-19.2] | -4.317 | -2.483 | 2.602 | yes | R |
| g_n_types | 2 | 2 | 0 [0.0] | 0 | 0 | 2 | yes |  |
| g_activity_trend | 2.307 | 0 | 2.307 [24.6] | 2.032 | 3.045 | 1.914 | yes | R |
| g_reversals | 1 | 0 | 1 [50.0] | 1 | 1 | 1 | yes | R |
| g_moran | 0.8303 | 0.0024 | 0.8279 [32.1] | 0.8595 | 0.5577 | 0.846 | yes | R |

## census:novelty:field-3  (COMBINATION EFFECT)

`CAF kC=2 kA=1 v=1.0 eta=1.0 fields=2 D=[0.02, 0.1] | C.MAJ(p=0.1) ; C.EMIT(b=1,f=0,p=1.0)` · view F0 · flagged True · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| f_mean | 501.8 | 501 | 0.8227 [0.0] | 4.762 | -27.83 | 539.2 | yes |  |
| f_std | 289.9 | 33.84 | 256.1 [39.5] | 273.7 | 693.1 | 300.8 | yes | R |
| f_change | 4.967 | 5 | -0.0332 [-0.2] | 0.1174 | -0.1513 | 5.616 | yes |  |
| f_n_blobs | 1.386 | 3.135 | -1.749 [-7.5] | -2.291 | -1.426 | 1.164 | yes | R |

## census:novelty:grid-10  (DECOMPOSES)

`CAF kC=2 kA=1 v=1.0 eta=1.0 fields=2 D=[0.02, 0.1] | C.MAJ(p=0.0003) ; A.CLIMB(f=0)` · view C · flagged True · verdict KNOWN · reproduced by rules [1]

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_n_domains | 6.092 | 6.366 | -0.2732 [-1.1] | -0.2915 | -0.8342 | 6.209 | yes |  |
| g_n_types | 2 | 2 | 0 [0.0] | 0 | 0 | 2.182 | yes |  |
| g_corrlen | 1 | 1 | 0 [0.0] | 0 | 1 | 1 | yes |  |
| g_type_entropy | 0.9998 | 0.9999 | -0.0001 [-0.0] | -0 | -0.001 | 1 | yes |  |
| g_activity_cv | 1.054 | 0 | 1.054 [49.9] | 0.5752 | 1.283 | 0.967 | yes | R |
| g_interface | 0.4366 | 0.4995 | -0.063 [-3.4] | -0.0663 | -0.1649 | 0.455 | yes | R |

## census:novelty:grid-12  (DECOMPOSES)

`CAF kC=2 kA=1 v=1.0 eta=1.0 fields=2 D=[0.02, 0.1] | C.COPY(p=0.001) ; A.CLIMB(f=0)` · view C · flagged False · verdict KNOWN · reproduced by rules [1]

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_n_domains | 6.099 | 6.379 | -0.2801 [-1.1] | -0.2892 | -0.5682 | 6.099 | yes |  |
| g_n_types | 2 | 2 | 0 [0.0] | 0 | 0 | 2 | yes |  |
| g_corrlen | 1 | 1 | 0 [0.0] | 0 | 0 | 1 | yes |  |
| g_reversals | 1 | 0 | 1 [50.0] | 1 | 1 | 1 | yes | R |
| g_type_entropy | 0.9995 | 0.9999 | -0.0003 [-0.0] | -0 | -0.0005 | 0.999 | yes |  |
| g_activity_cv | 0.4987 | 0 | 0.4987 [49.9] | 0.2464 | 0.5912 | 0.471 | yes |  |

## census:novelty:grid-15  (DECOMPOSES)

`C kC=2 | C.CONTAG(a=1,b=0,th=6,p=0.005) ; C.GAME(g=(1.9, 0.0)) ; C.KURA(K=0.0003)` · view C · flagged False · verdict KNOWN · reproduced by rules [2]

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_n_domains | 4.303 | 6.399 | -2.095 [-9.8] | -2.048 | -2.048 | 3.225 | no | R |
| g_corrlen | 2 | 1 | 1 [16.7] | 1 | 1 | 2.556 | yes | R |
| g_n_types | 2 | 2 | 0 [0.0] | 0 | 0 | 1.889 | yes |  |
| g_reversals | 1 | 0 | 1 [50.0] | 1 | 1 | 0.875 | yes | R |
| g_largest_domain | 0.6909 | 0.0453 | 0.6456 [34.7] | 0.6592 | 0.6564 | 0.849 | yes | R |
| g_moran | 0.5994 | -0.0078 | 0.6072 [32.0] | 0.6127 | 0.5993 | 0.452 | no | R |

## census:novelty:grid-11  (ARTIFACT)

`CAF kC=2 kA=1 v=1.0 eta=1.0 fields=2 D=[0.02, 0.1] | C.CYCLE(th=7,p=0.03) ; A.CLIMB(f=0)` · view C · flagged True · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_activity_trend | -1.914 | 0 | -1.914 [-50.0] | -1.883 | -0.9341 | -1.859 | yes |  |

## census:novelty:grid-13  (DECOMPOSES)

`CAF kC=4 kA=1 v=0.0 eta=2.0 fields=2 D=[0.02, 0.1] | C.COPY(p=1.0) ; C.EMIT(b=3,f=0,p=0.3)` · view C · flagged True · verdict KNOWN · reproduced by rules [1]

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_n_domains | 6.676 | 7.663 | -0.9878 [-3.4] | -1.072 | -1.764 | 6.946 | yes | R |
| g_n_types | 4 | 4 | 0 [0.0] | 0 | -1 | 3.933 | yes |  |
| g_corrlen | 1 | 1 | 0 [0.0] | 1 | 0 | 1.133 | yes |  |
| g_type_entropy | 0.9294 | 0.9998 | -0.0704 [-1.8] | -0.0369 | -0.4414 | 0.924 | yes |  |
| g_interface | 0.4608 | 0.7482 | -0.2874 [-11.9] | -0.2952 | -0.4568 | 0.545 | yes | R |
| g_reversals | 0.5206 | 0 | 0.5206 [50.1] | 0.5287 | 0.557 | 0.524 | yes | R |

## census:novelty:field-4  (DECOMPOSES)

`CAF kC=4 kA=1 v=0.0 eta=2.0 fields=2 D=[0.02, 0.1] | A.ALIGN(r=8.0) ; C.MAJ(p=0.02) ; C.SWAP(p=0.1) ; C.EMIT(b=3,f=0,p=0.3)` · view F0 · flagged False · verdict KNOWN · reproduced by rules [2, 4]

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| f_mean | 77.37 | 74.7 | 2.673 [0.7] | -0.641 | 25.07 | 88.94 | yes |  |
| f_std | 67.71 | 9.869 | 57.84 [37.3] | 56.3 | 166.2 | 56.15 | yes | R |
| f_n_blobs | 1.792 | 3.461 | -1.67 [-6.7] | -2.273 | -1.735 | 1.753 | yes | R |
| f_moran | 0.9912 | 0.9502 | 0.041 [1.1] | 0.0455 | 0.0187 | 0.99 | yes |  |
| f_change | 0.8195 | 0.7372 | 0.0823 [1.6] | 0.0331 | 0.1919 | 0.898 | yes |  |
| f_peak_sharp | 0.5362 | 0.3913 | 0.145 [1.9] | 0.2481 | 0.3998 | 0.604 | yes |  |

## review:001  (ARTIFACT)

`N kN=4 | N.COPY(p=1.0)` · view N · flagged True · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| n_components | 1.498 | 1.498 | 0 [0.0] | 0 | 0 | 1.551 | yes |  |
| n_giant | 0.9875 | 0.9875 | 0 [0.0] | 0 | 0 | 0.985 | yes |  |
| n_oscillation | 0.433 | 0 | 0.433 [16.7] | 0.6594 | 0.1676 | 0.4918 | yes |  |

## review:002  (SINGLE MECHANISM)

`C kC=1 | C[any] ALWAYS() -> PHASE_COUPLE(K=1.0) p=0.1` · view C.phase · flagged False · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| p_local_r | 0.971 | 0.4498 | 0.5212 [18.4] | 0.5242 | 0.5352 | 0.9728 | yes | R |
| p_r | 0.1126 | 0.013 | 0.0996 [12.1] | 0.1008 | 0.2261 | 0.1187 | yes | R |
| p_freq_spread | 0.0104 | 0.4987 | -0.4883 [-47.9] | -0.4913 | -1.403 | 0.01137 | yes | R |
| p_r_std | 0.0085 | 0.0067 | 0.0018 [0.9] | 0.0015 | 0.004 | 0.00699 | yes |  |
| p_vortex_density | 0.0051 | 0.3313 | -0.3262 [-48.7] | -0.3278 | -0.3325 | 0.00488 | yes | R |

## review:003  (SINGLE MECHANISM)

`C kC=3 | C.MAJ(p=1.0)` · view C · flagged False · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_corrlen | 6.5 | 1 | 5.5 [14.8] | 24 | -1 | 13.5 | yes |  |
| g_n_types | 3 | 3 | 0 [0.0] | 0 | -2 | 3 | yes |  |
| g_n_domains | 1.242 | 7.348 | -6.106 [-35.5] | -6.328 | -6.655 | 1.648 | yes | R |
| g_activity_trend | -1.444 | 0 | -1.444 [-24.8] | -1.494 | -4.536 | -1.2 | yes | R |
| g_moran | 0.878 | 0.0019 | 0.8761 [49.8] | 0.9616 | -0.0019 | 0.9351 | yes |  |
| g_reversals | 0.8808 | 0 | 0.8808 [50.0] | 0.8659 | 0.7308 | 0.8923 | yes | R |

## review:004  (SINGLE MECHANISM)

`N kN=3 | N.REWIRE(p=1.0)` · view N · flagged False · verdict KNOWN-VARIANT

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| n_components | 1.609 | 1.609 | 0 [0.0] | -0.532 | 0 | 1.647 | yes |  |
| n_type_entropy | 1.093 | 1.093 | 0 [0.0] | 0 | 0 | 1.094 | yes |  |
| n_type_modularity | 0.6615 | 0.0058 | 0.6556 [49.3] | 0.6662 | 0.6556 | 0.6629 | yes | R |
| n_edge_turnover | 0.6581 | 0 | 0.6581 [49.9] | 0.6661 | 0.6581 | 0.6499 | yes | R |
| n_cyclic_changes | -0.5 | -0.5 | 0 [0.0] | 0 | 0 | -0.5 | yes |  |
| n_degree_cv | 0.422 | 0.504 | -0.0819 [-4.4] | -0.068 | -0.0819 | 0.4117 | yes | R |

## review:005  (ARTIFACT)

`C kC=3 | C.SCHELL(th=0.2)` · view C · flagged False · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_n_domains | 7.165 | 7.352 | -0.1863 [-0.6] | -0.1798 | -0.1863 | 7.137 | yes |  |
| g_n_types | 3 | 3 | 0 [0.0] | 0 | 0 | 3 | yes |  |
| g_corrlen | 1 | 1 | 0 [0.0] | 0 | 0 | 1 | yes |  |
| g_type_entropy | 0.9998 | 0.9998 | 0 [0.0] | 0 | 0 | 0.9999 | yes |  |
| g_interface | 0.6259 | 0.6685 | -0.0426 [-1.6] | -0.0392 | -0.0426 | 0.6188 | yes |  |
| g_activity_trend | 0 | 0 | 0 [0.0] | 0 | 0 | -0.09394 | no |  |

## review:006  (SINGLE MECHANISM)

`N kN=1 | N.KURA(K=0.2)` · view N.phase · flagged False · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| p_r | 0.8611 | 0.064 | 0.7971 [36.2] | 0.8013 | 0.8008 | 0.869 | yes | R |
| p_freq_spread | 0.0688 | 0.5054 | -0.4366 [-16.5] | -0.4219 | -1.236 | 0.04109 | yes | R |
| p_r_std | 0.0062 | 0.0281 | -0.0219 [-5.6] | -0.0089 | -0.0225 | 0.004543 | yes |  |

## review:007  (DECOMPOSES)

`C kC=1 | C.KURA(K=1.0) ; C.KURA(K=1.0)` · view C.phase · flagged True · verdict KNOWN · reproduced by rules [1]

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| p_local_r | 0.9964 | 0.4484 | 0.548 [19.0] | 0.5435 | 0.5447 | 0.997 | yes | R |
| p_r | 0.7715 | 0.0138 | 0.7576 [28.3] | 0.3147 | 0.7775 | 0.8609 | yes |  |
| p_r_std | 0.0065 | 0.0073 | -0.0008 [-0.7] | 0.0143 | -0.007 | 0.00389 | yes |  |
| p_freq_spread | 0.0104 | 0.9948 | -0.9844 [-49.0] | -0.9768 | -1.795 | 0.00264 | yes | R |
| p_vortex_density | 0 | 0.3337 | -0.3337 [-49.8] | -0.333 | -0.3308 | 0 | yes | R |

## review:008  (SINGLE MECHANISM)

`C kC=3 | C.CYCLE(th=3,p=0.1)` · view C · flagged True · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_corrlen | 5 | 1 | 4 [33.3] | 3.5 | 3 | 4.667 | yes | R |
| g_n_domains | 4.102 | 7.355 | -3.252 [-14.2] | -3.378 | -3.228 | 4.215 | yes | R |
| g_n_types | 3 | 3 | 0 [0.0] | 0 | 0 | 3 | yes |  |
| g_type_entropy | 0.9988 | 0.9997 | -0.0009 [-0.0] | -0.0003 | -0.0005 | 0.9975 | yes |  |
| g_moran | 0.7923 | 0.0024 | 0.7899 [49.7] | 0.7744 | 0.7787 | 0.771 | yes | R |
| g_cyclic_changes | 0.4877 | -0.5 | 0.9877 [49.9] | 0.9878 | 0.9801 | 0.4877 | yes | R |

## review:009  (SINGLE MECHANISM)

`N kN=1 | N[any] ALWAYS() -> PHASE_COUPLE(K=0.2) p=0.1` · view N.phase · flagged True · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| p_r | 0.4459 | 0.0665 | 0.3794 [18.1] | 0.2922 | 0.7323 | 0.6454 | yes | R |
| p_r_std | 0.0269 | 0.028 | -0.0011 [-0.1] | 0.0276 | -0.0187 | 0.05675 | yes |  |
| p_freq_spread | 0.0251 | 0.5146 | -0.4895 [-45.3] | -0.4672 | -1.38 | 0.02483 | yes | R |

## review:010  (SINGLE MECHANISM)

`C kC=2 | C.COPY(p=0.005)` · view C · flagged False · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_n_domains | 5.712 | 6.374 | -0.6621 [-2.7] | -0.6757 | -0.8174 | 5.799 | yes |  |
| g_n_types | 2 | 2 | 0 [0.0] | 0 | 0 | 2 | yes |  |
| g_corrlen | 1 | 1 | 0 [0.0] | 0 | 1 | 1 | yes |  |
| g_reversals | 1 | 0 | 1 [50.0] | 1 | 1 | 1 | yes | R |
| g_type_entropy | 0.9996 | 0.9999 | -0.0002 [-0.0] | -0.0002 | -0.002 | 0.9992 | yes |  |
| g_interface | 0.3411 | 0.5015 | -0.1603 [-9.5] | -0.1583 | -0.2052 | 0.3434 | yes | R |

## review:011  (SINGLE MECHANISM)

`N kN=2 | N[0] ALWAYS() -> REWIRE_SAME() p=1.0` · view N · flagged True · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| n_components | 2.917 | 1.498 | 1.419 [4.4] | 1.506 | 1.419 | 2.512 | yes | R |

## review:012  (SINGLE MECHANISM)

`C kC=4 | C.CYCLE(th=1,p=0.3)` · view C · flagged True · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_n_domains | 7.37 | 7.665 | -0.2943 [-1.0] | -0.3072 | -0.3058 | 7.362 | yes |  |
| g_n_types | 4 | 4 | 0 [0.0] | 0 | 0 | 4 | yes |  |
| g_corrlen | 1 | 1 | 0 [0.0] | 0 | 0 | 1 | yes |  |
| g_type_entropy | 1 | 0.9998 | 0.0002 [0.0] | 0 | 0.0002 | 1 | yes |  |
| g_interface | 0.6521 | 0.7493 | -0.0972 [-3.5] | -0.0982 | -0.0961 | 0.6531 | yes | R |
| g_activity | 0.4174 | 0 | 0.4174 [50.3] | 0.4177 | 0.4177 | 0.4172 | yes | R |

## review:013  (SINGLE MECHANISM)

`A kA=1 v=1.0 eta=1.0 | A.ALIGN(r=0.5)` · view A · flagged False · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| a_step | 7.515 | 7.54 | -0.0247 [-0.1] | 0.6339 | 0.8307 | 7.12 | yes |  |
| a_density_cv | 1.169 | 0.7909 | 0.3786 [8.0] | 0.3961 | 0.3811 | 1.177 | yes | R |
| a_polar | 0.8539 | 0.0432 | 0.8107 [45.3] | 0.8277 | 0.8134 | 0.8554 | yes | R |
| a_nn_ratio | 0.8037 | 0.9981 | -0.1944 [-5.4] | -0.1981 | -0.1747 | 0.8206 | yes | R |
| a_local_align | 0.5571 | -0.0351 | 0.5922 [44.2] | 0.6046 | 0.6017 | 0.5528 | yes | R |
| a_ang_mom | 0.0814 | 0.0239 | 0.0575 [3.0] | 0.0231 | 0.0203 | 0.09767 | yes |  |

## review:014  (SINGLE MECHANISM)

`N kN=2 | N.CONTAG(a=0,b=1,th=3,p=0.1)` · view N · flagged False · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| n_components | 1.498 | 1.498 | 0 [0.0] | 0 | 0 | 1.521 | yes |  |
| n_giant | 0.9875 | 0.9875 | 0 [0.0] | 0 | 0 | 0.9867 | yes |  |
| n_degree_cv | 0.4956 | 0.4956 | 0 [0.0] | 0 | 0 | 0.5003 | yes |  |
| n_type_entropy | 0.427 | 0.6907 | -0.2637 [-11.8] | -0.2059 | -0.2637 | 0.4736 | yes | R |
| n_consensus_gain | 0.3075 | 0 | 0.3075 [49.6] | 0.305 | 0.3075 | 0.2767 | yes | R |
| n_oscillation | 0.0652 | 0 | 0.0652 [50.2] | 0.0688 | 0.0268 | 0.06346 | yes |  |

## review:015  (SINGLE MECHANISM)

`N kN=3 | N.CYCLE(th=2,p=0.1)` · view N · flagged True · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| n_components | 1.701 | 1.701 | 0 [0.0] | 0 | 0 | 1.33 | yes |  |
| n_giant | 0.9825 | 0.9825 | 0 [0.0] | 0 | 0 | 0.99 | yes |  |
| n_degree_cv | 0.5013 | 0.5013 | 0 [0.0] | 0 | 0 | 0.4705 | yes |  |
| n_cyclic_changes | 0.443 | -0.5 | 0.943 [49.9] | 0.9396 | 0.8619 | 0.4429 | yes | R |

## review:016  (SINGLE MECHANISM)

`CF kC=2 fields=1 D=[0.2] | C.GAME(g=(1.9, 0.0))` · view C · flagged True · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_n_domains | 4.35 | 6.34 | -1.991 [-9.3] | -2.121 | -2.044 | 4.337 | yes | R |
| g_corrlen | 2 | 1 | 1 [16.7] | 1 | 1 | 2.5 | yes | R |
| g_n_types | 2 | 2 | 0 [0.0] | 0 | 0 | 2 | yes |  |
| g_reversals | 1 | 0 | 1 [50.0] | 1 | 1 | 1 | yes | R |
| g_type_entropy | 0.8945 | 0.9999 | -0.1054 [-2.8] | -0.1011 | -0.1037 | 0.9018 | yes |  |
| g_largest_domain | 0.6847 | 0.0475 | 0.6372 [29.8] | 0.6621 | 0.6554 | 0.6738 | yes | R |

## review:017  (SINGLE MECHANISM)

`C kC=2 | C[any] ALWAYS() -> IMITATE_BEST(g=(1.5, 0.5)) p=0.1` · view C · flagged True · verdict KNOWN-VARIANT

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_n_domains | 6.217 | 6.373 | -0.1567 [-0.6] | -0.1949 | -0.2037 | 6.24 | yes |  |
| g_n_types | 2 | 2 | 0 [0.0] | 0 | 0 | 2 | yes |  |
| g_corrlen | 1 | 1 | 0 [0.0] | 0 | 0 | 1 | yes |  |
| g_reversals | 1 | 0 | 1 [50.0] | 1 | 1 | 1 | yes | R |
| g_type_entropy | 0.9954 | 0.9999 | -0.0045 [-0.1] | -0.0042 | -0.0043 | 0.9955 | yes |  |
| g_interface | 0.4448 | 0.501 | -0.0562 [-3.0] | -0.0533 | -0.055 | 0.4428 | yes |  |

## review:018  (COMBINATION EFFECT)

`N kN=3 | N.COPY(p=1.0) ; N.REWIRE(p=1.0)` · view N · flagged False · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| n_components | 1.609 | 1.609 | 0 [0.0] | -0.857 | 0 | 1.461 | yes |  |
| n_type_entropy | 1.073 | 1.093 | -0.0194 [-0.4] | -0.0126 | -0.0194 | 0.9997 | yes |  |
| n_edge_turnover | 0.921 | 0 | 0.921 [50.1] | 0.9198 | 0.921 | 0.8907 | yes | R |
| n_type_modularity | 0.6443 | -0.0062 | 0.6504 [37.0] | 0.658 | 0.6504 | 0.6115 | yes | R |
| n_degree_cv | 0.4733 | 0.4794 | -0.0062 [-0.3] | -0.02 | -0.0062 | 0.4995 | yes |  |
| n_giant | 0.4025 | 0.985 | -0.5825 [-21.0] | -0.5819 | -0.5825 | 0.4667 | yes | R |

## review:019  (COMBINATION EFFECT)

`N kN=3 | N.COPY(p=1.0) ; N.REWIRE(p=0.1)` · view N · flagged False · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| n_components | 1.386 | 1.386 | 0 [0.0] | -0.4204 | 0 | 1.304 | yes |  |
| n_giant | 0.99 | 0.9875 | 0.0025 [0.1] | 0.0081 | 0.0025 | 0.99 | yes |  |
| n_edge_turnover | 0.9695 | 0 | 0.9695 [50.0] | 0.9938 | 0.9695 | 0.9438 | yes | R |
| n_consensus_gain | 0.635 | 0 | 0.635 [50.0] | 0.6444 | 0.635 | 0.6483 | yes | R |
| n_degree_cv | 0.4887 | 0.4969 | -0.0082 [-0.3] | -0.0246 | -0.0082 | 0.4853 | yes |  |
| n_oscillation | 0.4639 | 0 | 0.4639 [8.6] | 0.6516 | 0.2284 | 0.3771 | yes |  |

## review:020  (DECOMPOSES)

`N kN=4 | N.COPY(p=0.1) ; N.KURA(K=0.1)` · view N · flagged False · verdict KNOWN · reproduced by rules [1]

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| n_components | 1.498 | 1.498 | 0 [0.0] | 0 | 0 | 1.752 | yes |  |
| n_giant | 0.9875 | 0.9875 | 0 [0.0] | 0 | 0 | 0.98 | yes |  |
| n_oscillation | 0.6285 | 0 | 0.6285 [47.6] | 0.8523 | 0.636 | 0.7221 | yes | R |
| n_consensus_gain | 0.4275 | 0 | 0.4275 [17.7] | 0.2562 | 0.7 | 0.5167 | yes | R |
| n_degree_cv | 0.4991 | 0.4991 | 0 [0.0] | 0 | 0 | 0.5121 | yes |  |
| n_type_entropy | 0.6073 | 1.38 | -0.7731 [-19.3] | -0.3587 | -1.317 | 0.3917 | yes |  |

## review:021  (COMBINATION EFFECT)

`C kC=2 | C.COPY(p=1.0) ; C.GAME(g=(1.5, 0.5))` · view C · flagged True · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_n_domains | 4.997 | 6.367 | -1.37 [-6.0] | -1.432 | -1.286 | 5.004 | yes | R |
| g_corrlen | 2 | 1 | 1 [16.7] | 1 | 1 | 2 | yes | R |
| g_n_types | 2 | 2 | 0 [0.0] | 0 | 0 | 2 | yes |  |
| g_reversals | 1 | 0 | 1 [50.0] | 1 | 1 | 1 | yes | R |
| g_type_entropy | 0.9902 | 0.9998 | -0.0096 [-0.2] | -0.009 | -0.0092 | 0.9893 | yes |  |
| g_activity | 0.4623 | 0 | 0.4623 [50.2] | 0.4621 | 0.4619 | 0.4599 | yes | R |

## review:022  (DECOMPOSES)

`C kC=2 | C.MAJ(p=1.0) ; C.SWAP(p=1.0)` · view C · flagged True · verdict KNOWN · reproduced by rules [1]

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_corrlen | 0 | 1 | -1 [-50.0] | 20 | -1 | 14 | no |  |
| g_n_types | 1 | 2 | -1 [-16.7] | 0 | -1 | 1.667 | no |  |
| g_n_domains | 0.6931 | 6.372 | -5.678 [-40.2] | -4.866 | -5.678 | 1.29 | no | R |
| g_reversals | 1 | 0 | 1 [50.0] | 1 | 1 | 1 | yes | R |

## review:023  (COMBINATION EFFECT)

`C kC=2 | C.GAME(g=(1.5, 0.5)) ; C.COPY(p=0.1)` · view C · flagged True · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_n_domains | 3.611 | 6.372 | -2.761 [-13.8] | -2.867 | -2.734 | 3.353 | yes | R |
| g_corrlen | 2 | 1 | 1 [16.7] | 1 | 1 | 2 | yes | R |
| g_n_types | 2 | 2 | 0 [0.0] | 0 | 0 | 2 | yes |  |
| g_reversals | 1 | 0 | 1 [50.0] | 1 | 1 | 1 | yes | R |
| g_largest_domain | 0.9397 | 0.0393 | 0.9004 [45.9] | 0.929 | 0.9053 | 0.9474 | yes | R |
| g_activity_trend | -0.625 | 0 | -0.625 [-50.0] | -0.5871 | -0.425 | -0.6036 | yes | R |

## review:024  (DECOMPOSES)

`C kC=2 | C.GAME(g=(1.5, 0.5)) ; C.SWAP(p=0.1)` · view C · flagged True · verdict KNOWN · reproduced by rules [1]

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_reversals | 1 | 0 | 1 [50.0] | 1 | 1 | 1 | yes | R |

## review:025  (COMBINATION EFFECT)

`C kC=3 | C.MAJ(p=1.0) ; C.SCHELL(th=0.3)` · view C · flagged True · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_n_types | 1 | 3 | -2 [-25.0] | -1 | -2 | 2.333 | no | R |

## review:026  (COMBINATION EFFECT)

`C kC=3 | C.SWAP(p=1.0) ; C.SCHELL(th=0.5)` · view C · flagged True · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_corrlen | 1 | 1 | 0 [0.0] | 9 | 2 | 11.33 | no |  |
| g_n_domains | 5.893 | 7.355 | -1.462 [-5.5] | -1.676 | -1.498 | 5.71 | yes | R |
| g_n_types | 3 | 3 | 0 [0.0] | 0 | 0 | 3 | yes |  |
| g_type_entropy | 0.9999 | 0.9999 | 0 [0.0] | 0 | 0 | 0.9996 | yes |  |
| g_reversals | 0.9417 | 0 | 0.9417 [50.1] | 0.9538 | 0.9591 | 0.9451 | yes | R |
| g_activity_cv | 0.7518 | 0 | 0.7518 [50.1] | 0.5731 | 0.7554 | 0.5554 | no | R |

## review:027  (COMBINATION EFFECT)

`C kC=3 | C.SWAP(p=0.1) ; C.SCHELL(th=0.5)` · view C · flagged True · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_n_domains | 6.086 | 7.347 | -1.262 [-4.7] | -1.315 | -1.272 | 6.079 | yes | R |
| g_n_types | 3 | 3 | 0 [0.0] | 0 | 0 | 3 | yes |  |
| g_corrlen | 1.5 | 1 | 0.5 [1.3] | 1 | 1 | 1.667 | yes |  |
| g_type_entropy | 0.9998 | 0.9998 | 0 [0.0] | 0 | 0 | 0.9997 | yes |  |
| g_reversals | 0.9531 | 0 | 0.9531 [49.9] | 0.9569 | 0.9706 | 0.955 | yes | R |
| g_moran | 0.3568 | 0.0006 | 0.3562 [21.2] | 0.3767 | 0.3786 | 0.37 | yes | R |

## review:028  (COMBINATION EFFECT)

`C kC=3 | C.SCHELL(th=0.5) ; C.COPY(p=1.0)` · view C · flagged True · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_n_domains | 5.861 | 7.337 | -1.476 [-5.6] | -1.635 | -1.514 | 5.926 | yes | R |
| g_n_types | 3 | 3 | 0 [0.0] | 0 | 0 | 3 | yes |  |
| g_corrlen | 3 | 1 | 2 [25.0] | 7 | 0 | 2.333 | yes |  |
| g_type_entropy | 0.9858 | 0.9999 | -0.0141 [-0.4] | -0.0124 | -0.3201 | 0.9305 | yes |  |
| g_reversals | 0.8921 | 0 | 0.8921 [50.1] | 0.9122 | 0.9071 | 0.8826 | yes | R |
| g_spec_peak | 0.3103 | 0 | 0.3103 [35.3] | 0.3828 | 0.3663 | 0.4181 | yes | R |

## review:029  (ARTIFACT)

`N kN=3 | N.MAJ(p=1.0) ; N.REWIRE(p=0.1)` · view N · flagged True · verdict KNOWN-VARIANT

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| n_components | 1.609 | 1.869 | -0.2594 [-0.7] | -0.2154 | -0.2594 | 1.439 | yes |  |
| n_giant | 0.95 | 0.9775 | -0.0275 [-0.7] | -0.4262 | -0.0275 | 0.6933 | no |  |

## review:030  (ARTIFACT)

`C kC=3 | C[any] NEXT_GE(th=1) -> SWAP_RAND() p=1.0` · view C · flagged True · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_n_types | 3 | 3 | 0 [0.0] | 0 | 0 | 3 | yes |  |
| g_type_entropy | 0.9997 | 0.9997 | 0 [0.0] | 0 | 0 | 0.9999 | yes |  |

## review:031  (SINGLE MECHANISM)

`C kC=3 | C[any] SAME_LT(th=0.5) -> MOVE_EMPTY() p=1.0` · view C · flagged False · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_activity_trend | -6.438 | 0 | -6.438 [-50.0] | -6.477 | -4.956 | -6.504 | yes | R |
| g_corrlen | 4 | 1 | 3 [30.0] | 3 | 3 | 4.667 | yes | R |
| g_n_domains | 3.959 | 7.351 | -3.393 [-15.0] | -3.683 | -3.393 | 3.859 | yes | R |
| g_n_types | 3 | 3 | 0 [0.0] | 0 | 0 | 3 | yes |  |
| g_type_entropy | 0.9999 | 0.9999 | 0 [0.0] | 0 | 0 | 0.9999 | yes |  |
| g_moran | 0.7911 | -0.0032 | 0.7943 [50.0] | 0.7989 | 0.7943 | 0.8044 | yes | R |

## review:032  (SINGLE MECHANISM)

`C kC=3 | C[1] SAME_LT(th=0.3) -> MOVE_EMPTY() p=0.1` · view C · flagged False · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_n_domains | 6.686 | 7.376 | -0.6894 [-2.5] | -0.7028 | -0.6894 | 6.669 | yes |  |
| g_activity_trend | -4.855 | 0 | -4.855 [-50.0] | -4.882 | -3.868 | -4.907 | yes | R |
| g_n_types | 3 | 3 | 0 [0.0] | 0 | 0 | 3 | yes |  |
| g_corrlen | 1 | 1 | 0 [0.0] | 0 | 0 | 1 | yes |  |
| g_reversals | 1 | 0 | 1 [50.0] | 1 | 1 | 1 | yes | R |
| g_type_entropy | 0.9999 | 0.9999 | 0 [0.0] | 0 | 0 | 0.9999 | yes |  |

## review:033  (SINGLE MECHANISM)

`C kC=4 | C[any] NEXT_GE(th=1) -> SET_NEXT() p=0.1` · view C · flagged False · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_n_domains | 7.383 | 7.67 | -0.2875 [-1.0] | -0.2831 | -0.2871 | 7.417 | yes |  |
| g_n_types | 4 | 4 | 0 [0.0] | 0 | 0 | 4 | yes |  |
| g_corrlen | 1 | 1 | 0 [0.0] | 0 | 0 | 1 | yes |  |
| g_type_entropy | 1 | 0.9998 | 0.0002 [0.0] | 0 | 0.0002 | 1 | yes |  |
| g_interface | 0.66 | 0.7498 | -0.0898 [-3.2] | -0.0904 | -0.0915 | 0.6631 | yes | R |
| g_cyclic_changes | 0.3556 | -0.3333 | 0.6889 [49.9] | 0.6892 | 0.3382 | 0.3539 | yes |  |

## review:034  (SINGLE MECHANISM)

`C kC=4 | C[any] NEXT_GE(th=1) -> COPY_RAND() p=1.0` · view C · flagged True · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_n_domains | 5.418 | 7.678 | -2.26 [-8.6] | -2.323 | -2.33 | 5.433 | yes | R |
| g_corrlen | 3 | 1 | 2 [25.0] | 2 | 1 | 4.333 | no | R |
| g_n_types | 4 | 4 | 0 [0.0] | 0 | 0 | 4 | yes |  |
| g_activity_trend | -0.7415 | 0 | -0.7415 [-47.8] | -0.9015 | -0.8283 | -1.031 | yes | R |
| g_type_entropy | 0.9636 | 0.9998 | -0.0362 [-0.9] | -0.0105 | -0.21 | 0.9287 | yes |  |
| g_moran | 0.5784 | -0.0017 | 0.58 [36.2] | 0.6014 | 0.5275 | 0.5946 | yes | R |

## review:035  (ARTIFACT)

`N kN=3 | N[any] NEXT_GE(th=1) -> COPY_RAND() p=0.1` · view N · flagged True · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| n_components | 1.609 | 1.609 | 0 [0.0] | 0 | 0 | 1.731 | yes |  |
| n_giant | 0.985 | 0.985 | 0 [0.0] | 0 | 0 | 0.9817 | yes |  |
| n_degree_cv | 0.491 | 0.491 | 0 [0.0] | 0 | 0 | 0.492 | yes |  |

## review:036  (SINGLE MECHANISM)

`N kN=4 | N[any] ALWAYS() -> COPY_RAND() p=0.2` · view N · flagged True · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| n_components | 1.701 | 1.701 | 0 [0.0] | 0 | 0 | 1.521 | yes |  |
| n_giant | 0.98 | 0.98 | 0 [0.0] | 0 | 0 | 0.985 | yes |  |
| n_oscillation | 0.7199 | 0 | 0.7199 [50.0] | 0.6997 | 0.3663 | 0.8078 | yes | R |
| n_consensus_gain | 0.7 | 0 | 0.7 [50.0] | 0.2687 | 0.71 | 0.6867 | yes |  |

## review:037  (COMBINATION EFFECT)

`C kC=2 | C.MAJ(p=0.1) ; C.GAME(g=(1.5, -0.5))` · view C · flagged True · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_n_domains | 3.676 | 6.378 | -2.701 [-13.4] | -2.967 | -2.728 | 3.631 | yes | R |
| g_corrlen | 2 | 1 | 1 [16.7] | 1 | 1 | 2.333 | yes | R |
| g_n_types | 2 | 2 | 0 [0.0] | 0 | 0 | 2 | yes |  |
| g_activity_trend | 2.124 | 0 | 2.124 [35.0] | 1.84 | 0.5468 | 1.486 | yes |  |
| g_reversals | 1 | 0 | 1 [50.0] | 1 | 1 | 1 | yes | R |
| g_type_entropy | 0.8378 | 0.9998 | -0.162 [-4.4] | -0.1681 | -0.1696 | 0.8321 | yes | R |

## review:038  (DECOMPOSES)

`C kC=4 | C.COPY(p=1.0) ; C.SCHELL(th=0.7)` · view C · flagged True · verdict KNOWN · reproduced by rules [1]

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_n_domains | 7.085 | 7.673 | -0.5881 [-2.0] | -0.2522 | -1.988 | 6.746 | yes |  |
| g_n_types | 4 | 4 | 0 [0.0] | 0 | -1 | 3.667 | yes |  |
| g_corrlen | 1 | 1 | 0 [0.0] | 0 | 4 | 2.667 | no |  |
| g_type_entropy | 0.9054 | 0.9998 | -0.0944 [-2.5] | -0.0281 | -0.2873 | 0.8703 | yes |  |
| g_spec_peak | 0.4684 | 0 | 0.4684 [18.3] | 0.682 | 0.6108 | 0.59 | yes | R |
| g_interface | 0.5504 | 0.7497 | -0.1993 [-7.7] | -0.0969 | -0.4398 | 0.4955 | yes |  |

## review:039  (DECOMPOSES)

`N kN=4 | N.CYCLE(th=3,p=1.0) ; N.COPY(p=0.1)` · view N · flagged False · verdict KNOWN · reproduced by rules [2]

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| n_components | 1.609 | 1.609 | 0 [0.0] | 0 | 0 | 1.866 | yes |  |
| n_giant | 0.985 | 0.985 | 0 [0.0] | 0 | 0 | 0.9767 | yes |  |
| n_oscillation | 0.5661 | 0 | 0.5661 [19.9] | 0.7324 | 0.4617 | 0.6412 | yes | R |
| n_consensus_gain | 0.65 | 0 | 0.65 [14.6] | 0.4587 | 0.71 | 0.5383 | yes | R |
| n_degree_cv | 0.4816 | 0.4816 | 0 [0.0] | 0 | 0 | 0.484 | yes |  |
| n_type_entropy | 0.1982 | 1.379 | -1.181 [-10.5] | -0.7516 | -1.319 | 0.4339 | yes | R |

## review:040  (ARTIFACT)

`N kN=2 | N.REWIRE(p=0.01) ; N.GAME(g=(1.9, 0.0))` · view N · flagged True · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| n_giant | 0.99 | 0.9775 | 0.0125 [0.3] | 0.0131 | 0.0175 | 0.9967 | yes |  |
| n_components | 1.386 | 1.869 | -0.4825 [-1.5] | -0.9159 | -0.7702 | 0.9635 | yes |  |

## review:041  (SINGLE MECHANISM)

`C kC=2 | C[any] CNT_GE(b=0,th=3) -> IMITATE_BEST(g=(1.5, 0.5)) p=1.0` · view C · flagged True · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_n_domains | 4.401 | 6.357 | -1.956 [-9.1] | -1.972 | -1.963 | 4.289 | yes | R |
| g_corrlen | 2 | 1 | 1 [16.7] | 1 | 1 | 2 | yes | R |
| g_n_types | 2 | 2 | 0 [0.0] | 0 | 0 | 2 | yes |  |
| g_reversals | 1 | 0 | 1 [50.0] | 1 | 1 | 1 | yes | R |
| g_largest_domain | 0.8939 | 0.0475 | 0.8464 [45.0] | 0.872 | 0.8457 | 0.897 | yes | R |
| g_moran | 0.4901 | 0.0062 | 0.4839 [25.1] | 0.5004 | 0.4872 | 0.5051 | yes | R |

## review:042  (SINGLE MECHANISM)

`C kC=2 | C[any] CNT_LE(b=0,th=3) -> IMITATE_BEST(g=(1.9, 0.0)) p=0.1` · view C · flagged True · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_n_domains | 4.382 | 6.346 | -1.964 [-9.2] | -2.202 | -1.951 | 4.334 | yes | R |
| g_corrlen | 2 | 1 | 1 [16.7] | 1 | 1 | 2 | yes | R |
| g_n_types | 2 | 2 | 0 [0.0] | 0 | 0 | 2 | yes |  |
| g_activity_trend | -1.084 | 0 | -1.084 [-49.9] | -1.128 | -1.077 | -1.093 | yes | R |
| g_reversals | 1 | 0 | 1 [50.0] | 1 | 1 | 1 | yes | R |
| g_type_entropy | 0.9994 | 0.9999 | -0.0005 [-0.0] | -0 | -0.0031 | 0.9997 | yes |  |

## review:043  (SINGLE MECHANISM)

`A kA=1 v=0.3 eta=0.1 | A[any] ALWAYS() -> HEAD_MEAN(b=any,r=0.5) p=1.0` · view A · flagged False · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| a_density_cv | 5.357 | 0.8195 | 4.538 [20.5] | 4.624 | 2.753 | 6.233 | yes | R |
| a_step | 5.997 | 5.992 | 0.0053 [0.0] | 0.0056 | 3.086 | 5.995 | yes |  |
| a_local_align | 0.8877 | -0.042 | 0.9296 [50.0] | 0.9526 | 0.883 | 0.8894 | yes | R |
| a_polar | 0.9889 | 0.0272 | 0.9617 [47.4] | 0.8028 | 0.9652 | 0.8439 | yes | R |
| a_ang_mom | 0.2113 | 0.0299 | 0.1814 [1.8] | 0.1828 | 0.0755 | 0.3267 | yes |  |
| a_nn_ratio | 0.2138 | 0.9958 | -0.7819 [-16.5] | -0.8273 | -0.6539 | 0.1769 | yes | R |

## review:044  (DECOMPOSES)

`C kC=2 | C.SWAP(p=1.0) ; C.CYCLE(th=3,p=0.1)` · view C · flagged True · verdict KNOWN · reproduced by rules [1]

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_corrlen | 0 | 1 | -1 [-50.0] | 23.5 | -1 | 20.5 | no |  |
| g_n_domains | 0.6931 | 6.378 | -5.685 [-40.2] | -1.688 | -5.685 | 4.605 | no |  |
| g_n_types | 1.5 | 2 | -0.5 [-1.3] | 0 | -1 | 2 | yes |  |
| g_activity_trend | -3.394 | 0 | -3.394 [-3.6] | -1.561 | -7.573 | -1.135 | yes |  |
| g_reversals | 1 | 0 | 1 [50.0] | 1 | 1 | 1 | yes | R |
| g_type_entropy | 0.2095 | 1 | -0.7905 [-6.2] | -0.1051 | -1 | 0.8362 | yes |  |

## review:045  (COMBINATION EFFECT)

`C kC=3 | C.COPY(p=1.0) ; C.CYCLE(th=3,p=1.0)` · view C · flagged True · verdict KNOWN-VARIANT

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_corrlen | 5 | 1 | 4 [33.3] | 4 | 4 | 4.667 | yes | R |
| g_n_domains | 3.526 | 7.352 | -3.826 [-17.6] | -4.125 | -3.841 | 3.434 | yes | R |
| g_n_types | 3 | 3 | 0 [0.0] | 0 | 0 | 3 | yes |  |
| g_type_entropy | 0.9998 | 0.9998 | 0.0001 [0.0] | 0 | 0.0001 | 0.9986 | yes |  |
| g_moran | 0.8219 | -0.0025 | 0.8244 [50.0] | 0.8243 | 0.8143 | 0.8201 | yes | R |
| g_reversals | 0.5017 | 0 | 0.5017 [50.2] | 0.5008 | 0.4277 | 0.5044 | yes | R |

## review:046  (COMBINATION EFFECT)

`C kC=3 | C.MAJ(p=1.0) ; C.CYCLE(th=1,p=1.0)` · view C · flagged True · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_n_domains | 3.97 | 7.356 | -3.385 [-14.9] | -3.735 | -3.464 | 3.975 | yes | R |
| g_corrlen | 3 | 1 | 2 [25.0] | 2 | 2 | 3 | yes | R |
| g_n_types | 3 | 3 | 0 [0.0] | 0 | 0 | 3 | yes |  |
| g_type_entropy | 0.9998 | 0.9998 | 0 [0.0] | -0 | 0.0001 | 0.9997 | yes |  |
| g_reversals | 0.9407 | 0 | 0.9407 [50.0] | 0.9382 | 0.8389 | 0.9404 | yes | R |
| g_moran | 0.7195 | -0.0023 | 0.7218 [44.0] | 0.7208 | 0.7241 | 0.7214 | yes | R |

## review:047  (COMBINATION EFFECT)

`C kC=3 | C.CYCLE(th=1,p=1.0) ; C.SCHELL(th=0.5)` · view C · flagged True · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_n_domains | 5.37 | 7.366 | -1.996 [-7.8] | -2.337 | -2.208 | 5.141 | yes | R |
| g_n_types | 3 | 3 | 0 [0.0] | 0 | 0 | 3 | yes |  |
| g_corrlen | 1 | 1 | 0 [0.0] | 0 | 0 | 1 | yes |  |
| g_largest_domain | 0.814 | 0.0076 | 0.8064 [49.2] | 0.8742 | 0.8568 | 0.8578 | yes | R |
| g_activity_trend | -0.6496 | 0 | -0.6496 [-12.0] | -0.8596 | -0.2554 | -0.8514 | yes |  |
| g_reversals | 0.5601 | 0 | 0.5601 [50.0] | 0.5577 | 0.6379 | 0.5599 | yes | R |

## review:048  (DECOMPOSES)

`C kC=3 | C.CYCLE(th=3,p=0.1) ; C.MAJ(p=1.0)` · view C · flagged True · verdict KNOWN · reproduced by rules [2]

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_corrlen | 0 | 1 | -1 [-50.0] | 24 | -1 | 3.667 | no |  |
| g_activity_trend | -5.442 | 0 | -5.442 [-50.0] | -1.428 | -4.902 | -2.861 | no |  |
| g_n_types | 1 | 3 | -2 [-25.0] | 0 | -2 | 2 | no |  |
| g_largest_domain | 1 | 0.0084 | 0.9916 [49.1] | 0.7286 | 0.9916 | 0.9644 | yes | R |
| g_reversals | 0.8076 | 0 | 0.8076 [49.9] | 0.8611 | 0.5718 | 0.8526 | yes | R |
| g_n_domains | 0.6931 | 7.346 | -6.653 [-41.4] | -6.704 | -6.653 | 0.8283 | yes | R |

## review:049  (COMBINATION EFFECT)

`C kC=3 | C.CYCLE(th=3,p=0.1) ; C.SCHELL(th=0.3)` · view C · flagged True · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_n_domains | 4.735 | 7.357 | -2.623 [-10.8] | -2.759 | -2.713 | 4.964 | yes | R |
| g_n_types | 3 | 3 | 0 [0.0] | 0 | 0 | 3 | yes |  |
| g_corrlen | 3 | 1 | 2 [25.0] | 2 | 2 | 2.667 | yes | R |
| g_type_entropy | 0.9132 | 0.9998 | -0.0866 [-2.3] | -0.0926 | -0.092 | 0.8962 | yes |  |
| g_moran | 0.6221 | -0.0056 | 0.6277 [37.1] | 0.6371 | 0.6277 | 0.5812 | yes | R |
| g_activity_trend | -0.407 | 0 | -0.407 [-25.4] | -0.3706 | -0.0898 | -0.4451 | yes |  |

## review:050  (COMBINATION EFFECT)

`C kC=4 | C.CYCLE(th=1,p=1.0) ; C.SCHELL(th=0.3)` · view C · flagged True · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_n_domains | 6.447 | 7.665 | -1.218 [-4.3] | -1.303 | -1.359 | 6.35 | yes | R |
| g_n_types | 4 | 4 | 0 [0.0] | 0 | 0 | 4 | yes |  |
| g_corrlen | 1 | 1 | 0 [0.0] | 0 | 0 | 2 | no |  |
| g_type_entropy | 0.9529 | 0.9997 | -0.0468 [-1.2] | -0.0476 | -0.044 | 0.9606 | yes |  |
| g_activity | 0.6131 | 0 | 0.6131 [49.8] | 0.6208 | 0.6229 | 0.6102 | yes | R |
| g_moran | 0.3197 | 0.0098 | 0.3099 [14.5] | 0.3401 | 0.3122 | 0.4407 | yes | R |

## review:051  (COMBINATION EFFECT)

`C kC=4 | C.SCHELL(th=0.5) ; C.CYCLE(th=1,p=1.0)` · view C · flagged True · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_n_domains | 6.244 | 7.669 | -1.425 [-5.1] | -1.404 | -1.316 | 6.326 | yes | R |
| g_n_types | 4 | 4 | 0 [0.0] | 0 | 0 | 4 | yes |  |
| g_corrlen | 1 | 1 | 0 [0.0] | 0.5 | 0 | 1.333 | yes |  |
| g_type_entropy | 0.9874 | 0.9997 | -0.0123 [-0.3] | -0.0141 | -0.018 | 0.9828 | yes |  |
| g_activity | 0.6728 | 0 | 0.6728 [49.8] | 0.6756 | 0.6735 | 0.6729 | yes | R |
| g_interface | 0.3755 | 0.7498 | -0.3743 [-16.6] | -0.3757 | -0.3808 | 0.3752 | yes | R |

## review:052  (COMBINATION EFFECT)

`N kN=4 | N.MAJ(p=1.0) ; N.REWIRE(p=0.03)` · view N · flagged False · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| n_components | 1.701 | 1.701 | 0 [0.0] | -0.0303 | 0 | 1.792 | yes |  |
| n_giant | 0.9725 | 0.9825 | -0.01 [-0.3] | -0.0194 | -0.01 | 0.9 | yes |  |
| n_consensus_gain | 0.6725 | 0 | 0.6725 [50.2] | 0.6931 | 0.6725 | 0.5983 | yes | R |
| n_degree_cv | 0.4883 | 0.4941 | -0.0058 [-0.1] | -0.0061 | -0.0058 | 0.4948 | yes |  |
| n_type_entropy | 0.1175 | 1.378 | -1.26 [-42.1] | -1.236 | -1.26 | 0.2705 | yes | R |
| n_edge_turnover | 0.1948 | 0 | 0.1948 [6.8] | 0.2416 | 0.1948 | 0.1551 | yes | R |

## review:053  (COMBINATION EFFECT)

`N kN=4 | N.CYCLE(th=1,p=1.0) ; N.REWIRE(p=1.0)` · view N · flagged True · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| n_components | 1.869 | 1.609 | 0.2594 [0.8] | -0.6675 | 0.2594 | 1.962 | yes |  |
| n_edge_turnover | 0.982 | 0 | 0.982 [50.1] | 0.9941 | 0.982 | 0.9812 | yes | R |
| n_giant | 0.8625 | 0.985 | -0.1225 [-3.3] | 0.0112 | -0.1225 | 0.9067 | yes |  |
| n_consensus_gain | 0.5925 | 0 | 0.5925 [49.8] | 0.0269 | 0.5925 | 0.6267 | yes |  |
| n_oscillation | 0.5207 | 0 | 0.5207 [3.4] | 0.3567 | 0.2706 | 0.463 | yes | R |
| n_degree_cv | 0.4626 | 0.4925 | -0.0299 [-1.6] | -0.0419 | -0.0299 | 0.455 | yes |  |

## review:054  (COMBINATION EFFECT)

`C kC=2 | C.MAJ(p=0.003) ; C.GAME(g=(1.3, 0.0))` · view C · flagged True · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_activity_cv | 2.477 | 0 | 2.477 [21.7] | 1.231 | 3.667 | 3.471 | yes |  |
| g_n_domains | 3.384 | 6.372 | -2.987 [-15.3] | -3.038 | -4.762 | 3.257 | yes | R |
| g_activity_trend | -2.717 | 0 | -2.717 [-50.0] | -2.858 | -2.399 | -3.203 | yes | R |
| g_corrlen | 2 | 1 | 1 [16.7] | 1 | 1 | 2 | yes | R |
| g_n_types | 2 | 2 | 0 [0.0] | 0 | -1 | 2 | yes |  |
| g_reversals | 1 | 0 | 1 [50.0] | 1 | 1 | 1 | yes | R |

## review:055  (DECOMPOSES)

`C kC=2 | C.CYCLE(th=7,p=0.1) ; C.GAME(g=(1.5, 0.5))` · view C · flagged True · verdict KNOWN · reproduced by rules [2]

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_n_domains | 4.29 | 6.356 | -2.066 [-9.7] | -2.018 | -2.032 | 4.317 | yes | R |
| g_corrlen | 2 | 1 | 1 [16.7] | 1 | 1 | 2 | yes | R |
| g_n_types | 2 | 2 | 0 [0.0] | 0 | 0 | 2 | yes |  |
| g_reversals | 1 | 0 | 1 [50.0] | 1 | 1 | 1 | yes | R |
| g_largest_domain | 0.9082 | 0.0542 | 0.854 [44.5] | 0.8814 | 0.8541 | 0.9045 | yes | R |
| g_activity_trend | 0.596 | 0 | 0.596 [50.1] | 0.7494 | 0.1573 | 0.5306 | yes |  |

## review:056  (COMBINATION EFFECT)

`C kC=2 | C.GAME(g=(1.5, 0.5)) ; C.SWAP(p=0.05)` · view C · flagged True · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_n_domains | 4.477 | 6.378 | -1.9 [-8.8] | -1.892 | -1.895 | 4.582 | yes | R |
| g_corrlen | 2 | 1 | 1 [16.7] | 1 | 1 | 2 | yes | R |
| g_n_types | 2 | 2 | 0 [0.0] | 0 | 0 | 2 | yes |  |
| g_reversals | 1 | 0 | 1 [50.0] | 1 | 1 | 1 | yes | R |
| g_largest_domain | 0.8746 | 0.0454 | 0.8292 [45.1] | 0.8449 | 0.813 | 0.8311 | yes | R |
| g_type_entropy | 0.5675 | 0.9999 | -0.4324 [-13.8] | -0.4297 | -0.4257 | 0.5776 | yes | R |

## review:057  (COMBINATION EFFECT)

`C kC=3 | C.FLIP(a=1,b=0,p=0.1) ; C.SCHELL(th=0.5)` · view C · flagged True · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_n_domains | 5.468 | 6.361 | -0.8932 [-3.8] | -0.9016 | -0.8932 | 5.503 | yes | R |
| g_activity_trend | -3.781 | -3.359 | -0.4227 [-3.0] | -0.498 | -0.2987 | -3.818 | yes |  |
| g_n_types | 2 | 2 | 0 [0.0] | 0 | 0 | 2 | yes |  |
| g_corrlen | 1 | 1 | 0 [0.0] | 0 | 0 | 1 | yes |  |
| g_reversals | 0.7064 | 0 | 0.7064 [50.1] | 0.7069 | 0.6711 | 0.7236 | yes | R |
| g_largest_domain | 0.5842 | 0.6194 | -0.0352 [-1.5] | -0.1044 | -0.0352 | 0.5813 | yes |  |

## review:058  (COMBINATION EFFECT)

`C kC=3 | C.FLIP(a=1,b=0,p=0.1) ; C.SCHELL(th=0.3)` · view C · flagged True · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_n_domains | 5.939 | 6.356 | -0.4173 [-1.7] | -0.4481 | -0.4173 | 5.837 | yes |  |
| g_activity_trend | -4.003 | -3.37 | -0.6337 [-4.3] | -0.5777 | -0.5706 | -3.932 | yes | R |
| g_n_types | 2 | 2 | 0 [0.0] | 0 | 0 | 2 | yes |  |
| g_corrlen | 1 | 1 | 0 [0.0] | 0 | 0 | 1 | yes |  |
| g_reversals | 0.8364 | 0 | 0.8364 [50.1] | 0.8356 | 0.6714 | 0.8073 | yes | R |
| g_type_entropy | 0.5803 | 0.5803 | 0 [0.0] | 0 | 0 | 0.5744 | yes |  |

## review:059  (ARTIFACT)

`N kN=2 | N.MAJ(p=1.0) ; N.REWIRE(p=0.5)` · view N · flagged True · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| n_components | 1.701 | 1.609 | 0.0912 [0.3] | -0.3483 | 0.0912 | 1.386 | yes |  |

## review:060  (SINGLE MECHANISM)

`CA kC=4 kA=3 v=0.3 eta=0.3 | C.CYCLE(th=3,p=1.0)` · view C · flagged True · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_corrlen | 10 | 1 | 9 [40.9] | 10 | 7 | 6.333 | yes | R |
| g_n_types | 4 | 4 | 0 [0.0] | 0 | 0 | 3.333 | yes |  |
| g_n_domains | 3.577 | 7.669 | -4.092 [-18.2] | -4.361 | -4.208 | 3.274 | yes | R |
| g_activity_trend | -0.5756 | 0 | -0.5756 [-18.4] | -0.6441 | -0.5527 | -2.441 | no | R |
| g_moran | 0.8744 | -0.0026 | 0.877 [50.1] | 0.8691 | 0.8555 | 0.8327 | yes | R |
| g_type_entropy | 0.9976 | 0.9997 | -0.002 [-0.1] | -0.0023 | -0.0133 | 0.7958 | no |  |

## review:061  (DECOMPOSES)

`C kC=2 | C.FLIP(a=1,b=0,p=0.03) ; C.GAME(g=(1.9, 0.0))` · view C · flagged True · verdict KNOWN · reproduced by rules [2]

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_n_domains | 4.204 | 0.6931 | 3.511 [35.9] | 4.781 | 3.376 | 4.144 | yes | R |
| g_corrlen | 2.5 | 0 | 2.5 [6.7] | 3 | 3 | 2.667 | yes | R |
| g_n_types | 2 | 1 | 1 [16.7] | 1 | 1 | 2 | yes | R |
| g_activity_trend | 0.2495 | -4.078 | 4.327 [50.0] | 4.448 | 2.755 | 1.102 | no | R |
| g_reversals | 1 | 0 | 1 [50.0] | 1 | 1 | 1 | yes | R |
| g_largest_domain | 0.7924 | 1 | -0.2076 [-5.8] | -0.2029 | -0.2129 | 0.7992 | yes | R |

## review:062  (COMBINATION EFFECT)

`C kC=3 | C.CYCLE(th=3,p=0.3) ; C.SWAP(p=1.0)` · view C · flagged True · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_corrlen | 6 | 1 | 5 [35.7] | 5 | 4 | 6 | yes | R |
| g_n_domains | 5.405 | 7.345 | -1.94 [-7.6] | -1.715 | -1.557 | 5.483 | yes | R |
| g_n_types | 3 | 3 | 0 [0.0] | 0 | 0 | 3 | yes |  |
| g_type_entropy | 0.9987 | 0.9998 | -0.0011 [-0.0] | -0 | 0 | 0.9998 | yes |  |
| g_moran | 0.6708 | 0.0083 | 0.6625 [40.9] | 0.6181 | 0.5373 | 0.638 | yes | R |
| g_cyclic_changes | 0.3924 | -0.5 | 0.8924 [50.1] | 0.8981 | 0.519 | 0.409 | yes | R |

## review:063  (COMBINATION EFFECT)

`C kC=4 | C.MAJ(p=0.005) ; C.SCHELL(th=0.7)` · view C · flagged True · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_n_domains | 6.353 | 7.666 | -1.314 [-4.7] | -1.107 | -5.181 | 5.941 | yes | R |
| g_n_types | 4 | 4 | 0 [0.0] | 0 | -1.5 | 4 | yes |  |
| g_corrlen | 2.5 | 1 | 1.5 [4.0] | 0 | 8 | 2.333 | yes |  |
| g_activity_trend | -0.6684 | 0 | -0.6684 [-49.9] | -0.3032 | -7.387 | -1.437 | yes |  |
| g_type_entropy | 0.8484 | 0.9998 | -0.1513 [-4.1] | -0.0853 | -0.5177 | 0.8184 | yes | R |
| g_spec_peak | 0.6579 | 0 | 0.6579 [49.8] | 0.6041 | 0.6592 | 0.6632 | yes | R |

## review:064  (COMBINATION EFFECT)

`C kC=4 | C.CYCLE(th=2,p=1.0) ; C.SWAP(p=1.0)` · view C · flagged True · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_n_domains | 4.219 | 7.672 | -3.453 [-14.5] | -2.993 | -3.903 | 4.931 | yes | R |
| g_corrlen | 4.5 | 1 | 3.5 [9.4] | 3.5 | 4.5 | 4.333 | yes | R |
| g_n_types | 4 | 4 | 0 [0.0] | 0 | 0 | 4 | yes |  |
| g_type_entropy | 0.9977 | 0.9997 | -0.002 [-0.1] | -0.0005 | -0 | 0.9989 | yes |  |
| g_moran | 0.7982 | 0.002 | 0.7962 [49.8] | 0.7116 | 0.8267 | 0.6906 | yes | R |
| g_activity | 0.4775 | 0 | 0.4775 [50.3] | 0.4759 | 0.4701 | 0.5237 | yes | R |

## review:065  (ARTIFACT)

`C kC=4 | C.SCHELL(th=0.5) ; C.SAND(p=0.001)` · view C · flagged True · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_n_domains | 7.062 | 7.664 | -0.6014 [-2.0] | -0.4157 | -0.5155 | 7.012 | yes |  |

## review:066  (COMBINATION EFFECT)

`C kC=2 | C.CYCLE(th=7,p=0.03) ; C.GAME(g=(1.5, 0.5))` · view C · flagged True · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_n_domains | 4.086 | 6.346 | -2.26 [-10.8] | -2.075 | -2.119 | 4.209 | yes | R |
| g_corrlen | 2 | 1 | 1 [16.7] | 1 | 1 | 2 | yes | R |
| g_n_types | 2 | 2 | 0 [0.0] | 0 | 0 | 2 | yes |  |
| g_activity_trend | 1.522 | 0 | 1.522 [24.4] | 1.68 | 0.8214 | 1.547 | yes | R |
| g_reversals | 1 | 0 | 1 [50.0] | 1 | 1 | 1 | yes | R |
| g_largest_domain | 0.9237 | 0.0446 | 0.8792 [45.3] | 0.8918 | 0.8721 | 0.9163 | yes | R |

## review:067  (COMBINATION EFFECT)

`C kC=2 | C.FLIP(a=1,b=0,p=0.01) ; C.GAME(g=(1.3, 0.0))` · view C · flagged True · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_n_domains | 3.178 | 0.6931 | 2.485 [32.1] | 3.749 | 2.274 | 3.147 | yes | R |
| g_corrlen | 2 | 0 | 2 [50.0] | 2.5 | 2 | 2.667 | yes | R |
| g_n_types | 2 | 1 | 1 [16.7] | 1 | 1 | 2 | yes | R |
| g_reversals | 1 | 0 | 1 [50.0] | 1 | 1 | 1 | yes | R |
| g_largest_domain | 0.9308 | 1 | -0.0692 [-1.8] | -0.0727 | -0.0867 | 0.9351 | yes |  |
| g_activity_trend | -0.798 | -3.925 | 3.127 [26.1] | 3.346 | 3.189 | -0.7485 | yes | R |

## review:068  (COMBINATION EFFECT)

`C kC=3 | C.CONTAG(a=0,b=1,th=3,p=1.0) ; C.SCHELL(th=0.3)` · view C · flagged True · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_n_types | 2 | 3 | -1 [-10.0] | -1 | -1 | 2 | yes | R |

## review:069  (COMBINATION EFFECT)

`C kC=3 | C.CYCLE(th=1,p=0.1) ; C.FLIP(a=0,b=1,p=0.1)` · view C · flagged True · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_n_domains | 6.917 | 6.338 | 0.579 [2.2] | 0.5692 | 0.5854 | 6.942 | yes |  |
| g_n_types | 3 | 2 | 1 [10.0] | 1 | 1 | 3 | yes | R |
| g_corrlen | 1 | 1 | 0 [0.0] | 0 | 0 | 1 | yes |  |
| g_type_entropy | 0.8868 | 0.5824 | 0.3045 [10.4] | 0.3063 | 0.3057 | 0.8871 | yes | R |
| g_interface | 0.5269 | 0.4447 | 0.0822 [4.2] | 0.0817 | 0.0837 | 0.5275 | yes | R |
| g_reversals | 0.4959 | 0 | 0.4959 [50.1] | 0.4951 | 0.5953 | 0.495 | yes | R |

## review:070  (DECOMPOSES)

`C kC=4 | C.CYCLE(th=3,p=1.0) ; C.KURA(K=0.3)` · view C · flagged True · verdict KNOWN · reproduced by rules [1]

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_corrlen | 10 | 1 | 9 [40.9] | 6 | 7.5 | 14 | yes | R |
| g_n_types | 4 | 4 | 0 [0.0] | 0 | 0 | 3.667 | yes |  |
| g_n_domains | 3.44 | 7.688 | -4.248 [-19.1] | -4.348 | -4.598 | 3.436 | yes | R |
| g_activity_trend | -0.6729 | 0 | -0.6729 [-34.5] | -0.535 | -2.084 | -1.431 | yes | R |
| g_type_entropy | 0.9847 | 0.9998 | -0.0152 [-0.4] | -0.0015 | -0.4803 | 0.7617 | no |  |
| g_moran | 0.8822 | -0.0048 | 0.887 [46.4] | 0.8526 | 0.8499 | 0.733 | yes | R |

## review:071  (ARTIFACT)

`N kN=3 | N.REWIRE(p=0.2) ; N.CYCLE(th=1,p=0.1)` · view N · flagged True · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| n_components | 1.099 | 1.701 | -0.602 [-1.2] | -1.131 | -0.3143 | 1.269 | yes |  |

## review:072  (ARTIFACT)

`N kN=2 | N.COPY(p=1.0) ; N.GAME(g=(1.9, 0.0)) ; N.REWIRE(p=1.0)` · view N · flagged True · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| n_components | 1.099 | 1.792 | -0.6931 [-2.0] | -0.8856 | -0.6931 | 1.099 | yes |  |
| n_giant | 0.985 | 0.9775 | 0.0075 [0.2] | 0.0069 | 0.0075 | 0.9567 | yes |  |

## review:073  (DECOMPOSES)

`C kC=2 | C.SWAP(p=0.3) ; C.CYCLE(th=4,p=0.01) ; C.GAME(g=(1.3, 0.0)) ; C.KURA(K=1.0)` · view C · flagged True · verdict KNOWN · reproduced by rules [3]

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_n_domains | 4.27 | 6.361 | -2.092 [-9.8] | -2.383 | -2.284 | 3.712 | yes | R |
| g_corrlen | 2 | 1 | 1 [16.7] | 1 | 1 | 3 | no | R |
| g_n_types | 2 | 2 | 0 [0.0] | 0 | 0 | 2 | yes |  |
| g_reversals | 1 | 0 | 1 [50.0] | 1 | 1 | 1 | yes | R |
| g_largest_domain | 0.7949 | 0.0371 | 0.7578 [45.4] | 0.8168 | 0.8298 | 0.9029 | yes | R |
| g_moran | 0.5749 | 0.0077 | 0.5672 [48.5] | 0.609 | 0.572 | 0.6077 | yes | R |

## qd:001  (SINGLE MECHANISM)

`AF kA=1 v=0.3 eta=1.0 fields=1 D=[0.05] | A.CLIMB(f=0)` · view A · flagged False · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| a_density_cv | 11.09 | 0.8207 | 10.27 [43.1] | 7.034 | 15.17 | 11.63 | yes | R |
| a_step | 0.2559 | 5.226 | -4.97 [-45.3] | -4.739 | -12.03 | 0.7122 | no | R |
| a_ang_mom | 0.0186 | 0.0187 | -0.0001 [-0.0] | 0.0236 | -0.0053 | 0.1001 | no |  |
| a_polar | 0.0574 | 0.0419 | 0.0155 [1.6] | 0.0365 | 0.0062 | 0.08836 | yes |  |
| a_polar_std | 0.0333 | 0.0206 | 0.0127 [2.3] | 0.0311 | 0.0005 | 0.08358 | no |  |
| a_local_align | -0.0011 | -0.0328 | 0.0317 [2.3] | 0.0627 | 0.0183 | 0.06496 | no |  |

## qd:002  (SINGLE MECHANISM)

`CA kC=4 kA=1 v=0.0 eta=0.3 | A.WRITECELL()` · view C · flagged False · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_n_domains | 7.676 | 7.676 | 0 [0.0] | 0 | 0 | 7.666 | yes |  |
| g_n_types | 4 | 4 | 0 [0.0] | 0 | 0 | 4 | yes |  |
| g_corrlen | 1 | 1 | 0 [0.0] | 0 | 0 | 1 | yes |  |
| g_reversals | 1 | 0 | 1 [50.0] | 1 | 1 | 1 | yes | R |
| g_type_entropy | 0.9997 | 0.9996 | 0.0001 [0.0] | -0 | 0.0001 | 0.9998 | yes |  |
| g_interface | 0.7494 | 0.752 | -0.0025 [-0.1] | -0.0007 | -0.0025 | 0.7476 | yes |  |

## qd:003  (DECOMPOSES)

`C kC=2 | C.GAME(g=(1.5, 0.5)) ; C.GAME(g=(1.5, 0.5))` · view C · flagged False · verdict KNOWN · reproduced by rules [1]

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_n_domains | 4.771 | 6.365 | -1.594 [-7.2] | -1.587 | -1.594 | 4.795 | yes | R |
| g_corrlen | 2 | 1 | 1 [16.7] | 1 | 1 | 2 | yes | R |
| g_n_types | 2 | 2 | 0 [0.0] | 0 | 0 | 2 | yes |  |
| g_reversals | 1 | 0 | 1 [50.0] | 1 | 0 | 1 | yes |  |
| g_activity_trend | -0.3679 | 0 | -0.3679 [-49.7] | -0.2232 | 0 | -0.9664 | yes |  |
| g_largest_domain | 0.8285 | 0.0396 | 0.7889 [45.3] | 0.8038 | 0.7892 | 0.8078 | yes | R |

## qd:004  (COMBINATION EFFECT)

`C kC=2 | C.GAME(g=(1.5, 0.5)) ; C.MAJ(p=0.1)` · view C.transient · flagged True · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_activity_trend | -7.208 | ko absent | | ko absent | ko absent | -6.142 | yes | R |
| g_activity_cv | 0.5469 | ko absent | | ko absent | ko absent | 2.184 | yes | R |
| g_corrlen | 0 | ko absent | | ko absent | ko absent | 1 | no | R |
| g_n_types | 1 | ko absent | | ko absent | ko absent | 1 | yes | R |
| g_reversals | 1 | ko absent | | ko absent | ko absent | 1 | yes | R |
| g_largest_domain | 1 | ko absent | | ko absent | ko absent | 0.999 | yes | R |

## qd:005  (COMBINATION EFFECT)

`C kC=2 | C.GAME(g=(1.5, 0.5)) ; C.COPY(p=0.1)` · view C · flagged True · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_n_domains | 3.611 | 6.372 | -2.761 [-13.8] | -2.867 | -2.734 | 3.353 | yes | R |
| g_corrlen | 2 | 1 | 1 [16.7] | 1 | 1 | 2 | yes | R |
| g_n_types | 2 | 2 | 0 [0.0] | 0 | 0 | 2 | yes |  |
| g_reversals | 1 | 0 | 1 [50.0] | 1 | 1 | 1 | yes | R |
| g_largest_domain | 0.9397 | 0.0393 | 0.9004 [45.9] | 0.929 | 0.9053 | 0.9474 | yes | R |
| g_activity_trend | -0.625 | 0 | -0.625 [-50.0] | -0.5871 | -0.425 | -0.6036 | yes | R |

## qd:006  (DECOMPOSES)

`C kC=2 | C.GAME(g=(1.9, 0.0)) ; C.KURA(K=0.1)` · view C · flagged False · verdict KNOWN · reproduced by rules [1]

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_n_domains | 2.478 | 6.399 | -3.921 [-2.9] | -2.071 | -3.936 | 4.377 | yes |  |
| g_corrlen | 1 | 1 | 0 [0.0] | 1 | 0 | 2 | yes |  |
| g_n_types | 1.5 | 2 | -0.5 [-1.3] | 0 | -0.5 | 2 | yes |  |
| g_reversals | 1 | 0 | 1 [50.0] | 1 | 1 | 1 | yes | R |
| g_type_entropy | 0.45 | 1 | -0.5499 [-1.6] | -0.1 | -0.551 | 0.9045 | yes |  |
| g_largest_domain | 0.8577 | 0.0371 | 0.8206 [6.1] | 0.6552 | 0.8207 | 0.6888 | yes | R |

## qd:007  (DECOMPOSES)

`C kC=2 | C.GAME(g=(1.9, 0.0)) ; C.KURA(K=0.1)` · view C.phase · flagged True · verdict KNOWN · reproduced by rules [2]

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| p_local_r | 0.673 | 0.4508 | 0.2223 [9.9] | 0.219 | 0.2139 | 0.6715 | yes | R |
| p_r | 0.0538 | 0.0139 | 0.0399 [8.1] | 0.0165 | 0.042 | 0.03712 | yes |  |

## qd:008  (DECOMPOSES)

`N kN=2 | N.GAME(g=(1.5, 0.5)) ; N.GAME(g=(1.9, 0.0))` · view N · flagged True · verdict KNOWN · reproduced by rules [2]

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| n_components | 1.869 | 1.869 | 0 [0.0] | 0 | 0 | 1.657 | yes |  |
| n_giant | 0.9775 | 0.9775 | 0 [0.0] | 0 | 0 | 0.9833 | yes |  |
| n_degree_cv | 0.5022 | 0.5022 | 0 [0.0] | 0 | 0 | 0.4719 | yes |  |
| n_type_modularity | 0.1849 | -0.0057 | 0.1905 [8.8] | 0.1968 | 0.1904 | 0.2074 | yes | R |

## qd:009  (ARTIFACT)

`N kN=2 | N.GAME(g=(1.5, 0.5)) ; N.KURA(K=1.0)` · view N · flagged True · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| n_components | 1.946 | 1.946 | 0 [0.0] | 0 | 0 | 2.036 | yes |  |

## qd:010  (DECOMPOSES)

`N kN=2 | N.GAME(g=(1.9, 0.0)) ; N.KURA(K=0.1)` · view N · flagged True · verdict KNOWN · reproduced by rules [1]

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| n_components | 1.869 | 1.869 | 0 [0.0] | 0 | 0 | 1.461 | yes |  |
| n_giant | 0.975 | 0.975 | 0 [0.0] | 0 | 0 | 0.9883 | yes |  |
| n_degree_cv | 0.4916 | 0.4916 | 0 [0.0] | 0 | 0 | 0.4738 | yes |  |
| n_consensus_gain | 0.04 | 0 | 0.04 [10.8] | 0.0381 | 0.04 | 0.1867 | yes | R |
| n_oscillation | 0.09 | 0 | 0.09 [24.3] | 0.1658 | 0.0591 | 0.1613 | yes | R |

## qd:011  (DECOMPOSES)

`N kN=2 | N.GAME(g=(1.9, 0.0)) ; N.KURA(K=0.1)` · view N.phase · flagged False · verdict KNOWN · reproduced by rules [2]

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| p_freq_spread | 0.347 | 0.5007 | -0.1537 [-3.4] | -0.1448 | -0.5827 | 0.3225 | yes | R |
| p_r | 0.235 | 0.0609 | 0.1742 [10.2] | 0.0673 | 0.1739 | 0.2379 | yes |  |
| p_r_std | 0.082 | 0.0295 | 0.0525 [5.2] | 0.0187 | 0.0599 | 0.1053 | yes |  |

## qd:012  (COMBINATION EFFECT)

`C kC=2 | C.GAME(g=(1.5, 0.5)) ; C.SWAP(p=1.0)` · view C · flagged False · verdict KNOWN-VARIANT

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_n_domains | 5.201 | 6.364 | -1.163 [-5.0] | -1.241 | -1.078 | 5.191 | yes | R |
| g_n_types | 2 | 2 | 0 [0.0] | 0 | 0 | 2 | yes |  |
| g_corrlen | 1 | 1 | 0 [0.0] | 1 | 0 | 1 | yes |  |
| g_reversals | 1 | 0 | 1 [50.0] | 1 | 1 | 1 | yes | R |
| g_type_entropy | 1 | 0.9999 | 0.0001 [0.0] | -0.0001 | 0.0001 | 0.9999 | yes |  |
| g_activity | 0.4514 | 0 | 0.4514 [50.2] | 0.4542 | 0.4568 | 0.4573 | yes | R |

## qd:013  (COMBINATION EFFECT)

`C kC=2 | C.GAME(g=(1.5, -0.5)) ; C.MAJ(p=0.1)` · view C · flagged True · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_n_domains | 3.636 | 6.372 | -2.736 [-13.7] | -2.833 | -2.612 | 3.696 | yes | R |
| g_corrlen | 2 | 1 | 1 [16.7] | 1 | 1 | 2 | yes | R |
| g_n_types | 2 | 2 | 0 [0.0] | 0 | 0 | 2 | yes |  |
| g_activity_trend | 2.344 | 0 | 2.344 [50.0] | 2.799 | 0.731 | 1.882 | yes |  |
| g_reversals | 1 | 0 | 1 [50.0] | 1 | 1 | 1 | yes | R |
| g_type_entropy | 0.8205 | 0.9999 | -0.1794 [-4.9] | -0.1836 | -0.164 | 0.8271 | yes | R |

## qd:014  (ARTIFACT)

`N kN=2 | N.GAME(g=(1.9, 0.0)) ; N.REWIRE(p=0.01)` · view N · flagged True · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| n_components | 1.099 | 1.701 | -0.602 [-1.7] | -0.8491 | -0.602 | 1.155 | yes |  |

## qd:015  (COMBINATION EFFECT)

`N kN=2 | N.GAME(g=(1.9, 0.0)) ; N.REWIRE(p=0.003)` · view N · flagged True · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| n_components | 1.498 | 1.589 | -0.0912 [-0.3] | -0.318 | -0.2027 | 1.246 | yes |  |
| n_giant | 0.9875 | 0.985 | 0.0025 [0.1] | 0.0062 | 0.005 | 0.99 | yes |  |
| n_type_entropy | 0.4882 | 0.6916 | -0.2034 [-7.9] | -0.0049 | -0.5569 | 0.6889 | yes |  |
| n_degree_cv | 0.5186 | 0.502 | 0.0166 [0.5] | 0.0216 | 0.0362 | 0.5019 | yes |  |
| n_oscillation | 0.3337 | 0 | 0.3337 [18.5] | 0.3341 | 0.4021 | 0.4213 | yes | R |
| n_edge_turnover | 0.3355 | 0 | 0.3355 [28.2] | 0.3645 | 0.4734 | 0.3488 | yes | R |

## qd:016  (COMBINATION EFFECT)

`C kC=2 | C.GAME(g=(1.5, 0.5)) ; C.MAJ(p=0.01)` · view C · flagged True · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_activity_trend | -4.208 | 0 | -4.208 [-50.0] | -4.216 | -2.744 | -4.22 | yes | R |
| g_largest_domain | 1 | 0.056 | 0.944 [44.7] | 0.9802 | 0.944 | 1 | yes | R |
| g_n_types | 1 | 2 | -1 [-16.7] | -1 | -1 | 1 | yes | R |
| g_reversals | 1 | 0 | 1 [50.0] | 1 | 1 | 1 | yes | R |
| g_n_domains | 0.6931 | 6.355 | -5.662 [-40.2] | -7.01 | -5.662 | 0.6932 | yes | R |
| g_spec_peak | 0.0931 | 0 | 0.0931 [23.9] | 0.1118 | 0.0938 | 0.08861 | yes | R |

## qd:017  (DECOMPOSES)

`C kC=2 | C.GAME(g=(1.5, 0.5)) ; C.SWAP(p=0.03)` · view C · flagged True · verdict KNOWN · reproduced by rules [1]

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_n_domains | 4.5 | 6.351 | -1.851 [-8.5] | -1.953 | -1.891 | 4.339 | yes | R |
| g_corrlen | 2 | 1 | 1 [16.7] | 1 | 1 | 2 | yes | R |
| g_n_types | 2 | 2 | 0 [0.0] | 0 | 0 | 2 | yes |  |
| g_reversals | 1 | 0 | 1 [50.0] | 1 | 1 | 1 | yes | R |
| g_largest_domain | 0.8756 | 0.0454 | 0.8302 [45.1] | 0.8601 | 0.8243 | 0.8856 | yes | R |
| g_moran | 0.5464 | 0.0072 | 0.5392 [48.6] | 0.5502 | 0.5463 | 0.5595 | yes | R |

## qd:018  (DECOMPOSES)

`C kC=2 | C.GAME(g=(1.5, 0.5)) ; C.KURA(K=0.2)` · view C.phase · flagged True · verdict KNOWN · reproduced by rules [2]

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| p_local_r | 0.9073 | 0.4513 | 0.456 [16.8] | 0.46 | 0.4596 | 0.9069 | yes | R |
| p_r | 0.1269 | 0.0149 | 0.1119 [17.8] | 0.0462 | 0.1059 | 0.1149 | yes |  |

## qd:019  (DECOMPOSES)

`N kN=2 | N.GAME(g=(1.5, 0.5)) ; N.KURA(K=0.2)` · view N.phase · flagged False · verdict KNOWN · reproduced by rules [2]

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| p_r | 0.8423 | 0.062 | 0.7803 [43.1] | 0.7891 | 0.7761 | 0.8403 | yes | R |
| p_freq_spread | 0.0836 | 0.5051 | -0.4214 [-14.9] | -0.4155 | -1.226 | 0.07791 | yes | R |
| p_r_std | 0.0074 | 0.0255 | -0.0181 [-4.9] | -0.011 | -0.024 | 0.01081 | yes | R |

## qd:020  (DECOMPOSES)

`C kC=2 | C.GAME(g=(1.5, 0.5)) ; C.CYCLE(th=7,p=0.1)` · view C · flagged True · verdict KNOWN · reproduced by rules [1]

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_n_domains | 4.29 | 6.388 | -2.097 [-9.8] | -1.925 | -2.104 | 4.375 | yes | R |
| g_corrlen | 2 | 1 | 1 [16.7] | 1 | 1 | 2 | yes | R |
| g_n_types | 2 | 2 | 0 [0.0] | 0 | 0 | 2 | yes |  |
| g_reversals | 1 | 0 | 1 [50.0] | 1 | 1 | 1 | yes | R |
| g_largest_domain | 0.9111 | 0.0487 | 0.8624 [44.9] | 0.8728 | 0.8605 | 0.8978 | yes | R |
| g_activity_trend | 0.6782 | 0 | 0.6782 [49.9] | 0.8511 | 0.1665 | 0.6621 | yes |  |

## qd:021  (ARTIFACT)

`C kC=2 | C.GAME(g=(1.5, 0.5)) ; C.CYCLE(th=3,p=0.03)` · view C · flagged True · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_n_types | 2 | 2 | 0 [0.0] | 0 | -1 | 1.333 | no |  |

## qd:022  (ARTIFACT, flagged view did not recur on fresh seeds)

`C kC=2 | C.GAME(g=(1.5, -0.5)) ; C.SWAP(p=0.002)` · view C.transient · flagged False · verdict KNOWN-VARIANT

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_activity_cv | view absent | | | | | 5.296 |  | |
| g_largest_domain | view absent | | | | | 1 |  | |
| g_n_types | view absent | | | | | 1 |  | |
| g_reversals | view absent | | | | | 1 |  | |
| g_n_domains | view absent | | | | | 0.6932 |  | |
| g_spec_freq | view absent | | | | | 0.2732 |  | |

## qd:023  (COMBINATION EFFECT)

`C kC=2 | C.GAME(g=(1.5, -0.5)) ; C.COPY(p=0.002)` · view C · flagged False · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_corrlen | 0 | 1 | -1 [-50.0] | 4 | -1 | 4 | no |  |
| g_activity_cv | 0 | 0 | 0 [0.0] | 1.846 | 0 | 2.717 | no |  |
| g_n_types | 1 | 2 | -1 [-16.7] | -1 | -1 | 1.667 | no | R |
| g_n_domains | 0.6931 | 6.379 | -5.686 [-40.2] | -5.926 | -5.686 | 1.195 | no | R |
| g_reversals | 1 | 0 | 1 [50.0] | 1 | 1 | 1 | yes | R |
| g_largest_domain | 1 | 0.0425 | 0.9575 [46.0] | 0.9746 | 0.9575 | 0.9908 | yes | R |

## qd:024  (DECOMPOSES)

`C kC=2 | C.SWAP(p=0.001) ; C.GAME(g=(1.5, -0.5))` · view C · flagged True · verdict KNOWN · reproduced by rules [2]

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_corrlen | 0 | 1 | -1 [-50.0] | 2 | -1 | 3 | no |  |
| g_activity_cv | 0 | 0 | 0 [0.0] | 4.791 | 0 | 2.828 | no |  |
| g_n_domains | 0.6931 | 6.384 | -5.691 [-40.2] | -6.338 | -5.691 | 1.589 | no | R |
| g_n_types | 1 | 2 | -1 [-16.7] | -1 | -1 | 1.5 | no | R |
| g_largest_domain | 1 | 0.0399 | 0.9601 [46.2] | 0.9799 | 0.9601 | 0.9915 | yes | R |
| g_moran | 0 | -0.0043 | 0.0043 [0.9] | 0.6651 | 0.0043 | 0.6458 | no |  |

## qd:025  (DECOMPOSES)

`C kC=2 | C.GAME(g=(1.5, 0.5)) ; C.FLIP(a=0,b=1,p=0.0003)` · view C · flagged True · verdict KNOWN · reproduced by rules [1]

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_n_domains | 4.945 | 6.371 | -1.426 [-6.3] | -1.434 | -1.31 | 4.933 | yes | R |
| g_corrlen | 2 | 1 | 1 [16.7] | 1 | 1 | 2 | yes | R |
| g_n_types | 2 | 2 | 0 [0.0] | 0 | 0 | 2 | yes |  |
| g_reversals | 1 | 0 | 1 [50.0] | 1 | 1 | 1 | yes | R |
| g_type_entropy | 0.7987 | 0.9641 | -0.1654 [-4.7] | -0.1718 | -0.0026 | 0.7918 | yes |  |
| g_largest_domain | 0.7556 | 0.447 | 0.3086 [8.1] | 0.2157 | -0.0249 | 0.7453 | yes |  |

## qd:026  (DECOMPOSES)

`C kC=2 | C.GAME(g=(1.5, 0.5)) ; C.CYCLE(th=7,p=0.05)` · view C · flagged True · verdict KNOWN · reproduced by rules [1]

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_n_types | 2 | 2 | 0 [0.0] | 0 | 0 | 2 | yes |  |
| g_reversals | 1 | 0 | 1 [50.0] | 1 | 1 | 1 | yes | R |

## qd:027  (COMBINATION EFFECT)

`AF kA=1 v=0.3 eta=0.3 fields=1 D=[0.02] | A.CLIMB(f=0) ; A.ALIGN(r=1.0)` · view A · flagged True · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| a_density_cv | 8.449 | 0.796 | 7.653 [39.1] | 4.203 | 14.51 | 10.57 | yes | R |
| a_local_align | 0.9106 | -0.0351 | 0.9457 [46.4] | 0.942 | 0.9817 | 0.914 | yes | R |
| a_polar | 0.5581 | 0.0362 | 0.5219 [5.8] | 0.2265 | 0.7264 | 0.6157 | yes |  |
| a_ang_mom | 0.6226 | 0.0309 | 0.5917 [8.1] | 0.1389 | 0.0169 | 0.3549 | yes |  |
| a_step | 0.098 | 5.925 | -5.827 [-48.4] | -5.636 | -17.13 | 0.2632 | no | R |
| a_polar_std | 0.1204 | 0.0133 | 0.1071 [22.3] | 0.1067 | 0.1249 | 0.1959 | yes | R |

## qd:028  (ARTIFACT)

`CA kC=1 kA=3 v=0.3 eta=0.0 | A.ALIGN(r=1.0) ; A.WRITECELL()` · view A · flagged True · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| a_type_entropy | 1.097 | 1.097 | 0 [0.0] | 0 | 0 | 1.095 | yes |  |
| a_type_segregation | 0.3362 | 0.3352 | 0.001 [0.1] | -0.0019 | 0.0077 | 0.3302 | yes |  |

## qd:029  (ARTIFACT)

`CA kC=3 kA=1 v=0.3 eta=3.14159 | C.COPY(p=1.0) ; A.WRITECELL()` · view C · flagged True · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_n_domains | 7.106 | 7.346 | -0.2402 [-0.8] | -0.2415 | -0.2435 | 7.107 | yes |  |
| g_corrlen | 1 | 1 | 0 [0.0] | 0 | 0 | 1 | yes |  |

## qd:030  (ARTIFACT)

`CA kC=3 kA=2 v=0.0 eta=0.3 | A.WRITECELL() ; C.MAJ(p=1.0)` · view C · flagged True · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_n_types | 3 | 3 | 0 [0.0] | 0 | 0 | 3 | yes |  |

## qd:031  (COMBINATION EFFECT)

`C kC=2 | C.GAME(g=(1.5, -0.5)) ; C.CYCLE(th=6,p=0.03)` · view C · flagged True · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_corrlen | 3 | 1 | 2 [25.0] | 2 | 2 | 3 | yes | R |
| g_activity_cv | 1.739 | 0 | 1.739 [50.0] | 1.26 | 1.404 | 2.96 | yes | R |
| g_n_domains | 1.609 | 6.371 | -4.761 [-29.8] | -5.522 | -4.761 | 1.099 | yes | R |
| g_n_types | 1 | 2 | -1 [-16.7] | -1 | -1 | 1 | yes | R |
| g_reversals | 1 | 0 | 1 [50.0] | 1 | 1 | 1 | yes | R |
| g_largest_domain | 0.9904 | 0.0533 | 0.937 [44.8] | 0.9743 | 0.937 | 0.9968 | yes | R |

## qd:032  (DECOMPOSES)

`CA kC=4 kA=1 v=0.0 eta=0.3 | C.SCHELL(th=0.3) ; A.WRITECELL()` · view C · flagged True · verdict KNOWN · reproduced by rules [1]

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_reversals | 0.3696 | 0 | 0.3696 [49.9] | 0.3722 | 0.3777 | 0.3757 | yes | R |

## qd:033  (COMBINATION EFFECT)

`CAF kC=1 kA=1 v=0.3 eta=0.3 fields=1 D=[0.05] | A.CLIMB(f=0) ; A.TURNCELL()` · view A · flagged True · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| a_step | 5.913 | 5.925 | -0.012 [-0.1] | -0.0103 | 0.135 | 5.898 | yes |  |
| a_density_cv | 1.159 | 0.8089 | 0.3496 [6.3] | 0.4652 | 0.5886 | 1.257 | yes | R |
| a_nn_ratio | 0.826 | 0.9926 | -0.1666 [-4.6] | -0.2046 | -0.2699 | 0.7857 | yes | R |
| a_ang_mom | 0.3603 | 0.0193 | 0.341 [16.6] | 0.1588 | 0.5122 | 0.3198 | yes |  |
| a_local_align | 0.2566 | -0.0405 | 0.2971 [11.5] | 0.3044 | 0.3409 | 0.269 | yes | R |
| a_polar | 0.0547 | 0.0388 | 0.0159 [1.6] | 0.0197 | 0.0109 | 0.05804 | yes |  |

## qd:034  (DECOMPOSES)

`CA kC=3 kA=2 v=0.0 eta=1.0 | C.KURA(K=0.1) ; A.WRITECELL()` · view C.phase · flagged True · verdict KNOWN · reproduced by rules [1]

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| p_local_r | 0.6674 | 0.4485 | 0.2189 [9.8] | 0.2216 | 0.2142 | 0.6731 | yes | R |
| p_r | 0.041 | 0.0136 | 0.0274 [5.4] | 0.0176 | 0.0482 | 0.06746 | yes | R |

## qd:035  (DECOMPOSES)

`AF kA=3 v=2.0 eta=0.3 fields=1 D=[0.02] | A.CLIMB(f=0) ; A.ALIGN(r=3.0)` · view A · flagged True · verdict KNOWN · reproduced by rules [1]

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| a_type_entropy | 1.097 | 1.097 | 0 [0.0] | 0 | 0 | 1.098 | yes |  |
| a_step | 0.5651 | 32.02 | -31.46 [-48.3] | -38.91 | -25.55 | 0.7002 | yes | R |
| a_type_segregation | 0.3283 | 0.3374 | -0.0091 [-0.7] | -0.0002 | 0.0072 | 0.3231 | yes |  |

## qd:036  (COMBINATION EFFECT)

`AF kA=1 v=0.3 eta=0.3 fields=1 D=[0.02] | A.ALIGN(r=3.0) ; A.DEPOSIT(f=0,p=0.1)` · view F0 · flagged False · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| f_mean | 10.75 | 10.75 | 0 [0.0] | 0 | 0 | 10.75 | yes |  |
| f_std | 4.251 | 0.4394 | 3.812 [32.7] | 5.018 | 3.47 | 2.813 | yes | R |
| f_n_blobs | 1.099 | 4.025 | -2.926 [-5.5] | -3.395 | -2.627 | 1.461 | yes | R |
| f_moran | 0.9953 | 0.9121 | 0.0833 [2.2] | 0.0634 | 0.0596 | 0.9905 | yes |  |
| f_peak_sharp | 0.8775 | 0.348 | 0.5295 [6.8] | 0.2122 | 0.3841 | 0.8027 | yes |  |
| f_change | 0.1182 | 0.1178 | 0.0004 [0.1] | 0.0032 | 0.0317 | 0.1169 | yes |  |

## qd:037  (COMBINATION EFFECT)

`CA kC=3 kA=4 v=0.0 eta=1.0 | C.SCHELL(th=0.2) ; A.WRITECELL()` · view C · flagged True · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_n_domains | 6.081 | 7.366 | -1.285 [-4.8] | -1.292 | -1.288 | 6.073 | yes | R |
| g_n_types | 3 | 3 | 0 [0.0] | 0 | 0 | 3 | yes |  |
| g_corrlen | 1 | 1 | 0 [0.0] | 0 | 0 | 1 | yes |  |
| g_type_entropy | 0.7364 | 0.9998 | -0.2634 [-7.6] | -0.263 | -0.2624 | 0.7373 | yes | R |
| g_reversals | 0.5455 | 0 | 0.5455 [50.0] | 0.5525 | 0.5758 | 0.5472 | yes | R |
| g_interface | 0.3748 | 0.6701 | -0.2953 [-14.1] | -0.2916 | -0.2982 | 0.3788 | yes | R |

## qd:038  (DECOMPOSES)

`CA kC=3 kA=3 v=0.0 eta=0.0 | A.ALIGN(r=3.0) ; A.WRITECELL()` · view A · flagged True · verdict KNOWN · reproduced by rules [1]

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| a_nn_ratio | 0.9914 | 0.9914 | 0 [0.0] | 0 | 0 | 0.9998 | yes |  |
| a_density_cv | 0.799 | 0.799 | 0 [0.0] | 0 | 0 | 0.7566 | yes |  |
| a_type_segregation | 0.3353 | 0.3353 | 0 [0.0] | 0 | 0 | 0.3287 | yes |  |
| a_local_align | 0.2017 | -0.0323 | 0.234 [19.0] | 0.2667 | 0.234 | 0.1944 | yes | R |

## qd:039  (COMBINATION EFFECT)

`C kC=2 | C.GAME(g=(1.5, -0.5)) ; C.CONTAG(a=0,b=1,th=4,p=0.03)` · view C · flagged True · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_corrlen | 6.5 | 1 | 5.5 [14.8] | 10 | 7.5 | 6.5 | yes | R |
| g_n_types | 2 | 2 | 0 [0.0] | 0 | 0 | 2 | yes |  |
| g_n_domains | 1.099 | 6.404 | -5.305 [-35.4] | -5.574 | -5.161 | 1.099 | yes | R |
| g_activity_cv | 0.8042 | 0 | 0.8042 [50.0] | 0.5086 | 0.5431 | 1.041 | yes | R |
| g_reversals | 1 | 0 | 1 [50.0] | 1 | 1 | 1 | yes | R |
| g_largest_domain | 0.9661 | 0.0433 | 0.9227 [43.1] | 0.9293 | 0.8499 | 0.9706 | yes | R |

## qd:040  (COMBINATION EFFECT)

`C kC=2 | C.GAME(g=(1.3, 0.0)) ; C.CONTAG(a=1,b=0,th=2,p=0.03)` · view C · flagged True · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_n_domains | 3.135 | 6.374 | -3.239 [-17.0] | -3.365 | -3.543 | 3.137 | yes | R |
| g_corrlen | 3 | 1 | 2 [25.0] | 2 | 1.5 | 3 | yes | R |
| g_n_types | 2 | 2 | 0 [0.0] | 0 | 0 | 2 | yes |  |
| g_reversals | 1 | 0 | 1 [50.0] | 1 | 1 | 1 | yes | R |
| g_largest_domain | 0.9426 | 0.0424 | 0.9003 [45.7] | 0.9197 | 0.8978 | 0.9381 | yes | R |
| g_moran | 0.5727 | 0.0001 | 0.5726 [33.9] | 0.5692 | 0.5639 | 0.578 | yes | R |

## qd:041  (COMBINATION EFFECT)

`CAF kC=2 kA=1 v=0.3 eta=0.3 fields=1 D=[0.05] | A.CLIMB(f=0) ; A.WRITECELL()` · view C · flagged True · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_n_domains | 6.357 | 6.367 | -0.0104 [-0.0] | -0.013 | -0.0052 | 6.383 | yes |  |
| g_activity_trend | -4.511 | 0 | -4.511 [-50.0] | -4.409 | -3.596 | -4.272 | yes | R |
| g_n_types | 2 | 2 | 0 [0.0] | 0 | 0 | 2 | yes |  |
| g_activity_cv | 1.521 | 0 | 1.521 [27.7] | 0.8032 | 3.937 | 1.493 | yes | R |
| g_corrlen | 1 | 1 | 0 [0.0] | 0 | 0 | 1 | yes |  |
| g_reversals | 1 | 0 | 1 [50.0] | 1 | 1 | 1 | yes | R |

## qd:042  (DECOMPOSES)

`AF kA=1 v=0.3 eta=0.3 fields=1 D=[0.02] | A.ALIGN(r=0.5) ; A.DEPOSIT(f=0,p=0.1)` · view A · flagged False · verdict KNOWN · reproduced by rules [1]

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| a_step | 5.904 | 5.928 | -0.0237 [-0.1] | -0.0365 | 0.0811 | 5.865 | yes |  |
| a_density_cv | 2.209 | 0.7859 | 1.424 [17.4] | 0.9117 | 1.06 | 2.262 | yes | R |
| a_local_align | 0.4889 | -0.0296 | 0.5185 [42.5] | 0.5478 | 0.4795 | 0.4978 | yes | R |
| a_nn_ratio | 0.3733 | 1.016 | -0.6426 [-12.9] | -0.6341 | -0.561 | 0.3729 | yes | R |
| a_polar | 0.6381 | 0.0317 | 0.6064 [16.7] | 0.242 | 0.7669 | 0.3101 | yes |  |
| a_ang_mom | 0.1225 | 0.0298 | 0.0927 [5.0] | 0.047 | 0.0508 | 0.09465 | yes | R |

## qd:043  (DECOMPOSES)

`CA kC=3 kA=1 v=0.5 eta=0.0 | C.MAJ(p=0.001) ; A.TURNCELL()` · view C · flagged True · verdict KNOWN · reproduced by rules [1]

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_n_domains | 6.304 | 7.35 | -1.045 [-3.8] | -1.065 | -2.333 | 6.308 | yes | R |
| g_n_types | 3 | 3 | 0 [0.0] | 0 | 0 | 3 | yes |  |
| g_corrlen | 2 | 1 | 1 [16.7] | 1 | 2 | 2 | yes | R |
| g_type_entropy | 0.9993 | 0.9998 | -0.0005 [-0.0] | -0.0001 | -0.0019 | 0.9994 | yes |  |
| g_reversals | 0.6667 | 0 | 0.6667 [50.1] | 0.625 | 0.7183 | 0.8333 | yes | R |
| g_activity_trend | -0.5672 | 0 | -0.5672 [-50.2] | -0.7057 | -1.199 | -0.7058 | yes | R |

## qd:044  (COMBINATION EFFECT)

`C kC=2 | C.COPY(p=0.3) ; C.MAJ(p=0.2) ; C.GAME(g=(1.5, 0.5))` · view C · flagged True · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_n_domains | 3.795 | 6.387 | -2.591 [-12.7] | -2.69 | -2.832 | 3.385 | yes | R |
| g_corrlen | 2 | 1 | 1 [16.7] | 1 | 1 | 2.667 | yes | R |
| g_n_types | 2 | 2 | 0 [0.0] | 0 | 0 | 2 | yes |  |
| g_reversals | 1 | 0 | 1 [50.0] | 1 | 1 | 1 | yes | R |
| g_largest_domain | 0.9025 | 0.0444 | 0.858 [45.4] | 0.9005 | 0.8796 | 0.9308 | yes | R |
| g_moran | 0.6071 | -0.002 | 0.609 [49.9] | 0.6144 | 0.6174 | 0.6384 | yes | R |

## qd:045  (COMBINATION EFFECT)

`CF kC=4 fields=2 D=[0.05, 0.02] | C.COPY(p=0.1) ; C.SENSE(a=0,b=2,f=0,th=0.5)` · view C · flagged False · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_n_domains | 5.924 | 7.681 | -1.757 [-6.5] | -1.809 | -1.884 | 5.942 | yes | R |
| g_corrlen | 3 | 1 | 2 [25.0] | 2 | 2 | 3.667 | yes | R |
| g_n_types | 3 | 4 | -1 [-7.1] | -1 | -1 | 3 | yes | R |
| g_reversals | 0.7641 | 0 | 0.7641 [49.9] | 0.7688 | 0.7499 | 0.7723 | yes | R |
| g_type_entropy | 0.7429 | 0.9998 | -0.2568 [-7.4] | -0.26 | -0.258 | 0.7435 | yes | R |
| g_moran | 0.5183 | -0.0033 | 0.5217 [36.2] | 0.5062 | 0.5397 | 0.5449 | yes | R |

## qd:046  (COMBINATION EFFECT)

`CF kC=4 fields=2 D=[0.05, 0.02] | C.COPY(p=0.2) ; C.EMIT(b=0,f=1,p=1.0)` · view F1 · flagged False · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| f_mean | 232.7 | 250 | -17.36 [-1.8] | 13.53 | -30.23 | 274.2 | yes |  |
| f_std | 93.58 | 30.83 | 62.75 [11.9] | 101.3 | 141.4 | 127.1 | yes | R |
| f_change | 2.694 | 2.5 | 0.1937 [1.1] | 0.4767 | 0.2104 | 3.065 | yes |  |
| f_n_blobs | 1.609 | 3.651 | -2.041 [-8.1] | -2.335 | -1.497 | 1.561 | yes | R |
| f_moran | 0.9925 | 0.9439 | 0.0485 [1.3] | 0.0458 | 0.0219 | 0.9925 | yes |  |
| f_peak_sharp | 0.7451 | 0.3226 | 0.4224 [7.0] | 0.2276 | 0.3556 | 0.755 | yes | R |

## qd:047  (ARTIFACT)

`CAF kC=1 kA=3 v=0.1 eta=0.1 fields=1 D=[0.05] | F.FEED(f=0,F=0.03) ; A.CLIMB(f=0)` · view A · flagged True · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| a_density_cv | 6.163 | 0.8069 | 5.356 [38.4] | 2.863 | 0.055 | 8.153 | yes |  |
| a_step | 1.997 | 1.997 | -0 [-0.0] | -0.0261 | -0.0006 | 1.997 | yes |  |
| a_type_entropy | 1.096 | 1.096 | 0 [0.0] | 0 | 0 | 1.096 | yes |  |
| a_local_align | 0.8227 | -0.0321 | 0.8549 [50.0] | 0.7447 | 0.0062 | 0.8289 | yes |  |
| a_ang_mom | 0.1562 | 0.0288 | 0.1274 [2.8] | 0.0612 | 0.0269 | 0.4325 | yes |  |
| a_type_segregation | 0.3355 | 0.3318 | 0.0038 [0.3] | 0.0007 | -0.005 | 0.3308 | yes |  |

## qd:048  (DECOMPOSES)

`CA kC=3 kA=1 v=0.3 eta=3.14159 | A.ATTRACT(b=any,r=3.0) ; C.COPY(p=1.0) ; A.TURNCELL()` · view C · flagged True · verdict KNOWN · reproduced by rules [2]

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_n_domains | 6.434 | 7.357 | -0.9238 [-3.3] | -0.9242 | -1.548 | 6.288 | yes | R |
| g_corrlen | 1 | 1 | 0 [0.0] | 0 | 0 | 1 | yes |  |

## qd:049  (ARTIFACT)

`CA kC=3 kA=2 v=0.0 eta=0.0 | C.CYCLE(th=6,p=0.003) ; A.TURNCELL()` · view C · flagged True · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_n_domains | 7.278 | 7.34 | -0.0623 [-0.2] | -0.0647 | -0.0682 | 7.293 | yes |  |
| g_activity_cv | 3.549 | 0 | 3.549 [12.4] | 1.866 | 0 | 3.42 | yes |  |
| g_n_types | 3 | 3 | 0 [0.0] | 0 | 0 | 3 | yes |  |
| g_corrlen | 1 | 1 | 0 [0.0] | 0 | 0 | 1 | yes |  |
| g_type_entropy | 0.9998 | 0.9999 | -0.0001 [-0.0] | -0 | -0.0001 | 0.9998 | yes |  |
| g_interface | 0.6355 | 0.6635 | -0.0279 [-1.1] | -0.0295 | -0.0313 | 0.6417 | yes |  |

## qd:050  (COMBINATION EFFECT)

`CF kC=4 fields=2 D=[0.05, 0.02] | C.SENSE(a=0,b=3,f=1,th=0.05) ; C.SCHELL(th=0.3)` · view C · flagged True · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_n_domains | 5.595 | 7.667 | -2.072 [-7.8] | -1.191 | -2.074 | 5.537 | yes | R |
| g_n_types | 4 | 4 | 0 [0.0] | 0 | 0 | 3.333 | yes |  |
| g_activity_trend | -2.332 | 0 | -2.332 [-10.8] | -3.679 | -1.39 | -3.316 | yes | R |
| g_corrlen | 2 | 1 | 1 [16.7] | 0 | 1 | 2 | yes |  |
| g_type_entropy | 0.8523 | 0.9998 | -0.1475 [-4.0] | -0.0061 | -0.1475 | 0.801 | yes |  |
| g_reversals | 0.767 | 0 | 0.767 [48.5] | 0.4686 | 0.6482 | 0.6846 | yes | R |

## qd:051  (ARTIFACT)

`CF kC=3 fields=2 D=[0.05, 0.02] | C.SCHELL(th=0.2) ; C.SENSE(a=0,b=2,f=1,th=0.05)` · view C · flagged True · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_n_domains | 6.902 | 7.352 | -0.4496 [-1.6] | -0.2425 | -0.4496 | 6.813 | yes |  |

## qd:052  (COMBINATION EFFECT)

`CF kC=2 fields=2 D=[0.02, 0.01] | F.AUTOCAT(f=0,g=1,p=1.0) ; C.SENSE(a=0,b=1,f=1,th=0.8)` · view C · flagged True · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_largest_domain | 1 | 0.0405 | 0.9595 [46.1] | 0.9829 | 0.9595 | 1 | yes | R |
| g_n_types | 1 | 2 | -1 [-16.7] | -1 | -1 | 1 | yes | R |

## qd:053  (COMBINATION EFFECT)

`CF kC=4 fields=2 D=[0.05, 0.02] | C.MAJ(p=0.001) ; C.SENSE(a=0,b=1,f=0,th=0.5)` · view C · flagged False · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_n_domains | 6.271 | 7.672 | -1.401 [-5.0] | -1.426 | -2.723 | 6.245 | yes | R |
| g_n_types | 3 | 4 | -1 [-7.1] | -1 | -1 | 3 | yes | R |
| g_corrlen | 1 | 1 | 0 [0.0] | 0 | 1 | 1 | yes |  |
| g_activity_trend | -0.6227 | 0 | -0.6227 [-49.8] | -0.5981 | -1.272 | -0.7312 | yes | R |
| g_type_entropy | 0.6487 | 0.9998 | -0.3511 [-10.6] | -0.3465 | -0.5449 | 0.6639 | yes | R |
| g_activity_cv | 0.7292 | 0 | 0.7292 [49.9] | 0.3466 | 1.096 | 0.6297 | yes |  |

## qd:054  (DECOMPOSES)

`CAF kC=4 kA=2 v=0.3 eta=0.3 fields=1 D=[0.05] | A.CLIMB(f=0) ; C.CYCLE(th=3,p=0.05)` · view C · flagged True · verdict KNOWN · reproduced by rules [2]

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_n_domains | 4.871 | 7.678 | -2.806 [-11.2] | -2.849 | -3.387 | 4.814 | yes | R |
| g_corrlen | 5 | 1 | 4 [33.3] | 4.5 | 6 | 4.333 | yes | R |
| g_n_types | 4 | 4 | 0 [0.0] | 0 | 0 | 4 | yes |  |
| g_type_entropy | 0.9795 | 0.9997 | -0.0202 [-0.5] | -0.0089 | -0.0152 | 0.9483 | yes |  |
| g_moran | 0.7419 | -0.0063 | 0.7481 [49.9] | 0.754 | 0.8271 | 0.7261 | yes | R |
| g_cyclic_changes | 0.653 | -0.3333 | 0.9863 [50.1] | 0.9861 | 0.9773 | 0.6534 | yes | R |

## qd:055  (DECOMPOSES)

`CA kC=3 kA=3 v=0.1 eta=0.3 | C.SCHELL(th=0.7) ; C.MAJ(p=1.0) ; A.WRITECELL()` · view C · flagged True · verdict KNOWN · reproduced by rules [2, 3]

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_n_domains | 5.746 | 7.359 | -1.613 [-6.2] | -1.592 | -1.607 | 5.724 | yes | R |
| g_n_types | 3 | 3 | 0 [0.0] | 0 | -1 | 2.333 | no |  |
| g_corrlen | 1 | 1 | 0 [0.0] | 0 | 0 | 1 | yes |  |
| g_largest_domain | 0.7854 | 0.0079 | 0.7775 [48.9] | 0.7639 | 0.8948 | 0.9043 | yes | R |
| g_reversals | 0.8437 | 0 | 0.8437 [49.9] | 0.8582 | 0.8538 | 0.8823 | yes | R |
| g_activity_trend | -0.2336 | 0 | -0.2336 [-16.5] | -0.2988 | -0.5192 | -0.5338 | yes | R |

## qd:056  (COMBINATION EFFECT)

`CA kC=3 kA=3 v=0.1 eta=0.1 | A.ALIGN(r=1.0) ; C.MAJ(p=1.0) ; A.WRITECELL()` · view C · flagged True · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_reversals | 0.5807 | 0 | 0.5807 [50.1] | 0.5789 | 0.71 | 0.5852 | yes | R |

## qd:057  (COMBINATION EFFECT)

`CA kC=3 kA=1 v=0.5 eta=0.0 | A.REPEL(b=any,r=1.0) ; C.SCHELL(th=0.4) ; A.WRITECELL()` · view C · flagged True · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_n_domains | 6.304 | 7.357 | -1.052 [-3.9] | -1.068 | -1.064 | 6.302 | yes | R |
| g_n_types | 3 | 3 | 0 [0.0] | 0 | 0 | 3 | yes |  |
| g_corrlen | 1 | 1 | 0 [0.0] | 0 | 0 | 1 | yes |  |
| g_type_entropy | 0.8884 | 0.9999 | -0.1115 [-2.9] | -0.1111 | -0.1117 | 0.888 | yes |  |
| g_reversals | 0.7029 | 0 | 0.7029 [49.9] | 0.7068 | 0.7418 | 0.7065 | yes | R |
| g_largest_domain | 0.3966 | 0.0077 | 0.3889 [48.0] | 0.3625 | 0.4049 | 0.4021 | yes | R |

## qd:058  (DECOMPOSES)

`CF kC=4 fields=2 D=[0.02, 0.01] | C.MAJ(p=0.005) ; C.EMIT(b=0,f=0,p=0.002)` · view C · flagged False · verdict KNOWN · reproduced by rules [1]

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_n_domains | 4.824 | 7.681 | -2.857 [-11.4] | -2.92 | -3.711 | 4.817 | yes | R |
| g_n_types | 4 | 4 | 0 [0.0] | 0 | 0 | 4 | yes |  |
| g_corrlen | 3 | 1 | 2 [25.0] | 2 | 4 | 3.333 | yes | R |
| g_activity_trend | -1.659 | 0 | -1.659 [-50.0] | -1.627 | -1.796 | -1.59 | yes | R |
| g_type_entropy | 0.9938 | 0.9997 | -0.006 [-0.2] | -0.0012 | -0.014 | 0.9865 | yes |  |
| g_moran | 0.7506 | -0.0074 | 0.758 [49.9] | 0.7513 | 0.8426 | 0.7432 | yes | R |

## qd:059  (DECOMPOSES)

`CAF kC=4 kA=1 v=2.0 eta=0.3 fields=1 D=[0.02] | C.CYCLE(th=1,p=0.2) ; A.TURNCELL()` · view C · flagged False · verdict KNOWN · reproduced by rules [1]

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_n_domains | 7.37 | 7.663 | -0.2928 [-1.0] | -0.2903 | -0.2909 | 7.376 | yes |  |
| g_n_types | 4 | 4 | 0 [0.0] | 0 | 0 | 4 | yes |  |
| g_corrlen | 1 | 1 | 0 [0.0] | 0 | 0 | 1 | yes |  |
| g_type_entropy | 1 | 0.9998 | 0.0002 [0.0] | 0.0001 | 0.0002 | 1 | yes |  |
| g_interface | 0.6588 | 0.7485 | -0.0897 [-3.2] | -0.0928 | -0.0908 | 0.658 | yes | R |
| g_reversals | 0.2935 | 0 | 0.2935 [49.7] | 0.2932 | 0.3356 | 0.2942 | yes | R |

## qd:060  (COMBINATION EFFECT)

`CAF kC=2 kA=1 v=0.3 eta=3.14159 fields=1 D=[0.05] | A.CLIMB(f=0) ; A.WRITECELL() ; C.GAME(g=(1.5, 0.5))` · view C · flagged True · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_n_domains | 4.466 | 6.348 | -1.882 [-8.7] | -1.917 | -1.894 | 4.527 | yes | R |
| g_activity_trend | -2.666 | 0 | -2.666 [-31.5] | -2.941 | -2.776 | -2.677 | yes | R |
| g_corrlen | 2 | 1 | 1 [16.7] | 1 | 1 | 2 | yes | R |
| g_n_types | 2 | 2 | 0 [0.0] | 0 | 0 | 2 | yes |  |
| g_reversals | 1 | 0 | 1 [50.0] | 1 | 1 | 1 | yes | R |
| g_largest_domain | 0.864 | 0.0569 | 0.8071 [43.9] | 0.8517 | 0.8192 | 0.8648 | yes | R |

## qd:061  (DECOMPOSES)

`C kC=2 | C.FLIP(a=1,b=0,p=0.2) ; C.KURA(K=1.0) ; C.GAME(g=(0.5, 0.5)) ; C.MAJ(p=1.0)` · view C.phase · flagged False · verdict KNOWN · reproduced by rules [2]

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| p_local_r | 0.993 | 0.451 | 0.542 [18.8] | 0.5407 | 0.544 | 0.9944 | yes | R |
| p_r | 0.5739 | 0.0132 | 0.5607 [9.1] | 0.2535 | 0.6843 | 0.6943 | yes |  |
| p_r_std | 0.0128 | 0.0073 | 0.0055 [2.5] | 0.0072 | 0.0182 | 0.02142 | yes |  |
| p_freq_spread | 0.0128 | 0.4969 | -0.4841 [-47.5] | -0.4857 | -1.396 | 0.01343 | yes | R |
| p_vortex_density | 0.0005 | 0.333 | -0.3325 [-49.6] | -0.3314 | -0.3335 | 0.00049 | yes | R |

## qd:062  (DECOMPOSES)

`CF kC=2 fields=2 D=[0.2, 0.2] | C.CONTAG(a=0,b=1,th=3,p=0.005) ; C.EMIT(b=1,f=0,p=0.001)` · view C · flagged False · verdict KNOWN · reproduced by rules [1]

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_n_domains | 2.862 | 6.371 | -3.509 [-19.0] | -3.591 | -5.677 | 2.908 | yes | R |
| g_activity_trend | -2.917 | 0 | -2.917 [-50.0] | -2.885 | -2.986 | -2.826 | yes | R |
| g_activity_cv | 1.739 | 0 | 1.739 [50.0] | 0.9221 | 0 | 1.629 | yes |  |
| g_n_types | 1 | 2 | -1 [-16.7] | -1 | -1 | 1.333 | yes | R |
| g_corrlen | 1 | 1 | 0 [0.0] | 0 | -1 | 1 | yes |  |
| g_largest_domain | 0.9961 | 0.0425 | 0.9536 [45.8] | 0.9803 | 0.9575 | 0.996 | yes | R |

## qd:063  (DECOMPOSES)

`CA kC=2 kA=3 v=0.5 eta=0.0 | C.KURA(K=0.2) ; A.WRITECELL() ; C.MAJ(p=1.0)` · view C · flagged False · verdict KNOWN · reproduced by rules [2, 3]

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_corrlen | 9 | 1 | 8 [10.8] | 19.5 | 0 | 12 | yes |  |
| g_n_domains | 2.639 | 6.378 | -3.739 [-20.7] | -3.863 | -3.893 | 2.267 | yes | R |
| g_n_types | 2 | 2 | 0 [0.0] | 0 | -1 | 1.667 | yes |  |
| g_activity_trend | -0.8304 | 0 | -0.8304 [-14.7] | -0.8496 | -1.32 | -1.152 | yes | R |
| g_reversals | 1 | 0 | 1 [50.0] | 1 | 1 | 1 | yes | R |
| g_largest_domain | 0.8854 | 0.0538 | 0.8315 [23.2] | 0.6581 | 0.9415 | 0.7956 | yes | R |

## qd:064  (COMBINATION EFFECT)

`CAF kC=2 kA=3 v=0.3 eta=0.3 fields=1 D=[0.05] | C.CONTAG(a=1,b=0,th=4,p=1.0) ; C.EMIT(b=0,f=0,p=1.0)` · view F0 · flagged False · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| f_mean | 1001 | 501.1 | 499.6 [16.6] | 500.5 | 1499 | 1001 | yes | R |
| f_change | 10 | 5.001 | 4.999 [16.7] | 5.008 | 15 | 10 | yes | R |
| f_n_blobs | 0.6931 | 2.398 | -1.705 [-12.6] | -2.057 | -1.008 | 1.269 | no | R |
| f_moran | 0.9958 | 0.9747 | 0.0211 [0.5] | 0.018 | 0.0107 | 0.996 | yes |  |
| f_peak_sharp | 0.9412 | 0.5377 | 0.4036 [5.3] | 0.1821 | 0.2374 | 0.9224 | yes |  |
| f_std | 0.035 | 21.72 | -21.69 [-49.8] | -24.21 | -32.38 | 0.03195 | yes | R |

## qd:065  (COMBINATION EFFECT)

`CAF kC=2 kA=4 v=2.0 eta=0.0 fields=2 D=[0.02, 0.05] | A.CLIMB(f=0) ; F.AUTOCAT(f=0,g=1,p=1.0)` · view A · flagged False · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| a_step | 27.98 | 24 | 3.982 [3.8] | 0 | 10.79 | 28.19 | yes |  |
| a_density_cv | 2.002 | 0.8207 | 1.181 [18.2] | 1.184 | 1.151 | 2.025 | yes | R |
| a_type_entropy | 1.384 | 1.384 | 0 [0.0] | 0 | 0 | 1.385 | yes |  |
| a_local_align | 0.4133 | -0.0494 | 0.4627 [24.4] | 0.6058 | 0.4377 | 0.4286 | yes | R |
| a_nn_ratio | 0.4609 | 0.9979 | -0.537 [-18.4] | -0.7118 | -0.5709 | 0.3903 | yes | R |
| a_type_segregation | 0.2466 | 0.2523 | -0.0056 [-0.6] | 0.0008 | -0.0015 | 0.2544 | yes |  |

## qd:066  (ARTIFACT)

`CF kC=4 fields=2 D=[0.05, 0.02] | C.SCHELL(th=0.2) ; C.EMIT(b=3,f=1,p=1.0) ; F.DECAY(f=1,d=0.1)` · view F1 · flagged True · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| f_mean | 2.239 | 2.239 | -0 [-0.0] | -0 | -0 | 2.265 | yes |  |

## qd:067  (DECOMPOSES)

`C kC=2 | C.FLIP(a=1,b=0,p=0.5) ; C.MAJ(p=0.5) ; C.KURA(K=1.0) ; C.GAME(g=(1.5, 0.5))` · view C · flagged True · verdict KNOWN · reproduced by rules [2]

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_n_domains | 0.6931 | 0.6931 | 0 [0.0] | 2.818 | 0 | 2.291 | no |  |
| g_n_types | 1.5 | 1 | 0.5 [1.3] | 1 | 0 | 2 | yes |  |
| g_corrlen | 0 | 0 | 0 [0.0] | 3 | 0 | 1.667 | no |  |
| g_activity_trend | -4.337 | -0.5954 | -3.742 [-1.6] | -0.7433 | -6.737 | -1.387 | yes |  |
| g_reversals | 1 | 0 | 1 [50.0] | 1 | 1 | 1 | yes | R |
| g_largest_domain | 1 | 1 | 0 [0.0] | -0.0407 | 0 | 0.9529 | yes |  |

## qd:068  (DECOMPOSES)

`CF kC=4 fields=2 D=[0.05, 0.02] | C.SCHELL(th=0.2) ; C.EMIT(b=2,f=1,p=1.0) ; C.MAJ(p=0.002)` · view C · flagged True · verdict KNOWN · reproduced by rules [1, 3]

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_n_domains | 5.539 | 7.666 | -2.127 [-8.1] | -2.16 | -3.206 | 5.553 | yes | R |
| g_n_types | 4 | 4 | 0 [0.0] | 0 | 0 | 4 | yes |  |
| g_corrlen | 2 | 1 | 1 [16.7] | 1 | 3 | 2 | yes | R |
| g_activity_trend | -1.448 | 0 | -1.448 [-49.9] | -1.376 | -1.615 | -1.392 | yes | R |
| g_type_entropy | 0.9743 | 0.9998 | -0.0255 [-0.6] | -0.0271 | -0.0956 | 0.9732 | yes |  |
| g_activity_cv | 0.6579 | 0 | 0.6579 [49.8] | 0.3349 | 0.9825 | 0.694 | yes | R |

## qd:069  (DECOMPOSES)

`CF kC=4 fields=2 D=[0.05, 0.02] | C.SCHELL(th=0.2) ; C.EMIT(b=2,f=1,p=1.0) ; C.SWAP(p=0.02)` · view C · flagged True · verdict KNOWN · reproduced by rules [1, 3]

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_n_domains | 6.468 | 7.667 | -1.2 [-4.2] | -1.193 | -1.357 | 6.579 | yes | R |
| g_corrlen | 1 | 1 | 0 [0.0] | 0 | 0 | 1 | yes |  |
| g_type_entropy | 0.9998 | 0.9998 | 0 [0.0] | 0 | 0 | 0.9999 | yes |  |
| g_activity_trend | -0.512 | 0 | -0.512 [-50.2] | -0.5114 | -0.5483 | -0.515 | yes | R |

## qd:070  (DECOMPOSES)

`C kC=2 | C.GAME(g=(1.5, -0.5)) ; C.FLIP(a=1,b=0,p=0.003) ; C.MAJ(p=0.01) ; C.COPY(p=0.001)` · view C · flagged True · verdict KNOWN · reproduced by rules [1, 2]

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_corrlen | 5 | 1 | 4 [5.4] | 8.5 | 2 | 7 | yes | R |
| g_n_domains | 3.234 | 4.58 | -1.346 [-8.6] | -1.521 | 3.277 | 3.312 | yes |  |
| g_n_types | 2 | 2 | 0 [0.0] | 0 | 1 | 2 | yes |  |
| g_reversals | 1 | 0 | 1 [50.0] | 1 | 1 | 1 | yes | R |
| g_type_entropy | 0.7134 | 0.2487 | 0.4647 [21.4] | 0.568 | 0.9138 | 0.8756 | yes | R |
| g_moran | 0.7537 | 0.0003 | 0.7533 [37.7] | 0.803 | 0.5574 | 0.7783 | yes | R |

## qd:071  (COMBINATION EFFECT)

`CF kC=4 fields=2 D=[0.05, 0.02] | F.AUTOCAT(f=0,g=1,p=0.5) ; C.SCHELL(th=0.2) ; C.EMIT(b=2,f=1,p=1.0)` · view C · flagged False · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_n_domains | 7.106 | 7.665 | -0.5588 [-1.9] | -0.5392 | -0.5588 | 7.141 | yes |  |
| g_n_types | 4 | 4 | 0 [0.0] | 0 | 0 | 4 | yes |  |
| g_activity_trend | -1.503 | 0 | -1.503 [-49.9] | -1.535 | -1.503 | -1.556 | yes | R |
| g_corrlen | 1 | 1 | 0 [0.0] | 0 | 0 | 1 | yes |  |
| g_type_entropy | 0.9998 | 0.9998 | 0 [0.0] | 0 | 0 | 0.9998 | yes |  |
| g_interface | 0.6146 | 0.7487 | -0.1341 [-4.9] | -0.1302 | -0.1341 | 0.6195 | yes | R |

## qd:072  (DECOMPOSES)

`CAF kC=4 kA=4 v=0.1 eta=0.1 fields=2 D=[0.1, 0.05] | C.SCHELL(th=0.7) ; A.CLIMB(f=0) ; C.MAJ(p=0.1)` · view C · flagged True · verdict KNOWN · reproduced by rules [1, 3]

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_corrlen | 10.5 | 1 | 9.5 [25.6] | 9.5 | 10 | 12.67 | yes | R |
| g_activity_trend | -4.305 | 0 | -4.305 [-48.6] | -4.09 | -3.952 | -4.002 | yes | R |
| g_n_types | 3 | 4 | -1 [-7.1] | -1 | -2 | 2.667 | yes | R |
| g_n_domains | 1.869 | 7.661 | -5.792 [-30.4] | -5.825 | -6.418 | 1.859 | yes | R |
| g_moran | 0.92 | 0.004 | 0.916 [49.5] | 0.9341 | 0.9404 | 0.9382 | yes | R |
| g_largest_domain | 0.5668 | 0.0046 | 0.5621 [37.2] | 0.236 | 0.7499 | 0.6177 | yes |  |

## qd:073  (DECOMPOSES)

`CAF kC=4 kA=4 v=2.0 eta=0.3 fields=2 D=[0.05, 0.2] | C.SCHELL(th=0.2) ; C.KURA(K=0.001) ; A.WRITECELL()` · view C · flagged True · verdict KNOWN · reproduced by rules [1, 2]

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_n_domains | 7.319 | 7.668 | -0.3492 [-1.2] | -0.3588 | -0.3549 | 7.318 | yes |  |
| g_n_types | 4 | 4 | 0 [0.0] | 0 | 0 | 4 | yes |  |
| g_corrlen | 1 | 1 | 0 [0.0] | 0 | 0 | 1 | yes |  |
| g_type_entropy | 1 | 0.9999 | 0.0001 [0.0] | 0 | 0.0001 | 1 | yes |  |
| g_interface | 0.6447 | 0.7501 | -0.1054 [-3.8] | -0.1055 | -0.1049 | 0.645 | yes | R |
| g_reversals | 0.3466 | 0 | 0.3466 [50.2] | 0.3467 | 0.3288 | 0.3461 | yes | R |

## qd:074  (DECOMPOSES)

`CAF kC=4 kA=1 v=0.0 eta=0.1 fields=2 D=[0.2, 0.01] | C.CYCLE(th=1,p=0.05) ; A.DEPOSIT(f=1,p=0.5)` · view C · flagged False · verdict KNOWN · reproduced by rules [1]

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_n_domains | 7.385 | 7.666 | -0.2817 [-0.9] | -0.28 | -0.2708 | 7.42 | yes |  |
| g_n_types | 4 | 4 | 0 [0.0] | 0 | 0 | 4 | yes |  |
| g_corrlen | 1 | 1 | 0 [0.0] | 0 | 0 | 1 | yes |  |
| g_type_entropy | 0.9999 | 0.9999 | 0 [0.0] | 0 | 0.0001 | 1 | yes |  |
| g_interface | 0.6599 | 0.7484 | -0.0885 [-3.1] | -0.0882 | -0.0873 | 0.6618 | yes | R |
| g_cyclic_changes | 0.5012 | -0.3333 | 0.8346 [50.0] | 0.8347 | 0.5619 | 0.5015 | yes | R |

## qd:075  (ARTIFACT)

`CF kC=4 fields=2 D=[0.05, 0.02] | C.SCHELL(th=0.2) ; C.FLIP(a=1,b=2,p=0.2) ; C.EMIT(b=2,f=1,p=1.0)` · view C · flagged True · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_n_types | 3 | 3 | 0 [0.0] | 0 | 0 | 3 | yes |  |

## qd:076  (DECOMPOSES)

`C kC=2 | C.CYCLE(th=4,p=0.003) ; C.KURA(K=0.05) ; C.GAME(g=(1.9, 0.0)) ; C.GAME(g=(0.5, -0.5)) ; C.MAJ(p=0.01)` · view C · flagged False · verdict KNOWN · reproduced by rules [3, 4, 5]

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_corrlen | 0 | 1 | -1 [-50.0] | 0 | -1 | 11 | no |  |
| g_activity_trend | -4.885 | 0 | -4.885 [-50.0] | -4.746 | -0.7209 | -4.123 | yes |  |
| g_n_types | 1 | 2 | -1 [-16.7] | -1 | -1 | 1.333 | yes | R |
| g_activity_cv | 0 | 0 | 0 [0.0] | 2.202 | 0 | 1.275 | no |  |
| g_n_domains | 0.6931 | 6.384 | -5.691 [-40.2] | -6.628 | -5.691 | 1.195 | no | R |
| g_reversals | 1 | 0 | 1 [50.0] | 1 | 1 | 1 | yes | R |

## qd:077  (ARTIFACT)

`CF kC=4 fields=2 D=[0.05, 0.02] | C.SCHELL(th=0.2) ; C.CONTAG(a=0,b=1,th=7,p=0.001) ; C.EMIT(b=2,f=1,p=0.3)` · view C · flagged True · verdict KNOWN

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_n_domains | 7.116 | 7.662 | -0.5469 [-1.9] | -0.5768 | -0.5493 | 7.09 | yes |  |

## qd:078  (DECOMPOSES)

`CAF kC=4 kA=4 v=0.1 eta=0.1 fields=2 D=[0.1, 0.05] | C.CONTAG(a=3,b=1,th=7,p=1.0) ; C.SCHELL(th=0.7) ; A.CLIMB(f=0)` · view C · flagged False · verdict KNOWN · reproduced by rules [1, 2]

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_n_domains | 4.46 | 7.672 | -3.212 [-13.2] | -3.288 | -3.224 | 4.459 | yes | R |
| g_n_types | 4 | 4 | 0 [0.0] | 0 | 0 | 4 | yes |  |
| g_corrlen | 3 | 1 | 2 [25.0] | 3 | 2 | 3 | yes | R |
| g_activity_trend | -1.777 | 0 | -1.777 [-50.1] | -1.715 | -8.462 | -1.703 | yes | R |
| g_activity_cv | 1.01 | 0 | 1.01 [50.0] | 0.9779 | 0 | 0.9764 | yes |  |
| g_type_entropy | 0.8355 | 0.9998 | -0.1642 [-4.5] | -0.1622 | -0.1711 | 0.8365 | yes | R |

## qd:079  (DECOMPOSES)

`CAF kC=4 kA=4 v=0.1 eta=0.1 fields=2 D=[0.1, 0.05] | C.SCHELL(th=0.7) ; C.CYCLE(th=4,p=0.2) ; C.EMIT(b=2,f=0,p=0.2)` · view C · flagged True · verdict KNOWN · reproduced by rules [1, 2]

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_activity_trend | -7.416 | 0 | -7.416 [-50.0] | -7.386 | -6.114 | -7.286 | yes | R |

## qd:080  (DECOMPOSES)

`CA kC=3 kA=1 v=2.0 eta=0.0 | A.WRITECELL() ; A.QUORUM(th=4,r=2.0) ; C.CYCLE(th=2,p=0.2) ; C.CYCLE(th=2,p=1.0)` · view C · flagged True · verdict KNOWN · reproduced by rules [1, 3, 4]

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_n_domains | 5.724 | 7.336 | -1.612 [-6.2] | -1.745 | -1.579 | 5.773 | yes | R |
| g_n_types | 3 | 3 | 0 [0.0] | 0 | 0 | 3 | yes |  |
| g_corrlen | 2 | 1 | 1 [16.7] | 1 | 1 | 2 | yes | R |
| g_type_entropy | 0.9999 | 0.9998 | 0.0001 [0.0] | -0 | 0.0002 | 0.9999 | yes |  |
| g_activity | 0.7437 | 0 | 0.7437 [49.9] | 0.7108 | 0.7407 | 0.7444 | yes | R |
| g_moran | 0.5441 | 0.0045 | 0.5395 [38.5] | 0.5543 | 0.5319 | 0.5417 | yes | R |

## qd:081  (DECOMPOSES)

`CA kC=4 kA=4 v=0.0 eta=0.1 | C.KURA(K=1.0) ; A.ATTRACT(b=1,r=3.0) ; C.SCHELL(th=0.4) ; C.CONTAG(a=2,b=3,th=1,p=1.0) ; A.WRITECELL()` · view C · flagged False · verdict KNOWN-VARIANT · reproduced by rules [3, 4, 5]

| feature | real | knock-out | e [e/s] | 4x | 3x | card | rep | R |
|---|---|---|---|---|---|---|---|---|
| g_n_domains | 5.822 | 7.668 | -1.846 [-6.8] | -1.843 | -1.861 | 5.825 | yes | R |
| g_n_types | 4 | 4 | 0 [0.0] | 0 | 0 | 4 | yes |  |
| g_corrlen | 3 | 1 | 2 [25.0] | 2 | 2 | 2.667 | yes | R |
| g_largest_domain | 0.6907 | 0.0044 | 0.6863 [49.4] | 0.7005 | 0.6874 | 0.7159 | yes | R |
| g_type_entropy | 0.5862 | 0.9998 | -0.4136 [-13.0] | -0.4244 | -0.414 | 0.566 | yes | R |
| g_moran | 0.5395 | 0.0022 | 0.5373 [49.8] | 0.5149 | 0.531 | 0.5164 | yes | R |

