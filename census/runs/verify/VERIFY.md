# Simulation verification of literature verdicts (VERIFICATION_PLAN.md)

214 phenomena classified; Stage 3 candidates: 84

| batch | class | n |
|---|---|---|
| census | ARTIFACT | 14 |
| census | COMBINATION EFFECT | 12 |
| census | DECOMPOSES | 15 |
| census | SINGLE MECHANISM | 19 |
| qd | ARTIFACT | 13 |
| qd | COMBINATION EFFECT | 27 |
| qd | DECOMPOSES | 39 |
| qd | SINGLE MECHANISM | 2 |
| review | ARTIFACT | 10 |
| review | COMBINATION EFFECT | 28 |
| review | DECOMPOSES | 12 |
| review | SINGLE MECHANISM | 23 |

| phenomenon | rep | flagged | class | robust features | reproduced by rules |
|---|---|---|---|---|---|
| census:g1:lattice-phase-1 | `C kC=1 | C.KURA(K=0.1)` | False | SINGLE MECHANISM | p_local_r, p_freq_spread, p_vortex_density |  |
| census:g1:network-phase-1 | `N kN=1 | N.KURA(K=0.1)` | False | SINGLE MECHANISM | p_freq_spread, p_r, p_r_std |  |
| census:g1:network-1 | `N kN=2 | N.GAME(g=(1.9, 0.0))` | False | SINGLE MECHANISM | n_type_modularity |  |
| census:g1:grid-1 | `C kC=2 | C.COPY(p=1.0)` | False | SINGLE MECHANISM | g_reversals, g_largest_domain |  |
| census:g2:lattice-phase-2 | `C kC=1 | C[any] ALWAYS() -> PHASE_COUPLE(K=0.1) p=0.1` | False | SINGLE MECHANISM | p_local_r, p_vortex_density, p_freq_spread |  |
| census:g2:network-phase-1 | `N kN=1 | N[any] ALWAYS() -> PHASE_COUPLE(K=0.1) p=1.0` | False | SINGLE MECHANISM | p_r, p_freq_spread, p_r_std |  |
| census:g1:grid-2 | `C kC=2 | C.GAME(g=(1.5, -0.5))` | False | SINGLE MECHANISM | g_corrlen, g_n_domains, g_n_types, g_largest_domain, g_moran |  |
| census:g1:network-6 | `N kN=3 | N.CYCLE(th=3,p=0.1)` | False | SINGLE MECHANISM | n_cyclic_changes |  |
| census:g1:network-15 | `N kN=3 | N.COPY(p=0.05)` | False | SINGLE MECHANISM | n_oscillation, n_consensus_gain |  |
| census:g1:grid-4 | `C kC=4 | C.CYCLE(th=1,p=1.0)` | False | SINGLE MECHANISM | g_activity |  |
| census:g1:network-12 | `N kN=2 | N.REWIRE(p=0.01)` | False | SINGLE MECHANISM | n_type_modularity, n_edge_turnover |  |
| census:g1:agents-1 | `A kA=1 v=1.0 eta=0.3 | A.ATTRACT(b=any,r=1.0)` | True | SINGLE MECHANISM | a_density_cv, a_step, a_nn_ratio, a_local_align |  |
| census:g1:agents-2 | `A kA=1 v=1.0 eta=0.3 | A.REPEL(b=any,r=3.0)` | True | SINGLE MECHANISM | a_step |  |
| census:g2:network-1 | `N kN=2 | N[any] ALWAYS() -> COPY_RAND() p=0.1` | False | SINGLE MECHANISM | n_oscillation |  |
| census:g2:grid-1 | `C kC=2 | C[any] ALWAYS() -> COPY_RAND() p=1.0` | False | SINGLE MECHANISM | g_reversals, g_largest_domain |  |
| census:g1:agents-3 | `A kA=1 v=1.0 eta=3.14159 | A.ALIGN(r=1.0)` | True | SINGLE MECHANISM | a_density_cv, a_local_align |  |
| census:g1:grid-6 | `C kC=2 | C.CONTAG(a=0,b=1,th=1,p=0.1)` | True | SINGLE MECHANISM | g_largest_domain, g_n_domains |  |
| census:g1:agents-8 | `AF kA=1 v=1.0 eta=0.3 fields=1 D=[0.05] | A.CLIMB(f=0)` | True | SINGLE MECHANISM | a_step |  |
| census:g1:grid-11 | `C kC=2 | C.COPY(p=1.0) ; C.MAJ(p=1.0)` | True | ARTIFACT (flagged view did not recur on fresh seeds) |  |  |
| census:g1:network-phase-2 | `N kN=3 | N.KURA(K=1.0) ; N.REWIRE(p=1.0)` | True | COMBINATION EFFECT | p_r_std, p_freq_spread |  |
| census:g1:grid-17 | `C kC=4 | C.MAJ(p=1.0) ; C.MAJ(p=0.1)` | True | ARTIFACT (flagged view did not recur on fresh seeds) |  |  |
| census:g2:agents-1 | `A kA=1 v=1.0 eta=0.3 | A[any] ALWAYS() -> HEAD_TO(b=any,r=1.0) p=1.0` | False | SINGLE MECHANISM | a_density_cv, a_step, a_nn_ratio, a_local_align |  |
| census:novelty:network-2 | `N kN=2 | N.GAME(g=(1.9, 0.0)) ; N.MAJ(p=0.1)` | False | COMBINATION EFFECT | n_activity, n_type_modularity |  |
| census:novelty:grid-1 | `C kC=2 | C.GAME(g=(1.5, 0.5)) ; C.KURA(K=1.0)` | False | DECOMPOSES | g_n_domains, g_corrlen, g_reversals, g_largest_domain | [0] |
| census:novelty:network-phase-1 | `N kN=2 | N.GAME(g=(1.9, 0.0)) ; N.KURA(K=0.1)` | False | DECOMPOSES | p_freq_spread | [1] |
| census:g1:network-7 | `N kN=3 | N.COPY(p=0.1) ; N.REWIRE(p=0.03)` | True | ARTIFACT |  |  |
| census:novelty:grid-3 | `CA kC=3 kA=3 v=0.0 eta=3.14159 | A.WRITECELL()` | True | ARTIFACT |  |  |
| census:g1:field-1 | `F fields=2 D=[0.01, 0.2] | F.AUTOCAT(f=0,g=1,p=0.05)` | True | ARTIFACT |  |  |
| census:g1:network-14 | `N kN=2 | N.REWIRE(p=0.1) ; N.CYCLE(th=1,p=0.1)` | True | ARTIFACT |  |  |
| census:g1:grid-25 | `C kC=4 | C.SWAP(p=1.0) ; C.CYCLE(th=3,p=0.1)` | True | COMBINATION EFFECT | g_n_domains, g_reversals, g_activity_cv |  |
| census:g1:network-10 | `N kN=4 | N.REWIRE(p=0.01) ; N.MAJ(p=1.0)` | False | DECOMPOSES | n_type_entropy, n_consensus_gain | [1] |
| census:g2:field-1 | `F fields=2 D=[0.01, 0.2] | F[any] ALWAYS() -> AUTOCAT(f=0,g=1,p=0.05) p=1.0` | False | ARTIFACT |  |  |
| census:novelty:network-3 | `N kN=2 | N.REWIRE(p=0.003) ; N.GAME(g=(1.9, 0.0))` | True | DECOMPOSES | n_components | [0] |
| census:g1:grid-19 | `C kC=2 | C.MAJ(p=0.003) ; C.GAME(g=(1.5, -0.5))` | True | COMBINATION EFFECT | g_corrlen, g_activity_cv, g_largest_domain |  |
| census:novelty:grid-4 | `C kC=2 | C.GAME(g=(1.3, 0.0)) ; C.COPY(p=0.03)` | True | COMBINATION EFFECT | g_n_types, g_largest_domain |  |
| census:novelty:lattice-phase-1 | `C kC=2 | C.KURA(K=0.01) ; C.GAME(g=(0.5, -0.5))` | False | ARTIFACT |  |  |
| census:g1:field-2 | `F fields=2 D=[0.2, 0.05] | F.FEED(f=0,F=1.0) ; F.AUTOCAT(f=0,g=1,p=1.0)` | True | ARTIFACT |  |  |
| census:g1:network-13 | `N kN=2 | N.KURA(K=0.03) ; N.REWIRE(p=0.01)` | False | DECOMPOSES | n_type_modularity, n_edge_turnover | [1] |
| census:g1:grid-29 | `C kC=4 | C.CYCLE(th=3,p=0.3) ; C.MAJ(p=1.0)` | True | ARTIFACT (flagged view did not recur on fresh seeds) |  |  |
| census:g1:grid-26 | `C kC=4 | C.SAND(p=0.002) ; C.SCHELL(th=0.5)` | True | ARTIFACT |  |  |
| census:novelty:agents-6 | `AF kA=1 v=1.0 eta=1.0 fields=1 D=[0.05] | F.FEED(f=0,F=0.1) ; A.CLIMB(f=0)` | False | ARTIFACT |  |  |
| census:novelty:agents-17 | `AF kA=1 v=1.0 eta=1.0 fields=1 D=[0.05] | A.CLIMB(f=0) ; A.REPEL(b=any,r=2.0)` | False | COMBINATION EFFECT | a_step, a_density_cv, a_nn_ratio, a_local_align |  |
| census:novelty:network-5 | `N kN=2 | N.KURA(K=0.1) ; N.COPY(p=1.0) ; N.GAME(g=(1.5, 0.5))` | False | DECOMPOSES | n_activity, n_consensus_gain | [1, 2] |
| census:novelty:agents-8 | `AF kA=1 v=1.0 eta=1.0 fields=1 D=[0.05] | A.CLIMB(f=0) ; A.ALIGN(r=2.0)` | True | COMBINATION EFFECT | a_density_cv, a_step, a_local_align |  |
| census:novelty:network-4 | `N kN=2 | N.COPY(p=1.0) ; N.REWIRE(p=0.3) ; N.GAME(g=(1.9, 0.0))` | True | DECOMPOSES | n_components, n_edge_turnover, n_oscillation | [0, 1] |
| census:novelty:field-1 | `AF kA=1 v=1.0 eta=1.0 fields=1 D=[0.02] | A.ALIGN(r=8.0) ; A.DEPOSIT(f=0,p=0.1)` | False | COMBINATION EFFECT | f_std, f_n_blobs, f_peak_sharp |  |
| census:novelty:agents-14 | `AF kA=2 v=0.5 eta=1.0 fields=1 D=[0.05] | A.ALIGN(r=1.0) ; A.DEPOSIT(f=0,p=1.0)` | True | DECOMPOSES | a_density_cv, a_nn_ratio, a_local_align | [0] |
| census:novelty:field-2 | `AF kA=1 v=1.0 eta=1.0 fields=1 D=[0.05] | A.ATTRACT(b=any,r=0.5) ; A.DEPOSIT(f=0,p=0.2)` | False | ARTIFACT |  |  |
| census:novelty:lattice-phase-2 | `C kC=2 | C.CYCLE(th=4,p=0.01) ; C.GAME(g=(1.5, 0.5)) ; C.KURA(K=0.1)` | False | DECOMPOSES | p_local_r, p_freq_spread, p_vortex_density, p_r | [2] |
| census:novelty:agents-16 | `CAF kC=1 kA=4 v=1.0 eta=0.3 fields=2 D=[0.02, 0.02] | A.CLIMB(f=1) ; A.TURNCELL()` | True | COMBINATION EFFECT | a_density_cv |  |
| census:novelty:grid-7 | `CAF kC=2 kA=1 v=1.0 eta=1.0 fields=2 D=[0.02, 0.1] | C.MAJ(p=0.1) ; A.WRITECELL()` | True | COMBINATION EFFECT | g_reversals, g_interface |  |
| census:novelty:grid-8 | `CAF kC=2 kA=1 v=1.0 eta=1.0 fields=2 D=[0.02, 0.1] | A.WRITECELL() ; A.CLIMB(f=0)` | True | COMBINATION EFFECT | g_activity_trend, g_reversals |  |
| census:novelty:grid-9 | `C kC=2 | C.CYCLE(th=4,p=0.01) ; C.GAME(g=(1.5, -0.5)) ; C.KURA(K=0.0003)` | False | DECOMPOSES | g_n_domains, g_activity_trend, g_reversals, g_moran | [0, 1] |
| census:novelty:field-3 | `CAF kC=2 kA=1 v=1.0 eta=1.0 fields=2 D=[0.02, 0.1] | C.MAJ(p=0.1) ; C.EMIT(b=1,f=0,p=1.0)` | True | COMBINATION EFFECT | f_std, f_n_blobs |  |
| census:novelty:grid-10 | `CAF kC=2 kA=1 v=1.0 eta=1.0 fields=2 D=[0.02, 0.1] | C.MAJ(p=0.0003) ; A.CLIMB(f=0)` | True | DECOMPOSES | g_activity_cv, g_interface | [0] |
| census:novelty:grid-12 | `CAF kC=2 kA=1 v=1.0 eta=1.0 fields=2 D=[0.02, 0.1] | C.COPY(p=0.001) ; A.CLIMB(f=0)` | False | DECOMPOSES | g_reversals | [0] |
| census:novelty:grid-15 | `C kC=2 | C.CONTAG(a=1,b=0,th=6,p=0.005) ; C.GAME(g=(1.9, 0.0)) ; C.KURA(K=0.0003)` | False | DECOMPOSES | g_n_domains, g_corrlen, g_reversals, g_largest_domain, g_moran | [1] |
| census:novelty:grid-11 | `CAF kC=2 kA=1 v=1.0 eta=1.0 fields=2 D=[0.02, 0.1] | C.CYCLE(th=7,p=0.03) ; A.CLIMB(f=0)` | True | ARTIFACT |  |  |
| census:novelty:grid-13 | `CAF kC=4 kA=1 v=0.0 eta=2.0 fields=2 D=[0.02, 0.1] | C.COPY(p=1.0) ; C.EMIT(b=3,f=0,p=0.3)` | True | DECOMPOSES | g_n_domains, g_interface, g_reversals | [0] |
| census:novelty:field-4 | `CAF kC=4 kA=1 v=0.0 eta=2.0 fields=2 D=[0.02, 0.1] | A.ALIGN(r=8.0) ; C.MAJ(p=0.02) ; C.SWAP(p=0.1) ; C.EMIT(b=3,f=0,p=0.3)` | False | DECOMPOSES | f_std, f_n_blobs | [1, 3] |
| review:001 | `N kN=4 | N.COPY(p=1.0)` | True | ARTIFACT |  |  |
| review:002 | `C kC=1 | C[any] ALWAYS() -> PHASE_COUPLE(K=1.0) p=0.1` | False | SINGLE MECHANISM | p_local_r, p_r, p_freq_spread, p_vortex_density |  |
| review:003 | `C kC=3 | C.MAJ(p=1.0)` | False | SINGLE MECHANISM | g_n_domains, g_activity_trend, g_reversals |  |
| review:004 | `N kN=3 | N.REWIRE(p=1.0)` | False | SINGLE MECHANISM | n_type_modularity, n_edge_turnover, n_degree_cv |  |
| review:005 | `C kC=3 | C.SCHELL(th=0.2)` | False | ARTIFACT |  |  |
| review:006 | `N kN=1 | N.KURA(K=0.2)` | False | SINGLE MECHANISM | p_r, p_freq_spread |  |
| review:007 | `C kC=1 | C.KURA(K=1.0) ; C.KURA(K=1.0)` | True | DECOMPOSES | p_local_r, p_freq_spread, p_vortex_density | [0] |
| review:008 | `C kC=3 | C.CYCLE(th=3,p=0.1)` | True | SINGLE MECHANISM | g_corrlen, g_n_domains, g_moran, g_cyclic_changes |  |
| review:009 | `N kN=1 | N[any] ALWAYS() -> PHASE_COUPLE(K=0.2) p=0.1` | True | SINGLE MECHANISM | p_r, p_freq_spread |  |
| review:010 | `C kC=2 | C.COPY(p=0.005)` | False | SINGLE MECHANISM | g_reversals, g_interface |  |
| review:011 | `N kN=2 | N[0] ALWAYS() -> REWIRE_SAME() p=1.0` | True | SINGLE MECHANISM | n_components |  |
| review:012 | `C kC=4 | C.CYCLE(th=1,p=0.3)` | True | SINGLE MECHANISM | g_interface, g_activity |  |
| review:013 | `A kA=1 v=1.0 eta=1.0 | A.ALIGN(r=0.5)` | False | SINGLE MECHANISM | a_density_cv, a_polar, a_nn_ratio, a_local_align |  |
| review:014 | `N kN=2 | N.CONTAG(a=0,b=1,th=3,p=0.1)` | False | SINGLE MECHANISM | n_type_entropy, n_consensus_gain |  |
| review:015 | `N kN=3 | N.CYCLE(th=2,p=0.1)` | True | SINGLE MECHANISM | n_cyclic_changes |  |
| review:016 | `CF kC=2 fields=1 D=[0.2] | C.GAME(g=(1.9, 0.0))` | True | SINGLE MECHANISM | g_n_domains, g_corrlen, g_reversals, g_largest_domain |  |
| review:017 | `C kC=2 | C[any] ALWAYS() -> IMITATE_BEST(g=(1.5, 0.5)) p=0.1` | True | SINGLE MECHANISM | g_reversals |  |
| review:018 | `N kN=3 | N.COPY(p=1.0) ; N.REWIRE(p=1.0)` | False | COMBINATION EFFECT | n_edge_turnover, n_type_modularity, n_giant |  |
| review:019 | `N kN=3 | N.COPY(p=1.0) ; N.REWIRE(p=0.1)` | False | COMBINATION EFFECT | n_edge_turnover, n_consensus_gain |  |
| review:020 | `N kN=4 | N.COPY(p=0.1) ; N.KURA(K=0.1)` | False | DECOMPOSES | n_oscillation, n_consensus_gain | [0] |
| review:021 | `C kC=2 | C.COPY(p=1.0) ; C.GAME(g=(1.5, 0.5))` | True | COMBINATION EFFECT | g_n_domains, g_corrlen, g_reversals, g_activity |  |
| review:022 | `C kC=2 | C.MAJ(p=1.0) ; C.SWAP(p=1.0)` | True | DECOMPOSES | g_n_domains, g_reversals | [0] |
| review:023 | `C kC=2 | C.GAME(g=(1.5, 0.5)) ; C.COPY(p=0.1)` | True | COMBINATION EFFECT | g_n_domains, g_corrlen, g_reversals, g_largest_domain, g_activity_trend |  |
| review:024 | `C kC=2 | C.GAME(g=(1.5, 0.5)) ; C.SWAP(p=0.1)` | True | DECOMPOSES | g_reversals | [0] |
| review:025 | `C kC=3 | C.MAJ(p=1.0) ; C.SCHELL(th=0.3)` | True | COMBINATION EFFECT | g_n_types |  |
| review:026 | `C kC=3 | C.SWAP(p=1.0) ; C.SCHELL(th=0.5)` | True | COMBINATION EFFECT | g_n_domains, g_reversals, g_activity_cv |  |
| review:027 | `C kC=3 | C.SWAP(p=0.1) ; C.SCHELL(th=0.5)` | True | COMBINATION EFFECT | g_n_domains, g_reversals, g_moran |  |
| review:028 | `C kC=3 | C.SCHELL(th=0.5) ; C.COPY(p=1.0)` | True | COMBINATION EFFECT | g_n_domains, g_reversals, g_spec_peak |  |
| review:029 | `N kN=3 | N.MAJ(p=1.0) ; N.REWIRE(p=0.1)` | True | ARTIFACT |  |  |
| review:030 | `C kC=3 | C[any] NEXT_GE(th=1) -> SWAP_RAND() p=1.0` | True | ARTIFACT |  |  |
| review:031 | `C kC=3 | C[any] SAME_LT(th=0.5) -> MOVE_EMPTY() p=1.0` | False | SINGLE MECHANISM | g_activity_trend, g_corrlen, g_n_domains, g_moran |  |
| review:032 | `C kC=3 | C[1] SAME_LT(th=0.3) -> MOVE_EMPTY() p=0.1` | False | SINGLE MECHANISM | g_activity_trend, g_reversals |  |
| review:033 | `C kC=4 | C[any] NEXT_GE(th=1) -> SET_NEXT() p=0.1` | False | SINGLE MECHANISM | g_interface |  |
| review:034 | `C kC=4 | C[any] NEXT_GE(th=1) -> COPY_RAND() p=1.0` | True | SINGLE MECHANISM | g_n_domains, g_corrlen, g_activity_trend, g_moran |  |
| review:035 | `N kN=3 | N[any] NEXT_GE(th=1) -> COPY_RAND() p=0.1` | True | ARTIFACT |  |  |
| review:036 | `N kN=4 | N[any] ALWAYS() -> COPY_RAND() p=0.2` | True | SINGLE MECHANISM | n_oscillation |  |
| review:037 | `C kC=2 | C.MAJ(p=0.1) ; C.GAME(g=(1.5, -0.5))` | True | COMBINATION EFFECT | g_n_domains, g_corrlen, g_reversals, g_type_entropy |  |
| review:038 | `C kC=4 | C.COPY(p=1.0) ; C.SCHELL(th=0.7)` | True | DECOMPOSES | g_spec_peak | [0] |
| review:039 | `N kN=4 | N.CYCLE(th=3,p=1.0) ; N.COPY(p=0.1)` | False | DECOMPOSES | n_oscillation, n_consensus_gain, n_type_entropy | [1] |
| review:040 | `N kN=2 | N.REWIRE(p=0.01) ; N.GAME(g=(1.9, 0.0))` | True | ARTIFACT |  |  |
| review:041 | `C kC=2 | C[any] CNT_GE(b=0,th=3) -> IMITATE_BEST(g=(1.5, 0.5)) p=1.0` | True | SINGLE MECHANISM | g_n_domains, g_corrlen, g_reversals, g_largest_domain, g_moran |  |
| review:042 | `C kC=2 | C[any] CNT_LE(b=0,th=3) -> IMITATE_BEST(g=(1.9, 0.0)) p=0.1` | True | SINGLE MECHANISM | g_n_domains, g_corrlen, g_activity_trend, g_reversals |  |
| review:043 | `A kA=1 v=0.3 eta=0.1 | A[any] ALWAYS() -> HEAD_MEAN(b=any,r=0.5) p=1.0` | False | SINGLE MECHANISM | a_density_cv, a_local_align, a_polar, a_nn_ratio |  |
| review:044 | `C kC=2 | C.SWAP(p=1.0) ; C.CYCLE(th=3,p=0.1)` | True | DECOMPOSES | g_reversals | [0] |
| review:045 | `C kC=3 | C.COPY(p=1.0) ; C.CYCLE(th=3,p=1.0)` | True | COMBINATION EFFECT | g_corrlen, g_n_domains, g_moran, g_reversals |  |
| review:046 | `C kC=3 | C.MAJ(p=1.0) ; C.CYCLE(th=1,p=1.0)` | True | COMBINATION EFFECT | g_n_domains, g_corrlen, g_reversals, g_moran |  |
| review:047 | `C kC=3 | C.CYCLE(th=1,p=1.0) ; C.SCHELL(th=0.5)` | True | COMBINATION EFFECT | g_n_domains, g_largest_domain, g_reversals |  |
| review:048 | `C kC=3 | C.CYCLE(th=3,p=0.1) ; C.MAJ(p=1.0)` | True | DECOMPOSES | g_largest_domain, g_reversals, g_n_domains | [1] |
| review:049 | `C kC=3 | C.CYCLE(th=3,p=0.1) ; C.SCHELL(th=0.3)` | True | COMBINATION EFFECT | g_n_domains, g_corrlen, g_moran |  |
| review:050 | `C kC=4 | C.CYCLE(th=1,p=1.0) ; C.SCHELL(th=0.3)` | True | COMBINATION EFFECT | g_n_domains, g_activity, g_moran |  |
| review:051 | `C kC=4 | C.SCHELL(th=0.5) ; C.CYCLE(th=1,p=1.0)` | True | COMBINATION EFFECT | g_n_domains, g_activity, g_interface |  |
| review:052 | `N kN=4 | N.MAJ(p=1.0) ; N.REWIRE(p=0.03)` | False | COMBINATION EFFECT | n_consensus_gain, n_type_entropy, n_edge_turnover |  |
| review:053 | `N kN=4 | N.CYCLE(th=1,p=1.0) ; N.REWIRE(p=1.0)` | True | COMBINATION EFFECT | n_edge_turnover, n_oscillation |  |
| review:054 | `C kC=2 | C.MAJ(p=0.003) ; C.GAME(g=(1.3, 0.0))` | True | COMBINATION EFFECT | g_n_domains, g_activity_trend, g_corrlen, g_reversals |  |
| review:055 | `C kC=2 | C.CYCLE(th=7,p=0.1) ; C.GAME(g=(1.5, 0.5))` | True | DECOMPOSES | g_n_domains, g_corrlen, g_reversals, g_largest_domain | [1] |
| review:056 | `C kC=2 | C.GAME(g=(1.5, 0.5)) ; C.SWAP(p=0.05)` | True | COMBINATION EFFECT | g_n_domains, g_corrlen, g_reversals, g_largest_domain, g_type_entropy |  |
| review:057 | `C kC=3 | C.FLIP(a=1,b=0,p=0.1) ; C.SCHELL(th=0.5)` | True | COMBINATION EFFECT | g_n_domains, g_reversals |  |
| review:058 | `C kC=3 | C.FLIP(a=1,b=0,p=0.1) ; C.SCHELL(th=0.3)` | True | COMBINATION EFFECT | g_activity_trend, g_reversals |  |
| review:059 | `N kN=2 | N.MAJ(p=1.0) ; N.REWIRE(p=0.5)` | True | ARTIFACT |  |  |
| review:060 | `CA kC=4 kA=3 v=0.3 eta=0.3 | C.CYCLE(th=3,p=1.0)` | True | SINGLE MECHANISM | g_corrlen, g_n_domains, g_activity_trend, g_moran |  |
| review:061 | `C kC=2 | C.FLIP(a=1,b=0,p=0.03) ; C.GAME(g=(1.9, 0.0))` | True | DECOMPOSES | g_n_domains, g_corrlen, g_n_types, g_activity_trend, g_reversals, g_largest_domain | [1] |
| review:062 | `C kC=3 | C.CYCLE(th=3,p=0.3) ; C.SWAP(p=1.0)` | True | COMBINATION EFFECT | g_corrlen, g_n_domains, g_moran, g_cyclic_changes |  |
| review:063 | `C kC=4 | C.MAJ(p=0.005) ; C.SCHELL(th=0.7)` | True | COMBINATION EFFECT | g_n_domains, g_type_entropy, g_spec_peak |  |
| review:064 | `C kC=4 | C.CYCLE(th=2,p=1.0) ; C.SWAP(p=1.0)` | True | COMBINATION EFFECT | g_n_domains, g_corrlen, g_moran, g_activity |  |
| review:065 | `C kC=4 | C.SCHELL(th=0.5) ; C.SAND(p=0.001)` | True | ARTIFACT |  |  |
| review:066 | `C kC=2 | C.CYCLE(th=7,p=0.03) ; C.GAME(g=(1.5, 0.5))` | True | COMBINATION EFFECT | g_n_domains, g_corrlen, g_activity_trend, g_reversals, g_largest_domain |  |
| review:067 | `C kC=2 | C.FLIP(a=1,b=0,p=0.01) ; C.GAME(g=(1.3, 0.0))` | True | COMBINATION EFFECT | g_n_domains, g_corrlen, g_n_types, g_reversals, g_activity_trend |  |
| review:068 | `C kC=3 | C.CONTAG(a=0,b=1,th=3,p=1.0) ; C.SCHELL(th=0.3)` | True | COMBINATION EFFECT | g_n_types |  |
| review:069 | `C kC=3 | C.CYCLE(th=1,p=0.1) ; C.FLIP(a=0,b=1,p=0.1)` | True | COMBINATION EFFECT | g_n_types, g_type_entropy, g_interface, g_reversals |  |
| review:070 | `C kC=4 | C.CYCLE(th=3,p=1.0) ; C.KURA(K=0.3)` | True | DECOMPOSES | g_corrlen, g_n_domains, g_activity_trend, g_moran | [0] |
| review:071 | `N kN=3 | N.REWIRE(p=0.2) ; N.CYCLE(th=1,p=0.1)` | True | ARTIFACT |  |  |
| review:072 | `N kN=2 | N.COPY(p=1.0) ; N.GAME(g=(1.9, 0.0)) ; N.REWIRE(p=1.0)` | True | ARTIFACT |  |  |
| review:073 | `C kC=2 | C.SWAP(p=0.3) ; C.CYCLE(th=4,p=0.01) ; C.GAME(g=(1.3, 0.0)) ; C.KURA(K=1.0)` | True | DECOMPOSES | g_n_domains, g_corrlen, g_reversals, g_largest_domain, g_moran | [2] |
| qd:001 | `AF kA=1 v=0.3 eta=1.0 fields=1 D=[0.05] | A.CLIMB(f=0)` | False | SINGLE MECHANISM | a_density_cv, a_step |  |
| qd:002 | `CA kC=4 kA=1 v=0.0 eta=0.3 | A.WRITECELL()` | False | SINGLE MECHANISM | g_reversals |  |
| qd:003 | `C kC=2 | C.GAME(g=(1.5, 0.5)) ; C.GAME(g=(1.5, 0.5))` | False | DECOMPOSES | g_n_domains, g_corrlen, g_largest_domain | [0] |
| qd:004 | `C kC=2 | C.GAME(g=(1.5, 0.5)) ; C.MAJ(p=0.1)` | True | COMBINATION EFFECT | g_activity_trend, g_activity_cv, g_corrlen, g_n_types, g_reversals, g_largest_domain |  |
| qd:005 | `C kC=2 | C.GAME(g=(1.5, 0.5)) ; C.COPY(p=0.1)` | True | COMBINATION EFFECT | g_n_domains, g_corrlen, g_reversals, g_largest_domain, g_activity_trend |  |
| qd:006 | `C kC=2 | C.GAME(g=(1.9, 0.0)) ; C.KURA(K=0.1)` | False | DECOMPOSES | g_reversals, g_largest_domain | [0] |
| qd:007 | `C kC=2 | C.GAME(g=(1.9, 0.0)) ; C.KURA(K=0.1)` | True | DECOMPOSES | p_local_r | [1] |
| qd:008 | `N kN=2 | N.GAME(g=(1.5, 0.5)) ; N.GAME(g=(1.9, 0.0))` | True | DECOMPOSES | n_type_modularity | [1] |
| qd:009 | `N kN=2 | N.GAME(g=(1.5, 0.5)) ; N.KURA(K=1.0)` | True | ARTIFACT |  |  |
| qd:010 | `N kN=2 | N.GAME(g=(1.9, 0.0)) ; N.KURA(K=0.1)` | True | DECOMPOSES | n_consensus_gain, n_oscillation | [0] |
| qd:011 | `N kN=2 | N.GAME(g=(1.9, 0.0)) ; N.KURA(K=0.1)` | False | DECOMPOSES | p_freq_spread | [1] |
| qd:012 | `C kC=2 | C.GAME(g=(1.5, 0.5)) ; C.SWAP(p=1.0)` | False | COMBINATION EFFECT | g_n_domains, g_reversals, g_activity |  |
| qd:013 | `C kC=2 | C.GAME(g=(1.5, -0.5)) ; C.MAJ(p=0.1)` | True | COMBINATION EFFECT | g_n_domains, g_corrlen, g_reversals, g_type_entropy |  |
| qd:014 | `N kN=2 | N.GAME(g=(1.9, 0.0)) ; N.REWIRE(p=0.01)` | True | ARTIFACT |  |  |
| qd:015 | `N kN=2 | N.GAME(g=(1.9, 0.0)) ; N.REWIRE(p=0.003)` | True | COMBINATION EFFECT | n_oscillation, n_edge_turnover |  |
| qd:016 | `C kC=2 | C.GAME(g=(1.5, 0.5)) ; C.MAJ(p=0.01)` | True | COMBINATION EFFECT | g_activity_trend, g_largest_domain, g_n_types, g_reversals, g_n_domains, g_spec_peak |  |
| qd:017 | `C kC=2 | C.GAME(g=(1.5, 0.5)) ; C.SWAP(p=0.03)` | True | DECOMPOSES | g_n_domains, g_corrlen, g_reversals, g_largest_domain, g_moran | [0] |
| qd:018 | `C kC=2 | C.GAME(g=(1.5, 0.5)) ; C.KURA(K=0.2)` | True | DECOMPOSES | p_local_r | [1] |
| qd:019 | `N kN=2 | N.GAME(g=(1.5, 0.5)) ; N.KURA(K=0.2)` | False | DECOMPOSES | p_r, p_freq_spread, p_r_std | [1] |
| qd:020 | `C kC=2 | C.GAME(g=(1.5, 0.5)) ; C.CYCLE(th=7,p=0.1)` | True | DECOMPOSES | g_n_domains, g_corrlen, g_reversals, g_largest_domain | [0] |
| qd:021 | `C kC=2 | C.GAME(g=(1.5, 0.5)) ; C.CYCLE(th=3,p=0.03)` | True | ARTIFACT |  |  |
| qd:022 | `C kC=2 | C.GAME(g=(1.5, -0.5)) ; C.SWAP(p=0.002)` | False | ARTIFACT (flagged view did not recur on fresh seeds) |  |  |
| qd:023 | `C kC=2 | C.GAME(g=(1.5, -0.5)) ; C.COPY(p=0.002)` | False | COMBINATION EFFECT | g_n_types, g_n_domains, g_reversals, g_largest_domain |  |
| qd:024 | `C kC=2 | C.SWAP(p=0.001) ; C.GAME(g=(1.5, -0.5))` | True | DECOMPOSES | g_n_domains, g_n_types, g_largest_domain | [1] |
| qd:025 | `C kC=2 | C.GAME(g=(1.5, 0.5)) ; C.FLIP(a=0,b=1,p=0.0003)` | True | DECOMPOSES | g_n_domains, g_corrlen, g_reversals | [0] |
| qd:026 | `C kC=2 | C.GAME(g=(1.5, 0.5)) ; C.CYCLE(th=7,p=0.05)` | True | DECOMPOSES | g_reversals | [0] |
| qd:027 | `AF kA=1 v=0.3 eta=0.3 fields=1 D=[0.02] | A.CLIMB(f=0) ; A.ALIGN(r=1.0)` | True | COMBINATION EFFECT | a_density_cv, a_local_align, a_step, a_polar_std |  |
| qd:028 | `CA kC=1 kA=3 v=0.3 eta=0.0 | A.ALIGN(r=1.0) ; A.WRITECELL()` | True | ARTIFACT |  |  |
| qd:029 | `CA kC=3 kA=1 v=0.3 eta=3.14159 | C.COPY(p=1.0) ; A.WRITECELL()` | True | ARTIFACT |  |  |
| qd:030 | `CA kC=3 kA=2 v=0.0 eta=0.3 | A.WRITECELL() ; C.MAJ(p=1.0)` | True | ARTIFACT |  |  |
| qd:031 | `C kC=2 | C.GAME(g=(1.5, -0.5)) ; C.CYCLE(th=6,p=0.03)` | True | COMBINATION EFFECT | g_corrlen, g_activity_cv, g_n_domains, g_n_types, g_reversals, g_largest_domain |  |
| qd:032 | `CA kC=4 kA=1 v=0.0 eta=0.3 | C.SCHELL(th=0.3) ; A.WRITECELL()` | True | DECOMPOSES | g_reversals | [0] |
| qd:033 | `CAF kC=1 kA=1 v=0.3 eta=0.3 fields=1 D=[0.05] | A.CLIMB(f=0) ; A.TURNCELL()` | True | COMBINATION EFFECT | a_density_cv, a_nn_ratio, a_local_align |  |
| qd:034 | `CA kC=3 kA=2 v=0.0 eta=1.0 | C.KURA(K=0.1) ; A.WRITECELL()` | True | DECOMPOSES | p_local_r, p_r | [0] |
| qd:035 | `AF kA=3 v=2.0 eta=0.3 fields=1 D=[0.02] | A.CLIMB(f=0) ; A.ALIGN(r=3.0)` | True | DECOMPOSES | a_step | [0] |
| qd:036 | `AF kA=1 v=0.3 eta=0.3 fields=1 D=[0.02] | A.ALIGN(r=3.0) ; A.DEPOSIT(f=0,p=0.1)` | False | COMBINATION EFFECT | f_std, f_n_blobs |  |
| qd:037 | `CA kC=3 kA=4 v=0.0 eta=1.0 | C.SCHELL(th=0.2) ; A.WRITECELL()` | True | COMBINATION EFFECT | g_n_domains, g_type_entropy, g_reversals, g_interface |  |
| qd:038 | `CA kC=3 kA=3 v=0.0 eta=0.0 | A.ALIGN(r=3.0) ; A.WRITECELL()` | True | DECOMPOSES | a_local_align | [0] |
| qd:039 | `C kC=2 | C.GAME(g=(1.5, -0.5)) ; C.CONTAG(a=0,b=1,th=4,p=0.03)` | True | COMBINATION EFFECT | g_corrlen, g_n_domains, g_activity_cv, g_reversals, g_largest_domain |  |
| qd:040 | `C kC=2 | C.GAME(g=(1.3, 0.0)) ; C.CONTAG(a=1,b=0,th=2,p=0.03)` | True | COMBINATION EFFECT | g_n_domains, g_corrlen, g_reversals, g_largest_domain, g_moran |  |
| qd:041 | `CAF kC=2 kA=1 v=0.3 eta=0.3 fields=1 D=[0.05] | A.CLIMB(f=0) ; A.WRITECELL()` | True | COMBINATION EFFECT | g_activity_trend, g_activity_cv, g_reversals |  |
| qd:042 | `AF kA=1 v=0.3 eta=0.3 fields=1 D=[0.02] | A.ALIGN(r=0.5) ; A.DEPOSIT(f=0,p=0.1)` | False | DECOMPOSES | a_density_cv, a_local_align, a_nn_ratio, a_ang_mom | [0] |
| qd:043 | `CA kC=3 kA=1 v=0.5 eta=0.0 | C.MAJ(p=0.001) ; A.TURNCELL()` | True | DECOMPOSES | g_n_domains, g_corrlen, g_reversals, g_activity_trend | [0] |
| qd:044 | `C kC=2 | C.COPY(p=0.3) ; C.MAJ(p=0.2) ; C.GAME(g=(1.5, 0.5))` | True | COMBINATION EFFECT | g_n_domains, g_corrlen, g_reversals, g_largest_domain, g_moran |  |
| qd:045 | `CF kC=4 fields=2 D=[0.05, 0.02] | C.COPY(p=0.1) ; C.SENSE(a=0,b=2,f=0,th=0.5)` | False | COMBINATION EFFECT | g_n_domains, g_corrlen, g_n_types, g_reversals, g_type_entropy, g_moran |  |
| qd:046 | `CF kC=4 fields=2 D=[0.05, 0.02] | C.COPY(p=0.2) ; C.EMIT(b=0,f=1,p=1.0)` | False | COMBINATION EFFECT | f_std, f_n_blobs, f_peak_sharp |  |
| qd:047 | `CAF kC=1 kA=3 v=0.1 eta=0.1 fields=1 D=[0.05] | F.FEED(f=0,F=0.03) ; A.CLIMB(f=0)` | True | ARTIFACT |  |  |
| qd:048 | `CA kC=3 kA=1 v=0.3 eta=3.14159 | A.ATTRACT(b=any,r=3.0) ; C.COPY(p=1.0) ; A.TURNCELL()` | True | DECOMPOSES | g_n_domains | [1] |
| qd:049 | `CA kC=3 kA=2 v=0.0 eta=0.0 | C.CYCLE(th=6,p=0.003) ; A.TURNCELL()` | True | ARTIFACT |  |  |
| qd:050 | `CF kC=4 fields=2 D=[0.05, 0.02] | C.SENSE(a=0,b=3,f=1,th=0.05) ; C.SCHELL(th=0.3)` | True | COMBINATION EFFECT | g_n_domains, g_activity_trend, g_reversals |  |
| qd:051 | `CF kC=3 fields=2 D=[0.05, 0.02] | C.SCHELL(th=0.2) ; C.SENSE(a=0,b=2,f=1,th=0.05)` | True | ARTIFACT |  |  |
| qd:052 | `CF kC=2 fields=2 D=[0.02, 0.01] | F.AUTOCAT(f=0,g=1,p=1.0) ; C.SENSE(a=0,b=1,f=1,th=0.8)` | True | COMBINATION EFFECT | g_largest_domain, g_n_types |  |
| qd:053 | `CF kC=4 fields=2 D=[0.05, 0.02] | C.MAJ(p=0.001) ; C.SENSE(a=0,b=1,f=0,th=0.5)` | False | COMBINATION EFFECT | g_n_domains, g_n_types, g_activity_trend, g_type_entropy |  |
| qd:054 | `CAF kC=4 kA=2 v=0.3 eta=0.3 fields=1 D=[0.05] | A.CLIMB(f=0) ; C.CYCLE(th=3,p=0.05)` | True | DECOMPOSES | g_n_domains, g_corrlen, g_moran, g_cyclic_changes | [1] |
| qd:055 | `CA kC=3 kA=3 v=0.1 eta=0.3 | C.SCHELL(th=0.7) ; C.MAJ(p=1.0) ; A.WRITECELL()` | True | DECOMPOSES | g_n_domains, g_largest_domain, g_reversals, g_activity_trend | [1, 2] |
| qd:056 | `CA kC=3 kA=3 v=0.1 eta=0.1 | A.ALIGN(r=1.0) ; C.MAJ(p=1.0) ; A.WRITECELL()` | True | COMBINATION EFFECT | g_reversals |  |
| qd:057 | `CA kC=3 kA=1 v=0.5 eta=0.0 | A.REPEL(b=any,r=1.0) ; C.SCHELL(th=0.4) ; A.WRITECELL()` | True | COMBINATION EFFECT | g_n_domains, g_reversals, g_largest_domain |  |
| qd:058 | `CF kC=4 fields=2 D=[0.02, 0.01] | C.MAJ(p=0.005) ; C.EMIT(b=0,f=0,p=0.002)` | False | DECOMPOSES | g_n_domains, g_corrlen, g_activity_trend, g_moran | [0] |
| qd:059 | `CAF kC=4 kA=1 v=2.0 eta=0.3 fields=1 D=[0.02] | C.CYCLE(th=1,p=0.2) ; A.TURNCELL()` | False | DECOMPOSES | g_interface, g_reversals | [0] |
| qd:060 | `CAF kC=2 kA=1 v=0.3 eta=3.14159 fields=1 D=[0.05] | A.CLIMB(f=0) ; A.WRITECELL() ; C.GAME(g=(1.5, 0.5))` | True | COMBINATION EFFECT | g_n_domains, g_activity_trend, g_corrlen, g_reversals, g_largest_domain |  |
| qd:061 | `C kC=2 | C.FLIP(a=1,b=0,p=0.2) ; C.KURA(K=1.0) ; C.GAME(g=(0.5, 0.5)) ; C.MAJ(p=1.0)` | False | DECOMPOSES | p_local_r, p_freq_spread, p_vortex_density | [1] |
| qd:062 | `CF kC=2 fields=2 D=[0.2, 0.2] | C.CONTAG(a=0,b=1,th=3,p=0.005) ; C.EMIT(b=1,f=0,p=0.001)` | False | DECOMPOSES | g_n_domains, g_activity_trend, g_n_types, g_largest_domain | [0] |
| qd:063 | `CA kC=2 kA=3 v=0.5 eta=0.0 | C.KURA(K=0.2) ; A.WRITECELL() ; C.MAJ(p=1.0)` | False | DECOMPOSES | g_n_domains, g_activity_trend, g_reversals, g_largest_domain | [1, 2] |
| qd:064 | `CAF kC=2 kA=3 v=0.3 eta=0.3 fields=1 D=[0.05] | C.CONTAG(a=1,b=0,th=4,p=1.0) ; C.EMIT(b=0,f=0,p=1.0)` | False | COMBINATION EFFECT | f_mean, f_change, f_n_blobs, f_std |  |
| qd:065 | `CAF kC=2 kA=4 v=2.0 eta=0.0 fields=2 D=[0.02, 0.05] | A.CLIMB(f=0) ; F.AUTOCAT(f=0,g=1,p=1.0)` | False | COMBINATION EFFECT | a_density_cv, a_local_align, a_nn_ratio |  |
| qd:066 | `CF kC=4 fields=2 D=[0.05, 0.02] | C.SCHELL(th=0.2) ; C.EMIT(b=3,f=1,p=1.0) ; F.DECAY(f=1,d=0.1)` | True | ARTIFACT |  |  |
| qd:067 | `C kC=2 | C.FLIP(a=1,b=0,p=0.5) ; C.MAJ(p=0.5) ; C.KURA(K=1.0) ; C.GAME(g=(1.5, 0.5))` | True | DECOMPOSES | g_reversals | [1] |
| qd:068 | `CF kC=4 fields=2 D=[0.05, 0.02] | C.SCHELL(th=0.2) ; C.EMIT(b=2,f=1,p=1.0) ; C.MAJ(p=0.002)` | True | DECOMPOSES | g_n_domains, g_corrlen, g_activity_trend, g_activity_cv | [0, 2] |
| qd:069 | `CF kC=4 fields=2 D=[0.05, 0.02] | C.SCHELL(th=0.2) ; C.EMIT(b=2,f=1,p=1.0) ; C.SWAP(p=0.02)` | True | DECOMPOSES | g_n_domains, g_activity_trend | [0, 2] |
| qd:070 | `C kC=2 | C.GAME(g=(1.5, -0.5)) ; C.FLIP(a=1,b=0,p=0.003) ; C.MAJ(p=0.01) ; C.COPY(p=0.001)` | True | DECOMPOSES | g_corrlen, g_reversals, g_type_entropy, g_moran | [0, 1] |
| qd:071 | `CF kC=4 fields=2 D=[0.05, 0.02] | F.AUTOCAT(f=0,g=1,p=0.5) ; C.SCHELL(th=0.2) ; C.EMIT(b=2,f=1,p=1.0)` | False | COMBINATION EFFECT | g_activity_trend, g_interface |  |
| qd:072 | `CAF kC=4 kA=4 v=0.1 eta=0.1 fields=2 D=[0.1, 0.05] | C.SCHELL(th=0.7) ; A.CLIMB(f=0) ; C.MAJ(p=0.1)` | True | DECOMPOSES | g_corrlen, g_activity_trend, g_n_types, g_n_domains, g_moran | [0, 2] |
| qd:073 | `CAF kC=4 kA=4 v=2.0 eta=0.3 fields=2 D=[0.05, 0.2] | C.SCHELL(th=0.2) ; C.KURA(K=0.001) ; A.WRITECELL()` | True | DECOMPOSES | g_interface, g_reversals | [0, 1] |
| qd:074 | `CAF kC=4 kA=1 v=0.0 eta=0.1 fields=2 D=[0.2, 0.01] | C.CYCLE(th=1,p=0.05) ; A.DEPOSIT(f=1,p=0.5)` | False | DECOMPOSES | g_interface, g_cyclic_changes | [0] |
| qd:075 | `CF kC=4 fields=2 D=[0.05, 0.02] | C.SCHELL(th=0.2) ; C.FLIP(a=1,b=2,p=0.2) ; C.EMIT(b=2,f=1,p=1.0)` | True | ARTIFACT |  |  |
| qd:076 | `C kC=2 | C.CYCLE(th=4,p=0.003) ; C.KURA(K=0.05) ; C.GAME(g=(1.9, 0.0)) ; C.GAME(g=(0.5, -0.5)) ; C.MAJ(p=0.01)` | False | DECOMPOSES | g_n_types, g_n_domains, g_reversals | [2, 3, 4] |
| qd:077 | `CF kC=4 fields=2 D=[0.05, 0.02] | C.SCHELL(th=0.2) ; C.CONTAG(a=0,b=1,th=7,p=0.001) ; C.EMIT(b=2,f=1,p=0.3)` | True | ARTIFACT |  |  |
| qd:078 | `CAF kC=4 kA=4 v=0.1 eta=0.1 fields=2 D=[0.1, 0.05] | C.CONTAG(a=3,b=1,th=7,p=1.0) ; C.SCHELL(th=0.7) ; A.CLIMB(f=0)` | False | DECOMPOSES | g_n_domains, g_corrlen, g_activity_trend, g_type_entropy | [0, 1] |
| qd:079 | `CAF kC=4 kA=4 v=0.1 eta=0.1 fields=2 D=[0.1, 0.05] | C.SCHELL(th=0.7) ; C.CYCLE(th=4,p=0.2) ; C.EMIT(b=2,f=0,p=0.2)` | True | DECOMPOSES | g_activity_trend | [0, 1] |
| qd:080 | `CA kC=3 kA=1 v=2.0 eta=0.0 | A.WRITECELL() ; A.QUORUM(th=4,r=2.0) ; C.CYCLE(th=2,p=0.2) ; C.CYCLE(th=2,p=1.0)` | True | DECOMPOSES | g_n_domains, g_corrlen, g_activity, g_moran | [0, 2, 3] |
| qd:081 | `CA kC=4 kA=4 v=0.0 eta=0.1 | C.KURA(K=1.0) ; A.ATTRACT(b=1,r=3.0) ; C.SCHELL(th=0.4) ; C.CONTAG(a=2,b=3,th=1,p=1.0) ; A.WRITECELL()` | False | DECOMPOSES | g_n_domains, g_corrlen, g_largest_domain, g_type_entropy, g_moran | [2, 3, 4] |
