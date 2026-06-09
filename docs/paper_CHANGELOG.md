# Paper Draft Changelog

Per-sprint mechanical changes to docs/paper_section*_draft.md files. Updates
applied by the orchestrator chain. Voice/framing changes flagged for chat-led
review at the next paper-review checkpoint.

## Sprint 74 (2026-06-09)
- epc/models/stochastic_resonance.py: new model — bistable double-well (Gammaitoni 1998 / Collins 1995) + threshold unit (T1b). Multi-trial noise-sweep design: n_trials independent runs per noise level.
- epc/detectors/p26_stochastic_resonance.py: new detector — coherent response |⟨x·signal⟩|, time-shuffle null at peak noise level, three-tier inverted-U detection.
- epc/orchestration.py: +1 model (bistable_double_well, noise_sweep_timeseries), +1 detector (P26). Counts 25 models × 24 detectors, 106 compatible pairs. New substrate: noise_sweep_timeseries.
- tests/test_stochastic_resonance_p26_e2e.py: 12 tests. Determinism, DEFINITIVE on canonical double-well, suprathreshold negative control, T1a observation bundle, threshold unit T1b.
- tests/test_cross_model.py: +2 T1b cross-model tests (P26 threshold unit generalization).
- analysis/reproductions/p26_collins.py: dim1 reproduction — gain 0.855 (> 0.05 tolerance), interior peak, DEFINITIVE.
- analysis/reproductions/p26_multiseed.py: dim2 — 20-seed campaign.
- docs/methods_notes/p26_methods.md: dim3 — coherent response metric, time-shuffle null, multi-trial design, design decisions.
- docs/observation_schema.md: +1 noise-sweep-timeseries bundle schema (P26).
- docs/depth_gap.md: +P26 row (dim1-3 PASS, dim4 pending). Implemented count 23→24.
- REPLICATION_NOTES.md: Sprint 74 P26 section added.
- docs/paper_section4_draft.md: §4.26 P26 stochastic resonance.
- docs/paper_CHANGELOG.md: this entry.

## Sprint 73 (2026-06-09)
- epc/phase2a/synthetic.py: added `scalar_timeseries` format to all 10 Class A generators. Uncontrolled drift trajectories (dx = perturbation·dt + noise) produce growth_ratio >> 2.0 → P24 rejects at screening.
- epc/phase2a/detector_invariance.py: +P24 invariance flags (permutation_invariant=True, time_shuffle_invariant=True). First run confirmed time_shuffled FP at DEFINITIVE; flag corrected to True (deviation integral is order-invariant for constant dt).
- epc/phase2a/catalog.py: +1 substrate (P24_proportional_homeostat) + generator + PATTERN_TO_SUBSTRATE_ID entry + scalar_timeseries format adapter.
- epc/phase2a/structured.py: +2 scalar_timeseries supplements (passive_ou_decay, uncontrolled_random_walk_scalar).
- epc/phase2a/failed_regimes/p24_homeostasis.py: new — Class C failed regimes (gain_zero_drift, no_perturbation).
- epc/phase2a/panel.py: scalar_timeseries format dispatch in Class A kwargs.
- analysis/run_phase2a_panel.py: P24 panel wiring (build_p24_positives, make_p24_detector_fn, run_p24).
- analysis/outputs/p24_phase2a_panel.json: Phase-2a v1.2 panel — TNR=1.000, Cohen's d=+inf, verdict=PASS.
- docs/depth_gap.md: P24 dim4 pending→PASS, grade GAP→AT-DEPTH. AT-DEPTH count 21→22, gap count 2→1.
- REPLICATION_NOTES.md: Sprint 73 P24 Phase-2a dim4 section added.
- docs/paper_section4_draft.md: §4.28 Phase-2a panel results appended.
- docs/paper_CHANGELOG.md: this entry.

## Sprint 72 (2026-06-09)
- epc/models/homeostasis.py: new model — Ashby (1956) proportional + integral homeostat. Scalar regulated variable with negative-feedback control under external perturbation.
- epc/detectors/p24_homeostasis.py: new detector — deviation-integral vs surrogate uncontrolled null, growth-ratio screening, T1a observation-bundle adapter.
- epc/orchestration.py: +1 model (proportional_homeostat, scalar_timeseries), +1 detector (P24). Counts 24 models × 23 detectors, 105 compatible pairs. New substrate: scalar_timeseries.
- tests/test_homeostasis_p24_e2e.py: 13 tests. Determinism, DEFINITIVE on canonical gain=5, gain=0 negative control, T1a observation bundle, pulse perturbation, metadata interaction.
- tests/test_cross_model.py: +2 T1b cross-model tests (P24 integral controller generalization).
- analysis/outputs/p24_homeostasis_reproduction.json: dim1 reproduction — deviation ratio 0.0027 (< 0.30 tolerance).
- analysis/outputs/p24_multiseed.json: dim2 — 20-seed campaign, deviation integral CV=0.8%, all DEFINITIVE.
- docs/methods_notes/p24_methods.md: dim3 — surrogate null, T1a contract, proportional vs integral control.
- docs/observation_schema.md: new — T1a scalar-regulated-variable bundle schema.
- docs/depth_gap.md: +P24 row (dim1–3 PASS, dim4 pending). Implemented count 22→23, gap count 1→2.
- docs/paper_section4_draft.md: §4.28 stub — P24 characterization.
- docs/paper_CHANGELOG.md: this entry.

## Sprint 70 (2026-06-08)
- epc/phase2a/failed_regimes/p19_informed_minority.py: new — Class C failed regimes (5 rho_zero + 5 bias_zero).
- epc/phase2a/catalog.py: +1 substrate (P19_informed_minority) + generator + PATTERN_TO_SUBSTRATE_ID entry.
- analysis/run_phase2a_panel.py: P19 panel wiring — detect_p19 wrapper with early-window leadership content prerequisite (Couzin 2005: informed minority must lead during convergence).
- analysis/outputs/p19_phase2a_panel.json: Phase-2a v1.2 panel — TNR=0.960, Cohen's d=5.418, verdict=PASS.
- docs/depth_gap.md: P19 dim4 pending→PASS, grade GAP→AT-DEPTH. AT-DEPTH count 20→21, gap count 2→1. Completes Milestone B Wave 1.
- REPLICATION_NOTES.md: Sprint 70 P19 Phase-2a dim4 section added.
- docs/paper_CHANGELOG.md: this entry.

## Sprint 69 (2026-06-08)
- epc/models/informed_minority.py: new model — Couzin et al. (2005) informed-minority flocking. Vicsek-style alignment + ω-weighted bias toward preferred direction for informed fraction ρ.
- epc/detectors/p19_emergent_leadership.py: new detector — group directional accuracy + label-shuffle influence asymmetry (directional pull) + guidance efficacy.
- epc/orchestration.py: +1 model (informed_minority, continuous_2d), +1 detector (P19). Counts 23 models × 22 detectors, 104 compatible pairs.
- tests/test_informed_minority_p19_e2e.py: 19 tests. Determinism, DEFINITIVE on canonical ρ=0.1, ρ=0 negative control, registry.
- tests/test_transfer_matrix_counts.py: EXPECTED updated for 23×22 registry.
- analysis/reproductions/p19_couzin2005.py: dim1 reproduction — accuracy rises with ρ (Spearman ρ=1.0); ρ=0 chance-level; ρ=0.025 near-perfect. All 5 tolerances PASS.
- analysis/reproductions/p19_multiseed.py: dim2 — 20 seeds at ρ=0.1, accuracy = 1.000 ± 0.000 (CV=0%).
- analysis/outputs/p19_couzin2005_reproduction.json: dim1 output.
- analysis/outputs/p19_multiseed.json: dim2 output.
- docs/methods_notes/p19_methods.md: dim3 — Vicsek+informed bias dynamics, label-shuffle null, pull metric, P5/P17/P18 distinctness.
- docs/depth_gap.md: P19 row added (dim1-3 PASS, dim4 pending). Patterns audited 21→22; gap count 1→2 (P12 + P19).

## Sprint 68 (2026-06-08)
- epc/phase2a/failed_regimes/p17_collective_sensing.py: new — Class C failed regimes (5 social_off + 5 field_too_strong).
- epc/phase2a/detector_invariance.py: P17 permutation_invariant False→True (group CI is CoM-based, agent-index invariant).
- epc/phase2a/catalog.py: +1 substrate (P17_collective_sensing) + generator + PATTERN_TO_SUBSTRATE_ID entry.
- analysis/run_phase2a_panel.py: P17 panel wiring — history-based CI adapter with 3 literature-grounded prerequisites (field_samples, individual SNR, social cohesion).
- analysis/outputs/p17_phase2a_panel.json: Phase-2a v1.2 panel — TNR=1.000, Cohen's d=11.117, verdict=PASS.
- docs/depth_gap.md: P17 dim4 pending→PASS, grade GAP→AT-DEPTH. AT-DEPTH count 19→20, gap count 2→1.
- REPLICATION_NOTES.md: Sprint 68 P17 Phase-2a section added.
- docs/paper_CHANGELOG.md: this entry.

## Sprint 67 (2026-06-08)
- epc/models/collective_sensing.py: new model — Berdahl et al. (2013) collective gradient sensing. Speed-modulation mechanism in periodic 2D domain with Gaussian scalar field.
- epc/detectors/p17_collective_sensing.py: new detector — group-size-scaling CI test (N-sweep + α=0 null).
- epc/orchestration.py: +1 model (collective_sensing, continuous_2d), +1 detector (P17). Counts 22 models × 21 detectors, 95 compatible pairs.
- tests/test_collective_sensing_p17_e2e.py: 12 tests (10 fast + 2 slow). Determinism, N-scaling, negative controls, registration.
- tests/test_transfer_matrix_counts.py: EXPECTED updated for 22×21 registry.
- analysis/reproductions/p17_berdahl2013.py: dim1 reproduction — CI slope vs log(N) = 0.133; N=1 at chance; N=50 positive. All tolerances PASS.
- analysis/reproductions/p17_multiseed.py: dim2 — 20 seeds, CI = 0.394 ± 0.130 (CV=33%), 100% positive.
- analysis/outputs/p17_berdahl2013_reproduction.json: dim1 output.
- analysis/outputs/p17_multiseed.json: dim2 output.
- docs/methods_notes/p17_methods.md: dim3 — speed-modulation mechanism, SNR∝√N, α=0 null, P5 distinctness.
- docs/depth_gap.md: P17 row added (dim1-3 PASS, dim4 pending). Patterns audited 20→21; gap count 1→2.

## Sprint 63 (2026-05-30)
- analysis/reproductions/p12_reichenbach2007.py: rewritten for L=200, zero-padded FFT ring-peak estimator (replaces ACF first-zero). 5 M values, 15 seeds each, T_eq=2500 gen.
- analysis/outputs/p12_reichenbach2007_reproduction.json: overwritten (sprint=63). Slope=0.161, R²=0.792, FAIL. Accepted as documented finite-size measurement limitation.
- REPLICATION_NOTES.md: Sprint 63 Dim1 final-attempt section added with per-M table and accepted-limitation verdict.
- docs/depth_gap.md: P12 dim1 notes updated with Sprint 63 result. C2/C3 carry-forwards reclassified as closed-as-documented-limitation.
- §4.11 paper_section4_draft.md: Sprint 63 L=200 FFT paragraph appended. Accepted-limitation conclusion documented.
- §6 paper_section6_draft.md: Sprint 63 paragraph added. AT-DEPTH count 18/19 unchanged.

## Sprint 62 (2026-05-30)
- epc/phase2a/failed_regimes/p8_nagel_schreckenberg.py: Class C regime correction — density range linspace(0.05,0.20,10) → linspace(0.02,0.07,10). Original densities ρ≥0.1167 above empirical jamming onset ρ_c≈0.10 at v_max=5, p=0.3 (brief-author error).
- analysis/outputs/p8_phase2a_panel.json: re-run — overall TNR 0.714→1.000, fai 0.400→1.000, Cohen's d 1.772→+inf. Verdict PARTIAL→PASS.
- REPLICATION_NOTES.md: Sprint 62 P8 panel re-run section added with jamming-onset calibration table.
- docs/depth_gap.md: P8 dim4 PARTIAL→PASS, grade GAP→AT-DEPTH. AT-DEPTH count 17→18.
- §4.14 paper_section4_draft.md: Sprint 62 panel re-run paragraph appended. P8 AT-DEPTH.
- §6.11 paper_section6_draft.md: Sprint 62 paragraph added. AT-DEPTH count 17→18.
- Carry-forward CLOSED: C-p8-class-c-near-onset.

## Sprint 61 (2026-05-30)
- epc/detectors/p1_aggregation.py: added multi-cluster prerequisite (Schelling 1971) — rejects substrates where every non-empty type forms a single contiguous block (monotonic spatial partition). Uses scipy.ndimage.label, 8-connected. Gradient rejected; canonical Schelling positive (10–20 components/type) unaffected.
- epc/phase2a/failed_regimes/p1_schelling.py: Class C regime correction — threshold range linspace(0.05,0.25,10) → linspace(0.01,0.10,10), grid_size 32→50. Original thresholds 0.161–0.250 above empirical critical threshold ≈0.13 at density=0.9 (brief-author error).
- analysis/outputs/p1_phase2a_panel.json: re-run — overall TNR 0.731→1.000, syn 0.889→1.000, fai 0.400→1.000, Cohen's d 1.740→+inf. Verdict PARTIAL→PASS.
- tests/test_p1_multi_cluster_prereq.py: 2 regression tests — gradient rejected, canonical positive not regressed.
- REPLICATION_NOTES.md: Sprint 61 P1 panel re-run section added.
- docs/depth_gap.md: P1 dim4 PARTIAL→PASS, grade GAP→AT-DEPTH. AT-DEPTH count 16→17.
- §3.5 paper_section3_draft.md: multi-cluster prerequisite paragraph added (alongside P11/P22/P27 prereqs).
- §4.6 paper_section4_draft.md: Sprint 61 panel re-run paragraph appended. P1 AT-DEPTH.
- §6.11 paper_section6_draft.md: Sprint 61 paragraph added. AT-DEPTH count 16→17.
- Carry-forwards CLOSED: C-p1-linear-gradient-fp, C-p1-class-c-subthreshold-fp.

## Sprint 60 (2026-05-30)
- analysis/p2_multiseed.py: new 20-seed multi-seed script for P2 MIPS (two_phase_score + density-speed r at Pe=100, φ=0.5, N=800).
- analysis/p21_multiseed.py: new 20-seed multi-seed script for P21 HK opinion (cluster_count at ε=0.20, N=100).
- analysis/p22_multiseed.py: new 20-seed multi-seed script for P22 SIR wavefront speed (paper-exact Datta-Acharyya CA, L=200).
- analysis/outputs/p2_multiseed.json: P2 aggregate — score=0.1134±0.0790 (CV=69.7%), r=−0.9585±0.0196 (CV=2.1%).
- analysis/outputs/p21_multiseed.json: P21 aggregate — cluster_count=1.90±0.31 (CV=16.2%, median=2).
- analysis/outputs/p22_multiseed.json: P22 aggregate — speed=0.4606±0.0163 (CV=3.5%, 19/20 valid).
- REPLICATION_NOTES.md: P2, P21, P22 Dim2 Multi-Seed Extension sections added. dim2 PARTIAL→PASS for all three.
- docs/depth_gap.md: P2, P21, P22 dim2→PASS, grade GAP→AT-DEPTH. AT-DEPTH count 13→16. Sprint 60 finding added.
- §4.15 (P2) paper_section4_draft.md: appended Sprint 60 multi-seed dispersion paragraph. P2 AT-DEPTH.
- §4.9 (P21) paper_section4_draft.md: appended Sprint 60 multi-seed dispersion paragraph. P21 AT-DEPTH.
- §4.10 (P22) paper_section4_draft.md: appended Sprint 60 multi-seed dispersion paragraph. P22 AT-DEPTH.
- §6 paper_section6_draft.md: updated AT-DEPTH count 13→16, Sprint 60 paragraph added. P2, P21, P22 listed.

## Sprint 59 (2026-05-29)
- analysis/reproductions/p12_reichenbach2007.py: rewritten for near-M_c dense sweep — M ∈ [2e-4, 5e-4] (7 points, 30 seeds each), fit_mask restricting slope fit to M ∈ [2e-4, 4.5e-4] with n_valid ≥ 15, updated JSON schema (sprint=59, m_c, fit_M_values, n_fit_points, near_extinction, n_fit_points_pass, fit_exclusion_note). spiral_wavelength() and run_one_seed() unchanged.
- analysis/outputs/p12_reichenbach2007_reproduction.json: overwritten — sprint=59, slope=0.244, R²=0.918, n_fit_points=6, overall_pass=False. Sprint 58 result superseded.
- tests/test_p12_reichenbach2007_reproduction.py: rewritten — B.1 (slow: slope/R²/overall_pass/n_fit_points), B.2 (slow: per-point rel_error for non-extinction points), B.3 (non-slow: sprint=59 JSON schema).
- REPLICATION_NOTES.md P12/Sprint 59 subsection: added with per-M table (7 points), slope=0.244, R²=0.918, FAIL verdict, and finite-size ACF bias diagnosis.
- docs/depth_gap.md: P12 dim1 note appended with Sprint 59 result (slope=0.244, FAIL, ACF finite-size diagnosis); Sprint 59 finding added to Aggregate findings; C2/C3 carry-forwards updated with finite-size bias path forward. AT-DEPTH count unchanged: 13/19.
- §4.11 (P12/RPS) paper_section4_draft.md: appended Sprint 59 paragraph after Sprint 58 paragraph. Documents slope=0.244, R²=0.918, finite-size ACF compression, path forward (L≥200).
- §6 paper_section6_draft.md: Sprint 59 paragraph added after Sprint 58 paragraph. AT-DEPTH count unchanged at 13/19.

## Sprint 58 (2026-05-29)
- analysis/reproductions/p12_reichenbach2007.py: extended in-place — wide M sweep [1e-5, 5e-4] (7 points, 15 seeds, T_eq=1000 gen, parallel Pool, updated JSON schema with log_log_slope/r_squared/per_point/sprint=58, unit verification block).
- analysis/outputs/p12_reichenbach2007_reproduction.json: overwritten — sprint=58, slope=0.107, R²=0.769, overall_pass=False. Sprint 54 narrow-sweep result superseded.
- tests/test_p12_reichenbach2007_reproduction.py: new file — B.1 (slow: slope/R²/overall_pass), B.2 (slow: per-point rel_error/n_valid), B.3 (non-slow: sprint==58 superseded).
- REPLICATION_NOTES.md P12/Sprint 58 subsection: added with per-M table (7 points), slope=0.107, R²=0.769, FAIL verdict, and detailed diagnosis (formula breaks down at M ≪ M_c; flat λ ~42–44 at M ≤ 5e-5; valid test range near M_c only).
- docs/depth_gap.md: P12 dim1 note appended with Sprint 58 wide-sweep result (slope=0.107, FAIL, diagnosis); Sprint 58 finding added to Aggregate findings; C2/C3 carry-forwards updated with new diagnosis and near-M_c sweep recommendation. AT-DEPTH count unchanged: 13/19.
- §4.12 (P12/RPS) paper_section4_draft.md: appended Sprint 58 paragraph after Sprint 54 paragraph. Documents slope=0.107, formula breakdown at M ≪ M_c, valid range near M_c, path forward.
- §6 paper_section6_draft.md: Sprint 58 paragraph added after Sprint 57 paragraph. AT-DEPTH count unchanged at 13/19.

## Sprint 57 (2026-05-29)
- docs/methods_notes/p2_methods.md: new methods note for P2 (MIPS / ABP). Documents Fily-Marchetti overdamped Langevin equations, two_phase_coexistence_score primary metric, Hartigan-dip unusability (ADR 44), FFT structure-factor context, burn-in / nucleation-lag requirements, mechanistic-null metadata flags (ADR 43), known limitations.
- docs/methods_notes/p21_methods.md: new methods note for P21 (polarization / Hegselmann-Krause). Documents synchronous HK update rule, L∞ convergence criterion (1e-8 default / 1e-6 reproduction), sorted-gap cluster counting (gap = ε/2), ε_c boundary zone, C-p21-time-shuffled-fp context, known limitations.
- docs/methods_notes/p22_methods.md: new methods note for P22 (information cascade / SIR CA). Documents S/I/R encoding (0/1/2), synchronous update, independent-neighbours infection probability, irreversibility prereq (Sprint 41), percolation threshold context (p_c ≈ 0.038 Moore at q=0.1), model difference between Sprint 51 reproduction and epc.models.sir_epidemic, known limitations.
- REPLICATION_NOTES.md P2/Sprint 52 section: appended "Dim3 Methods Note — Sprint 57" with coverage summary and `docs/methods_notes/p2_methods.md` reference. dim3 PARTIAL→PASS verdict.
- REPLICATION_NOTES.md P21/Sprint 53 section: appended "Dim3 Methods Note — Sprint 57" with coverage summary and `docs/methods_notes/p21_methods.md` reference. dim3 PARTIAL→PASS verdict.
- REPLICATION_NOTES.md P22/Sprint 51 section: appended "Dim3 Methods Note — Sprint 57" with coverage summary and `docs/methods_notes/p22_methods.md` reference. dim3 PARTIAL→PASS verdict.
- docs/depth_gap.md: P2 row dim3 PARTIAL→PASS; P21 row dim3 PARTIAL→PASS (effort L→M); P22 row dim3 PARTIAL→PASS (effort L→M); Sprint 57 finding added; C5 carry-forward CLOSED. AT-DEPTH count unchanged: 13/19 (all three still have dim2 PARTIAL).
- §3.4: added 3-sentence methods-note overview referencing p2/p21/p22_methods.md.
- §4.2 (GH/P13): added cross-reference to p22_methods.md §5 for P13/SIR boundary test.
- §4.9 (P21 HK): added cross-reference to p21_methods.md for implementation details.
- §4.10 (P22 SIR): added cross-reference to p22_methods.md for implementation details.
- §4.15 (P2 ABP): added cross-reference to p2_methods.md for implementation details.
- §6: Sprint 57 paragraph added; AT-DEPTH count 13/19 unchanged; C5 CLOSED noted.

## Sprint 56 (2026-05-29)
- analysis/p6_multiseed.py: new multi-seed analysis script for P6 D'Orsogna milling. N=100, random init, warmup=2500 steps, measure=500 steps, 20 seeds (seeds 100–119).
- analysis/p12_multiseed.py: new multi-seed analysis script for P12 spatial RPS. L=100, M=1e-4, T_eq=500 gen, T_measure=200 gen, 20 seeds (seeds 100–119).
- analysis/outputs/p6_multiseed.json: |L| = 0.9818 ± 0.0301 (CV = 3.1%) across 20 seeds. All seeds form stable mills (min |L| = 0.884).
- analysis/outputs/p12_multiseed.json: λ = 52.1 ± 10.4 (CV = 20.0%) across 20 seeds. All 20 seeds valid.
- REPLICATION_NOTES.md P6 section: appended "Dim2 Multi-seed Extension — Sprint 56" with per-seed |L| table and aggregate statistics. dim2 PARTIAL→PASS verdict.
- REPLICATION_NOTES.md P12 section: appended "Dim2 Multi-seed Extension — Sprint 56" with per-seed λ table and aggregate statistics. dim2 PARTIAL→PASS verdict.
- docs/depth_gap.md: P6 row dim2 PARTIAL→PASS, grade GAP→AT-DEPTH; P12 row dim2 PARTIAL→PASS, grade stays GAP (dim1 still PARTIAL); AT-DEPTH count 12→13; Gap count 7→6; Sprint 56 finding added.
- §4.4 (P6 D'Orsogna milling): appended "Multi-seed robustness (Sprint 56)" paragraph with |L| ± std, CV.
- §4.11 (P12 spatial RPS): appended "Multi-seed robustness at fixed M (Sprint 56)" paragraph with λ ± std, CV.
- §6 aggregate: Sprint 56 paragraph added; AT-DEPTH count 12→13/19; aggregate sentence updated to include P6.

## Sprint 55 (2026-05-29)
- analysis/p14_multiseed.py: new multi-seed analysis script for P14 BTW sandpile. L=32, n_drive=30,000, n_burn=3,000, 20 seeds (seeds 100–119).
- analysis/outputs/p14_multiseed.json: τ = 1.2914 ± 0.0012 (CV = 0.09%) across 20 seeds. All seeds detect SOC at CONFIRMATION tier.
- REPLICATION_NOTES.md P14 (BTW) section: appended "Dim2 Multi-seed Extension — Sprint 55" with per-seed τ table and aggregate statistics. dim2 PARTIAL→PASS verdict.
- docs/depth_gap.md: P14 row dim2 PARTIAL→PASS, grade GAP→AT-DEPTH; AT-DEPTH count 11→12; Gap count 8→7; Sprint 55 finding added; C4 carry-forward CLOSED.
- §4.7 (P14 BTW sandpile): appended "Multi-seed robustness (Sprint 55)" paragraph with τ ± std, CV, range.
- §6 aggregate: Sprint 55 paragraph added; AT-DEPTH count 11→12/19; aggregate sentence updated.

## Sprint 54 (2026-05-28)
- analysis/reproductions/p12_reichenbach2007.py: reproduction script for Reichenbach-Mobilia-Frey (2007) Fig. 2c spiral wavelength λ ~ M^(1/2). L=100, σ=μ=1, M ∈ {3e-4, 4e-4, 5e-4}, 10 seeds, radial ACF estimator, T_eq=500, T_meas=200, stride=20.
- analysis/outputs/p12_reichenbach2007_reproduction.json: passes_tolerance=False; measured slope=0.366 (target 0.5, tolerance [0.4, 0.6]); wavelengths qualitatively match formula (within 15%); tolerance failure due to narrow M range (1.67×) with ~10% per-point variance.
- REPLICATION_NOTES.md P12 section: updated Open Item #1, appended "Dim1 Reproduction — Sprint 54" with parameter table, per-M results, slope result, and PARTIAL verdict with carry-forward reasoning.
- docs/depth_gap.md: P12 row dim1 note updated (Sprint 54 attempt, slope 0.37, outside tolerance); Sprint 54 finding added; C2/C3 carry-forwards updated. dim1 PARTIAL count unchanged: 1. AT-DEPTH count unchanged: 11/19.
- §4.11 (P12 RPS): replaced "We did not replicate..." with "Numerical reproduction (Sprint 54)" paragraph including measured-vs-published table and tolerance verdict (outside [0.4, 0.6]).
- §3.6 Sprint 54: new subsection — dim1 reproduction attempt description, per-M results, cumulative dim1 table (P12 row added as OUTSIDE tolerance).
- §6 aggregate: Sprint 54 paragraph added; AT-DEPTH count unchanged at 11/19.

## Sprint 53 (2026-05-26)
- analysis/reproductions/p21_hegselmann2002.py: new reproduction script for Hegselmann & Krause (2002) Fig. 2 cluster-count vs ε curve. N=100, uniform IC, synchronous averaging, convergence tol=1e-6, 20 seeds per ε, 8 ε points (0.10–0.50).
- analysis/outputs/p21_hegselmann2002_reproduction.json: all 8 ε points PASS; median cluster counts match published ranges; ε=0.25 boundary zone documented (14/20 consensus, 6/20 two-cluster).
- REPLICATION_NOTES.md P21 section: appended "Dim1 Reproduction — Sprint 53" with parameter table, per-ε results table, boundary-zone footnote, and PASS verdict.
- docs/depth_gap.md: P21 row dim1 PARTIAL→PASS; dim1 PARTIAL count 2→1 (P12 only); Sprint 53 finding added; C2 carry-forward updated (P21 closed, P12 only remains). AT-DEPTH count unchanged: 11/19.
- §4.9 (P21 HK): appended "Numerical reproduction (Sprint 53)" paragraph with tolerance table.
- §3.6 Sprint 53: new subsection — dim1 reproduction table (cumulative through Sprint 53), P21 closure description, boundary-zone footnote.
- §6 aggregate: Sprint 53 paragraph added; AT-DEPTH count unchanged at 11/19 (P21 dims 2–3 still PARTIAL).

## Sprint 52 (2026-05-25)
- analysis/reproductions/p2_filymarchetti2012.py: new reproduction script for Fily & Marchetti (2012) PRL 108, 235702 Fig. 2 canonical MIPS state. N=800, φ=0.5, Pe=100 (MIPS) and Pe=5 (thermal), 5 seeds × 2500 steps.
- analysis/outputs/p2_filymarchetti2012_reproduction.json: passes_tolerance=True; Pe=100 two_phase_score=0.1237±0.077 (≥0.10 PASS); Pe=100 Pearson r=−0.958±0.020 (|r|≥0.70 PASS); Pe=5 score=0.052<0.08 PASS.
- REPLICATION_NOTES.md §Sprint 16: appended "Dim1 Reproduction — Sprint 52" section with parameter table, per-seed results, and PASS verdict for all three tolerance checks.
- docs/depth_gap.md: P2 row dim1 PARTIAL→PASS; dim1 PARTIAL count 3→2 (P12, P21 remain); Sprint 52 finding added; C2 carry-forward updated (P2 removed). AT-DEPTH count unchanged: 11/19.
- §4.15 (P2 MIPS): appended "Numerical reproduction (Sprint 52)" paragraph with measured-vs-published table.
- §3.6 Sprint 52: new subsection — dim1 reproduction table (cumulative through Sprint 52), P2 closure description.
- §6 aggregate: Sprint 52 paragraph added; AT-DEPTH count unchanged at 11/19 (P2 dims 2–3 still PARTIAL).

## Sprint 51 (2026-05-25)
- analysis/reproductions/p22_dattaacharyya2005.py: new reproduction script for Datta & Acharyya (2021) §3.1.1/Fig.11 wavefront speed. Implements paper's exact fixed-duration CA (t_τ=4, p2=0.10 re-infection) inline.
- analysis/outputs/p22_dattaacharyya2005_reproduction.json: reproduction output; passes_tolerance=True; measured speed=0.4612±0.0164 vs published 0.4405±0.0008 (rel error 4.7%).
- REPLICATION_NOTES.md SIR section: Open Item #1 closed; appended "Dim1 Reproduction — Sprint 51" section with parameter table, per-seed results, and PASS verdict.
- docs/depth_gap.md: P22 row dim1 PARTIAL→PASS; dim1 PARTIAL count 4→3; Sprint 51 finding added; C2 carry-forward updated (P22 removed). AT-DEPTH count unchanged: 11/19.
- §4.10 (P22 SIR): appended "Numerical reproduction (Sprint 51)" paragraph with measured-vs-published table.
- §3.6 Sprint 51: new subsection — dim1 reproduction table (cumulative through Sprint 51), P22 closure description.
- §6.11 aggregate: Sprint 51 paragraph added; AT-DEPTH count unchanged at 11/19 (P22 dims 2–3 still PARTIAL).

## Sprint 50 (2026-05-25)
- analysis/reproductions/p11_mobilia2007_fig2.py: new reproduction script for Mobilia-Georgiev-Täuber (2007) O(1/L) amplitude scaling law and coexistence verification.
- analysis/outputs/p11_mobilia2007_reproduction.json: reproduction output; all_passed=True; scaling exponent=−0.967 (R²=0.990); coexistence confirmed (FFT p2m=48.9).
- REPLICATION_NOTES.md §Sprint 11: appended "Dim1 Reproduction — Sprint 50" section with parameter table, scaling-law data, and tolerance verdict (PASS).
- docs/depth_gap.md: P11 row dim1 PARTIAL→PASS, grade GAP→AT-DEPTH; aggregate At-depth 10→11, Gap 9→8; Sprint 50 finding note added; C2 carry-forward updated (P11 removed from dim1 PARTIAL list).
- §4.12 (P11): appended "Numerical reproduction (Sprint 50)" paragraph: scaling exponent −0.967 vs published −1.0 (3.3% error), R²=0.990; large MF deviation documented as expected published finding.
- §3.6 Sprint 50: new subsection — dim1 reproduction table, P11 closure description, AT-DEPTH 11/19 milestone.
- §6.11 aggregate: Sprint 50 paragraph added; AT-DEPTH count updated 10→11; opening sentence updated to 11/19 AT-DEPTH patterns.

## Sprint 41 (2026-05-24)
- epc/detectors/p22_information_cascade.py: added Sprint 41 irreversibility prerequisite guard — `detect()` override + `_check_irreversibility_prereq()` helper. Short-circuits with `detected=False, confidence=0.0` if any cell has a backward (decreasing) state transition. Literature anchor: Datta & Acharyya (2021), Mobilia-Georgiev-Täuber (2007), Reichenbach (2007).
- tests/test_sir_p22_e2e.py: added 2 regression tests: `test_p22_short_circuits_on_lv_substrate` (guard fires on LV), `test_p22_still_fires_on_sir_canonical` (no regression on SIR DEFINITIVE).
- analysis/outputs/p22_phase2a_panel.json: re-run; TNR=1.000 (syn=1.000, cat=1.000, fai=1.000), Cohen's d=+∞, verdict=PASS. Class B false positives on LV+RPS eliminated.
- §3.5: appended P22 irreversibility-prerequisite paragraph alongside existing P11 + P27 guards.
- §4.10 P22 SIR: appended Phase-2a panel re-run paragraph (Sprint 41 PASS: TNR=1.000, d=+∞; C-p22-class-b-cascade-overlap CLOSED).
- §6.11 aggregate: updated Sprint 40/41 narrative; AT-DEPTH count unchanged at 6 (P22 dim4→PASS but dims 1–3 still PARTIAL).
- docs/depth_gap.md: P22 row dim4 PARTIAL→PASS; aggregate Sprint 41 finding note added. AT-DEPTH count remains 6/19.
- REPLICATION_NOTES.md: appended Phase-2a Panel Result (v1.2) Sprint 41 re-run section to SIR replication notes.

## Sprint 40 (2026-05-23)
- epc/phase2a/failed_regimes/p22_sir.py: corrected infection_prob range from [0.05, 0.18] to [0.005, 0.030] (all below Moore p_c≈0.038); updated docstring and description.
- epc/detectors/p27_spatial_reciprocity.py: added Sprint 40 prerequisite guard — short-circuits if `coop_fraction` absent from history; prevents out-of-domain fires on generic lattice_2d substrates.
- analysis/run_phase2a_panel.py: _augment_history_p27 converted to pass-through (no longer adds synthetic coop_fraction to generic grid histories).
- tests/test_nowak_may_p27_e2e.py: added 2 new regression tests (short-circuit without coop_fraction; no regression on Nowak-May canonical).
- §3.5: appended P27 observable-prerequisite paragraph alongside existing P11 discussion.
- §4.8 P27 Nowak-May: appended Phase-2a panel re-run paragraph (Sprint 40 PASS: TNR=1.000, d=2.950; C-p27-panel-screening-leak CLOSED; C-p27-time-shuffle-invariance VALIDATED).
- §4.10 P22 SIR: appended Phase-2a panel re-run paragraph (Sprint 40 PARTIAL: TNR=0.889, d=2.981; Class C fixed fai TNR 0.000→1.000; catalog FPs remain; C-p22-class-c-above-percolation CLOSED).
- §6.11 aggregate: updated AT-DEPTH count 5→6 (P27 advances); updated 5→6 AT-DEPTH patterns list; updated Sprint 39 finding note to Sprint 40.
- docs/depth_gap.md: P27 row dim4 PARTIAL→PASS, grade GAP→AT-DEPTH; P22 row dim4 note updated (Sprint 40 PARTIAL, TNR=0.889); aggregate AT-DEPTH 5→6, GAP 14→13.

## Sprint 39 (2026-05-23)
- analysis/run_phase2a_panel.py: added run_p27() + run_p22() + _augment_history_p27() adapter; failed-regime imports for p27_nowak_may and p22_sir.
- epc/phase2a/failed_regimes/p27_nowak_may.py: 10 extinction regimes at b ∈ linspace(2.0, 2.5, 10) (Sprint 39 Class C for P27).
- epc/phase2a/failed_regimes/p22_sir.py: 10 sub-percolation regimes at infection_prob ∈ linspace(0.05, 0.18, 10) (Sprint 39 Class C for P22).
- §4.8 P27 Nowak-May: appended Phase-2a panel FAIL paragraph (TNR=0.500, d=0.198; C-p27-panel-screening-leak carry-forward).
- §4.10 P22 SIR: appended Phase-2a panel PARTIAL paragraph (TNR=0.519, d=1.094; C-p22-class-c-above-percolation carry-forward).
- §6.11 aggregate: updated AT-DEPTH count note (remains 5/19; Sprint 39 panels did not add new AT-DEPTH patterns).

## Sprint 38 (2026-05-23)
- epc/phase2a/catalog.py: added _gen_p8_nagel_schreckenberg; updated PATTERN_TO_SUBSTRATE_ID for P8 lattice_1d Class B coverage.
- epc/phase2a/structured.py: added independent_lane_traffic + reverse_sorted_sequence lattice_1d B' supplements (C-supplements progress: lattice_1d now closed).
- 4 new tests in tests/test_phase2a_panel.py.

## Sprint 37 (2026-05-23)
- epc/phase2a/catalog.py: added _gen_p2_active_brownian + _gen_p6_dorsogna; updated PATTERN_TO_SUBSTRATE_ID for continuous_2d Class B coverage.
- epc/phase2a/structured.py: added uncorrelated_random_walks + independent_brownian_motion continuous_2d B' supplements (C-supplements progress).
- 6 new tests in tests/test_phase2a_panel.py.

## Sprint 36 (2026-05-22, paper catchup Sprints 27–35)
- §3.7 added: Phase-2a standard negative panel methodology (motivation,
  panel composition, PASS criterion + invariance flags, spec evolution as
  methodology contribution, reproducibility note)
- §4.3 P15 GoL: appended v1.1 panel result paragraph (PASS, TNR 1.000,
  Cohen's d 8.282; AT-DEPTH confirmed by panel)
- §4.5 P9 Kuramoto: appended v1.0 PARTIAL → v1.1 PARTIAL → v1.2
  PASS-with-weakness panel result paragraphs (TNR 0.952, Cohen's d 4.781,
  Class A weak due to constant_field carry-forward)
- §4.7 P14 BTW sandpile: appended v1.2 panel result paragraph (PASS,
  TNR 0.960, Cohen's d 10.585, Class C borderline at p_diss=0.350)
- §4.19 P10 chimera: appended Sprint 27 deferred Phase 1n multi-seed
  (A, β) phase boundary update — smooth basin-volume gradient finding,
  table of basin fractions across 5 β values, reconciliation with Sprint 26
- §4.20 P18 voter: appended v1.1 panel result paragraph (PASS, TNR 1.000,
  Cohen's d +inf; AT-DEPTH confirmed by panel)
- §6.11 added: aggregate AT-DEPTH count 5 / 19, GAP 14 / 19, transfer-matrix
  figure lock confirmation, note on 4-of-5 AT-DEPTH via panel

## Sprint 37 (placeholder — first orchestrator-driven sprint will write its own entry)

## Sprint 44 (2026-05-25)

- epc/phase2a/failed_regimes/p12_rps.py: NEW; Class C 10 high-mobility extinction regimes (mobility ∈ linspace(5e-3, 5e-2, 10)).
- epc/phase2a/failed_regimes/p13_gh.py: NEW; Class C 10 low-density init regimes (density ∈ linspace(0.01, 0.10, 10)).
- epc/phase2a/failed_regimes/__init__.py: p12_rps + p13_gh registered.
- analysis/run_phase2a_panel.py: `build_p12_positives()`, `make_p12_detector_fn()`, `run_p12()`, `build_p13_positives()`, `make_p13_detector_fn()`, `run_p13()` added.
- analysis/outputs/p12_phase2a_panel.json: NEW; TNR=1.000, Cohen's d=+∞, PASS. P12 dim4 PARTIAL→PASS.
- analysis/outputs/p13_phase2a_panel.json: NEW; TNR=1.000, Cohen's d=+∞, PASS. P13 advances to AT-DEPTH.
- §4.22 added: P12 RPS Phase-2a panel v1.2 PASS + P13 GH Phase-2a panel v1.2 PASS.
- §6: updated AT-DEPTH count 7→8/19; Sprint 44 narrative added.
- docs/depth_gap.md: P12 row dim4 PARTIAL→PASS (grade stays GAP: dim1+dim2 PARTIAL); P13 row dim4 PARTIAL→PASS, grade GAP→AT-DEPTH; aggregate counts updated 7→8 AT-DEPTH, 12→11 GAP.

## Sprint 43 (2026-05-25)

- epc/detectors/p1_aggregation.py: type-constancy guard extended to CONFIRMATION tier (CV threshold 0.01); literature anchor: Schelling (1971). C-p1-class-b-lattice2d-fp CLOSED.
- epc/phase2a/structured.py: 2 lattice_2d_continuous supplements added (`smooth_random_field`, `sinusoidal_traveling_wave`). C-lattice_2d_continuous-substrate-undercount CLOSED.
- epc/phase2a/failed_regimes/p3_gray_scott.py: P3 Class C declared N/A (non-Turing regimes rejected at field_std prerequisite).
- epc/phase2a/detector_invariance.py: P3 reclassified as `time_shuffle_invariant=True` (spatial FFT per frame; temporal order irrelevant).
- analysis/run_phase2a_panel.py: `make_p3_detector_fn()` added with `stability_stride=5`; `run_p3()` added.
- analysis/outputs/p1_phase2a_panel.json: re-run Sprint 43; TNR=0.704 (cat TNR 1.000 ↑ from 0.571), Cohen's d=1.624, PARTIAL.
- analysis/outputs/p3_phase2a_panel.json: NEW; TNR=1.000, Cohen's d=+∞, PASS. P3 advances to AT-DEPTH.
- §3.5: appended P1 type-constancy guard extension paragraph (Sprint 43, Schelling 1971).
- §4.6 P1 Schelling: appended Phase-2a panel re-run paragraph (Sprint 43 PARTIAL: cat 1.000, overall 0.704).
- §4.13 P3 Gray-Scott: appended Phase-2a panel result paragraph (Sprint 43 PASS: TNR=1.000, d=+∞; P3 AT-DEPTH).
- §6: updated AT-DEPTH count 6→7/19; Sprint 43 narrative added.
- docs/depth_gap.md: P3 row dim4 PARTIAL→PASS, grade GAP→AT-DEPTH; P1 notes updated (cat TNR ↑); aggregate counts updated 6→7 AT-DEPTH, 13→12 GAP.

## Sprint 42

- §4.21 added: P1 Schelling Phase-2a panel v1.2 PARTIAL (TNR=0.593, Cohen's d=1.298);
  4 carry-forwards opened; P3 Gray-Scott panel paused (lattice_2d_continuous undercount).
- §6 appended: Sprint 42 AT-DEPTH count unchanged 6/19; P1 PARTIAL finding noted.

## Sprint 45 (2026-05-25)

- §4.23 added: P11 Lotka-Volterra Phase-2a panel v1.2 PASS (TNR=1.000, Cohen's d=+∞);
  all 5 positives DEFINITIVE (0.900); P11 dim4 PARTIAL→PASS; grade remains GAP (dim1).
- §6: Sprint 45 AT-DEPTH count unchanged 8/19; P11 dim4 PASS finding noted.
- analysis/outputs/p11_phase2a_panel.json: new file (Sprint 45).

## Sprint 46 (2026-05-25)

- §4.24 added: Sprint 46 continuous_2d dim4 batch panel results.
  - P5 Vicsek: PASS-with-weakness (TNR=0.957, Cohen's d=4.987); all 5 positives DEFINITIVE;
    `time_shuffled` FP at DEFINITIVE (φ is per-frame temporal-order invariant); P5 dim4 PARTIAL→PASS;
    all 4 dims PASS → P5 advances to AT-DEPTH.
  - P2 ABP/MIPS: PASS (TNR=0.958, Cohen's d=3.401); 3/5 positives DEFINITIVE; `permutation_shuffled`
    FP at SCREENING (two_phase_score is spatial-distribution invariant; flag not auto-set per brief);
    P2 dim4 PARTIAL→PASS; grade remains GAP.
  - P6 D'Orsogna/milling: PASS (TNR=0.958, Cohen's d=5.087); all 5 DEFINITIVE; `time_shuffled` FP
    at DEFINITIVE (|L| per frame is temporal-order invariant); P6 dim4 PARTIAL→PASS; grade remains GAP.
- §6: Sprint 46 AT-DEPTH count 9/19 (+1: P5 advances to AT-DEPTH).
- analysis/outputs/p5_phase2a_panel.json: new file (Sprint 46).
- analysis/outputs/p2_phase2a_panel.json: new file (Sprint 46).
- analysis/outputs/p6_phase2a_panel.json: new file (Sprint 46).

## Sprint 47 (2026-05-25)

- §4.25 added: Sprint 47 lattice_1d + oscillator dim4 batch panel results.
  - P8 Nagel-Schreckenberg: PARTIAL (TNR=0.652, Cohen's d=1.751); 2 Class A FPs
    (permutation_shuffled, time_shuffled — stopped_fraction is time-average-invariant;
    carry-forwards C-p8-perm-shuffled-fp, C-p8-time-shuffle-fp); 6 Class C near-onset
    FPs (rho ≥ 0.1167 overlaps jamming transition; carry-forward C-p8-class-c-near-onset);
    P8 dim4 remains PARTIAL; escalate.
  - P10 chimera states: PASS (TNR=0.957, Cohen's d=9.679); all 5 positives DEFINITIVE;
    1 Class A FP (permutation_shuffled at SCREENING; carry-forward C-p10-perm-shuffled-fp);
    Class C (10 ordinary Kuramoto above K_c) all rejected; P10 dim4 PARTIAL→PASS;
    all 4 dims PASS → P10 advances to AT-DEPTH. AT-DEPTH count: 10/19.
- §6: Sprint 47 AT-DEPTH count 10/19 (+1: P10 advances to AT-DEPTH).
- analysis/outputs/p8_phase2a_panel.json: new file (Sprint 47).
- analysis/outputs/p10_phase2a_panel.json: new file (Sprint 47).

## Sprint 48 (2026-05-25)
- §4.26: P21 (Hegselmann-Krause polarization) Phase-2a panel v1.2 PARTIAL (TNR=0.913, d=4.543):
    all 5 positives DEFINITIVE; Class A 2 FPs (permutation_shuffled at CONFIRMATION,
    time_shuffled at SCREENING); Class B+C TNR=1.000. `opinions` detector_format added
    to panel harness. P21 dim4 remains PARTIAL; AT-DEPTH count unchanged: 10/19.
- §6: Sprint 48 AT-DEPTH count 10/19 (unchanged; P21 dim4 PARTIAL, escalate).
- analysis/outputs/p21_phase2a_panel.json: new file (Sprint 48).

## Sprint 49 (2026-05-25)
- §3.5 (new): Invariance-flag batch update rationale (Sprint 49). Documents mathematical grounding
    for six new flag assignments (P1 time_shuffle, P2 perm, P5 time_shuffle, P6 perm+time_shuffle,
    P8 perm+time_shuffle, P21 perm) and P10 exception (adapter artifact vs invariance).
- §4.27 (new): Sprint 49 panel re-run results for P1, P2, P5, P6, P8, P21:
    - P1: TNR 0.704→0.731, syn 0.800→0.889; `time_shuffled` SKIPPED; PARTIAL (C-p1-linear-gradient-fp + C-p1-class-c-subthreshold-fp persist).
    - P2: TNR 0.958→1.000, syn 0.900→1.000; `permutation_shuffled` SKIPPED; clean PASS. C-p2-perm-shuffled-fp CLOSED.
    - P5: TNR 0.957→1.000, syn 0.889→1.000; both substrates SKIPPED; clean PASS (was PASS-with-weakness). C-p5-time-shuffle-fp CLOSED.
    - P6: TNR 0.958→1.000, syn 0.900→1.000; both substrates SKIPPED; clean PASS. C-p6-time-shuffle-fp CLOSED.
    - P8: TNR 0.652→0.714, syn 0.800→1.000; both substrates SKIPPED; PARTIAL (C-p8-class-c-near-onset persists). C-p8-perm-shuffled-fp + C-p8-time-shuffle-fp CLOSED.
    - P21: TNR 0.913→0.955, syn 0.800→0.889; `permutation_shuffled` SKIPPED; PASS-with-weakness (advances from PARTIAL). C-p21-perm-shuffled-fp CLOSED; C-p21-time-shuffled-fp remains.
- §6: Sprint 49 AT-DEPTH count 10/19 (unchanged; P21 dim4 PARTIAL→PASS but dims 1–3 PARTIAL).
- analysis/outputs/p1_phase2a_panel.json: updated (Sprint 49 re-run).
- analysis/outputs/p2_phase2a_panel.json: updated (Sprint 49 re-run).
- analysis/outputs/p5_phase2a_panel.json: updated (Sprint 49 re-run).
- analysis/outputs/p6_phase2a_panel.json: updated (Sprint 49 re-run).
- analysis/outputs/p8_phase2a_panel.json: updated (Sprint 49 re-run).
- analysis/outputs/p21_phase2a_panel.json: updated (Sprint 49 re-run).

## Sprint 66 (2026-06-07)
- §4.21: P7 (lane formation) Phase-2a panel v1.2 PASS-with-weakness added:
    overall TNR=0.955, Cohen's d=6.932. Content prerequisite (counterflow requires
    ≥2 opposing populations, Helbing & Molnár 1995). Class C: 8 weak-repulsion
    (A ∈ [0.1, 0.8]) + 2 single-population regimes. Catalog mates P2/P5/P6
    rejected at prerequisite (no labels). C-p7-time-shuffled-fp carry-forward opened.
- §6: Sprint 66 AT-DEPTH count 19/20 (+1: P7 advances to AT-DEPTH).
- analysis/outputs/p7_phase2a_panel.json: new file (Sprint 66).
