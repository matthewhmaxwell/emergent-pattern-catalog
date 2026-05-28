# Paper Draft Changelog

Per-sprint mechanical changes to docs/paper_section*_draft.md files. Updates
applied by the orchestrator chain. Voice/framing changes flagged for chat-led
review at the next paper-review checkpoint.

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
