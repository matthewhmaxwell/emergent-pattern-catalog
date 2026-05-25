# Paper Draft Changelog

Per-sprint mechanical changes to docs/paper_section*_draft.md files. Updates
applied by the orchestrator chain. Voice/framing changes flagged for chat-led
review at the next paper-review checkpoint.

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
