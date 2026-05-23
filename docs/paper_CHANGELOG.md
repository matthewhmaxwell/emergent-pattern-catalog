# Paper Draft Changelog

Per-sprint mechanical changes to docs/paper_section*_draft.md files. Updates
applied by the orchestrator chain. Voice/framing changes flagged for chat-led
review at the next paper-review checkpoint.

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
