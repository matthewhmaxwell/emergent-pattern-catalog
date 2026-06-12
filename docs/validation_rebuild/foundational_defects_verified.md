# Foundational Detector Defects — Re-Verified (19/19 CONFIRMED)

**Date:** 2026-06-10 · **Branch:** `validation-rebuild` · **Method:** 4 independent
read-only agents re-checked each prior-audit claim against shipped code + live
behavioral probes on the VPS. **All 19 confirmed; none refuted.**

This is the *foundational* layer — does the detector faithfully detect its
canonical phenomenon. It is prior to (and blocks) the negative-panel / effect-size
QA. A pattern here is broken at the detector level, not just the panel level.

## The firmed-up list

### Group A — detector computes the WRONG metric (named canonical metric absent)

| P | spec's canonical metric | what the code actually computes | fix direction |
|--|--|--|--|
| **P8** | backward density-wave speed / fundamental-diagram breakdown / spacetime bands | `stopped_fraction` P(v=0) scalar + per-car jam-lifetime tail | add a spatial backward-propagation / fundamental-diagram metric as primary |
| **P13** | Transfer Entropy ≈0 across wave collisions (vs P15) | wavefront-speed CV + a **hardcoded model-name string lookup** for the P15 boundary (`"Full TE deferred"`) | **wire in the already-built `P13P15Discriminator` TE** (exists, never called by detector/panel) |
| **P15** | Transfer Entropy >0 across collisions / signal routing | GoL-rule replay reproducibility + outcome diversity (no TE module imported) | compute boundary-conditioned TE across collisions (TE machinery already in repo) |
| **P21** | persistent multimodality / between-cluster distance / final variance | "Hartigan's dip" that is really a KS uniformity distance (**inverted**: fires on Gaussians, passes on Uniform); real discrimination done by a sorted-gap `n_clusters>=2` gate | implement a genuine multimodality test; stop mislabeling |
| **P31** | per-agent distance-based delayed-gratification vs random-walk null | a global monotonicity-backtracking scalar; PASS only obtained by **removing the algorithm-identity control** (then DG re-encodes the algo label); headline ΔR²=+0.645 **does not reproduce** (sklearn not even installed) — honest result is FAIL | per-agent DG index vs random-walk null, keep the algo control (fails honestly) |

### Group B — detector FIRES ON INPUT IT SHOULDN'T (behavioral, probe-confirmed)

| P | probe result | fix direction |
|--|--|--|
| **P17** | **ignores its `history` argument entirely** — re-runs the model from metadata; garbage vs real history give **bit-identical** output; fires on any substrate | compute the CI/N-scaling signature from the supplied `history`, not a metadata re-sim |
| **P9** | fires on a **frozen constant phase field** (zero dynamics/coupling): `detected=True, r_mean=1.0` | gate on a coupling/entrainment signature, not the final-state order-parameter alone |
| **P24** | fires **DEFINITIVE** on 3 feedback-free lookalikes (passive relaxation to wrong target, perturbation-independent sine, step-and-hold) | require regulated return to set-point + measured perturbation→response coupling; drop the metadata flag |
| **P26** | fires **DEFINITIVE** on a non-SR deterministic tent with constant noise; the one guard fails open | test the SR mechanism (noise sweep, barrier-crossing/SNR peak) measured from input |
| **P4** | a **plain random walk** → DEFINITIVE conf 0.90 (indistinguishable from territoriality); only a metadata label demotes it. dim2 is really **15/20**, not the claimed 20/20 | discriminating scent-mediated-exclusion metric, no metadata label peek |
| **P19** | a **leaderless flock** (random informed-mask, declared dir = flock's own heading) → 20/20 detected | screen on influence-asymmetry (informed-vs-naive early alignment / TE), not group accuracy |
| **P23** | variance-below-baseline **alone** fires CONFIRMATION; canonical anti-persistence not enforced (detector's own comment claims both) | require the negative-autocorrelation AND-clause at confirmation |

### Group C — reproduction / manufactured-panel / tier defects

| P | finding | fix direction |
|--|--|--|
| **P6** | committed substrate config (ring init, 200 steps) → **tier=none**; DEFINITIVE only from the panel builder's hardcoded 3000-step driver; positive is a planted ring, not emergent | make the committed config the real positive; gate emergence on random-IC self-organization |
| **P7** | metric = end-state label-position correlation, not lane *formation*; dim1 `phi_init`=0.498 is ~3× the shipped model's real ~0.17 | add a dynamics/spontaneity prerequisite; regenerate dim1 from the shipped model |
| **P18** | phantom final-state "consensus" metric (real metric is time-ordered Spearman); PASS manufactured by reclassifying substrate `lattice→network` and swapping the 3 lattice neighbors for a trivially-rejected P21 mate | fix the invariance label/flag; restore lattice Class B; re-run honestly |
| **P20** | `step_r2` sharpness discriminator is **dead code** (computed, never gates); time-shuffled positive still DEFINITIVE; Class B has zero catalog mates | make R²-sharpness actually gate; stop reconstructing hysteresis from static tags |
| **P28** | **no panel exists**; mechanism gate reads 4 self-reported metadata booleans; flipping one flag on the same trajectory flips CONFIRMATION→DEFINITIVE | derive mechanism from trajectory (exchange conservation/symmetry); add Bouchaud-Mézard lookalike; build a real panel |
| **P30** | `closure_fraction` fires on random **non-membrane** scatter (=1.0); `association_score` is mere co-location; none of self-production / topological closure / gradient / self-repair implemented | measure real self-production + closed-boundary topology + maintained gradient + self-repair |
| **P32** | **DEFINITIVE mathematically unreachable** (null-p floor 1/200=0.005, gate requires `<0.005`); efficiency = within-run late-vs-early delta, not cross-condition | fix permutation floor or gate; compute efficiency vs a non-specialized baseline run |

## The deeper pattern — and what it does NOT mean

Several of these are **not "a detector with a bug"** but **"a detector that never
worked, with validation constructed to pass"**: P31 (PASS required deleting the
control), P18 (panel PASS via substrate reclassification), P6 (positive only via a
separate hardcoded driver), P32 (DEFINITIVE unreachable by construction), P17
(ignores input).

This is a statement about **our implementation**, NOT about the science. Every one
of the 32 is a **canonical, literature-validated emergent phenomenon** (Schelling
segregation, Vicsek flocking, SIR, Lotka–Volterra, RPS cycling, Kuramoto sync, BTW
sandpile SOC, Hopfield, D'Orsogna milling, the voter model, collective sensing, …)
with established detection criteria and, in most cases, reference implementations
available (literature + GitHub). **None of these get demoted or pulled.** The
phenomenon is real and detectable; our detector cheated or measured the wrong
thing. The fix for each is to implement the literature's canonical detection metric
faithfully, on the actual substrate, using the published method + reference code as
ground truth — then confirm it fires on the real positive and rejects genuine
lookalikes. The bar is **"faithful to the established science,"** not "passes," and
never "pull if hard."

## Fix order (user: easiest/highest-confidence first)

1. **P13** — the canonical TE discriminator already exists in-repo; just wire it in. Highest-confidence, smallest change.
2. **P15** — same TE machinery; compute TE across collisions.
3. **P8** — add the spatial backward-wave / fundamental-diagram metric.
Then Group B (behavioral), then Group C. None are pulled — each is fixed to the
literature's canonical metric, with tier/depth labels corrected to what the
detector honestly reaches.

## Fix log

### P13 — RESOLVED (commit 7e6c235, validation-rebuild)
Four layered defects, all fixed:
1. **Fake metric → real:** the P15 exclusion now computes boundary-conditioned
   **transfer entropy on the substrate** (`P13P15Discriminator`, was a
   `model_name in (...)` string lookup with a "TE deferred" comment). TE not above
   a spatial-shuffle null (p≈1.0, shuffle-invariant) ⇒ excitable waves ⇒ P13;
   above ⇒ P15. Verified: Greenberg-Hastings clears as P13, Game-of-Life does not.
2. **Cosmetic → gating:** the exclusion was computed *after* the tier was fixed and
   never demoted it. DEFINITIVE now genuinely requires `classification=="P13"`.
3. **Positive couldn't reach its tier:** canonical positive lengthened 300→900
   steps so it accrues ≥50 spiral rotations and reaches CONFIRMATION (legitimate —
   GH spirals genuinely persist; not threshold-gaming).
4. **Panel null under-resolved:** `run_p13` overrode the detector's 199-run null
   with 99, flooring `null_p` at 0.01 so CONFIRMATION (`<0.01`) was unreachable in
   the panel. Restored to 199.

Honest outcome: P13 reaches **CONFIRMATION** with the real TE discriminator
operative; panel TNR=1.0 (all 7 catalog mates + synthetics rejected). The
committed **DEFINITIVE** claim was never reachable — it requires `null_p<0.001`
but a 199-run null floors at 0.005. **Cross-cutting bug shared with P32:** the
DEFINITIVE `null_p<0.001` gate is unreachable unless the null uses ≥1000 runs;
to be fixed catalog-wide (raise null resolution or relax the gate to the null's
resolution) rather than per-pattern.

### P15 — RESOLVED (commit c4a1cbc, validation-rebuild)
Same TE machinery as P13, opposite condition. The P13 boundary in P15's
`_check_exclusions` was a bare `"inconclusive"` placeholder; the named canonical
metric (Transfer Entropy across collisions, >0 for P15) was never computed.
Fixed: added the `_p13_te` helper (boundary-TE via `P13P15Discriminator`,
substrate-based, memoized), made `_check_exclusions` report the real result, and
gated DEFINITIVE on `classification=="P15_candidate"` (directed information flow).
Verified — GoL canonical positive reaches **DEFINITIVE** (conf 0.85) with the TE
showing P15_candidate (p=0.02, ratio ~13–15×) and P13 excluded; a Greenberg-
Hastings substrate is rejected (detected=False). Panel: positive 0.71, TNR=1.0,
verdict PASS (d=8.28). Unlike P13, P15's DEFINITIVE is reachable (its gate is
rep+diversity+n_variations+TE, not the broken `null_p<0.001`).

### P8 — RESOLVED (commit 1bd91a2, validation-rebuild)
The detector discriminated on `stopped_fraction` = P(v=0) only, which is **0.61 in
gridlock-saturation** → DEFINITIVE false positive (high density ≠ emergent jamming;
the `jam_lifetime` gate doesn't catch it). Added the canonical **phase-coexistence**
metric: `min(stopped_fraction, free_fraction)` where free = P(v=v_max). Emergent
stop-and-go requires a jammed phase AND a free-flow phase at once (Sugiyama 2008 —
the NS jamming transition is a density phase separation); free-flow has no jam
phase, gridlock has no free phase. Gated screening + confirmation on it
(threshold 0.03). Calibrated against the real distribution: canonical positive
(ρ=0.30) coexistence 0.067–0.074 across 5 seeds vs saturation ≤0.013 (ρ≥0.45) — a
clean ~5× gap. Verified — jam ρ=0.30 → **DEFINITIVE**, free-flow + saturation
(ρ=0.55, 0.70) **rejected**. Panel: positive 0.9, catalog_tnr=1.0,
**failed_regime_tnr=1.0** (saturation negatives now rejected), overall_tnr=1.0.

### P17 — RESOLVED (commit 72c9d35) — behavioral cluster
Detector **ignored its `history` argument entirely** — re-ran CollectiveSensingModel
from metadata, firing on any substrate (bit-identical output for garbage vs real
input). Rewrote `detect()` to work from the substrate via three measured signatures:
**chemotactic index** to the true peak (climb), **cohesion** (group dispersion <
0.10 — rejects social_off, which disperses), and **emergence** (individual SNR <
0.5 — rejects field_too_strong, trivially sensable alone). All three at screening,
calibrated vs the real regimes. Decisive test passes: same metadata + different
substrates → different verdicts. Panel: positive 0.78, TNR=1.0, mates rejected at
METRIC stage.

### P9 — RESOLVED (commit 80275cb) — behavioral cluster
Fired on a **frozen constant phase field** (r=1) because the gate was `r_mean>0.7`
alone. Added the emergence signature: synchronization must arise from an
**incoherent start** — `r_init` (min order parameter over the transient) < 0.4 and
a rise > 0.3. Frozen field (r_init=1) and sub-critical (r_final low) now rejected.
Also lengthened the positive (n_T_osc≈103) and bumped null-runs 99→199 (the
**cross-cutting p-floor bug**, which here blocked CONFIRMATION's p<0.01). Verified:
frozen + sub-critical rejected, positive→confirmation/definitive; panel TNR=1.0.

### P24 — RESOLVED (commit f2d4044) — behavioral cluster
Fired DEFINITIVE on **feedback-free bounded series** (sine, constant offset,
relaxation to a wrong target) because it only tested bounded-vs-linear-ramp and a
self-reported `has_active_feedback` flag. Now requires **measured** homeostasis:
screening needs regulation (ss_deviation < 30% of perturbation amplitude — the
`ss_deviation` that was computed but unused), confirmation needs measured active
feedback (`restoring_slope < 0` from regressing dx/dt on x−setpoint). Dropped the
metadata flag. Verified: positive→DEFINITIVE, all four feedback-free lookalikes
rejected even with `has_active_feedback=True`; panel TNR=1.0.

### P26 — RESOLVED (commit 96dab30) — behavioral cluster
Fired DEFINITIVE on a non-SR inverted-U via a `has_subthreshold_signal` metadata
flag that "fails open" + a shape proxy. Now verifies the SR **mechanism** from the
substrate: confirmation requires **noise-driven** (per-level `std(diff(x))`
Spearman-correlates with the labelled noise sweep > 0.7 — rejects a fake
constant-noise sweep whose inverted-U comes from a ramped signal) AND **measured
subthreshold** (zero-noise response < 30% of peak). Dropped the metadata flag.
Verified: positive noise-corr=1.0, subthreshold-ratio=0.06 → DEFINITIVE; panel TNR=1.0.

### P4 — DEFERRED (task #35) — needs a movement-causality metric, NOT shipped
P4's metrics (exclusivity, pairwise overlap, occ_scent_corr) **fundamentally cannot
separate territoriality from a plain random walk**: a random walk scores *equal or
higher* exclusivity (sparse agents with cumulative occupancy separate regardless of
scent), and `occ_scent_corr` is noise (RW −0.07 vs positive −0.005). Confirmed across
N=4..25. The faithful discriminator is scent-mediated **movement causality** (does an
agent turn away from foreign scent?), which the 100-step cumulative-occupancy
snapshots don't expose. Needs a finer-grained substrate (per-step move vs foreign
scent) or a cleverer causal statistic — a gate-threshold tweak would just relabel the
same fake test. Returning after P19/P23.

### P19 — RESOLVED (commit 80407b9) — behavioral cluster
Screened on group directional accuracy (`cos(group_heading − preferred) > 0.3`),
which a **leaderless flock trivially passes** if handed its own emergent heading.
The correct causal metric (`_compute_influence_asymmetry`: informed-vs-naive pull
with a label-shuffle null) already existed but was only used at confirmation. Moved
it BEFORE screening and gated screening on a **significant** informed-minority pull
(p<0.01); confirmation adds consistent leadership (pull_fraction>0.5). Verified:
leaderless flocks (rho_zero, bias_zero) all rejected; positive PASS; panel TNR=1.0.

### P23 — RESOLVED (commit 094ca57) — behavioral cluster
Confirmation fired on variance-below-baseline ALONE (an `OR`), so a clamped Gaussian
(low variance, no coordination) fired DEFINITIVE. The audit suggested requiring
anti-persistence (`p_ac1<0.01`), but the canonical positives **don't show strong
anti-persistence at this config** (ac1≈−0.02, p_ac1>0.05) — that fix would reject
them. The signature this regime *does* exhibit is **adaptation**: attendance variance
falls over time as agents learn to anti-coordinate (positives 0.54-0.67 late/early;
a stationary clamped Gaussian ~1.0), measured on the full series (the learning
transient is in the burn-in). Confirmation now needs variance_below + p_sv<0.01 +
adaptation_ratio<0.85. Verified: clamped Gaussian rejected, positives→DEFINITIVE;
panel TNR=1.0.

---

## Behavioral cluster summary (7 patterns)
**Fixed (6):** P17, P9, P24, P26, P19, P23 — all panel TNR=1.0.
**Deferred (1):** P4 (task #35) — metrics can't separate territoriality from a random
walk; needs a movement-causality metric.
**Recurring theme:** several "canonical positives" only weakly exhibit their textbook
signature at the committed config (P9 short run, P13 unreachable tier, P23 weak
anti-persistence, P4 too-sparse) — the honest discriminator is often a *different but
real* feature (emergence-from-incoherence, adaptation, etc.) calibrated to what the
substrate actually shows.
