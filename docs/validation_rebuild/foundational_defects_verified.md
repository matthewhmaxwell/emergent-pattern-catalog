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
Then Group B (behavioral), then Group C (several likely demote/pull).
