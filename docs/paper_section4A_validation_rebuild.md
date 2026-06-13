# Section 4A: Validation — Rebuild and Honest Status

> **This section supersedes the per-pattern PASS claims and all reported
> Cohen's *d* values in the earlier drafts of §4 and §6.** Those drafts were
> written during an initial validation pass that a subsequent audit found to be
> substantially unsound. The findings below are the corrected, faithful state.

## 4A.1 Why a rebuild

After the 32-pattern catalog was implemented, an independent multi-agent code
audit of the detector suite and its validation harness returned a large number
of confirmed findings and **zero patterns without defects**. Two classes of
problem dominated:

1. **A vacuous effect-size gate.** The Phase-2a discrimination panel reported a
   Cohen's *d* computed over the **discrete tier-confidence score**, not over
   the continuous canonical metric. Because a correctly-discriminating detector
   assigns a constant high confidence to positives and constant zero to rejected
   negatives, the within-group variance is ~0 and the standardized effect size
   is undefined; the harness reported this as **±∞**, which then auto-passed the
   effect-size gate. **15 of 32 patterns reported *d* = ∞.** Every numeric *d*
   in the earlier §4/§6 drafts (e.g. "*d* = 8.282", "*d* = 10.585", "*d* = +∞")
   originates from this artifact and should be disregarded.

2. **Detectors that did not faithfully compute their stated metric.** Examples
   the audit confirmed and we reproduced: a detector that re-simulated from
   metadata instead of reading the substrate (P17); tier gates keyed on
   **self-reported metadata booleans** that a single flag-flip could change
   (P28); a "closure" metric that returned 1.0 on random scatter (P30); a
   "sharpness" discriminator that never ran because an amplitude guard rejected
   negatives first (P20); panels whose **PASS depended on substrate
   reclassification or a hardcoded positive driver** (P18, P6); and
   significance gates set **below the permutation null's resolution floor**, so
   the strongest tier was mathematically unreachable (P1, P12, P13, P31; also
   P18/P22/P27 at the panel level).

The science of each pattern is not in question — all 32 are canonical,
literature-validated emergent phenomena. The defects were in **our detector
code and validation harness**. The rebuild fixes them to the standard "faithful
to the established science," not "passes a panel."

## 4A.2 What changed in the methodology

- **Effect size is reported honestly.** `cohens_d` now returns **NaN** for
  degenerate constant-score comparisons instead of ±∞, and the panel emits a
  distinct verdict, **`TNR-PASS-EFFECT-UNDEFINED`**, when discrimination is
  perfect on a discrete tier score. Discrimination is therefore reported by the
  **true-negative rate (TNR)** against three negative classes — Class A
  synthetic nulls, Class B catalog-mate look-alikes, Class C failed regimes —
  not by a standardized effect size. A continuous-metric effect size is reported
  only where the metric is continuous (e.g. P32 *d* = 2.64; P27 *d* = 2.78).

- **Detection counts at the screening tier.** The panel scores a negative as a
  true negative only if the detector returns *not detected*. A look-alike that
  fires even at the weakest (screening) tier is a false positive. This forces
  the discriminating signal into the **screening gate** rather than letting it
  sit decoratively at confirmation — the change that drove the faithful fixes to
  P20 (sharpness), P28 (resource conservation), and P18 (domain coarsening).

- **Mechanisms are derived from the trajectory, not asserted.** Tier gates that
  read self-reported metadata booleans were replaced with quantities computed
  from the substrate (P28 conservation + emergence + monotonic growth; P24
  measured restoring feedback; P26 noise-driven coherence; P32 cross-condition
  demand efficiency).

- **Tier gates respect the null's resolution.** A permutation/surrogate *p* from
  *N* draws floors at 1/(*N*+1). The strongest tier now requires *p* ≤ that
  floor; *N* is chosen so the floor is ≤ 0.005 (≥ 199 draws), or ≤ 0.001 where a
  stronger claim is made (≥ 999 draws, e.g. P1, P32).

- **Canonical positives must show the phenomenon emerge.** Positives that began
  pre-ordered (P6 hand-planted mill, P7 pre-segregated lanes, P28 frame-1 already
  condensed) were replaced with runs that **self-organize from a disordered
  start**, so the detector is tested on emergence rather than on a frozen
  end-state.

## 4A.3 Honest results

Discrimination panels run three negative classes and report TNR. With the
corrections above, the verified state is:

| Status | Patterns | Notes |
|---|---|---|
| **Faithfully rebuilt + verified, TNR = 1.0** | P6, P7, P8, P9, P13, P15, P17, P18, P19, P20, P21, P23, P24, P26, P28, P32 | Detector recomputes the canonical metric from the substrate; positives self-organize; look-alikes + nulls + failed regimes all rejected. |
| **Tier-reachability corrected, TNR = 1.0** | P1, P12, P22, P27 | Null-floor fix made the intended tier reachable; P1/P12/P22 positives reach DEFINITIVE. |
| **Validated, reduced claim** | P21 (PASS-with-weakness: one synthetic static-bimodal null persists) | Reported honestly, not papered over. |
| **Carried forward from the original validation** | P2, P3, P5, P10, P11, P14, P16, P25, P29 | Panels pass on TNR; these were **not** re-audited at the rebuild's depth and the honest-effect-size caveat (§4A.1) applies to any *d* previously reported for them. |
| **Validated by a separate reproducible test** | P31 | The single-run detector honestly caps at CONFIRMATION; the DEFINITIVE evidence is the multi-run non-redundancy test (per-agent delayed-gratification carries variance beyond the algorithm label; ΔR² = +0.018…+0.033, vanishing under timing-shuffle). |
| **Honest depth = CONFIRMATION (not DEFINITIVE)** | P13, P31 | The null-floor correction removed unreachable DEFINITIVE gates; these patterns honestly top out at CONFIRMATION on the single-run detector. |
| **Deferred — genuine limitation, not a pass** | **P4**, **P30** | See §4A.4. Their earlier "PASS" was an artifact of the broken metrics and is **withdrawn**. |

The dominant verdict, `TNR-PASS-EFFECT-UNDEFINED`, should be read as: *the
detector perfectly separated the canonical positive from every synthetic null,
catalog look-alike, and failed regime, and a standardized effect size is
undefined because the tier score is discrete* — not as the inflated *d* the
earlier drafts reported.

## 4A.4 Honest limitations

- **P4 (territoriality).** The exclusivity / overlap / occupancy-scent-correlation
  metrics cannot separate scent-mediated territoriality from a plain random walk
  (a random walk scores equally or more exclusive; the scent correlation is
  noise). A faithful detector needs a scent-mediated **movement-causality** test
  on a finer per-step substrate. Deferred; the pattern is real, the current
  detector is not faithful, so no depth grade is claimed.

- **P30 (autopoiesis).** Across a parameter sweep (single catalyst, radial-spring
  strength up to 8, small *dt*, small box, balanced production/decay), the model
  does **not** form a thin membrane — links sit at ~1.5–2× the equilibrium radius
  with radial spread indistinguishable from random scatter. The detector cannot
  faithfully read a membrane the model does not produce; a model-level revision
  (capped production, dominant radial confinement) is required first. Deferred;
  the earlier "PASS" (broken closure metric) is withdrawn.

- **P27 (spatial reciprocity).** Passes on TNR, but the canonical positive is
  seed-fragile (2 of 5 seeds undetected). The discrimination is sound; the
  positive regime needs strengthening for a robust depth claim.

- **Effect size.** Where the tier score is discrete and discrimination is perfect,
  the standardized effect size is genuinely **undefined**, and we report it as
  such rather than as a large or infinite number.

## 4A.5 Reproducibility

All detector code, models, failed-regime substrates, panels, and the standalone
discrimination tests are on the `validation-rebuild` branch. Per-pattern fix
history with before/after metrics is in
`docs/validation_rebuild/foundational_defects_verified.md`. The P28 mechanism
discrimination + metadata-flip-invariance test is
`analysis/validation_rebuild/p28_mechanism_discrimination.py`; the P31
non-redundancy test is `analysis/validation_rebuild/p31_non_redundancy_test.py`.
