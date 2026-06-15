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

- **Effect size is reported honestly.** `cohens_d` over the discrete
  tier-confidence score returns **NaN** for degenerate constant-score
  comparisons instead of ±∞ (the panel labels these `TNR-PASS-EFFECT-UNDEFINED`).
  The reported effect size is instead computed over each pattern's **continuous
  canonical discriminating metric** — the scalar its screening gate keys on (the
  literature's statistic for the phenomenon) — as Cohen's *d* of the canonical
  positive vs. the pooled negative panel, oriented so a positive value means
  correct separation. Discrimination is reported by both the **true-negative
  rate (TNR)** against three negative classes — Class A synthetic nulls, Class B
  catalog-mate look-alikes, Class C failed regimes — and this continuous *d*
  (full per-pattern table in §4A.3). The discredited ±∞-over-the-discrete-score
  is never reported.

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

A full regression sweep of all 31 Phase-2a discrimination panels (P1–P30 and
P32; there is no separate P31 panel) at `validation-rebuild` HEAD `352f5e8`
returns **31 / 31 panels with true-negative rate (TNR) = 1.0** against all three
negative classes — Class A synthetic nulls, Class B catalog-mate look-alikes,
Class C failed regimes — with **zero errors**. Every canonical positive is
separated from every negative; **no false positives remain anywhere in the
catalog.**

**The re-audit is complete.** Nine patterns that the original validation pass had
not re-examined at the rebuild's depth (P2, P3, P5, P10, P11, P14, P16, P25, P29)
were re-audited under the same standard. Seven needed a faithful detector fix,
and in every case the root cause was the same one that dominated the foundational
rebuild — a discriminating statistic sitting decoratively at the confirmation
tier instead of in the screening gate, where the panel actually counts it (the
panel scores any detection, including the weakest screening tier, as a false
positive). The fixes: P2 (density-dependent motility, the MIPS signature), P5
(failed regimes regenerated as genuinely disordered, so the test is real), P10
(a degenerate chimera autocorrelation that returned 1.0 on flat input), P11
(instrumented and verified faithful), P14 (MLE-vs-log-binned slope consistency,
the Clauset signature that separates SOC from a dissipative cutoff), and P29
(trail temporal-reinforcement). P16 and P25 were re-verified as already faithful.
P21, reported as “PASS-with-weakness” in the prior draft, was likewise corrected
(emergence required at the screening gate) and is now TNR = 1.0.

**P3 (Turing wavelength) warrants a specific note** because its original metric
was not merely lenient but *wrong*: peak-to-mean of the radial FFT answers only
“is there a dominant wavelength?”, which any imposed sinusoid, oriented stripe,
or travelling wave satisfies — a static sinusoid reached DEFINITIVE. The faithful
detector now requires the two signatures that *define* a Turing pattern, both
gated at screening: **stationarity** (a Turing pattern is a standing wave;
spirals, travelling and target waves are rejected) and **spectral isotropy** (a
diffusion-driven instability selects a wave *number* |k\*| — an isotropic ring,
or finite band — not a single wave *vector*; imposed periodicity is rejected).
The adversarial negatives that expose this had been swapped out for degenerate
ones; they are restored as a populated Class C of five non-Turing periodic
fields (static sinusoid, diagonal stripes, rotating spiral, target waves,
travelling plane wave), all now rejected.

Effect size is reported over each pattern's **continuous canonical discriminating
metric** — the scalar its screening gate keys on (the literature's statistic for
the phenomenon) — as Cohen's *d* of the canonical positive (5 seeds) vs. the
pooled negative panel, oriented so a positive value means correct separation.
**All 31 patterns carry a finite, correctly-oriented continuous effect size**
(range 0.5–53.3, median 3.2); none is saturated or undefined. This retires the
discrete tier-confidence proxy that the earlier drafts reported as undefined for
23 patterns.

| Pattern | discriminating metric | Cohen's *d* |
|---|---|---|
| P1 — similarity-driven aggregation | Moran's I | 36.7 |
| P2 — activity-induced phase separation (MIPS) | two-phase coexistence score | 2.6 |
| P3 — Turing wavelength | Turing discriminant ¹ | 5.3 |
| P4 — territoriality | scent-avoidance ratio | 11.7 |
| P5 — flocking (alignment) | polar order φ | 3.3 |
| P6 — milling / vortex | angular momentum \|L\| | 42.2 |
| P7 — lane formation | lane order parameter | 4.7 |
| P8 — traffic jamming | stopped fraction | 1.0 |
| P9 — synchronization | Kuramoto r | 3.9 |
| P10 — chimera state | position–velocity autocorr. | 4.3 |
| P11 — predator–prey oscillation | quarter-period anti-corr. ρ | 53.3 |
| P12 — cyclic dominance | intransitivity score | 9.7 |
| P13 — excitable spiral / target waves | wavefront-speed CV | 1.9 |
| P14 — self-organized criticality | slope consistency ¹ | 0.8 |
| P15 — persistent propagating computation | distinct outcomes | 0.9 |
| P16 — associative memory | completion accuracy | 7.7 |
| P17 — collective sensing | chemotactic index | 7.3 |
| P18 — consensus / coarsening | early-Moran Spearman | 2.5 |
| P19 — emergent leadership | group accuracy | 2.5 |
| P20 — quorum sensing | step-response R² | 2.4 |
| P21 — polarization | between-camp distance | 4.7 |
| P22 — information cascade | Moran's I over time | 4.3 |
| P23 — anti-coordination | scaled variance (Fano) | 0.5 |
| P24 — homeostatic regulation | deviation integral | 6.0 |
| P25 — equifinality | convergence variance ratio | 1.3 |
| P26 — stochastic resonance | peak performance | 1.8 |
| P27 — spatial reciprocity | cooperation fraction | 2.6 |
| P28 — wealth condensation | final Gini | 3.2 |
| P29 — trail / network formation | weight–distance corr. | 1.7 |
| P30 — autopoiesis (boundary formation) | membrane radial CV | 2.2 |
| P32 — emergent specialization | role-entropy decline | 1.6 |

¹ Derived gating scalar: P3 = min(angular_entropy/0.55, field_stationarity/0.95),
the binding margin of the 2-D Turing gate (a single FFT scalar is insufficient —
peak-to-mean alone gives *d* = −0.15, since the adversaries are equally peaked);
P14 = \|τ_MLE − τ_logbin\|, the SOC slope-consistency the gate keys on. For
lower-is-positive metrics (P4, P11, P13, P14, P23, P25, P30) *d* is oriented so a
positive value denotes correct separation. The two weakest (P23 = 0.5, P14 =
0.8) reflect genuinely subtle single-scalar effects — sub-Poissonian variance
suppression and the SOC slope gap over a diverse avalanche-null panel — and still
discriminate at TNR = 1.0.

How each pattern reached TNR = 1.0:

| How fixed | Patterns |
|---|---|
| Faithfully rebuilt: detector recomputes the canonical metric from the substrate; positives self-organize from a disordered start; look-alikes + nulls + failed regimes all rejected | P3, P4, P6, P7, P8, P9, P13, P15, P17, P18, P19, P20, P21, P23, P24, P26, P28, P30, P32 |
| Re-audited this pass: discriminator moved into the screening gate / made faithful | P2, P5, P10, P11, P14, P29 |
| Re-verified, already faithful | P16, P25 |
| Tier-reachability corrected (null-floor fix made the intended tier reachable) | P1, P12, P22, P27 |

The panel still emits the verdict label `TNR-PASS-EFFECT-UNDEFINED` for the
patterns whose *discrete tier score* is degenerate, but this is now a property of
the discrete proxy only — the **continuous effect size above is the reported
*d***, and it is finite and correctly oriented for every pattern. Neither the
inflated nor the infinite *d* of the earlier drafts is used anywhere.

## 4A.4 Honest limitations

- **Effect size — RESOLVED.** Previously the principal remaining task: the
  reported *d* was over the discrete tier-confidence score, undefined at perfect
  separation for 23 patterns. A per-pattern registry of each detector's
  *continuous* canonical discriminating metric (`epc/phase2a/continuous_metrics.py`)
  now drives a faithful continuous *d* (positives vs. pooled negatives, oriented),
  finite for all 31 patterns (§4A.3). The discrete `cohens_d_positive_vs_panel`
  is retained alongside for provenance; the continuous *d* is the reported value.

- **P4 (territoriality) — RESOLVED.** Initially deferred: the cumulative
  exclusivity / overlap / occupancy-scent-correlation metrics cannot separate
  scent-mediated territoriality from a plain random walk. The faithful detector
  uses a scent-mediated **movement-causality** test on a per-step substrate
  (`avoidance_ratio`: foreign-scent at the chosen cell vs. the mean over
  available neighbours), which scores ≪1 for territorial movement and ≈1 for a
  scent-blind walk on identical scent dynamics. TNR = 1.0, continuous *d* = 12.8.

- **P30 (autopoiesis) — RESOLVED.** Initially deferred: the model produced a
  diffuse cloud (links at ~1.5–2× the equilibrium radius, radial spread
  indistinguishable from random scatter) and the “closure” metric returned 1.0
  on that scatter. A model-level fix (a link-production cap so a thin shell can
  form) plus a radial-tightness detector now yields a genuine self-maintained
  membrane that is distinguished from a random cloud and from mechanism-off
  regimes. TNR = 1.0.

- **P27 (spatial reciprocity).** Passes on TNR with a continuous *d* = 2.8, but
  the canonical positive is seed-fragile (2 of 5 seeds undetected). The
  discrimination is sound; the positive regime needs strengthening for a robust
  depth claim. P2 and P29 carry a milder version of the same seed-fragility note.
  These are positive-robustness calibrations, not false positives.

- **Honest depth caps.** P13 and P31 honestly top out at CONFIRMATION on the
  single-run detector (the null-floor correction removed unreachable DEFINITIVE
  gates). P31's DEFINITIVE evidence is the separate multi-run non-redundancy test
  (per-agent delayed-gratification carries variance beyond the algorithm label;
  ΔR² = +0.018…+0.033, vanishing under timing-shuffle).

## 4A.5 Reproducibility

All detector code, models, failed-regime substrates, panels, and the standalone
discrimination tests are on the `validation-rebuild` branch. Per-pattern fix
history with before/after metrics is in
`docs/validation_rebuild/foundational_defects_verified.md`. The P28 mechanism
discrimination + metadata-flip-invariance test is
`analysis/validation_rebuild/p28_mechanism_discrimination.py`; the P31
non-redundancy test is `analysis/validation_rebuild/p31_non_redundancy_test.py`.
