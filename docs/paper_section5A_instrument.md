# Section 5A: The Catalog as an Instrument

> This section establishes the catalog's **outward** use: pointing the detector
> battery at systems it was not built on, and asking "what emergent behavior does
> this exhibit, with what calibrated confidence?" It builds on the within-family
> transfer matrix of §5 and the faithful detectors of §4A. All results are at
> `validation-rebuild` HEAD `68f105b`.

## 5A.1 From taxonomy to instrument

A validated catalog of 32 patterns is a taxonomy. The artifact we actually want
is a **measurement instrument**: a detector battery that can be pointed at an
arbitrary system — an independent agent-based model, a multi-agent AI system, or
a real time series — and return a calibrated emergent-pattern profile rather than
a yes/no on a single pre-specified pattern. This section reports the four steps
that turn the taxonomy into that instrument and the validation of each:
a cross-detector **calibration layer** (5A.2), a substrate-agnostic
**generic-emergence / "none-of-the-above"** primitive (5A.3), the battery's
**self-recognition** and **cross-model** behavior (5A.4–5A.5), an
**out-of-distribution validation suite** on held-out systems (5A.6), and an
**external demonstration** on a real multi-agent LLM swarm (5A.7).

The instrument's output for one observation is produced by a single entry point
(`analysis/battery_profile.profile_observation`): it runs every detector on the
observation, extracts each pattern's continuous canonical metric (the scalar of
§4A.3), calibrates it, ranks the patterns, and combines the result with the
generic-emergence indicator into a **three-way verdict**:

- **MATCH** — a specific detector fires through its full multi-criterion gate;
- **EMERGENT-UNCLASSIFIED** — generic emergence is high but no detector fires
  (the novelty / discovery signal);
- **NO-EMERGENCE** — neither a detector nor the emergence indicator trips.

## 5A.2 Calibration: a cross-comparable confidence (T2a)

A detector's native output is a tier (screening / confirmation / definitive) and
a tier-capped confidence, tuned on that pattern's own positives and negatives.
These are not comparable across patterns: a "confirmation" from one detector and
a "confirmation" from another encode different evidence. To rank patterns against
one another for an unknown system, each detector's continuous canonical metric
(§4A.3) is mapped through a **calibrator** built from the rebuild's faithful
reference distributions — the canonical positives and the pooled negative panel.

For an observation's oriented metric value, the calibrator reports where it sits
among genuine positives, a conformal tail probability against the negatives, and
an **AUC-style separation score** in [0, 1] — the probability that the value is
more pattern-like than a random negative. This separation score is the single
cross-detector-comparable instrument confidence. Detectors whose substrate the
observation does not match reject at prerequisite and calibrate as *absent*,
which is the correct battery behavior: only the patterns actually present light
up. (Implementation: `epc/phase2a/calibration.py`, `epc/phase2a/battery.py`.)

## 5A.3 Generic emergence: the none-of-the-above primitive (T2b)

A measurement instrument must be able to say "something self-organizing is
happening here, and it matches none of my 32 detectors." Without this, every
observation is forced onto the nearest catalog label, and the instrument can
neither abstain honestly nor flag a candidate novel pattern. The
generic-emergence indicator (`epc/phase2a/emergence.py`) supplies it.

The principle is the detectors' own null logic lifted to a generic structural
statistic: the late-window observation is scored as emergent when it is both (a)
more structured than an order-destroying shuffle of itself (a null-excess
*z*-score) **and** (b) more structured than its own early window (an
order-gain term). The structure statistic is chosen by data kind — field
autocorrelation (Moran's *I*), point clustering, phase order (Kuramoto *r*),
vector concentration, and an **orientation channel** (polar and nematic order
parameters from velocity/heading data, added during the external demonstration of
5A.7 to capture apolar/nematic order that a clustering measure cannot see). The
orientation channel is combined by taking the stronger of the spatial and
orientation scores; because no non-emergent null carries velocity data, it can
only raise a true-emergent score and introduces no false-novelty.

The indicator covers field/grid autocorrelation, point clustering, and phase
order well (positive-vs-null separation E ≈ 0.75–0.98 against E ≈ 0.02–0.07 for
random nulls). It has **characterized gaps** — morphologies whose late-frame
structure is neither point-clustering, field-autocorrelation, nor phase/orientation
order (rotational milling, banded lanes, transient waves that decay to a uniform
final state, and distributional concentration). On those a low emergence score is
*inconclusive*, not "no emergence," and the caller is expected to treat it as
such. These gaps are reported, not hidden (5A.8).

## 5A.4 Self-recognition: the battery is not a lookup table

The first test of an instrument is that it recognizes each catalog pattern from
its own canonical positive **without firing the wrong detector**. Running every
canonical positive through the full battery yields a confusion matrix in which
**30 of 31 patterns rank their own pattern #1**, **29 of 31 are firm self-matches**
(the detector fires: 23 at definitive, 3 at confirmation, 3 at screening), and
cross-fire is near-nil. The one non-self-match (P13, excitable waves, edged at
screening by P18) and the screening-only self-matches (P22, P27, P29) are exactly
the seed-fragile / honest-depth-capped patterns flagged in §4A.4; they matter for
the operating-point choice in 5A.6. (`docs/validation_rebuild/battery_confusion_matrix.md`.)

## 5A.5 Cross-model recognition: the phenomenon, not the implementation (T1b)

A detector that recognizes its *phenomenon* — rather than its *native
implementation* — should still rank its pattern #1 on an independent second
implementation it was never tuned on. For each pattern with an available
alternative model, a positive was built from the alternative and run through the
full battery. **All 7 rank the correct pattern #1; 6 of 7 are firm MATCHes:** P16
on a Boolean gene-regulatory network (vs. Hopfield), P26 on a threshold unit (vs.
double-well), P29 on a Physarum model (vs. ant trail), P20 on a fraction-threshold
model (vs. autoinducer), P25 on a multi-basin GRN (vs. canalized landscape), and
P23 on El Farol (vs. Minority Game). The seventh, P24 on an integral controller
(vs. proportional), ranks P24 #1 but stays below its proportional-tuned firing
threshold — a genuine, surfaced generalization limit (a controller-agnostic
restoring-feedback signature would close it). (§5;
`docs/validation_rebuild/cross_model_generalization.md`.)

## 5A.6 Out-of-distribution validation (T2c)

Self-recognition and cross-model recognition test the battery on systems drawn
from, or adjacent to, the catalog. The instrument claim requires testing it on
systems the detectors were **not built on**. The OOD suite
(`analysis/t2c_ood_suite.py`) runs the battery, blind, over sixteen held-out
systems in four arms:

- **recognition** — independent implementations of catalog phenomena (the §5A.5
  set), which should MATCH the correct pattern;
- **novelty** — genuinely emergent systems whose phenomenon is *not* in the
  catalog (diffusion-limited aggregation, Keller-Segel chemotactic aggregation,
  an active-nematic rod swarm, Eden growth), which should be
  EMERGENT-UNCLASSIFIED, not a false MATCH;
- **null** — non-emergent systems (uniform noise, bounded diffusion, uncoupled
  oscillators, a frozen field), which should be NO-EMERGENCE;
- **probe** — characterized-limit demonstrations, reported but not scored (5A.8).

**Operating point.** A screening-tier firing is exploratory evidence, appropriate
for self-recognition but too weak to assert a MATCH on an unfamiliar system: the
non-specific P18 voter-consensus screening gate, for example, fires on any growing
contiguous region. The battery therefore exposes a `match_min_tier` parameter.
Its default ("screening") preserves native behavior and leaves the confusion
matrix and the public gallery unchanged; for unknown systems the recommended
setting is "confirmation," which folds screening-only firings into the
emergence-driven verdict (surfaced as a `demoted_match`). It is **not** made the
global default precisely because P22, P27, and P29 self-recognize only at
screening (5A.4); a forced global gate would regress native self-recognition.

At the recommended ≥confirmation operating point, over three seeds:

| metric | result |
|---|---|
| recognition recall | 0.83 (6/6 rank the correct pattern #1; 5/6 fire) |
| novelty correct-abstention | 1.0 (4/4 EMERGENT-UNCLASSIFIED) |
| null specificity | 1.0 (12/12 NO-EMERGENCE) |
| **false-MATCH rate** | **0.0** |

The conservative gate is what delivers the zero false-MATCH rate: at the native
(screening) gate two near-miss novelties (DLA, Eden growth) trip the P18 voter
screening gate, and the false-MATCH rate is 0.21. (Report:
`docs/validation_rebuild/t2c_ood_validation.md`.)

## 5A.7 External demonstration: a multi-agent LLM swarm (T2d)

The strongest test of the instrument is a genuinely external multi-agent system
whose emergent behavior is not designed to match a catalog model. We constructed
one from independent large-language-model agents
(`analysis/t2d_llm_swarm.py`): twenty-four agents move on a 2-D torus, and each
agent, every round, receives **only its local egocentric neighborhood** as text —
the bearings, distances, and relative headings of the agents within its sensing
radius — and returns a single turn. No agent ever observes the global state, so
any global order is emergent from local LLM decisions rather than imposed. Three
swarms differ only in the one-line local rule, and each resulting trajectory was
profiled by the battery **blind** to the condition, at the recommended
≥confirmation operating point.

| swarm (local rule) | what emerged | instrument verdict |
|---|---|---|
| alignment ("match neighbors' heading") | polar order 0.31 → 0.97 | **MATCH P5 (flocking)** — φ = 0.96, Cohen's *d* = 32, null *p* = 0.005 |
| segregation ("toward your own type; stay when surrounded by it") | same-type fraction 0.58 → 1.0 | **EMERGENT-UNCLASSIFIED** — emergence 1.0 |
| random ("ignore neighbors, turn randomly") | polar order ≈ 0.15 | **NO-EMERGENCE** — emergence 0.17 |

The single demonstration exercises **all three instrument verdicts on a real LLM
multi-agent system** and discriminates the regimes by emergence kind
(orientational order for alignment, positional clustering for segregation, neither
for random). Two honest qualifications:

1. *The flocking detection is real but reaches screening, not confirmation, on
   this short trace.* P5's confirmation gate requires roughly ten system
   crossing-times of data; the 91-round run provides fewer. The signal itself is
   unambiguous (φ = 0.96, *d* = 32, *p* = 0.005); under the conservative gate the
   verdict is reported as strong emergence with P5 as the demoted top candidate.
   This is the same conservatism that produced the zero false-MATCH rate in 5A.6,
   applied honestly to an unfamiliar system; a longer run lifts P5 to confirmation.

2. *The segregation swarm is EMERGENT-UNCLASSIFIED by substrate, not by failure.*
   The catalog's aggregation detector (P1, Schelling) requires a **lattice**
   substrate; this swarm is in **continuous space**, so P1 rejects at prerequisite
   (`substrate_mismatch`). A spatial-binning observation adapter was tried, but a
   24-agent grid is too sparse to trigger P1. The honest reading is the one the
   instrument gives: the swarm self-organizes unmistakably (emergence 1.0) and is
   flagged as such, while no pattern is named because none matches the substrate —
   the none-of-the-above signal of 5A.3 working as intended. A denser swarm or a
   continuous-space aggregation detector would convert this to a firm MATCH P1.
   (Report: `docs/validation_rebuild/t2d_external_demo.md`.)

## 5A.8 Honest limitations of the instrument

- **Three characterized OOD limits, reported as probes, not scored.** (i) *P3 ↔
  phase ordering:* a fully-coarsened Allen-Cahn domain field is observationally
  indistinguishable from a coarse Gray-Scott Turing pattern by any single-field
  statistic (both stationary, isotropic, low-wavenumber, bimodal); the
  distinguishing physics — an intrinsic vs. a box-scaled wavelength — requires
  grid-size-invariance testing, which re-runs the model and is outside
  observation-only scope. (ii) *Structure without growth:* a pre-formed static
  structure scores as emergent because the order-gain term cannot see growth that
  already happened before the observation began. (iii) *Connectivity, not
  autocorrelation:* a percolation cluster at threshold is missed because its
  signature is graph connectivity rather than spatial autocorrelation. Each is a
  known boundary with a named resolving test, not a silent failure.

- **Calibration reference sizes.** The conformal tail probability floors at
  ~1/(|N|+1) for the current reference-set sizes, too coarse for hard *p* ≤ 0.05
  gates; the calibrated separation score is the headline figure and the verdict
  tiers are operating-point placeholders pending a held-out precision/recall sweep.

- **Emergence-indicator coverage gaps.** Rotational, banded, transient-decaying,
  and distributional morphologies (5A.3) are not yet covered; a low emergence
  score on those is inconclusive, and a universal indicator needs a battery of
  generic structure measures with appropriate nulls.

- **Substrate adaptation.** Continuous-space and graph systems require an
  observation adapter to a detector's native substrate; the LLM-swarm aggregation
  case (5A.7) shows both the need and a case where the naive adapter was too
  sparse at the available agent count.

## 5A.9 Reproducibility

The calibration layer, generic-emergence indicator, battery keystone, OOD suite,
held-out systems, and the LLM-swarm demonstration are all on the
`validation-rebuild` branch: `epc/phase2a/calibration.py`,
`epc/phase2a/emergence.py`, `epc/phase2a/battery.py`,
`analysis/battery_profile.py`, `analysis/t2c_ood_suite.py`,
`analysis/t2c_systems.py`, `analysis/t2d_llm_swarm.py`, and
`analysis/t2d_profile.py`. The instrument reports are
`docs/validation_rebuild/{battery_confusion_matrix,cross_model_generalization,t2c_ood_validation,t2d_external_demo}.md`.
The remaining roadmap item is a public command-line tool
(`epc-detect <observation-bundle>`) that wraps `profile_observation` at the
≥confirmation operating point.
