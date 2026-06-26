# Section 5B: The lens battery and the genuine-novelty frontier (Ring 2)

Section 5.8.1 established that the ontology's empty cells are often predictive, but every
confirmed prediction landed on a phenomenon already in the literature — coverage discovery
(novel-to-catalog), not novel-to-science. The decisive test is a prediction, or a search,
that lands on a phenomenon *absent* from the literature. Reaching for it requires both a
different search — sampling continuous generative families rather than reading a discrete
ontology — and a richer instrument, because a novelty search can only flag as novel what its
lenses cannot already name. This section reports a first, fully audited pass at that frontier.
Its headline result is negative and, we argue, correct: zero novel-to-science findings, with
the contribution lying in the auditable method and the instrument improvements the search
drove.

## 5B.1 The discovery ceiling is lens comprehensiveness

The reach of an instrument-steered novelty search is bounded by the breadth and
interpretability of its lens set: a system is a lead only when the instrument judges it
structured yet has no category for it, so a missing lens manufactures false leads and a
mislabeling hides real ones. This framing is also the cleanest statement of how the present
work differs from foundation-model approaches to artificial-life search (notably ASAL; Kumar
et al. 2024), which steer by a single, opaque, perception-aligned vision-model embedding. That
embedding measures human-visual distinctiveness, is not validated against a ground-truth
emergence taxonomy, and cannot say *what* a flagged system is. The battery below is the
opposite design point: many lenses, each interpretable, each validated to a discrimination
standard, each naming a specific structural axis — plus a literature gate before any novelty
claim. The two approaches are complementary, not competing; we make no race of it.

## 5B.2 A model-free bridge

Between the named lenses (which assert *what* a structure is) and the open frontier (structure
with no name) sits a lens-agnostic bridge: a tripwire that fires when an observation is
COMPLEX yet UNCLASSIFIED. Complexity is read from two substrate-free measures — Martin–
Plastino–Rosso statistical complexity and a PID-based causal-emergence estimate (Ψ_CE) — and
"unclassified" means no named lens and no catalogued detector recognizes it. The thresholds
are null-calibrated against noise, random-walk, and random-graph baselines. Out-of-distribution
exposure during the search (Section 5B.4) hardened the bridge twice: a surrogate structure-gate
(the raw permutation-entropy complexity is finite-size-biased upward on short iid series, so
the complex verdict now also requires the series to be less random than its time-shuffle), and
a dead-state gate (a system that collapses to a trivial constant final state has a complex
transient but no sustained structure, and is no longer counted complex). Both gates were
verified to leave the calibrated baseline unchanged.

## 5B.3 An imported, validated lens battery

Five lenses were imported from adjacent fields, each held to the same bar the catalog's
detectors meet — *earn your place*: separate genuine structure on a positive control, stay
quiet on the null set, and not false-fire on the mundane leads the search itself produces.

- **structure_factor** — characteristic length scale (a peak in the spatial power spectrum;
  condensed matter).
- **persistent_homology** — topology (loops/voids), with a Vietoris–Rips path for point clouds
  and a superlevel-set path for fields.
- **graph_structure** — interaction-network topology (degree heterogeneity, modularity,
  clustering), the substrate class spatial and field lenses cannot ingest.
- **directed_info_flow** — directed information transfer (who drives whom), distinct from
  coupling magnitude and from spatial propagation.
- **fractal_dimension** — scale-free spatial structure, admitted in a scoped form: box-counting
  dimension alone conflates self-similarity with density, so the discriminator is lacunarity
  (gliding-box heterogeneity), and the scope is sparse self-similar aggregates, with dense
  near-critical fractals left to a future multifractal treatment.

A sixth candidate (recurrence-quantification determinism) was retired as redundant: its target
is already covered by the spectral and synergy channels, and its standalone statistic could not
be separated from a smooth random walk even with surrogate correction. The recurring lesson —
the same one that shaped the catalog's discriminators — is that the signature is always a
*specific* structural feature, never a shared absolute: the most-persistent single loop rather
than the sum of loop persistences; lacunarity rather than the bare dimension; flow
*directionality* rather than flow magnitude. The admitted lenses, the named emergence indicator,
and the model-free bridge are composed into one substrate-aware descriptor: each lens self-
guards off-substrate, so the fingerprint is sparse and the set of lenses that fire is itself the
coordinate subspace an observation occupies.

## 5B.4 An audited first pass: zero novelty, two instrument gains

The descriptor was run over 690 configurations spanning four generative families —
particle-life (signed interaction matrices), Lenia (continuous cellular automata, the ASAL
search space), Gray–Scott reaction–diffusion (the Pearson (F,k) atlas), and Kuramoto
oscillator lattices — clustering each observation by its fingerprint and surfacing the
tripwire's complex-and-unclassified leads.

**Zero genuine-novelty leads survived vetting.** Three families produced none. Lenia produced a
handful that vetting traced to dying configurations whose death transient is temporally complex
but whose final state is trivially dead — closed by the dead-state gate. Kuramoto produced
seven that vetting identified as twisted and phase-wave states: locally phase-coherent but
globally incoherent, imposed by the initial condition and persisting — known oscillator physics,
not novelty. The honest prior — that genuine novelty is rare — held.

The reliable deliverable is therefore not a discovery but an **audited frontier**: an
interpretable, disclosed map of the explored emergence-space, in which we can state exactly
what was searched, along which named axes each region was characterized, and why nothing
qualified. This auditability is the substantive difference from an opaque-embedding search.

The search also performed its second function — surfacing the instrument's own gaps. The
Kuramoto twisted states exposed a real blind spot: the emergence indicator's phase channel keys
on global coherence and so misses locally-ordered, globally-incoherent phase fields. This is the
same local-versus-global distinction that the active-nematic pattern (P33) forced for
orientation fields, now recurring for oscillators. We closed it with a local-phase-order channel
(mean neighborhood coherence on the phase grid, scored against a phase-shuffle null), a
net-widening of exactly the kind the catalog's earlier coverage rounds performed, and re-validated
the modified indicator with no regression (out-of-distribution null-specificity 1.0, recognition
and novelty rates unchanged; blind-spot recall 17/17). The seven leads now classify correctly.
The search thus realizes a self-improving loop — search, surface a blind spot, widen the net,
re-validate — driven not by hand-specification but by the instrument's own failures on
out-of-distribution input.

## 5B.5 Honest status and what remains

No phenomenon novel to science was found in this pass, and on the stated prior none was expected.
The contribution is threefold: an interpretable, validated, substrate-aware lens battery and
model-free bridge; an auditable novelty-search engine that produces a disclosed map rather than a
claim; and a demonstrated self-improvement loop in which the search hardens the instrument. What
remains is breadth — more generative families, finer parameter grids, larger samples — with the
literature-novelty gate reserved for any future complex-and-unclassified survivor. A
foundation-model signal would be reconsidered as one more imported lens only if it demonstrably
clears the same bar every lens here had to clear; reputation does not earn a lens its place,
evidence does.
