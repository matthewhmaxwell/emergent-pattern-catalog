# Ring 2 — Strategy (living document)

**Status: design phase.** The LENS framework below is settled (co-developed). Sourcing,
the literature-novelty gate, and the success bar are still under discussion — flagged
OPEN at the end. Nothing here commits us to building yet.

## The core honesty

Ring 2 hunts for emergent behaviors not in the catalog — ideally not in the literature.
The central problem is **completeness**, and it has two dual faces:

- **Lens completeness** — can we *see* a novel phenomenon?
- **Source completeness** — can we *produce* one?

Neither can be *proven* (unknown-unknowns). So we do not claim completeness; we **bound,
audit, and backstop** it. The strongest honest claim Ring 2 can ever make is "novel
relative to an audited, disclosed coordinate system," never "novel, full stop."

## Lens taxonomy (settled)

A *lens* is a detection channel; together the lenses form the emergence descriptor the
search navigates by. Sorted by the provenance / maturity of each lens:

### Tier 1 — grid-rooted lenses
Already in the instrument. Each earned its place by discriminating a real catalogued
pattern from its look-alikes at TNR 1.0 *in our own validation* — evidence-based in the
strong sense (demonstrated here, not merely literature-attested). Current set: morphology
/ clustering, orientation (polar + nematic), Ψ_CE synergy, temporal (spectral-peak),
heavy-tail / SOC, network (modularity / scale-free), ψ6 hexatic-gain, Pareto-tail
(heavy-tail-wealth), coexistence-oscillation — plus the 36 per-pattern detectors.

### Tier 2 — other evidence-based lenses, including from adjacent fields
Established order parameters / collectivity signatures from fields that study their own
emergence — to be **imported, not invented**. The mine:
- condensed matter — order parameters, broken-symmetry measures
- neuroscience — integrated information (Φ), criticality / avalanche scaling
- information theory — transfer entropy, predictive information / excess entropy
- applied topology — persistent homology / TDA
- nonlinear dynamics — Lyapunov exponents, recurrence quantification, attractor dimension
- network science — percolation / connectivity transitions, assortativity
- ecology — diversity–stability measures

**Discipline:** evidence-based-*elsewhere* ≠ trustworthy-*here*. Every import clears the
same bar the grid lenses did — discriminates real structure, does not false-fire on the
null set — before it joins the trusted descriptor.

### Tier 2-special — the model-free BRIDGE (the completeness tripwire)
A distinguished sub-class of Tier 2: **lens-agnostic complexity measures** — compression-
based complexity / logical depth, predictive information / excess entropy, causal
emergence (Ψ_CE; the MPR-complexity estimator is already built). These do **not name** a
structure. They report: *"there is real structure here, beyond randomness, that none of
the specific lenses caught."* That signal — high on the model-free measure, silent on
every named channel — is simultaneously:
1. the **prime novelty lead** (the instrument is looking at something it has no category for), and
2. the **trigger** that a Tier-3 lens needs inventing.

This is how the instrument reports its *own* incompleteness, and roughly *where*. It is the
hinge between "exhaust the known lenses" and "invent a new one."

### Tier 3 — novel lenses we invent
Not speculative: the discovery phase did this three times — ψ6 → P35, Pareto-tail → P36,
coexistence-oscillation → P37 — each forced by a phenomenon the existing lenses could not
see. Trigger: the Tier-2-special tripwire. Once invented and validated to the standard
bar, a Tier-3 lens graduates into Tier 1.

## What this replaces

The unanswerable "are the lenses complete?" becomes a **posture**:
> exhaust the evidence-based lenses (Tier 1 + Tier 2); keep the model-free tripwire armed
> for structure we cannot yet name; invent a Tier-3 lens when it trips; and report what we
> remain blind to.

An ever-expanding, audited coverage frontier — never a completeness claim. The blind-spot
audit (probe library → recall) is the instrument that measures the frontier and is itself
generalized for Ring 2 into a literature-spanning coverage test.

## Related work & positioning

ASAL — *Automating the Search for Artificial Life with Foundation Models* (Kumar et al.,
MIT / Sakana AI / OpenAI / IDSIA; *Artificial Life*, 2024) — automates ALife discovery
across substrates (Boids, Particle Life, Game of Life, Lenia, Neural CA) using a
**vision-language foundation-model representation** as its behaviour descriptor, in three
modes: target-prompt search, open-ended novelty in FM-representation space, and diversity
illumination.

How this work differs (stated plainly, not as competition): our descriptor is a set of
**validated, interpretable, emergence-grounded lenses** — each a discriminating detector
held to a TNR-1.0 bar against look-alikes — plus the model-free completeness backstop.
That lets us say *what* emerged, whether it is *real* (vs. an artifact or a look-alike),
and — via the literature gate — whether it is novel, in an audited, disclosed coordinate
system. ASAL's foundation-model descriptor is broad but opaque and perception-aligned
(it measures human-visual distinctiveness, not emergent mechanism) and is not validated
against a ground-truth emergence taxonomy. Different approach, different aim; complementary,
not competing — we do not frame this as a race.

**FM-as-lens policy.** A foundation-model signal is *not* part of our lens process. It may
be reconsidered as one more Tier-2 import only if it demonstrably clears the same bar every
lens must — discriminates real structure, does not false-fire on the null set, and
measurably improves coverage. Reputation does not earn a lens its place; evidence does.

## Lens ledger (running — the Ring-2 import sprint)

Every candidate lens must *earn its place*: discriminate genuine structure on a positive
control, stay quiet on the null set, and not false-fire on the mundane particle-life leads.
Status is recorded here so the lens library is an auditable artifact, not commit-message lore.

| Lens | Family / axis | Status | Discriminator & evidence |
|---|---|---|---|
| `structure_factor` | condensed-matter / characteristic length scale | **ADMITTED** (877cb91) | principal S(k>0) peak prominence; patterned 20–166 vs nulls ~6 |
| `persistent_homology` | topology / loops + voids | **ADMITTED** (this sprint) | `h1_max` (single dominant H1 loop); ring ctrl 0.58, vortex 0.43 vs all non-loop ≤0.12, gap +0.31, thr ~0.27. **NOT** `h1_total` (sum is noise-confounded). Positions-substrate only; field extension via sublevel-set filtration deferred. |
| `graph_structure` | network topology / interaction graphs | **ADMITTED** (this sprint) | two axes vs random-graph null: `degree_cv` (hubs) scale-free 0.89–1.00 vs random ~0.4 (gap +0.45); `modularity` (communities) 0.54–0.58 vs random ~0.35 floor (gap +0.14); clustering disambiguates small-world. Covers the network substrate class positions/field lenses can't touch. |
| `directed_info_flow` | causal direction / info-flow asymmetry | **ADMITTED** (this sprint) | directed transfer entropy among component series. `directionality` (net asym / magnitude): cascade 1.72–1.74 vs symmetric mesh 0.05–0.13 (gap +1.60); `mean_te` coupled ~0.06–0.19 vs independent ~0.004. New axis: who-drives-whom, independent of coupling magnitude. Gate: directionality meaningful only above the coupling floor; needs T≥60. |
| `fractal_dimension` | scale-free spatial structure (box-counting) | **DEFERRED** (this sprint) | box-counting D conflates self-similarity with density/boundary: random gas D~1.57 ≈ Sierpinski D~1.585, percolation D~1.83 ≈ uniform disk D~1.90. D is a shared absolute, not a signature. Needs lacunarity / multifractal spectrum. Works as a coarse strong-fractal flag (dla D~1.35) meanwhile. |
| `recurrence` (RQA) | nonlinear-dynamics / determinism | **DEFERRED** (877cb91) | RQA determinism confounded by smooth-stochastic trajectories (null_walk 0.54); needs RR-gating + phase-randomized surrogates before re-test |
| `novelty_tripwire` | model-free bridge (Tier-2-special) | **ARMED** (aab2e95) | fires on COMPLEX (MPR-C ≥0.16 ∧ Ψ_CE ≥0.05) ∧ UNCLASSIFIED; baseline-validated quiet (0/3 nulls, 17/17 classified) |

Next candidates (venv unblocks all): directed transfer entropy (beyond the global TE already
present), networkx graph-structure measures, optimal-transport drift. Then multi-family
substrates (Lenia / reaction-diffusion / coupled oscillators) to exercise the broadened set.

## OPEN — still under discussion (do NOT finalize)

- **Sourcing the systems under test** — generative families vs. program-space grammar vs.
  imported corpora; provenance + reproducibility; the streetlight risk of writing our own.
- **The literature-novelty gate** — how to establish novel-to-science despite the
  absence-of-evidence problem; how much rigor before we'd stake a claim.
- **The success bar** — novel-to-catalog vs. novel-regime-of-known-mechanism vs.
  novel-to-science, and what the deliverable is even if no genuine novelty turns up.
