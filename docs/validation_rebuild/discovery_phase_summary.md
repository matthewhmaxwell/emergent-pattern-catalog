# Discovery Phase — Growing the Catalog with Its Own Instrument

**Honest framing (read first).** The three patterns added in this phase — P33 active
nematic, P34 adaptive-network fragmentation, P35 entropy-driven crystallization — are
**established emergent phenomena from the scientific literature**. They are novel **to
this catalog**, *not* novel to science. This phase is **coverage discovery**: using the
validated instrument (and the catalog's own ontology) to find known-but-uncatalogued
emergence and give each a *validated discriminating detector*. The contribution is
(1) the detectors, which did not previously exist; (2) the demonstration that the
instrument can predict gaps and self-grow the catalog; (3) the detection-channel
widening required to even see one of them. Genuine novelty-to-science is a separate,
harder target (Ring 2, not yet attempted; the honest prior there is that it is rare).

## The discovery loop

Each pattern was taken through the same cycle, to the same Phase-2a bar the original 32
were held to:

  find → verify → catalog → recognize

- **find**: surface an EMERGENT-UNCLASSIFIED candidate (the instrument detects emergence
  but no catalog detector fires).
- **verify**: build a minimal model + a discriminating detector that fires on the target
  and rejects look-alikes; validate TNR (true-negative rate) = 1.0 + a continuous effect
  size, across seeds and finite-size where applicable.
- **catalog**: promote to `epc/{models,metrics,detectors}` + register in the battery.
- **recognize**: `profile_observation` now returns MATCH for the new pattern, and
  nothing else false-matches it.

Two search strategies ("rings") fed the loop:

- **Ring 0 — known-but-uncatalogued sweep.** Run a library of canonical emergent systems
  the catalog doesn't contain through the instrument; vet the EMERGENT-UNCLASSIFIED hits.
  (`analysis/discovery/sweep.py`, `ring0_systems.py`.) → P33.
- **Ring 1 — dimensional gaps (periodic-table logic).** Map the 32 patterns onto the
  11-dimension emergence ontology (`docs/ontology_v0_4.md`); the EMPTY cells are
  predictions of where a pattern *could* live. (`analysis/discovery/coverage_map.py`.)
  → P34 (empty `substrate=evolving_net`), P35 (empty `interact=none_entropy`).

A recurring methodological finding: **the clean discriminator is almost always a
*gain* (emergence from a disordered start), not an absolute value** — because
look-alikes share the absolute (random close packing has local ψ6≈0.5; a static modular
graph has high modularity Q; a coarsening field has a dominant wavelength). Measuring the
change from the initial state is what isolates the phenomenon. This holds across P3
(wavelength growth), P34 (modularity gain), and P35 (ψ6 gain).

## The three patterns

### P33 — Active nematic order with ±1/2 topological defects
- **Established science:** active nematics / active turbulence. Sanchez et al. 2012
  (*Nature* 491:431, microtubule–kinesin); Doostmohammadi et al. 2018 (*Nat Commun*
  review). Decades of active-matter literature.
- **Distinct from catalog:** P5 flocking is *polar* (φ high); P6 milling is a single
  +1 vortex. Active nematics are *apolar* (φ≈0) with motile ±1/2 (half-integer) defects.
- **Discriminator:** coherent half-integer defect density = ±1/2 defect density gated by
  local nematic order high, polar order low, angular momentum low.
- **Validation:** positives 5/5 definitive; TNR 20/20 = 1.000; continuous d = 4.81;
  real P5 Vicsek rejected; finite-size definitive at G = 64/96/128.

### P34 — Adaptive-network fragmentation
- **Established science:** adaptive/coevolutionary networks; the adaptive-voter
  fragmentation transition. Holme & Newman 2006 (*PRE* 74:056108); Vazquez et al. 2008;
  Gross & Blasius 2008 review (*J R Soc Interface*).
- **Distinct from catalog:** P18 voter, P21 polarization, P22 cascade all run on FIXED
  topology. Here the *topology coevolves* with node states and fragments into
  disconnected same-opinion communities — the emergence is in the network structure.
- **Discriminator:** coevolution-driven modularity gain = Q_late − Q_early, gated by
  final modularity high, started-connected (low initial Q), giant-component
  fragmentation, opinion–topology assortativity.
- **Validation:** positives 5/5 definitive; TNR 25/25 = 1.000; continuous d = 16.74.
  (Regime: rewiring probability φ = 0.99, above the fragmentation transition.)

### P35 — Entropy-driven crystallization
- **Established science:** the Alder transition — hard/repulsive disks crystallize from
  entropy (excluded volume) alone, with no attraction. Alder & Wainwright 1957
  (*J Chem Phys* 27:1208).
- **Distinct from catalog:** P1 is type segregation, P2 is active motility. Here ordering
  is passive, identical-agent, repulsion-only. The signature is hexatic
  (bond-orientational) order, which a clustering measure misses (a crystal has
  near-uniform density).
- **Discriminator:** hexatic-order emergence = local-ψ6 *gain* from a disordered fluid,
  gated by reaching a hexatic floor. Absolute ψ6 ≈ 0.5 is shared with random close
  packing / aggregation, so only the gain is specific.
- **Validation:** positives 5/5 definitive; TNR 30/30 = 1.000; continuous d = 8.31
  (vs dilute gas, static crystal, Keller–Segel clumps, active nematic, noise, walk).
  *Caveat:* the minimal soft-disk sim produces polycrystalline order; a single-crystal
  demonstration would need event-driven Monte Carlo. The ψ6-gain discriminator is
  robust regardless.

## Instrument improvements driven by discovery

- **P3 sharpening (Ring-0 finding).** The sweep showed P3 Turing over-claims on
  Cahn–Hilliard conserved coarsening (a coarsening field is observationally a
  stationary, isotropic, low-wavenumber periodic structure). Fix: an
  intrinsic-wavelength test — a Turing pattern locks its wavelength; coarsening's grows
  monotonically. P3 re-validated TNR 1.0, positives 5/5; cahn 0/3 false-match. Known
  residual limit (documented): non-conserved Allen–Cahn box-saturation needs grid-size
  invariance (out of observation-only scope).
- **Round-3 ψ6 channel.** P35 was initially a detection *blind spot* (no
  bond-orientational channel). Added a ψ6-gain channel to the generic emergence
  indicator (additive — does not fire on flocking/clumps/random). Blind-spot audit
  recall 15/15, 0 null false-positives; T2c null-spec 1.0, false-MATCH 0.0 — no
  regression. The net now covers the hexatic family for future systems too.

## Predictive validation of the ontology  (★ write up in the next paper pass)

The headline scientific result of the discovery phase — flagged for §5.8 / §7.5.

A descriptive taxonomy's empty cells could be empty for trivial reasons (physically
unrealizable combinations, or artifacts of the classification scheme), carrying no
predictive content. Ring 1 is a first test of whether this ontology's gaps are instead
*predictive*:

1. The 11-dimension ontology was constructed to *describe* the 32 catalog patterns
   (post-hoc, descriptive).
2. Coverage analysis (`coverage_map.py`) identified two ENTIRELY EMPTY ontology values:
   `substrate = evolving network` and `interaction = none / entropy-driven`.
3. We treated each empty cell as a *prediction* — "an emergent pattern should exist
   here" — and built a canonical minimal system in the cell.
4. Both cells turned out to host **established, independently-discovered emergent
   phenomena** the catalog had simply missed: adaptive-network fragmentation
   (Holme & Newman 2006) and entropy-driven crystallization (the Alder transition,
   1957). Each then passed the full Phase-2a discrimination bar (TNR 1.0).

This is the **periodic-table property**: the ontology's gaps predict *real* phenomena,
not artifacts — evidence that the dimensions carve emergence at meaningful joints, and
that the instrument can grow its own catalog by coverage analysis. **Two predictions
made; two confirmed.**

Honesty guardrails the paper claim MUST carry:
- This validates the method *against the existing literature* — the predictions were
  confirmed against phenomena already known to science (coverage discovery), not a
  discovery of anything novel-to-science.
- n = 2 confirmed predictions is suggestive, not conclusive; the ontology axes are
  descriptive tags, not a proven generative axis (unlike atomic number in Mendeleev's
  table). More empty cells should be tested to strengthen the claim.
- The decisive future test — a gap prediction that lands on a phenomenon ABSENT from
  the literature (genuine novelty), via generative-family novelty search — is the open
  frontier, with the honest prior that genuine novelty is rare.

## State

Catalog **32 → 35**; battery **34 detectors**; all validated, committed, pushed
(branch `validation-rebuild`). The discovery engine is demonstrated across five
substrates (orientation field, reaction–diffusion field, evolving graph, point packing,
plus the catalog's existing types) and two search strategies (sweep + ontology gaps).

**Not yet done:** Ring 1 combination-cell gaps; Ring 2 (composition / generative-family
novelty search), which is where genuine novelty-to-science would be hunted — with the
honest expectation that it is rare.
