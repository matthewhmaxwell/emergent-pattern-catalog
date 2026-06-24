# Ring 1 — Combination-Cell Investigation (closure)

**Status: CLOSED.** Catalog 35 → 37. This document records the outcome of probing the
plausible two-way *combination* cells of the 11-dimension emergence ontology, after
first-order coverage was completed.

## First-order coverage: complete

With P33 (active nematic), P34 (adaptive-network fragmentation), and P35 (entropy
crystallization) added, `analysis/discovery/coverage_map.py` confirms **every one of
the 11 first-order ontology values is now occupied** — zero single-value gaps.
(`interaction = none/entropy` → P35; `substrate = evolving-network` → P34.)

## The four curated combination cells

Coverage analysis flagged four physically-plausible, empty two-way cells. Each was
treated as a *prediction* and probed to the same standard the catalog is held to.

| Cell | Predicted phenomenon | Outcome |
|---|---|---|
| `update=resource_exchange × homog=heterogeneous` | heterogeneous kinetic exchange | ✅ **P36** — Pareto wealth tail |
| `interact=field × conflict=competitive` | field-mediated resource competition | ✅ **P37** — Huisman–Weissing chaos |
| `substrate=fixed_net × update=reproduction` | network reciprocity | ⊘ **evidenced bound** (statistical signature) |
| `driving=ext_forced × feedback=nontransitive` | driven cyclic dominance | ⊘ **non-prediction** (probed) |

### ✅ P36 — Heterogeneous kinetic exchange (Pareto wealth tail)
Chatterjee–Chakrabarti–Manna (2004). Per-agent saving propensity λ_i ~ U[0,1) →
power-law (Pareto) wealth tail, P(w)~w⁻². Discriminator = Pareto-over-exponential
tail advantage (Gini alone is shared with P28 condensation; the *tail shape* is the
specific signature). **TNR 25/25 = 1.0, continuous d = 13.87, 5/5 definitive.**
Round-4 net-widening: a heavy-tail(pareto-wealth) generic-emergence channel.

### ✅ P37 — Field-mediated resource competition (Huisman–Weissing chaos)
Huisman & Weissing (1999/2001). Species competing through shared depletable resources
coexist via sustained oscillation/chaos beyond the competitive-exclusion limit.
Discriminator = multi-species coexistence + sustained abundance oscillation; the
stable-coexistence control (same coexistence, settles to a fixed point) is rejected,
so the signature is the non-equilibrium dynamics, not coexistence. **TNR 15/15 = 1.0,
continuous d = 15.09, 5/5 definitive.** Round-5 net-widening: a coexistence-oscillation channel.

### ⊘ `fixed_net × reproduction` — network reciprocity: evidenced bound
Network reciprocity (Santos–Pacheco 2005; Ohtsuki–Nowak 2006) is a *real* phenomenon
for this cell, but it does **not** yield a clean within-single-observation detector at
the catalog's bar. Five configurations were tested — synchronous Fermi (cooperation
collapses), asynchronous proportional imitation (defector-hubs, degree–coop corr ≈ −0.27,
bistable), long runs (unchanged), death-birth (fixation to all-C), and DB+mutation
across the full b/c range (coexistence but degree-**indifferent**, corr ≈ 0). The
cooperator-hub signature — the only within-observation structure that would distinguish
it from P27 (lattice spatial reciprocity) — never robustly appears.

The reason is principled: this cell's phenomena (network reciprocity, evolutionary
amplifiers, epidemics on scale-free graphs) have **statistical / comparative
signatures** — a cooperation *level* enhanced *relative to a baseline*, fixation
*probabilities* — not the single-observation structural fingerprints that P36 (a tail
shape) and P37 (sustained oscillation) leave. A faithful detector would require a
different, statistical instrument format. Recorded as an honest bound; model WIP kept
at `epc/models/network_games.py` (untracked).

### ⊘ `ext_forced × nontransitive` — driven cyclic dominance: non-prediction
Probed via externally-forced spatial rock-paper-scissors (cyclic external drive on the
P12 substrate), `analysis/discovery/cell4_driven_cyclic_probe.py`. The driven system
shows the **same field-type emergence** as the autonomous cyclic-dominance pattern
(emergence score ≈ 0.9, `kind=field`), with no distinct new structural or temporal
signature registered by the instrument (the forcing entrains/perturbs but does not
create a novel pattern). Within the probed scope, this cell is a non-prediction.

## Scientific conclusion

Of the four plausible combination cells, **two host clean, validated emergent patterns
(P36, P37) and two do not** (one a real phenomenon outside the observation-only format,
one a non-prediction). This is the honest, useful result: the ontology's gaps are
**often, but not always, predictive within the instrument's format.** That nuance
*strengthens* the periodic-table claim — the dimensions carve emergence at meaningful
joints often enough to be predictive, while the bounds delimit the claim and
characterize the instrument's scope: it detects *single-observation structural*
emergence, not statistical/comparative phenomena.

A recurring secondary finding: **each new pattern revealed a generic-net blind spot**
(P36 heavy-tail-wealth, P37 resource-competition-oscillation), closed by an additive
gated channel (round-4, round-5). The catalog's coverage analysis surfaces *detection*
gaps, not only *pattern* gaps.

## State after closure
Catalog **37 patterns**, battery **36 detectors**. T2c OOD null-spec 1.0 / STRICT
false-MATCH 0.0 / nov 1.0; blind-spot recall 17/17; no regression. Branch
`validation-rebuild` == `main`. Next frontier: Ring 2 (generative-family novelty search).
