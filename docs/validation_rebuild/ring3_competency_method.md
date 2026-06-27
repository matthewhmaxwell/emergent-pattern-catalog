# Ring 3 — Interventional competency detection (the agent-observer)

**Status: METHOD VALIDATED. Competency found so far = navigation (a KNOWN form). No new
science yet.** This documents a pivot, not a discovery.

## Why this exists
Rings 0–2 hunted by *passive structural observation* (a battery of pattern-derived lenses) and
found 0 genuine novelty across 440+ configs. The honest diagnosis, reached by going back to first
principles: behavior worth cataloguing is **competency** — a system *doing more than its code*,
in the sense of pursuing/achieving an outcome it was not programmed for. Michael Levin's position
is explicit: **you cannot infer goal-directedness from passive observation; you must do
perturbative experiments** — the signature is William James's *"same ends through different
means."* Our entire prior approach was structurally blind to this. Hence the inversion:

> Detect **behavior itself** (competency) *interventionally and generically*, FIRST; only then
> compare to the catalog (same / variant / new) and check cross-substrate recurrence (Levin's
> "Platonic" forms). The catalog stops being the detector and becomes only the *namer*.

## The detector: the agent-observer
An automated Levin-scientist — a closed loop, with the burden of proof set against itself:
1. **Observe** the system untouched; hypothesize what (if anything) it is *achieving*.
2. **Design a discriminating barrier** — a *novel* obstacle such that genuine competency reaches
   the same end a different way, and mere mechanism breaks. State the prediction under both
   hypotheses first.
3. **Intervene, judge, then DEBUNK** — actively try to explain the apparent competency away as
   mechanism, re-equilibration, or relabeled spec. Competency is accepted only if it survives.

Default verdict is *no competency*. Convergence-from-many-starts is **not** competency (a ball
rolls to the bottom of a bowl). The decisive trap is a barrier whose only escape requires moving
**away** from the apparent goal.

## What was built and validated (files in analysis/ring3_competency/)
- `levin_barrier_probe.py` — falsifiable toy: same-ends-different-means under a novel barrier
  separates competency from mechanism (clean on navigation; flat on 1-D sorting, which taught us
  competency is multi-faced and barrier design is the crux).
- `agent_observer_probe.py` + `agent_observer_toy_systems.py` — LLM-as-scientist, black-box, on
  three opaque systems. Returned **A=competent, B=mechanism, C=illusion** — correctly, and
  **debunked the seductive convergence illusion (C)** by inventing the away-from-goal trap.
- `plife.py` — particle-life substrate with intervention hooks. The LLM scientist, **cold and with
  no hints**, correctly called System 1 *self-organising / re-equilibration (NOT competency)* —
  refusing to be fooled by clusters that re-form after a scatter — and System 2 a *gas*. Neither
  earned competency.
- `embodied_funnel.py` / `stateful_funnel.py` — **program search** (vary the rule, not parameters):
  generate simple algorithms → surplus filter → competency probe → generalization debunk.
  - Memoryless reactive rules: **true 0/4000** solve the away-from-goal trap; a hand-built reactive
    navigator *also* fails it → that 0 is a **representation wall**, not rarity. Committed detours
    need memory.
  - Stateful finite-state controllers: a hand Pledge/DFS control solves the trap (test is passable).
    Random search surfaced **seed 2959 (1/6000)**, which initially looked competent — BUT this was
    **RETRACTED**: the catalog-classification broad-battery (open-field reach + re-route + away-from-goal
    + generalise, from a FULL start set) showed 2959 is **fragile / start-specific** — it fails
    open-field reach from corner starts (0,0),(12,0). The earlier "generalization debunk" used a
    narrow favorable (upper-left) start set, which over-fit. Honest result: **0/6000 random programs
    are verified robustly-competent.** The only core-competent navigator is the **hand-built Pledge**
    (designed, not discovered; and even it fails the serpentine). The broad functional battery is the
    real verifier — a narrow debunk is insufficient (it caught this over-fit across sessions).
- **Directed search** (`evolve_navigation.py`) — the method that actually works. Random sampling: 0
  verified. Evolution on a FIXED training set: train 0.92 but held-out fail (Goodhart). Evolution +
  **domain randomization** (fresh random conditions every generation, so nothing is memorisable):
  converges, and the best survivor is a **genuinely generalising navigator — 0.71 on 123 SOLVABLE
  novel multi-wall mazes it never saw** (vs 0 for random, vs 2959 which failed open-field corners).
  NOTE: its narrow held-out signature read 1.0; the broad 123-maze debunk corrected it to 0.71
  (partial, not robustly general) — the 4th time a narrow test over-read and the broad battery fixed
  it. So: directed-search FINDS genuine generalising competency; broad verification KEEPS IT HONEST;
  the 0.71→robust gap is an engineering knob (more states / longer search / richer sensing). The
  competency is still navigation = a KNOWN form (not new science).

## Honest findings
1. The pipeline works end-to-end and **self-corrects**: it caught its own false positives twice
   (weak-barrier 27→0; then the representation wall). The self-debunking *is* the differentiator
   from a passive embedding (ASAL).
2. **Memory is required** for the interesting competencies (committed detours / delayed
   gratification). A search over memoryless rules is doubly inadequate.
3. The competency *target* (navigation/re-routing) is a **KNOWN form** → method validated, **not
   novelty-to-science.** No new phenomenon discovered — and no robustly-competent *random program*
   discovered either (the one candidate, 2959, was retracted by the broad battery). The only
   competent navigator on record is hand-built.
4. Catalog-classification (analysis/ring3_competency/catalog_classify.py): a FORM = functional
   signature, not name/substrate; tiers = new-form / known-form-NEW-INGRESSION / dup. Per the
   Platonic refinement, a known competency in a *structurally-distant algorithm* is a genuine data
   point (another portal), and the set of ingressions reveals the form's GOVERNING DYNAMICS (e.g.
   navigation requires internal state + goal/obstacle sensing). Its broad battery doubles as the
   strongest competency verifier we have.
4. Cost: the LLM agent-observer is a per-system experimental loop (minutes); for a *known*
   competency-face the barrier-test is cheap and programmatic. Random program sampling is sparse
   (1/6000) → **directed search** (evolve toward competency) is what scales.

## Next
- **Catalog-classification of verified hits**: known competency (like 2959's navigation) = same/
  variant; the prize is a verified competency whose form the catalog lacks = candidate new, then
  check **cross-substrate recurrence** (Platonic).
- Richer schemas + more competency-faces (memory, anticipation, generalisation, regulation).
- Directed search over stateful program spaces to make competency findable at scale.
