# Ring 3 — Synthesis: the structure of the emergent-competency space

*Consolidation of the interventional competency search (registry: 21 competencies; branch `validation-rebuild`).*

## 1. What we set out to do
Detect **competency** in emergent systems — a system "doing more than its code," in Levin's operational sense
(established by *perturbation*; signature = *same ends by different means*) — and honestly classify each instance
as **catalogued** (known to science) or **novel** (new-to-science). The standing goal was to find something
genuinely new; the standing method was to never claim novelty until it survived a literature gate.

## 2. The instrument (method)
An **agent-observer, interventional** detector with a **default-skeptical** verdict: a behavior is presumed to be
a known mechanism unless it survives *debunking*. Debunking = targeted ablations that should collapse the behavior
iff the claimed competency is real (memory-wipe, channel-scramble, blind-partner, feedback-ablation, anonymize,
no-artifact, provably-symmetric controls, train/held-out splits, culture-buffer removal). Two methodological
refinements the search itself forced:
- **Fair baselines, not eval-time ablations.** Ablating a cue at eval on a policy trained *with* it produces
  degenerate behavior that over-states dependence. The fix is a separately-*trained* baseline (used for the
  memoryless control in ToM #20 and the Condorcet control in the swarm probe — both corrected inflated readings).
- **Provably-symmetric controls** are the gold standard: when a single frame is symmetric *by construction*
  (ToM #20), the baseline sits at chance mathematically, so any above-chance performance is unambiguous.

The instrument is *self-debunking*: it caught and corrected **~9 of our own over-claims** across the arc (broken
memoryless baselines, reward-design artifacts, focal-point escapes, confounded ablations). This honesty tooling
is itself a contribution — most emergence claims in the wild are over-claimed for exactly these reasons.

## 3. The catalog (21 competencies, 4 levels)
Every one **gated to KNOWN** → the honest bound is **0 new-to-science**.
- **Object (world-facing):** navigation, memory, counting/accumulation, sequencing, delayed-gratification,
  regulation, tool-use/instrumental-construction.
- **Meta (rule-facing):** rule-inference (1-shot), rule-tracking (continual), conditional-selection,
  abstraction/systematic-generalization (partial), predictive intention-reading / theory-of-mind (#20).
- **Multi-agent (other-facing):** communication (Lewis), role-division, contest-resolution (Bourgeois),
  fair-alternation, count-communication (compositional), stigmergy, reputation/indirect-reciprocity,
  distributed-info integration.
- **Cross-generational (culture-facing):** cumulative-culture ratchet (#21) — the first competency realized
  across lifetimes rather than within one.

## 4. Three structural regularities (the Track-2 payoff so far)
The search did not find a new competency, but it produced three candidate **regularities** about the *structure*
of the competency space. Each is composed of known pieces; the candidate contribution is the *synthesis*.

**(a) The coordination-media COST HIERARCHY + COMMUNICATION-NICHE theorem.**
Emergent coordination uses the *cheapest available* medium, in a fixed cost order:
`environmental focal-point < observation < coordinate focal-point < communication`.
Established by an 8-rung *cascade* (fingerprint the medium in use → remove it → re-run). The sharp claim:
**symmetric coordination never forces communication** — a cheaper focal point always exists; communication is
forced *only* by information ASYMMETRY (single-source → referential; distributed → composed). This UNIFIES the
whole multi-agent catalog (every genuine communication case had asymmetry; every forced-symmetric case used a
focal point). It holds across DOMAINS (space + time), REGIMES (2/3/25 agents), and KINDS (coordination, social
choice, collective computation). Corollary: **the primitive-media set is CLOSED within a lifetime** — which
*explains* the 0-new-to-science bound (competencies are minimal compositions of a closed, catalogued basis).

**(b) SHORTCUT-CASCADE DEPTH** (a competency-complexity measure).
To force a target competency you must adversarially remove *every cheaper mechanism*. The **number** of shortcuts
that must be closed is a countable measure of how high a competency sits on the mechanism ladder. Theory-of-mind
(#20) required closing **six** successive single-frame shortcuts before velocity-integration emerged; simple
competencies require zero. This turns "richer competency" from a vibe into a number.

**(c) The TEMPORAL-SCALE axis.**
Competencies stratify by the timescale over which they are realized: within-episode → within-lifetime →
across-generations. The cumulative-culture ratchet (#21) is the first at the across-generations scale, and the
only probe whose mechanism does not reduce within a lifetime. Refined completeness statement: **the primitive set
is closed within a lifetime; cross-generational accumulation is a new SCALE, not a new primitive.**

## 5. The unifying law: MINIMAL SUFFICIENT MECHANISM
The through-line behind all of the above, replicated ~a dozen times (evolutionary and gradient learners; 2- to
25-agent; spatial/temporal/cultural): **a learner realizes the cheapest mechanism the task's structure permits.**
Consequences: (i) emergent competencies decompose into a small closed set of catalogued primitives; (ii) the
minimal mechanism is always something science already named → 0 new-to-science; (iii) forcing a richer mechanism
requires adversarially removing every cheaper one (→ the cascade method and the depth measure); (iv) composition
appears only under *sufficiency pressure* (when no single minimal mechanism suffices — the compositional-comms
and distributed-integration results).

## 6. Honest scope + what a genuinely-new result would require
Everything here is at CPU-scale toy RL. The bound and the completeness argument are honest *for that substrate*.
A genuinely new-to-science competency would require a **new primitive** — which, on this evidence, needs a
materially different substrate (LLM-agent societies, real embodiment, open-ended artifact/tool evolution), not a
bigger version of the same toy. That is a directional claim about *where* to look, earned empirically, not an
a-priori dismissal.

## 7. The candidate contribution (for the field)
Not a new competency, but a **theory + method for the emergent-competency space**:
1. A **self-debunking interventional instrument** for honestly gating emergence claims (most are over-claimed).
2. The **coordination-media cost hierarchy + communication-niche theorem** (a predictive account of *which*
   coordination mechanism emerges = which media the environment leaves available, in cost order).
3. The **closed-basis / minimal-mechanism explanation** for why emergent-competency search returns known results,
   with two quantitative handles (shortcut-cascade depth; the temporal-scale axis).
Track 2 will hard-gate these against the literature to determine which is genuinely under-named.
