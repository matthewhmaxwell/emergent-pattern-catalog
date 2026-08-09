# Section 5C: From structure to competency — interventional discovery (Ring 3)

Section 5B's frontier search was passive: it read structure from recorded trajectories and asked
whether any region was complex yet unnameable. It found nothing novel, and in re-examining why, we
arrived at a sharper question about what the catalog is *for*. The patterns worth cataloguing are
not merely structured states; they are **competencies** — cases where a system does more than its
specification, achieving or pursuing an outcome it was not programmed to. This is the sense in
which the project's motivating example (Zhang, Goldstein & Levin 2024) is striking: a sorting
algorithm is interesting precisely when it exhibits a problem-solving competency its swap-rule does
not contain. Competency, however, is invisible to the Section 5B instrument, and not by accident.
This section reports the consequence: a different kind of detector, an interventional one, and the
map of competency-forms it produced. Its headline result is, again, negative and we argue correct —
zero competencies novel to science across thirty catalogued faces — with the contribution lying in a
validated, self-debunking method and a governing-dynamics map.

## 5C.1 Why passive observation cannot see competency

Levin's position is explicit and, on reflection, unavoidable: goal-directedness cannot be inferred
from passive observation; it must be established by perturbative experiment, and its signature is
William James's "same ends through different means" — when the normal route to an outcome is
blocked, a competent system reaches the outcome a different way. A recorded trajectory shows what a
system *did*, never what it *would have done* under obstruction, so the entire Section 5B apparatus —
however many lenses — is structurally blind to it. This also dissolves a circularity in the passive
search: a structural lens declares novelty relative to *known structures*, whereas a competency is
defined relative to the system's *own specification* (it does more than its code), a per-system,
non-circular baseline.

## 5C.2 The agent-observer

The detector is therefore an automated experimenter — a closed loop run with the burden of proof
against itself: observe the system; hypothesize what, if anything, it is achieving; design a
*novel barrier* such that genuine competency reaches the same end a different way while mere
mechanism breaks; predict both outcomes; intervene; and then actively try to explain the apparent
competency away as mechanism, re-equilibration, or relabeled specification. Only competency that
survives debunking is accepted; the default verdict is mechanism, and convergence-from-many-starts
is explicitly *not* competency (a ball reaches the bottom of a bowl from anywhere). Validated
black-box on a controlled triad, the loop returned competent / mechanism / illusion correctly,
including debunking a seductive convergence illusion by inventing the away-from-goal trap that a
greedy follower cannot pass. Run cold on a particle-life substrate with no hints, it correctly
declined the tempting case — clusters that re-form after a scatter — as passive re-equilibration
rather than competency, refusing the self-healing mirage.

The loop's defining property is that it is built to refute its own positives, and it repeatedly did.
In the foundational work it caught its own errors five separate times (a too-weak barrier
manufacturing false positives; a representation wall mistaken for rarity; an over-fit to favorable
test conditions; a narrow held-out over-reading a generalization rate; and an emergent-strategy
mirage). The competency-map and scale phases that follow added several more self-corrections —
memoryless baselines that were secretly recurrent, a reward design that made non-participation
optimal, two "untrainable" verdicts that were under-training artifacts (§5C.4), an eval-time
ablation that overstated a swarm's dependence on observation until a fair separately-trained
baseline corrected it, and a hard threshold that mislabelled a partial generalization as
"memorization." That self-debunking is the substantive difference from a passive perception-aligned
embedding (ASAL; Kumar et al. 2024): the method is interventional rather than observational,
interpretable rather than opaque, and adversarial toward its own positives. It is, in effect, the
interventional successor to the classic ALife emergence test (Ronald, Sipper & Capcarrère 1999,
"Design, Observation, Surprise!"), replacing observer-*surprise* — which over-fires and is
observer-psychological — with a *survives-a-collapse-ablation* criterion, and adding a
catalogued-versus-novel verdict that the surprise test does not provide.

## 5C.3 Discovering competency by program search

Because the catalog's patterns are all simple algorithms, the search is over the *space of simple
algorithms* (the rule), not the parameters of one model — new behaviors live in new algorithms, not
new settings of a known one. The funnel is: generate simple programs, cheaply filter for surplus
(non-trivial activity), then competency-probe the survivors. Two findings fixed the design. First,
a *representation* result: memoryless reactive rules cannot solve a barrier that requires committing
to a detour away from the goal — and neither can a hand-built reactive navigator — so the
controlled zero is a representation wall, not competency rarity. The interesting competencies
require internal state. Second, a *search* result: random sampling finds verified competency at
essentially zero rate, and a fixed training set is Goodharted (high train score, failed held-out);
only directed search with domain randomization — fresh conditions every generation, so nothing is
memorizable — discovers competency that survives a broad held-out battery on conditions never
trained on. Two competencies were so discovered: a generalizing navigator (partial, 0.71 on novel
mazes it never saw) and a robust working memory (cue-dependent action across unseen delays, with a
memoryless baseline pinned at chance). Both are genuine; both, by the literature gate, are *known*
competency-kinds.

## 5C.4 A governing-dynamics map of competency-forms

Iterating the detector across object-level, meta-level, multi-agent, and cross-generational
substrates produced a catalog of **thirty competency faces** (plus three *learner boundaries*, §5C.5),
each fingerprinted with an ablation battery — a channel is *load-bearing* only if scrambling or
removing it collapses the behavior — and each literature-gated to a known competency-kind. The full
catalog with per-face diagnostics is the supplementary `canonical_catalog.json`; here we report the
two axes it is organized by.

The first axis is **state-requirement**. Navigation and memory require internal state (a memoryless
policy scores at zero or chance), whereas delayed gratification (cash out at the peak of an
unknown-height reward), regulation (hold a setpoint against disturbance), and niche construction
(durably modify the world with no reward for doing so) are achieved *reactively* — adding state buys
little or nothing. A pre-registered hypothesis that "control/feedback" would be its own state-mode
was *falsified*: regulation came back reactive against both stochastic and structured disturbances.
State-requirement is thus a discriminating axis, not a universal property.

The second axis, for the state-requiring forms, is **state-mode** — what the held state is *for*:

- *commitment* — holding a detour against the goal-pull (navigation);
- *storage* — holding a cue across distractors (memory), or a partner's history (reputation / social
  memory);
- *accumulation* — a running tally toward a threshold (counting);
- *belief over a latent* — an inferred hidden rule (rule-inference, rule-tracking) or an
  other-agent model (predictive intention-reading / a theory-of-mind component; non-stationary
  co-player adaptation);
- *internal-state* — reading one's own energy, velocity, or morphology (metabolic-commons foraging,
  continuous-physics momentum control, morphology-conditioned specialization);
- *externalized* — coordination state written into the world rather than the agent (stigmergy, and
  its n-agent form stigmergic coordination; the cumulative-culture ratchet stores its state in a
  persistent artifact that outlives any individual);
- **joint-conditioning** *(new mode)* — conditioning on two simultaneously-necessary channels at
  once (compositional attention: own-morphology *and* an opponent's strategy-switch). This face
  emerges only under an attention (Transformer) policy; every recurrent variant carried one channel
  and dropped the other. It is the map's first genuinely two-channel state-mode, and it is reachable
  by changing the *learner architecture* rather than the task (its boundary complement is LB3, §5C.5).

Two honest-status notes qualify the map. First, several multi-agent faces are not new *modes* but
new *media* of coordination — a shared symbol channel (communication), direct mutual observation
(role-division, contest-resolution), and the environment itself (stigmergy) — and across the
multi-agent probes coordination reliably rides the *cheapest available medium*, an observation we
report as a proposed regularity rather than a validated law (§5C.6). Second, the reachable-competent
set is a **budget-dependent lower bound**: two faces (momentum control, morphology specialization)
were first logged "untrainable" at a modest optimization budget and both train cleanly at a larger
one, so short-run negatives are not trustworthy negatives and the map's boundaries move with
training budget. Each face was a genuine test that *refined* the map rather than confirming it —
the correct behavior of real data.

## 5C.5 The open-ended frontier: environment and learner-architecture bottlenecks

A known competency is recognized by its functional signature; a *new* one can only appear if the
detector can name a competency it was not given in advance. We first probed this open-ended naming in
a multi-affordance survival world (forage regrowing food, avoid hazards, manage depleting energy),
evolving survivors and asking the agent-observer to name, with no pre-specification, what the best
survivor does. The naming mechanism worked and *debunked* the tempting "patch-rotation to exploit
regrowth" reading via four interventions, finding instead a mundane reactive forager whose apparent
strategy was an epiphenomenon of greedy foraging in a regrowing cluster. No new competency emerged,
and the probe localized why: a task solvable by mundane means yields mundane solutions.

We then built and ran the route §5C names as most likely to reach a novel competency —
environment/agent co-evolution — as Minimal Criterion Coevolution (Brant & Stanley 2017, with
resource limitation): two populations, environments as evolvable genomes and recurrent agents with a
policy-gradient inner loop, each reproducing only on a binary minimal criterion, with no reward or
objective. Over four hundred generations the populations turned over continuously and produced many
competent behaviors, but a fair shared-channel re-screen against the catalog found **zero
novel-to-catalog** — an in-loop screen first flagged 1356 "candidates," all a support-mismatch
artifact of comparing sparse co-evolution fingerprints against the dense catalog, which the honest
re-screen resolved to existing behaviors one by one. Removing the human from objective-specification
did not escape the reachable set.

This refines the section's central claim. The limiting factor on discovering a new competency is not
the detector, and not the environment alone: it is the **environment and the learner architecture
jointly**. The cleanest demonstration is a *learner boundary* we mapped directly. In a task requiring
two simultaneously-necessary own-state channels, a recurrent (GRU) policy holds roughly one and
drops the other — an auxiliary loss can steer *which* channel it keeps but not make it hold both —
whereas holding the world fixed and swapping the recurrent policy for an attention (Transformer) one
lets both be necessary at once (this is LB3, the boundary complement of competency #30). The
instrument, in other words, is part of what determines the reachable phenomena. We record three such
boundaries — reciprocal cooperation in a social dilemma, resolved not by more capacity but by a
learner that shapes its co-learner's gradient (LOLA; Foerster et al. 2018); periodic
synchronization, resolved by inter-agent coupling (Kuramoto dynamics) rather than more training; and
the one-channel budget above — each a limit of the *learner*, gated known, and each pointing to a
richer learner class, not a richer objective or world, as the un-crossed lever toward novelty.

## 5C.6 Honest status

No competency novel to science was found across thirty catalogued faces, and on the stated prior none
was expected from these substrates. We take care not to overclaim the *reason*: we do not assert the
space of competency-primitives is closed — that is unfalsifiable at this sample size — only that
minimal-sufficient mechanisms are exactly the region of behavior-space science mapped first, so a
search that (like any gradient or evolutionary learner) settles on the *cheapest sufficient*
mechanism is structurally biased toward rediscovering known primitives. That tendency is itself a
consolidation, into the multi-agent competency setting, of established results — shortcut learning
(Geirhos et al. 2020), simplicity bias (Shah et al. 2020), specification gaming (Krakovna et al.),
and the compression–expressivity account of when structure emerges (Kirby et al. 2015) — and it is
bounded to the toy scale we tested; open-ended-evolution theory (Lehman & Stanley 2011) holds that
the minimal basis is *escaped* by climbing complexity, which is consistent with our reading that a
richer learner or world, not this one, is where novelty would live.

The contribution is, as in Section 5B, methodological, and it is the durable result of the work: an
interventional, self-debunking competency detector (the agent-observer), grounded in the requirement
that goal-directedness be established by perturbation, positioned as the interventional successor to
the ALife design-observation-surprise test and distinguished from the emergent-abilities-as-mirage
critique (Schaeffer et al. 2023) and the pitfalls-of-measuring-emergent-communication line (Lowe,
Foerster et al. 2019) by acting on systems, naming what it finds, and being built to refute its own
positives; a directed-search discovery funnel that genuinely finds verified competency in evolved
simple programs where passive search and random sampling find none; and a governing-dynamics map of
thirty competency-forms organized by state-requirement and state-mode. Three further regularities
the corpus suggests — a cost hierarchy over coordination media (communication is forced only by
information asymmetry), a "shortcut-cascade depth" measure of how many cheaper mechanisms must be
foreclosed to force a target competency, and a temporal-scale axis (within-episode, within-lifetime,
across-generations) — we report explicitly as **proposed, not validated**: each is composed of known
pieces and would need a dedicated experiment (a dialable information-structure ladder; a
well-definedness proof; genuine multi-timescale runs) before any novelty claim. What remains, and the
route we judge most likely to reach a genuinely new competency, is a richer learner class inside
open-ended environment/agent co-evolution, with the same literature gate reserved for any off-map
survivor.
