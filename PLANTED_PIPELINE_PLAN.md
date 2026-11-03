# Planted behaviors through the whole census pipeline (plan, frozen before any planted model is simulated)

**Status:** post-registration, exploratory, added 2026-11-03. Owner approval: "run the planted behaviors through the
whole pipeline" (2026-11-03). This file is committed, and its hash recorded in git, before any planted model is
implemented or simulated. It extends PLANTED_CONTROLS_PLAN.md (commit 4fcf65f), whose "What this can and cannot show"
says the stages before the gate were not tested for these behaviors. Nothing registered for the census (e5hkf) or
fixed by the verification plan is changed; frozen census code is imported unchanged.

## Why
The gate-only test measured whether the literature gate can say "not documented" when handed a card written from the
paper. The census never hands the gate such a card. A behavior has to be simulated, pass the screen and the interaction
knock-out, escape the namer, be described by the census card builder from the census fingerprint, pass the gate, and
then pass simulation verification and the source review. This test runs the planted novel behaviors through those
stages and records where each one is lost.

## Behaviors and scope
The 20 N cards of the planted controls (frozen at de277c8). A behavior is **in engine scope** if its model's own rules
and parameters can run unchanged in a world the census engine observes: a 2D square periodic lattice of discrete
types or phases, a 2D periodic box of moving particles with positions and headings, a graph (any graph, including the
complete graph) with node types or phases, or 2D periodic continuous fields. Only system size and run length are set
by the census scale below. Fixed now:
- **In scope (13):** N01, N02, N03, N04, N07, N08, N11, N12, N13, N14, N15, N17, N19.
- **Out of scope (7):** N05 (one-dimensional ring), N06 and N18 (three-component spin vectors), N09 and N16 (three-
  dimensional space), N10 (aperiodic tiling), N20 (triangular lattice with open edges started from imposed defect
  textures). These count as not reachable by the census engine: a coverage limit, reported separately from the
  sensitivity of the stages.
- N17 is a census grammar-1 program as it stands: `A kA=1 v=0.5 eta=2.0 | A.ALIGN(r=1.0)` (16 bits, census index
  2358; heading noise uniform in [-1, 1] rad = eta*pi with eta = 0.318, inside the paper's window 0.28 to 0.35). It is
  re-simulated with the frozen interpreter on its census seed. The other 12 are new simulator modules ("plug-ins")
  under `census/planted_models/`.

## Census scale (all plug-ins)
- Lattices and fields: 64 x 64 sites at the paper's own grid spacing (the domain shrinks, the pattern scale in sites
  is kept). Particles: 400 at the paper's density (box side set by the density). Graphs: 200 nodes; parameters the
  paper defines per population size (such as noise of order 1/N) keep the paper's scaling. Particle models that
  create or remove particles (N11, N14) cannot fill the census particle view, which needs a fixed number of particles;
  they are recorded as the coarse-grained particle density on a 64 x 64 grid over the paper's domain (a field view),
  with the paper's domain if it holds at most 10,000 particles, otherwise a domain holding about 400 at close packing.
- 3 seeds, 1,000 census steps, recorded at the census intervals (frames every 2 steps for lattices and particles, 5
  for phases, 10 for fields). One census step is one update of every rule for discrete-time models, one sweep (one
  attempted update per entity) for Monte Carlo and random-sequential models, and one unit of the paper's
  dimensionless time for continuous-time models, with the integrator's sub-steps inside it. Exceptions fixed now:
  N01 one free period; N11 one mean lifetime at the density-independent death rate; N15 one time unit (200 steps of
  0.005).
- Initial conditions follow the paper (random where the paper's are random).
- Compute cap: if the census-scale run (3 seeds, real and knock-out) would take more than 60 CPU-minutes, the time per
  census step is reduced to the largest that fits, and the change is recorded in the module before any pipeline run.
- Each model runs at one **headline setting** inside the range where its paper reports the new behavior (listed
  below). The census tested single programs; a sweep is not part of this test.
- Particle headings: the propulsion direction for self-propelled particles; the body axis for passive anisotropic
  particles (N04 ellipse axis, N08 dumbbell axis, N14 ellipse axis).

## Interaction knock-out and rules (fixed per model)
The census principle: remove every rule through which an entity's update reads another entity; keep one-way rules
(an entity's own propulsion, noise, decay, feed, firing, emission or deposition, an external drive) and the header
(sizes, speeds, noise levels, diffusion of fields). Where a coupling is a strength, it is set to zero. Verification
sub-programs remove interaction rules only (one-way rules and the header are kept in every sub-program); this differs
from the census verification, whose programs had at most four rules, and is fixed now because several plug-ins have
five to seven terms. Interaction rules are marked (I).

| Card | View(s) | Headline setting | Rules ((I) = interaction) | Knock-out |
|---|---|---|---|---|
| N01 | graph phases + node types (complete graph) | linear charging, kick k/N | 1 charge and fire; 2 pulse kick with absorption (I) | kick = 0 |
| N02 | fields u, v | a = 0.015, b = 10, D1 = 0.01, D2 = 1000 D1, c1 = c, c2 = -c, c = 0.17 | 1 decay of u; 2 feed and decay of v; 3 conversion a v -> u (I); 4 autocatalysis u^2 v (I); 5 u drifts along grad v (I); 6 v drifts along grad u (I) | 3-6 removed |
| N03 | particles | 400 disks, area fraction 0.05, screening 0.20, turning strength inside the paper's bubbly window | 1 propulsion; 2 chemical attraction force (I); 3 chemo-repulsive torque (I); 4 rotational noise; 5 contact repulsion (I) | 2, 3, 5 removed |
| N04 | particles | density 0.87, e_a/e_p = 33, T = 6 (restored order) | 1 orientation-dependent pair potential (I); 2 transverse odd pair force (I); 3 thermal noise | 1, 2 removed |
| N07 | field (mass per site) | keep 40% / 60% shares as on the card | 1 sender hops its remainder one step in a random direction; 2 receiver pools, keeps its share and passes the rest one step to its left (I) | 2 removed (receiver keeps what it gets) |
| N08 | particles (dumbbell centres, axis) | area fraction 0.65, Pe = 5, T = 0.01, paper's mass and friction | 1 spring inside a dumbbell (header); 2 constant torque; 3 inertia, friction, noise; 4 bead contact repulsion between dumbbells (I) | 4 removed |
| N11 | density field | g = 3.14e-5 (stripes) | 1 birth with Gaussian dispersal; 2 density-independent death; 3 two-body competition (I); 4 three-body competition (I) | 3, 4 removed (population growth beyond 50x the start is a blow-up and stops the run, as in the census) |
| N12 | lattice phases + cell types (k = 1) | beta = 1 | 1 long-range coupling J(r) ~ 1/r^4 (I); 2 thermal noise | J = 0 |
| N13 | graph node types (complete graph) | the paper's trimodal setting of eta and shock rate | 1 noise flip; 2 copy a random other agent (I); 3 population-wide shock to mean opinion 0.3 (external) | 2 removed |
| N14 | density field | T = 0.145, mu inside the paper's kagome window | 1 hard-core exclusion (I); 2 patch bonds (I); 3 exchange with the reservoir; 4 translation and rotation moves | 1, 2 removed |
| N15 | particles | 400 at density 1 (box 20), R = 5, cone half-angle pi/8, noise 0.08 | 1 propulsion with small inertia; 2 steering toward visible particles (I); 3 rotational noise | 2 removed |
| N17 | particles | census program idx 2358 | census interpreter | census knock-out |
| N19 | fields b, a, n | the paper's parameters (card) | 1 b drifts up grad a (I); 2 b drifts up grad n (I); 3 a produced by b; 4 a decays; 5 n consumed by b; 6 n resupplied (diffusion = header) | 1, 2 removed |

Production by cells and consumption by cells write to a field without reading another entity's state beyond the
site's own cell density, as census emission and deposition do; they count as one-way.

## Implementation and fidelity (before the code freeze)
- Each plug-in states the paper's model equations and parameters in its docstring, with the section of the paper they
  come from, and implements: the simulation at a given size and length; the knock-out and sub-programs by rule list;
  output arrays in the census format, so the frozen screen, fingerprint, knock-out decision, textbook measures, namer
  and verification classifier run on them unchanged.
- **Fidelity check:** each plug-in reproduces one observation from its paper at the headline setting, with the paper's
  own measure, at the largest size and length that fit 2 CPU-hours (a few models listed below also need a second
  setting for comparison). The check, its tolerance and its result are recorded. A failed or infeasible check does
  not remove the model; its results are reported as "implementation not confirmed". Fixed checks: N01 probability of
  full synchrony at N = 4, 5, 6 (paper 0.887, 0.99968, 0.906); N02 pattern amplitude larger at c = 0.17 than at
  c = 0.01 and absent at 0.49; N03 cluster fraction and cavities inside clusters; N04 bond-orientational order at the
  headline above the passive system at the same temperature; N07 density correlations decaying close to r^-4 (fitted
  exponent 3.5 to 4.5) where the achiral version does not; N08 two coexisting densities at Pe = 5 and none at Pe = 0;
  N11 stripes at the headline and hexagonal spots at g = 3.27e-6; N12 squared magnetization decreasing slowly with L
  at beta = 1 as the paper reports; N13 three peaks of the mean-opinion distribution, one near 0.3; N14 ring order
  after the quench; N15 a non-coarsening stream network with fractal dimension near 1.5; N17 not testable within the
  budget (the paper's lanes need boxes of side ~1,200); N19 one dense aggregate coexisting with travelling bands.
- **Presence at census scale:** the same measure is computed on the census-scale real runs and reported next to the
  pipeline outcome, so that a loss caused by the small census world is not charged to a pipeline stage. N01's new
  behavior is a property of many starts across N and cannot appear in one run; it is recorded as "not observable in a
  census run".
- Plug-ins are written by Claude agents from the papers and checked rule by rule against the paper by a separate
  agent before the freeze. Fidelity runs read only the paper's measure. No pipeline output (screen, knock-out, name,
  card, gate) is computed or read before all modules are committed (the **code freeze**).

## Pipeline (frozen census code; same order as the census)
1. **Flag** (census runner): screen in >= 2 of 3 seeds and interaction-driven by the frozen knock-out decision, per
   view; fingerprints and every applicable textbook measure stored for flagged views.
2. **Name** (census triage, library `census/reflib/census_v3`): each flagged view is NAMED, KNOWN-BEHAVIOUR-ATYPICAL-LOOK
   or UNNAMED. NAMED views never reached the gate in the census; the behavior is lost there (its name recorded).
3. **Card:** UNNAMED views get the census literature card (rule in plain words; parameter variants "(none distinct)";
   the six largest-magnitude mean fingerprint features with the frozen glosses; "appears across 1 distinct rule
   programs"); atypical-look views get the post-run review card. The rule text is the frozen planted card's rule text,
   edited only where the simulated size, length or setting differs, and checked against the planted leak list before
   any gate call.
4. **Gate:** census protocol as in PLANTED_CONTROLS_PLAN.md (first pass claude-sonnet-5-5, 10 searches; strict pass
   claude-opus-5-5, 12 searches, on every first-pass KNOWN-VARIANT; for atypical-look cards the post-run review's own
   strict-pass trigger and prompt), 3 independent runs per card, target blocked. Block list: the planted list
   (`census.litcheck_planted.BASE_BLOCK`, journal and author-page domains, without shared hosts such as github.com),
   plus arXiv mirrors and aggregators, fixed now: pith.science (the host behind every URL leak in the planted runs),
   ar5iv.org, papers.cool, emergentmind.com, chatpaper.com, synthical.com, scholar.archive.org, paperdigest.org,
   catalyzex.com, arxiv-sanity-lite.com, connectedpapers.com, academia.edu. Leaks are detected and adjudicated as in
   the planted controls and leave the denominator.
5. **Verification** (VERIFICATION_PLAN.md Stage 2, frozen `verify.classify_one`): 10 fresh seeds of the program and its
   knock-out at census settings, 4x area (lattice and field side doubled; particles x4 at the same density; nodes x4)
   and 3x length; sub-programs at census settings; class ARTIFACT, DECOMPOSES, COMBINATION EFFECT or SINGLE MECHANISM
   from the card features (or the features a mismatch-flagged gate verdict names). Every planted behavior that
   reached a gate verdict is verified, as in the census. A condition that cannot finish within 24 CPU-hours per model
   is recorded as not run; the frozen classifier then treats that feature as not persistent, and such cases are
   reported separately.
6. **Stage 3** for COMBINATION EFFECT and mismatch-flagged SINGLE MECHANISM: (a) the strict prompt with the simulation
   evidence (`verify_stage3` blind card format), target blocked as in step 4; (b) one Claude research agent per card
   checks the cited sources with the brief of the census source review (VERIFICATION_PLAN Stage 3 (b), V4), told to
   count only sources dated before the target paper's first version, without being told why. Labels DOCUMENTED,
   CLOSELY RELATED, UNDOCUMENTED CANDIDATE.

## Scoring (fixed now)
- **Gate stage, per card:** correct = LITERATURE-SILENT, or KNOWN-VARIANT whose reply states that the observed behavior
  on the card is not described (decided by a blinded reading, with the quoted sentence; one reader, as in the planted
  controls). Card result = majority over its valid runs.
- **End-to-end, per behavior (primary):** survives (a) on the registered census path if a flagged, unnamed or
  atypical-look view's card gets final label LITERATURE-SILENT in the majority of its valid runs; (b) on the
  verification path if it ends Stage 3 as UNDOCUMENTED CANDIDATE. Reported: per-stage pass counts, the stage where each
  behavior is lost, end-to-end shares over the 20 cards and over the 13 in scope with exact 95% intervals, and the
  lenient gate rule (KNOWN-VARIANT read as correct) alongside.
- **Descriptive:** whether any card feature measures the paper's new behavior (read against the paper's measure);
  the gate result for the pipeline card next to the hand-written card's result in the planted controls; presence at
  census scale; fidelity.
- **Novelty dispute:** recomputed with the planted controls' disputed set (N14, N16, N17, N19 treated as disputed).
- Thresholds are not changed after any result is seen. A rule that proves unworkable when first run is logged as a
  deviation, with its reason, before any affected result is read.

## What this can and cannot show
It shows where real behaviors that are new by an external standard are lost in the census pipeline, for the behaviors
the census engine can host, at census scale. It cannot separate pipeline sensitivity from scale for behaviors that do
not appear in a 64 x 64 or 400-particle world (presence is reported); 12 of the 13 in-scope models run on new code
written for this test, not the frozen interpreter; implementations, readings and reviews are done by Claude agents of
the same family as the gate; one program per behavior is tested, where the census grouped many; and the authors'
novelty claims are taken as stated, as in the planted controls.
