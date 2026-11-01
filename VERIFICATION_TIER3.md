# Tier-3 vetting of the undocumented candidate qd:060 (written before any Tier-3 run)

**Status:** post-registration, exploratory. VERIFICATION_PLAN.md says any undocumented candidate gets the registered
Tier-3 checks; the census preregistration (§4, item 4) lists them: replication (10 seeds, two system sizes), rule
knock-outs, the interventional debunking battery (null and provably symmetric controls), a second independent
literature pass, and owner review. This file fixes how each is done for qd:060 before any of them is run. Owner
approval: the verification plan ("do it now", all phenomena) and deviation V4 (review by Claude, 2026-11-01).

## The candidate
qd:060, representative `CAF kC=2 kA=1 v=0.3 eta=3.14159 fields=1 D=[0.05] | A.CLIMB(f=0) ; A.WRITECELL() ;
C.GAME(g=(1.5, 0.5))`: 400 agents turn up the gradient of a diffusing field (initially 1.0 with 20 random 3x3 patches
at 0.5 and 1% noise, no source or decay), each agent advances the type of the cell under it every step, and cells play
a snowdrift game with imitate-the-best. Robust features (Stage 2): g_n_domains, g_activity_trend, g_corrlen,
g_reversals, g_largest_domain. No sub-program reproduces all five. Stage 3 (b) found no source for the behaviour or a
near relative.

## Already done (Stage 2)
1. Replication: 10 fresh seeds at census settings, 4x area and 3x length; all five features robust.
2. Rule knock-outs: every proper sub-program at census settings (in the card). Rule 3 alone gives the ordering
   (largest domain 0.82) and rules 1+2 the steep activity decline (-3.53); rules 2+3 give a disordered lattice
   (largest domain 0.34).

## The debunking hypothesis to test
H0 (decomposition into documented pieces): agents climbing a field that only decays converge on its maxima
(chemotactic trapping), so their cell flips become confined to few cells; the snowdrift game with imitate-the-best
then orders the rest of the lattice as it does alone. The activity decline is the agents' convergence; the ordering is
the game under sparse, localized perturbation.

## Controls (census settings, 10 seeds, the verification seeds of qd:060, frozen census fingerprint, view C)
- **N (null, registered):** the full program with `null_shuffle=True` (spatial structure destroyed after every step).
- **S (provably symmetric):** the full program with the field set exactly uniform (1.0 everywhere, no patches, no
  noise). The gradient is then exactly zero everywhere at every step, so the climbing rule never changes a heading.
- **P (pinned flippers, the debunking control):** the full program with agent speed 0 (agents never move), with the
  number of agents set to the median number of distinct cells occupied by agents at the last frame of the full
  program (from the same seeds), placed uniformly at random as in the census.

## Decision rules (fixed now)
For each robust feature f, a control *reproduces* f if |median(full) - median(control)| <= 3 s (s pooled over the two
seed sets as in the plan), the same test as sub-program reproduction.
- N and S are expected to reproduce no robust feature. If S reproduces all five, the gradient is not needed and the
  Stage 2 sub-program result is contradicted; that is reported as a failure of the harness, not a finding.
- If P reproduces the ordering features (g_largest_domain, g_n_domains, g_corrlen, g_reversals) and sub-program 1+2
  accounts for the activity decline (its activity trend has the full program's sign and at least half its size), H0
  is supported: the final label becomes CLOSELY RELATED (a decomposition into documented pieces: chemotactic trapping
  and a spatial game under localized perturbation, each documented, their coupling not), and the candidate is closed.
- Otherwise the candidate stays an UNDOCUMENTED CANDIDATE.

## Second independent literature pass (run in every case, after the controls)
A fresh research agent gets only the rule in plain words, the measured features, the sub-program and control results,
and the strict question (is this behaviour, or a near relative, documented?). It does not see the Stage 3 replies or
the Stage 3 (b) review. Its label is combined with the control result: a source that documents the behaviour makes it DOCUMENTED,
a near relative makes it CLOSELY RELATED, and only if the pass finds neither and H0 is not supported does the
candidate remain an UNDOCUMENTED CANDIDATE.

## Owner review
Covered by deviation V4: the owner has asked Claude to review. Whatever the outcome, qd:060 is reported as a candidate
or a closed candidate, never as a discovery.
