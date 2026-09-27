# EPC Emergence Census — Design v1 (for review; nothing built yet)

**Owner:** Matt Maxwell · **Drafted:** 2026-09-27 · **Status:** awaiting approval
**Goal (the year-long one):** find genuinely novel emergent behavior in simple algorithms — and, whatever happens, leave a
contribution that is meaningful on its own.

**Decisions locked in the design conversation (2026-09-27):** success = map + novelty · scope = fixed rules + simple
adaptation · budget = weeks on the existing VPS · rigor = pilot, then pre-register on OSF · search = hybrid (exhaustive
census + novelty arm) · delivery = daily digest only (no immediate alerts).

---

## 1. Contribution (framed to survive zero discoveries)
1. **Emergence-complexity table.** For every catalog phenomenon: the **shortest program** that produces it (under two
   independent grammars, with their rank correlation), its **robustness** (the fraction of all quantized parameter settings
   of that minimal program's rule structure that still produce it), and its **necessary ingredients** (primitive knockouts,
   after Ma et al. 2009). This gives the catalog a *generative* ordering principle — the "atomic number" the periodic table
   has been missing — and answers the original "trends among them all" question non-circularly (the census counts every
   program, not just admitted ones).
2. **Census statistics.** Prevalence of each phenomenon by program length; rate of "unclassified" emergence vs a null base rate.
3. **Catalog expansion.** Every known-but-uncatalogued phenomenon the census surfaces gets a validated detector.
4. **Novelty dossier.** Candidates that pass every pre-registered gate — possibly empty, reported either way.

## 2. Positioning (from the 2026-09-27 literature check)
- **Done elsewhere:** simplest-first enumeration of abstract rule systems (Wolfram, *NKS* 2002); automated novelty search
  over *parameters* of fixed families with image-embedding novelty (ASAL, Kumar et al., *Artificial Life* 2025 — closest
  prior work); minimal structure for *one* phenomenon at a time (Ma et al., *Cell* 2009; Scholes et al., *Cell Systems* 2019).
- **Open:** a *cross-substrate structural* census with a *validated mechanistic* filter; a *catalog-wide* ordering by
  shortest generating program.
- **Distinctions we must keep:** search program **structure**, not parameters of known families; mechanistic detectors,
  not CLIP; concentrate novelty effort on **ingredient intersections** (where new phases historically appear: MIPS 2008,
  non-reciprocal transitions 2021). Complementary to, not competing with, symmetry-based universality classes.

## 3. Grammars
- **G1 (high-level primitives).** Substrates: lattice cells, particles on a 2-D torus, network nodes, up to two diffusing
  fields with reaction terms. Entity state: discrete type (2–4 values) plus heading/phase where relevant. Primitive
  families: *sense* (neighbor counts by type, mean heading, field value/gradient, own payoff), *move* (toward/away, align,
  density-dependent speed), *state* (copy, threshold switch, cyclic dominance, reproduce/die, conserved transfer),
  *field* (deposit, react, diffuse/decay), *network* (rewire), *adapt* (imitate-best, win-stay-lose-shift with 2×2 payoffs).
- **G2 (low-level primitives).** Same capabilities built from finer operations (e.g., "align" = sense-mean-heading +
  set-heading). Exists to test grammar dependence.
- **Length.** An explicit prefix-free code charges bits for every choice; parameters are quantized
  (probabilities {0.05, 0.2, 0.5, 0.9}, thresholds {1–4}, radii {1, 2}, diffusion {0.05, 0.2}, decay {0.01, 0.1}) and
  charged in bits. Program length = total bits.
- **Deduplication.** Canonical forms: drop no-op rules, merge equivalent rule orders, collapse state relabelings.

## 4. Simulation
Fixed per substrate: lattice 64×64; 400 particles; 200-node network; fields on a 64×64 grid; 1,000 steps; 3 seeds per
program. The pilot may adjust sizes/steps for throughput; final values are frozen at registration.

## 5. Search (a ~3-week run, compute split as below)
- **Census, G1 (~50% ≈ 1.5 weeks).** Exhaustive over every canonical G1 program of length ≤ **L**, where L is the largest
  length whose measured pilot throughput fits this share. Output: proven-shortest (within G1, ≤ L) program per
  phenomenon; phenomena not reached are reported as "> L". Follow-ups on each minimal program: robustness sweep and
  ingredient knockout.
- **Census, G2 (~15%).** Exhaustive over every canonical G2 program of length ≤ **L₂**, chosen the same way for this share.
  Rank correlation is computed on phenomena found in both grammars; if the ordering is unstable, that instability is
  reported as a result.
- **Novelty arm (~35%).** Quality-diversity (MAP-Elites) over programs longer than L that **couple at least two of**
  {particles, lattice, network, fields, adaptation}; seeded from the census's highest-emergence programs; behavior space =
  the frozen 21-lens Ring-2 fingerprint.

## 6. Filter pipeline + catalog-expansion loop
1. **Screen** (every program): discard dead, frozen, or noise-level runs via the existing generic-emergence channels.
2. **Known filter:** the 37-detector battery at the **confirmation** tier. MATCH → phenomenon recorded in the map.
3. **Unclassified:** emergent in ≥ 2 of 3 seeds and no MATCH. Never called "novel".
4. **Triage:** unclassified programs are clustered by **agglomerative clustering on the standardized 21-lens
   fingerprint** (linkage and distance threshold set in the pilot, frozen at registration). Each cluster gets a
   **literature check on its specific behavior** by an independent research agent with web search and verified citations
   (the protocol used for the P1 literature pass).
   - Known but uncatalogued → **catalog expansion**: new detector + Phase-2a validation panel (negative controls,
     TNR 1.0) → the known filter is **re-run over every previously unclassified program** → phenomenon enters the map.
   - Literature-silent → **novelty candidate**.
5. **Vetting (novelty candidates):** replication (more seeds, two system sizes), ingredient knockout, interventional
   debunking (null and provably-symmetric controls), independent literature pass. Only then any claim.
- **Base-rate control (pre-registered):** the whole pipeline also runs on **null programs** (non-interacting entities;
  shuffled rules). A candidate's fingerprint must fall outside the null programs' unclassified region (threshold frozen at
  registration).

## 7. Validation before anything counts
- **Recovery benchmark:** each textbook model — Vicsek flocking, Schelling segregation, BTW sandpile, Kuramoto
  synchronization, Gray-Scott Turing patterns, May-Leonard rock-paper-scissors spirals, MIPS, co-evolving-voter
  fragmentation, plus every other catalog entry the grammar can reach — is first **hand-encoded in G1**. The census must
  then produce that phenomenon at a length **no greater than the hand encoding**, and the filter must recognize it.
  Publish the confusion matrix on these positives and the false-positive rate on null programs. Failures are fixed **in
  the pilot, before freeze**. Entries outside a spatial/network/field grammar (e.g., sorting algorithms, Hopfield memory)
  are labeled "out of grammar scope", not forced.

## 8. Visibility — the daily digest
- **Channel:** a daily digest only (no immediate alerts). Everything is logged on the VPS as it happens, so a missed day
  loses nothing — the next digest catches up.
- **Contents, in priority order:** Tier 4 *vetted* findings and Tier 3 *novelty candidates* flagged at the top; then
  Tier 2 *catalog expansions* and Tier 1 *unclassified clusters*; then census progress (programs run, current length,
  phenomena found, the emergence-complexity table so far).
- **Delivery:** a scheduled daily check-in task reads the VPS status, writes the digest, and refreshes a **private
  dashboard page**.
- **Candidate card** (every Tier 3): animation; the program in plain words; nearest known phenomena and why it is not
  them; knockout results. The owner's judgment is part of the gate.

## 9. Process and timeline (~5 weeks)
- **Phase 0 (≈ 1 day):** VPS disk is at 95% → *archive* old run artifacts to the DO Space (list shown to owner before
  anything moves; nothing deleted without approval).
- **Phase 1 — pilot (3–5 days, not counted):** build G1/G2, simulator, pipeline, digest + dashboard; recovery benchmark;
  null base rate; throughput → L, L₂; small census slice to debug.
- **Phase 2 — freeze + OSF registration (≈ 1 day):** see §11.
- **Phase 3 — full run (≈ 3 weeks):** niced, resumable, checkpointed tmux jobs as matthewhmaxwell.
- **Phase 4 — analysis + write-up (≈ 1 week).** Separate paper from the instrument paper (the instrument is this paper's filter).

## 10. Stopping rules
- Census: stops when every canonical program of length ≤ L (G1) and ≤ L₂ (G2) has run.
- Novelty arm: stops at its compute share, or when no new archive cell is filled for N consecutive generations
  (N set in the pilot, frozen at registration).

## 11. Frozen at registration (checklist)
Grammar specs + code hashes · bit-coding scheme · L, L₂ · simulation sizes/steps/seeds · screen thresholds · battery tier
(confirmation) · unclassified definition · clustering linkage + threshold · literature-check protocol · null programs +
base-rate threshold · catalog-expansion procedure (incl. the full re-filter) · vetting steps · stopping rules (incl. N) ·
analysis plan (rank correlation, robustness, knockouts) · digest tiers.

## 12. Risks and honest priors
- **Novelty prior is low** (40 years of mining small rule spaces; random programs mostly produce simple outputs —
  Dingle et al. 2018). The design's value does not depend on a discovery.
- **Grammar dependence** of description length is the main scientific risk → two grammars + rank correlation; instability
  is reportable.
- **Detector coverage:** 37 detectors cover a fraction of named phenomena (chimera states, traveling bands, active
  turbulence, chiral phases, explosive synchronization…) → most unclassified hits will be known; the expansion loop
  exists for exactly this, and its cost (≈ hours–days per new detector) is the main schedule risk.
- **Compute:** a shared 8-core VPS; L is chosen from measured throughput rather than hoped for.
