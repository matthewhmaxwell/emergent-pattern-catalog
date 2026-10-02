# EPC Emergence Census — Implementation Plan (Phase 1 pilot → freeze)

**Design:** `DESIGN.md` v1 (approved by owner 2026-09-27: "proceed").
**Code:** git branch `emergence-census` — VPS worktree `/home/matthewhmaxwell/epc-census` (off `validation-rebuild`), package
`census/`. Local mirror of docs: `~/VPS-staging/projects/epc-emergence-census/`. All VPS jobs run as matthewhmaxwell in
tmux, niced.

## Phase 0 — disk archive: NOT NEEDED
VPS root disk 89% (18 GB free, 2026-09-27); all EPC work ≈ 16 MB; census footprint estimated < 3 GB (per-program summary
rows; raw trajectories kept only for flagged programs). Nothing is moved. Revisit only if free space < 8 GB.

## Phase 1 — pilot (3–5 days; nothing here counts as census data)
Each task ends with a check that can be verified on its own.

| # | Task | Done when |
|---|------|-----------|
| T1 | **G1 grammar spec** (`GRAMMAR_G1.md`): substrates, primitives, parameter sets, prefix-free bit code | spec reviewed against the recovery list (every benchmark model expressible) |
| T2 | **Program core** (`census/program.py`, `census/coding.py`): AST, encoder/decoder, canonicalizer, length-ordered enumerator | round-trip encode→decode on 10k random programs; count-by-length table; canonical dedup tests |
| T3 | **Interpreter** (`census/sim/`): lattice, particles, network, fields, adaptation; seeds vectorized; deterministic | unit test per primitive; hand-coded Vicsek/Schelling/voter/sandpile/Kuramoto/Gray-Scott/RPS/MIPS reproduce textbook behaviour visually + numerically |
| T4 | **Battery adapter** (`census/filter.py`): simulator output → generic emergence screen → 37-detector battery at confirmation tier → 21-lens fingerprint | adapter runs every applicable detector on each substrate's output without error |
| T5 | **Recovery benchmark**: hand encodings in G1 of every reachable catalog entry + textbook model; filter must MATCH each; enumerator must reach each at ≤ its hand length | confusion matrix on positives; failures fixed here (grammar, adapter, or labelled "out of grammar scope") |
| T6 | **Null programs + base rate**: non-interacting and shuffled-rule programs through the full pipeline | unclassified rate + null fingerprint region (threshold proposal for freeze) |
| T7 | **Throughput** → choose L (G1, ~50% of 3-week budget) and L₂ (G2, ~15%) | measured programs/s per core on the census mix; L, L₂ proposed |
| T8 | **G2 grammar** (low-level primitives) + enumerator + hand encodings of the benchmark | benchmark models expressible; count-by-length table |
| T9 | **Triage**: agglomerative clustering on standardized 21-lens fingerprints; candidate cards; literature-check agent protocol (P1 protocol) | linkage + distance threshold proposed from pilot data |
| T10 | **Novelty arm**: MAP-Elites over programs coupling ≥ 2 of {particles, lattice, network, fields, adaptation}, behaviour space = 21-lens | runs on a pilot budget; archive fills |
| T11 | **Visibility**: daily digest (tiers 1–4, candidate cards) + private dashboard page; scheduled daily check-in | one digest produced end-to-end from a pilot slice |
| T12 | **Pilot slice**: full pipeline over every G1 program of length ≤ a small ℓ | runs clean; resumable; checkpointed |
| T13 | **Freeze**: `PREREG_census_v1.md` (DESIGN §11 checklist), code hashes; owner reads and approves; OSF draft → register | owner approval; registration id recorded |

## Phase 2 — freeze + OSF (≈1 day) → Phase 3 — full run (≈3 weeks) → Phase 4 — analysis + write-up (≈1 week)
As in DESIGN §9. Registration visibility is the owner's call at T13 (the Ring-3 prereg was public; the OSF project node stays private).

## Standing rules
- Nothing from Phase 1 counts as census data; pilot outputs live in `census/pilot/` and are labelled as such.
- Any post-hoc analysis is labelled post-hoc. Every probe/detector gets a negative control (lesson of Ring-3 prereg D4).
- Literature pre-check before any novelty claim (lesson of P1).
- Commits only on `emergence-census`; push to origin.

## Pilot log (Phase 1; nothing here is census data)
**2026-09-27 (day 1)** — branch `emergence-census`, commits 3da21b2 → e7e10a8.
- **T1/T2 done (G1 v0):** 31 mechanism templates over layers C/A/N/F (+ CA, CF, AF, CAF couplings); prefix-free
  code with coarse-to-fine parameter grids; canonical dedup (type relabeling under each rule's implicit semantics:
  cyclic / empty-0 / parity / fixed). Round-trip verified on all programs <= 17 bits. Count grows ~2.1x per bit
  (<= 17 bits: 8,746). Textbook models: voter/majority/Schelling 11, Nowak-May 12, cyclic CA / GH 13, Vicsek 14,
  sandpile 15, MIPS 16, RPS agents 17, co-evolving voter 19, Langton ants 22, chemotaxis 28, Gray-Scott 29 bits.
- **T3 done:** vectorized interpreter (3 seeds batched); agents via periodic k-d tree (40 s -> 2-3 s per program);
  BTW with open boundary; quorum = slow-down (stop-rule jammed into an absorbing state). Per-key recording
  (grid/agents every 2 steps, fields/rewiring graphs every 10); runs never stopped early (absorbed state is data).
- **T4/T5 findings (battery as known-filter):**
  - Detectors are tuned to canonical scales/run lengths ("run length < 5 tau", "< 5 T_cross", "post-burn-in < 300").
    Thinning frames broke them -> battery now gets full resolution.
  - `powerlaw` was missing from the EPC venv -> P14 (SOC) could never fire anywhere. Installed (declared dependency).
  - Recognized so far: voter p=0.3/0.1 -> P18, Schelling -> P1 (definitive). Not recognized: synchronous voter
    (misses P18 screening narrowly), majority coarsening, cyclic-CA spirals (P12 is built for May-Leonard with empty
    sites), GH at k=3 (degenerate; canonical GH has 8 states), Vicsek at our density/speed, MIPS (dilute world),
    Kuramoto at K=0.1 (too weak), co-evolving voter (P34 gates), Nowak-May (P27 needs coop_fraction).
  - Battery FALSE MATCH seen: `C kC=3 | C.MAJ(p=0.1)` -> P12 (cyclic dominance) at confirmation.
  - Cost: per-program battery ~37 s mean (P1 alone 70 s at 999 permutations) -> infeasible at census scale.
- **Design revision proposed (needs owner sign-off before freeze): CLUSTER-FIRST pipeline.** Every program: sim +
  screen + cheap fingerprint (~2.2 s incl. sim, 3 seeds) -> cluster emergent views into behaviour classes per
  substrate family -> run the validated battery (full frames) on each class's 3 shortest representatives -> label
  classes; unlabelled classes -> literature check. Emergence-complexity table = shortest program per labelled class.
  Rationale: battery cost per class not per program; robust to the battery's OOD recall gaps (classes get labelled
  from their clearest members + literature). Per-program battery kept as a reference mode (`runner --battery`).
- **Throughput (cluster-first, shortest 40 programs):** sim 1.8 s + filter 0.45 s per program (3 seeds). Network
  screen was 8 s (generic_emergence recomputing modularity on static graphs) -> screen static networks on node
  states only. Re-measure on the full mix before choosing L (first estimate: L ~ 24 bits).
- **Running now:** reference slice (per-program battery, <= 13 bits, 430 programs) + cluster-first slice (<= 16
  bits, 4,134 programs). Next: triage on cf16 (clusters + battery on representatives), compare labels with the
  reference slice, null base rate (`runner --null`, same programs), G2 counts/benchmark lengths on the VPS.

**2026-09-28 — owner decisions**
- Cluster-first pipeline APPROVED ("the grouping approach makes sense").
- VALIDATION GATE required before any scale-up (owner: "even 1% wrong over a million = disaster"). v1 FAILED (naming
  2/10; one circularity in my negative design found and fixed in v2). v2 = non-circular held-out negatives (>= 300)
  + positive families.
- NAMING = REFERENCE LIBRARY APPROVED ("a better detector, using our catalog as examples"): run each known catalog
  phenomenon many times in the census world, record fingerprints, name groups by nearest reference within an
  acceptance radius; error rates measured on held-out examples + negatives. Catalog detectors kept as a second
  opinion where they work (P1, P14, P18 confirmed).
- ADDED TO PLAN: census finds of catalog patterns shown on the website (per catalog page: simplest census rule
  that produces it + animation). Deploy only with owner OK.
- ADDED (owner, 2026-09-28): record EVERY census program that produces a catalogued behaviour, not just the simplest
  -> per catalog pattern, a "routes" page: all distinct mechanisms (rule combinations, substrates) that produce it,
  which ingredients all routes share (necessary) vs which vary. Unexpected routes to a known behaviour get their own
  literature check (a new mechanism for a known behaviour can itself be a finding).
- **Anti-overfitting rule (2026-09-28):** features/namer are being tuned on the validation set (rounds 1-7). Once a
  round passes, a CONFIRMATION run on fresh, never-seen data (new seed sets for every library variant + held-out
  negatives with new seeds) must also pass, with identical code and criteria, before anything is frozen.

**2026-09-28 — VALIDATION GATE PASSED** (round 9 dev + fresh-data confirmation; code 25a99b7, reports committed).
Confirmation: recall 477/480, 0 false flags (24 interacting-trivial + 400 interaction-free), 0 names not shown,
0 names surviving knock-out, primary exact 86%, novelty-risk (withheld class given a primary name) 0/477.
Open before scaling: (1) expand interacting-but-trivial negatives to ~300 (0/24 only bounds the rate at ~12%);
(2) merge names that share a textbook measure into one behaviour name (e.g. "domain coarsening"), mechanism = route;
(3) re-measure throughput with the knock-out (+1 sim per flagged program) -> choose L; (4) freeze doc + OSF.
**2026-09-28 (later) — GATE PASSED ON FRESH DATA WITH EXTENDED NEGATIVES (confirm3; code 4d0b566).**
Path: confirmation 1 passed but only 24 interacting-trivial negatives -> expanded to 312 -> 1/312 then 20/312 flagged
(noise-sensitive / side-effect knock-out comparisons) -> knock-out redefined: the DETECTED EMERGENCE SIGNAL must drop
(seed-averaged evidence vs spread) OR a textbook ORDER measure must be clearly higher with interaction -> confirm3:
0/312 + 0/400 flagged (false-flag rate <= ~1%), recall 468/472, 0 unshown names, 0 names surviving knock-out, 88%
exact primary, novelty-risk 1/468 (genuinely co-present). Each rule change was followed by a fresh-seed re-run.
NEXT: (1) merge names sharing a textbook measure into behaviour names; (2) throughput with knock-out -> L;
(3) freeze doc (PREREG_census_v1.md) for owner approval -> OSF; (4) Phase 3 run.

**2026-10-01 — items 1-2 done (owner: "yes")**
- BEHAVIOUR NAMES (census/namer.py): library classes sharing a textbook measure merged for display ("domain
  coarsening" = voter + majority); mechanism = the program's rule set (route). Display layer only — confirm3 report
  regenerated IDENTICALLY (criteria + per-class table). Runner now runs the validated pipeline (screen + knock-out +
  textbook measures; program-list cache); triage builds the behaviour map with ROUTES; digest updated.
- THROUGHPUT (600 random programs <= 22 bits, 5 workers, validated pipeline): 11.9 s/program (A 24, C 16, AF/CAF 13,
  CA 12, CF 9, N 5, F 3), 50% flagged, 0 errors, 0 fingerprint errors. Mix at this length is 61% agent-containing.
  Budget 10.5 days x 5 workers = 1,260 worker-hours -> **L = 22 fits (289,798 programs, 986 worker-hours, ~8.2 days);
  L = 23 needs ~15.6 days.** G2 at L2 = 21 (98,082 programs) ~2.7 days. Enumeration to 24 bits takes ~1.7 h -> cached.
- First 22-bit sample map: 9 known behaviours with multiple routes each (domain coarsening 19 routes, flocking 9,
  lattice phase locking 10); 153 unnamed views in 15 classes (agent clustering, chemotaxis, ... -> literature check /
  library expansion); 64 known-behaviour-atypical-look (review).
- VPS disk hit 100% on 2026-10-01 (not the census: 20 MB); owner cleared savepoints -> 88%.
- BEFORE FREEZE: final gate re-run on the frozen code (regression), novelty arm + G2 runner on the validated
  pipeline, library choice for the census, PREREG_census_v1.md for owner approval -> OSF.
