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
