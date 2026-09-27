# Literature-check protocol for unlabelled behaviour classes (pilot draft v0; frozen at registration)

Used at triage step 4 (DESIGN §6) for every behaviour class the battery could not label. Same protocol as the
Ring-3 P1 literature pass (2026-09-27): an **independent** research agent (fresh context, web search, no access to
the census codebase or to our hopes for the result) receives one class card and must classify it with verified
citations. The agent is told explicitly that "known" is the expected, fine outcome.

## Input: the class card
- **Program in plain words** — the shortest program of the class, translated from G1 (e.g. "64x64 torus lattice,
  3 cell types; each step, every cell adopts the type held by the majority of its 8 neighbours with probability 0.1").
- **Also in the class** — up to 5 other member programs in plain words (to show what varies without changing the
  behaviour).
- **What was observed** — the fingerprint features that stand out (robust z > 2) translated to words (e.g. "large
  single-type domains, Moran's I rising over time, activity decaying"), the generic-lens kind, and the battery's
  near-misses (detectors that fired below confirmation).
- **Nearest catalogued patterns** — the 3 catalog entries whose reference fingerprints are closest.

## Instructions to the agent (verbatim)
> You are checking whether a behaviour produced by a very simple simulated rule system is already described in the
> scientific literature. Most such behaviours are known; "known" is the expected answer and a perfectly good one.
> 1. Identify the model family this rule system belongs to (statistical physics, cellular automata, active matter,
>    opinion dynamics, evolutionary game theory, synchronization, reaction–diffusion, network science, ...).
> 2. Search for the specific behaviour described (not just the model family): is this behaviour, in this or an
>    equivalent system, reported? Give the earliest clear source and one review if one exists.
> 3. Verify every citation (title, authors, year, venue) by opening it. Never cite from memory alone.
> 4. Answer with exactly one verdict:
>    - **KNOWN** — the behaviour is described for this or an equivalent system (cite).
>    - **KNOWN-VARIANT** — a close relative is described; say precisely what differs and whether the difference is
>      likely to matter.
>    - **LITERATURE-SILENT** — after a genuine search, nothing describes it; list the searches you ran.
> 5. Say how confident you are (low / medium / high) and why.

## Output handling
- **KNOWN / KNOWN-VARIANT and in the catalog's scope** -> catalog expansion (new detector + Phase-2a validation
  panel with negative controls, TNR 1.0) -> the known filter is re-run over every previously unlabelled class.
- **LITERATURE-SILENT** -> novelty candidate (Tier 3) -> vetting: replication (more seeds, two system sizes),
  ingredient knockout, interventional debunking (null and provably symmetric controls), a SECOND independent
  literature pass. Only then any claim.
- Every verdict, its citations and the searches are stored with the class id (`litcheck/<class-id>.json`).
