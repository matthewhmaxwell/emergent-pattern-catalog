EPC Emergence Census — pre-registration v1 (Matt Maxwell)

We list every simple rule-system up to a fixed size, run each one, and ask three questions in a fixed order:
(1) Did anything emerge? A pattern counts only if it disappears when the system's interactions are removed
(an interaction knock-out on the same seeds). (2) Is it a behaviour we already know? It is named only if it looks
like our reference examples of that behaviour AND passes that behaviour's textbook measure; otherwise it stays
unnamed. (3) If unnamed, is it known to science? An independent literature check decides; only a behaviour the
literature is silent about becomes a novelty candidate, which must then survive 10-seed replication at two system
sizes, rule knock-outs, null and symmetric controls, a second literature pass, and owner review.

Search space (frozen): grammar G1 — 31 mechanism-level rule templates on lattice cells, moving agents, network
nodes and diffusing fields, and their couplings, with an exact prefix-free length in bits. Exhaustive census of
every G1 program of length <= 22 bits (289,798 programs, shortest first). Second grammar G2 (finer parts: subject,
condition, action, rate), exhaustive to <= 21 bits (98,082 programs), to test whether "simplest" depends on the
grammar. Novelty arm: CVT-MAP-Elites over coupled G1 programs, 150,000 evaluations or 20,000 without a new niche.
Simulation: 64x64 worlds, 400 agents, 200-node networks, 1,000 steps, 3 seeds per program, seeds = hash of the
program's bits.

Promised results whatever happens: a map of known behaviours (the simplest program producing each, in both
grammars, every other route to it, robustness and ingredient knock-outs), Spearman correlation of G1 vs G2
minimal lengths, census statistics, every unnamed group with its literature verdict, every candidate with its
vetting outcome, a post-run false-flag audit (placebo test on 3,000 census programs), and all deviations.
A discovery is possible but not promised.

Validation before registration (frozen code, fresh seeds): gate passed on four separate seed sets — 474/476
real examples flagged, 0/1,243 trivial interacting programs flagged, 0 names the run does not show, 0/474 named
behaviours surviving the knock-out, 87% exact primary name (the rest sent to review, none wrongly named);
placebo test (interaction-free runs compared with interaction-free runs, any flag false by construction):
0/6,642 views. Failures found and fixed on the way are listed in the document (section 6).

Frozen artifacts: EPC repository, branch emergence-census, commit 52cee0f; SHA-256 of every census module, the
reference library and both program lists in CODE_census_v1.sha256; this document's hash in PREREG_census_v1.sha256.
Deviations are logged in docs/census/DEVIATIONS.md with date, reason, and whether made before or after seeing data.
