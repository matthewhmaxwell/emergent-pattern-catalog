# Per-program census records

`census_v1_g1/` (G1, every program of up to 22 bits, 289,798 programs) and `census_v1_g2/` (G2, every program of up to
21 bits, 98,082 programs) hold one JSON line per program, xz-compressed (`results_*.jsonl.xz`), as written by the
frozen census runner, plus the triage output (`triage.json.xz`: names, review queue, unnamed groups and the behavior
map) and the run digests (`DIGEST.md`, `TRIAGE.md`, `manifest.json`). Decompress with `xz -dk <file>`.

Each record gives the program (`bits`, `len`, `prog`, `layers`, `templates`, `grammar`), timing, whether the program
is interaction-free, and per view (`views`): `flagged`, `em_seeds` (seeds emergent under the screen), `driven`
(interaction-driven under the knock-out) and per-seed screen evidence. G2's `results_0..4` include 1,206 ERROR
records for the sandpile programs that crashed before Deviation 1; `results_fix.jsonl.xz` holds their re-runs on the
original seeds (see docs/census/DEVIATIONS.md).

`python -m census.paper_tables` recomputes the screen, knock-out and flag counts by grammar and by kind of view from
these files. Every run can also be regenerated exactly from its program, since the seeds are a hash of the program's
bit string.
